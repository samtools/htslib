#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "kstring.h"
#include "bgzf.h"
#include "vcf.h"

#include "khash.h"
KHASH_MAP_INIT_STR(vdict, vcf_keyinfo_t)
typedef khash_t(vdict) vdict_t;

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int vcf_verbose = 3; // 1: error; 2: warning; 3: message; 4: progress; 5: debugging; >=10: pure debugging
uint8_t vcf_type_size[] = { 0, 1, 2, 4, 8, 4, 8, 1, 1, 0, 1, 2, 4, 1, 0, 0 };
static vcf_keyinfo_t vcf_keyinfo_def = { { 15, 15, 15 }, -1, -1, -1, -1 };

/******************
 * Basic routines *
 ******************/

static inline void align_mem(kstring_t *s)
{
	if (s->l&7) {
		uint64_t zero = 0;
		int l = ((s->l + 7)>>3<<3) - s->l;
		kputsn((char*)&zero, l, s);
	}
}

/*************
 * Basic I/O *
 *************/

vcfFile *vcf_open(const char *fn, const char *mode, const char *fn_ref)
{
	const char *p;
	vcfFile *fp;
	fp = calloc(1, sizeof(vcfFile));
	for (p = mode; *p; ++p) {
		if (*p == 'w') fp->is_write = 1;
		else if (*p == 'b') fp->is_bin = 1;
	}
	if (fp->is_bin) {
		if (fp->is_write) fp->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_dopen(fileno(stdout), mode);
		else fp->fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
	} else {
		if (fp->is_write) {
			fp->fp = strcmp(fn, "-")? fopen(fn, "rb") : stdout;
		} else {
			gzFile gzfp;
			gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
			if (gzfp) fp->fp = ks_init(gzfp);
			if (fn_ref) fp->fn_ref = strdup(fn_ref);
		}
	}
	if (fp->fp == 0) {
		if (vcf_verbose >= 2)
			fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);
		free(fp->fn_ref); free(fp);
		return 0;
	}
	return fp;
}

void vcf_close(vcfFile *fp)
{
	if (!fp->is_bin) {
		free(fp->line.s);
		if (!fp->is_write) {
			gzFile gzfp = ((kstream_t*)fp->fp)->f;
			ks_destroy(fp->fp);
			gzclose(gzfp);
			free(fp->fn_ref);
		} else fclose(fp->fp);
	} else bgzf_close(fp->fp);
	free(fp);
}

/*********************
 * VCF header parser *
 *********************/

// return: positive => contig; zero => INFO/FILTER/FORMAT; negative => error or skipped
int vcf_hdr_parse_line2(const char *str, uint32_t *info, int *id_beg, int *id_end)
{
	const char *p, *q;
	int ctype;
	int type = -1; // Type
	int num = -1; // Number
	int var = -1; // A, G, ., or fixed
	int ctg_len = -1;

	if (*str != '#' && str[1] != '#') return -1;
	*id_beg = *id_end = *info = -1;
	p = str + 2;
	for (q = p; *q && *q != '='; ++q); // FIXME: do we need to check spaces?
	if (*q == 0) return -2;
	if (q - p == 4 && strncmp(p, "INFO", 4) == 0) ctype = VCF_DT_INFO;
	else if (q - p == 6 && strncmp(p, "FILTER", 6) == 0) ctype = VCF_DT_FLT;
	else if (q - p == 6 && strncmp(p, "FORMAT", 6) == 0) ctype = VCF_DT_FMT;
	else if (q - p == 6 && strncmp(p, "contig", 6) == 0) ctype = VCF_DT_CTG;
	else return -3;
	for (; *q && *q != '<'; ++q);
	if (*q == 0) return -3;
	p = q + 1; // now p points to the first character following '<'
	while (*p && *p != '>') {
		int which = 0;
		char *tmp;
		const char *val;
		for (q = p; *q && *q != '='; ++q);
		if (*q == 0) break;
		if (q - p == 2 && strncmp(p, "ID", 2) == 0) which = 1; // ID
		else if (q - p == 4 && strncmp(p, "Type", 4) == 0) which = 2; // Number
		else if (q - p == 6 && strncmp(p, "Number", 6) == 0) which = 3; // Type
		else if (q - p == 6 && strncmp(p, "length", 6) == 0) which = 4; // length
		val = q + 1;
		if (*val == '"') { // quoted string
			for (q = val + 1; *q && *q != '"'; ++q)
				if (*q == '\\' && *(q+1) != 0) ++q;
			if (*q != '"') return -4; // open double quotation mark
			p = q + 1;
			if (*p == ',') ++p;
			continue;
		}
		for (q = val; *q && *q != ',' && *q != '>'; ++q); // parse val
		if (which == 1) {
			*id_beg = val - str; *id_end = q - str;
		} else if (which == 2) {
			if (q - val == 7 && strncmp(val, "Integer", 7) == 0) type = VCF_TP_INT;
			else if (q - val == 5 && strncmp(val, "Float", 5) == 0) type = VCF_TP_REAL;
			else if (q - val == 6 && strncmp(val, "String", 6) == 0) type = VCF_TP_STR;
		} else if (which == 3) {
			if (*val == 'A') var = VCF_VTP_A;
			else if (*val == 'G') var = VCF_VTP_G;
			else if (isdigit(*val)) var = VCF_VTP_FIXED, num = strtol(val, &tmp, 10);
			else var = VCF_VTP_VAR;
			if (var != VCF_VTP_FIXED) num = 0xfffff;
		} else if (which == 4) {
			if (isdigit(*val)) ctg_len = strtol(val, &tmp, 10);
		}
		p = q + 1;
	}
	if (ctype == VCF_DT_CTG) {
		if (ctg_len > 0) return ctg_len;
		else return -5;
	} else {
		if (ctype == VCF_DT_FLT) num = 0;
		if (type == VCF_TP_FLAG) {
			if (num != 0 && vcf_verbose >= 2)
				fprintf(stderr, "[W::%s] ignore NUmber for a Flag\n", __func__);
			num = 0, var = VCF_VTP_FIXED; // if Flag VCF type, force to change num to 0
		}
		if (num == 0) type = VCF_TP_FLAG, var = VCF_VTP_FIXED; // conversely, if num==0, force the type to Flag
		if (*id_beg < 0 || type < 0 || num < 0 || var < 0) return -5; // missing information
		*info = (uint32_t)num<<12 | var<<8 | type<<4 | ctype;
		//printf("%d, %s, %d, %d, [%d,%d]\n", ctype, vcf_type_name[type], var, num, *id_beg, *id_end);
		return 0;
	}
}

vcf_hdr_t *vcf_hdr_init(void)
{
	vcf_hdr_t *h;
	h = calloc(1, sizeof(vcf_hdr_t));
	h->dict = kh_init(vdict);
	return h;
}

void vcf_hdr_destroy(vcf_hdr_t *h)
{
	khint_t k;
	vdict_t *d = (vdict_t*)h->dict;
	for (k = kh_begin(d); k != kh_end(d); ++k)
		if (kh_exist(d, k)) free((char*)kh_key(d, k));
	kh_destroy(vdict, d);
	free(h->mem.s); free(h->text); free(h->key); free(h->r2k); free(h->s2k);
	free(h);
}

int vcf_hdr_parse1(vcf_hdr_t *h, const char *str)
{
	khint_t k;
	vdict_t *d = (vdict_t*)h->dict;
	if (*str != '#') return -1;
	if (str[1] == '#') {
		uint32_t info;
		int len, ret, id_beg, id_end;
		char *s;

		len = vcf_hdr_parse_line2(str, &info, &id_beg, &id_end);
		if (len < 0) return -1;
		s = malloc(id_end - id_beg + 1);
		strncpy(s, str + id_beg, id_end - id_beg);
		s[id_end - id_beg] = 0;
		k = kh_put(vdict, d, s, &ret);
		if (ret == 0) { // present
			if (len > 0) {
				if (kh_val(d, k).rlen > 0) {
					if (vcf_verbose >= 2) // check contigs with identical name
						fprintf(stderr, "[W::%s] Duplicated contig name '%s'. Skipped.\n", __func__, s);
				} else {
					kh_val(d, k).rid = h->n_ref++;
					kh_val(d, k).rlen = len;
				}
			} else kh_val(d, k).info[info&0xf] = info;
			free(s);
		} else { // absent
			kh_val(d, k) = vcf_keyinfo_def; // set default values
			if (len > 0) { // contig
				kh_val(d, k).rid = h->n_ref++;
				kh_val(d, k).rlen = len;
			} else kh_val(d, k).info[info&0xf] = info;
			kh_val(d, k).kid = h->n_key++;
		}
	} else {
		int i = 0;
		const char *p, *q;
		// check if "PASS" is in the dictionary
		k = kh_get(vdict, d, "PASS");
		if (k == kh_end(d)) vcf_hdr_parse1(h, "##FILTER=<ID=PASS>"); // if not, add it
		// add samples
		for (p = q = str;; ++q) {
			int ret;
			if (*q != '\t' && *q != 0) continue;
			if (++i > 9) {
				char *s;
				s = malloc(q - p + 1);
				strncpy(s, p, q - p);
				s[q - p] = 0;
				k = kh_put(vdict, d, s, &ret);
				if (ret == 0) { // present
					if (kh_val(d, k).sid >= 0) {
						if (vcf_verbose >= 2)
							fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
					} else kh_val(d, k).sid = h->n_sample++;
					free(s);
				} else { // absent
					kh_val(d, k) = vcf_keyinfo_def;
					kh_val(d, k).sid = h->n_sample++;
					kh_val(d, k).kid = h->n_key++;
				}
			}
			if (*q == 0) break;
			p = q + 1;
		}
	}
	return 0;
}

int vcf_hdr_sync(vcf_hdr_t *h)
{
	khint_t k;
	vdict_t *d = (vdict_t*)h->dict;
	h->key = malloc(h->n_key * sizeof(vcf_keypair_t));
	h->r2k = malloc(h->n_ref * sizeof(int));
	h->s2k = malloc(h->n_sample * sizeof(int));
	for (k = kh_begin(d); k != kh_end(d); ++k) {
		int i;
		if (!kh_exist(d, k)) continue;
		i = kh_val(d, k).kid;
		assert(i < h->n_key);
		h->key[i].key = kh_key(d, k);
		h->key[i].info = &kh_val(d, k);
		if (kh_val(d, k).rid >= 0) h->r2k[kh_val(d, k).rid] = i;
		if (kh_val(d, k).sid >= 0) h->s2k[kh_val(d, k).sid] = i;
	}
	return 0;
}

int vcf_hdr_parse(vcf_hdr_t *h)
{
	char *p, *q;
	for (p = q = h->text;; ++q) {
		int c;
		if (*q != '\n' && *q != 0) continue;
		c = *q; *q = 0;
		vcf_hdr_parse1(h, p);
		*q = c;
		if (*q == 0) break;
		p = q + 1;
	}
	vcf_hdr_sync(h);
	return 0;
}

/******************
 * VCF header I/O *
 ******************/

vcf_hdr_t *vcf_hdr_read(vcfFile *fp)
{
	vcf_hdr_t *h;
	if (fp->is_write) return 0;
	h = vcf_hdr_init();
	if (fp->is_bin) {
		uint8_t magic[4];
		bgzf_read(fp->fp, magic, 4);
		bgzf_read(fp->fp, &h->l_text, 4);
		h->text = malloc(h->l_text);
		bgzf_read(fp->fp, h->text, h->l_text);
	} else {
		int dret;
		kstring_t txt, *s = &fp->line;
		txt.l = txt.m = 0; txt.s = 0;
		if (fp->fn_ref) {
			gzFile f;
			kstream_t *ks;
			f = gzopen(fp->fn_ref, "r");
			ks = ks_init(f);
			while (ks_getuntil(ks, 0, s, &dret) >= 0) {
				int c;
				txt.l = 0;
				kputs("##contig=<ID=", &txt); kputs(s->s, &txt);
				ks_getuntil(ks, 0, s, &dret);
				kputs(",length=", &txt); kputw(atol(s->s), &txt); kputc('>', &txt);
				vcf_hdr_parse1(h, txt.s);
				if (dret != '\n')
					while ((c = ks_getc(ks)) != '\n' && c != -1);
			}
			ks_destroy(ks);
			gzclose(f);
		}
		txt.l = 0;
		while (ks_getuntil(fp->fp, KS_SEP_LINE, s, &dret) >= 0) {
			if (s->l == 0) continue;
			if (s->s[0] != '#') {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[E::%s] no sample line\n", __func__);
				free(txt.s);
				vcf_hdr_destroy(h);
				return 0;
			}
			kputsn(s->s, s->l, &txt);
			if (s->s[1] != '#') break;
			kputc('\n', &txt);
		}
		h->l_text = txt.l + 1; // including NULL
		h->text = txt.s;
	}
	if (vcf_hdr_parse(h) != 0) { // FIXME: vcf_hdr_parse() always returns 0
		vcf_hdr_destroy(h);
		return 0;
	} else return h;
}

void vcf_hdr_write(vcfFile *fp, const vcf_hdr_t *h)
{
	if (fp->is_bin) {
		bgzf_write(fp->fp, "BCF\2", 4);
		bgzf_write(fp->fp, &h->l_text, 4);
		bgzf_write(fp->fp, h->text, h->l_text);
	} else {
		fwrite(h->text, 1, h->l_text, fp->fp);
		fputc('\n', fp->fp);
	}
}

/*******************
 * Typed value I/O *
 *******************/

void vcf_enc_int(kstring_t *s, int n, int32_t *a, int wsize)
{
	int32_t max = INT32_MIN + 1, min = INT32_MAX;
	int i;
	if (n == 0) vcf_enc_size(s, 0, VCF_RT_INT8);
	else if (n == 1) vcf_enc_int1(s, a[0]);
	else {
		if (wsize <= 0) wsize = n;
		for (i = 0; i < n; ++i) {
			if (a[i] == INT32_MIN) continue;
			if (max < a[i]) max = a[i];
			if (min > a[i]) min = a[i];
		}
		if (max <= INT8_MAX && min > INT8_MIN) {
			vcf_enc_size(s, wsize, VCF_RT_INT8);
			for (i = 0; i < n; ++i)
				kputc(a[i] == INT32_MIN? INT8_MIN : a[i], s);
		} else if (max <= INT16_MAX && min > INT16_MIN) {
			vcf_enc_size(s, wsize, VCF_RT_INT16);
			for (i = 0; i < n; ++i) {
				int16_t x = a[i] == INT32_MIN? INT16_MIN : a[i];
				kputsn((char*)&x, 2, s);
			}
		} else {
			vcf_enc_size(s, wsize, VCF_RT_INT32);
			for (i = 0; i < n; ++i) {
				int32_t x = a[i] == INT32_MIN? INT32_MIN : a[i];
				kputsn((char*)&x, 4, s);
			}
		}
	}
}

void vcf_enc_float(kstring_t *s, int n, float *a)
{
	vcf_enc_size(s, n, VCF_RT_FLOAT);
	kputsn((char*)a, 4 * n, s);
}

void vcf_fmt_array(kstring_t *s, int n, int type, void *data)
{
	int j = 0;
	if (type == VCF_RT_INT8) {
		int8_t *p = (int8_t*)data;
		for (j = 0; j < n && *p != INT8_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == VCF_RT_INT16) {
		int16_t *p = (int16_t*)data;
		for (j = 0; j < n && *p != INT16_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == VCF_RT_INT32) {
		int32_t *p = (int32_t*)data;
		for (j = 0; j < n && *p != INT32_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == VCF_RT_FLOAT) {
		float *p = (float*)data;
		for (j = 0; j < n && *(int32_t*)p != 0x7F800001; ++j, ++p) {
			if (j) kputc(',', s);
			ksprintf(s, "%g", *p);
		}
	} else if (type == VCF_RT_CHAR) {
		char *p = (char*)data;
		for (j = 0; j < n && *p; ++j, ++p)
			kputc(*p, s);
	}
	if (n && j == 0) kputc('.', s);
}

/****************************
 * Parsing VCF record lines *
 ****************************/

vcf1_t *vcf_init1()
{
	vcf1_t *v;
	v = calloc(1, sizeof(vcf1_t));
	return v;
}

void vcf_destroy1(vcf1_t *v)
{
	free(v->shared.s); free(v->indiv.s);
	free(v);
}

typedef struct {
	int key, size;
	int max_l, max_m, offset;
	uint32_t y;
	uint8_t *buf;
} fmt_aux_t;

int vcf_parse1(kstring_t *s, const vcf_hdr_t *h, vcf1_t *v)
{
	int i = 0;
	char *p, *q, *r, *t;
	fmt_aux_t *fmt;
	vdict_t *d = (vdict_t*)h->dict;
	kstring_t *str, *mem = (kstring_t*)&h->mem;
	khint_t k;
	ks_tokaux_t aux;
	uint16_t n_fmt = 0;

	mem->l = v->shared.l = v->indiv.l = 0;
	str = &v->shared;
	for (p = kstrtok(s->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
		q = (char*)aux.p;
		*q = 0;
		if (i == 0) { // CHROM
			k = kh_get(vdict, d, p);
			if (k == kh_end(d) || kh_val(d, k).rid < 0) {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[W::%s] can't find '%s' in the sequence dictionary\n", __func__, p);
				return 0;
			} else v->rid = kh_val(d, k).rid;
		} else if (i == 1) { // POS
			v->pos = atoi(p) - 1;
		} else if (i == 2) { // ID
			if (strcmp(p, ".")) kputsn(p, q - p, str);
			kputc(0, str);
		} else if (i == 3) { // REF
			kputsn(p, q - p + 1, str); // +1 to include NULL
		} else if (i == 4) { // ALT
			uint16_t n_alt;
			if (strcmp(p, ".")) {
				for (r = p, n_alt = 1; r != q; ++r)
					if (*r == ',') *r = 0, ++n_alt;
			} else n_alt = 0;
			kputsn((char*)&n_alt, 2, str);
			if (n_alt) kputsn(p, q - p + 1, str);
		} else if (i == 5) { // QUAL
			if (strcmp(p, ".")) v->qual = atof(p);
			else *(uint32_t*)&v->qual = 0x7F800001;
		} else if (i == 6) { // FILTER
			if (strcmp(p, ".")) {
				int32_t *a;
				int n_flt = 1, i;
				ks_tokaux_t aux1;
				// count the number of filters
				if (*(q-1) == ';') *(q-1) = 0;
				for (r = p; *r; ++r)
					if (*r == ';') ++n_flt;
				a = alloca(n_flt * 4);
				// add filters
				for (t = kstrtok(p, ";", &aux1), i = 0; t; t = kstrtok(0, 0, &aux1)) {
					*(char*)aux1.p = 0;
					k = kh_get(vdict, d, t);
					if (k == kh_end(d)) { // not defined
						if (vcf_verbose >= 2) fprintf(stderr, "[W::%s] undefined FILTER '%s'\n", __func__, t);
					} else a[i++] = kh_val(d, k).kid;
				}
				n_flt = i;
				vcf_enc_int(str, n_flt, a, -1);
			} else vcf_enc_int(str, 0, 0, -1);
		} else if (i == 7) { // INFO
			char *key;
			int n_info = 0, o_info = str->l;
			kputsn("\0\0", 2, str); // place holder for n_info
			if (strcmp(p, ".")) {
				if (*(q-1) == ';') *(q-1) = 0;
				for (r = key = p;; ++r) {
					int c;
					char *val, *end;
					if (*r != ';' && *r != '=' && *r != 0) continue;
					val = end = 0;
					c = *r; *r = 0;
					if (c == '=') {
						val = r + 1;
						for (end = val; *end != ';' && *end != 0; ++end);
						c = *end; *end = 0;
					} else end = r;
					k = kh_get(vdict, d, key);
					if (k == kh_end(d) || kh_val(d, k).info[VCF_DT_INFO] == 15) { // not defined in the header
						if (vcf_verbose >= 2) fprintf(stderr, "[W::%s] undefined INFO '%s'\n", __func__, key);
					} else { // defined in the header
						uint32_t y = kh_val(d, k).info[VCF_DT_INFO];
						if ((y>>4&0xf) == VCF_TP_FLAG || val == 0) { // a flag defined in the dict or without value
							vcf_enc_int1(str, kh_val(d, k).kid);
							++n_info;
							if (val != 0 && vcf_verbose >= 2)
								fprintf(stderr, "[W::%s] INFO '%s' is defined as a flag in the header but has a value '%s' in VCF; value skipped\n", __func__, key, val);
							if (val == 0 && (y>>4&0xf) != VCF_TP_FLAG && vcf_verbose >= 2)
								fprintf(stderr, "[W::%s] INFO '%s' takes at least a value, but no value is found\n", __func__, key);
						} else if ((y>>4&0xf) == VCF_TP_STR) { // a string
							vcf_enc_int1(str, kh_val(d, k).kid);
							vcf_enc_size(str, 1, VCF_RT_CSTR);
							kputsn(val, end - val + 1, str); // +1 to include NULL
							++n_info;
						} else { // an integer or float value or array
							int i, n_val = 1;
							char *t;
							for (t = val; *t; ++t)
								if (*t == ',') ++n_val;
							++n_info;
							vcf_enc_int1(str, kh_val(d, k).kid);
							if ((y>>4&0xf) == VCF_TP_INT) {
								int32_t *z;
								z = alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtol(t, &t, 10);
								vcf_enc_int(str, n_val, z, -1);
							} else if ((y>>4&0xf) == VCF_TP_REAL) {
								float *z;
								z = alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtod(t, &t);
								vcf_enc_float(str, n_val, z);
							}
						}
					}
					if (c == 0) break;
					r = end;
					key = r + 1;
				}
			}
			*(uint16_t*)(str->s + o_info) = n_info;
		} else if (i == 8 && h->n_sample > 0) { // FORMAT
			int j, l, m;
			ks_tokaux_t aux1;
			// count the number of format fields
			for (r = p, n_fmt = 1; *r; ++r)
				if (*r == ':') ++n_fmt;
			fmt = alloca(n_fmt * sizeof(fmt_aux_t));
			// get format information from the dictionary
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				*(char*)aux1.p = 0;
				k = kh_get(vdict, d, t);
				if (k == kh_end(d) || kh_val(d, k).info[VCF_DT_FMT] == 15) {
					if (vcf_verbose >= 2)
						fprintf(stderr, "[W::%s] FORMAT '%s' is not defined in the header\n", __func__, t);
					n_fmt = 0;
					break;
				} else {
					fmt[j].max_l = fmt[j].max_m = 0;
					fmt[j].key = kh_val(d, k).kid;
					fmt[j].y = h->key[fmt[j].key].info->info[VCF_DT_FMT];
				}
			}
			// compute max
			for (r = q + 1, j = 0, m = l = 1;; ++r, ++l) {
				if (*r == ':' || *r == '\t' || *r == '\0') { // end of a sample
					if (fmt[j].max_m < m) fmt[j].max_m = m;
					if (fmt[j].max_l < l) fmt[j].max_l = l;
					l = 0, m = 1;
					if (*r != ':') j = 0;
					else ++j;
				} else if (*r == ',') ++m;
				if (*r == 0) break;
			}
			// allocate memory for arrays
			for (j = 0; j < n_fmt; ++j) {
				fmt_aux_t *f = &fmt[j];
				if ((f->y>>4&0xf) == VCF_TP_STR) {
					f->size = f->max_l;
				} else if ((f->y>>4&0xf) == VCF_TP_REAL || (f->y>>4&0xf) == VCF_TP_INT) {
					f->size = f->max_m << 2;
				} else abort(); // I do not know how to do with Flag in the genotype fields
				align_mem(mem);
				f->offset = mem->l;
				ks_resize(mem, mem->l + h->n_sample * f->size);
				mem->l += h->n_sample * f->size;
			}
			for (j = 0; j < n_fmt; ++j)
				fmt[j].buf = (uint8_t*)mem->s + fmt[j].offset;
		} else if (i >= 9 && h->n_sample > 0) {
			int j, l;
			ks_tokaux_t aux1;
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				fmt_aux_t *z = &fmt[j];
				*(char*)aux1.p = 0;
				if ((z->y>>4&0xf) == VCF_TP_STR) {
					r = (char*)z->buf + z->size * (i - 9);
					for (l = 0; l < aux1.p - t; ++l) r[l] = t[l];
					for (; l != z->size; ++l) r[l] = 0;
				} else if ((z->y>>4&0xf) == VCF_TP_INT) {
					int32_t *x = (int32_t*)(z->buf + z->size * (i - 9));
					for (r = t, l = 0; r < aux1.p; ++r) {
						if (*r == '.') x[l++] = INT32_MIN, ++r;
						else x[l++] = strtol(r, &r, 10);
					}
					for (; l != z->size>>2; ++l) x[l] = INT32_MIN;
				} else if ((z->y>>4&0xf) == VCF_TP_REAL) {
					float *x = (float*)(z->buf + z->size * (i - 9));
					for (r = t, l = 0; r < aux1.p; ++r) {
						if (*r == '.' && !isdigit(r[1])) *(int32_t*)&x[l++] = 0x7F800001;
						else x[l++] = strtod(r, &r);
					}
					for (; l != z->size>>2; ++l) *(int32_t*)(x+l) = 0x7F800001;
				}
			}
		}
	}
	if (h->n_sample > 0) { // write individual genotype information
		str = &v->indiv;
		kputsn((char*)&n_fmt, 2, str);
		for (i = 0; i < n_fmt; ++i) {
			fmt_aux_t *z = &fmt[i];
			vcf_enc_int1(str, z->key);
			if ((z->y>>4&0xf) == VCF_TP_STR) {
				vcf_enc_size(str, z->size, VCF_RT_CHAR);
				kputsn((char*)z->buf, z->size * h->n_sample, str);
			} else if ((z->y>>4&0xf) == VCF_TP_INT) {
				vcf_enc_int(str, (z->size>>2) * h->n_sample, (int32_t*)z->buf, z->size>>2);
			} else {
				vcf_enc_size(str, z->size>>2, VCF_RT_FLOAT);
				kputsn((char*)z->buf, z->size * h->n_sample, str);
			}
		}
	}
	return 0;
}

int vcf_read1(vcfFile *fp, const vcf_hdr_t *h, vcf1_t *v)
{
	if (fp->is_bin) {
	} else {
		int ret, dret;
		ret = ks_getuntil(fp->fp, KS_SEP_LINE, &fp->line, &dret);
		if (ret < 0) return -1;
		ret = vcf_parse1(&fp->line, h, v);
	}
	return 0;
}

/**************************
 * Print VCF record lines *
 **************************/

typedef struct {
	int key, type, n, size;
	uint8_t *p;
} fmt_daux_t;

int vcf_format1(const vcf_hdr_t *h, const vcf1_t *v, kstring_t *s)
{
	uint8_t *ptr = (uint8_t*)v->shared.s;
	int i, l;
	s->l = 0;
	kputs(h->key[h->r2k[v->rid]].key, s); kputc('\t', s); // CHROM
	kputw(v->pos + 1, s); kputc('\t', s); // POS
	if (*ptr) { // ID
		l = strlen((char*)ptr);
		kputsn((char*)ptr, l, s);
		kputc('\t', s);
		ptr += l;
	} else kputsn(".\t", 2, s);
	++ptr;
	if (*ptr) { // REF
		l = strlen((char*)ptr);
		kputsn((char*)ptr, l, s);
		kputc('\t', s);
		ptr += l;
	} else kputsn(".\t", 2, s);
	++ptr;
	l = *(uint16_t*)ptr; // ALT
	ptr += 2;
	if (l) { // n_alt != 0
		for (i = 0; i < l; ++i) {
			int t = strlen((char*)ptr);
			kputsn((char*)ptr, t, s);
			ptr += t + 1;
		}
		kputc('\t', s);
	} else kputsn(".\t", 2, s);
	if (*(uint32_t*)&v->qual == 0x7F800001) kputsn(".\t", 2, s); // QUAL
	else ksprintf(s, "%g\t", v->qual);
	if (*ptr>>4) { // FILTER
		int32_t x, y;
		int type, i;
		x = vcf_dec_size(ptr, &ptr, &type);
		for (i = 0; i < x; ++i) {
			if (i) kputc(';', s);
			y = vcf_dec_int1(ptr, type, &ptr);
			kputs(h->key[y].key, s);
		}
		kputc('\t', s);
	} else {
		kputsn(".\t", 2, s);
		++ptr;
	}
	l = *(uint16_t*)ptr; // INFO
	ptr += 2;
	if (*(uint16_t*)ptr) { // INFO
		int i, n_info = l, type;
		for (i = 0; i < n_info; ++i) {
			int32_t x;
			uint32_t y;
			if (i) kputc(';', s);
			x = vcf_dec_typed_int1(ptr, &ptr);
			kputs(h->key[x].key, s);
			y = h->key[x].info->info[VCF_DT_INFO];
			if ((y>>4&0xf) != VCF_TP_FLAG) {
				kputc('=', s);
				if ((y>>4&0xf) == VCF_TP_STR) {
					++ptr; // skip the typing byte
					l = strlen((char*)ptr);
					kputsn((char*)ptr, l, s);
					ptr += l + 1;
				} else {
					x = vcf_dec_size(ptr, &ptr, &type);
					vcf_fmt_array(s, x, type, ptr);
					ptr += vcf_type_size[type] * x;
				}
			} else continue;
		}
	} else kputc('.', s);
	// FORMAT and individual information
	ptr = (uint8_t*)v->indiv.s;
	l = *(uint16_t*)ptr;
	ptr += 2;
	if (h->n_sample && l) { // FORMAT
		int i, j, n_fmt = l;
		fmt_daux_t *faux;
		faux = alloca(n_fmt * sizeof(fmt_daux_t));
		kputc('\t', s);
		for (i = 0; i < n_fmt; ++i) {
			fmt_daux_t *f = &faux[i];
			f->key = vcf_dec_typed_int1(ptr, &ptr);
			f->n = vcf_dec_size(ptr, &ptr, &f->type);
			f->size = vcf_type_size[f->type] * f->n;
			f->p = ptr;
			ptr += h->n_sample * f->size;
			if (i) kputc(':', s);
			kputs(h->key[f->key].key, s);
		}
		for (j = 0; j < h->n_sample; ++j) {
			kputc('\t', s);
			for (i = 0; i < n_fmt; ++i) {
				fmt_daux_t *f = &faux[i];
				if (i) kputc(':', s);
				vcf_fmt_array(s, f->n, f->type, f->p + j * f->size);
			}
		}
	}
	return 0;
}

int vcf_write1(vcfFile *fp, const vcf_hdr_t *h, const vcf1_t *v)
{
	if (fp->is_bin) {
	} else {
		vcf_format1(h, v, &fp->line);
		fwrite(fp->line.s, 1, fp->line.l, fp->fp);
		fputc('\n', fp->fp);
	}
	return 0;
}

int main(int argc, char *argv[])
{
	vcf_hdr_t *h;
	vcfFile *in, *out;
	vcf1_t *v;

	if (argc == 1) {
		fprintf(stderr, "Usage: bcf2ls <in.vcf>\n");
		return 1;
	}
	in = vcf_open(argv[1], "r", argc > 2? argv[2] : 0);
	out = vcf_open("-", "w", 0);
	h = vcf_hdr_read(in);
	vcf_hdr_write(out, h);
	v = vcf_init1();
	while (vcf_read1(in, h, v) >= 0) vcf_write1(out, h, v);
	vcf_destroy1(v);
	vcf_hdr_destroy(h);
	vcf_close(out);
	vcf_close(in);
	return 0;
}
