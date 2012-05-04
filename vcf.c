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
KHASH_MAP_INIT_STR(vdict, vcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int vcf_verbose = 3; // 1: error; 2: warning; 3: message; 4: progress; 5: debugging; >=10: pure debugging
uint32_t vcf_missing_float = 0x7F800001;
uint8_t vcf_type_shift[] = { 0, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static vcf_idinfo_t vcf_idinfo_def = { { 15, 15, 15 }, -1 };

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
	if (q - p == 4 && strncmp(p, "INFO", 4) == 0) ctype = VCF_HL_INFO;
	else if (q - p == 6 && strncmp(p, "FILTER", 6) == 0) ctype = VCF_HL_FLT;
	else if (q - p == 6 && strncmp(p, "FORMAT", 6) == 0) ctype = VCF_HL_FMT;
	else if (q - p == 6 && strncmp(p, "contig", 6) == 0) ctype = VCF_HL_CTG;
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
			if (q - val == 7 && strncmp(val, "Integer", 7) == 0) type = VCF_HT_INT;
			else if (q - val == 5 && strncmp(val, "Float", 5) == 0) type = VCF_HT_REAL;
			else if (q - val == 6 && strncmp(val, "String", 6) == 0) type = VCF_HT_STR;
			else if (q - val == 4 && strncmp(val, "Flag", 6) == 0) type = VCF_HT_FLAG;
		} else if (which == 3) {
			if (*val == 'A') var = VCF_VL_A;
			else if (*val == 'G') var = VCF_VL_G;
			else if (isdigit(*val)) var = VCF_VL_FIXED, num = strtol(val, &tmp, 10);
			else var = VCF_VL_VAR;
			if (var != VCF_VL_FIXED) num = 0xfffff;
		} else if (which == 4) {
			if (isdigit(*val)) ctg_len = strtol(val, &tmp, 10);
		}
		p = q + 1;
	}
	if (ctype == VCF_HL_CTG) {
		if (ctg_len > 0) return ctg_len;
		else return -5;
	} else {
		if (ctype == VCF_HL_FLT) num = 0;
		if (type == VCF_HT_FLAG) {
			if (num != 0 && vcf_verbose >= 2)
				fprintf(stderr, "[W::%s] ignore Number for a Flag\n", __func__);
			num = 0, var = VCF_VL_FIXED; // if Flag VCF type, force to change num to 0
		}
		if (num == 0) type = VCF_HT_FLAG, var = VCF_VL_FIXED; // conversely, if num==0, force the type to Flag
		if (*id_beg < 0 || type < 0 || num < 0 || var < 0) return -5; // missing information
		*info = (uint32_t)num<<12 | var<<8 | type<<4 | ctype;
		//printf("%d, %s, %d, %d, [%d,%d]\n", ctype, vcf_type_name[type], var, num, *id_beg, *id_end);
		return 0;
	}
}

vcf_hdr_t *vcf_hdr_init(void)
{
	int i;
	vcf_hdr_t *h;
	h = calloc(1, sizeof(vcf_hdr_t));
	for (i = 0; i < 3; ++i)
		h->dict[i] = kh_init(vdict);
	return h;
}

void vcf_hdr_destroy(vcf_hdr_t *h)
{
	int i;
	khint_t k;
	for (i = 0; i < 3; ++i) {
		vdict_t *d = (vdict_t*)h->dict[i];
		for (k = kh_begin(d); k != kh_end(d); ++k)
			if (kh_exist(d, k)) free((char*)kh_key(d, k));
		kh_destroy(vdict, d);
		free(h->id[i]);
	}
	free(h->mem.s); free(h->text);
	free(h);
}

int vcf_hdr_parse1(vcf_hdr_t *h, const char *str)
{
	khint_t k;
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
		if (len > 0) { // a contig line
			vdict_t *d = (vdict_t*)h->dict[VCF_DT_CTG];
			k = kh_put(vdict, d, s, &ret);
			if (ret == 0) { // a contig line
				if (vcf_verbose >= 2)
					fprintf(stderr, "[W::%s] Duplicated contig name '%s'. Skipped.\n", __func__, s);
			} else {
				kh_val(d, k) = vcf_idinfo_def;
				kh_val(d, k).id = kh_size(d) - 1;
				kh_val(d, k).info[0] = len;
			}
		} else { // a FILTER/INFO/FORMAT line
			vdict_t *d = (vdict_t*)h->dict[VCF_DT_ID];
			k = kh_put(vdict, d, s, &ret);
			if (ret) { // absent from the dict
				kh_val(d, k) = vcf_idinfo_def;
				kh_val(d, k).info[info&0xf] = info;
				kh_val(d, k).id = kh_size(d) - 1;
			} else kh_val(d, k).info[info&0xf] = info;
		}
	} else {
		int i = 0;
		const char *p, *q;
		vdict_t *d = (vdict_t*)h->dict[VCF_DT_ID];
		// check if "PASS" is in the dictionary
		k = kh_get(vdict, d, "PASS");
		if (k == kh_end(d)) vcf_hdr_parse1(h, "##FILTER=<ID=PASS>"); // if not, add it; this is a recursion
		// add samples
		d = (vdict_t*)h->dict[VCF_DT_SAMPLE];
		for (p = q = str;; ++q) {
			int ret;
			if (*q != '\t' && *q != 0) continue;
			if (++i > 9) {
				char *s;
				s = malloc(q - p + 1);
				strncpy(s, p, q - p);
				s[q - p] = 0;
				k = kh_put(vdict, d, s, &ret);
				if (ret) { // absent
					kh_val(d, k) = vcf_idinfo_def;
					kh_val(d, k).id = kh_size(d) - 1;
				} else {
					if (vcf_verbose >= 2)
						fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
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
	int i;
	for (i = 0; i < 3; ++i) {
		khint_t k;
		vdict_t *d = (vdict_t*)h->dict[i];
		h->n[i]  = kh_size(d);
		h->id[i] = malloc(kh_size(d) * sizeof(vcf_idpair_t));
		for (k = kh_begin(d); k != kh_end(d); ++k) {
			if (!kh_exist(d, k)) continue;
			h->id[i][kh_val(d, k).id].id   = kh_key(d, k);
			h->id[i][kh_val(d, k).id].info = &kh_val(d, k);
		}
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
		while (ks_getuntil(fp->fp, KS_SEP_LINE, s, &dret) >= 0) {
			if (s->l == 0) continue;
			if (s->s[0] != '#') {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[E::%s] no sample line\n", __func__);
				free(txt.s);
				vcf_hdr_destroy(h);
				return 0;
			}
			if (s->s[1] != '#' && fp->fn_ref) { // insert contigs here
				gzFile f;
				kstream_t *ks;
				kstring_t tmp;
				tmp.l = tmp.m = 0; tmp.s = 0;
				f = gzopen(fp->fn_ref, "r");
				ks = ks_init(f);
				while (ks_getuntil(ks, 0, &tmp, &dret) >= 0) {
					int c;
					kputs("##contig=<ID=", &txt); kputs(tmp.s, &txt);
					ks_getuntil(ks, 0, &tmp, &dret);
					kputs(",length=", &txt); kputw(atol(tmp.s), &txt);
					kputsn(">\n", 2, &txt);
					if (dret != '\n')
						while ((c = ks_getc(ks)) != '\n' && c != -1); // skip the rest of the line
				}
				free(tmp.s);
				ks_destroy(ks);
				gzclose(f);
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

void vcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
{
	int32_t max = INT32_MIN + 1, min = INT32_MAX;
	int i;
	if (n == 0) vcf_enc_size(s, 0, VCF_BT_NULL);
	else if (n == 1) vcf_enc_int1(s, a[0]);
	else {
		if (wsize <= 0) wsize = n;
		for (i = 0; i < n; ++i) {
			if (a[i] == INT32_MIN) continue;
			if (max < a[i]) max = a[i];
			if (min > a[i]) min = a[i];
		}
		if (max <= INT8_MAX && min > INT8_MIN) {
			vcf_enc_size(s, wsize, VCF_BT_INT8);
			for (i = 0; i < n; ++i)
				kputc(a[i] == INT32_MIN? INT8_MIN : a[i], s);
		} else if (max <= INT16_MAX && min > INT16_MIN) {
			vcf_enc_size(s, wsize, VCF_BT_INT16);
			for (i = 0; i < n; ++i) {
				int16_t x = a[i] == INT32_MIN? INT16_MIN : a[i];
				kputsn((char*)&x, 2, s);
			}
		} else {
			vcf_enc_size(s, wsize, VCF_BT_INT32);
			for (i = 0; i < n; ++i) {
				int32_t x = a[i] == INT32_MIN? INT32_MIN : a[i];
				kputsn((char*)&x, 4, s);
			}
		}
	}
}

void vcf_enc_vfloat(kstring_t *s, int n, float *a)
{
	vcf_enc_size(s, n, VCF_BT_FLOAT);
	kputsn((char*)a, n << 2, s);
}

void vcf_enc_vchar(kstring_t *s, int l, char *a)
{
	vcf_enc_size(s, l, VCF_BT_CHAR);
	kputsn(a, l, s);
}

void vcf_fmt_array(kstring_t *s, int n, int type, void *data)
{
	int j = 0;
	if (type == VCF_BT_INT8) {
		int8_t *p = (int8_t*)data;
		for (j = 0; j < n && *p != INT8_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == VCF_BT_CHAR) {
		char *p = (char*)data;
		for (j = 0; j < n && *p; ++j, ++p) kputc(*p, s);
	} else if (type == VCF_BT_INT32) {
		int32_t *p = (int32_t*)data;
		for (j = 0; j < n && *p != INT32_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == VCF_BT_FLOAT) {
		float *p = (float*)data;
		for (j = 0; j < n && *(int32_t*)p != 0x7F800001; ++j, ++p) {
			if (j) kputc(',', s);
			ksprintf(s, "%g", *p);
		}
	} else if (type == VCF_BT_INT16) {
		int16_t *p = (int16_t*)data;
		for (j = 0; j < n && *p != INT16_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
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

static inline void align_mem(kstring_t *s)
{
	if (s->l&7) {
		uint64_t zero = 0;
		int l = ((s->l + 7)>>3<<3) - s->l;
		kputsn((char*)&zero, l, s);
	}
}

int vcf_parse1(kstring_t *s, const vcf_hdr_t *h, vcf1_t *v)
{
	int i = 0;
	char *p, *q, *r, *t;
	fmt_aux_t *fmt = 0;
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
			vdict_t *d = (vdict_t*)h->dict[VCF_DT_CTG];
			k = kh_get(vdict, d, p);
			if (k == kh_end(d)) {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[W::%s] can't find '%s' in the sequence dictionary\n", __func__, p);
				return 0;
			} else v->rid = kh_val(d, k).id;
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
			else memcpy(&v->qual, &vcf_missing_float, 4);
		} else if (i == 6) { // FILTER
			if (strcmp(p, ".")) {
				int32_t *a;
				int n_flt = 1, i;
				ks_tokaux_t aux1;
				vdict_t *d = (vdict_t*)h->dict[VCF_DT_ID];
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
					} else a[i++] = kh_val(d, k).id;
				}
				n_flt = i;
				vcf_enc_vint(str, n_flt, a, -1);
			} else vcf_enc_vint(str, 0, 0, -1);
		} else if (i == 7) { // INFO
			char *key;
			int n_info = 0, o_info = str->l;
			vdict_t *d = (vdict_t*)h->dict[VCF_DT_ID];
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
					if (k == kh_end(d) || kh_val(d, k).info[VCF_HL_INFO] == 15) { // not defined in the header
						if (vcf_verbose >= 2) fprintf(stderr, "[W::%s] undefined INFO '%s'\n", __func__, key);
					} else { // defined in the header
						uint32_t y = kh_val(d, k).info[VCF_HL_INFO];
						++n_info;
						vcf_enc_int1(str, kh_val(d, k).id);
						if (val == 0) {
							vcf_enc_size(str, 0, VCF_BT_NULL);
						} else if ((y>>4&0xf) == VCF_HT_FLAG || (y>>4&0xf) == VCF_HT_STR) { // if Flag has a value, treat it as a string
							vcf_enc_vchar(str, end - val, val);
						} else { // int/float value/array
							int i, n_val;
							char *t;
							for (t = val, n_val = 1; *t; ++t) // count the number of values
								if (*t == ',') ++n_val;
							if ((y>>4&0xf) == VCF_HT_INT) {
								int32_t *z;
								z = alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtol(t, &t, 10);
								vcf_enc_vint(str, n_val, z, -1);
							} else if ((y>>4&0xf) == VCF_HT_REAL) {
								float *z;
								z = alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtod(t, &t);
								vcf_enc_vfloat(str, n_val, z);
							}
						}
					}
					if (c == 0) break;
					r = end;
					key = r + 1;
				}
			}
			*(uint16_t*)(str->s + o_info) = n_info;
		} else if (i == 8 && h->n[VCF_DT_SAMPLE] > 0) { // FORMAT
			int j, l, m;
			ks_tokaux_t aux1;
			vdict_t *d = (vdict_t*)h->dict[VCF_DT_ID];
			// count the number of format fields
			for (r = p, n_fmt = 1; *r; ++r)
				if (*r == ':') ++n_fmt;
			fmt = alloca(n_fmt * sizeof(fmt_aux_t));
			// get format information from the dictionary
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				*(char*)aux1.p = 0;
				k = kh_get(vdict, d, t);
				if (k == kh_end(d) || kh_val(d, k).info[VCF_HL_FMT] == 15) {
					if (vcf_verbose >= 2)
						fprintf(stderr, "[W::%s] FORMAT '%s' is not defined in the header\n", __func__, t);
					n_fmt = 0;
					break;
				} else {
					fmt[j].max_l = fmt[j].max_m = 0;
					fmt[j].key = kh_val(d, k).id;
					fmt[j].y = h->id[0][fmt[j].key].info->info[VCF_HL_FMT];
				}
			}
			// compute max
			for (r = q + 1, j = 0, m = l = 1;; ++r, ++l) {
				if (*r == ':' || *r == '\t' || *r == '\0') { // end of a sample
					if (fmt[j].max_m < m) fmt[j].max_m = m;
					if (fmt[j].max_l < l - 1) fmt[j].max_l = l - 1;
					l = 0, m = 1;
					if (*r != ':') j = 0;
					else ++j;
				} else if (*r == ',') ++m;
				if (*r == 0) break;
			}
			// allocate memory for arrays
			for (j = 0; j < n_fmt; ++j) {
				fmt_aux_t *f = &fmt[j];
				if ((f->y>>4&0xf) == VCF_HT_STR) {
					f->size = f->max_l;
				} else if ((f->y>>4&0xf) == VCF_HT_REAL || (f->y>>4&0xf) == VCF_HT_INT) {
					f->size = f->max_m << 2;
				} else abort(); // I do not know how to do with Flag in the genotype fields
				align_mem(mem);
				f->offset = mem->l;
				ks_resize(mem, mem->l + h->n[VCF_DT_SAMPLE] * f->size);
				mem->l += h->n[VCF_DT_SAMPLE] * f->size;
			}
			for (j = 0; j < n_fmt; ++j)
				fmt[j].buf = (uint8_t*)mem->s + fmt[j].offset;
		} else if (i >= 9 && h->n[VCF_DT_SAMPLE] > 0) {
			int j, l;
			ks_tokaux_t aux1;
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				fmt_aux_t *z = &fmt[j];
				*(char*)aux1.p = 0;
				if ((z->y>>4&0xf) == VCF_HT_STR) {
					r = (char*)z->buf + z->size * (i - 9);
					for (l = 0; l < aux1.p - t; ++l) r[l] = t[l];
					for (; l != z->size; ++l) r[l] = 0;
				} else if ((z->y>>4&0xf) == VCF_HT_INT) {
					int32_t *x = (int32_t*)(z->buf + z->size * (i - 9));
					for (r = t, l = 0; r < aux1.p; ++r) {
						if (*r == '.') x[l++] = INT32_MIN, ++r;
						else x[l++] = strtol(r, &r, 10);
					}
					for (; l != z->size>>2; ++l) x[l] = INT32_MIN;
				} else if ((z->y>>4&0xf) == VCF_HT_REAL) {
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
	// write individual genotype information
	str = &v->indiv;
	kputsn((char*)&n_fmt, 2, str);
	if (h->n[VCF_DT_SAMPLE] > 0) {
		for (i = 0; i < n_fmt; ++i) {
			fmt_aux_t *z = &fmt[i];
			vcf_enc_int1(str, z->key);
			if ((z->y>>4&0xf) == VCF_HT_STR) {
				vcf_enc_size(str, z->size, VCF_BT_CHAR);
				kputsn((char*)z->buf, z->size * h->n[VCF_DT_SAMPLE], str);
			} else if ((z->y>>4&0xf) == VCF_HT_INT) {
				vcf_enc_vint(str, (z->size>>2) * h->n[VCF_DT_SAMPLE], (int32_t*)z->buf, z->size>>2);
			} else {
				vcf_enc_size(str, z->size>>2, VCF_BT_FLOAT);
				kputsn((char*)z->buf, z->size * h->n[VCF_DT_SAMPLE], str);
			}
		}
	}
	return 0;
}

int vcf_read1(vcfFile *fp, const vcf_hdr_t *h, vcf1_t *v)
{
	if (fp->is_bin) {
		uint32_t x[4];
		int ret;
		if ((ret = bgzf_read(fp->fp, x, 16)) != 16) {
			if (ret == 0) return -1;
			return -2;
		}
		ks_resize(&v->shared, x[0] - 12);
		v->shared.l = x[0] - 12;
		v->rid = x[1]; v->pos = x[2];
		memcpy(&v->qual, &x[3], 4);
		bgzf_read(fp->fp, v->shared.s, v->shared.l);
		bgzf_read(fp->fp, x, 4);
		ks_resize(&v->indiv, x[0]);
		v->indiv.l = x[0];
		bgzf_read(fp->fp, v->indiv.s, v->indiv.l);
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
	kputs(h->id[VCF_DT_CTG][v->rid].id, s); kputc('\t', s); // CHROM
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
			if (i) kputc(',', s);
			kputsn((char*)ptr, t, s);
			ptr += t + 1;
		}
		kputc('\t', s);
	} else kputsn(".\t", 2, s);
	if (memcmp(&v->qual, &vcf_missing_float, 4) == 0) kputsn(".\t", 2, s); // QUAL
	else ksprintf(s, "%g\t", v->qual);
	if (*ptr>>4) { // FILTER
		int32_t x, y;
		int type, i;
		x = vcf_dec_size(ptr, &ptr, &type);
		for (i = 0; i < x; ++i) {
			if (i) kputc(';', s);
			y = vcf_dec_int1(ptr, type, &ptr);
			kputs(h->id[VCF_DT_ID][y].id, s);
		}
		kputc('\t', s);
	} else {
		kputsn(".\t", 2, s);
		++ptr;
	}
	l = *(uint16_t*)ptr; // INFO
	ptr += 2;
	if (l) {
		int i, n_info = l, type;
		for (i = 0; i < n_info; ++i) {
			int32_t x;
			if (i) kputc(';', s);
			x = vcf_dec_typed_int1(ptr, &ptr);
			kputs(h->id[VCF_DT_ID][x].id, s);
			if (*ptr>>4) { // more than zero element
				kputc('=', s);
				x = vcf_dec_size(ptr, &ptr, &type);
				vcf_fmt_array(s, x, type, ptr);
				ptr += x << vcf_type_shift[type];
			} else ++ptr;
		}
	} else kputc('.', s);
	// FORMAT and individual information
	ptr = (uint8_t*)v->indiv.s;
	l = *(uint16_t*)ptr;
	ptr += 2;
	if (h->n[VCF_DT_SAMPLE] && l) { // FORMAT
		int i, j, n_fmt = l;
		fmt_daux_t *faux;
		faux = alloca(n_fmt * sizeof(fmt_daux_t));
		kputc('\t', s);
		for (i = 0; i < n_fmt; ++i) {
			fmt_daux_t *f = &faux[i];
			f->key = vcf_dec_typed_int1(ptr, &ptr);
			f->n = vcf_dec_size(ptr, &ptr, &f->type);
			f->size = f->n << vcf_type_shift[f->type];
			f->p = ptr;
			ptr += h->n[VCF_DT_SAMPLE] * f->size;
			if (i) kputc(':', s);
			kputs(h->id[VCF_DT_ID][f->key].id, s);
		}
		for (j = 0; j < h->n[VCF_DT_SAMPLE]; ++j) {
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
		uint32_t x[4];
		x[0] = 12 + v->shared.l;
		x[1] = v->rid; x[2] = v->pos;
		memcpy(&x[3], &v->qual, 4);
		bgzf_write(fp->fp, x, 16);
		bgzf_write(fp->fp, v->shared.s, v->shared.l);
		x[0] = v->indiv.l;
		bgzf_write(fp->fp, x, 4);
		bgzf_write(fp->fp, v->indiv.s, v->indiv.l);
	} else {
		vcf_format1(h, v, &fp->line);
		fwrite(fp->line.s, 1, fp->line.l, fp->fp);
		fputc('\n', fp->fp);
	}
	return 0;
}
