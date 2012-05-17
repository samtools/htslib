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
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

uint32_t bcf_missing_float = 0x7F800001;
uint8_t bcf_type_shift[] = { 0, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static bcf_idinfo_t bcf_idinfo_def = { { 15, 15, 15 }, -1 };

/*************************
 *** VCF header parser ***
 *************************/

// return: positive => contig; zero => INFO/FILTER/FORMAT; negative => error or skipped
int bcf_hdr_parse_line2(const char *str, uint32_t *info, int *id_beg, int *id_end)
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
	if (q - p == 4 && strncmp(p, "INFO", 4) == 0) ctype = BCF_HL_INFO;
	else if (q - p == 6 && strncmp(p, "FILTER", 6) == 0) ctype = BCF_HL_FLT;
	else if (q - p == 6 && strncmp(p, "FORMAT", 6) == 0) ctype = BCF_HL_FMT;
	else if (q - p == 6 && strncmp(p, "contig", 6) == 0) ctype = BCF_HL_CTG;
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
			if (q - val == 7 && strncmp(val, "Integer", 7) == 0) type = BCF_HT_INT;
			else if (q - val == 5 && strncmp(val, "Float", 5) == 0) type = BCF_HT_REAL;
			else if (q - val == 6 && strncmp(val, "String", 6) == 0) type = BCF_HT_STR;
			else if (q - val == 4 && strncmp(val, "Flag", 6) == 0) type = BCF_HT_FLAG;
		} else if (which == 3) {
			if (*val == 'A') var = BCF_VL_A;
			else if (*val == 'G') var = BCF_VL_G;
			else if (isdigit(*val)) var = BCF_VL_FIXED, num = strtol(val, &tmp, 10);
			else var = BCF_VL_VAR;
			if (var != BCF_VL_FIXED) num = 0xfffff;
		} else if (which == 4) {
			if (isdigit(*val)) ctg_len = strtol(val, &tmp, 10);
		}
		p = q + 1;
	}
	if (ctype == BCF_HL_CTG) {
		if (ctg_len > 0) return ctg_len;
		else return -5;
	} else {
		if (ctype == BCF_HL_FLT) num = 0;
		if (type == BCF_HT_FLAG) {
			if (num != 0 && hts_verbose >= 2)
				fprintf(stderr, "[W::%s] ignore Number for a Flag\n", __func__);
			num = 0, var = BCF_VL_FIXED; // if Flag VCF type, force to change num to 0
		}
		if (num == 0) type = BCF_HT_FLAG, var = BCF_VL_FIXED; // conversely, if num==0, force the type to Flag
		if (*id_beg < 0 || type < 0 || num < 0 || var < 0) return -5; // missing information
		*info = (uint32_t)num<<12 | var<<8 | type<<4 | ctype;
		//printf("%d, %s, %d, %d, [%d,%d]\n", ctype, bcf_type_name[type], var, num, *id_beg, *id_end);
		return 0;
	}
}

int bcf_hdr_parse1(bcf_hdr_t *h, const char *str)
{
	khint_t k;
	if (*str != '#') return -1;
	if (str[1] == '#') {
		uint32_t info;
		int len, ret, id_beg, id_end;
		char *s;

		len = bcf_hdr_parse_line2(str, &info, &id_beg, &id_end);
		if (len < 0) return -1;
		s = (char*)malloc(id_end - id_beg + 1);
		strncpy(s, str + id_beg, id_end - id_beg);
		s[id_end - id_beg] = 0;
		if (len > 0) { // a contig line
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
			k = kh_put(vdict, d, s, &ret);
			if (ret == 0) {
				if (hts_verbose >= 2)
					fprintf(stderr, "[W::%s] Duplicated contig name '%s'. Skipped.\n", __func__, s);
				free(s);
			} else {
				kh_val(d, k) = bcf_idinfo_def;
				kh_val(d, k).id = kh_size(d) - 1;
				kh_val(d, k).info[0] = len;
			}
		} else { // a FILTER/INFO/FORMAT line
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
			k = kh_put(vdict, d, s, &ret);
			if (ret) { // absent from the dict
				kh_val(d, k) = bcf_idinfo_def;
				kh_val(d, k).info[info&0xf] = info;
				kh_val(d, k).id = kh_size(d) - 1;
			} else {
				kh_val(d, k).info[info&0xf] = info;
				free(s);
			}
		}
	} else {
		int i = 0;
		const char *p, *q;
		vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
		// check if "PASS" is in the dictionary
		k = kh_get(vdict, d, "PASS");
		if (k == kh_end(d)) bcf_hdr_parse1(h, "##FILTER=<ID=PASS>"); // if not, add it; this is a recursion
		// add samples
		d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
		for (p = q = str;; ++q) {
			int ret;
			if (*q != '\t' && *q != 0) continue;
			if (++i > 9) {
				char *s;
				s = (char*)malloc(q - p + 1);
				strncpy(s, p, q - p);
				s[q - p] = 0;
				k = kh_put(vdict, d, s, &ret);
				if (ret) { // absent
					kh_val(d, k) = bcf_idinfo_def;
					kh_val(d, k).id = kh_size(d) - 1;
				} else {
					if (hts_verbose >= 2)
						fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
				}
			}
			if (*q == 0) break;
			p = q + 1;
		}
	}
	return 0;
}

int bcf_hdr_sync(bcf_hdr_t *h)
{
	int i;
	for (i = 0; i < 3; ++i) {
		khint_t k;
		vdict_t *d = (vdict_t*)h->dict[i];
		h->n[i]  = kh_size(d);
		h->id[i] = (bcf_idpair_t*)malloc(kh_size(d) * sizeof(bcf_idpair_t));
		for (k = kh_begin(d); k != kh_end(d); ++k) {
			if (!kh_exist(d, k)) continue;
			h->id[i][kh_val(d, k).id].key = kh_key(d, k);
			h->id[i][kh_val(d, k).id].val = &kh_val(d, k);
		}
	}
	return 0;
}

int bcf_hdr_parse(bcf_hdr_t *h)
{
	char *p, *q;
	for (p = q = h->text;; ++q) {
		int c;
		if (*q != '\n' && *q != 0) continue;
		c = *q; *q = 0;
		bcf_hdr_parse1(h, p);
		*q = c;
		if (*q == 0) break;
		p = q + 1;
	}
	bcf_hdr_sync(h);
	return 0;
}

/**********************
 *** BCF header I/O ***
 **********************/

bcf_hdr_t *bcf_hdr_init(void)
{
	int i;
	bcf_hdr_t *h;
	h = (bcf_hdr_t*)calloc(1, sizeof(bcf_hdr_t));
	for (i = 0; i < 3; ++i)
		h->dict[i] = kh_init(vdict);
	return h;
}

void bcf_hdr_destroy(bcf_hdr_t *h)
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

bcf_hdr_t *bcf_hdr_read(BGZF *fp)
{
	uint8_t magic[4];
	bcf_hdr_t *h;
	h = bcf_hdr_init();
	bgzf_read(fp, magic, 4);
	bgzf_read(fp, &h->l_text, 4);
	h->text = (char*)malloc(h->l_text);
	bgzf_read(fp, h->text, h->l_text);
	bcf_hdr_parse(h);
	return h;
}

void bcf_hdr_write(BGZF *fp, const bcf_hdr_t *h)
{
	bgzf_write(fp, "BCF\2", 4);
	bgzf_write(fp, &h->l_text, 4);
	bgzf_write(fp, h->text, h->l_text);
}

/********************
 *** BCF site I/O ***
 ********************/

bcf1_t *bcf_init1()
{
	bcf1_t *v;
	v = (bcf1_t*)calloc(1, sizeof(bcf1_t));
	return v;
}

void bcf_destroy1(bcf1_t *v)
{
	free(v->shared.s); free(v->indiv.s);
	free(v);
}

int bcf_read1(BGZF *fp, bcf1_t *v)
{
	uint32_t x[8];
	int ret;
	if ((ret = bgzf_read(fp, x, 32)) != 32) {
		if (ret == 0) return -1;
		return -2;
	}
	ks_resize(&v->shared, x[0]);
	ks_resize(&v->indiv, x[1]);
	memcpy(v, x + 2, 24);
	v->shared.l = x[0], v->indiv.l = x[1];
	bgzf_read(fp, v->shared.s, v->shared.l);
	bgzf_read(fp, v->indiv.s, v->indiv.l);
	return 0;
}

int bcf_write1(BGZF *fp, const bcf1_t *v)
{
	uint32_t x[8];
	x[0] = v->shared.l;
	x[1] = v->indiv.l;
	memcpy(x + 2, v, 24);
	bgzf_write(fp, x, 32);
	bgzf_write(fp, v->shared.s, v->shared.l);
	bgzf_write(fp, v->indiv.s, v->indiv.l);
	return 0;
}

/**********************
 *** VCF header I/O ***
 **********************/

bcf_hdr_t *vcf_hdr_read(htsFile *fp)
{
	if (!fp->is_bin) {
		kstring_t txt, *s = &fp->line;
		bcf_hdr_t *h;
		h = bcf_hdr_init();
		txt.l = txt.m = 0; txt.s = 0;
		while (hts_getline(fp, KS_SEP_LINE, s) >= 0) {
			if (s->l == 0) continue;
			if (s->s[0] != '#') {
				if (hts_verbose >= 2)
					fprintf(stderr, "[E::%s] no sample line\n", __func__);
				free(txt.s);
				bcf_hdr_destroy(h);
				return 0;
			}
			if (s->s[1] != '#' && fp->fn_aux) { // insert contigs here
				int dret;
				gzFile f;
				kstream_t *ks;
				kstring_t tmp;
				tmp.l = tmp.m = 0; tmp.s = 0;
				f = gzopen(fp->fn_aux, "r");
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
		bcf_hdr_parse(h);
		return h;
	} else return bcf_hdr_read((BGZF*)fp->fp);
}

void vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h)
{
	if (!fp->is_bin) {
		fwrite(h->text, 1, h->l_text, (FILE*)fp->fp);
		fputc('\n', (FILE*)fp->fp);
	} else bcf_hdr_write((BGZF*)fp->fp, h);
}

/***********************
 *** Typed value I/O ***
 ***********************/

void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
{
	int32_t max = INT32_MIN + 1, min = INT32_MAX;
	int i;
	if (n == 0) bcf_enc_size(s, 0, BCF_BT_NULL);
	else if (n == 1) bcf_enc_int1(s, a[0]);
	else {
		if (wsize <= 0) wsize = n;
		for (i = 0; i < n; ++i) {
			if (a[i] == INT32_MIN) continue;
			if (max < a[i]) max = a[i];
			if (min > a[i]) min = a[i];
		}
		if (max <= INT8_MAX && min > INT8_MIN) {
			bcf_enc_size(s, wsize, BCF_BT_INT8);
			for (i = 0; i < n; ++i)
				kputc(a[i] == INT32_MIN? INT8_MIN : a[i], s);
		} else if (max <= INT16_MAX && min > INT16_MIN) {
			bcf_enc_size(s, wsize, BCF_BT_INT16);
			for (i = 0; i < n; ++i) {
				int16_t x = a[i] == INT32_MIN? INT16_MIN : a[i];
				kputsn((char*)&x, 2, s);
			}
		} else {
			bcf_enc_size(s, wsize, BCF_BT_INT32);
			for (i = 0; i < n; ++i) {
				int32_t x = a[i] == INT32_MIN? INT32_MIN : a[i];
				kputsn((char*)&x, 4, s);
			}
		}
	}
}

void bcf_enc_vfloat(kstring_t *s, int n, float *a)
{
	bcf_enc_size(s, n, BCF_BT_FLOAT);
	kputsn((char*)a, n << 2, s);
}

void bcf_enc_vchar(kstring_t *s, int l, char *a)
{
	bcf_enc_size(s, l, BCF_BT_CHAR);
	kputsn(a, l, s);
}

void bcf_fmt_array(kstring_t *s, int n, int type, void *data)
{
	int j = 0;
	if (type == BCF_BT_INT8) {
		int8_t *p = (int8_t*)data;
		for (j = 0; j < n && *p != INT8_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == BCF_BT_CHAR) {
		char *p = (char*)data;
		for (j = 0; j < n && *p; ++j, ++p) kputc(*p, s);
	} else if (type == BCF_BT_INT32) {
		int32_t *p = (int32_t*)data;
		for (j = 0; j < n && *p != INT32_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	} else if (type == BCF_BT_FLOAT) {
		float *p = (float*)data;
		for (j = 0; j < n && *(int32_t*)p != 0x7F800001; ++j, ++p) {
			if (j) kputc(',', s);
			ksprintf(s, "%g", *p);
		}
	} else if (type == BCF_BT_INT16) {
		int16_t *p = (int16_t*)data;
		for (j = 0; j < n && *p != INT16_MIN; ++j, ++p) {
			if (j) kputc(',', s);
			kputw(*p, s);
		}
	}
	if (n && j == 0) kputc('.', s);
}

uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr)
{
	int x, type;
	x = bcf_dec_size(ptr, &ptr, &type);
	bcf_fmt_array(s, x, type, ptr);
	return ptr + (x << bcf_type_shift[type]);
}

/********************
 *** VCF site I/O ***
 ********************/

typedef struct {
	int key, max_m, size, offset;
	uint32_t is_gt:1, max_g:15, max_l:16;
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

int vcf_parse1(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
{
	int i = 0;
	char *p, *q, *r, *t;
	fmt_aux_t *fmt = 0;
	kstring_t *str, *mem = (kstring_t*)&h->mem;
	khint_t k;
	ks_tokaux_t aux;

	mem->l = v->shared.l = v->indiv.l = 0;
	str = &v->shared;
	v->n_fmt = 0;
	memset(&aux, 0, sizeof(ks_tokaux_t));
	for (p = kstrtok(s->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
		q = (char*)aux.p;
		*q = 0;
		if (i == 0) { // CHROM
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
			k = kh_get(vdict, d, p);
			if (k == kh_end(d)) {
				if (hts_verbose >= 2)
					fprintf(stderr, "[W::%s] can't find '%s' in the sequence dictionary\n", __func__, p);
				return 0;
			} else v->rid = kh_val(d, k).id;
		} else if (i == 1) { // POS
			v->pos = atoi(p) - 1;
		} else if (i == 2) { // ID
			if (strcmp(p, ".")) bcf_enc_vchar(str, q - p, p);
			else bcf_enc_size(str, 0, BCF_BT_CHAR);
		} else if (i == 3) { // REF
			bcf_enc_vchar(str, q - p, p);
			v->n_allele = 1, v->rlen = q - p;
		} else if (i == 4) { // ALT
			if (strcmp(p, ".")) {
				for (r = t = p;; ++r) {
					if (*r == ',' || *r == 0) {
						bcf_enc_vchar(str, r - t, t);
						t = r + 1;
						++v->n_allele;
					}
					if (r == q) break;
				}
			}
		} else if (i == 5) { // QUAL
			if (strcmp(p, ".")) v->qual = atof(p);
			else memcpy(&v->qual, &bcf_missing_float, 4);
		} else if (i == 6) { // FILTER
			if (strcmp(p, ".")) {
				int32_t *a;
				int n_flt = 1, i;
				ks_tokaux_t aux1;
				vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
				// count the number of filters
				if (*(q-1) == ';') *(q-1) = 0;
				for (r = p; *r; ++r)
					if (*r == ';') ++n_flt;
				a = (int32_t*)alloca(n_flt * 4);
				// add filters
				for (t = kstrtok(p, ";", &aux1), i = 0; t; t = kstrtok(0, 0, &aux1)) {
					*(char*)aux1.p = 0;
					k = kh_get(vdict, d, t);
					if (k == kh_end(d)) { // not defined
						if (hts_verbose >= 2) fprintf(stderr, "[W::%s] undefined FILTER '%s'\n", __func__, t);
					} else a[i++] = kh_val(d, k).id;
				}
				n_flt = i;
				bcf_enc_vint(str, n_flt, a, -1);
			} else bcf_enc_vint(str, 0, 0, -1);
		} else if (i == 7) { // INFO
			char *key;
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
			v->n_info = 0;
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
					if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_INFO] == 15) { // not defined in the header
						if (hts_verbose >= 2) fprintf(stderr, "[W::%s] undefined INFO '%s'\n", __func__, key);
					} else { // defined in the header
						uint32_t y = kh_val(d, k).info[BCF_HL_INFO];
						++v->n_info;
						bcf_enc_int1(str, kh_val(d, k).id);
						if (val == 0) {
							bcf_enc_size(str, 0, BCF_BT_NULL);
						} else if ((y>>4&0xf) == BCF_HT_FLAG || (y>>4&0xf) == BCF_HT_STR) { // if Flag has a value, treat it as a string
							bcf_enc_vchar(str, end - val, val);
						} else { // int/float value/array
							int i, n_val;
							char *t;
							for (t = val, n_val = 1; *t; ++t) // count the number of values
								if (*t == ',') ++n_val;
							if ((y>>4&0xf) == BCF_HT_INT) {
								int32_t *z;
								z = (int32_t*)alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtol(t, &t, 10);
								bcf_enc_vint(str, n_val, z, -1);
								if (strcmp(key, "END") == 0) v->rlen = z[0] - v->pos;
							} else if ((y>>4&0xf) == BCF_HT_REAL) {
								float *z;
								z = (float*)alloca(n_val<<2);
								for (i = 0, t = val; i < n_val; ++i, ++t)
									z[i] = strtod(t, &t);
								bcf_enc_vfloat(str, n_val, z);
							}
						}
					}
					if (c == 0) break;
					r = end;
					key = r + 1;
				}
			}
		} else if (i == 8) { // FORMAT
			int j, l, m, g;
			ks_tokaux_t aux1;
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
			char *end = s->s + s->l;
			// count the number of format fields
			for (r = p, v->n_fmt = 1; *r; ++r)
				if (*r == ':') ++v->n_fmt;
			fmt = (fmt_aux_t*)alloca(v->n_fmt * sizeof(fmt_aux_t));
			// get format information from the dictionary
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				*(char*)aux1.p = 0;
				k = kh_get(vdict, d, t);
				if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_FMT] == 15) {
					if (hts_verbose >= 2)
						fprintf(stderr, "[W::%s] FORMAT '%s' is not defined in the header\n", __func__, t);
					v->n_fmt = 0;
					break;
				} else {
					fmt[j].max_l = fmt[j].max_m = fmt[j].max_g = 0;
					fmt[j].key = kh_val(d, k).id;
					fmt[j].is_gt = !strcmp(t, "GT");
					fmt[j].y = h->id[0][fmt[j].key].val->info[BCF_HL_FMT];
				}
			}
			// compute max
			for (r = q + 1, j = 0, m = l = g = 1, v->n_sample = 0;; ++r, ++l) {
				if (*r == '\t') *r = 0;
				if (*r == ':' || *r == '\0') { // end of a sample
					if (fmt[j].max_m < m) fmt[j].max_m = m;
					if (fmt[j].max_l < l - 1) fmt[j].max_l = l - 1;
					if (fmt[j].is_gt && fmt[j].max_g < g) fmt[j].max_g = g;
					l = 0, m = g = 1;
					if (*r) ++j;
					else j = 0, ++v->n_sample;
				} else if (*r == ',') ++m;
				else if (*r == '|' || *r == '/') ++g;
				if (r == end) break;
			}
			// allocate memory for arrays
			for (j = 0; j < v->n_fmt; ++j) {
				fmt_aux_t *f = &fmt[j];
				if ((f->y>>4&0xf) == BCF_HT_STR) {
					f->size = f->is_gt? f->max_g << 2 : f->max_l;
				} else if ((f->y>>4&0xf) == BCF_HT_REAL || (f->y>>4&0xf) == BCF_HT_INT) {
					f->size = f->max_m << 2;
				} else abort(); // I do not know how to do with Flag in the genotype fields
				align_mem(mem);
				f->offset = mem->l;
				ks_resize(mem, mem->l + v->n_sample * f->size);
				mem->l += v->n_sample * f->size;
			}
			for (j = 0; j < v->n_fmt; ++j)
				fmt[j].buf = (uint8_t*)mem->s + fmt[j].offset;
			// fill the sample fields; at beginning of the loop, t points to the first char of a format
			for (t = q + 1, j = m = 0;; ++t) { // j: fmt id, m: sample id
				fmt_aux_t *z = &fmt[j];
				if ((z->y>>4&0xf) == BCF_HT_STR) {
					if (z->is_gt) { // genotypes
						int32_t is_phased = 0, *x = (int32_t*)(z->buf + z->size * m);
						for (l = 0;; ++t) {
							if (*t == '.') ++t, x[l++] = is_phased;
							else x[l++] = (strtol(t, &t, 10) + 1) << 1 | is_phased;
							is_phased = (*t == '|');
							if (*t == ':' || *t == 0) break;
						}
						for (; l != z->size>>2; ++l) x[l] = INT32_MIN;
					} else {
						char *x = (char*)z->buf + z->size * m;
						for (r = t, l = 0; *t != ':' && *t; ++t) x[l++] = *t;
						for (; l != z->size; ++l) x[l] = 0;
					}
				} else if ((z->y>>4&0xf) == BCF_HT_INT) {
					int32_t *x = (int32_t*)(z->buf + z->size * m);
					for (l = 0;; ++t) {
						if (*t == '.') x[l++] = INT32_MIN, ++t; // ++t to skip "."
						else x[l++] = strtol(t, &t, 10);
						if (*t == ':' || *t == 0) break;
					}
					for (; l != z->size>>2; ++l) x[l] = INT32_MIN;
				} else if ((z->y>>4&0xf) == BCF_HT_REAL) {
					float *x = (float*)(z->buf + z->size * m);
					for (l = 0;; ++t) {
						if (*t == '.' && !isdigit(t[1])) *(int32_t*)&x[l++] = 0x7F800001, ++t; // ++t to skip "."
						else x[l++] = strtod(t, &t);
						if (*t == ':' || *t == 0) break;
					}
					for (; l != z->size>>2; ++l) *(int32_t*)(x+l) = 0x7F800001;
				} else abort();
				if (*t == 0) {
					if (t == end) break;
					++m, j = 0;
				} else if (*t == ':') ++j;
			}
			break;
		}
	}
	// write individual genotype information
	str = &v->indiv;
	if (v->n_sample > 0) {
		for (i = 0; i < v->n_fmt; ++i) {
			fmt_aux_t *z = &fmt[i];
			bcf_enc_int1(str, z->key);
			if ((z->y>>4&0xf) == BCF_HT_STR && !z->is_gt) {
				bcf_enc_size(str, z->size, BCF_BT_CHAR);
				kputsn((char*)z->buf, z->size * v->n_sample, str);
			} else if ((z->y>>4&0xf) == BCF_HT_INT || z->is_gt) {
				bcf_enc_vint(str, (z->size>>2) * v->n_sample, (int32_t*)z->buf, z->size>>2);
			} else {
				bcf_enc_size(str, z->size>>2, BCF_BT_FLOAT);
				kputsn((char*)z->buf, z->size * v->n_sample, str);
			}
		}
	}
	return 0;
}

int vcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
	if (!fp->is_bin) {
		int ret;
		ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
		if (ret < 0) return -1;
		ret = vcf_parse1(&fp->line, h, v);
		return 0;
	} else return bcf_read1((BGZF*)fp->fp, v);
}

uint8_t *bcf_unpack_fmt_core(uint8_t *ptr, int n_sample, int n_fmt, bcf_fmt_t *fmt)
{
	int i;
	for (i = 0; i < n_fmt; ++i) {
		bcf_fmt_t *f = &fmt[i];
		f->id = bcf_dec_typed_int1(ptr, &ptr);
		f->n = bcf_dec_size(ptr, &ptr, &f->type);
		f->size = f->n << bcf_type_shift[f->type];
		f->p = ptr;
		ptr += n_sample * f->size;
	}
	return ptr;
}

int vcf_format1(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)
{
	uint8_t *ptr = (uint8_t*)v->shared.s;
	int i;
	s->l = 0;
	kputs(h->id[BCF_DT_CTG][v->rid].key, s); kputc('\t', s); // CHROM
	kputw(v->pos + 1, s); kputc('\t', s); // POS
	// ID
	ptr = bcf_fmt_sized_array(s, ptr);
	kputc('\t', s);
	if (v->n_allele) { // REF and ALT
		for (i = 0; i < v->n_allele; ++i) {
			if (i) kputc(i == 1? '\t' : ',', s);
			ptr = bcf_fmt_sized_array(s, ptr);
		}
		if (v->n_allele == 1) kputsn("\t.\t", 3, s);
		else kputc('\t', s);
	} else kputsn(".\t.\t", 4, s);
	if (memcmp(&v->qual, &bcf_missing_float, 4) == 0) kputsn(".\t", 2, s); // QUAL
	else ksprintf(s, "%g\t", v->qual);
	if (*ptr>>4) { // FILTER
		int32_t x, y;
		int type, i;
		x = bcf_dec_size(ptr, &ptr, &type);
		for (i = 0; i < x; ++i) {
			if (i) kputc(';', s);
			y = bcf_dec_int1(ptr, type, &ptr);
			kputs(h->id[BCF_DT_ID][y].key, s);
		}
		kputc('\t', s);
	} else {
		kputsn(".\t", 2, s);
		++ptr;
	}
	if (v->n_info) {
		for (i = 0; i < (int)v->n_info; ++i) {
			int32_t x;
			if (i) kputc(';', s);
			x = bcf_dec_typed_int1(ptr, &ptr);
			kputs(h->id[BCF_DT_ID][x].key, s);
			if (*ptr>>4) { // more than zero element
				kputc('=', s);
				ptr = bcf_fmt_sized_array(s, ptr);
			} else ++ptr; // skip 0
		}
	} else kputc('.', s);
	// FORMAT and individual information
	ptr = (uint8_t*)v->indiv.s;
	if (v->n_sample && v->n_fmt) { // FORMAT
		int i, j, l, gt_i = -1;
		bcf_fmt_t *fmt;
		fmt = (bcf_fmt_t*)alloca(v->n_fmt * sizeof(bcf_fmt_t));
		ptr = bcf_unpack_fmt_core(ptr, v->n_sample, v->n_fmt, fmt);
		for (i = 0; i < (int)v->n_fmt; ++i) {
			kputc(i? ':' : '\t', s);
			kputs(h->id[BCF_DT_ID][fmt[i].id].key, s);
			if (strcmp(h->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0) gt_i = i;
		}
		for (j = 0; j < v->n_sample; ++j) {
			kputc('\t', s);
			for (i = 0; i < (int)v->n_fmt; ++i) {
				bcf_fmt_t *f = &fmt[i];
				if (i) kputc(':', s);
				if (gt_i == i) {
					int8_t *x = (int8_t*)(f->p + j * f->size); // FIXME: does not work with n_alt >= 64
					for (l = 0; l < f->n && x[l] != INT8_MIN; ++l) {
						if (l) kputc("/|"[x[l]&1], s);
						if (x[l]>>1) kputw((x[l]>>1) - 1, s);
						else kputc('.', s);
					}
					if (l == 0) kputc('.', s);
				} else bcf_fmt_array(s, f->n, f->type, f->p + j * f->size);
			}
		}
	}
	return 0;
}

int vcf_write1(htsFile *fp, const bcf_hdr_t *h, const bcf1_t *v)
{
	if (!fp->is_bin) {
		vcf_format1(h, v, &fp->line);
		fwrite(fp->line.s, 1, fp->line.l, (FILE*)fp->fp);
		fputc('\n', (FILE*)fp->fp);
	} else return bcf_write1((BGZF*)fp->fp, v);
	return 0;
}

/************************
 * Data access routines *
 ************************/

int vcf_id2int(const bcf_hdr_t *h, int which, const char *id)
{
	khint_t k;
	vdict_t *d = (vdict_t*)h->dict[which];
	k = kh_get(vdict, d, id);
	return k == kh_end(d)? -1 : kh_val(d, k).id;
}

bcf_fmt_t *vcf_unpack_fmt(const bcf_hdr_t *h, const bcf1_t *v)
{
	bcf_fmt_t *fmt;
	if (v->n_fmt == 0) return 0;
	fmt = (bcf_fmt_t*)malloc(v->n_fmt * sizeof(bcf_fmt_t));
	bcf_unpack_fmt_core((uint8_t*)v->indiv.s, v->n_sample, v->n_fmt, fmt);
	return fmt;
}
