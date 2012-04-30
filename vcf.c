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
char *vcf_col_name[] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" };
char *vcf_type_name[] = { "Flag", "Integer", "Float", "String" };
uint8_t vcf_type_size[] = { 0, 1, 2, 4, 8, 4, 8, 1, 1, 0, 1, 2, 4, 1, 0, 0 };
static vcf_keyinfo_t vcf_keyinfo_def = { { 15, 15, 15 }, -1, -1, -1, -1 };

/******************
 * Basic routines *
 ******************/

static inline void align_mem(kstring_t *s)
{
	if (s->l&7) {
		uint64_t zero = 0;
		int l = (s->l + 7)>>3<<3;
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
		fp->buf = calloc(1, sizeof(kstring_t));
		if (fp->is_write) {
			fp->fp = strcmp(fn, "-")? fopen(fn, "rb") : fdopen(fileno(stdout), "rb");
		} else {
			gzFile gzfp;
			gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
			if (gzfp) fp->fp = ks_init(gzfp);
		}
	}
	if (fp->fp == 0) {
		if (vcf_verbose >= 2)
			fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);
		free(fp);
		return 0;
	}
	return fp;
}

void vcf_close(vcfFile *fp)
{
	if (!fp->is_bin) {
		kstring_t *s = (kstring_t*)fp->buf;
		free(s->s); free(s);
		if (!fp->is_write) {
			gzFile gzfp = ((kstream_t*)fp->fp)->f;
			ks_destroy(fp->fp);
			gzclose(gzfp);
		} else fclose(fp->fp);
	} else bgzf_close(fp->fp);
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
	free(h->key); free(h->r2k); free(h->s2k);
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
	if (fp->is_bin) {
	} else {
		int dret;
		kstring_t txt, *s = (kstring_t*)fp->buf;
		txt.l = txt.m = 0; txt.s = 0;
		while (ks_getuntil(fp->fp, KS_SEP_LINE, s, &dret) >= 0) {
			if (s->l == 0) continue;
			if (s->s[0] != '#') {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[E::%s] no sample line\n", __func__);
				free(txt.s);
				return 0;
			}
			kputsn(s->s, s->l, &txt);
			if (s->s[1] != '#') break;
			else kputc('\n', &txt);
		}
		h = vcf_hdr_init();
		h->l_text = txt.l + 1; // including NULL
		h->text = txt.s;
	}
	if (vcf_hdr_parse(h) != 0) { // FIXME: vcf_hdr_parse() always returns 0
		vcf_hdr_destroy(h);
		return 0;
	} else return h;
}

/*******************
 * Typed value I/O *
 *******************/

void vcf_enc_int(kstring_t *s, int n, long *a)
{
	long max = LONG_MIN + 1, min = LONG_MAX;
	int i;
	if (n == 0) vcf_enc_size(s, 0, VCF_RT_INT8);
	else if (n == 1) vcf_enc_int1(s, a[0]);
	else {
		for (i = 0; i < n; ++i) {
			if (a[i] == LONG_MIN) continue;
			if (max < a[i]) max = a[i];
			if (min > a[i]) min = a[i];
		}
		if (max <= INT8_MAX && min > INT8_MIN) {
			vcf_enc_size(s, n, VCF_RT_INT8);
			for (i = 0; i < n; ++i)
				kputc(a[i] == LONG_MIN? INT8_MIN : a[i], s);
		} else if (max <= INT16_MAX && min > INT16_MIN) {
			vcf_enc_size(s, n, VCF_RT_INT16);
			for (i = 0; i < n; ++i) {
				int16_t x = a[i] == LONG_MIN? INT16_MIN : a[i];
				kputsn((char*)&x, 2, s);
			}
		} else {
			vcf_enc_size(s, n, VCF_RT_INT32);
			for (i = 0; i < n; ++i) {
				int32_t x = a[i] == LONG_MIN? INT32_MIN : a[i];
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

/******************
 * VCF record I/O *
 ******************/

vcf1_t *vcf_init1()
{
	vcf1_t *v;
	v = calloc(1, sizeof(vcf1_t));
	return v;
}

int vcf_parse1(kstring_t *s, const vcf_hdr_t *h, vcf1_t *v)
{
	int ret, i = 0, ori_l;
	char *p, *q, *r;
	vdict_t *d = (vdict_t*)h->dict;
	kstring_t str;
	khint_t k;

	str.l = v->l_str = 0; str.m = v->m_str; str.s = v->str;
	kputc(0, s); align_mem(s); ori_l = s->l; // the end of s will be used as a temporary buffer
	if (ret < 0) return -1;
	for (p = q = s->s;; ++q) {
		int c;
		if (*q != 0 && *q != '\t') continue;
		c = *q; *q = 0;
		if (i == 0) { // CHROM
			k = kh_get(vdict, d, p);
			if (k == kh_end(d) || kh_val(d, k).rid < 0) {
				if (vcf_verbose >= 2)
					fprintf(stderr, "[W::%s] can't find '%s' in the sequence dictionary\n", __func__, p);
				return 0;
			} else v->rid = kh_val(d, k).rid;
		} else if (i == 1) { // POS
			v->pos = atoi(p);
		} else if (i == 2) { // ID
			if (strcmp(p, ".")) kputsn(p, q - p, &str);
			kputc(0, &str);
		} else if (i == 3) { // REF
			v->o_ref = str.l;
			kputsn(p, q - p + 1, &str); // +1 to include NULL
		} else if (i == 4) { // ALT
			v->o_alt = str.l;
			if (strcmp(p, ".")) {
				v->n_alt = 1;
				for (r = p; r != q; ++r)
					if (*r == ',') *r = 0, ++v->n_alt;
			} else v->n_alt = 0;
			kputsn((char*)&v->n_alt, 2, &str);
			if (v->n_alt) kputsn(p, q - p + 1, &str);
		} else if (i == 5) { // QUAL
			if (strcmp(p, ".")) {
				v->qual = atof(p);
				kputsn((char*)&v->qual, 4, &str);
			} else *(uint32_t*)&v->qual = 0x7F800001;
		} else if (i == 6) { // FILTER
			v->o_flt = str.l;
			if (strcmp(p, ".")) {
				char *t;
				long *a;
				int n_flt = 1, i;
				for (r = p; *r; ++r)
					if (*r == ';') ++n_flt;
				ks_resize(s, ori_l + n_flt * sizeof(long)); // reuse the end of s as a buffer
				a = (long*)(s->s + ori_l);
				for (r = t = p, i = 0;; ++r) {
					int c;
					if (*r != ';' && *r != 0) continue;
					c = *r; *r = 0;
					k = kh_get(vdict, d, t);
					if (k == kh_end(d)) { // not defined
						if (vcf_verbose >= 2) fprintf(stderr, "[W::%s] undefined FILTER '%s'\n", __func__, t);
					} else a[i++] = kh_val(d, k).kid;
					*r = c;
					if (*r == 0 || r[1] == 0) break; // r[1] if the last char is ';'
					t = r + 1;
				}
				vcf_enc_int(&str, n_flt, a);
			} else vcf_enc_int(&str, 0, 0);
		} else if (i == 7) { // INFO
			char *key;
			int n_info = 0;
			v->o_info = str.l;
			kputsn("\0\0", 2, &str); // place holder for n_info
			if (strcmp(p, ".")) {
				// fill info
				printf("*** %s\n", p);
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
					printf("%s, %s\n", key, val? val : ".");
					k = kh_get(vdict, d, key);
					if (k == kh_end(d) || kh_val(d, k).info[VCF_DT_INFO] == 15) { // not defined in the header
						if (vcf_verbose >= 2) fprintf(stderr, "[W::%s] undefined INFO '%s'\n", __func__, key);
					} else { // defined in the header
						uint32_t y = kh_val(d, k).info[VCF_DT_INFO];
						if ((y>>4&0xf) == VCF_TP_FLAG) { // a flag defined in the dict
							if (val != 0 && vcf_verbose >= 2)
								fprintf(stderr, "[W::%s] INFO '%s' is defined as a flag in the header but has a value '%s' in VCF; value skipped\n", __func__, key, val);
							else {
								vcf_enc_int1(&str, kh_val(d, k).kid);
								++n_info;
							}
						} else if (val == 0) { // if not Flag, there must be value(s)
							if (vcf_verbose >= 2)
								fprintf(stderr, "[W::%s] INFO '%s' takes at least a value, but no value is found\n", __func__, key);
						} else if ((y>>4&0xf) == VCF_TP_STR) { // a string
							vcf_enc_int1(&str, kh_val(d, k).kid);
							vcf_enc_size(&str, 1, VCF_RT_CSTR);
							kputsn(val, end - val + 1, &str); // +1 to include NULL
							++n_info;
						} else { // an integer or float value or array
							int n_rec_val = 1, n_hdr_val;
							char *t;
							for (t = val; *t; ++t)
								if (*t == ',') ++n_rec_val;
							n_hdr_val = vcf_hdr_n_val(y, v->n_alt);
							if (n_hdr_val == -1) n_hdr_val = n_rec_val; // -1 if variable size
							else if (n_hdr_val != n_rec_val) { // check if number of values agrees with the header
								if (vcf_verbose >= 2)
									fprintf(stderr, "[W::%s] INFO '%s' takes %d value(s) but the header requires %d value(s); skipped\n", __func__, key, n_rec_val, n_hdr_val);
							} else {
								int i;
								++n_info;
								vcf_enc_int1(&str, kh_val(d, k).kid);
								ks_resize(s, ori_l + n_rec_val * 8);
								if ((y>>4&0xf) == VCF_TP_INT) {
									long *z = (long*)(s->s + ori_l);
									for (i = 0, t = val; i < n_rec_val; ++i, ++t)
										z[i] = strtol(t, &t, 10);
									vcf_enc_int(&str, n_rec_val, z);
								} else if ((y>>4&0xf) == VCF_TP_REAL) {
									float *z = (float*)(s->s + ori_l);
									for (i = 0, t = val; i < n_rec_val; ++i, ++t)
										z[i] = strtod(t, &t);
									vcf_enc_float(&str, n_rec_val, z);
								}
							}
						}
					}
					if (c == 0) break;
					r = end;
					key = r + 1;
				}
			}
			*(uint16_t*)(str.s + v->o_info) = n_info;
		}
		++i;
		*q = c;
		if (*q == 0) break;
		p = q + 1;
	}
	v->l_str = str.l; v->m_str = str.m; v->str = str.s;
	{
		str.l = str.m = 0; str.s = 0;
		vcf_format1(h, v, &str);
		puts(str.s);
		free(str.s);
	}
	return 0;
}

int vcf_read1(vcfFile *fp, const vcf_hdr_t *h, vcf1_t *v)
{
	if (fp->is_bin) {
	} else {
		int ret, dret;
		ret = ks_getuntil(fp->fp, KS_SEP_LINE, fp->buf, &dret);
		if (ret < 0) return -1;
		ret = vcf_parse1(fp->buf, h, v);
	}
	return 0;
}

int vcf_format1(const vcf_hdr_t *h, const vcf1_t *v, kstring_t *s)
{
	const char *p, *q;
	int i;

	s->l = 0;
	kputs(h->key[h->r2k[v->rid]].key, s); kputc('\t', s); // CHROM
	kputw(v->pos + 1, s); kputc('\t', s); // POS
	if (v->str[0]) { // ID
		kputsn(v->str, v->o_ref - 1, s);
		kputc('\t', s);
	} else kputsn(".\t", 2, s);
	if (v->str[v->o_ref]) { // REF
		kputsn(v->str + v->o_ref, v->o_alt - v->o_ref - 1, s);
		kputc('\t', s);
	} else kputsn(".\t", 2, s);
	if (v->n_alt) { // ALT
		int i;
		for (p = q = v->str + v->o_alt + 2, i = 0;; ++q)
			if (*q == 0) {
				if (i) kputc(',', s);
				kputsn(p, q - p, s);
				p = q + 1;
				if (++i == v->n_alt) break;
			}
		kputc('\t', s);
	} else kputsn(".\t", 2, s);
	if (*(uint32_t*)&v->qual == 0x7F800001) kputsn(".\t", 2, s); // QUAL
	else ksprintf(s, "%g\t", v->qual);
	{ // FILTER
		long x;
		x = vcf_dec_size((uint8_t*)v->str + v->o_flt);
		if (x>>4) {
			int size, o;
			size = vcf_type_size[v->str[v->o_flt]&0xf];
			for (i = 0, o = v->o_flt + (x&0xf); i < x>>4; ++i, o += size) {
				long y;
				if (i) kputc(';', s);
				if (size == 1) y = *(int8_t*)(v->str + o);
				else if (size == 2) y = *(int16_t*)(v->str + o);
				else y = *(int32_t*)(v->str + o);
				assert(y < h->n_key);
				kputs(h->key[y].key, s);
			}
			kputc('\t', s);
		} else kputsn(".\t", 2, s);
	}
	{ // INFO
		int n_info;
		uint8_t *q, *p = (uint8_t*)v->str + v->o_info + 2;
		n_info = *(uint16_t*)(v->str + v->o_info);
		for (i = 0; i < n_info; ++i) {
			long x;
			uint32_t y;
			if (i) kputc(';', s);
			x = vcf_dec_int1(p);
			p += 1 + vcf_type_size[*p&0xf];
			kputs(h->key[x].key, s);
			y = h->key[x].info->info[VCF_DT_INFO];
			if ((y>>4&0xf) == VCF_TP_FLAG) continue;
			else if ((y>>4&0xf) == VCF_TP_STR) {
				kputc('=', s);
				for (q = p; *q; ++q);
				kputsn((char*)p, q - p, s);
				p = q + 1;
			} else {
				int j, size;
				kputc('=', s);
				x = vcf_dec_size(p);
				size = vcf_type_size[*p&0xf];
				p += x&0xf;
				if ((y>>4&0xf) == VCF_TP_REAL) {
					for (j = 0; j < x>>4; ++j, p += 4) {
						if (j) kputc(',', s);
						ksprintf(s, "%g", *(float*)p);
					}
				} else if ((y>>4&0xf) == VCF_TP_INT) {
					if (size == 1) {
						for (j = 0; j < x>>4; ++j, p += size) {
							if (j) kputc(',', s);
							ksprintf(s, "%d", *(int8_t*)p);
						}
					} else if (size == 2) {
						for (j = 0; j < x>>4; ++j, p += size) {
							if (j) kputc(',', s);
							ksprintf(s, "%d", *(int16_t*)p);
						}
					} else {
						for (j = 0; j < x>>4; ++j, p += size) {
							if (j) kputc(',', s);
							ksprintf(s, "%d", *(int32_t*)p);
						}
					}
				}
			}
		}
	}
	return 0;
}

int main()
{
	vcf_hdr_t *h;
	/*
	h = vcf_hdr_init();
	vcf_hdr_parse1(h, "##INFO=<Number=A,Type=Integer,Description=\"gatca,agct=gact,\",ID=MQ>");
	vcf_hdr_parse1(h, "##contig=<ID=NA00001,length=62435964,assembly=B36>");
	vcf_hdr_parse1(h, "##contig=<ID=20,length=62435964,assembly=B36>");
	vcf_hdr_parse1(h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003");
	vcf_hdr_sync(h);
	vcf_hdr_destroy(h);
	*/

	vcfFile *fp;
	vcf1_t *v;
	fp = vcf_open("ex2.vcf", "r", 0);
	h = vcf_hdr_read(fp);
	v = vcf_init1();
	while (vcf_read1(fp, h, v) >= 0) {
	}
	vcf_hdr_destroy(h);
	return 0;
}
