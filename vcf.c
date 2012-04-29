#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "kstring.h"
#include "bgzf.h"
#include "vcf.h"

#include "khash.h"
KHASH_MAP_INIT_STR(vdict, vcf_keyinfo_t)
typedef khash_t(vdict) vdict_t;

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define VCF_DT_FLT  0
#define VCF_DT_INFO 1
#define VCF_DT_FMT  2
#define VCF_DT_CTG  3

#define VCF_TP_BOOL 0
#define VCF_TP_INT  1
#define VCF_TP_REAL 2
#define VCF_TP_STR  3

#define VCF_VTP_FIXED 0
#define VCF_VTP_VAR   1
#define VCF_VTP_A     2
#define VCF_VTP_G     3

int vcf_verbose = 3; // 1: error; 2: warning; 3: message; 4: progress; 5: debugging; >=10: pure debugging
char *vcf_col_name[] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" };
char *vcf_type_name[] = { "Flag", "Integer", "Float", "String" };
static vcf_keyinfo_t vcf_keyinfo_def = { { 15, 15, 15 }, -1, -1, -1, -1 };

/******************
 * Basic routines *
 ******************/

static inline uint8_t *expand_mem(uint8_t *mem, int32_t *max, int32_t new_len)
{
	if (new_len > *max) {
		*max = new_len;
		kroundup32(*max);
		mem = (uint8_t*)realloc(mem, *max);
	}
	return mem;
}

static inline uint8_t *append_text(uint8_t *mem, int32_t *len, int32_t *max, const char *text, int text_len)
{
	mem = expand_mem(mem, max, (*len) + text_len + 1); // +1 to include '\0'
	memcpy(mem + (*len), text, text_len);
	*len += text_len;
	mem[*len] = 0;
	return mem;
}

/********************
 * Parse VCF header *
 ********************/

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
		if (num == 0) type = VCF_TP_BOOL, var = VCF_VTP_FIXED;
		if (*id_beg < 0 || type < 0 || num < 0 || var < 0) return -5; // missing information
		*info = (uint32_t)num<<20 | var<<8 | type<<4 | ctype;
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

int main()
{
	vcf_hdr_t *h;
	h = vcf_hdr_init();
	vcf_hdr_parse1(h, "##INFO=<Number=A,Type=Integer,Description=\"gatca,agct=gact,\",ID=MQ>");
	vcf_hdr_parse1(h, "##contig=<ID=NA00001,length=62435964,assembly=B36>");
	vcf_hdr_parse1(h, "##contig=<ID=20,length=62435964,assembly=B36>");
	vcf_hdr_parse1(h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003");
	vcf_hdr_sync(h);
	vcf_hdr_destroy(h);

	vcfFile *fp;
	fp = vcf_open("t.vcf", "r", 0);
	h = vcf_hdr_read(fp);
	vcf_hdr_destroy(h);
	return 0;
}
