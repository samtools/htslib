#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "kstring.h"
#include "bgzf.h"
#include "vcf.h"

typedef struct {
	uint32_t info[3]; // Number:20, var:4, Type:4, ColType:4
} vcf_keyinfo_t;

#include "khash.h"
KHASH_MAP_INIT_STR(vdict, vcf_keyinfo_t)

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define VCF_DT_CTG  0
#define VCF_DT_FLT  1
#define VCF_DT_INFO 2
#define VCF_DT_FMT  3

#define VCF_TP_BOOL 0
#define VCF_TP_INT  1
#define VCF_TP_REAL 2
#define VCF_TP_STR  3

#define VCF_VTP_FIXED 0
#define VCF_VTP_VAR   1
#define VCF_VTP_A     2
#define VCF_VTP_G     3

char *vcf_col_name[] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" };
char *vcf_type_name[] = { "Flag", "Integer", "Float", "String" };

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

int vcf_parse_hline2(const char *str, uint32_t *info, int *id_beg, int *id_end)
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

int main()
{
	int id_beg, id_end;
	uint32_t rst;
	vcf_parse_hline2("##INFO=<Number=A,Type=Integer,Description=\"gatca,agct=gact,\",ID=MQ>", &rst, &id_beg, &id_end);
	vcf_parse_hline2("##contig=<ID=20,length=62435964,assembly=B36>", &rst, &id_beg, &id_end);
	return 0;
}
