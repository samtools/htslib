#ifndef VCF_H
#define VCF_H

#include <stdint.h>

/*******************
 * VCF file struct *
 *******************/

typedef struct {
	uint32_t is_bin:1, is_write:1, dummy:30;
	void *buf; // a string buffer, only for VCF
	void *fp;  // file pointer; actual type depending on is_bin and is_write
} vcfFile;

/*****************
 * Header struct *
 *****************/

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

typedef struct {
	uint32_t info[3]; // Number:20, var:4, Type:4, ColType:4
	int kid, rid, sid, rlen; // key-ID, ref-ID, sample-ID, ref-length
} vcf_keyinfo_t;

typedef struct {
	const char *key;
	const vcf_keyinfo_t *info;
} vcf_keypair_t;

typedef struct {
	int32_t n_ref, n_sample, n_key, l_text;
	vcf_keypair_t *key;
	int *r2k, *s2k;
	char *text;
	void *dict;
} vcf_hdr_t;

/**************
 * VCF record *
 **************/

#define VCF_RT_INT8		1
#define VCF_RT_INT16	2
#define VCF_RT_INT32	3
#define VCF_RT_INT64	4
#define VCF_RT_FLOAT	5
#define VCF_RT_BOOL		8
#define VCF_RT_CSTR		9
#define VCF_RT_UINT8	13

typedef struct {
	uint32_t key:28, type:3, is_vec:1;
	union {
		int32_t x;
		float f;
	} x;
} vcf_info_t;

typedef struct {
	int32_t rid; // CHROM
	int32_t pos; // POS
	int32_t end; // end position
	float qual;  // QUAL
	uint16_t n_alt, n_flt, n_info, n_fmt;
	int l_str, m_str;
	int o_ref, o_alt, o_flt, o_info; // offsets in str
	int *alt, *flt;
	char *str;
} vcf1_t;

/*******
 * API *
 *******/

#ifdef __cplusplus
extern "C" {
#endif

	vcfFile *vcf_open(const char *fn, const char *mode, const char *fn_ref);
	void vcf_close(vcfFile *fp);
	vcf_hdr_t *vcf_hdr_read(vcfFile *fp);
	void vcf_hdr_destroy(vcf_hdr_t *h);

	vcf1_t *vcf_init1(void);
	void vcf_destroy1(vcf1_t *v);
	int vcf_read1(vcfFile *fp, const vcf_hdr_t *h, vcf1_t *v);
	void vcf_write1(vcfFile *fp, const vcf1_t *v);

#ifdef __cplusplus
}
#endif

#endif
