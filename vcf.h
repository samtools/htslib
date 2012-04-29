#ifndef VCF_H
#define VCF_H

#include <stdint.h>

typedef struct {
	uint32_t is_vcf:1, is_write:1, dummy:30;
	union {
		void *fp_bcf;
		void *fp_vcf;
	} x;
} vcfFile;

typedef struct {
	uint32_t info[3]; // Number:20, var:4, Type:4, ColType:4
	int kid, rid, sid, rlen; // key-ID, ref-ID, sample-ID, ref-length
} vcf_keyinfo_t;

typedef struct {
	const char *key;
	const vcf_keyinfo_t *info;
} vcf_keypair_t;

typedef struct {
	int32_t n_ref, n_sample, n_key;
	int32_t l_str, m_str;
	int32_t l_text;
	vcf_keypair_t *key;
	int *r2k, *s2k;
	char *text;
	void *dict;
} vcf_hdr_t;

typedef struct {
	int32_t rid; // CHROM
	int32_t pos; // POS
	int32_t end; // end position
	float qual;  // QUAL
	int32_t n_info, n_fmt;
	int32_t l_mem, m_mem;
	char *id, *ref, *alt; // point to mem
	uint8_t *flt; // point to mem
	uint8_t *mem;
} vcf1_t;

#ifdef __cplusplus
extern "C" {
#endif

	vcfFile *vcf_open(const char *fn, const char *mode, const char *fn_ref);
	void vcf_close(vcfFile *fp);
	vcf_hdr_t *vcf_hdr_read(vcfFile *fp);
	void vcf_hdr_destroy(vcf_hdr_t *h);

	vcf1_t *vcf_init1(void);
	void vcf_destroy1(vcf1_t *v);
	void vcf_read1(vcfFile *fp, vcf_hdr_t *h, vcf1_t *v);
	void vcf_write1(vcfFile *fp, const vcf1_t *v);

#ifdef __cplusplus
}
#endif

#endif
