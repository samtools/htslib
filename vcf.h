#ifndef VCF_H
#define VCF_H

#include <stdint.h>
#include <limits.h>
#include <assert.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*******************
 * VCF file struct *
 *******************/

typedef struct {
	uint32_t is_bin:1, is_write:1, dummy:30;
	kstring_t line;
	char *fn_ref; // external reference sequence dictionary
	void *fp; // file pointer; actual type depending on is_bin and is_write
} vcfFile;

/*****************
 * Header struct *
 *****************/

#define VCF_HL_FLT  0 // header line
#define VCF_HL_INFO 1
#define VCF_HL_FMT  2
#define VCF_HL_CTG  3

#define VCF_HT_FLAG 0 // header type
#define VCF_HT_INT  1
#define VCF_HT_REAL 2
#define VCF_HT_STR  3

#define VCF_VL_FIXED 0 // variable length
#define VCF_VL_VAR   1
#define VCF_VL_A     2
#define VCF_VL_G     3

/* === Dictionary ===

   The header keeps three dictonaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
   in the "contig" lines and the last keeps the sample names. vcf_hdr_t::dict[]
   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   vcf_idinfo_t struct. vcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. vcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/

#define VCF_DT_ID		0 // dictionary type
#define VCF_DT_CTG		1
#define VCF_DT_SAMPLE	2

typedef struct {
	uint32_t info[3]; // Number:20, var:4, Type:4, ColType:4
	int id;
} vcf_idinfo_t;

typedef struct {
	const char *key;
	const vcf_idinfo_t *val;
} vcf_idpair_t;

typedef struct {
	int32_t l_text, n[3];
	vcf_idpair_t *id[3];
	void *dict[3]; // ID dictionary, contig dict and sample dict
	char *text;
	kstring_t mem;
} vcf_hdr_t;

extern uint8_t vcf_type_shift[];

#define vcf_hdr_n_val(x, _n_alt) (((x)>>8&0xf) == VCF_VTP_FIXED? (x)>>12 : ((x)>>8&0xf) == VCF_VTP_A? (_n_alt) : ((x)>>8&0xf) == VCF_VTP_G? -2 : -1)

/**************
 * VCF record *
 **************/

#define VCF_BT_NULL		0
#define VCF_BT_INT8		1
#define VCF_BT_INT16	2
#define VCF_BT_INT32	3
#define VCF_BT_FLOAT	5
#define VCF_BT_CHAR		7

typedef struct {
	int32_t rid; // CHROM
	int32_t pos; // POS
	float qual;  // QUAL
	kstring_t shared, indiv;
} vcf1_t;

typedef struct {
	int id, n, type, size;
	uint8_t *p;
} vcf_fmt_t;

/*******
 * API *
 *******/

#ifdef __cplusplus
extern "C" {
#endif

	vcfFile *vcf_open(const char *fn, const char *mode, const char *fn_ref);
	void vcf_close(vcfFile *fp);
	vcf_hdr_t *vcf_hdr_read(vcfFile *fp);
	void vcf_hdr_write(vcfFile *fp, const vcf_hdr_t *h);
	void vcf_hdr_destroy(vcf_hdr_t *h);

	vcf1_t *vcf_init1(void);
	void vcf_destroy1(vcf1_t *v);
	int vcf_read1(vcfFile *fp, const vcf_hdr_t *h, vcf1_t *v);
	int vcf_format1(const vcf_hdr_t *h, const vcf1_t *v, kstring_t *s);
	int vcf_write1(vcfFile *fp, const vcf_hdr_t *h, const vcf1_t *v);

	int vcf_id2int(const vcf_hdr_t *h, int which, const char *id);
	vcf_fmt_t *vcf_unpack_fmt(const vcf_hdr_t *h, const vcf1_t *v, int *n_fmt);

#ifdef __cplusplus
}
#endif

/*******************
 * Typed value I/O *
 *******************/

#include "kstring.h"

static inline void vcf_enc_size(kstring_t *s, int size, int type)
{
	if (size >= 15) {
		kputc(15<<4|type, s);
		if (size >= 128) {
			int16_t x = size;
			assert(size <= 32767);
			kputc(1<<4|VCF_BT_INT16, s);
			kputsn((char*)&x, 2, s);
		} else {
			kputc(1<<4|VCF_BT_INT8, s);
			kputc(size, s);
		}
	} else kputc(size<<4|type, s);
}

static inline int vcf_enc_inttype(long x)
{
	if (x <= INT8_MAX && x > INT8_MIN) return VCF_BT_INT8;
	if (x <= INT16_MAX && x > INT16_MIN) return VCF_BT_INT16;
	return VCF_BT_INT32;
}

static inline void vcf_enc_int1(kstring_t *s, int32_t x)
{
	if (x == INT32_MIN) {
		vcf_enc_size(s, 1, VCF_BT_INT8);
		kputc(INT8_MIN, s);
	} else if (x <= INT8_MAX && x > INT8_MIN) {
		vcf_enc_size(s, 1, VCF_BT_INT8);
		kputc(x, s);
	} else if (x <= INT16_MAX && x > INT16_MIN) {
		int16_t z = x;
		vcf_enc_size(s, 1, VCF_BT_INT16);
		kputsn((char*)&z, 2, s);
	} else {
		int32_t z = x;
		vcf_enc_size(s, 1, VCF_BT_INT32);
		kputsn((char*)&z, 4, s);
	}
}

static inline int32_t vcf_dec_int1(const uint8_t *p, int type, uint8_t **q)
{
	if (type == VCF_BT_INT8) {
		*q = (uint8_t*)p + 1;
		return *(int8_t*)p;
	} else if (type == VCF_BT_INT16) {
		*q = (uint8_t*)p + 2;
		return *(int16_t*)p;
	} else {
		*q = (uint8_t*)p + 4;
		return *(int32_t*)p;
	}
}

static inline int32_t vcf_dec_typed_int1(const uint8_t *p, uint8_t **q)
{
	return vcf_dec_int1(p + 1, *p&0xf, q);
}

static inline int32_t vcf_dec_size(const uint8_t *p, uint8_t **q, int *type)
{
	*type = *p & 0xf;
	if (*p>>4 != 15) {
		*q = (uint8_t*)p + 1;
		return *p>>4;
	} else return vcf_dec_typed_int1(p + 1, q);
}

#endif
