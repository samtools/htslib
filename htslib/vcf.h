#ifndef BCF_H
#define BCF_H

#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "bgzf.h"
#include "hts.h"

/*****************
 * Header struct *
 *****************/

#define BCF_HL_FLT  0 // header line
#define BCF_HL_INFO 1
#define BCF_HL_FMT  2
#define BCF_HL_CTG  3

#define BCF_HT_FLAG 0 // header type
#define BCF_HT_INT  1
#define BCF_HT_REAL 2
#define BCF_HT_STR  3

#define BCF_VL_FIXED 0 // variable length
#define BCF_VL_VAR   1
#define BCF_VL_A     2
#define BCF_VL_G     3

/* === Dictionary ===

   The header keeps three dictonaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
   in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/

#define BCF_DT_ID		0 // dictionary type
#define BCF_DT_CTG		1
#define BCF_DT_SAMPLE	2

typedef struct {
	uint32_t info[3]; // for each number => Number:20, var:4, Type:4, ColType:4
	int id;
} bcf_idinfo_t;

typedef struct {
	const char *key;
	const bcf_idinfo_t *val;
} bcf_idpair_t;

typedef struct {
	int32_t l_text, n[3];
	bcf_idpair_t *id[3];
	void *dict[3]; // ID dictionary, contig dict and sample dict
	char *text;
	kstring_t mem;
} bcf_hdr_t;

extern uint8_t bcf_type_shift[];

/**************
 * VCF record *
 **************/

#define BCF_BT_NULL		0
#define BCF_BT_INT8		1
#define BCF_BT_INT16	2
#define BCF_BT_INT32	3
#define BCF_BT_FLOAT	5
#define BCF_BT_CHAR		7

typedef struct {
	int id, n, type, size;
	uint8_t *p;
} bcf_fmt_t;

typedef struct {
	int key, type, len;
	union {
		int32_t i;
		float f;
	} v1; // only set if $len==1
	uint8_t *vptr;
} bcf_info_t;

typedef struct {
	int m_fmt, m_info, m_str, m_allele, m_flt, n_flt;
	char *id, **allele;
	int *flt;
	bcf_fmt_t *fmt;
	bcf_info_t *info;
} bcf_dec_t;

typedef struct {
	int32_t rid;  // CHROM
	int32_t pos;  // POS
	int32_t rlen; // length of REF
	float qual;   // QUAL
	uint32_t n_info:16, n_allele:16;
	uint32_t n_fmt:8, n_sample:24;
	kstring_t shared, indiv;
	bcf_dec_t d; // lazy evaluation: $d is filled by bcf_unpack()
} bcf1_t;

/*******
 * API *
 *******/

#ifdef __cplusplus
extern "C" {
#endif

	/***************
	 *** BCF I/O ***
	 ***************/

	bcf_hdr_t *bcf_hdr_read(BGZF *fp);
	void bcf_hdr_write(BGZF *fp, const bcf_hdr_t *h);
	void bcf_hdr_destroy(bcf_hdr_t *h);

	bcf1_t *bcf_init1();
	void bcf_destroy1(bcf1_t *v);
	int bcf_read1(BGZF *fp, bcf1_t *v);
	int bcf_write1(BGZF *fp, const bcf1_t *v);
	int bcf_readrec(BGZF *fp, void *null, bcf1_t *v, int *tid, int *beg, int *end);

	int bcf_id2int(const bcf_hdr_t *h, int which, const char *id);
	int bcf_name2id(const bcf_hdr_t *h, const char *id);
	bcf_fmt_t *bcf_unpack_fmt(const bcf_hdr_t *h, const bcf1_t *v);
    int bcf_unpack(bcf1_t *b);
	uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr);

	/*****************
	 *** BCF index ***
	 *****************/

	#define bcf_itr_destroy(iter) hts_itr_destroy(iter)
	#define bcf_itr_queryi(idx, tid, beg, end) hts_itr_query((idx), (tid), (beg), (end))
	#define bcf_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bcf_name2id), (hdr))
	#define bcf_itr_next(fp, itr, r) hts_itr_next((fp), (itr), (r), (hts_readrec_f)(bcf_readrec), 0)
	#define bcf_index_load(fn) hts_idx_load(fn, HTS_FMT_CSI)

	int bcf_index_build(const char *fn, int min_shift);

	/***************
	 *** VCF I/O ***
	 ***************/
	
	bcf_hdr_t *vcf_hdr_read(htsFile *fp);
	void vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h);

	int vcf_parse1(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v);
	int vcf_format1(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);
	int vcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);
	int vcf_write1(htsFile *fp, const bcf_hdr_t *h, const bcf1_t *v);

	bcf_hdr_t *vcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);
	int vcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap);

#ifdef __cplusplus
}
#endif

/*******************
 * Typed value I/O *
 *******************/

#include "kstring.h"

static inline void bcf_enc_size(kstring_t *s, int size, int type)
{
	if (size >= 15) {
		kputc(15<<4|type, s);
		if (size >= 128) {
			if (size >= 32768) {
				int32_t x = size;
				kputc(1<<4|BCF_BT_INT32, s);
				kputsn((char*)&x, 4, s);
			} else {
				int16_t x = size;
				kputc(1<<4|BCF_BT_INT16, s);
				kputsn((char*)&x, 2, s);
			}
		} else {
			kputc(1<<4|BCF_BT_INT8, s);
			kputc(size, s);
		}
	} else kputc(size<<4|type, s);
}

static inline int bcf_enc_inttype(long x)
{
	if (x <= INT8_MAX && x > INT8_MIN) return BCF_BT_INT8;
	if (x <= INT16_MAX && x > INT16_MIN) return BCF_BT_INT16;
	return BCF_BT_INT32;
}

static inline void bcf_enc_int1(kstring_t *s, int32_t x)
{
	if (x == INT32_MIN) {
		bcf_enc_size(s, 1, BCF_BT_INT8);
		kputc(INT8_MIN, s);
	} else if (x <= INT8_MAX && x > INT8_MIN) {
		bcf_enc_size(s, 1, BCF_BT_INT8);
		kputc(x, s);
	} else if (x <= INT16_MAX && x > INT16_MIN) {
		int16_t z = x;
		bcf_enc_size(s, 1, BCF_BT_INT16);
		kputsn((char*)&z, 2, s);
	} else {
		int32_t z = x;
		bcf_enc_size(s, 1, BCF_BT_INT32);
		kputsn((char*)&z, 4, s);
	}
}

static inline int32_t bcf_dec_int1(const uint8_t *p, int type, uint8_t **q)
{
	if (type == BCF_BT_INT8) {
		*q = (uint8_t*)p + 1;
		return *(int8_t*)p;
	} else if (type == BCF_BT_INT16) {
		*q = (uint8_t*)p + 2;
		return *(int16_t*)p;
	} else {
		*q = (uint8_t*)p + 4;
		return *(int32_t*)p;
	}
}

static inline int32_t bcf_dec_typed_int1(const uint8_t *p, uint8_t **q)
{
	return bcf_dec_int1(p + 1, *p&0xf, q);
}

static inline int32_t bcf_dec_size(const uint8_t *p, uint8_t **q, int *type)
{
	*type = *p & 0xf;
	if (*p>>4 != 15) {
		*q = (uint8_t*)p + 1;
		return *p>>4;
	} else return bcf_dec_typed_int1(p + 1, q);
}

#endif
