#ifndef SAM_H
#define SAM_H

#include <stdint.h>

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
 * SAM file struct *
 *******************/

typedef struct { // identical to vcfFile
	uint32_t is_bin:1, is_write:1, dummy:30;
	kstring_t line;
	char *fn_ref; // external reference sequence dictionary
	void *fp; // file pointer; actual type depending on is_bin and is_write
} samFile;

/******************
 * SAM/BAM header *
 ******************/

typedef struct {
	int32_t n_targets, is_be;
	uint32_t l_text, n_text;
	uint32_t *target_len;
	char **target_name;
	char *text;
	void *dict;
} sam_hdr_t;

/************************
 * CIGAR related macros *
 ************************/

#define SAM_CMATCH      0
#define SAM_CINS        1
#define SAM_CDEL        2
#define SAM_CREF_SKIP   3
#define SAM_CSOFT_CLIP  4
#define SAM_CHARD_CLIP  5
#define SAM_CPAD        6
#define SAM_CEQUAL      7
#define SAM_CDIFF       8
#define SAM_CBACK       9

#define SAM_CIGAR_STR  "MIDNSHP=XB"
#define SAM_CIGAR_TYPE 0x3C1A7

#define sam_cigar_op(c) ((c)&SAM_CIGAR_MASK)
#define sam_cigar_oplen(c) ((c)>>SAM_CIGAR_SHIFT)
#define sam_cigar_opchr(c) (SAM_CIGAR_STR[sam_cigar_op(c)])
#define sam_cigar_gen(l, o) ((l)<<SAM_CIGAR_SHIFT|(o))
#define sam_cigar_type(o) (SAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*********************
 * Alignment records *
 *********************/

typedef struct {
	int32_t rid;
	int32_t pos;
	uint32_t bin:16, qual:8, l_qname:8;
	uint32_t flag:16, n_cigar:16;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t tlen;
} sam1_core_t;

typedef struct {
	sam1_core_t core;
	int l_aux, l_data, m_data;
	uint8_t *data;
} sam1_t;

#define sam_get_strand(b) (((b)->core.flag&SAM_FREVERSE) != 0)
#define sam_get_mstrand(b) (((b)->core.flag&SAM_FMREVERSE) != 0)
#define sam_get_qname(b) ((char*)((b)->data))
#define sam1_seq(b) ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define sam1_qual(b) ((b)->data + ((b)->core.n_cigar<<1) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define sam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

#endif
