#ifndef SAM_H
#define SAM_H

#include <stdint.h>
#include "hts.h"

/******************
 * SAM/BAM header *
 ******************/

typedef struct {
	int32_t n_targets;
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

#define SAM_CIGAR_STR   "MIDNSHP=XB"
#define SAM_CIGAR_SHIFT 4
#define SAM_CIGAR_MASK  0xf
#define SAM_CIGAR_TYPE  0x3C1A7

#define sam_cigar_op(c) ((c)&SAM_CIGAR_MASK)
#define sam_cigar_oplen(c) ((c)>>SAM_CIGAR_SHIFT)
#define sam_cigar_opchr(c) (SAM_CIGAR_STR[sam_cigar_op(c)])
#define sam_cigar_gen(l, o) ((l)<<SAM_CIGAR_SHIFT|(o))
#define sam_cigar_type(o) (SAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*********************
 * Alignment records *
 *********************/

typedef struct {
	int32_t tid;
	int32_t pos;
	uint32_t bin:16, qual:8, l_qname:8;
	uint32_t flag:16, n_cigar:16;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} sam1_core_t;

typedef struct {
	sam1_core_t core;
	int l_aux, l_data, m_data;
	uint8_t *data;
} sam1_t;

#define sam_get_strand(b) (((b)->core.flag&SAM_FREVERSE) != 0)
#define sam_get_mstrand(b) (((b)->core.flag&SAM_FMREVERSE) != 0)
#define sam_get_qname(b) ((char*)(b)->data)
#define sam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define sam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define sam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define sam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
#define sam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

#ifdef __cplusplus
extern "C" {
#endif

	void sam_hdr_destroy(sam_hdr_t *h);
	sam_hdr_t *sam_hdr_read(htsFile *fp);
	int sam_hdr_write(htsFile *fp, const sam_hdr_t *h);

	sam1_t *sam_init1(void);
	void sam_destroy1(sam1_t *b);
	int sam_read1(htsFile *fp, const sam_hdr_t *h, sam1_t *b);
	int sam_format1(const sam_hdr_t *h, const sam1_t *b, kstring_t *str);
	int sam_write1(htsFile *fp, const sam_hdr_t *h, const sam1_t *b);

#ifdef __cplusplus
}
#endif

#endif
