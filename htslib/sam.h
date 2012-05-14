#ifndef SAM_H
#define SAM_H

#include <stdint.h>
#include "hts.h"

/******************
 * SAM/BAM header *
 ******************/

typedef struct {
	int32_t n_targets, has_SQ;
	uint32_t l_text;
	uint32_t *target_len;
	uint8_t *cigar_tab;
	char **target_name;
	char *text;
	void *sdict;
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

#define SAM_FPAIRED        1
#define SAM_FPROPER_PAIR   2
#define SAM_FUNMAP         4
#define SAM_FMUNMAP        8
#define SAM_FREVERSE      16
#define SAM_FMREVERSE     32
#define SAM_FREAD1        64
#define SAM_FREAD2       128
#define SAM_FSECONDARY   256
#define SAM_FQCFAIL      512
#define SAM_FDUP        1024

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
	int l_data, m_data;
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

/************************************
 * Alias of hts types and functions *
 ************************************/

typedef htsFile samFile;
typedef hts_idx_t sam_idx_t;

#define sam_open(fn, mode, fnaux) hts_open(fn, mode, fnaux)
#define sam_close(fp) hts_close(fp)
#define sam_iter_queryi(idx, tid, beg, end) hts_iter_query(idx, tid, beg, end)

/**********************
 * Exported functions *
 **********************/

#ifdef __cplusplus
extern "C" {
#endif

	void sam_hdr_destroy(sam_hdr_t *h);
	sam_hdr_t *sam_hdr_read(htsFile *fp);
	int sam_hdr_write(htsFile *fp, const sam_hdr_t *h);
	int sam_get_tid(sam_hdr_t *h, const char *ref);

	sam1_t *sam_init1(void);
	void sam_destroy1(sam1_t *b);
	int sam_read1(htsFile *fp, sam_hdr_t *h, sam1_t *b);
	int sam_format1(const sam_hdr_t *h, const sam1_t *b, kstring_t *str);
	int sam_write1(htsFile *fp, const sam_hdr_t *h, const sam1_t *b);

	hts_idx_t *sam_index(htsFile *fp);
	int sam_index_build(const char *fn, const char *_fnidx);
	hts_idx_t *sam_index_load_local(const char *fnidx);
	hts_iter_t *sam_iter_querys(hts_idx_t *idx, sam_hdr_t *h, const char *reg);
	int sam_iter_read(samFile *fp, hts_iter_t *iter, sam1_t *b);

#ifdef __cplusplus
}
#endif

#endif
