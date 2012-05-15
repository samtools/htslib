#ifndef BAM_H
#define BAM_H

#include <stdint.h>
#include "hts.h"

/******************
 * SAM/BAM header *
 ******************/

typedef struct {
	int32_t n_targets;
	uint32_t l_text;
	uint32_t *target_len;
	uint8_t *cigar_tab;
	char **target_name;
	char *text;
	void *sdict;
} bam_hdr_t;

/************************
 * CIGAR related macros *
 ************************/

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024

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
} bam1_core_t;

typedef struct {
	bam1_core_t core;
	int l_data, m_data;
	uint8_t *data;
} bam1_t;

#define bam_get_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam_get_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define bam_get_qname(b) ((char*)(b)->data)
#define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

/**********************
 * Exported functions *
 **********************/

#ifdef __cplusplus
extern "C" {
#endif

	/***************
	 *** BAM I/O ***
	 ***************/

	bam_hdr_t *bam_hdr_read(void *fp);
	int bam_hdr_write(void *fp, const bam_hdr_t *h);
	void bam_hdr_destroy(bam_hdr_t *h);
	int bam_get_tid(bam_hdr_t *h, const char *ref);

	bam1_t *bam_init1(void);
	void bam_destroy1(bam1_t *b);
	int bam_read1(void *fp, bam1_t *b);
	int bam_write1(void *fp, const bam1_t *b);

	/********************
	 *** BAM indexing ***
	 ********************/

	#define bam_iter_queryi(idx, tid, beg, end) hts_iter_query(idx, tid, beg, end)
	#define bam_iter_destroy(iter) hts_iter_destroy(iter)

	typedef hts_idx_t bam_idx_t;
	typedef hts_iter_t bam_iter_t;

	bam_idx_t *bam_index(void *fp);
	int bam_index_build(const char *fn, const char *_fnidx);
	bam_idx_t *bam_index_load_local(const char *fnidx);
	bam_iter_t *bam_iter_querys(hts_idx_t *idx, bam_hdr_t *h, const char *reg);
	int bam_iter_read(void *fp, bam_iter_t *iter, bam1_t *b);

	/***************
	 *** SAM I/O ***
	 ***************/

	#define sam_open(fn, mode, fnaux) hts_open(fn, mode, fnaux)
	#define sam_close(fp) hts_close(fp)

	typedef htsFile samFile;
	bam_hdr_t *sam_hdr_parse(int l_text, const char *text);
	bam_hdr_t *sam_hdr_read(samFile *fp);
	int sam_hdr_write(samFile *fp, const bam_hdr_t *h);

	int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b);
	int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);
	int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
	int sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b);

#ifdef __cplusplus
}
#endif

#endif
