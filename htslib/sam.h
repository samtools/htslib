#ifndef BAM_H
#define BAM_H

#include <stdint.h>
#include "hts.h"

/**********************
 *** SAM/BAM header ***
 **********************/

/*! @typedef
 @abstract Structure for the alignment header.
 @field n_targets   number of reference sequences
 @field l_text      length of the plain text in the header
 @field target_len  lengths of the referene sequences
 @field target_name names of the reference sequences
 @field text        plain text
 @field sdict       header dictionary
 */

typedef struct {
	int32_t n_targets, ignore_sam_err;
	uint32_t l_text;
	uint32_t *target_len;
	uint8_t *cigar_tab;
	char **target_name;
	char *text;
	void *sdict;
} bam_hdr_t;

/****************************
 *** CIGAR related macros ***
 ****************************/

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

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

/*************************
 *** Alignment records ***
 *************************/

/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 */
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

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
 
 @discussion Notes:
 
 1. qname is zero tailing and core.l_qname includes the tailing '\0'.
 2. l_qseq is calculated from the total length of an alignment block
 on reading or from CIGAR.
 3. cigar data is encoded 4 bytes per CIGAR operation.
 4. seq is nybble-encoded according to bam_nt16_table.
 */
typedef struct {
	bam1_core_t core;
	int l_data, m_data;
	uint8_t *data;
#ifndef BAM_NO_ID
	uint64_t id;
#endif
} bam1_t;

#define bam_is_rev(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam_is_mrev(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define bam_get_qname(b) ((char*)(b)->data)
#define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
#define bam_get_l_aux(b) ((b)->l_data - ((b)->core.n_cigar<<2) - (b)->core.l_qname - (b)->core.l_qseq - (((b)->core.l_qseq + 1)>>1))
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

/**************************
 *** Exported functions ***
 **************************/

#ifdef __cplusplus
extern "C" {
#endif

	/***************
	 *** BAM I/O ***
	 ***************/

	bam_hdr_t *bam_hdr_init();
	bam_hdr_t *bam_hdr_read(BGZF *fp);
	int bam_hdr_write(BGZF *fp, const bam_hdr_t *h);
	void bam_hdr_destroy(bam_hdr_t *h);
	int bam_name2id(bam_hdr_t *h, const char *ref);

	bam1_t *bam_init1(void);
	void bam_destroy1(bam1_t *b);
	int bam_read1(BGZF *fp, bam1_t *b);
	int bam_write1(BGZF *fp, const bam1_t *b);
	bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc);

	/* Internal helper function used by bam_itr_next() */
	int bam_readrec(BGZF *fp, void *null, bam1_t *b, int *tid, int *beg, int *end);

	int bam_cigar2qlen(int n_cigar, const uint32_t *cigar);
	int bam_cigar2rlen(int n_cigar, const uint32_t *cigar);

	/*!
	  @abstract Calculate the rightmost base position of an alignment on the
	  reference genome.

	  @param  b  pointer to an alignment
	  @return    the coordinate of the first base after the alignment, 0-based

	  @discussion For a mapped read, this is just b->core.pos + bam_cigar2rlen.
	  For an unmapped read (either according to its flags or if it has no cigar
	  string), we return b->core.pos + 1 by convention.
	*/
	int32_t bam_endpos(const bam1_t *b);

	/********************
	 *** BAM indexing ***
	 ********************/

	#define bam_itr_destroy(iter) hts_itr_destroy(iter)
	#define bam_itr_queryi(idx, tid, beg, end) hts_itr_query(idx, tid, beg, end)
	#define bam_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bam_name2id), (hdr))
	#define bam_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), (hts_readrec_f)(bam_readrec), 0)
	#define bam_index_load(fn) hts_idx_load((fn), HTS_FMT_BAI)

	int bam_index_build(const char *fn, int min_shift);

	/***************
	 *** SAM I/O ***
	 ***************/

	#define sam_open(fn, mode) (hts_open((fn), (mode)))
	#define sam_close(fp) hts_close(fp)

	typedef htsFile samFile;
	bam_hdr_t *sam_hdr_parse(int l_text, const char *text);
	bam_hdr_t *sam_hdr_read(samFile *fp);
	int sam_hdr_write(samFile *fp, const bam_hdr_t *h);

	int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b);
	int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);
	int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
	int sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b);

	/*************************************
	 *** Manipulating auxiliary fields ***
	 *************************************/

	uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);
	int32_t bam_aux2i(const uint8_t *s);
	double bam_aux2f(const uint8_t *s);
	char bam_aux2A(const uint8_t *s);
	char *bam_aux2Z(const uint8_t *s);

	void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data);
	int bam_aux_del(bam1_t *b, uint8_t *s);

#ifdef __cplusplus
}
#endif

/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

typedef struct {
	bam1_t *b;
	int32_t qpos;
	int indel, level;
	uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

struct __bam_plp_t;
typedef struct __bam_plp_t *bam_plp_t;

struct __bam_mplp_t;
typedef struct __bam_mplp_t *bam_mplp_t;

#ifdef __cplusplus
extern "C" {
#endif

	bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
	void bam_plp_destroy(bam_plp_t iter);
	int bam_plp_push(bam_plp_t iter, const bam1_t *b);
	const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
	const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
	void bam_plp_set_mask(bam_plp_t iter, int mask);
	void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
	void bam_plp_reset(bam_plp_t iter);

	bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
	void bam_mplp_destroy(bam_mplp_t iter);
	void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
	int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);

#ifdef __cplusplus
}
#endif

#endif // ~!defined(BAM_NO_PILEUP)

#endif
