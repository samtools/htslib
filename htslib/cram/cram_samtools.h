#ifndef _CRAM_SAMTOOLS_H_
#define _CRAM_SAMTOOLS_H_

/* Samtools compatible API */
#define bam_blk_size(b)  ((b)->l_data)
#define bam_set_blk_size(b,v) ((b)->data_len = (v))

#define bam_ref(b)       (b)->core.tid
#define bam_pos(b)       (b)->core.pos
#define bam_mate_pos(b)  (b)->core.mpos
#define bam_mate_ref(b)  (b)->core.mtid
#define bam_ins_size(b)  (b)->core.isize
#define bam_seq_len(b)   (b)->core.l_qseq
#define bam_cigar_len(b) (b)->core.n_cigar
#define bam_flag(b)      (b)->core.flag
#define bam_bin(b)       (b)->core.bin
#define bam_map_qual(b)  (b)->core.qual
#define bam_name_len(b)  (b)->core.l_qname
#define bam_name(b)      bam_get_qname((b))
#define bam_qual(b)      bam_get_qual((b))
#define bam_seq(b)       bam_get_seq((b))
#define bam_cigar(b)     bam_get_cigar((b))
#define bam_aux(b)       bam_get_aux((b))

#define bam_dup(b)       bam_copy1(bam_init1(), (b))

#define bam_free(b)      bam_destroy1((b))

#define bam_reg2bin(beg,end) hts_reg2bin((beg),(end),14,5)

#include "sam.h"

enum cigar_op {
    BAM_CMATCH_=BAM_CMATCH,
    BAM_CINS_=BAM_CINS,
    BAM_CDEL_=BAM_CDEL,
    BAM_CREF_SKIP_=BAM_CREF_SKIP,
    BAM_CSOFT_CLIP_=BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP_=BAM_CHARD_CLIP,
    BAM_CPAD_=BAM_CPAD,
    BAM_CBASE_MATCH=BAM_CEQUAL,
    BAM_CBASE_MISMATCH=BAM_CDIFF
};

typedef bam1_t bam_seq_t;

// Sorry this is a bit loopy
#include "cram/cram.h"

int cram_get_bam1_seq(cram_fd *fd, bam1_t *b);
bam_hdr_t *cram_header_to_bam(SAM_hdr *h);
SAM_hdr *bam_header_to_cram(bam_hdr_t *h);

int bam_construct_seq(bam_seq_t **bp, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end,        // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual);

#endif /* _CRAM_SAMTOOLS_H_ */
