#ifndef TBX_H
#define TBX_H

#include "hts.h"

#define TBX_MAX_SHIFT 31

#define TBX_GENERIC 0
#define TBX_SAM     1
#define TBX_VCF     2
#define TBX_UCSC    0x10000

typedef struct {
	int32_t preset;
	int32_t sc, bc, ec; // seq col., beg col. and end col.
	int32_t meta_char, line_skip;
} tbx_conf_t;

typedef struct {
	tbx_conf_t conf;
	hts_idx_t *idx;
	void *dict;
} tbx_t;

extern tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;

#ifdef __cplusplus
extern "C" {
#endif

	#define tbx_itr_destroy(iter) hts_itr_destroy(iter)
	#define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end))
	#define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx))
	#define tbx_itr_next(fp, tbx, itr, r) hts_itr_next((fp), (itr), (r), (hts_readrec_f)(tbx_readrec), (tbx))

	int tbx_name2id(tbx_t *tbx, const char *ss);
	int tbx_readrec(BGZF *fp, tbx_t *tbx, kstring_t *s, int *tid, int *beg, int *end);
	int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);
	tbx_t *tbx_index_load(const char *fn);
	void tbx_destroy(tbx_t *tbx);

#ifdef __cplusplus
}
#endif

#endif
