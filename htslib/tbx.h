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
	int n, m;
	tbx_conf_t conf;
	hts_idx_t *idx;
	void *dict;
} tbx_t;

#ifdef __cplusplus
extern "C" {
#endif

	int tbx_name2id(tbx_t *tbx, const char *ss);
	int tbx_readrec(BGZF *fp, tbx_t *tbx, kstring_t *s, int *tid, int *beg, int *end);
	int tbx_index_build(const char *fn, const char *_fnidx, int min_shift, const tbx_conf_t *conf);
	tbx_t *tbx_index_load(const char *fn);

#define tbx_itr_querys(tbx, itr, s) hts_iter_querys((itr), (s), (hts_name2id_f)(tbx_name2id), (tbx))
#define tbx_itr_next(fp, tbx, itr, r) hts_iter_next((fp), (itr), (r), (hts_readrec_f)(tbx_readrec), (tbx))

#ifdef __cplusplus
}
#endif

#endif
