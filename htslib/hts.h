#ifndef HTS_H
#define HTS_H

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

/************
 * File I/O *
 ************/

typedef struct {
	uint32_t is_bin:1, is_write:1, is_be:1, dummy:29;
	int64_t lineno;
	kstring_t line;
	char *fn_aux;
	void *fp; // file pointer; actual type depending on is_bin and is_write
} htsFile;

/**********************
 * Exported functions *
 **********************/

extern int hts_verbose;
extern unsigned char seq_nt16_table[256];

#ifdef __cplusplus
extern "C" {
#endif

	htsFile *hts_open(const char *fn, const char *mode, const char *fn_aux);
	void hts_close(htsFile *fp);
	int hts_getline(htsFile *fp, int delimiter, kstring_t *str);

#ifdef __cplusplus
}
#endif

/************
 * Indexing *
 ************/

#define HTS_IDX_NOCOOR (-1)
#define HTS_IDX_START  (-2)

struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
	uint64_t u, v;
} hts_pair64_t;

typedef struct {
	uint32_t from_first:1, finished:1, dummy:30;
	int tid, beg, end, n_off, i;
	uint64_t curr_off;
	hts_pair64_t *off;
} hts_iter_t;

#ifdef __cplusplus
extern "C" {
#endif

	hts_idx_t *hts_idx_init(int n, uint64_t offset0);
	void hts_idx_destroy(hts_idx_t *idx);
	int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int bin, int is_mapped);
	void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);
	void hts_idx_save(const hts_idx_t *idx, void *fp, int is_bgzf);
	hts_idx_t *hts_idx_load(void *fp, int is_bgzf);
	const char *hts_parse_reg(const char *s, int *beg, int *end);
	hts_iter_t *hts_iter_query(const hts_idx_t *idx, int tid, int beg, int end);
	void hts_iter_destroy(hts_iter_t *iter);

#ifdef __cplusplus
}
#endif

static inline int hts_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

#endif
