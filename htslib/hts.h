#ifndef HTS_H
#define HTS_H

#include <stdint.h>
#include "bgzf.h"

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
extern char seq_nt16_str[];

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
	char *hts_idx_getfn(const char *fn, const char *ext);

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

/**************
 * Endianness *
 **************/

static inline int ed_is_big()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
	*(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
	*(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
	*(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
	return x;
}

#endif
