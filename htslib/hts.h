#ifndef HTS_H
#define HTS_H

#define HTS_VERSION "r199"

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

/**
 * hts_expand()  - expands memory block pointed to by $ptr;
 * hts_expand0()   the latter sets the newly allocated part to 0.
 *
 * @param n     requested number of elements of type type_t
 * @param m     size of memory allocated
 */
#define hts_expand(type_t, n, m, ptr) if ((n) > (m)) { \
		(m) = (n); kroundup32(m); \
		(ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
	}
#define hts_expand0(type_t, n, m, ptr) if ((n) > (m)) { \
		int t = (m); (m) = (n); kroundup32(m); \
		(ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
        memset((ptr)+t,0,sizeof(type_t)*((m)-t)); \
	}

/************
 * File I/O *
 ************/

typedef struct {
	uint32_t is_bin:1, is_write:1, is_be:1, dummy:29;
	int64_t lineno;
	kstring_t line;
	char *fn, *fn_aux;
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
	char **hts_readlines(const char *fn, int *_n);

#ifdef __cplusplus
}
#endif

/************
 * Indexing *
 ************/

#define HTS_IDX_NOCOOR (-2)
#define HTS_IDX_START  (-3)
#define HTS_IDX_REST   (-4)

#define HTS_FMT_CSI 0
#define HTS_FMT_BAI 1
#define HTS_FMT_TBI 2

struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
	uint64_t u, v;
} hts_pair64_t;

typedef struct {
	uint32_t read_rest:1, finished:1, dummy:29;
	int tid, beg, end, n_off, i;
	uint64_t curr_off;
	hts_pair64_t *off;
	struct {
		int n, m;
		int *a;
	} bins;
} hts_itr_t;

#ifdef __cplusplus
extern "C" {
#endif

	#define hts_bin_first(l) (((1<<(((l)<<1) + (l))) - 1) / 7)
	#define hts_bin_parent(l) (((l) - 1) >> 3)

	hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
	void hts_idx_destroy(hts_idx_t *idx);
	int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped);
	void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);

	void hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt);
	hts_idx_t *hts_idx_load(const char *fn, int fmt);

	uint8_t *hts_idx_get_meta(hts_idx_t *idx, int *l_meta);
	void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy);

	const char *hts_parse_reg(const char *s, int *beg, int *end);
	hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end);
	void hts_itr_destroy(hts_itr_t *iter);

	typedef int (*hts_readrec_f)(BGZF*, void*, void*, int*, int*, int*);
	typedef int (*hts_name2id_f)(void*, const char*);

	hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr);
	int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, hts_readrec_f readrec, void *hdr);

    /**
     * file_type() - Convenience function to determine file type
     * @fname: the file name
     *
     * Returns one of the IS_* defines.
     *
     * This function was added in order to avoid the need for excessive command
     * line switches. Note that in the current implementation only the file name is
     * checked, looking at magic string is also a possibility.
     */
    #define IS_VCF    1
    #define IS_VCF_GZ 2
    #define IS_BCF    4
    int file_type(const char *fname);


#ifdef __cplusplus
}
#endif

static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
{
	int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
	for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
		if (beg>>s == end>>s) return t + (beg>>s);
	return 0;
}

static inline int hts_bin_bot(int bin, int n_lvls)
{
	int l, b;
	for (l = 0, b = bin; b; ++l, b = hts_bin_parent(b)); // compute the level of bin
	return (bin - hts_bin_first(l)) << (n_lvls - l) * 3;
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
