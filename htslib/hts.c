#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bgzf.h"
#include "endian.h"
#include "hts.h"

#include "kseq.h"
KSTREAM_INIT2(, gzFile, gzread, 16384)

#include "khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = 3;

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

htsFile *hts_open(const char *fn, const char *mode, const char *fn_aux)
{
	htsFile *fp;
	fp = (htsFile*)calloc(1, sizeof(htsFile));
	fp->is_be = ed_is_big();
	if (strchr(mode, 'w')) fp->is_write = 1;
	if (strchr(mode, 'b')) fp->is_bin = 1;
	if (fp->is_bin) {
		if (fp->is_write) fp->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_dopen(fileno(stdout), mode);
		else fp->fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
	} else {
		if (!fp->is_write) {
			gzFile gzfp;
			gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
			if (gzfp) fp->fp = ks_init(gzfp);
			if (fn_aux) fp->fn_aux = strdup(fn_aux);
		} else fp->fp = strcmp(fn, "-")? fopen(fn, "rb") : stdout;
	}
	if (fp->fp == 0) {
		if (hts_verbose >= 2)
			fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);
		free(fp->fn_aux); free(fp);
		return 0;
	}
	return fp;
}

void hts_close(htsFile *fp)
{
	if (!fp->is_bin) {
		free(fp->line.s);
		if (!fp->is_write) {
			gzFile gzfp = ((kstream_t*)fp->fp)->f;
			ks_destroy((kstream_t*)fp->fp);
			gzclose(gzfp);
			free(fp->fn_aux);
		} else fclose((FILE*)fp->fp);
	} else bgzf_close((BGZF*)fp->fp);
	free(fp);
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
	int ret, dret;
	ret = ks_getuntil((kstream_t*)fp->fp, delimiter, str, &dret);
	++fp->lineno;
	return ret;
}

/************
 * Indexing *
 ************/

#ifndef HTS_NO_INDEX

#define IDX_MIN_CHUNK_GAP 32768
#define IDX_LIDX_SHIFT    14
#define IDX_MAX_BIN       37450

typedef struct {
	uint64_t u, v;
} pair64_t;

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bins_t;

#include "khash.h"
KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} lidx_t;

struct __hts_index_t {
	int n, m;
	bidx_t **bidx;
	lidx_t *lidx;
	struct {
		uint32_t last_bin, save_bin;
		int last_coor, last_tid, save_tid, finished;
		uint64_t last_off, save_off, offset0;
		uint64_t off_beg, off_end;
		uint64_t n_no_coor, n_mapped, n_unmapped;
	} z; // keep internal states
};

static inline void insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
	khint_t k;
	bins_t *l;
	int absent;
	k = kh_put(bin, b, bin, &absent);
	l = &kh_value(b, k);
	if (absent) {
		l->m = 1; l->n = 0;
		l->list = (pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (pair64_t*)realloc(l->list, l->m * 16);
	}
	l->list[l->n].u = beg;
	l->list[l->n++].v = end;
}

static inline uint64_t insert_to_l(lidx_t *l, int _beg, int _end, uint64_t offset)
{
	int i, beg, end;
	beg = _beg >> IDX_LIDX_SHIFT;
	end = (_end - 1) >> IDX_LIDX_SHIFT;
	if (l->m < end + 1) {
		int old_m = l->m;
		l->m = end + 1;
		kroundup32(l->m);
		l->offset = (uint64_t*)realloc(l->offset, l->m * 8);
		memset(l->offset + old_m, 0, 8 * (l->m - old_m));
	}
	if (beg == end) { // to save a loop in this case
		if (l->offset[beg] == 0) l->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (l->offset[i] == 0) l->offset[i] = offset;
	}
	if (l->n < end + 1) l->n = end + 1;
	return (uint64_t)beg<<32 | end;
}

hts_index_t *hts_idx_init(uint64_t offset0)
{
	hts_index_t *idx;
	idx = (hts_index_t*)calloc(1, sizeof(hts_index_t));
	idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
	idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
	idx->z.last_coor = 0xffffffffu;
	idx->z.offset0 = (uint64_t)-1;
	return idx;
}

void hts_idx_finish(hts_index_t *idx)
{
	int i;
	if (idx->z.finished) return; // do not run this function multiple times
	if (idx->z.save_tid >= 0) {
		insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, idx->z.last_off);
		insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.off_beg, idx->z.last_off);
		insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.n_mapped, idx->z.n_unmapped);
	}
	for (i = 0; i < idx->n; ++i) {
		bidx_t *bidx = idx->bidx[i];
		lidx_t *lidx = &idx->lidx[i];
		khint_t k;
		int l, m;
		// merge adjacent blocks that start from the same BGZF block
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
			bins_t *p;
			if (!kh_exist(bidx, k)) continue;
			p = &kh_value(bidx, k);
			for (l = 1, m = 0; l < p->n; ++l) {
				if (p->list[m].v>>16 == p->list[l].u>>16) p->list[m].v = p->list[l].v; // if in the same BGZF block, merge them
				else p->list[++m] = p->list[l];
			}
			p->n = m + 1;
		}
		// fill missing values in the linear index
		for (l = 1; l < lidx->n; ++l)
			if (lidx->offset[l] == 0)
				lidx->offset[l] = lidx->offset[l-1];
	}
	// I forgot the purpose of the following block. It is in tabix and probably for a bug fix.
	if (idx->z.offset0 != (uint64_t)-1 && idx->n && idx->lidx[0].offset) {
		int beg = idx->z.offset0 >> 32, end = idx->z.offset0 & 0xfffffffu;
		for (i = beg; i < end; ++i) idx->lidx[0].offset[i] = 0;
	}
	idx->z.finished = 1;
}

int hts_idx_push(hts_index_t *idx, int tid, int beg, int end, uint64_t offset, int bin, int is_mapped)
{
	if (tid < 0) ++idx->z.n_no_coor;
	if (idx->z.finished) return 0;
	if (idx->z.last_tid < tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
		idx->z.last_tid = tid;
		idx->z.last_bin = 0xffffffffu;
	} else if ((uint32_t)idx->z.last_tid > (uint32_t)tid) {
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted chromosomes\n", __func__);
		return -1;
	} else if (tid >= 0 && idx->z.last_coor > beg) {
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted positions\n", __func__);
		return -1;
	}
	if (tid >= 0 && is_mapped) {
		uint64_t ret;
		ret = insert_to_l(&idx->lidx[tid], beg, end, idx->z.last_off);
		if (idx->z.last_off == 0) idx->z.offset0 = ret;
	}
	if (bin < 0) bin = hts_reg2bin(beg, end); // compute bin if this has not been done
	if (idx->z.last_bin != bin) { // then possibly write the binning index
		if (idx->z.save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
			insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, idx->z.last_off);
		if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // keep meta information
			idx->z.off_end = idx->z.last_off;
			insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.off_beg, idx->z.off_end);
			insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.n_mapped, idx->z.n_unmapped);
			idx->z.n_mapped = idx->z.n_unmapped = 0;
			idx->z.off_beg = idx->z.off_end;
		}
		idx->z.save_off = idx->z.last_off;
		idx->z.save_bin = idx->z.last_bin = bin;
		idx->z.save_tid = tid;
		if (tid < 0) {
			idx->z.last_off = offset;
			hts_idx_finish(idx);
			return 0;
		}
	}
	if (is_mapped) ++idx->z.n_mapped;
	else ++idx->z.n_unmapped;
	idx->z.last_off = offset;
	idx->z.last_coor = beg;
	return 0;
}

#endif
