#ifndef HTS_NO_INDEX

#include <stdio.h>
#include "endian.h"
#include "bgzf.h"
#include "hts.h"

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

struct __hts_idx_t {
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

hts_idx_t *hts_idx_init(uint64_t offset0)
{
	hts_idx_t *idx;
	idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
	idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
	idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
	idx->z.last_coor = 0xffffffffu;
	idx->z.offset0 = (uint64_t)-1;
	return idx;
}

void hts_idx_finish(hts_idx_t *idx)
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

int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int bin, int is_mapped)
{
	if (tid < 0) ++idx->z.n_no_coor;
	if (idx->z.finished) return 0;
	if (idx->z.last_tid < tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
		idx->z.last_tid = tid;
		idx->z.last_bin = 0xffffffffu;
	} else if ((uint32_t)idx->z.last_tid > (uint32_t)tid) { // test if chromosomes are out of order
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted chromosomes\n", __func__);
		return -1;
	} else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted positions\n", __func__);
		return -1;
	}
	if (tid >= 0 && is_mapped) {
		uint64_t ret;
		ret = insert_to_l(&idx->lidx[tid], beg, end, idx->z.last_off); // last_off points to the start of the current record
		if (idx->z.last_off == 0) idx->z.offset0 = ret; // I forgot the purpose of offset0
	}
	if (bin < 0) bin = hts_reg2bin(beg, end); // compute bin if this has not been done
	if (idx->z.last_bin != bin) { // then possibly write the binning index
		if (idx->z.save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
			insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, idx->z.last_off);
		if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
			idx->z.off_end = idx->z.last_off;
			insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.off_beg, idx->z.off_end);
			insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.n_mapped, idx->z.n_unmapped);
			idx->z.n_mapped = idx->z.n_unmapped = 0;
			idx->z.off_beg = idx->z.off_end;
		}
		idx->z.save_off = idx->z.last_off;
		idx->z.save_bin = idx->z.last_bin = bin;
		idx->z.save_tid = tid;
		if (tid < 0) { // come to the end of the records having coordinates
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

void hts_idx_destroy(hts_idx_t *idx)
{
	khint_t k;
	int i;
	if (idx == 0) return;
	for (i = 0; i < idx->n; ++i) {
		bidx_t *bidx = idx->bidx[i];
		lidx_t *lidx = &idx->lidx[i];
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
			if (kh_exist(bidx, k))
				free(kh_value(bidx, k).list);
		kh_destroy(bin, bidx);
		free(lidx->offset);
	}
	free(idx->bidx); free(idx->lidx);
	free(idx);
}

static inline long idx_write(int is_bgzf, void *fp, const void *buf, long l)
{
	if (is_bgzf) return bgzf_write((BGZF*)fp, buf, l);
	else return (long)fwrite(buf, 1, l, (FILE*)fp);
}

void hts_idx_save(const hts_idx_t *idx, void *fp, int is_bgzf)
{
	int32_t i, size, is_be;
	is_be = ed_is_big();
	if (is_be) {
		uint32_t x = idx->n;
		idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
	} else idx_write(is_bgzf, fp, &idx->n, 4);
	for (i = 0; i < idx->n; ++i) {
		khint_t k;
		bidx_t *bidx = idx->bidx[i];
		lidx_t *lidx = &idx->lidx[i];
		// write binning index
		size = kh_size(bidx);
		if (is_be) { // big endian
			uint32_t x = size;
			idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
		} else idx_write(is_bgzf, fp, &size, 4);
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
			if (kh_exist(bidx, k)) {
				bins_t *p = &kh_value(bidx, k);
				if (is_be) { // big endian
					uint32_t x;
					x = kh_key(bidx, k); idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
					x = p->n; idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
					for (x = 0; (int)x < p->n; ++x) {
						ed_swap_8p(&p->list[x].u);
						ed_swap_8p(&p->list[x].v);
					}
					idx_write(is_bgzf, fp, p->list, 16 * p->n);
					for (x = 0; (int)x < p->n; ++x) {
						ed_swap_8p(&p->list[x].u);
						ed_swap_8p(&p->list[x].v);
					}
				} else {
					idx_write(is_bgzf, fp, &kh_key(bidx, k), 4);
					idx_write(is_bgzf, fp, &p->n, 4);
					idx_write(is_bgzf, fp, p->list, p->n << 4);
				}
			}
		}
		// write linear index
		if (is_be) {
			int32_t x = lidx->n;
			idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
			for (x = 0; x < lidx->n; ++x) ed_swap_8p(&lidx->offset[x]);
			idx_write(is_bgzf, fp, lidx->offset, lidx->n << 3);
			for (x = 0; x < lidx->n; ++x) ed_swap_8p(&lidx->offset[x]);
		} else {
			idx_write(is_bgzf, fp, &lidx->n, 4);
			idx_write(is_bgzf, fp, lidx->offset, lidx->n << 3);
		}
	}
}

#endif // ~!defined(HTS_NO_INDEX)
