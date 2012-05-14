#ifndef HTS_NO_INDEX

#include <stdio.h>
#include "endian.h"
#include "bgzf.h"
#include "hts.h"

#define IDX_LIDX_SHIFT    14
#define IDX_MAX_BIN       37450

#define pair64_lt(a,b) ((a).u < (b).u)

#include "ksort.h"
KSORT_INIT(_off, hts_pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	hts_pair64_t *list;
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
	uint64_t n_no_coor;
	bidx_t **bidx;
	lidx_t *lidx;
	struct {
		uint32_t last_bin, save_bin;
		int last_coor, last_tid, save_tid, finished;
		uint64_t last_off, save_off, offset0;
		uint64_t off_beg, off_end;
		uint64_t n_mapped, n_unmapped;
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
		l->list = (hts_pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (hts_pair64_t*)realloc(l->list, l->m * 16);
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

hts_idx_t *hts_idx_init(int n, uint64_t offset0)
{
	hts_idx_t *idx;
	idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
	idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
	idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
	idx->z.last_coor = 0xffffffffu;
	idx->z.offset0 = (uint64_t)-1;
	if (n) {
		idx->n = idx->m = n;
		idx->bidx = (bidx_t**)calloc(n, sizeof(void*));
		idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
	}
	return idx;
}

void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
	int i;
	if (idx->z.finished) return; // do not run this function multiple times
	if (idx->z.save_tid >= 0) {
		insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
		insert_to_b(idx->bidx[idx->z.save_tid], IDX_MAX_BIN, idx->z.off_beg, final_offset);
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
	if (tid >= idx->m) { // enlarge the index
		int32_t oldm = idx->m;
		idx->m = idx->m? idx->m<<1 : 2;
		idx->bidx = (bidx_t**)realloc(idx->bidx, sizeof(void*));
		idx->lidx = (lidx_t*) realloc(idx->lidx, sizeof(lidx_t));
		memset(idx->bidx[oldm],  0, (idx->m - oldm) * sizeof(void*));
		memset(&idx->lidx[oldm], 0, (idx->m - oldm) * sizeof(lidx_t));
	}
	if (tid < 0) ++idx->n_no_coor;
	if (idx->z.finished) return 0;
	if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
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
			hts_idx_finish(idx, offset);
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
	for (i = 0; i < idx->m; ++i) {
		bidx_t *bidx = idx->bidx[i];
		free(idx->lidx[i].offset);
		if (bidx == 0) continue;
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
			if (kh_exist(bidx, k))
				free(kh_value(bidx, k).list);
		kh_destroy(bin, bidx);
	}
	free(idx->bidx); free(idx->lidx);
	free(idx);
}

static inline long idx_read(int is_bgzf, void *fp, void *buf, long l)
{
	if (is_bgzf) return bgzf_read((BGZF*)fp, buf, l);
	else return (long)fread(buf, 1, l, (FILE*)fp);
}

static inline long idx_write(int is_bgzf, void *fp, const void *buf, long l)
{
	if (is_bgzf) return bgzf_write((BGZF*)fp, buf, l);
	else return (long)fwrite(buf, 1, l, (FILE*)fp);
}

static inline void swap_bins(bins_t *p)
{
	int i;
	for (i = 0; i < p->n; ++i) {
		ed_swap_8p(&p->list[i].u);
		ed_swap_8p(&p->list[i].v);
	}
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
					swap_bins(p);
					idx_write(is_bgzf, fp, p->list, 16 * p->n);
					swap_bins(p);
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
	if (is_be) { // write the number of reads without coordinates
		uint64_t x = idx->n_no_coor;
		idx_write(is_bgzf, fp, &x, 8);
	} else idx_write(is_bgzf, fp, &idx->n_no_coor, 8);
}

hts_idx_t *hts_idx_load(void *fp, int is_bgzf)
{
	int32_t i, n, is_be;
	hts_idx_t *idx;

	is_be = ed_is_big();
	idx_read(is_bgzf, fp, &n, 4);
	if (is_be) ed_swap_4p(&n);
	idx = hts_idx_init(n, 0);
	for (i = 0; i < idx->n; ++i) {
		bidx_t *h;
		lidx_t *l = &idx->lidx[i];
		uint32_t key;
		int j, absent;
		bins_t *p;
		h = idx->bidx[i] = kh_init(bin);
		// load binning index
		idx_read(is_bgzf, fp, &n, 4);
		if (is_be) ed_swap_4p(&n);
		for (j = 0; j < n; ++j) {
			khint_t k;
			idx_read(is_bgzf, fp, &key, 4);
			if (is_be) ed_swap_4p(&key);
			k = kh_put(bin, h, key, &absent);
			p = &kh_val(h, k);
			idx_read(is_bgzf, fp, &p->n, 4);
			if (is_be) ed_swap_4p(&p->n);
			p->m = p->n;
			p->list = (hts_pair64_t*)malloc(p->m * 16);
			idx_read(is_bgzf, fp, p->list, p->n<<4);
			if (is_be) swap_bins(p);
		}
		// load linear index
		idx_read(is_bgzf, fp, &l->n, 4);
		if (is_be) ed_swap_4p(&l->n);
		l->m = l->n;
		l->offset = (uint64_t*)malloc(l->n << 3);
		idx_read(is_bgzf, fp, l->offset, l->n << 3);
		if (is_be) for (j = 0; j < l->n; ++j) ed_swap_8p(&l->offset[j]);
	}
	if (idx_read(is_bgzf, fp, &idx->n_no_coor, 8) != 8) idx->n_no_coor = 0;
	if (is_be) ed_swap_8p(&idx->n_no_coor);
	return idx;
}


static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[IDX_MAX_BIN])
{
	int i = 0, k;
	if (beg >= end) return 0;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}

hts_iter_t *hts_iter_query(const hts_idx_t *idx, int tid, int beg, int end)
{
	uint16_t *bins;
	int i, n_bins, n_off, l;
	hts_pair64_t *off;
	khint_t k;
	bidx_t *bidx;
	uint64_t min_off;
	hts_iter_t *iter = 0;

	if (beg < 0) beg = 0;
	if (end < beg) return 0;
	iter = (hts_iter_t*)calloc(1, sizeof(hts_iter_t));
	iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;

	if ((bidx = idx->bidx[tid]) == 0) return 0;
	bins = (uint16_t*)alloca(IDX_MAX_BIN<<1);
	n_bins = reg2bins(beg, end, bins);
	if (idx->lidx[tid].n > 0) {
		min_off = (beg>>IDX_LIDX_SHIFT >= idx->lidx[tid].n)? idx->lidx[tid].offset[idx->lidx[tid].n-1]
			: idx->lidx[tid].offset[beg>>IDX_LIDX_SHIFT];
		if (min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
			int n = beg>>IDX_LIDX_SHIFT;
			if (n > idx->lidx[tid].n) n = idx->lidx[tid].n;
			for (i = n - 1; i >= 0; --i)
				if (idx->lidx[tid].offset[i] != 0) break;
			if (i >= 0) min_off = idx->lidx[tid].offset[i];
		}
	} else min_off = 0; // tabix 0.1.2 may produce such index files
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(bin, bidx, bins[i])) != kh_end(bidx))
			n_off += kh_value(bidx, k).n;
	}
	if (n_off == 0) return iter;
	off = (hts_pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(bin, bidx, bins[i])) != kh_end(bidx)) {
			int j;
			bins_t *p = &kh_value(bidx, k);
			for (j = 0; j < p->n; ++j)
				if (p->list[j].v > min_off) off[n_off++] = p->list[j];
		}
	}
	if (n_off == 0) {
		free(off); return iter;
	}
	ks_introsort(_off, n_off, off);
	// resolve completely contained adjacent blocks
	for (i = 1, l = 0; i < n_off; ++i)
		if (off[l].v < off[i].v) off[++l] = off[i];
	n_off = l + 1;
	// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
	for (i = 1; i < n_off; ++i)
		if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
	// merge adjacent blocks
	for (i = 1, l = 0; i < n_off; ++i) {
		if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
		else off[++l] = off[i];
	}
	n_off = l + 1;
	iter->n_off = n_off; iter->off = off;
	return iter;
}

void hts_iter_destroy(hts_iter_t *iter)
{
	if (iter) { free(iter->off); free(iter); }
}

#endif // ~!defined(HTS_NO_INDEX)
