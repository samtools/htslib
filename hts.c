#include <zlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "cram/cram.h"
#include "hfile.h"
#include "version.h"

#include "htslib/kseq.h"
#define KS_BGZF 1
#if KS_BGZF
    // pd3 todo: gzread() in BGZF
    KSTREAM_INIT2(, BGZF*, bgzf_read, 65536)
#else
    KSTREAM_INIT2(, gzFile, gzread, 16384)
#endif

#include "htslib/khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = 3;

const char *hts_version()
{
	return HTS_VERSION;
}

const unsigned char seq_nt16_table[256] = {
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

const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

/**********************
 *** Basic file I/O ***
 **********************/

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
	// Typically at most a couple of hundred bytes of input are required
	// to get a few bytes of output from inflate(), so hopefully this buffer
	// size suffices in general.
	unsigned char buffer[512];
	z_stream zs;
	ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

	if (npeek < 0) return 0;

	zs.zalloc = NULL;
	zs.zfree = NULL;
	zs.next_in = buffer;
	zs.avail_in = npeek;
	zs.next_out = dest;
	zs.avail_out = destsize;
	if (inflateInit2(&zs, 31) != Z_OK) return 0;

	while (zs.total_out < destsize)
		if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

	destsize = zs.total_out;
	inflateEnd(&zs);

	return destsize;
}

// Returns whether the block contains any control characters, i.e.,
// characters less than SPACE other than whitespace etc (ASCII BEL..CR).
static int is_binary(unsigned char *s, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++)
		if (s[i] < 0x07 || (s[i] >= 0x0e && s[i] < 0x20)) return 1;
	return 0;
}

htsFile *hts_open(const char *fn, const char *mode)
{
	htsFile *fp = NULL;
	hFILE *hfile = hopen(fn, mode);
	if (hfile == NULL) goto error;

	fp = (htsFile*)calloc(1, sizeof(htsFile));
	if (fp == NULL) goto error;

	fp->fn = strdup(fn);
	fp->is_be = ed_is_big();

	if (strchr(mode, 'r')) {
		unsigned char s[18];
		if (hpeek(hfile, s, 6) == 6 && memcmp(s, "CRAM", 4) == 0 &&
			s[4] >= 1 && s[4] <= 2 && s[5] <= 1) {
			fp->is_cram = 1;
		}
		else if (hpeek(hfile, s, 18) == 18 && s[0] == 0x1f && s[1] == 0x8b &&
				 (s[3] & 4) && memcmp(&s[12], "BC\2\0", 4) == 0) {
			// The stream is BGZF-compressed.  Decompress a few bytes to see
			// whether it's in a binary format (e.g., BAM or BCF, starting
			// with four bytes of magic including a control character) or is
			// a bgzipped SAM or VCF text file.
			fp->is_compressed = 1;
			if (is_binary(s, decompress_peek(hfile, s, 4))) fp->is_bin = 1;
			else fp->is_kstream = 1;
		}
		else if (hpeek(hfile, s, 2) == 2 && s[0] == 0x1f && s[1] == 0x8b) {
			// Plain GZIP header... so a gzipped text file.
			fp->is_compressed = 1;
			fp->is_kstream = 1;
		}
		else if (hpeek(hfile, s, 4) == 4 && is_binary(s, 4)) {
			// Binary format, but in a raw non-compressed form.
			fp->is_bin = 1;
		}
		else {
			fp->is_kstream = 1;
		}
	}
	else if (strchr(mode, 'w')) {
		fp->is_write = 1;
		if (strchr(mode, 'b')) fp->is_bin = 1;
		if (strchr(mode, 'c')) fp->is_cram = 1;
		if (strchr(mode, 'z')) fp->is_compressed = 1;
		else if (strchr(mode, 'u')) fp->is_compressed = 0;
		else fp->is_compressed = 2;    // not set, default behaviour
	}
	else goto error;

	if (fp->is_bin || (fp->is_write && fp->is_compressed==1)) {
		fp->fp.bgzf = bgzf_hopen(hfile, mode);
		if (fp->fp.bgzf == NULL) goto error;
	}
	else if (fp->is_cram) {
		fp->fp.cram = cram_dopen(hfile, fn, mode);
		if (fp->fp.cram == NULL) goto error;
	}
	else if (fp->is_kstream) {
	#if KS_BGZF
		BGZF *gzfp = bgzf_hopen(hfile, mode);
	#else
		// TODO Implement gzip hFILE adaptor
		hclose(hfile); // This won't work, especially for stdin
		gzFile gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	#endif
		if (gzfp) fp->fp.voidp = ks_init(gzfp);
		else goto error;
	}
	else {
		fp->fp.hfile = hfile;
	}

	return fp;

error:
	if (hts_verbose >= 2)
		fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

	if (hfile)
		hclose_abruptly(hfile);

	if (fp) {
		free(fp->fn);
		free(fp->fn_aux);
		free(fp);
	}
	return NULL;
}

int hts_close(htsFile *fp)
{
	int ret, save;

	if (fp->is_bin || (fp->is_write && fp->is_compressed==1)) {
		ret = bgzf_close(fp->fp.bgzf);
	} else if (fp->is_cram) {
		ret = cram_close(fp->fp.cram);
		if (!fp->is_write) {
			switch (cram_eof(fp->fp.cram)) {
			case 0:
				fprintf(stderr, "[E::%s] Failed to decode sequence.\n", __func__);
				return -1;
			case 2:
				fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
				break;
			default: /* case 1, expected EOF */
				break;
			}
		}
	} else if (fp->is_kstream) {
	#if KS_BGZF
		BGZF *gzfp = ((kstream_t*)fp->fp.voidp)->f;
		ret = bgzf_close(gzfp);
	#else
		gzFile gzfp = ((kstream_t*)fp->fp.voidp)->f;
		ret = gzclose(gzfp);
	#endif
		ks_destroy((kstream_t*)fp->fp.voidp);
	} else {
		ret = hclose(fp->fp.hfile);
	}

	save = errno;
	free(fp->fn);
	free(fp->fn_aux);
    free(fp->line.s);
	free(fp);
	errno = save;
	return ret;
}

int hts_set_fai_filename(htsFile *fp, const char *fn_aux)
{
	free(fp->fn_aux);
	if (fn_aux) {
		fp->fn_aux = strdup(fn_aux);
		if (fp->fn_aux == NULL) return -1;
	}
	else fp->fn_aux = NULL;

	return 0;
}

// For VCF/BCF backward sweeper. Not exposing these functions because their
// future is uncertain. Things will probably have to change with hFILE...
BGZF *hts_get_bgzfp(htsFile *fp)
{
    if ( fp->is_bin )
        return fp->fp.bgzf;
    else
        return ((kstream_t*)fp->fp.voidp)->f;
}
int hts_useek(htsFile *fp, long uoffset, int where)
{
    if ( fp->is_bin )
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
    {
        ks_rewind((kstream_t*)fp->fp.voidp);
        ((kstream_t*)fp->fp.voidp)->seek_pos = uoffset;
        return bgzf_useek(((kstream_t*)fp->fp.voidp)->f, uoffset, where);
    }
}
long hts_utell(htsFile *fp)
{
    if ( fp->is_bin )
        return bgzf_utell(fp->fp.bgzf);
    else
        return ((kstream_t*)fp->fp.voidp)->seek_pos;
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
	int ret, dret;
	ret = ks_getuntil((kstream_t*)fp->fp.voidp, delimiter, str, &dret);
	++fp->lineno;
	return ret;
}

char **hts_readlist(const char *string, int *n)
{
    if ( string[0]==':' ) 
        return hts_readlines(string+1,n);
    char *str = (char*) malloc(strlen(string)+2);
    if ( !str ) { *n = -1; return NULL; }
    str[0] = ':';
    strcpy(str+1, string);
    char **ret = hts_readlines(str,n);
    free(str);
    return ret;
}

char **hts_readlines(const char *fn, int *_n)
{
	int m = 0, n = 0, dret;
	char **s = 0;
#if KS_BGZF
	BGZF *fp = bgzf_open(fn, "r");
#else
	gzFile fp = gzopen(fn, "r");
#endif
	if ( fp ) { // read from file
		kstream_t *ks;
		kstring_t str;
		str.s = 0; str.l = str.m = 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			if (str.l == 0) continue;
			if (m == n) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(char*));
			}
			s[n++] = strdup(str.s);
		}
		ks_destroy(ks);
        #if KS_BGZF
            bgzf_close(fp);
        #else
		    gzclose(fp);
        #endif
		s = (char**)realloc(s, n * sizeof(char*));
		free(str.s);
	} else if (*fn == ':') { // read from string
		const char *q, *p;
		for (q = p = fn + 1;; ++p)
			if (*p == ',' || *p == 0) {
				if (m == n) {
					m = m? m<<1 : 16;
					s = (char**)realloc(s, m * sizeof(char*));
				}
				s[n] = (char*)calloc(p - q + 1, 1);
				strncpy(s[n++], q, p - q);
				q = p + 1;
				if (*p == 0) break;
			}
	} else return 0;
	s = (char**)realloc(s, n * sizeof(char*));
	*_n = n;
	return s;
}

int hts_file_type(const char *fname)
{
    int len = strlen(fname);
    if ( !strcasecmp(".vcf.gz",fname+len-7) ) return FT_VCF_GZ;
    if ( !strcasecmp(".vcf",fname+len-4) ) return FT_VCF;
    if ( !strcasecmp(".bcf",fname+len-4) ) return FT_BCF_GZ;
    if ( !strcmp("-",fname) ) return FT_STDIN;
    // ... etc

    int fd = open(fname, O_RDONLY);
    if ( !fd ) return 0;

    uint8_t magic[5];
    if ( read(fd,magic,2)!=2 ) { close(fd); return 0; }
    if ( !strncmp((char*)magic,"##",2) ) { close(fd); return FT_VCF; }
    if ( !strncmp((char*)magic,"BCF",3) ) { close(fd); return FT_BCF; }
    close(fd);

    if ( magic[0]==0x1f && magic[1]==0x8b ) // compressed
    {
        BGZF *fp = bgzf_open(fname, "r");
        if ( !fp ) return 0;
        if ( bgzf_read(fp, magic, 3)!=3 ) { bgzf_close(fp); return 0; }
        bgzf_close(fp);
        if ( !strncmp((char*)magic,"##",2) ) return FT_VCF;
        if ( !strncmp((char*)magic,"BCF",3) ) return FT_BCF_GZ;
    }
    return 0;
}

/****************
 *** Indexing ***
 ****************/

#define HTS_MIN_MARKER_DIST 0x10000

#define pair64_lt(a,b) ((a).u < (b).u)

#include "htslib/ksort.h"
KSORT_INIT(_off, hts_pair64_t, pair64_lt)

typedef struct {
	int32_t m, n;
	uint64_t loff;
	hts_pair64_t *list;
} bins_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
	int fmt, min_shift, n_lvls, n_bins;
	uint32_t l_meta;
	int32_t n, m;
	uint64_t n_no_coor;
	bidx_t **bidx;
	lidx_t *lidx;
	uint8_t *meta;
	struct {
		uint32_t last_bin, save_bin;
		int last_coor, last_tid, save_tid, finished;
		uint64_t last_off, save_off;
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

static inline void insert_to_l(lidx_t *l, int64_t _beg, int64_t _end, uint64_t offset, int min_shift)
{
	int i, beg, end;
	beg = _beg >> min_shift;
	end = (_end - 1) >> min_shift;
	if (l->m < end + 1) {
		int old_m = l->m;
		l->m = end + 1;
		kroundup32(l->m);
		l->offset = (uint64_t*)realloc(l->offset, l->m * 8);
		memset(l->offset + old_m, 0xff, 8 * (l->m - old_m)); // fill l->offset with (uint64_t)-1
	}
	if (beg == end) { // to save a loop in this case
		if (l->offset[beg] == (uint64_t)-1) l->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (l->offset[i] == (uint64_t)-1) l->offset[i] = offset;
	}
	if (l->n < end + 1) l->n = end + 1;
}

hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
{
	hts_idx_t *idx;
	idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
	idx->fmt = fmt;
	idx->min_shift = min_shift;
	idx->n_lvls = n_lvls;
	idx->n_bins = ((1<<(3 * n_lvls + 3)) - 1) / 7;
	idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
	idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
	idx->z.last_coor = 0xffffffffu;
	if (n) {
		idx->n = idx->m = n;
		idx->bidx = (bidx_t**)calloc(n, sizeof(bidx_t*));
		idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
	}
	return idx;
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
	bidx_t *bidx = idx->bidx[i];
	lidx_t *lidx = &idx->lidx[i];
	khint_t k;
	int l;
	uint64_t offset0 = 0;
	if (bidx) {
		k = kh_get(bin, bidx, idx->n_bins + 1);
		if (k != kh_end(bidx))
			offset0 = kh_val(bidx, k).list[0].u;
		for (l = 0; l < lidx->n && lidx->offset[l] == (uint64_t)-1; ++l)
			lidx->offset[l] = offset0;
	} else l = 1;
	for (; l < lidx->n; ++l) // fill missing values
		if (lidx->offset[l] == (uint64_t)-1)
			lidx->offset[l] = lidx->offset[l-1];
	if (bidx == 0) return;
	for (k = kh_begin(bidx); k != kh_end(bidx); ++k) // set loff
		if (kh_exist(bidx, k))
			kh_val(bidx, k).loff = kh_key(bidx, k) < idx->n_bins? lidx->offset[hts_bin_bot(kh_key(bidx, k), idx->n_lvls)] : 0;
	if (free_lidx) {
		free(lidx->offset);
		lidx->m = lidx->n = 0;
		lidx->offset = 0;
	}
}

static void compress_binning(hts_idx_t *idx, int i)
{
	bidx_t *bidx = idx->bidx[i];
	khint_t k;
	int l, m;
	if (bidx == 0) return;
	// merge a bin to its parent if the bin is too small
	for (l = idx->n_lvls; l > 0; --l) {
		unsigned start = hts_bin_first(l);
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
			bins_t *p, *q;
			if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins || kh_key(bidx, k) < start) continue;
			p = &kh_value(bidx, k);
			if (l < idx->n_lvls && p->n > 1) ks_introsort(_off, p->n, p->list);
			if ((p->list[p->n - 1].v>>16) - (p->list[0].u>>16) < HTS_MIN_MARKER_DIST) {
				khint_t kp;
				kp = kh_get(bin, bidx, hts_bin_parent(kh_key(bidx, k)));
				if (kp == kh_end(bidx)) continue;
				q = &kh_val(bidx, kp);
				if (q->n + p->n > q->m) {
					q->m = q->n + p->n;
					kroundup32(q->m);
					q->list = (hts_pair64_t*)realloc(q->list, q->m * 16);
				}
				memcpy(q->list + q->n, p->list, p->n * 16);
				q->n += p->n;
				free(p->list);
				kh_del(bin, bidx, k);
			}
		}
	}
	k = kh_get(bin, bidx, 0);
	if (k != kh_end(bidx)) ks_introsort(_off, kh_val(bidx, k).n, kh_val(bidx, k).list);
	// merge adjacent chunks that start from the same BGZF block
	for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
		bins_t *p;
		if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins) continue;
		p = &kh_value(bidx, k);
		for (l = 1, m = 0; l < p->n; ++l) {
			if (p->list[m].v>>16 >= p->list[l].u>>16) {
				if (p->list[m].v < p->list[l].v) p->list[m].v = p->list[l].v;
			} else p->list[++m] = p->list[l];
		}
		p->n = m + 1;
	}
}

void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
	int i;
	if (idx == NULL || idx->z.finished) return; // do not run this function on an empty index or multiple times
	if (idx->z.save_tid >= 0) {
		insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
		insert_to_b(idx->bidx[idx->z.save_tid], idx->n_bins + 1, idx->z.off_beg, final_offset);
		insert_to_b(idx->bidx[idx->z.save_tid], idx->n_bins + 1, idx->z.n_mapped, idx->z.n_unmapped);
	}
	for (i = 0; i < idx->n; ++i) {
		update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
		compress_binning(idx, i);
	}
	idx->z.finished = 1;
}

int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
{
	int bin;
	if (tid >= idx->m) { // enlarge the index
		int32_t oldm = idx->m;
		idx->m = idx->m? idx->m<<1 : 2;
		idx->bidx = (bidx_t**)realloc(idx->bidx, idx->m * sizeof(bidx_t*));
		idx->lidx = (lidx_t*) realloc(idx->lidx, idx->m * sizeof(lidx_t));
		memset(&idx->bidx[oldm], 0, (idx->m - oldm) * sizeof(bidx_t*));
		memset(&idx->lidx[oldm], 0, (idx->m - oldm) * sizeof(lidx_t));
	}
	if (idx->n < tid + 1) idx->n = tid + 1;
	if (tid < 0) ++idx->n_no_coor;
	if (idx->z.finished) return 0;
	if (idx->z.last_tid != tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
		idx->z.last_tid = tid;
		idx->z.last_bin = 0xffffffffu;
    /*
        // This seems a too prohibitive condition: When concatenating VCFs, the
        // sequences may come in a different order.

        } else if ((uint32_t)idx->z.last_tid > (uint32_t)tid) { // test if chromosomes are out of order
            if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted chromosomes\n", __func__);
            return -1;
    */
        if (idx->bidx[tid] != 0)
        {
            if (hts_verbose >= 1) fprintf(stderr, "[E::%s] chromosome blocks not continuous\n", __func__);
            return -1;
        }
	} else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted positions\n", __func__);
		return -1;
	}
	if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
	if (tid >= 0 && is_mapped)
		insert_to_l(&idx->lidx[tid], beg, end, idx->z.last_off, idx->min_shift); // last_off points to the start of the current record
	bin = hts_reg2bin(beg, end, idx->min_shift, idx->n_lvls);
	if ((int)idx->z.last_bin != bin) { // then possibly write the binning index
		if (idx->z.save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
			insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, idx->z.last_off);
		if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
			idx->z.off_end = idx->z.last_off;
			insert_to_b(idx->bidx[idx->z.save_tid], idx->n_bins + 1, idx->z.off_beg, idx->z.off_end);
			insert_to_b(idx->bidx[idx->z.save_tid], idx->n_bins + 1, idx->z.n_mapped, idx->z.n_unmapped);
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
	// For HTS_FMT_CRAI, idx actually points to a different type -- see sam.c
	if (idx->fmt == HTS_FMT_CRAI) { free(idx); return; }

	for (i = 0; i < idx->m; ++i) {
		bidx_t *bidx = idx->bidx[i];
		free(idx->lidx[i].offset);
		if (bidx == 0) continue;
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
			if (kh_exist(bidx, k))
				free(kh_value(bidx, k).list);
		kh_destroy(bin, bidx);
	}
	free(idx->bidx); free(idx->lidx); free(idx->meta);
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

static void hts_idx_save_core(const hts_idx_t *idx, void *fp, int fmt)
{
	int32_t i, size, is_be;
	int is_bgzf = (fmt != HTS_FMT_BAI);
	is_be = ed_is_big();
	if (is_be) {
		uint32_t x = idx->n;
		idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
	} else idx_write(is_bgzf, fp, &idx->n, 4);
	if (fmt == HTS_FMT_TBI && idx->l_meta) idx_write(is_bgzf, fp, idx->meta, idx->l_meta);
	for (i = 0; i < idx->n; ++i) {
		khint_t k;
		bidx_t *bidx = idx->bidx[i];
		lidx_t *lidx = &idx->lidx[i];
		// write binning index
		size = bidx? kh_size(bidx) : 0;
		if (is_be) { // big endian
			uint32_t x = size;
			idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
		} else idx_write(is_bgzf, fp, &size, 4);
		if (bidx == 0) goto write_lidx;
		for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
			bins_t *p;
			if (!kh_exist(bidx, k)) continue;
			p = &kh_value(bidx, k);
			if (is_be) { // big endian
				uint32_t x;
				x = kh_key(bidx, k); idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
				if (fmt == HTS_FMT_CSI) {
					uint64_t y = kh_val(bidx, k).loff;
					idx_write(is_bgzf, fp, ed_swap_4p(&y), 8);
				}
				x = p->n; idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
				swap_bins(p);
				idx_write(is_bgzf, fp, p->list, 16 * p->n);
				swap_bins(p);
			} else {
				idx_write(is_bgzf, fp, &kh_key(bidx, k), 4);
				if (fmt == HTS_FMT_CSI) idx_write(is_bgzf, fp, &kh_val(bidx, k).loff, 8);
				//int j;for(j=0;j<p->n;++j)fprintf(stderr,"%d,%llx,%d,%llx:%llx\n",kh_key(bidx,k),kh_val(bidx, k).loff,j,p->list[j].u,p->list[j].v);
				idx_write(is_bgzf, fp, &p->n, 4);
				idx_write(is_bgzf, fp, p->list, p->n << 4);
			}
		}
write_lidx:
		if (fmt != HTS_FMT_CSI) {
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
	if (is_be) { // write the number of reads without coordinates
		uint64_t x = idx->n_no_coor;
		idx_write(is_bgzf, fp, &x, 8);
	} else idx_write(is_bgzf, fp, &idx->n_no_coor, 8);
}

void hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
{
	char *fnidx;
	fnidx = (char*)calloc(1, strlen(fn) + 5);
	strcpy(fnidx, fn);
	if (fmt == HTS_FMT_CSI) {
		BGZF *fp;
		uint32_t x[3];
		int is_be, i;
		is_be = ed_is_big();
		fp = bgzf_open(strcat(fnidx, ".csi"), "w");
		bgzf_write(fp, "CSI\1", 4);
		x[0] = idx->min_shift; x[1] = idx->n_lvls; x[2] = idx->l_meta;
		if (is_be) {
			for (i = 0; i < 3; ++i)
				bgzf_write(fp, ed_swap_4p(&x[i]), 4);
		} else bgzf_write(fp, &x, 12);
		if (idx->l_meta) bgzf_write(fp, idx->meta, idx->l_meta);
		hts_idx_save_core(idx, fp, HTS_FMT_CSI);
		bgzf_close(fp);
	} else if (fmt == HTS_FMT_TBI) {
		BGZF *fp;
		fp = bgzf_open(strcat(fnidx, ".tbi"), "w");
		bgzf_write(fp, "TBI\1", 4);
		hts_idx_save_core(idx, fp, HTS_FMT_TBI);
		bgzf_close(fp);
	} else if (fmt == HTS_FMT_BAI) {
		FILE *fp;
		fp = fopen(strcat(fnidx, ".bai"), "w");
		fwrite("BAI\1", 1, 4, fp);
		hts_idx_save_core(idx, fp, HTS_FMT_BAI);
		fclose(fp);
	} else abort();
	free(fnidx);
}

static void hts_idx_load_core(hts_idx_t *idx, void *fp, int fmt)
{
	int32_t i, n, is_be;
	int is_bgzf = (fmt != HTS_FMT_BAI);
	is_be = ed_is_big();
	for (i = 0; i < idx->n; ++i) {
		bidx_t *h;
		lidx_t *l = &idx->lidx[i];
		uint32_t key;
		int j, absent;
		bins_t *p;
		h = idx->bidx[i] = kh_init(bin);
		idx_read(is_bgzf, fp, &n, 4);
		if (is_be) ed_swap_4p(&n);
		for (j = 0; j < n; ++j) {
			khint_t k;
			idx_read(is_bgzf, fp, &key, 4);
			if (is_be) ed_swap_4p(&key);
			k = kh_put(bin, h, key, &absent);
			p = &kh_val(h, k);
			if (fmt == HTS_FMT_CSI) {
				idx_read(is_bgzf, fp, &p->loff, 8);
				if (is_be) ed_swap_8p(&p->loff);
			} else p->loff = 0;
			idx_read(is_bgzf, fp, &p->n, 4);
			if (is_be) ed_swap_4p(&p->n);
			p->m = p->n;
			p->list = (hts_pair64_t*)malloc(p->m * 16);
			idx_read(is_bgzf, fp, p->list, p->n<<4);
			if (is_be) swap_bins(p);
		}
		if (fmt != HTS_FMT_CSI) { // load linear index
			int j;
			idx_read(is_bgzf, fp, &l->n, 4);
			if (is_be) ed_swap_4p(&l->n);
			l->m = l->n;
			l->offset = (uint64_t*)malloc(l->n << 3);
			idx_read(is_bgzf, fp, l->offset, l->n << 3);
			if (is_be) for (j = 0; j < l->n; ++j) ed_swap_8p(&l->offset[j]);
			for (j = 1; j < l->n; ++j) // fill missing values; may happen given older samtools and tabix
				if (l->offset[j] == 0) l->offset[j] = l->offset[j-1];
			update_loff(idx, i, 1);
		}
	}
	if (idx_read(is_bgzf, fp, &idx->n_no_coor, 8) != 8) idx->n_no_coor = 0;
	if (is_be) ed_swap_8p(&idx->n_no_coor);
}

hts_idx_t *hts_idx_load_local(const char *fn, int fmt)
{
	uint8_t magic[4];
	int i, is_be;
	hts_idx_t *idx;
	is_be = ed_is_big();
	if (fmt == HTS_FMT_CSI) {
		BGZF *fp;
		uint32_t x[3], n;
		uint8_t *meta = 0;
		if ((fp = bgzf_open(fn, "r")) == 0) return NULL;
		if ( bgzf_read(fp, magic, 4) != 4 ) return NULL;
		if ( bgzf_read(fp, x, 12) != 12 ) return NULL;
		if (is_be) for (i = 0; i < 3; ++i) ed_swap_4p(&x[i]);
		if (x[2]) {
			meta = (uint8_t*)malloc(x[2]);
			if ( bgzf_read(fp, meta, x[2]) != x[2] ) return NULL;
		}
		if ( bgzf_read(fp, &n, 4) != 4 ) return NULL;
		if (is_be) ed_swap_4p(&n);
		idx = hts_idx_init(n, fmt, 0, x[0], x[1]);
		idx->l_meta = x[2];
		idx->meta = meta;
		hts_idx_load_core(idx, fp, HTS_FMT_CSI);
		bgzf_close(fp);
	} else if (fmt == HTS_FMT_TBI) {
		BGZF *fp;
		uint32_t x[8];
		if ((fp = bgzf_open(fn, "r")) == 0) return NULL;
		if ( bgzf_read(fp, magic, 4) != 4 ) return NULL;
		if ( bgzf_read(fp, x, 32) != 32 ) return NULL;
		if (is_be) for (i = 0; i < 8; ++i) ed_swap_4p(&x[i]);
		idx = hts_idx_init(x[0], fmt, 0, 14, 5);
		idx->l_meta = 28 + x[7];
		idx->meta = (uint8_t*)malloc(idx->l_meta);
		memcpy(idx->meta, &x[1], 28);
		bgzf_read(fp, idx->meta + 28, x[7]);
		hts_idx_load_core(idx, fp, HTS_FMT_TBI);
		bgzf_close(fp);
	} else if (fmt == HTS_FMT_BAI) {
		uint32_t n;
		FILE *fp;
		if ((fp = fopen(fn, "rb")) == 0) return NULL;
		if ( fread(magic, 1, 4, fp) != 4 ) return NULL;
		if ( fread(&n, 4, 1, fp) != 1 ) return NULL;
		if (is_be) ed_swap_4p(&n);
		idx = hts_idx_init(n, fmt, 0, 14, 5);
		hts_idx_load_core(idx, fp, HTS_FMT_BAI);
		fclose(fp);
	} else abort();
	return idx;
}

void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)
{
	if (idx->meta) free(idx->meta);
	idx->l_meta = l_meta;
	if (is_copy) {
		idx->meta = (uint8_t*)malloc(l_meta);
		memcpy(idx->meta, meta, l_meta);
	} else idx->meta = meta;
}

uint8_t *hts_idx_get_meta(hts_idx_t *idx, int *l_meta)
{
	*l_meta = idx->l_meta;
	return idx->meta;
}

const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)
{
    if ( !idx->n )
    {
        *n = 0;
        return NULL;
    }

    int tid = 0, i;
    const char **names = (const char**) calloc(idx->n,sizeof(const char*));
    for (i=0; i<idx->n; i++)
    {
        bidx_t *bidx = idx->bidx[i];
        if ( !bidx ) continue;
        names[tid++] = getid(hdr,i);
    }
    *n = tid;
    return names;
}

/****************
 *** Iterator ***
 ****************/

static inline int reg2bins(int64_t beg, int64_t end, hts_itr_t *itr, int min_shift, int n_lvls)
{
	int l, t, s = min_shift + (n_lvls<<1) + n_lvls;
	if (beg >= end) return 0;
	if (end >= 1LL<<s) end = 1LL<<s;
	for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
		int b, e, n, i;
		b = t + (beg>>s); e = t + (end>>s); n = e - b + 1;
		if (itr->bins.n + n > itr->bins.m) {
			itr->bins.m = itr->bins.n + n;
			kroundup32(itr->bins.m);
			itr->bins.a = (int*)realloc(itr->bins.a, sizeof(int) * itr->bins.m);
		}
		for (i = b; i <= e; ++i) itr->bins.a[itr->bins.n++] = i;
	}
	return itr->bins.n;
}

hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
{
	int i, n_off, l, bin;
	hts_pair64_t *off;
	khint_t k;
	bidx_t *bidx;
	uint64_t min_off;
	hts_itr_t *iter = 0;

	if (tid < 0) {
		uint64_t off0 = (uint64_t)-1;
		khint_t k;
		if (tid == HTS_IDX_START) {
            if ( idx->n <=0 ) return 0;
            // Find the smallest offset, note that sequence ids may not be ordered sequentially
            for (i=0; i<idx->n; i++)
            {
                bidx = idx->bidx[i];
                k = kh_get(bin, bidx, idx->n_bins + 1);
                if (k == kh_end(bidx)) continue;
                if ( off0 > kh_val(bidx, k).list[0].u ) off0 = kh_val(bidx, k).list[0].u;
            }
		} else if (tid == HTS_IDX_NOCOOR) {
			if (idx->n > 0) {
				bidx = idx->bidx[idx->n - 1];
				k = kh_get(bin, bidx, idx->n_bins + 1);
				if (k == kh_end(bidx)) return 0;
				off0 = kh_val(bidx, k).list[0].v;
			} else return 0;
		} else off0 = 0;
		if (off0 != (uint64_t)-1) {
			iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
			iter->read_rest = 1;
			iter->curr_off = off0;
			iter->readrec = readrec;
			return iter;
		} else return 0;
	}
	if (beg < 0) beg = 0;
	if (end < beg) return 0;
	if ((bidx = idx->bidx[tid]) == 0) return 0;

	iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
	iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;
	iter->readrec = readrec;

	// compute min_off
	bin = hts_bin_first(idx->n_lvls) + (beg>>idx->min_shift);
	do {
		int first;
		k = kh_get(bin, bidx, bin);
		if (k != kh_end(bidx)) break;
		first = (hts_bin_parent(bin)<<3) + 1;
		if (bin > first) --bin;
		else bin = hts_bin_parent(bin);
	} while (bin);
	if (bin == 0) k = kh_get(bin, bidx, bin);
	min_off = k != kh_end(bidx)? kh_val(bidx, k).loff : 0;
	// retrieve bins
	reg2bins(beg, end, iter, idx->min_shift, idx->n_lvls);
	for (i = n_off = 0; i < iter->bins.n; ++i)
		if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx))
			n_off += kh_value(bidx, k).n;
	if (n_off == 0) return iter;
	off = (hts_pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < iter->bins.n; ++i) {
		if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx)) {
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

void hts_itr_destroy(hts_itr_t *iter)
{
	if (iter) { free(iter->off); free(iter->bins.a); free(iter); }
}

const char *hts_parse_reg(const char *s, int *beg, int *end)
{
	int i, k, l, name_end;
	*beg = *end = -1;
	name_end = l = strlen(s);
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
	}
	// parse the interval
	if (name_end < l) {
		char *tmp;
		tmp = (char*)alloca(l - name_end + 1);
		for (i = name_end + 1, k = 0; i < l; ++i)
			if (s[i] != ',') tmp[k++] = s[i];
		tmp[k] = 0;
		if ((*beg = strtol(tmp, &tmp, 10) - 1) < 0) *beg = 0;
		*end = *tmp? strtol(tmp + 1, &tmp, 10) : 1<<29;
		if (*beg > *end) name_end = l;
	}
	if (name_end == l) *beg = 0, *end = 1<<29;
	return s + name_end;
}

hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)
{
	int tid, beg, end;
	char *q, *tmp;
	if (strcmp(reg, ".") == 0)
		return itr_query(idx, HTS_IDX_START, 0, 1<<29, readrec);
	else if (strcmp(reg, "*") != 0) {
		q = (char*)hts_parse_reg(reg, &beg, &end);
		tmp = (char*)alloca(q - reg + 1);
		strncpy(tmp, reg, q - reg);
		tmp[q - reg] = 0;
		if ((tid = getid(hdr, tmp)) < 0)
			tid = getid(hdr, reg);
		if (tid < 0) return 0;
		return itr_query(idx, tid, beg, end, readrec);
	} else return itr_query(idx, HTS_IDX_NOCOOR, 0, 0, readrec);
}

int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
{
	int ret, tid, beg, end;
	if (iter && iter->finished) return -1;
	if (iter->read_rest) {
		if (iter->curr_off) { // seek to the start
			bgzf_seek(fp, iter->curr_off, SEEK_SET);
			iter->curr_off = 0; // only seek once
		}
		ret = iter->readrec(fp, data, r, &tid, &beg, &end);
		if (ret < 0) iter->finished = 1;
		return ret;
	}
	if (iter->off == 0) return -1;
	for (;;) {
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
			if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
				bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bgzf_tell(fp);
			}
			++iter->i;
		}
		if ((ret = iter->readrec(fp, data, r, &tid, &beg, &end)) >= 0) {
			iter->curr_off = bgzf_tell(fp);
			if (tid != iter->tid || beg >= iter->end) { // no need to proceed
				ret = -1; break;
			} else if (end > iter->beg && iter->end > beg) return ret;
		} else break; // end of file or error
	}
	iter->finished = 1;
	return ret;
}

/**********************
 *** Retrieve index ***
 **********************/

static char *test_and_fetch(const char *fn)
{
	FILE *fp;
	if (strstr(fn, "ftp://") == fn || strstr(fn, "http://") == fn) {
#ifdef _USE_KETFILE
		const int buf_size = 1 * 1024 * 1024;
		knetFile *fp_remote;
		uint8_t *buf;
		const char *p;
		for (p = fn + strlen(fn) - 1; p >= url; --p)
			if (*p == '/') break;
		++p; // p now points to the local file name
		if ((fp_remote = knet_open(fn, "r")) == 0) {
			if (hts_verbose >= 1) fprintf(stderr, "[E::%s] fail to open remote file\n", __func__);
			return 0;
		}
		if ((fp = fopen(fn, "w")) == 0) {
			if (hts_verbose >= 1) fprintf(stderr, "[E::%s] fail to create file in the working directory\n", __func__);
			knet_close(fp_remote);
			return 0;
		}
		if (hts_verbose >= 3) fprintf(stderr, "[M::%s] downloading file '%s' to local directory\n", __func__, fn);
		buf = (uint8_t*)calloc(buf_size, 1);
		while ((l = knet_read(fp_remote, buf, buf_size)) != 0) fwrite(buf, 1, l, fp);
		free(buf);
		fclose(fp);
		knet_close(fp_remote);
		return (char*)p;
#else
		return 0;
#endif
	} else {
		if ((fp = fopen(fn, "rb")) == 0) return 0;
		fclose(fp);
		return (char*)fn;
	}
}

char *hts_idx_getfn(const char *fn, const char *ext)
{
	int i, l_fn, l_ext;
	char *fnidx, *ret;
	l_fn = strlen(fn); l_ext = strlen(ext);
	fnidx = (char*)calloc(l_fn + l_ext + 1, 1);
	strcpy(fnidx, fn); strcpy(fnidx + l_fn, ext);
	if ((ret = test_and_fetch(fnidx)) == 0) {
		for (i = l_fn - 1; i > 0; --i)
			if (fnidx[i] == '.') break;
		strcpy(fnidx + i, ext);
		ret = test_and_fetch(fnidx);
	}
	if (ret == 0) {
		free(fnidx);
		return 0;
	}
	l_fn = strlen(ret);
	memmove(fnidx, ret, l_fn + 1);
	return fnidx;
}

hts_idx_t *hts_idx_load(const char *fn, int fmt)
{
	char *fnidx;
	hts_idx_t *idx;
	fnidx = hts_idx_getfn(fn, ".csi");
	if (fnidx) fmt = HTS_FMT_CSI;
	else fnidx = hts_idx_getfn(fn, fmt == HTS_FMT_BAI? ".bai" : ".tbi");
	if (fnidx == 0) return 0;
	idx = hts_idx_load_local(fnidx, fmt);
	free(fnidx);
	return idx;
}
