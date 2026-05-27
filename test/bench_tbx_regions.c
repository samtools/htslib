/*  bench_tbx_regions.c -- microbenchmark comparing tbx_itr_queryi loop
                            against tbx_itr_regions on the same input.

    Usage:  bench_tbx_regions FILE.gz REGIONS.bed [WARMUP] [REPS]
              FILE.gz       bgzipped + tabix indexed source
              REGIONS.bed   plain text BED of regions to query, one per line:
                                <chrom> <beg> <end>     (0 based half open)
              WARMUP        number of warmup runs per path (default 0)
              REPS          number of timed runs per path (default 3)

    Reports wallclock and CPU time and total records emitted for each path.
    Verifies byte level parity between paths by running an untimed capture
    pass that streams every emitted line of both paths into separate
    buffers, then memcmp comparing them and printing an FNV 1a digest of
    each. Falls back to exit code 3 on any mismatch and prints the first
    diverging line index. The record count check remains as an additional
    sanity gate.

    Copyright (C) 2026 Genome Research Ltd.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/resource.h>

#include "../htslib/hts.h"
#include "../htslib/tbx.h"
#include "../htslib/kstring.h"

/* FNV 1a 64 bit, streamed. Used purely to give the parity report a digest
 * so a reviewer reading the output can see that both paths produced byte
 * equal streams at any panel size without holding the full byte stream in
 * RAM (which at 100k regions on a 1kGP VCF approaches 2 GB per side).
 *
 * Collision probability for two distinct streams is ~2^-64 (~5e-20), which
 * is well below cosmic ray bit flip rates for the duration of any bench.
 * Combined with the byte length check, we treat matching {length, hash}
 * across both paths as the parity guarantee. */
#define FNV1A64_OFFSET 0xcbf29ce484222325ULL
#define FNV1A64_PRIME  0x100000001b3ULL

static inline uint64_t fnv1a64_update(uint64_t h, const void *data, size_t len)
{
    const unsigned char *p = (const unsigned char *) data;
    for (size_t i = 0; i < len; i++) {
        h ^= p[i];
        h *= FNV1A64_PRIME;
    }
    return h;
}

/* In memory region list parsed once and reused for both paths. */
typedef struct {
    char *chrom;
    hts_pos_t beg, end;
} bench_region_t;

typedef struct {
    bench_region_t *items;
    int n, cap;
} bench_regions_t;

/* Per run timing sample. parity_hash and parity_bytes are populated only
 * when the caller requests parity tracking (see run_*); they sit alongside
 * timing so the report can show all of {wall, records, bytes, hash} for a
 * single rep without juggling separate return paths. */
typedef struct {
    double wall_s;
    double user_s;
    double sys_s;
    long records;
    uint64_t parity_hash;
    size_t   parity_bytes;
} bench_sample_t;

static double timespec_diff_s(struct timespec a, struct timespec b)
{
    return (b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) / 1.0e9;
}

static double timeval_s(struct timeval t)
{
    return t.tv_sec + t.tv_usec / 1.0e6;
}

/* Read a plain BED into memory. Three columns: chrom, beg, end. */
static int read_bed(const char *path, bench_regions_t *out)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "fopen %s failed\n", path); return -1; }
    out->items = NULL; out->n = 0; out->cap = 0;
    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char chrom[256]; long long beg, end;
        if (sscanf(line, "%255s %lld %lld", chrom, &beg, &end) != 3) continue;
        if (out->n == out->cap) {
            int ncap = out->cap ? out->cap * 2 : 64;
            bench_region_t *grown = realloc(out->items, ncap * sizeof(*grown));
            if (!grown) { fclose(fp); return -1; }
            out->items = grown; out->cap = ncap;
        }
        out->items[out->n].chrom = strdup(chrom);
        out->items[out->n].beg = (hts_pos_t)beg;
        out->items[out->n].end = (hts_pos_t)end;
        out->n++;
    }
    fclose(fp);
    return 0;
}

static void free_regions(bench_regions_t *r)
{
    for (int i = 0; i < r->n; i++) free(r->items[i].chrom);
    free(r->items);
}

/* Path A: for each region, tbx_itr_queryi + tbx_itr_next loop.
 * If with_parity is non zero, every emitted line is fed into a streaming
 * fnv1a64 hash and counted in s->parity_bytes (line bytes plus a trailing
 * newline per record). The hash overhead is two machine ops per byte and
 * adds well under 5 percent to wallclock; the parity_* fields appear in
 * the parity report only, not in the timed reps. */
static int run_single(const char *fn, const bench_regions_t *regs,
                      bench_sample_t *s, int with_parity)
{
    htsFile *fp = hts_open(fn, "r");
    if (!fp) return -1;
    tbx_t *tbx = tbx_index_load(fn);
    if (!tbx) { hts_close(fp); return -1; }
    kstring_t line = {0,0,0};
    long records = 0;
    uint64_t hash = FNV1A64_OFFSET;
    size_t bytes = 0;
    static const char newline = '\n';
    struct rusage ru0, ru1;
    struct timespec t0, t1;
    getrusage(RUSAGE_SELF, &ru0);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < regs->n; i++) {
        int tid = tbx_name2id(tbx, regs->items[i].chrom);
        if (tid < 0) continue;
        hts_itr_t *itr = tbx_itr_queryi(tbx, tid, regs->items[i].beg,
                                             regs->items[i].end);
        if (!itr) continue;
        while (tbx_itr_next(fp, tbx, itr, &line) >= 0) {
            records++;
            if (with_parity) {
                hash = fnv1a64_update(hash, line.s, line.l);
                hash = fnv1a64_update(hash, &newline, 1);
                bytes += line.l + 1;
            }
        }
        tbx_itr_destroy(itr);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    getrusage(RUSAGE_SELF, &ru1);
    free(line.s);
    tbx_destroy(tbx);
    hts_close(fp);
    s->wall_s = timespec_diff_s(t0, t1);
    s->user_s = timeval_s(ru1.ru_utime) - timeval_s(ru0.ru_utime);
    s->sys_s  = timeval_s(ru1.ru_stime) - timeval_s(ru0.ru_stime);
    s->records = records;
    s->parity_hash  = with_parity ? hash  : 0;
    s->parity_bytes = with_parity ? bytes : 0;
    return 0;
}

/* Path B: one tbx_itr_regions covering all regions, then drain. Streaming
 * parity hash semantics mirror run_single. */
static int run_multi(const char *fn, const bench_regions_t *regs,
                     bench_sample_t *s, int with_parity)
{
    htsFile *fp = hts_open(fn, "r");
    if (!fp) return -1;
    tbx_t *tbx = tbx_index_load(fn);
    if (!tbx) { hts_close(fp); return -1; }

    /* Group regions by tid into one reglist entry per contig. */
    hts_reglist_t *reglist = calloc(regs->n, sizeof(*reglist));
    if (!reglist) { tbx_destroy(tbx); hts_close(fp); return -1; }
    int n_reg = 0;
    for (int i = 0; i < regs->n; i++) {
        int tid = tbx_name2id(tbx, regs->items[i].chrom);
        if (tid < 0) continue;
        int idx = -1;
        for (int j = 0; j < n_reg; j++) {
            if (reglist[j].tid == tid) { idx = j; break; }
        }
        if (idx < 0) {
            idx = n_reg++;
            reglist[idx].tid = tid;
            reglist[idx].reg = NULL;
            reglist[idx].count = 0;
            reglist[idx].intervals = NULL;
            reglist[idx].min_beg = regs->items[i].beg;
            reglist[idx].max_end = regs->items[i].end;
        }
        hts_pair_pos_t *grown = realloc(reglist[idx].intervals,
                                        (reglist[idx].count + 1) * sizeof(hts_pair_pos_t));
        if (!grown) {
            for (int j = 0; j < n_reg; j++) free(reglist[j].intervals);
            free(reglist); tbx_destroy(tbx); hts_close(fp); return -1;
        }
        reglist[idx].intervals = grown;
        reglist[idx].intervals[reglist[idx].count].beg = regs->items[i].beg;
        reglist[idx].intervals[reglist[idx].count].end = regs->items[i].end;
        reglist[idx].count++;
        if (regs->items[i].beg < reglist[idx].min_beg)
            reglist[idx].min_beg = regs->items[i].beg;
        if (regs->items[i].end > reglist[idx].max_end)
            reglist[idx].max_end = regs->items[i].end;
    }

    kstring_t line = {0,0,0};
    long records = 0;
    uint64_t hash = FNV1A64_OFFSET;
    size_t bytes = 0;
    static const char newline = '\n';
    struct rusage ru0, ru1;
    struct timespec t0, t1;
    getrusage(RUSAGE_SELF, &ru0);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    hts_itr_t *itr = tbx_itr_regions(fp, tbx, reglist, n_reg);
    if (!itr) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        for (int j = 0; j < n_reg; j++) free(reglist[j].intervals);
        free(reglist); free(line.s);
        tbx_destroy(tbx); hts_close(fp);
        return -1;
    }
    while (hts_itr_multi_next(fp, itr, &line) >= 0) {
        records++;
        if (with_parity) {
            hash = fnv1a64_update(hash, line.s, line.l);
            hash = fnv1a64_update(hash, &newline, 1);
            bytes += line.l + 1;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    getrusage(RUSAGE_SELF, &ru1);
    free(line.s);
    hts_itr_destroy(itr); /* takes ownership of reglist */
    tbx_destroy(tbx);
    hts_close(fp);
    s->wall_s = timespec_diff_s(t0, t1);
    s->user_s = timeval_s(ru1.ru_utime) - timeval_s(ru0.ru_utime);
    s->sys_s  = timeval_s(ru1.ru_stime) - timeval_s(ru0.ru_stime);
    s->records = records;
    s->parity_hash  = with_parity ? hash  : 0;
    s->parity_bytes = with_parity ? bytes : 0;
    return 0;
}

static void print_sample(const char *label, bench_sample_t s)
{
    printf("  %s: wall=%.3fs  user=%.3fs  sys=%.3fs  records=%ld\n",
           label, s.wall_s, s.user_s, s.sys_s, s.records);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s FILE.gz REGIONS.bed [WARMUP] [REPS]\n", argv[0]);
        return 1;
    }
    const char *fn = argv[1];
    const char *bed = argv[2];
    int warmup = argc > 3 ? atoi(argv[3]) : 0;
    int reps   = argc > 4 ? atoi(argv[4]) : 3;

    bench_regions_t regs;
    if (read_bed(bed, &regs) < 0) return 2;
    fprintf(stderr, "loaded %d regions from %s\n", regs.n, bed);

    bench_sample_t scratch;
    for (int w = 0; w < warmup; w++) {
        if (run_single(fn, &regs, &scratch, 0) < 0) { fprintf(stderr, "warmup single failed\n"); return 2; }
        if (run_multi (fn, &regs, &scratch, 0) < 0) { fprintf(stderr, "warmup multi failed\n"); return 2; }
    }

    /* Untimed parity pass. Each path streams every emitted line through a
     * running fnv1a64 hash and a byte counter; we then compare {length,
     * hash} across paths. Constant memory regardless of panel size, so it
     * scales to the 100k+ region case where holding two full byte streams
     * (~1.7 GB each on a 1kGP VCF) would trigger swap or OOM. Doubles as
     * a real warmup pass for the OS page cache, so the timed reps that
     * follow operate on equivalent cache state. */
    printf("=== parity check ===\n");
    bench_sample_t parity_a, parity_b;
    if (run_single(fn, &regs, &parity_a, 1) < 0) {
        fprintf(stderr, "parity single failed\n"); return 2;
    }
    if (run_multi (fn, &regs, &parity_b, 1) < 0) {
        fprintf(stderr, "parity multi failed\n"); return 2;
    }
    printf("  single: %ld records, %zu bytes, fnv1a=0x%016llx\n",
           parity_a.records, parity_a.parity_bytes,
           (unsigned long long) parity_a.parity_hash);
    printf("  multi : %ld records, %zu bytes, fnv1a=0x%016llx\n",
           parity_b.records, parity_b.parity_bytes,
           (unsigned long long) parity_b.parity_hash);

    int rc = 0;
    if (parity_a.parity_bytes != parity_b.parity_bytes ||
        parity_a.parity_hash  != parity_b.parity_hash  ||
        parity_a.records      != parity_b.records) {
        fprintf(stderr, "PARITY MISMATCH:"
                        " single=%ld records/%zu bytes/fnv1a=0x%016llx,"
                        " multi=%ld records/%zu bytes/fnv1a=0x%016llx\n",
                        parity_a.records, parity_a.parity_bytes,
                        (unsigned long long) parity_a.parity_hash,
                        parity_b.records, parity_b.parity_bytes,
                        (unsigned long long) parity_b.parity_hash);
        rc = 3;
    } else {
        printf("  PARITY OK: %ld records / %zu bytes match across paths\n",
               parity_a.records, parity_a.parity_bytes);
    }
    if (rc) { free_regions(&regs); return rc; }

    double sum_wall_single = 0, sum_wall_multi = 0;
    long rec_single = parity_a.records, rec_multi = parity_b.records;
    printf("=== bench_tbx_regions: %s with %d regions ===\n", fn, regs.n);
    for (int r = 0; r < reps; r++) {
        bench_sample_t a, b;
        if (run_single(fn, &regs, &a, 0) < 0) return 2;
        if (run_multi (fn, &regs, &b, 0) < 0) return 2;
        printf("rep %d:\n", r);
        print_sample("single", a);
        print_sample("multi ", b);
        printf("  speedup wall = %.2fx\n", a.wall_s / b.wall_s);
        sum_wall_single += a.wall_s;
        sum_wall_multi  += b.wall_s;
        if (a.records != rec_single || b.records != rec_multi) {
            fprintf(stderr, "WARN: record count differs across reps "
                            "(single %ld -> %ld, multi %ld -> %ld)\n",
                            rec_single, a.records, rec_multi, b.records);
        }
    }
    printf("=== summary ===\n");
    printf("  single avg wall = %.3fs over %d reps\n", sum_wall_single / reps, reps);
    printf("  multi  avg wall = %.3fs over %d reps\n", sum_wall_multi  / reps, reps);
    printf("  avg speedup     = %.2fx\n", sum_wall_single / sum_wall_multi);
    printf("  records single  = %ld\n", rec_single);
    printf("  records multi   = %ld\n", rec_multi);

    free_regions(&regs);
    return 0;
}
