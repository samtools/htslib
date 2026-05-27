/*  test_tbx_multi.c -- correctness test for tbx_itr_regions.

    Compares the records returned by the multi-region path (tbx_itr_regions
    plus hts_itr_multi_next) against the records returned by walking the
    equivalent regions one at a time through the single-region path
    (tbx_itr_queryi plus tbx_itr_next plus tbx_readrec).

    If both paths return byte-identical line sets, the new function is
    correct by construction relative to the well-trodden single-region
    callback. Otherwise the test reports the difference and exits non-zero.

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

#include "../htslib/hts.h"
#include "../htslib/tbx.h"
#include "../htslib/kstring.h"

/* Regions queried by both paths. Designed to exercise the multi-region code
 * path under several shapes: point queries, wide range hitting multiple
 * records, empty region returning nothing, and cross-contig spread.
 * The fixture contains:
 *   1\t99\t100\tregion_a    (1-based pos 100)
 *   1\t149\t150\tregion_b   (pos 150)
 *   1\t199\t200\tregion_c   (pos 200)
 *   1\t299\t300\tregion_d   (pos 300)
 *   2\t99\t100\tregion_e
 *   2\t199\t200\tregion_f
 */
static const struct { const char *chrom; hts_pos_t beg, end; } QUERIES[] = {
    { "1",   99, 100 },   /* point query, hits region_a */
    { "1",  140, 210 },   /* wide range, hits region_b and region_c */
    { "1", 5000, 6000 },  /* empty region, no records */
    { "1",  290, 310 },   /* point-ish query, hits region_d */
    { "2",    0, 1000 },  /* spans full contig 2 range, hits region_e and region_f */
};
static const int N_QUERIES = sizeof(QUERIES) / sizeof(QUERIES[0]);

/* Append a line to a kstring. */
static void append_line(kstring_t *out, const char *line)
{
    kputs(line, out);
    kputc('\n', out);
}

/* Single-region path: one tbx_itr_queryi per query, append results. */
static int run_single_region(const char *filename, kstring_t *out)
{
    htsFile *fp = hts_open(filename, "r");
    if (!fp) { fprintf(stderr, "single: hts_open failed: %s\n", filename); return -1; }
    tbx_t *tbx = tbx_index_load(filename);
    if (!tbx) { fprintf(stderr, "single: tbx_index_load failed\n"); hts_close(fp); return -1; }

    kstring_t s = {0, 0, 0};
    int rc = 0;
    for (int i = 0; i < N_QUERIES; i++) {
        int tid = tbx_name2id(tbx, QUERIES[i].chrom);
        if (tid < 0) {
            fprintf(stderr, "single: no contig %s in fixture\n", QUERIES[i].chrom);
            rc = -1; break;
        }
        hts_itr_t *itr = tbx_itr_queryi(tbx, tid, QUERIES[i].beg, QUERIES[i].end);
        if (!itr) {
            fprintf(stderr, "single: tbx_itr_queryi failed for %s:%lld-%lld\n",
                    QUERIES[i].chrom, (long long)QUERIES[i].beg, (long long)QUERIES[i].end);
            rc = -1; break;
        }
        while (tbx_itr_next(fp, tbx, itr, &s) >= 0) {
            append_line(out, s.s);
        }
        tbx_itr_destroy(itr);
    }
    free(s.s);
    tbx_destroy(tbx);
    hts_close(fp);
    return rc;
}

/* Multi-region path: build one reglist, one iterator via tbx_itr_regions,
 * drain it. */
static int run_multi_region(const char *filename, kstring_t *out)
{
    htsFile *fp = hts_open(filename, "r");
    if (!fp) { fprintf(stderr, "multi: hts_open failed: %s\n", filename); return -1; }
    tbx_t *tbx = tbx_index_load(filename);
    if (!tbx) { fprintf(stderr, "multi: tbx_index_load failed\n"); hts_close(fp); return -1; }

    /* Build a reglist: group queries by contig. hts_reglist_t holds intervals
     * for ONE tid each, so one reglist per distinct chrom. */
    hts_reglist_t *reglist = calloc(N_QUERIES, sizeof(hts_reglist_t));
    if (!reglist) {
        tbx_destroy(tbx); hts_close(fp); return -1;
    }
    int n_reg = 0;
    for (int i = 0; i < N_QUERIES; i++) {
        int tid = tbx_name2id(tbx, QUERIES[i].chrom);
        if (tid < 0) {
            fprintf(stderr, "multi: no contig %s in fixture\n", QUERIES[i].chrom);
            free(reglist); tbx_destroy(tbx); hts_close(fp); return -1;
        }
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
            reglist[idx].min_beg = QUERIES[i].beg;
            reglist[idx].max_end = QUERIES[i].end;
        }
        /* realloc into a temp so the original buffer is not lost on NULL. */
        hts_pair_pos_t *grown = realloc(reglist[idx].intervals,
                                         (reglist[idx].count + 1) * sizeof(hts_pair_pos_t));
        if (!grown) {
            for (int j = 0; j < n_reg; j++) free(reglist[j].intervals);
            free(reglist); tbx_destroy(tbx); hts_close(fp); return -1;
        }
        reglist[idx].intervals = grown;
        reglist[idx].intervals[reglist[idx].count].beg = QUERIES[i].beg;
        reglist[idx].intervals[reglist[idx].count].end = QUERIES[i].end;
        reglist[idx].count++;
        if (QUERIES[i].beg < reglist[idx].min_beg) reglist[idx].min_beg = QUERIES[i].beg;
        if (QUERIES[i].end > reglist[idx].max_end) reglist[idx].max_end = QUERIES[i].end;
    }

    hts_itr_t *itr = tbx_itr_regions(fp, tbx, reglist, n_reg);
    if (!itr) {
        fprintf(stderr, "multi: tbx_itr_regions returned NULL\n");
        for (int j = 0; j < n_reg; j++) free(reglist[j].intervals);
        free(reglist);
        tbx_destroy(tbx); hts_close(fp); return -1;
    }

    kstring_t s = {0, 0, 0};
    while (hts_itr_multi_next(fp, itr, &s) >= 0) {
        append_line(out, s.s);
    }

    free(s.s);
    hts_itr_destroy(itr);
    /* Caller always owns the original reglist (tbx_itr_regions deep copies). */
    for (int j = 0; j < n_reg; j++) free(reglist[j].intervals);
    free(reglist);
    tbx_destroy(tbx);
    hts_close(fp);
    return 0;
}

/* Dedup test: tbx_itr_regions documents that overlapping and duplicate
 * intervals yield each matching record at most once. Verify by constructing
 * a reglist with three overlapping intervals that together cover the same
 * three records four times over (each record is in two intervals at least),
 * and asserting we still get exactly three records back. */
static int test_dedup_and_overlap(const char *filename)
{
    htsFile *fp = hts_open(filename, "r");
    if (!fp) return -1;
    tbx_t *tbx = tbx_index_load(filename);
    if (!tbx) { hts_close(fp); return -1; }

    int tid = tbx_name2id(tbx, "1");
    if (tid < 0) { tbx_destroy(tbx); hts_close(fp); return -1; }

    /* Three intervals on contig 1, deliberately overlapping. The fixture
     * has records at positions 100, 150, 200, 300 (1-based) on contig 1.
     *   [99, 200) covers records a (100), b (150), c (200-boundary)
     *   [120, 250) covers b, c
     *   [99, 100) covers a (full duplicate of the start of interval 1)
     * If dedup is working, we expect exactly 3 records: a, b, c. */
    hts_reglist_t *reg = calloc(1, sizeof(hts_reglist_t));
    if (!reg) { tbx_destroy(tbx); hts_close(fp); return -1; }
    reg[0].tid = tid;
    reg[0].count = 3;
    reg[0].intervals = calloc(3, sizeof(hts_pair_pos_t));
    if (!reg[0].intervals) { free(reg); tbx_destroy(tbx); hts_close(fp); return -1; }
    reg[0].intervals[0].beg = 99;  reg[0].intervals[0].end = 200;
    reg[0].intervals[1].beg = 120; reg[0].intervals[1].end = 250;
    reg[0].intervals[2].beg = 99;  reg[0].intervals[2].end = 100;
    reg[0].min_beg = 99; reg[0].max_end = 250;

    hts_itr_t *itr = tbx_itr_regions(fp, tbx, reg, 1);
    if (!itr) {
        fprintf(stderr, "FAIL: dedup test could not construct iterator\n");
        free(reg[0].intervals); free(reg);
        tbx_destroy(tbx); hts_close(fp);
        return -1;
    }

    kstring_t s = {0, 0, 0};
    int n = 0;
    while (hts_itr_multi_next(fp, itr, &s) >= 0) n++;
    free(s.s);
    hts_itr_destroy(itr);
    /* Caller always owns the original reglist (tbx_itr_regions deep copies). */
    free(reg[0].intervals); free(reg);
    tbx_destroy(tbx);
    hts_close(fp);

    /* Three records expected: a, b, c. If dedup were broken we'd get 5+. */
    if (n != 3) {
        fprintf(stderr, "FAIL: dedup expected 3 records, got %d\n", n);
        return -1;
    }
    fprintf(stderr, "OK: dedup yields exactly 3 records across overlapping intervals\n");
    return 0;
}

/* Failure-path test: confirm tbx_itr_regions returns NULL cleanly when given
 * bad inputs.
 *
 * Per tbx.h, the caller ALWAYS owns the original reglist regardless of
 * outcome (the function deep copies internally). This test exercises the
 * input-validation guard plus a successful-then-freed case, and the
 * ASan-instrumented build will catch any double-free or use-after-free. */
static int test_failure_paths(const char *filename)
{
    htsFile *fp = hts_open(filename, "r");
    if (!fp) return -1;
    tbx_t *tbx = tbx_index_load(filename);
    if (!tbx) { hts_close(fp); return -1; }

    hts_reglist_t *reg = calloc(1, sizeof(hts_reglist_t));
    if (!reg) { tbx_destroy(tbx); hts_close(fp); return -1; }
    reg[0].tid = tbx_name2id(tbx, "1");
    reg[0].count = 1;
    reg[0].intervals = calloc(1, sizeof(hts_pair_pos_t));
    if (!reg[0].intervals) { free(reg); tbx_destroy(tbx); hts_close(fp); return -1; }
    reg[0].intervals[0].beg = 99;
    reg[0].intervals[0].end = 100;
    reg[0].min_beg = 99;
    reg[0].max_end = 100;

    int rc = 0;
    hts_itr_t *itr;

    struct { const char *name; htsFile *fp; tbx_t *tbx; hts_reglist_t *reg; int count; } cases[] = {
        { "NULL htsFile",     NULL, tbx,  reg,  1 },
        { "NULL tbx",         fp,   NULL, reg,  1 },
        { "NULL fp+tbx",      NULL, NULL, reg,  1 },
        { "NULL reglist",     fp,   tbx,  NULL, 1 },
        { "count == 0",       fp,   tbx,  reg,  0 },
        { "negative count",   fp,   tbx,  reg, -1 },
    };
    const int n_cases = sizeof(cases) / sizeof(cases[0]);

    for (int i = 0; i < n_cases; i++) {
        itr = tbx_itr_regions(cases[i].fp, cases[i].tbx, cases[i].reg, cases[i].count);
        if (itr != NULL) {
            fprintf(stderr, "FAIL: %s should return NULL, got non-NULL\n", cases[i].name);
            hts_itr_destroy(itr);
            rc = 1;
        } else {
            fprintf(stderr, "OK: %s returns NULL\n", cases[i].name);
        }
    }

    /* Test the .reg name-string resolution paths added by the deep-copy
     * refactor. Both should succeed in constructing an iterator: one
     * resolves to a known tid, one warns and resolves to a negative tid
     * (the iterator still constructs but yields no records for that
     * entry). Caller still owns reg afterwards. */
    {
        hts_reglist_t *named = calloc(1, sizeof(*named));
        named[0].reg = "1";
        named[0].count = 1;
        named[0].intervals = calloc(1, sizeof(hts_pair_pos_t));
        named[0].intervals[0].beg = 99;
        named[0].intervals[0].end = 100;
        named[0].min_beg = 99; named[0].max_end = 100;

        itr = tbx_itr_regions(fp, tbx, named, 1);
        if (!itr) {
            fprintf(stderr, "FAIL: .reg=\"1\" should resolve and succeed\n");
            rc = 1;
        } else {
            fprintf(stderr, "OK: .reg=\"1\" resolved via tbx_name2id\n");
            hts_itr_destroy(itr);
        }
        free(named[0].intervals); free(named);
    }
    {
        hts_reglist_t *named = calloc(1, sizeof(*named));
        named[0].reg = "nonexistent_contig";
        named[0].count = 1;
        named[0].intervals = calloc(1, sizeof(hts_pair_pos_t));
        named[0].intervals[0].beg = 99;
        named[0].intervals[0].end = 100;
        named[0].min_beg = 99; named[0].max_end = 100;

        itr = tbx_itr_regions(fp, tbx, named, 1);
        if (!itr) {
            fprintf(stderr, "OK: unknown .reg name returns NULL (acceptable)\n");
        } else {
            fprintf(stderr, "OK: unknown .reg name resolved to negative tid"
                            " and iterator constructed (will yield nothing)\n");
            hts_itr_destroy(itr);
        }
        free(named[0].intervals); free(named);
    }

    /* The ASan-instrumented build catches any double-free or
     * use-after-free; the deep-copy contract means this trailing free is
     * always safe regardless of which paths above hit. */
    free(reg[0].intervals);
    free(reg);

    tbx_destroy(tbx);
    hts_close(fp);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s file.tabix.gz\n", argv[0]);
        return 1;
    }
    const char *filename = argv[1];

    kstring_t single_out = {0, 0, 0};
    kstring_t multi_out  = {0, 0, 0};

    if (run_single_region(filename, &single_out) < 0) {
        fprintf(stderr, "single-region run failed\n");
        free(single_out.s); free(multi_out.s);
        return 2;
    }
    if (run_multi_region(filename, &multi_out) < 0) {
        fprintf(stderr, "multi-region run failed\n");
        free(single_out.s); free(multi_out.s);
        return 2;
    }

    fprintf(stderr, "single-region output (%zu bytes):\n%s", single_out.l, single_out.s);
    fprintf(stderr, "multi-region output  (%zu bytes):\n%s", multi_out.l, multi_out.s);

    int rc = 0;
    if (single_out.l != multi_out.l || memcmp(single_out.s, multi_out.s, single_out.l) != 0) {
        fprintf(stderr, "PARITY MISMATCH between single-region and multi-region outputs\n");
        rc = 1;
    } else {
        fprintf(stderr, "PARITY OK: %zu bytes identical\n", single_out.l);
    }

    /* Dedup coverage (only for BED fixture, which has the positions the
     * test expects; VCF/GFF fixtures share positions so it works for them
     * too, but the constants in the test are BED-centric). */
    if (test_dedup_and_overlap(filename) != 0) {
        fprintf(stderr, "DEDUP/OVERLAP TESTS failed\n");
        rc = 1;
    } else {
        fprintf(stderr, "DEDUP/OVERLAP TESTS OK\n");
    }

    /* Failure-path coverage */
    if (test_failure_paths(filename) != 0) {
        fprintf(stderr, "FAILURE-PATH TESTS failed\n");
        rc = 1;
    } else {
        fprintf(stderr, "FAILURE-PATH TESTS OK\n");
    }

    free(single_out.s);
    free(multi_out.s);
    return rc;
}
