/*  test/test_format_plan_cache.c -- FORMAT planner cache tests.

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
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../htslib/kstring.h"
#include "../htslib/vcf.h"

void hts_vcf_format_plan_stats(uint64_t *attempts, uint64_t *hits,
                               uint64_t *fallback, uint64_t *parsed_samples);

static void fail(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

#define check0(expr) do { if ((expr) != 0) fail("check failed: " #expr); } while (0)
#define check1(expr) do { if (!(expr)) fail("check failed: " #expr); } while (0)

static void parse_line(bcf_hdr_t *hdr, bcf1_t *rec, kstring_t *line,
                       const char *text)
{
    ks_clear(line);
    if (kputsn(text, strlen(text), line) < 0)
        fail("failed to build VCF line");
    check0(vcf_parse(line, hdr, rec));
}

static void check_x_values(bcf_hdr_t *hdr, bcf1_t *rec,
                           const int32_t *expected, int n_expected)
{
    int32_t *values = NULL;
    int n_values = 0, ret, i;

    ret = bcf_get_format_int32(hdr, rec, "X", &values, &n_values);
    if (ret != n_expected) {
        free(values);
        fail("unexpected X vector length");
    }
    for (i = 0; i < n_expected; i++) {
        if (values[i] != expected[i]) {
            free(values);
            fail("unexpected X value");
        }
    }
    free(values);
}

int main(void)
{
    static char header[] =
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=1>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"cache generation test\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
    static const int32_t x1[] = { 7 };
    static const int32_t x2[] = { 11, 13 };
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    kstring_t line = KS_INITIALIZE;
    uint64_t attempts = 0, hits = 0, fallback = 0, parsed_samples = 0;

    check0(setenv("HTS_VCF_FORMAT_PLAN", "1", 1));
    hdr = bcf_hdr_init("r");
    rec = bcf_init();
    check1(hdr);
    check1(rec);
    check0(bcf_hdr_parse(hdr, header));

    parse_line(hdr, rec, &line,
               "1\t1\t.\tA\tC\t.\tPASS\t.\tGT:X\t0/1:7");
    check_x_values(hdr, rec, x1, 1);

    /*
     * Rebuild the same FORMAT string against changed metadata.  A stale plan
     * would still think X is Number=1 and would either fall back or encode the
     * second row incorrectly.  The header-owned generation must force a fresh
     * compile, preserving both correctness and fast-path coverage.
     */
    bcf_hdr_remove(hdr, BCF_HL_FMT, "X");
    check0(bcf_hdr_append(hdr,
                          "##FORMAT=<ID=X,Number=2,Type=Integer,Description=\"cache generation test\">"));
    check0(bcf_hdr_sync(hdr));
    bcf_clear1(rec);
    parse_line(hdr, rec, &line,
               "1\t2\t.\tA\tC\t.\tPASS\t.\tGT:X\t0/1:11,13");
    check_x_values(hdr, rec, x2, 2);

    hts_vcf_format_plan_stats(&attempts, &hits, &fallback, &parsed_samples);
    if (attempts != 2 || hits != 2 || fallback != 0 || parsed_samples != 2) {
        fprintf(stderr,
                "unexpected planner stats: attempts=%" PRIu64
                " hits=%" PRIu64 " fallback=%" PRIu64
                " parsed_samples=%" PRIu64 "\n",
                attempts, hits, fallback, parsed_samples);
        return EXIT_FAILURE;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    free(line.s);
    return EXIT_SUCCESS;
}
