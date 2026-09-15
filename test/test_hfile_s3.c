/*  test/test_hfile_s3.c -- Integration tests for the hfile_s3 backend,
    exercised against a real S3-compatible server (MinIO in CI).

    Copyright (C) 2026 Broad Institute.

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

/*
 * This is not part of the default `make check` suite -- it needs a live
 * S3-compatible endpoint (see test/docker/minio/ and `make check-s3`).  It
 * self-skips (exit 0) unless HTS_S3_HOST is set, so it's safe to leave the
 * binary built without affecting contributors who don't have Docker/MinIO.
 *
 * See test/s3_io_inventory.md for the file-I/O survey this suite is based
 * on, and for which findings are exercised here vs. documented only.
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include "../htslib/hfile.h"
#include "../htslib/hts.h"
#include "../htslib/sam.h"
#include "../htslib/vcf.h"
#include "../htslib/tbx.h"
#include "../htslib/faidx.h"
#include "../htslib/bgzf.h"
#include "../htslib/cram.h"

static int failures = 0;
static char bucket[256];
static char run_id[64];

#define PASS(name) fprintf(stderr, "  PASS: %s\n", (name))
#define FAIL(name, ...) do { \
    fprintf(stderr, "  FAIL: %s: ", (name)); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    failures++; \
} while (0)

// A trivial "4M" CIGAR, needed because bam_set1() requires a CIGAR for any
// mapped record (tid >= 0).
static const uint32_t cigar_4m[1] = { bam_cigar_gen(4, BAM_CMATCH) };

// Builds "s3+http://<bucket>/<run_id>/<key>" into buf.
static void s3_url(char *buf, size_t bufsz, const char *key)
{
    snprintf(buf, bufsz, "s3+http://%s/%s/%s", bucket, run_id, key);
}

static void generate_bytes(unsigned char *buf, size_t n, unsigned seed)
{
    size_t i;
    for (i = 0; i < n; i++)
        buf[i] = (unsigned char)((i * 7 + seed) & 0xFF);
}

/* ---------------------------------------------------------------------- */
/* 1. Raw hFILE read/write/seek                                            */
/* ---------------------------------------------------------------------- */

static void test_hfile_roundtrip(void)
{
    const char *name = "hFILE roundtrip";
    char url[512];
    unsigned char wbuf[5000], rbuf[5000];
    hFILE *fp;
    ssize_t n;
    size_t total;

    s3_url(url, sizeof(url), "hfile_roundtrip.bin");
    generate_bytes(wbuf, sizeof(wbuf), 13);

    fp = hopen(url, "w");
    if (!fp) { FAIL(name, "hopen(w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, wbuf, sizeof(wbuf)) != (ssize_t)sizeof(wbuf)) {
        FAIL(name, "hwrite failed: %s", strerror(errno));
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "hclose(w) failed: %s", strerror(errno)); return; }

    fp = hopen(url, "r");
    if (!fp) { FAIL(name, "hopen(r) failed: %s", strerror(errno)); return; }
    total = 0;
    while (total < sizeof(rbuf) &&
           (n = hread(fp, rbuf + total, sizeof(rbuf) - total)) > 0)
        total += n;
    hclose_abruptly(fp);

    if (total != sizeof(wbuf))
        FAIL(name, "size mismatch: got %zu, expected %zu", total, sizeof(wbuf));
    else if (memcmp(wbuf, rbuf, sizeof(wbuf)) != 0)
        FAIL(name, "data mismatch");
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 2. BGZF write/read + .gzi sidecar index                                 */
/* ---------------------------------------------------------------------- */

static void test_bgzf_roundtrip(void)
{
    const char *name = "BGZF roundtrip + .gzi index";
    char url[512];
    unsigned char wbuf[200000], rbuf[200000];
    BGZF *bgzf;
    ssize_t n;
    size_t total;

    s3_url(url, sizeof(url), "bgzf_roundtrip.gz");
    generate_bytes(wbuf, sizeof(wbuf), 41);

    bgzf = bgzf_open(url, "w");
    if (!bgzf) { FAIL(name, "bgzf_open(w) failed: %s", strerror(errno)); return; }
    if (bgzf_write(bgzf, wbuf, sizeof(wbuf)) != (ssize_t)sizeof(wbuf)) {
        FAIL(name, "bgzf_write failed");
        bgzf_close(bgzf);
        return;
    }
    if (bgzf_close(bgzf) != 0) { FAIL(name, "bgzf_close(w) failed"); return; }

    // Plain read-back sanity check.
    bgzf = bgzf_open(url, "r");
    if (!bgzf) { FAIL(name, "bgzf_open(r) failed: %s", strerror(errno)); return; }
    total = 0;
    while (total < sizeof(rbuf) &&
           (n = bgzf_read(bgzf, rbuf + total, sizeof(rbuf) - total)) > 0)
        total += n;
    bgzf_close(bgzf);
    if (total != sizeof(wbuf) || memcmp(wbuf, rbuf, sizeof(wbuf)) != 0) {
        FAIL(name, "plain read-back mismatch (got %zu bytes)", total);
        return;
    }

    // Build the .gzi index by scanning on a read pass (mirrors fai_build3_core).
    bgzf = bgzf_open(url, "r");
    if (!bgzf) { FAIL(name, "bgzf_open(r) for indexing failed"); return; }
    if (bgzf_index_build_init(bgzf) != 0) {
        FAIL(name, "bgzf_index_build_init failed");
        bgzf_close(bgzf);
        return;
    }
    total = 0;
    while ((n = bgzf_read(bgzf, rbuf, sizeof(rbuf))) > 0)
        total += n;
    if (bgzf_index_dump(bgzf, url, ".gzi") < 0) {
        FAIL(name, "bgzf_index_dump failed: %s", strerror(errno));
        bgzf_close(bgzf);
        return;
    }
    bgzf_close(bgzf);

    // Load the index into a fresh handle and seek into the middle.
    bgzf = bgzf_open(url, "r");
    if (!bgzf) { FAIL(name, "bgzf_open(r) for seek test failed"); return; }
    if (bgzf_index_load(bgzf, url, ".gzi") < 0) {
        FAIL(name, "bgzf_index_load failed: %s", strerror(errno));
        bgzf_close(bgzf);
        return;
    }
    {
        off_t mid = sizeof(wbuf) / 2;
        size_t want = sizeof(wbuf) - mid;
        if (bgzf_useek(bgzf, mid, SEEK_SET) < 0) {
            FAIL(name, "bgzf_useek failed: %s", strerror(errno));
            bgzf_close(bgzf);
            return;
        }
        total = 0;
        while (total < want &&
               (n = bgzf_read(bgzf, rbuf + total, want - total)) > 0)
            total += n;
        bgzf_close(bgzf);

        if (total != want)
            FAIL(name, "post-seek size mismatch: got %zu, expected %zu", total, want);
        else if (memcmp(wbuf + mid, rbuf, want) != 0)
            FAIL(name, "post-seek data mismatch");
        else
            PASS(name);
    }
}

/* ---------------------------------------------------------------------- */
/* 3. SAM/BAM main path + index build/load                                 */
/* ---------------------------------------------------------------------- */

static void test_sam_bam_roundtrip(void)
{
    const char *name = "SAM/BAM roundtrip + index";
    char url[512];
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int tid, count, ret;

    s3_url(url, sizeof(url), "test.bam");

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wb");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (sam_hdr_write(out, hdr) < 0) {
        FAIL(name, "sam_hdr_write failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }

    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 100, 60, 1, cigar_4m, -1, -1, 0,
                    4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read1) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    ret = bam_set1(rec, 4, "read2", 0, 0, 500, 60, 1, cigar_4m, -1, -1, 0,
                    4, "TTTT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read2) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    if (sam_index_build3(url, NULL, 0, 0) < 0) {
        FAIL(name, "sam_index_build3 failed: %s", strerror(errno));
        return;
    }

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read failed"); hts_close(in); return; }
    idx = sam_index_load3(in, url, NULL, 0);
    if (!idx) {
        FAIL(name, "sam_index_load3 failed: %s", strerror(errno));
        sam_hdr_destroy(hdr2); hts_close(in);
        return;
    }

    tid = sam_hdr_name2tid(hdr2, "chr1");
    itr = sam_itr_queryi(idx, tid, 0, 1000);
    rec = bam_init1();
    count = 0;
    while ((ret = sam_itr_next(in, itr, rec)) >= 0)
        count++;

    hts_itr_destroy(itr);
    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr2);
    hts_close(in);

    if (count != 2)
        FAIL(name, "expected 2 records via index query, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 4. BCF main path + .csi index                                           */
/* ---------------------------------------------------------------------- */

static void test_bcf_roundtrip(void)
{
    const char *name = "BCF roundtrip + index";
    char url[512];
    htsFile *out, *in;
    bcf_hdr_t *hdr, *hdr2;
    bcf1_t *rec;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int count, ret;

    s3_url(url, sizeof(url), "test.bcf");

    hdr = bcf_hdr_init("w");
    if (!hdr) { FAIL(name, "bcf_hdr_init failed"); return; }
    if (bcf_hdr_append(hdr, "##contig=<ID=chr1,length=1000>") < 0) {
        FAIL(name, "bcf_hdr_append failed");
        bcf_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wb");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); bcf_hdr_destroy(hdr); return; }
    if (bcf_hdr_write(out, hdr) < 0) {
        FAIL(name, "bcf_hdr_write failed");
        hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }

    rec = bcf_init();
    rec->rid = bcf_hdr_name2id(hdr, "chr1");
    rec->pos = 99;
    rec->qual = 50;
    if (bcf_update_alleles_str(hdr, rec, "A,G") < 0 || bcf_write1(out, hdr, rec) < 0) {
        FAIL(name, "bcf_write1(rec1) failed");
        bcf_destroy(rec); hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }
    bcf_clear1(rec);
    rec->rid = bcf_hdr_name2id(hdr, "chr1");
    rec->pos = 199;
    rec->qual = 30;
    if (bcf_update_alleles_str(hdr, rec, "C,T") < 0 || bcf_write1(out, hdr, rec) < 0) {
        FAIL(name, "bcf_write1(rec2) failed");
        bcf_destroy(rec); hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }
    bcf_destroy(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); bcf_hdr_destroy(hdr); return; }
    bcf_hdr_destroy(hdr);

    if (bcf_index_build3(url, NULL, 14, 0) < 0) {
        FAIL(name, "bcf_index_build3 failed: %s", strerror(errno));
        return;
    }

    in = hts_open(url, "rb");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    hdr2 = bcf_hdr_read(in);
    if (!hdr2) { FAIL(name, "bcf_hdr_read failed"); hts_close(in); return; }
    idx = bcf_index_load3(url, NULL, 0);
    if (!idx) {
        FAIL(name, "bcf_index_load3 failed: %s", strerror(errno));
        bcf_hdr_destroy(hdr2); hts_close(in);
        return;
    }

    itr = bcf_itr_queryi(idx, bcf_hdr_name2id(hdr2, "chr1"), 0, 1000);
    rec = bcf_init();
    count = 0;
    while ((ret = bcf_itr_next(in, itr, rec)) >= 0)
        count++;

    hts_itr_destroy(itr);
    bcf_destroy(rec);
    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr2);
    hts_close(in);

    if (count != 2)
        FAIL(name, "expected 2 records via index query, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 5. Generic tabix path (BED), independent of BCF                         */
/* ---------------------------------------------------------------------- */

static void test_tabix_roundtrip(void)
{
    const char *name = "tabix roundtrip (BED)";
    char url[512];
    static const char bed[] =
        "chr1\t10\t20\tfeatureA\n"
        "chr1\t30\t40\tfeatureB\n";
    BGZF *bgzf;
    htsFile *in;
    tbx_t *tbx;
    hts_itr_t *itr;
    kstring_t str = KS_INITIALIZE;
    int tid, count, ret;

    s3_url(url, sizeof(url), "test.bed.gz");

    bgzf = bgzf_open(url, "w");
    if (!bgzf) { FAIL(name, "bgzf_open(w) failed: %s", strerror(errno)); return; }
    if (bgzf_write(bgzf, bed, sizeof(bed) - 1) != (ssize_t)(sizeof(bed) - 1)) {
        FAIL(name, "bgzf_write failed");
        bgzf_close(bgzf);
        return;
    }
    if (bgzf_close(bgzf) != 0) { FAIL(name, "bgzf_close(w) failed"); return; }

    if (tbx_index_build3(url, NULL, 0, 0, &tbx_conf_bed) < 0) {
        FAIL(name, "tbx_index_build3 failed: %s", strerror(errno));
        return;
    }

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    tbx = tbx_index_load3(url, NULL, 0);
    if (!tbx) {
        FAIL(name, "tbx_index_load3 failed: %s", strerror(errno));
        hts_close(in);
        return;
    }

    tid = tbx_name2id(tbx, "chr1");
    itr = tbx_itr_queryi(tbx, tid, 0, 100);
    count = 0;
    while ((ret = tbx_itr_next(in, tbx, itr, &str)) >= 0)
        count++;

    free(str.s);
    hts_itr_destroy(itr);
    tbx_destroy(tbx);
    hts_close(in);

    if (count != 2)
        FAIL(name, "expected 2 records via tabix query, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 6. FASTA/faidx path                                                     */
/* ---------------------------------------------------------------------- */

static void test_faidx_roundtrip(void)
{
    const char *name = "faidx roundtrip";
    char url[512];
    static const char fasta[] =
        ">chr1\nACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n";
    hFILE *fp;
    faidx_t *fai;
    char *seq;
    int len;

    s3_url(url, sizeof(url), "test.fa");

    fp = hopen(url, "w");
    if (!fp) { FAIL(name, "hopen(w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, fasta, sizeof(fasta) - 1) != (ssize_t)(sizeof(fasta) - 1)) {
        FAIL(name, "hwrite failed");
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "hclose(w) failed"); return; }

    if (fai_build3(url, NULL, NULL) < 0) {
        FAIL(name, "fai_build3 failed: %s", strerror(errno));
        return;
    }

    fai = fai_load3(url, NULL, NULL, 0);
    if (!fai) {
        FAIL(name, "fai_load3 failed: %s", strerror(errno));
        return;
    }

    seq = fai_fetch(fai, "chr1:1-10", &len);
    if (!seq) {
        FAIL(name, "fai_fetch failed");
    } else if (len != 10 || strcmp(seq, "ACGTACGTAC") != 0) {
        FAIL(name, "fai_fetch mismatch: len=%d seq=%s", len, seq);
    } else {
        PASS(name);
    }
    free(seq);
    fai_destroy(fai);
}

/* ---------------------------------------------------------------------- */
/* 7. CRAM main path, reference-free (CRAM_OPT_NO_REF)                     */
/* ---------------------------------------------------------------------- */

static void test_cram_roundtrip(void)
{
    const char *name = "CRAM roundtrip (no_ref)";
    char url[512];
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    int count, ret;

    s3_url(url, sizeof(url), "test.cram");

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wc");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (hts_set_opt(out, CRAM_OPT_NO_REF, 1) < 0) {
        FAIL(name, "hts_set_opt(CRAM_OPT_NO_REF) failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    if (sam_hdr_write(out, hdr) < 0) {
        FAIL(name, "sam_hdr_write failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }

    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 10, 60, 1, cigar_4m, -1, -1, 0,
                    4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read1) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    ret = bam_set1(rec, 4, "read2", 0, 0, 20, 60, 1, cigar_4m, -1, -1, 0,
                    4, "GGCC", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read2) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read failed"); hts_close(in); return; }

    rec = bam_init1();
    count = 0;
    while ((ret = sam_read1(in, hdr2, rec)) >= 0)
        count++;

    bam_destroy1(rec);
    sam_hdr_destroy(hdr2);
    hts_close(in);

    if (ret < -1)
        FAIL(name, "sam_read1 returned error %d", ret);
    else if (count != 2)
        FAIL(name, "expected 2 records, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 8. CRAM with an explicit S3-hosted reference (CRAM_OPT_REFERENCE)       */
/* ---------------------------------------------------------------------- */

static void test_cram_explicit_s3_reference(void)
{
    const char *name = "CRAM roundtrip (explicit S3 reference)";
    char cram_url[512], ref_url[512];
    static const char fasta[] =
        ">chr1\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    hFILE *fp;
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    int count, ret;

    s3_url(cram_url, sizeof(cram_url), "test_with_ref.cram");
    s3_url(ref_url, sizeof(ref_url), "ref.fa");

    // Upload the reference and build its .fai, all on S3 -- see
    // test_faidx_roundtrip for the same pattern.
    fp = hopen(ref_url, "w");
    if (!fp) { FAIL(name, "hopen(ref, w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, fasta, sizeof(fasta) - 1) != (ssize_t)(sizeof(fasta) - 1)) {
        FAIL(name, "hwrite(ref) failed");
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "hclose(ref, w) failed"); return; }
    if (fai_build3(ref_url, NULL, NULL) < 0) {
        FAIL(name, "fai_build3(ref) failed: %s", strerror(errno));
        return;
    }

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "64", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    out = hts_open(cram_url, "wc");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (hts_set_opt(out, CRAM_OPT_REFERENCE, ref_url) < 0) {
        FAIL(name, "hts_set_opt(CRAM_OPT_REFERENCE, S3 url) failed for write: %s",
             strerror(errno));
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    if (sam_hdr_write(out, hdr) < 0) {
        FAIL(name, "sam_hdr_write failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 5, 60, 1, cigar_4m, -1, -1, 0,
                    4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1 failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    // Read back, explicitly pointing the CRAM reference resolver at the
    // *S3* reference -- this proves cram_load_reference()/refs_load_fai()
    // (cram/cram_io.c) fetch reference bytes through hFILE when told to,
    // unlike a CRAM's embedded UR: tag which htslib refuses to follow for
    // non-file:// URLs (see test/s3_io_inventory.md).
    in = hts_open(cram_url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    if (hts_set_opt(in, CRAM_OPT_REFERENCE, ref_url) < 0) {
        FAIL(name, "hts_set_opt(CRAM_OPT_REFERENCE, S3 url) failed for read: %s",
             strerror(errno));
        hts_close(in);
        return;
    }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read failed"); hts_close(in); return; }

    rec = bam_init1();
    count = 0;
    while ((ret = sam_read1(in, hdr2, rec)) >= 0)
        count++;

    bam_destroy1(rec);
    sam_hdr_destroy(hdr2);
    hts_close(in);

    if (ret < -1)
        FAIL(name, "sam_read1 returned error %d", ret);
    else if (count != 1)
        FAIL(name, "expected 1 record, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 9. Bypass proof: local-CWD index fallback fires for a remote main file  */
/*    (hts_idx_check_local, hts.c) -- see test/s3_io_inventory.md category 1 */
/* ---------------------------------------------------------------------- */

static void test_index_local_fallback_bypass(void)
{
    const char *name = "bypass: local .bai fallback for S3 main file";
    char url[512];
    char local_bam[] = "/tmp/htslib_s3_bypass_XXXXXX";
    char local_bam_path[300], local_bai_path[310];
    char cwd_bai_path[300];
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int tid, count, ret, tmpfd;

    s3_url(url, sizeof(url), "bypass/local_idx.bam");

    // Build identical BAM content both on S3 (no index built there) and in
    // a local temp file (index built locally), then plant the local index
    // in CWD under the S3 file's basename so hts_idx_check_local() will
    // find it purely by local stat().
    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    tmpfd = mkstemp(local_bam);
    if (tmpfd < 0) { FAIL(name, "mkstemp failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    close(tmpfd);
    snprintf(local_bam_path, sizeof(local_bam_path), "%s", local_bam);

    out = hts_open(url, "wb");
    if (!out) { FAIL(name, "hts_open(S3 w) failed: %s", strerror(errno)); goto cleanup_hdr; }
    if (sam_hdr_write(out, hdr) < 0) { FAIL(name, "sam_hdr_write(S3) failed"); hts_close(out); goto cleanup_hdr; }
    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 100, 60, 1, cigar_4m, -1, -1, 0, 4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) { FAIL(name, "sam_write1(S3) failed"); bam_destroy1(rec); hts_close(out); goto cleanup_hdr; }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(S3 w) failed"); goto cleanup_hdr; }

    out = hts_open(local_bam_path, "wb");
    if (!out) { FAIL(name, "hts_open(local w) failed: %s", strerror(errno)); goto cleanup_hdr; }
    if (sam_hdr_write(out, hdr) < 0) { FAIL(name, "sam_hdr_write(local) failed"); hts_close(out); goto cleanup_hdr; }
    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 100, 60, 1, cigar_4m, -1, -1, 0, 4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) { FAIL(name, "sam_write1(local) failed"); bam_destroy1(rec); hts_close(out); goto cleanup_hdr; }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(local w) failed"); goto cleanup_hdr; }

    if (sam_index_build3(local_bam_path, NULL, 0, 0) < 0) {
        FAIL(name, "sam_index_build3(local) failed: %s", strerror(errno));
        goto cleanup_local;
    }
    snprintf(local_bai_path, sizeof(local_bai_path), "%s.bai", local_bam_path);

    // Plant it in CWD as "<basename-of-S3-file>.bai".
    snprintf(cwd_bai_path, sizeof(cwd_bai_path), "local_idx.bam.bai");
    {
        FILE *src = fopen(local_bai_path, "rb");
        FILE *dst = fopen(cwd_bai_path, "wb");
        char buf[65536];
        size_t n;
        if (!src || !dst) {
            FAIL(name, "copying local index into CWD failed: %s", strerror(errno));
            if (src) fclose(src);
            if (dst) fclose(dst);
            goto cleanup_local;
        }
        while ((n = fread(buf, 1, sizeof(buf), src)) > 0)
            fwrite(buf, 1, n, dst);
        fclose(src);
        fclose(dst);
    }

    // Sanity check: confirm S3 genuinely has no index of its own.
    {
        char idx_url[560];
        hFILE *probe;
        snprintf(idx_url, sizeof(idx_url), "%s.bai", url);
        probe = hopen(idx_url, "r");
        if (probe) {
            FAIL(name, "unexpected: S3 already has an index at %s", idx_url);
            hclose_abruptly(probe);
            unlink(cwd_bai_path);
            goto cleanup_local;
        }
    }

    // The actual bypass: open the S3 BAM, load its index -- it should
    // succeed by silently reading the local CWD copy, not by fetching (or
    // failing to fetch) anything from S3.
    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(S3 r) failed: %s", strerror(errno)); unlink(cwd_bai_path); goto cleanup_local; }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read(S3) failed"); hts_close(in); unlink(cwd_bai_path); goto cleanup_local; }
    idx = sam_index_load3(in, url, NULL, 0);
    if (!idx) {
        FAIL(name, "sam_index_load3 unexpectedly failed -- local fallback did not fire");
        sam_hdr_destroy(hdr2); hts_close(in); unlink(cwd_bai_path);
        goto cleanup_local;
    }

    tid = sam_hdr_name2tid(hdr2, "chr1");
    itr = sam_itr_queryi(idx, tid, 0, 1000);
    rec = bam_init1();
    count = 0;
    while ((ret = sam_itr_next(in, itr, rec)) >= 0)
        count++;
    hts_itr_destroy(itr);
    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr2);
    hts_close(in);
    unlink(cwd_bai_path);

    if (count != 1)
        FAIL(name, "index query via local fallback returned %d records, expected 1", count);
    else
        PASS(name);

cleanup_local:
    unlink(local_bam_path);
    unlink(local_bai_path);
cleanup_hdr:
    sam_hdr_destroy(hdr);
}

/* ---------------------------------------------------------------------- */
/* 10. Large file: exercises real multipart upload                        */
/* ---------------------------------------------------------------------- */

// hfile_s3.c's write-side part size defaults to MINIMUM_S3_WRITE_SIZE
// (5 MiB, hfile_s3.c:1498) and can only be raised via HTS_S3_PART_SIZE, never
// lowered -- 5 MiB is S3's own minimum multipart part size. s3_write()
// (hfile_s3.c:1728) uploads whatever is buffered as soon as it exceeds the
// part size, so three writes of ~6/6/2 MiB below force three separate
// upload_part() calls (two real parts plus a final short part), genuinely
// exercising InitiateMultipartUpload/UploadPart/CompleteMultipartUpload
// rather than a single-shot PUT. The read-back also exceeds the default 1
// MiB READ_PART_SIZE (hfile_s3.c:1968), so it's served by multiple ranged
// GETs too.
static void test_multipart_large_file(void)
{
    const char *name = "multipart large file (write + ranged read)";
    char url[512];
    static const size_t chunk1 = 6 * 1024 * 1024;
    static const size_t chunk2 = 6 * 1024 * 1024;
    static const size_t chunk3 = 2 * 1024 * 1024;
    const size_t total = chunk1 + chunk2 + chunk3;
    unsigned char *wbuf, *rbuf;
    hFILE *fp;
    ssize_t n;
    size_t off, got;

    s3_url(url, sizeof(url), "multipart_large.bin");

    wbuf = malloc(total);
    rbuf = malloc(total);
    if (!wbuf || !rbuf) {
        FAIL(name, "malloc failed for %zu-byte buffers", total);
        free(wbuf); free(rbuf);
        return;
    }
    generate_bytes(wbuf, total, 197);

    fp = hopen(url, "w");
    if (!fp) {
        FAIL(name, "hopen(w) failed: %s", strerror(errno));
        free(wbuf); free(rbuf);
        return;
    }

    off = 0;
    if (hwrite(fp, wbuf + off, chunk1) != (ssize_t)chunk1) {
        FAIL(name, "hwrite chunk1 failed: %s", strerror(errno));
        hclose_abruptly(fp); free(wbuf); free(rbuf);
        return;
    }
    off += chunk1;
    if (hwrite(fp, wbuf + off, chunk2) != (ssize_t)chunk2) {
        FAIL(name, "hwrite chunk2 failed: %s", strerror(errno));
        hclose_abruptly(fp); free(wbuf); free(rbuf);
        return;
    }
    off += chunk2;
    if (hwrite(fp, wbuf + off, chunk3) != (ssize_t)chunk3) {
        FAIL(name, "hwrite chunk3 failed: %s", strerror(errno));
        hclose_abruptly(fp); free(wbuf); free(rbuf);
        return;
    }

    if (hclose(fp) != 0) {
        FAIL(name, "hclose(w) failed (multipart completion): %s", strerror(errno));
        free(wbuf); free(rbuf);
        return;
    }

    fp = hopen(url, "r");
    if (!fp) {
        FAIL(name, "hopen(r) failed: %s", strerror(errno));
        free(wbuf); free(rbuf);
        return;
    }
    got = 0;
    while (got < total && (n = hread(fp, rbuf + got, total - got)) > 0)
        got += n;
    hclose_abruptly(fp);

    if (got != total)
        FAIL(name, "size mismatch: got %zu bytes, expected %zu", got, total);
    else if (memcmp(wbuf, rbuf, total) != 0)
        FAIL(name, "data mismatch after multipart roundtrip");
    else
        PASS(name);

    free(wbuf);
    free(rbuf);
}

/* ---------------------------------------------------------------------- */
/* 11. Unhappy paths                                                       */
/* ---------------------------------------------------------------------- */

static void test_unhappy_read_nonexistent_key(void)
{
    const char *name = "unhappy: read a key that was never written";
    char url[512];
    hFILE *fp;

    s3_url(url, sizeof(url), "does/not/exist.bin");
    fp = hopen(url, "r");
    if (fp) {
        FAIL(name, "hopen unexpectedly succeeded");
        hclose_abruptly(fp);
    } else if (errno != ENOENT) {
        FAIL(name, "expected ENOENT, got errno=%d (%s)", errno, strerror(errno));
    } else {
        PASS(name);
    }
}

static void test_unhappy_read_nonexistent_bucket(void)
{
    const char *name = "unhappy: read from a bucket that doesn't exist";
    char url[256];
    hFILE *fp;

    snprintf(url, sizeof(url), "s3+http://htslib-test-no-such-bucket/probe");
    fp = hopen(url, "r");
    if (fp) {
        FAIL(name, "hopen unexpectedly succeeded");
        hclose_abruptly(fp);
    } else if (errno != ENOENT) {
        FAIL(name, "expected ENOENT, got errno=%d (%s)", errno, strerror(errno));
    } else {
        PASS(name);
    }
}

static void test_unhappy_wrong_credentials(void)
{
    const char *name = "unhappy: wrong credentials";
    char url[512], bad_url[560];
    hFILE *fp;
    static const char data[] = "wrong-credentials-probe";

    // Write a real object with the (correct) ambient credentials...
    s3_url(url, sizeof(url), "wrong_creds_probe.bin");
    fp = hopen(url, "w");
    if (!fp) { FAIL(name, "setup: hopen(w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, data, sizeof(data) - 1) != (ssize_t)(sizeof(data) - 1)) {
        FAIL(name, "setup: hwrite failed: %s", strerror(errno));
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "setup: hclose failed: %s", strerror(errno)); return; }

    // ...then try to read it back with bogus credentials embedded directly
    // in the URL (s3[+SCHEME]://ID:SECRET@BUCKET/PATH, hfile_s3.c:559),
    // which overrides AWS_ACCESS_KEY_ID/AWS_SECRET_ACCESS_KEY for just this
    // one request without touching global env state.
    snprintf(bad_url, sizeof(bad_url), "s3+http://wrongid:wrongsecret@%s/%s/wrong_creds_probe.bin",
             bucket, run_id);
    fp = hopen(bad_url, "r");
    if (fp) {
        FAIL(name, "hopen unexpectedly succeeded with bad credentials");
        hclose_abruptly(fp);
    } else if (errno != EACCES && errno != EPERM) {
        FAIL(name, "expected EACCES/EPERM, got errno=%d (%s)", errno, strerror(errno));
    } else {
        PASS(name);
    }
}

static void test_unhappy_connection_refused(void)
{
    const char *name = "unhappy: connection refused (bad HTS_S3_HOST)";
    char url[512];
    hFILE *fp;
    const char *orig = getenv("HTS_S3_HOST");
    char *saved = orig ? strdup(orig) : NULL;

    // Port 1 is essentially never listening in a CI container -- fails fast
    // with ECONNREFUSED rather than hanging. Note HTS_RETRY_MAX isn't
    // relevant here: hfile_s3.c has no retry logic of its own to disable
    // (see test/s3_gaps.md) -- every request is a single curl_easy_perform().
    setenv("HTS_S3_HOST", "127.0.0.1:1", 1);

    s3_url(url, sizeof(url), "unreachable_probe.bin");
    fp = hopen(url, "r");

    if (saved) { setenv("HTS_S3_HOST", saved, 1); free(saved); }
    else unsetenv("HTS_S3_HOST");

    if (fp) {
        FAIL(name, "hopen unexpectedly succeeded against an unreachable host");
        hclose_abruptly(fp);
    } else {
        PASS(name);
    }
}

static void test_unhappy_hts_open_nonexistent(void)
{
    const char *name = "unhappy: hts_open() on a BAM that was never written";
    char url[512];
    htsFile *fp;

    s3_url(url, sizeof(url), "does/not/exist.bam");
    fp = hts_open(url, "r");
    if (fp) {
        FAIL(name, "hts_open unexpectedly succeeded");
        hts_close(fp);
    } else {
        PASS(name);
    }
}

/* ---------------------------------------------------------------------- */
/* 12. CRAM's own index format (.crai)                                     */
/* ---------------------------------------------------------------------- */

static void test_cram_index_roundtrip(void)
{
    const char *name = "CRAM .crai index build/load";
    char url[512];
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    hts_idx_t *idx;
    hts_itr_t *itr;
    int tid, count, ret;
    // cram_index_load() (cram/cram_index.c) always calls hts_idx_getfn(),
    // i.e. always passes HTS_IDX_SAVE_REMOTE -- unlike hts_idx_load3() for
    // BAI/CSI/TBI, it ignores the flags argument entirely and unconditionally
    // downloads+persists a local copy of the .crai (see test/s3_gaps.md).
    // We didn't ask for that; clean it up so it doesn't pollute the CWD.
    static const char local_crai[] = "test_indexed.cram.crai";

    s3_url(url, sizeof(url), "test_indexed.cram");

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wc");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (hts_set_opt(out, CRAM_OPT_NO_REF, 1) < 0) {
        FAIL(name, "hts_set_opt(CRAM_OPT_NO_REF) failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    if (sam_hdr_write(out, hdr) < 0) {
        FAIL(name, "sam_hdr_write failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 100, 60, 1, cigar_4m, -1, -1, 0, 4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read1) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    ret = bam_set1(rec, 4, "read2", 0, 0, 500, 60, 1, cigar_4m, -1, -1, 0, 4, "TTTT", "IIII", 0);
    if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
        FAIL(name, "sam_write1(read2) failed");
        bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    if (sam_index_build3(url, NULL, 0, 0) < 0) {
        FAIL(name, "sam_index_build3 (.crai) failed: %s", strerror(errno));
        return;
    }

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read failed"); hts_close(in); return; }
    idx = sam_index_load3(in, url, NULL, 0);
    unlink(local_crai);
    if (!idx) {
        FAIL(name, "sam_index_load3 (.crai) failed: %s", strerror(errno));
        sam_hdr_destroy(hdr2); hts_close(in);
        return;
    }

    // Query a region covering only read2, to prove the index is actually
    // being used for selective access rather than a full linear scan.
    tid = sam_hdr_name2tid(hdr2, "chr1");
    itr = sam_itr_queryi(idx, tid, 400, 1000);
    rec = bam_init1();
    count = 0;
    while ((ret = sam_itr_next(in, itr, rec)) >= 0)
        count++;

    hts_itr_destroy(itr);
    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr2);
    hts_close(in);

    if (count != 1)
        FAIL(name, "expected 1 record from indexed region query, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 13. SigV2 signing (read-only -- hfile_s3.c has no SigV2 write support)  */
/* ---------------------------------------------------------------------- */

static void test_sigv2_read(void)
{
    const char *name = "SigV2 read (HTS_S3_V2)";
    char url[512];
    hFILE *fp;
    static const char data[] = "sigv2-probe-data";
    unsigned char rbuf[64];
    ssize_t n;
    size_t got;

    // Written normally (default SigV4) -- hfile_s3.c's SigV2 path only
    // implements reads, see test/s3_gaps.md.
    s3_url(url, sizeof(url), "sigv2_probe.bin");
    fp = hopen(url, "w");
    if (!fp) { FAIL(name, "setup: hopen(w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, data, sizeof(data) - 1) != (ssize_t)(sizeof(data) - 1)) {
        FAIL(name, "setup: hwrite failed: %s", strerror(errno));
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "setup: hclose failed: %s", strerror(errno)); return; }

    setenv("HTS_S3_V2", "1", 1);
    fp = hopen(url, "r");
    got = 0;
    if (fp) {
        while (got < sizeof(rbuf) &&
               (n = hread(fp, rbuf + got, sizeof(rbuf) - got)) > 0)
            got += n;
        hclose_abruptly(fp);
    }
    unsetenv("HTS_S3_V2");

    if (!fp)
        FAIL(name, "hopen with HTS_S3_V2=1 failed: %s", strerror(errno));
    else if (got != sizeof(data) - 1 || memcmp(rbuf, data, got) != 0)
        FAIL(name, "data mismatch reading via SigV2");
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 14. Credentials sourced from a file, not env vars                       */
/* ---------------------------------------------------------------------- */

static void test_credentials_file_fallback(void)
{
    const char *name = "credentials file fallback (~/.aws/credentials style)";
    char cred_path[] = "/tmp/htslib_s3_creds_XXXXXX";
    char url[512];
    hFILE *fp;
    FILE *cf;
    int tmpfd;
    static const char data[] = "creds-file-probe";
    unsigned char rbuf[64];
    ssize_t n;
    size_t got;
    char *saved_id, *saved_secret, *saved_file;
    const char *v;

    tmpfd = mkstemp(cred_path);
    if (tmpfd < 0) { FAIL(name, "mkstemp failed: %s", strerror(errno)); return; }
    cf = fdopen(tmpfd, "w");
    if (!cf) { FAIL(name, "fdopen failed: %s", strerror(errno)); close(tmpfd); unlink(cred_path); return; }
    fprintf(cf, "[default]\n");
    fprintf(cf, "aws_access_key_id = %s\n", getenv("AWS_ACCESS_KEY_ID") ? getenv("AWS_ACCESS_KEY_ID") : "minioadmin");
    fprintf(cf, "aws_secret_access_key = %s\n", getenv("AWS_SECRET_ACCESS_KEY") ? getenv("AWS_SECRET_ACCESS_KEY") : "minioadmin");
    fclose(cf);

    v = getenv("AWS_ACCESS_KEY_ID");     saved_id = v ? strdup(v) : NULL;
    v = getenv("AWS_SECRET_ACCESS_KEY"); saved_secret = v ? strdup(v) : NULL;
    v = getenv("AWS_SHARED_CREDENTIALS_FILE"); saved_file = v ? strdup(v) : NULL;

    unsetenv("AWS_ACCESS_KEY_ID");
    unsetenv("AWS_SECRET_ACCESS_KEY");
    setenv("AWS_SHARED_CREDENTIALS_FILE", cred_path, 1);

    s3_url(url, sizeof(url), "creds_file_probe.bin");
    fp = hopen(url, "w");
    if (fp) {
        if (hwrite(fp, data, sizeof(data) - 1) != (ssize_t)(sizeof(data) - 1)) {
            hclose_abruptly(fp);
            fp = NULL;
        } else if (hclose(fp) != 0) {
            fp = NULL;
        }
    }

    got = 0;
    if (fp) {
        fp = hopen(url, "r");
        if (fp) {
            while (got < sizeof(rbuf) &&
                   (n = hread(fp, rbuf + got, sizeof(rbuf) - got)) > 0)
                got += n;
            hclose_abruptly(fp);
        }
    }

    if (saved_id) { setenv("AWS_ACCESS_KEY_ID", saved_id, 1); free(saved_id); }
    if (saved_secret) { setenv("AWS_SECRET_ACCESS_KEY", saved_secret, 1); free(saved_secret); }
    if (saved_file) { setenv("AWS_SHARED_CREDENTIALS_FILE", saved_file, 1); free(saved_file); }
    else unsetenv("AWS_SHARED_CREDENTIALS_FILE");
    unlink(cred_path);

    if (!fp)
        FAIL(name, "roundtrip via credentials file failed: %s", strerror(errno));
    else if (got != sizeof(data) - 1 || memcmp(rbuf, data, got) != 0)
        FAIL(name, "data mismatch reading back via credentials file");
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 15. CLI tools end-to-end (htsfile, tabix) against S3                    */
/* ---------------------------------------------------------------------- */

// Runs cmd via popen, captures combined-ish stdout into buf, returns the
// process exit status (as from pclose) or -1 on a popen-level failure.
static int run_cli(const char *cmd, char *buf, size_t bufsz)
{
    FILE *p = popen(cmd, "r");
    size_t got = 0;
    int status;

    if (!p)
        return -1;

    if (buf && bufsz) {
        got = fread(buf, 1, bufsz - 1, p);
        buf[got] = '\0';
    }
    status = pclose(p);
    return status;
}

static void test_cli_htsfile(void)
{
    const char *name = "CLI: htsfile identifies an S3 BAM";
    char url[512], cmd[600], out[512];
    htsFile *w;
    sam_hdr_t *hdr;
    bam1_t *rec;
    int ret, status;

    s3_url(url, sizeof(url), "cli/probe.bam");

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed"); sam_hdr_destroy(hdr); return;
    }
    w = hts_open(url, "wb");
    if (!w) { FAIL(name, "setup: hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (sam_hdr_write(w, hdr) < 0) { FAIL(name, "setup: sam_hdr_write failed"); hts_close(w); sam_hdr_destroy(hdr); return; }
    rec = bam_init1();
    ret = bam_set1(rec, 4, "read1", 0, 0, 10, 60, 1, cigar_4m, -1, -1, 0, 4, "ACGT", "IIII", 0);
    if (ret < 0 || sam_write1(w, hdr, rec) < 0) {
        FAIL(name, "setup: sam_write1 failed");
        bam_destroy1(rec); hts_close(w); sam_hdr_destroy(hdr);
        return;
    }
    bam_destroy1(rec);
    if (hts_close(w) < 0) { FAIL(name, "setup: hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    snprintf(cmd, sizeof(cmd), "./htsfile '%s' 2>&1", url);
    status = run_cli(cmd, out, sizeof(out));

    if (status < 0)
        FAIL(name, "popen failed: %s", strerror(errno));
    else if (status != 0)
        FAIL(name, "htsfile exited with status %d, output: %s", status, out);
    else if (!strstr(out, "BAM"))
        FAIL(name, "htsfile output did not mention BAM: %s", out);
    else
        PASS(name);
}

static void test_cli_tabix(void)
{
    const char *name = "CLI: tabix queries an S3 BED";
    char url[512], cmd[600], out[4096];
    static const char bed[] =
        "chr1\t10\t20\tfeatureA\n"
        "chr1\t30\t40\tfeatureB\n";
    BGZF *bgzf;
    int status;

    s3_url(url, sizeof(url), "cli/probe.bed.gz");

    bgzf = bgzf_open(url, "w");
    if (!bgzf) { FAIL(name, "setup: bgzf_open(w) failed: %s", strerror(errno)); return; }
    if (bgzf_write(bgzf, bed, sizeof(bed) - 1) != (ssize_t)(sizeof(bed) - 1)) {
        FAIL(name, "setup: bgzf_write failed");
        bgzf_close(bgzf);
        return;
    }
    if (bgzf_close(bgzf) != 0) { FAIL(name, "setup: bgzf_close failed"); return; }

    if (tbx_index_build3(url, NULL, 0, 0, &tbx_conf_bed) < 0) {
        FAIL(name, "setup: tbx_index_build3 failed: %s", strerror(errno));
        return;
    }

    // -D: don't cache the downloaded index locally (the tabix CLI defaults
    // to downloading+persisting a local copy for any remote main file, per
    // args.download_index in tabix.c -- real, but not what this test is
    // about, and we don't want to litter the CWD).
    snprintf(cmd, sizeof(cmd), "./tabix -D '%s' chr1:10-20 2>&1", url);
    status = run_cli(cmd, out, sizeof(out));

    if (status < 0)
        FAIL(name, "popen failed: %s", strerror(errno));
    else if (status != 0)
        FAIL(name, "tabix exited with status %d, output: %s", status, out);
    else if (!strstr(out, "featureA"))
        FAIL(name, "tabix output missing featureA: %s", out);
    else if (strstr(out, "featureB"))
        FAIL(name, "tabix output unexpectedly included featureB (out-of-range): %s", out);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 16. Threaded compression writing to an S3 sink                          */
/* ---------------------------------------------------------------------- */

static void test_threaded_write(void)
{
    const char *name = "threaded BGZF compression writing to S3";
    char url[512];
    htsFile *out, *in;
    sam_hdr_t *hdr, *hdr2;
    bam1_t *rec;
    int i, count, ret;
    char qname[16];

    s3_url(url, sizeof(url), "threaded.bam");

    hdr = sam_hdr_init();
    if (!hdr) { FAIL(name, "sam_hdr_init failed"); return; }
    if (sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "10000", NULL) < 0) {
        FAIL(name, "sam_hdr_add_line failed");
        sam_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wb");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); sam_hdr_destroy(hdr); return; }
    if (hts_set_threads(out, 4) < 0) {
        FAIL(name, "hts_set_threads failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }
    if (sam_hdr_write(out, hdr) < 0) {
        FAIL(name, "sam_hdr_write failed");
        hts_close(out); sam_hdr_destroy(hdr);
        return;
    }

    rec = bam_init1();
    for (i = 0; i < 200; i++) {
        snprintf(qname, sizeof(qname), "read%d", i);
        ret = bam_set1(rec, strlen(qname), qname, 0, 0, i * 10, 60, 1, cigar_4m,
                        -1, -1, 0, 4, "ACGT", "IIII", 0);
        if (ret < 0 || sam_write1(out, hdr, rec) < 0) {
            FAIL(name, "sam_write1 failed at record %d", i);
            bam_destroy1(rec); hts_close(out); sam_hdr_destroy(hdr);
            return;
        }
    }
    bam_destroy1(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); sam_hdr_destroy(hdr); return; }
    sam_hdr_destroy(hdr);

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    if (hts_set_threads(in, 4) < 0) {
        FAIL(name, "hts_set_threads(r) failed");
        hts_close(in);
        return;
    }
    hdr2 = sam_hdr_read(in);
    if (!hdr2) { FAIL(name, "sam_hdr_read failed"); hts_close(in); return; }

    rec = bam_init1();
    count = 0;
    while ((ret = sam_read1(in, hdr2, rec)) >= 0)
        count++;
    bam_destroy1(rec);
    sam_hdr_destroy(hdr2);
    hts_close(in);

    if (ret < -1)
        FAIL(name, "sam_read1 returned error %d", ret);
    else if (count != 200)
        FAIL(name, "expected 200 records, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 17. VCF text (not BCF) + tabix .tbi                                     */
/* ---------------------------------------------------------------------- */

static void test_vcf_text_tabix_roundtrip(void)
{
    const char *name = "VCF text roundtrip + tabix .tbi index";
    char url[512];
    htsFile *out, *in;
    bcf_hdr_t *hdr, *hdr2;
    bcf1_t *rec;
    tbx_t *tbx;
    hts_itr_t *itr;
    kstring_t str = KS_INITIALIZE;
    int count, ret;

    s3_url(url, sizeof(url), "test.vcf.gz");

    hdr = bcf_hdr_init("w");
    if (!hdr) { FAIL(name, "bcf_hdr_init failed"); return; }
    if (bcf_hdr_append(hdr, "##contig=<ID=chr1,length=1000>") < 0) {
        FAIL(name, "bcf_hdr_append failed");
        bcf_hdr_destroy(hdr);
        return;
    }

    out = hts_open(url, "wz");
    if (!out) { FAIL(name, "hts_open(w) failed: %s", strerror(errno)); bcf_hdr_destroy(hdr); return; }
    if (bcf_hdr_write(out, hdr) < 0) {
        FAIL(name, "bcf_hdr_write failed");
        hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }

    rec = bcf_init();
    rec->rid = bcf_hdr_name2id(hdr, "chr1");
    rec->pos = 99;
    rec->qual = 50;
    if (bcf_update_alleles_str(hdr, rec, "A,G") < 0 || bcf_write1(out, hdr, rec) < 0) {
        FAIL(name, "bcf_write1(rec1) failed");
        bcf_destroy(rec); hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }
    bcf_clear1(rec);
    rec->rid = bcf_hdr_name2id(hdr, "chr1");
    rec->pos = 199;
    rec->qual = 30;
    if (bcf_update_alleles_str(hdr, rec, "C,T") < 0 || bcf_write1(out, hdr, rec) < 0) {
        FAIL(name, "bcf_write1(rec2) failed");
        bcf_destroy(rec); hts_close(out); bcf_hdr_destroy(hdr);
        return;
    }
    bcf_destroy(rec);
    if (hts_close(out) < 0) { FAIL(name, "hts_close(w) failed"); bcf_hdr_destroy(hdr); return; }
    bcf_hdr_destroy(hdr);

    if (tbx_index_build3(url, NULL, 0, 0, &tbx_conf_vcf) < 0) {
        FAIL(name, "tbx_index_build3 failed: %s", strerror(errno));
        return;
    }

    in = hts_open(url, "r");
    if (!in) { FAIL(name, "hts_open(r) failed: %s", strerror(errno)); return; }
    hdr2 = bcf_hdr_read(in);
    if (!hdr2) { FAIL(name, "bcf_hdr_read failed"); hts_close(in); return; }
    tbx = tbx_index_load3(url, NULL, 0);
    if (!tbx) {
        FAIL(name, "tbx_index_load3 failed: %s", strerror(errno));
        bcf_hdr_destroy(hdr2); hts_close(in);
        return;
    }

    itr = tbx_itr_queryi(tbx, tbx_name2id(tbx, "chr1"), 0, 1000);
    count = 0;
    while ((ret = tbx_itr_next(in, tbx, itr, &str)) >= 0)
        count++;

    free(str.s);
    hts_itr_destroy(itr);
    tbx_destroy(tbx);
    bcf_hdr_destroy(hdr2);
    hts_close(in);

    if (count != 2)
        FAIL(name, "expected 2 records via tabix query on VCF text, got %d", count);
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 18. Key names requiring URL-escaping                                    */
/* ---------------------------------------------------------------------- */

static void test_key_url_encoding(void)
{
    const char *name = "key names requiring URL-escaping";
    char url[512];
    hFILE *fp;
    static const char data[] = "escaped-key-probe";
    unsigned char rbuf[64];
    ssize_t n;
    size_t got;

    // Space, '+', '#' and '%' all fall outside escape_path()'s unescaped
    // set (hfile_s3.c:425) and must be percent-encoded on every request.
    s3_url(url, sizeof(url), "weird key name/has space+plus#hash%percent.bin");

    fp = hopen(url, "w");
    if (!fp) { FAIL(name, "hopen(w) failed: %s", strerror(errno)); return; }
    if (hwrite(fp, data, sizeof(data) - 1) != (ssize_t)(sizeof(data) - 1)) {
        FAIL(name, "hwrite failed: %s", strerror(errno));
        hclose_abruptly(fp);
        return;
    }
    if (hclose(fp) != 0) { FAIL(name, "hclose(w) failed: %s", strerror(errno)); return; }

    fp = hopen(url, "r");
    if (!fp) { FAIL(name, "hopen(r) failed: %s", strerror(errno)); return; }
    got = 0;
    while (got < sizeof(rbuf) &&
           (n = hread(fp, rbuf + got, sizeof(rbuf) - got)) > 0)
        got += n;
    hclose_abruptly(fp);

    if (got != sizeof(data) - 1)
        FAIL(name, "size mismatch: got %zu, expected %zu", got, sizeof(data) - 1);
    else if (memcmp(rbuf, data, got) != 0)
        FAIL(name, "data mismatch");
    else
        PASS(name);
}

/* ---------------------------------------------------------------------- */
/* 19. HTS_S3_PART_SIZE / HTS_S3_READ_PART_SIZE overrides                  */
/* ---------------------------------------------------------------------- */

static void test_part_size_override(void)
{
    const char *name = "HTS_S3_PART_SIZE / HTS_S3_READ_PART_SIZE overrides";
    char url[512];
    const size_t total = 8 * 1024 * 1024; // 8 MiB: > the 6 MiB override below
    unsigned char *wbuf, *rbuf;
    hFILE *fp;
    ssize_t n;
    size_t got;

    s3_url(url, sizeof(url), "part_size_override.bin");

    wbuf = malloc(total);
    rbuf = malloc(total);
    if (!wbuf || !rbuf) {
        FAIL(name, "malloc failed");
        free(wbuf); free(rbuf);
        return;
    }
    generate_bytes(wbuf, total, 61);

    // HTS_S3_PART_SIZE is in MiB and can only raise the 5 MiB floor
    // (hfile_s3.c:2256-2260); HTS_S3_READ_PART_SIZE is also in MiB and sets
    // the ranged-GET chunk size outright (hfile_s3.c:2386-2389).
    setenv("HTS_S3_PART_SIZE", "6", 1);
    setenv("HTS_S3_READ_PART_SIZE", "2", 1);

    fp = hopen(url, "w");
    if (fp) {
        if (hwrite(fp, wbuf, total) != (ssize_t)total) {
            hclose_abruptly(fp);
            fp = NULL;
        } else if (hclose(fp) != 0) {
            fp = NULL;
        }
    }

    got = 0;
    if (fp) {
        fp = hopen(url, "r");
        if (fp) {
            while (got < total && (n = hread(fp, rbuf + got, total - got)) > 0)
                got += n;
            hclose_abruptly(fp);
        }
    }

    unsetenv("HTS_S3_PART_SIZE");
    unsetenv("HTS_S3_READ_PART_SIZE");

    if (!fp)
        FAIL(name, "roundtrip with overridden part sizes failed: %s", strerror(errno));
    else if (got != total)
        FAIL(name, "size mismatch: got %zu, expected %zu", got, total);
    else if (memcmp(wbuf, rbuf, total) != 0)
        FAIL(name, "data mismatch");
    else
        PASS(name);

    free(wbuf);
    free(rbuf);
}

int main(void)
{
    const char *host = getenv("HTS_S3_HOST");
    const char *b = getenv("HTSLIB_TEST_S3_BUCKET");
    pid_t pid = getpid();

    if (!host || !*host) {
        fprintf(stderr, "HTS_S3_HOST not set, skipping S3 integration tests "
                        "(see test/docker/minio/ and 'make check-s3')\n");
        return 0;
    }

    snprintf(bucket, sizeof(bucket), "%s", (b && *b) ? b : "htslib-test");
    snprintf(run_id, sizeof(run_id), "test-hfile-s3-%ld-%ld",
             (long)pid, (long)time(NULL));

    fprintf(stderr, "test_hfile_s3: bucket=%s host=%s run_id=%s\n",
            bucket, host, run_id);

    test_hfile_roundtrip();
    test_bgzf_roundtrip();
    test_sam_bam_roundtrip();
    test_bcf_roundtrip();
    test_tabix_roundtrip();
    test_faidx_roundtrip();
    test_cram_roundtrip();
    test_cram_explicit_s3_reference();
    test_index_local_fallback_bypass();
    test_multipart_large_file();
    test_unhappy_read_nonexistent_key();
    test_unhappy_read_nonexistent_bucket();
    test_unhappy_wrong_credentials();
    test_unhappy_connection_refused();
    test_unhappy_hts_open_nonexistent();
    test_cram_index_roundtrip();
    test_sigv2_read();
    test_credentials_file_fallback();
    test_cli_htsfile();
    test_cli_tabix();
    test_threaded_write();
    test_vcf_text_tabix_roundtrip();
    test_key_url_encoding();
    test_part_size_override();

    if (failures > 0) {
        fprintf(stderr, "%d test(s) FAILED\n", failures);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "All tests passed.\n");
    return EXIT_SUCCESS;
}
