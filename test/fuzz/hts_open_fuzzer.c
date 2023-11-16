/*  test/fuzz/hts_open_fuzzer.c -- Fuzz driver for hts_open.

    Copyright (C) 2018 Google LLC.
    Copyright (C) 2019-2020 Genome Research Ltd.

    Author: Markus Kusano <kusano@google.com>

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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../htslib/hfile.h"
#include "../../htslib/hts.h"
#include "../../htslib/sam.h"
#include "../../htslib/vcf.h"

static void hts_close_or_abort(htsFile* file) {
    if (hts_close(file) != 0) {
        abort();
    }
}

static void view_sam(htsFile *in) {
    if (!in) {
        return;
    }
    samFile *outS = sam_open("/dev/null", "w");
    samFile *outB = sam_open("/dev/null", "wb");
    samFile *outC = sam_open("/dev/null", "wc");
    if (!outS || !outB || !outC) {
        abort();
    }
    // Not critical if this doesn't work, but can test more if
    // we're in the right location.
    if (hts_set_fai_filename(outC, "../c2.fa") < 0) {
        static int warned = 0;
        if (!warned) {
            warned = 1;
            fprintf(stderr, "Warning couldn't find the c2.fa file\n");
        }
    }
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        hts_close_or_abort(outS);
        hts_close_or_abort(outB);
        hts_close(outC);
        return;
    }

    // This will force the header to be parsed.
    (void) sam_hdr_count_lines(hdr, "SQ");

    if (sam_hdr_write(outS, hdr) != 0 ||
        sam_hdr_write(outB, hdr) != 0 ||
        sam_hdr_write(outC, hdr) != 0) {
        sam_hdr_destroy(hdr);
        hts_close_or_abort(outS);
        hts_close_or_abort(outB);
        hts_close(outC);
        return;
    }
    bam1_t *b = bam_init1();
    if (b == NULL) {
        sam_hdr_destroy(hdr);
        hts_close_or_abort(outS);
        hts_close_or_abort(outB);
        hts_close(outC);
        return;
    }
    while (sam_read1(in, hdr, b) >= 0) {
        if (sam_write1(outS, hdr, b) < 0 ||
            sam_write1(outB, hdr, b) < 0 ||
            sam_write1(outC, hdr, b) < 0) {
            break;
        }
    }
    bam_destroy1(b);

    sam_hdr_destroy(hdr);
    hts_close_or_abort(outS);
    hts_close_or_abort(outB);
    hts_close(outC);
}

static void view_vcf(htsFile *in) {
    if (!in) {
        return;
    }
    vcfFile *outv = vcf_open("/dev/null", "w");
    vcfFile *outb = vcf_open("/dev/null", "wb");
    if (!outv || !outb) {
        abort();
    }
    bcf_hdr_t *hdr = bcf_hdr_read(in);
    if (hdr == NULL) {
        hts_close_or_abort(outv);
        hts_close_or_abort(outb);
        return;
    }

    if (bcf_hdr_write(outv, hdr) != 0 ||
        bcf_hdr_write(outb, hdr) != 0) {
        bcf_hdr_destroy(hdr);
        hts_close_or_abort(outv);
        hts_close_or_abort(outb);
    }
    bcf1_t *rec = bcf_init();
    if (rec == NULL) {
        bcf_hdr_destroy(hdr);
        hts_close_or_abort(outv);
        hts_close_or_abort(outb);
    }
    while (bcf_read(in, hdr, rec) >= 0) {
        if (bcf_write(outv, hdr, rec) < 0 ||
            bcf_write(outb, hdr, rec) < 0) {
            break;
        }
    }
    bcf_destroy(rec);

    bcf_hdr_destroy(hdr);
    hts_close_or_abort(outv);
    hts_close_or_abort(outb);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    hFILE *memfile;
    uint8_t *copy = malloc(size);
    if (copy == NULL) {
        abort();
    }
    memcpy(copy, data, size);
    // hopen does not take ownership of `copy`, but hts_hopen does.
    memfile = hopen("mem:", "rb:", copy, size);
    if (memfile == NULL) {
        free(copy);
        return 0;
    }

    htsFile *ht_file = hts_hopen(memfile, "data", "rb");
    if (ht_file == NULL) {
        if (hclose(memfile) != 0) {
            abort();
        }
        return 0;
    }
    switch (ht_file->format.category) {
        case sequence_data:
            view_sam(ht_file);
            break;
        case variant_data:
            view_vcf(ht_file);
            break;
        default:
            break;
    }
    hts_close(ht_file);
    return 0;
}
