/*  fuzz_expr.c -- Fuzz driver for hts_filter_init.

    Copyright (C) 2026 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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
#include "../../htslib/hts_expr.h"

#ifndef FUZZ_INPUT
#  define FUZZ_INPUT "expr.sam"
#endif

const char *SAM_DAT="data:,"
"@SQ\tSN:CHROMOSOME_I\tLN:1009800\n"
"SRR065390.14978392\t16\tCHROMOSOME_I\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tAS:i:-18\tXS:i:-18\tYT:Z:UU\tRG:Z:ID\tA!:A:!\tAc:A:c\tAC:A:C\tI0:i:0\tI1:i:1\tI2:i:127\tI3:i:128\tI4:i:255\tI5:i:256\tI6:i:32767\tI7:i:32768\tI8:i:65535\tI9:i:65536\tIA:i:2147483647\ti1:i:-1\ti2:i:-127\ti3:i:-128\ti4:i:-255\ti5:i:-256\ti6:i:-32767\ti7:i:-32768\ti8:i:-65535\ti9:i:-65536\tiA:i:-2147483647\tiB:i:-2147483648\tF0:f:-1\tF1:f:0\tF2:f:1\tF3:f:9.9e-19\tF4:f:-9.9e-19\tF5:f:9.9e+19\tF6:f:-9.9e+19\tH0:H:AA\tH1:H:dead00beef\tZ0:Z:space space\tZn:Z:Hn:H:\n";

static int opened = 0;
static bam1_t *b = NULL;
static sam_hdr_t *hdr = NULL;
static samFile *in = NULL;

// Possible race condition, but don't particularly care for the fuzzer.
// We could also construct a bam using bam_set1().
static void init_bam(void) {
    puts("INIT_BAM START");
    opened = 1;
    if (!(in = sam_open(SAM_DAT, "r")))
        abort();

    if (!(b = bam_init1()))
        abort();

    sam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == NULL)
        abort();

    if (sam_read1(in, hdr, b) < 0)
        abort();

    puts("INIT_BAM END");
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    // Create the filter
    char expr[8192];
    size_t len = strnlen((char *)data, size);
    if (len > 8191)
        len = 8191;
    memcpy(expr, data, len);
    expr[len] = 0;

    hts_filter_t *filt = hts_filter_init(expr);
    if (!filt)
        return 0;

    if (!opened)
        init_bam();

    int r = sam_passes_filter(hdr, b, filt);
    printf("r=%d\n", r);
    hts_filter_free(filt);

    return 0;
}
