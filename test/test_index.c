/*  test/test_index.c -- simple tool to build an index, for the test harness.

    Copyright (C) 2018 Genome Research Ltd.

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
#include <stdio.h>
#include <unistd.h>

#include "../htslib/sam.h"
#include "../htslib/vcf.h"

void HTS_NORETURN usage(FILE *fp) {
    fprintf(fp, "Usage: test_index [opts] in.{sam.gz,bam,cram}|in.{vcf.gz,bcf}\n\n");
    fprintf(fp, "  -b       Use BAI index (BAM, SAM)\n");
    fprintf(fp, "  -c       Use CSI index (BAM, SAM, VCF, BCF)\n");
    fprintf(fp, "  -t       Use TBI index (VCF) \n");
    fprintf(fp, "  -m bits  Adjust min_shift; implies CSI\n");
    fprintf(fp, "\nThe default index format is CSI for sam/bam/vcf/bcf and CRAI for crams\n");
    exit(fp == stderr ? 1 : 0);
}

int main(int argc, char **argv) {
    int c, min_shift = 14;

    while ((c = getopt(argc, argv, "bctm:")) >= 0) {
        switch (c) {
        case 't': case 'b': min_shift = 0; break;
        case 'c': min_shift = 14; break;
        case 'm': min_shift = atoi(optarg); break;
        case 'h': usage(stdout);
        default:  usage(stderr);
        }
    }

    if (optind >= argc) usage(stderr);

    htsFile *in = hts_open(argv[optind], "r");
    if (!in) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        exit(1);
    }

    int ret;
    if (in->format.format == sam ||
        in->format.format == bam ||
        in->format.format == cram) {
        ret = sam_index_build(argv[optind], min_shift);
    } else {
        ret = bcf_index_build(argv[optind], min_shift);
    }

    if (ret < 0) {
        fprintf(stderr, "Failed to build index for \"%s\"\n", argv[optind]);
        exit(1);
    }

    if (hts_close(in) < 0) {
        fprintf(stderr, "Error closing \"%s\"\n", argv[optind]);
        exit(1);
    }

    return 0;
}
