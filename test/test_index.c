/*  Derived from samtools bam_index.c -- index subcommand.

    Copyright (C) 2008-2011, 2013, 2014, 2015 Genome Research Ltd.
    Portions copyright (C) 2010 Broad Institute.
    Portions copyright (C) 2013 Peter Cock, The James Hutton Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define BAM_LIDX_SHIFT 14

static void index_usage(FILE *fp) {
    char *msg =
        "Usage: samtools index [-bc] [-m INT] <in.[bam|cram]> [out.index]\n"
        "Options:\n"
        "  -b      Generate BAI-format index for BAM files [default]\n"
        "  -c      Generate CSI-format index for BAM files\n"
        "  -m INT  Set minimum interval size for CSI indices to 2^INT [%d]\n";
    fprintf(fp, msg, BAM_LIDX_SHIFT);
}

int main(int argc, char *argv[]) {
    int csi = 0;
    int min_shift = BAM_LIDX_SHIFT;
    int c;

    while ((c = getopt(argc, argv, "bcm:")) >= 0) {
        switch (c) {
        case 'b': csi = 0; break;
        case 'c': csi = 1; break;
        case 'm': csi = 1; min_shift = atoi(optarg); break;
        default:
            index_usage(stderr);
            return 1;
        }
    }

    if (optind == argc) {
        index_usage(stdout);
        return 1;
    }
    
    int ret = bam_index_build(argv[optind], csi ? min_shift : 0);
    if (ret < 0) {
        return 1;
    }

    return ret;
}
