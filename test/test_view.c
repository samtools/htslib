/*  test/test_view.c -- simple view tool, purely for use in a test harness.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "cram/cram.h"

#include "htslib/sam.h"

int main(int argc, char *argv[])
{
    samFile *in;
    char *fn_ref = 0;
    int flag = 0, c, clevel = -1, ignore_sam_err = 0;
    char moder[8];
    bam_hdr_t *h;
    bam1_t *b;
    htsFile *out;
    char modew[800];
    int r = 0, exit_code = 0;
    hts_opt *in_opts = NULL, *out_opts = NULL;
    int nreads = 0;
    int extra_hdr_nuls = 0;
    int benchmark = 0;

    while ((c = getopt(argc, argv, "IbDCSl:t:i:o:N:BZ:")) >= 0) {
        switch (c) {
        case 'S': flag |= 1; break;
        case 'b': flag |= 2; break;
        case 'D': flag |= 4; break;
        case 'C': flag |= 8; break;
        case 'B': benchmark = 1; break;
        case 'l': clevel = atoi(optarg); flag |= 2; break;
        case 't': fn_ref = optarg; break;
        case 'I': ignore_sam_err = 1; break;
        case 'i': if (hts_opt_add(&in_opts,  optarg)) return 1; break;
        case 'o': if (hts_opt_add(&out_opts, optarg)) return 1; break;
        case 'N': nreads = atoi(optarg); break;
        case 'Z': extra_hdr_nuls = atoi(optarg); break;
        }
    }
    if (argc == optind) {
        fprintf(stderr, "Usage: samview [-bSCSIB] [-N num_reads] [-l level] [-o option=value] [-Z hdr_nuls] <in.bam>|<in.sam>|<in.cram> [region]\n");
        return 1;
    }
    strcpy(moder, "r");
    if (flag&4) strcat(moder, "c");
    else if ((flag&1) == 0) strcat(moder, "b");

    in = sam_open(argv[optind], moder);
    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    h = sam_hdr_read(in);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    h->ignore_sam_err = ignore_sam_err;
    if (extra_hdr_nuls) {
        char *new_text = realloc(h->text, h->l_text + extra_hdr_nuls);
        if (new_text == NULL) {
            fprintf(stderr, "Error reallocing header text\n");
            return EXIT_FAILURE;
        }
        h->text = new_text;
        memset(&h->text[h->l_text], 0, extra_hdr_nuls);
        h->l_text += extra_hdr_nuls;
    }

    b = bam_init1();

    strcpy(modew, "w");
    if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
    if (flag&8) strcat(modew, "c");
    else if (flag&2) strcat(modew, "b");
    out = hts_open("-", modew);
    if (out == NULL) {
        fprintf(stderr, "Error opening standard output\n");
        return EXIT_FAILURE;
    }

    /* CRAM output */
    if (flag & 8) {
        int ret;

        // Parse input header and use for CRAM output
        out->fp.cram->header = sam_hdr_parse_(h->text, h->l_text);

        // Create CRAM references arrays
        if (fn_ref)
            ret = cram_set_option(out->fp.cram, CRAM_OPT_REFERENCE, fn_ref);
        else
            // Attempt to fill out a cram->refs[] array from @SQ headers
            ret = cram_set_option(out->fp.cram, CRAM_OPT_REFERENCE, NULL);

        if (ret != 0)
            return EXIT_FAILURE;
    }

    // Process any options; currently cram only.
    if (hts_opt_apply(in, in_opts))
        return EXIT_FAILURE;
    hts_opt_free(in_opts);

    if (hts_opt_apply(out, out_opts))
        return EXIT_FAILURE;
    hts_opt_free(out_opts);

    if (!benchmark)
        sam_hdr_write(out, h);
    if (optind + 1 < argc && !(flag&1)) { // BAM input and has a region
        int i;
        hts_idx_t *idx;
        if ((idx = sam_index_load(in, argv[optind])) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            return 1;
        }
        for (i = optind + 1; i < argc; ++i) {
            hts_itr_t *iter;
            if ((iter = sam_itr_querys(idx, h, argv[i])) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
                continue;
            }
            while ((r = sam_itr_next(in, iter, b)) >= 0) {
                if (!benchmark && sam_write1(out, h, b) < 0) {
                    fprintf(stderr, "Error writing output.\n");
                    exit_code = 1;
                    break;
                }
                if (nreads && --nreads == 0)
                    break;
            }
            hts_itr_destroy(iter);
        }
        hts_idx_destroy(idx);
    } else while ((r = sam_read1(in, h, b)) >= 0) {
        if (!benchmark && sam_write1(out, h, b) < 0) {
            fprintf(stderr, "Error writing output.\n");
            exit_code = 1;
            break;
        }
        if (nreads && --nreads == 0)
            break;
    }

    if (r < -1) {
        fprintf(stderr, "Error parsing input.\n");
        exit_code = 1;
    }

    r = sam_close(out);
    if (r < 0) {
        fprintf(stderr, "Error closing output.\n");
        exit_code = 1;
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);

    r = sam_close(in);
    if (r < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = 1;
    }

    return exit_code;
}
