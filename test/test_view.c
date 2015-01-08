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

typedef struct hts_opt {
    enum cram_option opt;
    union {
        int i;
        char *s;
    } val;
    struct hts_opt *next;
} hts_opt;

/*
 * Parses arg and appends it to the option list.
 * Returns 0 on success;
 *        -1 on failure.
 */
int add_option(hts_opt **opts, char *arg) {
    hts_opt *o, *t;
    char *cp;

    if (!(cp = strchr(arg, '=')))
        cp = "1"; // assume boolean
    else
        *cp++ = 0;

    if (!(o =  malloc(sizeof(*o))))
        return -1;

    if (strcmp(arg, "DECODE_MD") == 0)
        o->opt = CRAM_OPT_DECODE_MD, o->val.i = atoi(cp);
    else if (strcmp(arg, "VERBOSITY") == 0)
        o->opt = CRAM_OPT_VERBOSITY, o->val.i = atoi(cp);
    else if (strcmp(arg, "SEQS_PER_SLICE") == 0)
        o->opt = CRAM_OPT_SEQS_PER_SLICE, o->val.i = atoi(cp);
    else if (strcmp(arg, "SLICES_PER_CONTAINER") == 0)
        o->opt = CRAM_OPT_SLICES_PER_CONTAINER, o->val.i = atoi(cp);
    else if (strcmp(arg, "EMBED_REF") == 0)
        o->opt = CRAM_OPT_EMBED_REF, o->val.i = atoi(cp);
    else if (strcmp(arg, "NO_REF") == 0)
        o->opt = CRAM_OPT_NO_REF, o->val.i = atoi(cp);
    else if (strcmp(arg, "IGNORE_MD5") == 0)
        o->opt = CRAM_OPT_IGNORE_MD5, o->val.i = atoi(cp);
    else if (strcmp(arg, "USE_BZIP2") == 0)
        o->opt = CRAM_OPT_USE_BZIP2, o->val.i = atoi(cp);
    else if (strcmp(arg, "USE_RANS") == 0)
        o->opt = CRAM_OPT_USE_RANS, o->val.i = atoi(cp);
    else if (strcmp(arg, "USE_LZMA") == 0)
        o->opt = CRAM_OPT_USE_LZMA, o->val.i = atoi(cp);
    else if (strcmp(arg, "REFERENCE") == 0)
        o->opt = CRAM_OPT_REFERENCE, o->val.s = cp;
    else if (strcmp(arg, "VERSION") == 0)
        o->opt = CRAM_OPT_VERSION, o->val.s =cp;
    else if (strcmp(arg, "MULTI_SEQ_PER_SLICE") == 0)
        o->opt = CRAM_OPT_MULTI_SEQ_PER_SLICE, o->val.i = atoi(cp);
    else if (strcmp(arg, "NTHREADS") == 0)
        o->opt = CRAM_OPT_NTHREADS, o->val.i = atoi(cp);
    else if (strcmp(arg, "REQUIRED_FIELDS") == 0)
        o->opt = CRAM_OPT_REQUIRED_FIELDS, o->val.i = strtol(cp, NULL, 0);
    else {
        fprintf(stderr, "Unknown option '%s'\n", arg);
        free(o);
        return -1;
    }

    o->next = NULL;

    if (*opts) {
        t = *opts;
        while (t->next)
            t = t->next;
        t->next = o;
    } else {
        *opts = o;
    }

    return 0;
}

int main(int argc, char *argv[])
{
    samFile *in;
    char *fn_ref = 0;
    int flag = 0, c, clevel = -1, ignore_sam_err = 0;
    char moder[8];
    bam_hdr_t *h;
    bam1_t *b;
    htsFile *out;
    char modew[8];
    int r = 0, exit_code = 0;
    hts_opt *in_opts = NULL, *out_opts = NULL, *last = NULL;

    while ((c = getopt(argc, argv, "IbDCSl:t:i:o:")) >= 0) {
        switch (c) {
        case 'S': flag |= 1; break;
        case 'b': flag |= 2; break;
        case 'D': flag |= 4; break;
        case 'C': flag |= 8; break;
        case 'l': clevel = atoi(optarg); flag |= 2; break;
        case 't': fn_ref = optarg; break;
        case 'I': ignore_sam_err = 1; break;
        case 'i': if (add_option(&in_opts,  optarg)) return 1; break;
        case 'o': if (add_option(&out_opts, optarg)) return 1; break;
        }
    }
    if (argc == optind) {
        fprintf(stderr, "Usage: samview [-bSCSI] [-l level] [-o option=value] <in.bam>|<in.sam>|<in.cram> [region]\n");
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
    h->ignore_sam_err = ignore_sam_err;
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
        // Parse input header and use for CRAM output
        out->fp.cram->header = sam_hdr_parse_(h->text, h->l_text);

        // Create CRAM references arrays
        if (fn_ref)
            cram_set_option(out->fp.cram, CRAM_OPT_REFERENCE, fn_ref);
        else
            // Attempt to fill out a cram->refs[] array from @SQ headers
            cram_set_option(out->fp.cram, CRAM_OPT_REFERENCE, NULL);
    }

    // Process any options; currently cram only.
    for (; in_opts;  in_opts = (last=in_opts)->next, free(last)) {
        hts_set_opt(in,  in_opts->opt,  in_opts->val);
        if (in_opts->opt == CRAM_OPT_REFERENCE)
            hts_set_opt(out,  in_opts->opt,  in_opts->val);
    }
    for (; out_opts;  out_opts = (last=out_opts)->next, free(last))
        hts_set_opt(out, out_opts->opt,  out_opts->val);

    sam_hdr_write(out, h);
    if (optind + 1 < argc && !(flag&1)) { // BAM input and has a region
        int i;
        hts_idx_t *idx;
        if ((idx = bam_index_load(argv[optind])) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            return 1;
        }
        for (i = optind + 1; i < argc; ++i) {
            hts_itr_t *iter;
            if ((iter = bam_itr_querys(idx, h, argv[i])) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
                continue;
            }
            while ((r = bam_itr_next(in, iter, b)) >= 0) {
                if (sam_write1(out, h, b) < 0) {
                    fprintf(stderr, "Error writing output.\n");
                    exit_code = 1;
                    break;
                }
            }
            hts_itr_destroy(iter);
        }
        hts_idx_destroy(idx);
    } else while ((r = sam_read1(in, h, b)) >= 0) {
        if (sam_write1(out, h, b) < 0) {
            fprintf(stderr, "Error writing output.\n");
            exit_code = 1;
            break;
        }
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
