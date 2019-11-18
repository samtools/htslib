/*  test/test_view.c -- simple view tool, purely for use in a test harness.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013-2019 Genome Research Ltd.

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

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>

#include "cram/cram.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts_log.h"

struct opts {
    char *fn_ref;
    int flag;
    int clevel;
    int ignore_sam_err;
    int nreads;
    int extra_hdr_nuls;
    int benchmark;
    int nthreads;
    int multi_reg;
    char *index;
    int min_shift;
};

enum test_op {
    READ_COMPRESSED    = 1,
    WRITE_BINARY_COMP  = 2, // eg bam, bcf
    READ_CRAM          = 4,
    WRITE_CRAM         = 8,
    WRITE_UNCOMPRESSED = 16,
    WRITE_COMPRESSED   = 32, // eg vcf.gz, sam.gz
};

int sam_loop(int argc, char **argv, int optind, struct opts *opts, htsFile *in, htsFile *out) {
    int r = 0;
    sam_hdr_t *h = NULL;
    hts_idx_t *idx = NULL;
    bam1_t *b = NULL;

    h = sam_hdr_read(in);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }
    h->ignore_sam_err = opts->ignore_sam_err;
    if (opts->extra_hdr_nuls > 0) {
        char *new_text = realloc(h->text, h->l_text + opts->extra_hdr_nuls);
        if (new_text == NULL) {
            fprintf(stderr, "Error reallocing header text\n");
            goto fail;
        }
        h->text = new_text;
        memset(&h->text[h->l_text], 0, opts->extra_hdr_nuls);
        h->l_text += opts->extra_hdr_nuls;
    }

    b = bam_init1();
    if (b == NULL) {
        fprintf(stderr, "Out of memory allocating BAM struct\n");
        goto fail;
    }

    /* CRAM output */
    if ((opts->flag & WRITE_CRAM) && opts->fn_ref) {
        // Create CRAM references arrays
        int ret = hts_set_fai_filename(out, opts->fn_ref);

        if (ret != 0)
            goto fail;
    }

    if (!opts->benchmark && sam_hdr_write(out, h) < 0) {
        fprintf(stderr, "Error writing output header.\n");
        goto fail;
    }

    if (opts->index) {
        if (sam_idx_init(out, h, opts->min_shift, opts->index) < 0) {
            fprintf(stderr, "Failed to initialise index\n");
            goto fail;
        }
    }

    if (optind + 1 < argc && !(opts->flag & READ_COMPRESSED)) { // BAM input and has a region
        int i;
        if ((idx = sam_index_load(in, argv[optind])) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
            goto fail;
        }
        if (opts->multi_reg) {
            hts_itr_t *iter = sam_itr_regarray(idx, h, &argv[optind + 1], argc - optind-1);
            if (!iter)
                goto fail;
            while ((r = sam_itr_next(in, iter, b)) >= 0) {
                if (!opts->benchmark && sam_write1(out, h, b) < 0) {
                    fprintf(stderr, "Error writing output.\n");
                    hts_itr_destroy(iter);
                    goto fail;
                }
                if (opts->nreads && --opts->nreads == 0)
                    break;
            }
            hts_itr_destroy(iter);
            if (r < -1) {
                fprintf(stderr, "Error reading input.\n");
                goto fail;
            }
        } else {
            for (i = optind + 1; i < argc; ++i) {
                hts_itr_t *iter;
                if ((iter = sam_itr_querys(idx, h, argv[i])) == 0) {
                    fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
                    goto fail;
                }
                while ((r = sam_itr_next(in, iter, b)) >= 0) {
                    if (!opts->benchmark && sam_write1(out, h, b) < 0) {
                        fprintf(stderr, "Error writing output.\n");
                        hts_itr_destroy(iter);
                        goto fail;
                    }
                    if (opts->nreads && --opts->nreads == 0)
                        break;
                }
                hts_itr_destroy(iter);
                if (r < -1) {
                    fprintf(stderr, "Error reading input.\n");
                    goto fail;
                }
            }
        }
        hts_idx_destroy(idx); idx = NULL;
    } else while ((r = sam_read1(in, h, b)) >= 0) {
        if (!opts->benchmark && sam_write1(out, h, b) < 0) {
            fprintf(stderr, "Error writing output.\n");
            goto fail;
        }
        if (opts->nreads && --opts->nreads == 0)
            break;
    }

    if (r < -1) {
        fprintf(stderr, "Error parsing input.\n");
        goto fail;
    }

    if (opts->index) {
        if (sam_idx_save(out) < 0) {
            fprintf(stderr, "Error saving index\n");
            goto fail;
        }
    }

    bam_destroy1(b);
    sam_hdr_destroy(h);

    return 0;
 fail:
    if (b) bam_destroy1(b);
    if (h) sam_hdr_destroy(h);
    if (idx) hts_idx_destroy(idx);

    return 1;
}

int vcf_loop(int argc, char **argv, int optind, struct opts *opts, htsFile *in, htsFile *out) {
    bcf_hdr_t *h = bcf_hdr_read(in);
    bcf1_t *b = bcf_init1();
    hts_idx_t *idx;
    int i, exit_code = 0, r = 0;

    if (!h)
        return 1;
    if (!b)
        return 1;

    if (!opts->benchmark && bcf_hdr_write(out, h) < 0)
        return 1;

    if (opts->index) {
        if (bcf_idx_init(out, h, opts->min_shift, opts->index) < 0) {
            fprintf(stderr, "Failed to initialise index\n");
            return 1;
        }
    }

    if (optind + 1 < argc) {
        // A series of regions.
        if ((idx = bcf_index_load(argv[optind])) == 0) {
            fprintf(stderr, "[E::%s] fail to load the BVCF index\n", __func__);
            return 1;
        }

        for (i = optind + 1; i < argc; i++) {
            hts_itr_t *iter;
            if ((iter = bcf_itr_querys(idx, h, argv[i])) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
                continue;
            }
            while ((r = bcf_itr_next(in, iter, b)) >= 0) {
                if (!opts->benchmark && bcf_write1(out, h, b) < 0) {
                    fprintf(stderr, "Error writing output.\n");
                    exit_code = 1;
                    break;
                }
                if (opts->nreads && --opts->nreads == 0)
                    break;
            }
            if (r < -1) {
                fprintf(stderr, "Error reading input.\n");
                exit_code = 1;
            }
            hts_itr_destroy(iter);
            if (exit_code != 0) break;
        }

        hts_idx_destroy(idx);

    } else {
        // Whole file
        while ((r = bcf_read1(in, h, b)) >= 0) {
            if (!opts->benchmark && bcf_write1(out, h, b) < 0) {
                fprintf(stderr, "Error writing output.\n");
                exit_code = 1;
                break;
            }
            if (opts->nreads && --opts->nreads == 0)
                break;
        }
        if (r < -1) {
            fprintf(stderr, "Error reading input.\n");
            exit_code = 1;
        }
    }

    if (exit_code == 0 && opts->index) {
        if (bcf_idx_save(out) < 0) {
            fprintf(stderr, "Error saving index\n");
            exit_code = 1;
        }
    }

    bcf_destroy1(b);
    bcf_hdr_destroy(h);
    return exit_code;
}

int main(int argc, char *argv[])
{
    htsFile *in, *out;
    char moder[8];
    char modew[800];
    int c, exit_code = EXIT_SUCCESS;
    hts_opt *in_opts = NULL, *out_opts = NULL;
    char *out_fn = "-";

    struct opts opts;
    opts.fn_ref = NULL;
    opts.flag = 0;
    opts.clevel = -1;
    opts.ignore_sam_err = 0;
    opts.nreads = 0;
    opts.extra_hdr_nuls = 0;
    opts.benchmark = 0;
    opts.nthreads = 0; // shared pool
    opts.multi_reg = 0;
    opts.index = NULL;
    opts.min_shift = 0;

    while ((c = getopt(argc, argv, "DSIt:i:bzCul:o:N:BZ:@:Mx:m:p:v")) >= 0) {
        switch (c) {
        case 'D': opts.flag |= READ_CRAM; break;
        case 'S': opts.flag |= READ_COMPRESSED; break;
        case 'I': opts.ignore_sam_err = 1; break;
        case 't': opts.fn_ref = optarg; break;
        case 'i': if (hts_opt_add(&in_opts, optarg)) return 1; break;
        case 'b': opts.flag |= WRITE_BINARY_COMP; break;
        case 'z': opts.flag |= WRITE_COMPRESSED; break;
        case 'C': opts.flag |= WRITE_CRAM; break;
        case 'u': opts.flag |= WRITE_UNCOMPRESSED; break; // eg u-BAM not SAM
        case 'l': opts.clevel = atoi(optarg); break;
        case 'o': if (hts_opt_add(&out_opts, optarg)) return 1; break;
        case 'N': opts.nreads = atoi(optarg); break;
        case 'B': opts.benchmark = 1; break;
        case 'Z': opts.extra_hdr_nuls = atoi(optarg); break;
        case 'M': opts.multi_reg = 1; break;
        case '@': opts.nthreads = atoi(optarg); break;
        case 'x': opts.index = optarg; break;
        case 'm': opts.min_shift = atoi(optarg); break;
        case 'p': out_fn = optarg; break;
        case 'v': hts_verbose++; break;
        }
    }
    if (argc == optind) {
        fprintf(stderr, "Usage: test_view [-DSI] [-t fn_ref] [-i option=value] [-bC] [-l level] [-o option=value] [-N num_reads] [-B] [-Z hdr_nuls] [-@ num_threads] [-x index_fn] [-m min_shift] [-p out] [-v] <in.bam>|<in.sam>|<in.cram> [region]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-D: read CRAM format (mode 'c')\n");
        fprintf(stderr, "-S: read compressed BCF, BAM, FAI (mode 'b')\n");
        fprintf(stderr, "-I: ignore SAM parsing errors\n");
        fprintf(stderr, "-t: fn_ref: load CRAM references from the specificed fasta file instead of @SQ headers when writing a CRAM file\n");
        fprintf(stderr, "-i: option=value: set an option for CRAM input\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-b: write binary compressed BCF, BAM, FAI (mode 'b')\n");
        fprintf(stderr, "-z: write text compressed VCF.gz, SAM.gz (mode 'z')\n");
        fprintf(stderr, "-C: write CRAM format (mode 'c')\n");
        fprintf(stderr, "-l 0-9: set zlib compression level\n");
        fprintf(stderr, "-o option=value: set an option for CRAM output\n");
        fprintf(stderr, "-N: num_reads: limit the output to the first num_reads reads\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "-B: enable benchmarking\n");
        fprintf(stderr, "-M: use hts_itr_multi iterator\n");
        fprintf(stderr, "-Z hdr_nuls: append specified number of null bytes to the SAM header\n");
        fprintf(stderr, "-@ num_threads: use thread pool with specified number of threads\n\n");
        fprintf(stderr, "-x fn: write index to fn\n");
        fprintf(stderr, "-m min_shift: specifies BAI/CSI bin size; 0 is BAI(BAM) or TBI(VCF), 14 is CSI default\n");
        fprintf(stderr, "-p out_fn: output to out_fn instead of stdout\n");
        fprintf(stderr, "-v: increase verbosity\n");
        fprintf(stderr, "The region list entries should be specified as 'reg:beg-end', with intervals of a region being disjunct and sorted by the starting coordinate.\n");
        return 1;
    }
    strcpy(moder, "r");
    if (opts.flag & READ_CRAM) strcat(moder, "c");
    else if ((opts.flag & READ_COMPRESSED) == 0) strcat(moder, "b");

    in = hts_open(argv[optind], moder);
    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    strcpy(modew, "w");
    if (opts.clevel >= 0 && opts.clevel <= 9) sprintf(modew + 1, "%d", opts.clevel);
    if (opts.flag & WRITE_CRAM) strcat(modew, "c");
    else if (opts.flag & WRITE_BINARY_COMP) strcat(modew, "b");
    else if (opts.flag & WRITE_COMPRESSED) strcat(modew, "z");
    else if (opts.flag & WRITE_UNCOMPRESSED) strcat(modew, "bu");
    out = hts_open(out_fn, modew);
    if (out == NULL) {
        fprintf(stderr, "Error opening standard output\n");
        return EXIT_FAILURE;
    }

    // Process any options; currently cram only.
    if (hts_opt_apply(in, in_opts))
        return EXIT_FAILURE;
    hts_opt_free(in_opts);

    if (hts_opt_apply(out, out_opts))
        return EXIT_FAILURE;
    hts_opt_free(out_opts);

    // Create and share the thread pool
    htsThreadPool p = {NULL, 0};
    if (opts.nthreads > 0) {
        p.pool = hts_tpool_init(opts.nthreads);
        if (!p.pool) {
            fprintf(stderr, "Error creating thread pool\n");
            exit_code = 1;
        } else {
            hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
            hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
        }
    }

    int ret;
    switch (hts_get_format(in)->category) {
    case sequence_data:
        ret = sam_loop(argc, argv, optind, &opts, in, out);
        break;

    case variant_data:
        ret = vcf_loop(argc, argv, optind, &opts, in, out);
        break;

    default:
        fprintf(stderr, "Unsupported or unknown category of data in input file\n");
        return EXIT_FAILURE;
    }

    if (ret != 0)
        exit_code = EXIT_FAILURE;

    ret = hts_close(out);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        exit_code = EXIT_FAILURE;
    }
    ret = hts_close(in);
    if (ret < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = EXIT_FAILURE;
    }

    if (p.pool)
        hts_tpool_destroy(p.pool);

    return exit_code;
}
