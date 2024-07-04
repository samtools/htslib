/* test/test_fadix.c -- Test faidx interfaces

    Copyright (C) 2022 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

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
#include <stdlib.h>
#include <unistd.h>

#include "../htslib/faidx.h"

int file_compare(const char *file1, const char *file2) {
    FILE *f1 = NULL;
    FILE *f2 = NULL;
    unsigned int lno = 1;
    size_t got1, got2, i;
    char buf1[1024], buf2[1024];
    int ret = -1;

    f1 = fopen(file1, "rb");
    if (!f1) {
        perror(file1);
        goto out;
    }
    f2 = fopen(file2, "rb");
    if (!f2) {
        perror(file2);
        goto out;
    }

    do {
        got1 = fread(buf1, 1, sizeof(buf1), f1);
        got2 = fread(buf2, 1, sizeof(buf2), f2);

        for (i = 0; i < got1 && i < got2 && buf1[i] == buf2[i]; i++)
            lno += (buf1[i] == '\n');
        if (i < got1 || i < got2) {
            fprintf(stderr, "%s and %s differ at line %u\n",
                    file1, file2, lno);
            goto out;
        }
    } while (got1 > 0 && got2 > 0);

    if (ferror(f1)) {
        perror(file1);
        goto out;
    }
    if (ferror(f2)) {
        perror(file2);
        goto out;
    }

    if (got1 > 0 || got2 > 0) {
        fprintf(stderr, "EOF on %s at line %u\n",
                got1 ? file2 : file1, lno);
        goto out;
    }

    ret = 0;
 out:
    if (f1) fclose(f1);
    if (f2) fclose(f2);
    return ret;
}

faidx_t * load_index(const char *fn, const char *fnfai, const char *fngzi,
                     int flags, enum fai_format_options format) {
    faidx_t *fai = fai_load3_format(fn, fnfai, fngzi, flags, format);
    if (!fai) {
        fprintf(stderr, "Failed: fai_load3(%s, %s, %s, %d, %d)\n",
                fn, fnfai ? fnfai : "NULL", fngzi ? fngzi : "NULL", flags,
                (int) format);
        return NULL;
    }
    return fai;
}

int do_retrieval(const char *fn, const char *fnfai, const char *fngzi,
                 int flags, enum fai_format_options format, const char *fnout,
                 const char *interface, int nreg, char **regions) {
    int i, use_64bit = 1, use_parse_reg = 0, use_adjust_reg = 0;
    faidx_t *fai = NULL;
    FILE *out = stdout;

    if (interface) {
        if (strcmp(interface, "fai_fetch") == 0) {
            use_64bit = 0;
        } else if (strcmp(interface, "faidx_fetch_seq") == 0) {
            use_64bit = 0;
            use_parse_reg = 1;
        } else if (strcmp(interface, "faidx_fetch_seq64") == 0
                   || strcmp(interface, "fai_parse_region") == 0) {
            use_parse_reg = 1;
        } else if (strcmp(interface, "fai_adjust_region") == 0) {
            use_parse_reg = 1;
            use_adjust_reg = 1;
        }
    }

    if (fnout) {
        out = fopen(fnout, "wb");
        if (!out) {
            perror(fnout);
            return -1;
        }
    }

    fai = load_index(fn, fnfai, fngzi, flags, format);
    if (!fai)
        goto fail;

    for (i = 0; i < nreg; i++) {
        hts_pos_t len = 0, pos, beg = 0, end = 0;
        int tid = 0;
        char *seq = NULL;
        size_t l;

        if (use_parse_reg) {
            const char *e = fai_parse_region(fai, regions[i],
                                             &tid, &beg, &end, 0);
            if (e == NULL) {
                fprintf(stderr, "Failed: "
                        "fai_parse_region(fai, %s, &tid, &beg, &end, 0)\n",
                        regions[i]);
                goto fail;
            }
            if (use_adjust_reg) {
                hts_pos_t orig_beg = beg, orig_end = end;
                int r = fai_adjust_region(fai, tid, &beg, &end);
                if (r < 0
                    || (((r & 1) != 0) ^ (beg != orig_beg))
                    || (((r & 2) != 0) ^ (end != orig_end))) {
                    fprintf(stderr, "Failed: fai_adjust_region(fai, %d, "
                            "%"PRIhts_pos", %"PRIhts_pos") returned %d\n"
                            "After: beg = %"PRIhts_pos" end = %"PRIhts_pos"\n",
                            tid, orig_beg, orig_end, r, beg, end);
                    goto fail;
                }
            }
            if (use_64bit) {
                seq = faidx_fetch_seq64(fai, faidx_iseq(fai, tid),
                                        beg, end - 1, &len);
            } else {
                int ilen = 0;
                seq = faidx_fetch_seq(fai, faidx_iseq(fai, tid),
                                      beg, end - 1, &ilen);
                len = ilen;
            }
            if (!seq) {
                fprintf(stderr, "Failed: faidx_fetch_seq%s(fai, %s, "
                        "%"PRIhts_pos", %"PRIhts_pos", &len)\n",
                        use_64bit ? "64" : "", faidx_iseq(fai, tid), beg, end);
                goto fail;
            }
        } else {
            if (use_64bit) {
                seq = fai_fetch64(fai, regions[i], &len);
            } else {
                int ilen = 0;
                seq = fai_fetch(fai, regions[i], &ilen);
                len = ilen;
            }
            if (!seq) {
                fprintf(stderr, "Failed: fai_fetch%s(fai, %s, &len)\n",
                        use_64bit ? "64" : "", regions[i]);
                goto fail;
            }
        }

        l = strlen(seq);
        fprintf(out, "%c%s length: %"PRIhts_pos"\n",
                format == FAI_FASTQ ? '@' : '>', regions[i], len);
        for (pos = 0; pos < l; pos += 50) {
            fprintf(out, "%.*s\n", 50, seq + pos);
        }
        free(seq);
        if (format == FAI_FASTQ) {
            hts_pos_t qual_len = 0;
            char *qual;
            if (use_parse_reg) {
                if (use_64bit) {
                    qual = faidx_fetch_qual64(fai, faidx_iseq(fai, tid),
                                              beg, end - 1, &qual_len);
                } else {
                    int ilen = 0;
                    qual = faidx_fetch_qual(fai, faidx_iseq(fai, tid),
                                            beg, end - 1, &ilen);
                    qual_len = ilen;
                }
            } else {
                if (use_64bit) {
                    qual = fai_fetchqual64(fai, regions[i], &qual_len);
                } else {
                    int ilen = 0;
                    qual = fai_fetchqual(fai, regions[i], &ilen);
                    qual_len = ilen;
                }
                if (!qual) {
                    fprintf(stderr, "Failed: fai_fetchqual64(fai, %s, &len)\n",
                            regions[i]);
                    goto fail;
                }
            }
            if (qual_len != len) {
                fprintf(stderr,
                        "Sequence and quality lengths differ for %s %s\n",
                        fn, regions[i]);
                free(qual);
                goto fail;
            }
            fprintf(out, "+\n");
            l = strlen(qual);
            for (pos = 0; pos < l; pos+=50) {
                fprintf(out, "%.*s\n", 50, qual + pos);
            }
            free(qual);
        }
    }

    fai_destroy(fai);

    if (fnout) {
        if (fclose(out) != 0) {
            perror(fnout);
            return -1;
        }
    }
    return 0;

 fail:
    if (fai)
        fai_destroy(fai);
    if (fnout)
        fclose(out);

    return -1;
}

int test_fai_line_length(const char *fn, const char *fnfai, const char *fngzi,
                         enum fai_format_options format, const char *expected,
                         const char *reg) {
    hts_pos_t found_len;
    faidx_t *fai = NULL;

    fai = load_index(fn, fnfai, fngzi, 0, format);
    if (!fai)
        return -1;

    found_len = fai_line_length(fai, reg);
    fai_destroy(fai);
    if (expected) {
        long long exp_len = strtoll(expected, NULL, 10);
        if (found_len != exp_len) {
            fprintf(stderr, "Unexpected result %"PRIhts_pos" from "
                    "fai_line_length, expected %s\n", found_len, expected);
            return -1;
        }
    } else {
        printf("%"PRIhts_pos"\n", found_len);
    }
    return 0;
}

int test_faidx_has_seq(const char *fn, const char *fnfai, const char *fngzi,
                       enum fai_format_options format, const char *expected,
                       const char *seq) {
    int res;
    faidx_t *fai = NULL;

    fai = load_index(fn, fnfai, fngzi, 0, format);
    if (!fai)
        return -1;

    res = faidx_has_seq(fai, seq);
    fai_destroy(fai);
    if (expected) {
        long exp_res = strtol(expected, NULL, 10);
        if (res != exp_res) {
            fprintf(stderr, "Unexpected result %d from faidx_has_seq(%s) "
                    "expected %s\n", res, seq, expected);
            return -1;
        }
    } else {
        printf("%d\n", res);
    }
    return 0;
}

int test_faidx_iseq(const char *fn, const char *fnfai, const char *fngzi,
                    enum fai_format_options format, const char *expected,
                    const char *index) {
    const char *found_name = NULL;
    int idx = atoi(index);
    faidx_t *fai = NULL;

    fai = load_index(fn, fnfai, fngzi, 0, format);
    if (!fai)
        return -1;

    found_name = faidx_iseq(fai, idx);

    if (expected) {
        if (!found_name || strcmp(found_name, expected) != 0) {
            fprintf(stderr, "Unexpected result %s from faidx_iseq(fai, %d), "
                    "expected %s\n", found_name ? found_name : "(null)",
                    idx, expected);
            fai_destroy(fai);
            return -1;
        }
    } else {
        printf("%s\n", found_name ? found_name : "(null)");
    }

    fai_destroy(fai);
    return 0;
}

int test_faidx_seq_len(const char *fn, const char *fnfai, const char *fngzi,
                       enum fai_format_options format, const char *expected,
                       const char *seq) {
    int found_len;
    faidx_t *fai = NULL;

    fai = load_index(fn, fnfai, fngzi, 0, format);
    if (!fai)
        return -1;

    found_len = faidx_seq_len(fai, seq);
    fai_destroy(fai);

    if (expected) {
        int exp_len = atoi(expected);
        if (found_len != exp_len) {
            fprintf(stderr, "Unexpected result %d from faidx_seq_len(fai, %s) "
                    "expected %s\n", found_len, seq, expected);
            return -1;
        }
    } else {
        printf("%d\n", found_len);
    }

    return 0;
}

int test_faidx_seq_len64(const char *fn, const char *fnfai, const char *fngzi,
                         enum fai_format_options format, const char *expected,
                         const char *seq) {
    hts_pos_t found_len;
    faidx_t *fai = NULL;

    fai = load_index(fn, fnfai, fngzi, 0, format);
    if (!fai)
        return -1;

    found_len = faidx_seq_len(fai, seq);
    fai_destroy(fai);

    if (expected) {
        long long exp_len = strtoll(expected, NULL, 10);
        if (found_len != exp_len) {
            fprintf(stderr, "Unexpected result %"PRIhts_pos
                    " from fai_seq_len64(fai, %s) expected %s\n",
                    found_len, seq, expected);
            return -1;
        }
    } else {
        printf("%"PRIhts_pos"\n", found_len);
    }

    return 0;
}

void usage(FILE *out, const char *arg0) {
    fprintf(out,
            "Usage: %s [-c] -i fasta/q [-f fai_file] [-g gzi_file] [-e expected_fai]\n"
            "       %s [-cQ] -i fasta/q [-f fai_file] [-g gzi_file] [region]\n"
            "       %s -t FUNC -i fasta/q [-f fai_file] [-g gzi_file] [-e expected] <PARAM>\n"
            "       %s -h\n",
            arg0, arg0, arg0, arg0);
}

void help(FILE *out, const char *arg0) {
    usage(out, arg0);
    fprintf(out,
            "Options:\n"
            "  -i FILE      Input file\n"
            "  -f FILE      Fasta/q index file name\n"
            "  -g FILE      Bgzip index file name\n"
            "  -o FILE      Output file name\n"
            "  -e FILE|STR  Expected output\n"
            "  -c           Set FAI_CREATE flag\n"
            "  -Q           Output fastq format\n"
            "  -t FUNC      Test function\n"
            "  -h           Print this help\n"
            "\n"
            "Expected output is compared to the FAI file in indexing mode;"
            " the output file\n"
            "in retrieval mode; "
            "expected output for various -t function tests.\n"
            "\n"
            "Unit tests (-t option):\n"
            "   fai_line_length, faidx_has_seq, faidx_iseq, faidx_seq_len, faidx_seq_len64\n"
            "In retrieval mode, -t can change the functions used to fetch data:\n"
            "   fai_fetch, fai_fetch64, faidx_fetch_seq, faidx_fetch_seq64,\n"
            "   fai_parse_region, fai_adjust_region\n"
            "\n");
}

int main(int argc, char **argv) {
    int opt;
    const char *fn = NULL;
    const char *fnout = NULL;
    const char *fnfai = NULL;
    const char *fngzi = NULL;
    const char *expected = NULL;
    const char *func = "";
    int flags = 0;
    enum fai_format_options format = FAI_FASTA;
    int res;

    while ((opt = getopt(argc, argv, "i:f:g:o:e:t:cQh")) > 0) {
        switch (opt) {
        case 'i':
            fn = optarg;
            break;
        case 'f':
            fnfai = optarg;
            break;
        case 'g':
            fngzi = optarg;
            break;
        case 'o':
            fnout = optarg;
            break;
        case 'e':
            expected = optarg;
            break;
        case 'c':
            flags |= FAI_CREATE;
            break;
        case 'Q':
            format = FAI_FASTQ;
            break;
        case 't':
            func = optarg;
            break;
        case 'h':
            help(stdout, argv[0]);
            return EXIT_SUCCESS;
        default:
            usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (!fn) {
        usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    if (optind == argc) {
        // Index building mode
        res = fai_build3(fn, fnfai, fngzi);
        if (res) {
            fprintf(stderr, "Failed: fai_build3(%s, %s, %s)\n",
                    fn, fnfai ? fnfai : "NULL", fngzi ? fngzi : "NULL");
        } else if (expected) {
            res = file_compare(fnfai, expected);
        }
    } else {
        if (strcmp(func, "fai_line_length") == 0) {
            res = test_fai_line_length(fn, fnfai, fngzi, format, expected,
                                       argv[optind]);
        } else if (strcmp(func, "faidx_has_seq") == 0) {
            res = test_faidx_has_seq(fn, fnfai, fngzi, format, expected,
                                     argv[optind]);
        } else if (strcmp(func, "faidx_iseq") == 0) {
            res = test_faidx_iseq(fn, fnfai, fngzi, format, expected,
                                  argv[optind]);
        } else if (strcmp(func, "faidx_seq_len") == 0) {
            res = test_faidx_seq_len(fn, fnfai, fngzi, format, expected,
                                     argv[optind]);
        } else if (strcmp(func, "faidx_seq_len64") == 0) {
            res = test_faidx_seq_len64(fn, fnfai, fngzi, format, expected,
                                       argv[optind]);
        } else {
            res =  do_retrieval(fn, fnfai, fngzi, flags, format, fnout,
                                func, argc - optind, &argv[optind]);
            if (res == 0 && fnout && expected) {
                res = file_compare(fnout, expected);
            }
        }
    }
    return res == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
