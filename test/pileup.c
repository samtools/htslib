/*  test/pileup.c -- simple pileup tester

    Copyright (C) 2014,2018-2019, 2024 Genome Research Ltd.

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

/*
  The output from this program isn't quite the same as that from
  `samtools mpileup`.  It doesn't print the reference base column,
  it puts brackets around insertion sequences to make them easier to spot
  and it writes empty brackets after a reported deletion.

  The output from `samtools mpileup` can be converted to the same format like
  this:

samtools mpileup -B -Q 0 in.bam | perl -lane \
  'pop(@F);
   splice(@F, 2, 1);
   $F[3] =~ s/\+(\d+)([ACGTN]+)/sprintf("+%d(%s)%s",$1,substr($2,0,$1),substr($2,$1))/ieg;
   $F[3] =~ s/\-(\d+)([ACGTN]+)/sprintf("-%d()%s",$1,substr($2,$1))/ieg;
   print join("\t", @F);'

 */

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>

#include "../htslib/sam.h"
#include "../htslib/kstring.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef struct ptest_t {
    const char *fname;
    samFile *fp;
    sam_hdr_t *fp_hdr;
} ptest_t;

static int readaln(void *data, bam1_t *b) {
    ptest_t *g = (ptest_t*)data;
    int ret;

    while (1) {
        ret = sam_read1(g->fp, g->fp_hdr, b);
        if (ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        break;
    }

    return ret;
}

static int print_pileup_seq(const bam_pileup1_t *p, int n) {
    kstring_t ks = { 0, 0, NULL };
    int i;

    for (i = 0; i < n; i++, p++) {
        uint8_t *seq = bam_get_seq(p->b);
        int del_len, is_rev = bam_is_rev(p->b);

        if (p->is_head)
            putchar('^'), putchar('!'+MIN(p->b->core.qual,93));

        if (p->is_del)
            putchar(p->is_refskip ? (is_rev ? '<' : '>') : '*');
        else {
            unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
            putchar(is_rev ? tolower(c) : toupper(c));
        }

        del_len = -p->indel;
        if (p->indel > 0) {
            int j, len = bam_plp_insertion(p, &ks, &del_len);
            if (len < 0) {
                perror("bam_plp_insertion");
                goto fail;
            }
            printf("%+d(", len);
            for (j = 0; j < len; j++)
                putchar(is_rev ?
                        tolower((uint8_t) ks.s[j]) :
                        toupper((uint8_t) ks.s[j]));
            putchar(')');
        }
        if (del_len > 0) {
            printf("-%d()", del_len);
        }
        if (p->is_tail)
            putchar('$');
    }
    free(ks.s);
    return 0;

 fail:
    free(ks.s);
    return -1;
}

static void print_pileup_qual(const bam_pileup1_t *p, int n) {
    int i;

    for (i = 0; i < n; i++, p++) {
        uint8_t *qual = bam_get_qual(p->b);
        uint8_t q = '~';
        if (p->qpos < p->b->core.l_qseq &&
            qual[p->qpos]+33 < '~')
            q = qual[p->qpos]+33;
        putchar(q);
    }
}

static int test_pileup(ptest_t *input) {
    bam_plp_t plp = NULL;
    const bam_pileup1_t *p;
    int tid, pos, n = 0;

    plp = bam_plp_init(readaln, input);
    if (!plp) {
        perror("bam_plp_init");
        goto fail;
    }
    while ((p = bam_plp_auto(plp, &tid, &pos, &n)) != 0) {
        if (tid < 0) break;
        if (tid >= input->fp_hdr->n_targets) {
            fprintf(stderr,
                    "bam_plp_auto returned tid %d >= header n_targets %d\n",
                    tid, input->fp_hdr->n_targets);
            goto fail;
        }

        printf("%s\t%d\t%d\t", input->fp_hdr->target_name[tid], pos+1, n);

        if (print_pileup_seq(p, n) < 0)
            goto fail;

        putchar('\t');
        print_pileup_qual(p, n);

        putchar('\n');
    }
    if (n < 0) {
        fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", input->fname);
        goto fail;
    }

    bam_plp_destroy(plp);
    return 0;

 fail:
    bam_plp_destroy(plp);
    return -1;
}

static int test_mpileup(ptest_t *input) {
    bam_mplp_t iter = NULL;
    const bam_pileup1_t *pileups[1] = { NULL };
    int n_plp[1] = { 0 };
    int tid, pos, n = 0;

    iter = bam_mplp_init(1, readaln, (void **) &input);
    if (!iter) {
        perror("bam_plp_init");
        goto fail;
    }
    if (bam_mplp_init_overlaps(iter) < 0) {
        perror("bam_mplp_init_overlaps");
        goto fail;
    }

    while ((n = bam_mplp_auto(iter, &tid, &pos, n_plp, pileups)) > 0) {
        if (tid < 0) break;
        if (tid >= input->fp_hdr->n_targets) {
            fprintf(stderr,
                    "bam_mplp_auto returned tid %d >= header n_targets %d\n",
                    tid, input->fp_hdr->n_targets);
            goto fail;
        }

        printf("%s\t%d\t%d\t", input->fp_hdr->target_name[tid], pos+1, n_plp[0]);

        if (print_pileup_seq(pileups[0], n_plp[0]) < 0)
            goto fail;

        putchar('\t');
        print_pileup_qual(pileups[0], n_plp[0]);

        putchar('\n');
    }
    if (n < 0) {
        fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", input->fname);
        goto fail;
    }

    bam_mplp_destroy(iter);
    return 0;

 fail:
    bam_mplp_destroy(iter);
    return -1;
}

int main(int argc, char **argv) {
    ptest_t g = { NULL, NULL, NULL };
    int use_mpileup = 0, opt;

    while ((opt = getopt(argc, argv, "m")) != -1) {
        switch (opt) {
        case 'm':
            use_mpileup = 1;
            break;
        default:
            fprintf(stderr, "Usage: %s [-m] <sorted.sam>\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Usage: %s [-m] <sorted.sam>\n", argv[0]);
        return EXIT_FAILURE;
    }

    g.fname = argv[optind];
    g.fp = sam_open(g.fname, "r");
    if (!g.fp) {
        fprintf(stderr, "Couldn't open \"%s\" : %s", g.fname, strerror(errno));
        goto fail;
    }
    g.fp_hdr = sam_hdr_read(g.fp);
    if (!g.fp_hdr) {
        fprintf(stderr, "Couldn't read header from \"%s\" : %s",
                g.fname, strerror(errno));
        goto fail;
    }

    if (use_mpileup) {
        if (test_mpileup(&g) < 0)
            goto fail;
    } else {
        if (test_pileup(&g) < 0)
            goto fail;
    }

    sam_hdr_destroy(g.fp_hdr);
    sam_close(g.fp);

    return EXIT_SUCCESS;

 fail:
    if (g.fp_hdr) sam_hdr_destroy(g.fp_hdr);
    if (g.fp) sam_close(g.fp);
    return EXIT_FAILURE;
}
