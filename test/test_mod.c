/*  test/test_mod.c -- testing of base modification functions

    Copyright (C) 2020-2021 Genome Research Ltd.

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
This tests multiple APIs.  The simplest is to parse the MM/ML tags with
bam_parse_basemod and then call bam_mods_at_next_pos once for each base in
the bam sequence to check for modifications.

Ie:

    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod(b, m); // b=bam1_t pointer
    hts_base_mod mods[5];
    for (i = 0; i < b->core.l_qseq; i++) {
        n = bam_mods_at_next_pos(b, m, mods, 5);
        for (j = 0; j < n && j < 5; j++) {
            // Report 'n'th mod at seq pos 'i'.
            // mods[j].modified_base holds the base mod itself, with
            // mods[j].canonical_base, mods[j].strand and mods[j].qual
            // also present in hts_base_mod struct.
            // ...
        }
    }
    hts_base_mod_state_free(m);

The extended mode has the same loop above, but calls bam_mods_query_type
to return additional meta-data including the strand, canonical base and
whether the base modification is recorded implicitly or explicitly:

            int ret = bam_mods_query_type(m, mods[j].modified_base,
                                          &m_strand, &m_implicit,
                                          &m_canonical);

Looping over every base in the sequence is not particularly efficient
however unless this fits your natural processing order.  The alternative
is to call bam_next_base_mod to iterate only over modified locations:

    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod(b, m); // b=bam1_t pointer
    hts_base_mod mods[5];
    while ((n=bam_next_basemod(b, m, mods, 5, &pos)) > 0) {
        for (j = 0; j < n && j < 5; j++) {
            // Report 'n'th mod at sequence position 'pos'
        }
    }
    hts_base_mod_state_free(m);

*/

#include <config.h>
#include <stdio.h>

#include "../htslib/sam.h"

static char *code(int id) {
    static char code[20];
    if (id > 0) {
        code[0] = id;
        code[1] = 0;
    } else {
        sprintf(code, "(%d)", -id);
    }

    return code;
}

int main(int argc, char **argv) {
    char out[1024] = {0};
    int extended = 0;

    if (argc > 1 && strcmp(argv[1], "-x") == 0) {
        extended = 1;
        argv++;
        argc--;
    }

    if (argc < 2)
        return 1;

    samFile *in = sam_open(argv[1], "r");
    if (!in)
        return 1;

    bam1_t *b = bam_init1();
    sam_hdr_t *h = sam_hdr_read(in);
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    if (!h || !b || !m)
        goto err;

    int r;
    while ((r = sam_read1(in, h, b)) >= 0) {
        if (bam_parse_basemod(b, m) < 0) {
            fprintf(stderr, "Failed to parse MM/ML aux tags\n");
            goto err;
        }

        // per-base iterator
        int i, j, n;
        hts_base_mod mods[5];
        for (i = 0; i < b->core.l_qseq; i++) {
            char line[8192], *lp = line;
            n = bam_mods_at_next_pos(b, m, mods, 5);
            lp += sprintf(lp, "%d\t%c\t",
                          i, seq_nt16_str[bam_seqi(bam_get_seq(b), i)]);
            for (j = 0; j < n && j < 5; j++) {
                if (extended) {
                    int m_strand, m_implicit;
                    char m_canonical;
                    int ret = bam_mods_query_type(m, mods[j].modified_base,
                                                  &m_strand, &m_implicit,
                                                  &m_canonical);
                    if (ret < 0 ||
                        m_canonical != mods[j].canonical_base ||
                        m_strand    != mods[j].strand)
                        goto err;
                    lp += sprintf(lp, "%c%c%s%c%d ",
                                  mods[j].canonical_base,
                                  "+-"[mods[j].strand],
                                  code(mods[j].modified_base),
                                  "?."[m_implicit],
                                  mods[j].qual);
                } else {
                    lp += sprintf(lp, "%c%c%s%d ",
                                  mods[j].canonical_base,
                                  "+-"[mods[j].strand],
                                  code(mods[j].modified_base),
                                  mods[j].qual);
                }
            }
            *lp++ = '\n';
            *lp++ = 0;

            if (argc > 1)
                printf("%s", line);
            else
                strcat(out, line);
        }

        if (argc > 1) puts("---");

        bam_parse_basemod(b, m);

        // List possible mod choices.
        int *all_mods;
        int all_mods_n = 0;
        all_mods = bam_mods_recorded(m, &all_mods_n);
        printf("Present:");
        for (i = 0; i < all_mods_n; i++)
            printf(all_mods[i] > 0 ? " %c" : " #%d", all_mods[i]);
        putchar('\n');

        int pos;
        while ((n=bam_next_basemod(b, m, mods, 5, &pos)) > 0) {
            char line[8192]={0}, *lp = line;
            lp += sprintf(lp, "%d\t%c\t", pos,
                          seq_nt16_str[bam_seqi(bam_get_seq(b), pos)]);
            for (j = 0; j < n && j < 5; j++) {
                lp += sprintf(lp, "%c%c%s%d ",
                              mods[j].canonical_base,
                              "+-"[mods[j].strand],
                              code(mods[j].modified_base),
                              mods[j].qual);
            }
            *lp++ = '\n';
            *lp++ = 0;

            if (argc > 1)
                printf("%s", line);
            else
                strcat(out, line);
        }
        if (n < 0)
            goto err;

        if (argc > 1) puts("\n===\n");
    }
    fflush(stdout);
    if (sam_close(in) != 0 || r < -1)
        goto err;

    bam_destroy1(b);
    sam_hdr_destroy(h);
    hts_base_mod_state_free(m);
    return 0;

 err:
    bam_destroy1(b);
    sam_hdr_destroy(h);
    hts_base_mod_state_free(m);
    return 1;
}
