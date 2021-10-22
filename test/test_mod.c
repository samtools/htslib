/*  test/test_mod.c -- testing of base modification functions

    Copyright (C) 2020 Genome Research Ltd.

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
            for (j = 0; j < n && j < 5; j++)
                lp += sprintf(lp, "%c%c%s%d ",
                              mods[j].canonical_base,
                              "+-"[mods[j].strand],
                              code(mods[j].modified_base),
                              mods[j].qual);
            *lp++ = '\n';
            *lp++ = 0;

            if (argc > 1)
                printf("%s", line);
            else
                strcat(out, line);
        }

        if (argc > 1) puts("---");

        bam_parse_basemod(b, m);

        int pos;
        while ((n=bam_next_basemod(b, m, mods, 5, &pos)) > 0) {
            char line[8192]={0}, *lp = line;
            lp += sprintf(lp, "%d\t%c\t", pos,
                          seq_nt16_str[bam_seqi(bam_get_seq(b), pos)]);
            for (j = 0; j < n && j < 5; j++)
                lp += sprintf(lp, "%c%c%s%d ",
                              mods[j].canonical_base,
                              "+-"[mods[j].strand],
                              code(mods[j].modified_base),
                              mods[j].qual);
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
