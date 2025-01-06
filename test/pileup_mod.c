/*  test/pileup_mod.c -- simple pileup tester with base modifications

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
#include <ctype.h>
#include <math.h>
#include "../htslib/sam.h"

typedef struct {
    samFile *fp;
    sam_hdr_t *h;
}  plp_dat;

static int readaln(void *data, bam1_t *b) {
    plp_dat *dat = (plp_dat *)data;
    return sam_read1(dat->fp, dat->h, b);
}

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

// No modification reporting.
// This is just a simple base-line for comparison against mod_pileup1 for
// performance testing.
void process_pileup(sam_hdr_t *h, const bam_pileup1_t *p,
                    int tid, int pos, int n) {
    kstring_t s = {0,0};
    printf("%s\t%d\t", sam_hdr_tid2name(h, tid), pos);
    int i;
    for (i = 0; i < n; i++, p++) {
        if (p->is_del) {
            putchar('*');
            continue;
        }

        uint8_t *seq = bam_get_seq(p->b);
        uint8_t *qual = bam_get_qual(p->b);
        unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
        putchar(c);
        kputc(MIN('~','!'+qual[p->qpos]), &s);
    }
    putchar('\t');
    puts(s.l ? s.s : "");

    free(s.s);
}

// Initialise and destroy the base modifier state data. This is called
// as each new read is added or removed from the pileups.
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    if (bam_parse_basemod(b, m) < 0) {
        hts_base_mod_state_free(m);
        return -1;
    }
    cd->p = m;
    return 0;
}

int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state_free(cd->p);
    return 0;
}

// Report a line of pileup, including base modifications inline with
// the sequence (including qualities), as [<strand><dir><qual>...]
void process_mod_pileup1(sam_hdr_t *h, const bam_pileup1_t *p,
                         int tid, int pos, int n) {
    kstring_t s = {0,0};
    printf("%s\t%d\t", sam_hdr_tid2name(h, tid), pos);
    int i;
    for (i = 0; i < n; i++, p++) {
        if (p->is_del) {
            putchar('*');
            continue;
        }

        uint8_t *seq = bam_get_seq(p->b);
        uint8_t *qual = bam_get_qual(p->b);
        unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
        putchar(c);
        kputc(MIN('~','!'+qual[p->qpos]), &s);

        // Simple mod detection; assumes at most 5 mods
        hts_base_mod_state *m = p->cd.p;
        hts_base_mod mod[5];
        int nm;
        if ((nm = bam_mods_at_qpos(p->b, p->qpos, m, mod, 5)) > 0) {
            int j;
            putchar('[');
            for (j = 0; j < nm && j < 5; j++) {
                if (mod[j].modified_base < 0)
                    // ChEBI
                    printf("%c(%d)%d", "+-"[mod[j].strand],
                           -mod[j].modified_base, mod[j].qual);
                else
                    printf("%c%c%d", "+-"[mod[j].strand],
                           mod[j].modified_base, mod[j].qual);
            }
            putchar(']');
        }
    }
    putchar('\t');
    puts(s.l ? s.s : "");

    free(s.s);
}

// Report a line of pileup, including base modifications.
// This replaces the base with the mod call (NB this can be confusing
// as both C and G can map to m depending on orientation).
// It also reports qualities in the QUAl column, remapped to
// phred scale as only one single mod is supported and hence extreme
// unlikely probabilities shouldn't be reported (although we don't
// scan to pick the highest).
void process_mod_pileup2(sam_hdr_t *h, const bam_pileup1_t *p,
                        int tid, int pos, int n) {
    kstring_t s = {0,0};
    printf("%s\t%d\t%d\t", sam_hdr_tid2name(h, tid), pos, n);
    int i;
    for (i = 0; i < n; i++, p++) {
        if (p->is_del) {
            putchar('*');
            continue;
        }

        uint8_t *seq = bam_get_seq(p->b);
        uint8_t *qual = bam_get_qual(p->b);
        unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];

        // Simple mod detection; assumes at most 2 non-ChEBI mods
        hts_base_mod_state *m = p->cd.p;
        int n, is_rev = bam_is_rev(p->b);
        hts_base_mod mod;
        char base;
        uint8_t q = qual[p->qpos];
        if ((n = bam_mods_at_qpos(p->b, p->qpos, m, &mod, 1)) > 0) {
            base = mod.modified_base;
            // base mod as phred scale
            q = -10 * log10(1-((mod.qual+0.5)/256)) + 0.5;
        } else {
            base = c;
        }

        // Case is inappropriate here as some mods (eg "a") are lc.
        // So we dim/bold them instead using ANSI escape codes.
        // It's a test script, so I'm not going to care about curses.
        if (is_rev) {
            printf("\033[2m%c\033[0m", base);
        } else {
            printf("\033[1m%c\033[0m", base);
        }
        kputc(MIN('~','!'+q), &s);
    }
    putchar('\t');
    puts(s.l ? s.s : "");

    free(s.s);
}

int main(int argc, char **argv) {
    int compact = 0;
    while (argc > 1 && strcmp(argv[1], "-c") == 0) {
        compact++;
        argc--;
        argv++;
    }

    samFile *in = sam_open(argc > 1 ? argv[1] : "-", "r");
    bam1_t *b = bam_init1();
    sam_hdr_t *h = sam_hdr_read(in);

    // Pileup iterator with constructor/destructor to parse base mod tags
    plp_dat dat = {
        .fp = in,
        .h = h,
    };
    bam_plp_t iter = bam_plp_init(readaln, &dat);
    bam_plp_constructor(iter, pileup_cd_create);
    bam_plp_destructor(iter, pileup_cd_destroy);

    const bam_pileup1_t *p;
    int tid, pos, n = 0;
    while ((p = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
        switch (compact) {
        case 0:
            process_mod_pileup1(h, p, tid, pos, n);
            break;
        case 1:
            process_mod_pileup2(h, p, tid, pos, n);
            break;
        default:
            process_pileup(h, p, tid, pos, n);
            break;
        }
    }
    bam_plp_destroy(iter);

    sam_close(in);
    bam_destroy1(b);
    sam_hdr_destroy(h);

    return n != 0;
}
