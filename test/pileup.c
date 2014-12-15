/*  test/pileup.c -- simple pileup tester

    Copyright (C) 2014 Genome Research Ltd.

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

#include <stdio.h>
#include "htslib/sam.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef struct ptest_t {
    samFile *fp;
    bam_hdr_t *fp_hdr;
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

int main(int argc, char **argv) {
    bam_plp_t plp;
    const bam_pileup1_t *p;
    int tid, pos, n;
    ptest_t g;

    if (argc != 2) {
	fprintf(stderr, "Usage: pileup foo.bam\n");
	return 1;
    }

    g.fp = sam_open(argv[1], "r");
    g.fp_hdr = sam_hdr_read(g.fp);

    plp = bam_plp_init(readaln, &g);
    while ((p = bam_plp_auto(plp, &tid, &pos, &n)) != 0) {
	int i;

        if (tid < 0) break;

	printf("%2d\t%6d\t", tid, pos+1);

	for (i = 0; i < n; i++, p++) {
            uint8_t *seq = bam_get_seq(p->b);
	    if (p->is_head)
		putchar('^'), putchar('!'+MIN(p->b->core.qual,93));

	    if (p->is_del)
		putchar('*');
	    else
		putchar(seq_nt16_str[bam_seqi(seq, p->qpos)]);

	    if (p->indel > 0) {
		int j;
		printf("%+d(", p->indel);
		for (j = 1; j<=p->indel; j++)
		    putchar(seq_nt16_str[bam_seqi(seq, p->qpos+j-p->is_del)]);
		putchar(')');
	    }
	    if (p->indel < 0) {
		printf("%+d(?)", p->indel);
	    }
	    if (p->is_tail)
		putchar('$');
	}

	putchar('\n');
    }
    bam_plp_destroy(plp);

    bam_hdr_destroy(g.fp_hdr);
    sam_close(g.fp);

    return 0;
}
