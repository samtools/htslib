#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "kstring.h"
#include "sam.h"

int main_bam2bed(int argc, char *argv[])
{
	samFile *in;
	bam_hdr_t *h;
	bam1_t *b;
	int c, is_sam = 0, is_ext = 0;
	kstring_t str;

	while ((c = getopt(argc, argv, "Se")) >= 0)
		if (c == 'S') is_sam = 1;
		else if (c == 'e') is_ext = 1;
	
	if (argc == optind) {
		fprintf(stderr, "Usage: bam2bed [-Se] <in.bam> | <in.sam>\n");
		return 1;
	}
	
	in = sam_open(argv[optind], is_sam? "r" : "rb", 0);
	h = sam_hdr_read(in);

	b = bam_init1();
	str.l = str.m = 0; str.s = 0;
	while (sam_read1(in, h, b) >= 0) {
		const bam1_core_t *c = &b->core;
		uint32_t *cigar = bam_get_cigar(b);
		int n[16], k;
		if (c->tid < 0 || (c->flag&BAM_FUNMAP) || c->n_cigar <= 0) continue;
		memset(n, 0, 16 * sizeof(int));
		for (k = 0; k < c->n_cigar; ++k)
			n[bam_cigar_op(cigar[k])] += bam_cigar_oplen(cigar[k]);
		str.l = 0;
		kputs(h->target_name[c->tid], &str);
		kputc('\t', &str); kputw(c->pos, &str);
		kputc('\t', &str); kputw(c->pos + n[0] + n[2] + n[3] + n[7] + n[8], &str);
		kputc('\t', &str); kputs(bam_get_qname(b), &str);
		kputc('\t', &str); kputw(c->qual, &str);
		kputc('\t', &str); kputc(bam_is_rev(b)? '-' : '+', &str);
		if (is_ext) {
			int qbeg = 0, qlen = n[0] + n[1] + n[7] + n[8];
			int op = bam_cigar_op(cigar[0]), nm = n[1] + n[2];
			uint8_t *s;
			if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
				qbeg = bam_cigar_oplen(cigar[0]);
			if (bam_is_rev(b)) qbeg = c->l_qseq - (qbeg + qlen);
			kputc('\t', &str); kputw(qbeg, &str);
			kputc('\t', &str); kputw(qbeg + qlen, &str);
			if ((s = bam_aux_get(b, "NM")) != 0) nm = bam_aux2i(s);
			ksprintf(&str, "\t%.6f", (double)nm / (n[0] + n[1] + n[2] + n[7] + n[8]));
		}
		puts(str.s);
	}
	free(str.s);
	bam_destroy1(b);

	bam_hdr_destroy(h);
	sam_close(in);
	return 0;
}
