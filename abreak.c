#include <stdlib.h>
#include <unistd.h>
#include "sam.h"
#include "kstring.h"
#include "ksort.h"

typedef struct {
	int is_bam, min_len, min_q, max_gap;
	float mask_level;
} cmdopt_t;

typedef struct {
	int64_t L, n_un, l_un, n_dropped;
	int n_b[5], n_bg[5];
	int n, m;
	int *len;
} stat_t;

typedef struct {
	int tid, pos, len, qlen, rlen, flag, mapq, qbeg, clip[2];
} aln_t;

#define aln_lt(a, b) ((a).tid < (b).tid || ((a).tid == (b).tid && (a).pos < (b).pos))
KSORT_INIT(s2b, aln_t, aln_lt)

#define intr_lt(a, b) ((a) > (b))
KSORT_INIT(intr, int, intr_lt)

static void count_break(int c[5], int n_aa, aln_t *aa, const cmdopt_t *o)
{
	int i, b[5] = { n_aa, 0, 0, 0, 0 };
	for (i = 0; i < n_aa; ++i) {
		aln_t *p = &aa[i];
		if (p->mapq < o->min_q) continue;
		++b[1];
		if (p->qlen >= 100) {
			++b[2];
			if (p->qlen >= 200) {
				++b[3];
				if (p->qlen >= 500) ++b[4];
			}
		}
	}
	for (i = 0; i < 5; ++i)
		if (b[i]) c[i] += b[i] - 1;
}

static void analyze_aln(int n_aa, aln_t *aa, stat_t *s, cmdopt_t *o)
{
	int n_tmp, i;
	aln_t *p, *tmp = 0;
	// special treatment of unmapped
	if (n_aa == 1 && (aa[0].flag&4)) {
		++s->n_un; s->l_un += aa[0].len;
		return;
	}
	tmp = (aln_t*)alloca(n_aa * sizeof(aln_t));
	// apply mask_level
	if (n_aa > 1) { // multi-part alignment; check if this is a BWA-SW problem
		aln_t *p, *q;
		for (n_tmp = 0, p = aa; p < aa + n_aa; ++p) {
			int dropped = 0;
			for (q = tmp; q < tmp + n_tmp; ++q) {
				int beg = p->qbeg > q->qbeg? p->qbeg : q->qbeg;
				int end = p->qbeg + p->qlen < q->qbeg + q->qlen? p->qbeg + p->qlen : q->qbeg + q->qlen;
				if (beg < end && (double)(end - beg) > p->qlen * o->mask_level) {
					dropped = 1;
					break;
				}
			}
			if (!dropped) tmp[n_tmp++] = *p;
			else ++s->n_dropped;
		}
		memcpy(aa, tmp, n_tmp * sizeof(aln_t));
		n_aa = n_tmp;
		count_break(s->n_b, n_aa, aa, o);
	}
	if (s->n + n_aa > s->m) {
		s->m = s->n + n_aa;
		kroundup32(s->m);
		s->len = (int*)realloc(s->len, sizeof(int) * s->m);
	}
	for (p = aa; p < aa + n_aa; ++p) s->len[s->n++] = p->qlen;
	if (n_aa) { // still multi-part
		ks_introsort(s2b, n_aa, aa);
		for (i = 1; i < n_aa; ++i) {
			aln_t *p = &aa[i], *q = &aa[i-1];
			if (p->tid == q->tid && (p->flag&16) == (q->flag&16)) {
				int gapr = p->pos - (q->pos + q->rlen);
				int gapq = p->clip[0] - (q->clip[0] + q->qlen);
				if (gapr < 0) gapr = -gapr;
				if (gapq < 0) gapq = -gapq;
				if (gapr < o->max_gap && gapq < o->max_gap) {
					p->qlen = p->clip[0] + p->qlen - q->clip[0]; p->clip[0] = q->clip[0];
					p->rlen = p->pos + p->rlen - q->pos; p->pos = q->pos;
					q->flag |= 4;
				}
			}
		}
		for (n_tmp = i = 0; i < n_aa; ++i)
			if ((aa[i].flag&4) == 0) tmp[n_tmp++] = aa[i];
		memcpy(aa, tmp, n_tmp * sizeof(aln_t));
		n_aa = n_tmp;
		count_break(s->n_bg, n_aa, aa, o);
	}
}

int main_abreak(int argc, char *argv[])
{
	cmdopt_t o;
	samFile *in;
	bam_hdr_t *h;
	int c, k, n_aa = 0, m_aa = 0;
	aln_t a, *aa = 0;
	kstring_t last, out;
	stat_t s;
	bam1_t *b;
	
	memset(&a, 0, sizeof(aln_t));
	memset(&s, 0, sizeof(stat_t));
	memset(&last, 0, sizeof(kstring_t));
	memset(&out, 0, sizeof(kstring_t));
	memset(&o, 0, sizeof(cmdopt_t));
	o.min_len = 150; o.min_q = 10; o.mask_level = 0.5; o.max_gap = 500;
	while ((c = getopt(argc, argv, "l:b")) >= 0)
		if (c == 'b') o.is_bam = 1;
		else if (c == 'l') o.min_len = atoi(optarg);
		else if (c == 'q') o.min_q = atoi(optarg);
		else if (c == 'm') o.mask_level = atof(optarg);
		else if (c == 'g') o.max_gap = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   htscmd abreak [options] <aln.sam>|<aln.bam>\n\n");
		fprintf(stderr, "Options: -b        assume the input is BAM (default is SAM)\n");
		fprintf(stderr, "         -l INT    exclude contigs shorter than INT [%d]\n", o.min_len);
		fprintf(stderr, "         -q INT    exclude alignments with mapQ below INT [%d]\n", o.min_q);
		fprintf(stderr, "         -m FLOAT  exclude alignments overlapping another long alignment by FLOAT fraction [%g]\n", o.mask_level);
		fprintf(stderr, "         -g INT    join alignments separated by a gap shorter than INT bp [%d]\n\n", o.max_gap);
		fprintf(stderr, "Note: recommended BWA-SW setting is '-b9 -q16 -r1 -w500'\n\n");
		return 1;
	}

	in = sam_open(argv[optind], o.is_bam? "rb" : "r", 0);
	h = sam_hdr_read(in);
	b = bam_init1();
	while (sam_read1(in, h, b) >= 0) {
		uint32_t *cigar = bam_get_cigar(b);
		if (last.s == 0 || strcmp(last.s, bam_get_qname(b))) {
			analyze_aln(n_aa, aa, &s, &o);
			last.l = 0;
			kputs(bam_get_qname(b), &last);
			n_aa = 0;
		}
		a.qlen = a.rlen = 0; a.clip[0] = a.clip[1] = 0;
		for (k = 0; k < b->core.n_cigar; ++k) {
			int op = bam_cigar_op(cigar[k]);
			int oplen = bam_cigar_oplen(cigar[k]);
			if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP) a.qlen += oplen;
			if (bam_cigar_type(op)&2) a.rlen += oplen;
			if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
				a.clip[k == 0? 0 : 1] = oplen;
		}
		a.len = a.qlen + a.clip[0] + a.clip[1];
		if (a.len == 0) a.len = b->core.l_qseq;
		a.tid = b->core.tid; a.pos = b->core.pos; a.flag = b->core.flag; a.mapq = b->core.qual;
		a.qbeg = a.clip[!!(a.flag&BAM_FREVERSE)];
		if (a.len >= o.min_len) {
			if (n_aa == m_aa) {
				m_aa = m_aa? m_aa<<1 : 8;
				aa = (aln_t*)realloc(aa, sizeof(aln_t) * m_aa);
			}
			aa[n_aa++] = a;
		}
	}
	analyze_aln(n_aa, aa, &s, &o);
	bam_hdr_destroy(h);
	sam_close(in);
	{
		uint64_t L = 0, tmp;
		int N50;
		ks_introsort(intr, s.n, s.len);
		for (k = 0; k < s.n; ++k) L += s.len[k];
		for (k = 0, tmp = 0; k < s.n; ++k)
			if ((tmp += s.len[k]) >= L/2) break;
		N50 = s.len[k];
		printf("Number of unmapped contigs: %ld\n", (long)s.n_un);
		printf("Total length of unmapped contigs: %ld\n", (long)s.l_un);
		printf("Number of alignments dropped due to excessive overlaps: %ld\n", (long)s.n_dropped);
		printf("Mapped contig bases: %ld\n", (long)L);
		printf("Mapped N50: %d\n", N50);
		printf("Number of break points: %d\n", s.n_b[0]);
		printf("Number of Q%d break points longer than (0,100,200,500)bp: (%d,%d,%d,%d)\n", o.min_q, s.n_b[1], s.n_b[2], s.n_b[3], s.n_b[4]);
		printf("Number of break points after patching gaps short than %dbp: %d\n", o.max_gap, s.n_bg[0]);
		printf("Number of Q%d break points longer than (0,100,200,500)bp after gap patching: (%d,%d,%d,%d)\n", o.min_q, s.n_bg[1], s.n_bg[2], s.n_bg[3], s.n_bg[4]);
	}
	return 0;
}
