// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include "sam.h"
#include "faidx.h"
#include "ksort.h"

const char *hts_parse_reg(const char *s, int *beg, int *end);

typedef struct {     // auxiliary data structure
	BGZF *fp;        // the file handler
	hts_itr_t *itr;  // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->itr? bam_itr_next(aux->fp, aux->itr, b) : bam_read1(aux->fp, b);
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
		else if (aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
	}
	return ret;
}

typedef struct {
	uint32_t is_skip:1, is_rev:1, b:14, q:16;
	int indel;
	uint64_t hash;
	uint64_t pos;
} allele_t;

#define allele_lt(a, b) ((a).indel < (b).indel || ((a).indel == (b).indel && (a).hash < (b).hash))
KSORT_INIT(allele, allele_t, allele_lt)

static inline allele_t pileup2allele(const bam_pileup1_t *p, int min_baseQ, uint64_t pos)
{
	allele_t a;
	int i;
	const uint8_t *seq = bam_get_seq(p->b);
	a.q = bam_get_qual(p->b)[p->qpos];
	a.is_rev = bam_is_rev(p->b);
	a.is_skip = (p->is_del || p->is_refskip || a.q < min_baseQ);
	a.indel = p->indel;
	a.b = a.hash = bam_seqi(seq, p->qpos);
	a.pos = pos;
	if (p->indel > 0)
		for (i = 0; i < p->indel; ++i)
			a.hash = (a.hash<<4) + a.hash + bam_seqi(seq, p->qpos + i + 1);
	return a;
}

static inline void print_allele(const bam_pileup1_t *p, int l_ref, const char *ref, int pos)
{
	const uint8_t *seq = bam_get_seq(p->b);
	int i;
	putchar(seq_nt16_str[bam_seqi(seq, p->qpos)]);
	if (p->indel > 0) {
		printf("+%d", p->indel);
		for (i = 1; i <= p->indel; ++i)
			putchar(seq_nt16_str[bam_seqi(seq, p->qpos + i)]);
	} else if (p->indel < 0) {
		printf("%d", p->indel);
		for (i = 1; i <= -p->indel; ++i)
			putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
	}
}

typedef struct {
	int tot_dp, max_dp, n_cnt, max_cnt;
	allele_t *a;
	int *cnt_b, *cnt_q;
} paux_t;

int main_pileup(int argc, char *argv[])
{
	int i, j, n, tid, beg, end, pos, last_tid, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, l_ref = 0, depth_only = 0, min_sum_q = 0;
	const bam_pileup1_t **plp;
	char *ref = 0, *reg = 0, *chr_end; // specified region
	faidx_t *fai = 0;
	bam_hdr_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	paux_t aux;
	bam_mplp_t mplp;

	// parse the command line
	while ((n = getopt(argc, argv, "r:q:Q:l:f:ds:")) >= 0) {
		if (n == 'f') fai = fai_load(optarg);
		else if (n == 'l') min_len = atoi(optarg); // minimum query length
		else if (n == 'r') reg = strdup(optarg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(optarg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(optarg);    // mapping quality threshold
		else if (n == 's') min_sum_q = atoi(optarg);
		else if (n == 'd') depth_only = 1;
	}
	if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   pileup [options] in1.bam [in2.bam [...]]\n\n");
        fprintf(stderr, "Options: -f FILE    reference genome [null]\n");
		fprintf(stderr, "         -l INT     minimum query length [%d]\n", min_len);
		fprintf(stderr, "         -q INT     minimum mapping quality [%d]\n", mapQ);
		fprintf(stderr, "         -Q INT     minimum base quality [%d]\n", baseQ);
		fprintf(stderr, "         -r STR     region [null]\n");
        fprintf(stderr, "\n");
		return 1;
	}

	// initialize the auxiliary data structures
	n = argc - optind; // the number of BAMs on the command line
	data = (aux_t**)calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	if (reg) {
		chr_end = (char*)hts_parse_reg(reg, &beg, &end);
		ref = fai? fai_fetch(fai, reg, &l_ref) : 0;
	} else chr_end = 0;

	// load the index or put the file position at the right place
	last_tid = -1;
	for (i = 0; i < n; ++i) {
		bam_hdr_t *htmp;
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = bgzf_open(argv[optind+i], "r"); // open BAM
		data[i]->min_mapQ = mapQ;                     // set the mapQ filter
		data[i]->min_len  = min_len;                  // set the qlen filter
		htmp = bam_hdr_read(data[i]->fp);             // read the BAM header
		if (i == 0 && chr_end) {
			char c = *chr_end;
			*chr_end = 0;
			last_tid = tid = bam_name2id(htmp, reg);
			*chr_end = c;
		}
		if (i) bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header
		else h = htmp; // keep the header of the 1st BAM
		if (tid >= 0) { // if a region is specified and parsed successfully
			hts_idx_t *idx = bam_index_load(argv[optind+i]); // load the index
			data[i]->itr = bam_itr_queryi(idx, tid, beg, end); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = (int*)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = (const bam_pileup1_t**)calloc(n, sizeof(const bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
	memset(&aux, 0, sizeof(paux_t));
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		for (i = aux.tot_dp = 0; i < n; ++i) aux.tot_dp += n_plp[i];
		if (aux.tot_dp == 0) continue; // well, this should not happen
		if (last_tid != tid && fai) {
			free(ref);
			ref = fai_fetch(fai, h->target_name[tid], &l_ref);
		}
		if (depth_only) { // only print read depth; no allele information
			fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
			if (ref == 0 || pos >= l_ref + beg) printf("\tN");
			else printf("\t%c", ref[pos - beg]);
			for (i = 0; i < n; ++i) { // base level filters have to go here
				int m = 0;
				for (j = 0; j < n_plp[i]; ++j) {
					const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
					if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
					else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
				}
				printf("\t%d", n_plp[i] - m); // this the depth to output
			}
		} else { // print alleles and allele counts
			int m = 0, n_alleles, k;
			allele_t *a;
			if (aux.tot_dp > aux.max_dp) {
				aux.max_dp = aux.tot_dp;
				kroundup32(aux.max_dp);
				aux.a = (allele_t*)realloc(aux.a, aux.max_dp * sizeof(allele_t));
			}
			a = aux.a;
			// collect alleles
			for (i = 0; i < n; ++i)
				for (j = 0; j < n_plp[i]; ++j) {
					a[m] = pileup2allele(&plp[i][j], baseQ, (uint64_t)i<<32 | j);
					if (!a[m].is_skip) ++m;
				}
			ks_introsort(allele, m, aux.a);
			// count alleles
			for (i = n_alleles = 1; i < m; ++i)
				if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash)
					++n_alleles;
			if (n_alleles == 0) continue;
			if (ref && pos - beg < l_ref && min_sum_q > 0) {
				int r = seq_nt16_table[(int)ref[pos - beg]], sum = 0;
				for (i = 0; i < m; ++i)
					if (a[i].indel != 0 || a[i].b != r)
						sum += a[i].q;
				if (sum < min_sum_q) continue;
			}
			// print
			fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
			if (ref == 0 || pos >= l_ref + beg) printf("\tN\t");
			else printf("\t%c\t", ref[pos - beg]);
			// print alleles
			print_allele(&plp[a[0].pos>>32][(uint32_t)a[0].pos], l_ref, ref, pos - beg);
			for (i = n_alleles = 1; i < m; ++i)
				if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) {
					putchar(',');
					print_allele(&plp[a[i].pos>>32][(uint32_t)a[i].pos], l_ref, ref, pos - beg);
					++n_alleles;
				}
			// collection per-BAM counts
			aux.n_cnt = n_alleles * n;
			if (aux.n_cnt > aux.max_cnt) {
				aux.max_cnt = aux.n_cnt;
				kroundup32(aux.max_cnt);
				aux.cnt_b = (int*)realloc(aux.cnt_b, aux.max_cnt * 2 * sizeof(int));
				aux.cnt_q = (int*)realloc(aux.cnt_q, aux.max_cnt * sizeof(int));
			}
			memset(aux.cnt_b, 0, aux.n_cnt * 2 * sizeof(int));
			memset(aux.cnt_q, 0, aux.n_cnt * sizeof(int));
			j = (a[0].pos>>32)*n_alleles;
			++aux.cnt_b[j]; aux.cnt_q[j] += a[0].q;
			for (i = 1, k = 0; i < m; ++i) {
				if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) ++k;
				j = (a[i].pos>>32)*n_alleles + k;
				++aux.cnt_b[j<<1|a[i].is_rev]; aux.cnt_q[j] += a[i].q;
			}
			// print counts
			for (i = m = 0; i < n; ++i, m += n_alleles) {
				putchar('\t');
				for (j = 0; j < n_alleles; ++j) {
					if (j) putchar(',');
					printf("%d:%d:%d", aux.cnt_b[(m+j)<<1], aux.cnt_b[(m+j)<<1|1], aux.cnt_q[m+j]);
				}
			}
		}
		putchar('\n');
		last_tid = tid;
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_hdr_destroy(h);
	for (i = 0; i < n; ++i) {
		bgzf_close(data[i]->fp);
		if (data[i]->itr) bam_itr_destroy(data[i]->itr);
		free(data[i]);
	}
	if (ref) free(ref);
	if (fai) fai_destroy(fai);
	free(aux.cnt_b); free(aux.cnt_q); free(aux.a);
	free(data); free(reg);
	return 0;
}
