// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "sam.h"
#include "faidx.h"

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

int main_pileup(int argc, char *argv[])
{
	int i, n, tid, beg, end, pos, last_tid, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, l_ref = 0;
	const bam_pileup1_t **plp;
	char *ref = 0, *reg = 0, *chr_end; // specified region
	faidx_t *fai = 0;
	bam_hdr_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;

	// parse the command line
	while ((n = getopt(argc, argv, "r:q:Q:l:f:")) >= 0) {
		if (n == 'f') fai = fai_load(optarg);
		else if (n == 'l') min_len = atoi(optarg); // minimum query length
		else if (n == 'r') reg = strdup(optarg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(optarg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(optarg);    // mapping quality threshold
	}
	if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   pileup [options] in1.bam [in2.bam [...]]\n");
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
		data[i] = calloc(1, sizeof(aux_t));
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
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		if (last_tid != tid && fai) {
			free(ref);
			ref = fai_fetch(fai, h->target_name[tid], &l_ref);
		}
		if (ref == 0 || pos >= l_ref + beg) printf("\tN");
		else printf("\t%c", ref[pos - beg]);
		for (i = 0; i < n; ++i) { // base level filters have to go here
			int j, m = 0;
			for (j = 0; j < n_plp[i]; ++j) {
				const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
				else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
			}
			printf("\t%d", n_plp[i] - m); // this the depth to output
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
	free(data); free(reg);
	return 0;
}
