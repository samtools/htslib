#include <stdio.h>

#include "htslib/sam.h"

int ntests = 0;
int nfailures = 0;

void check(const bam1_t *aln, const char *testname, const char *tag, int value)
{
    int32_t refvalue;
    uint8_t *aux = bam_aux_get(aln, tag);
    if (!aux) return;
    ntests++;
    refvalue = bam_aux2i(aux);
    if (value != refvalue) {
        fprintf(stderr, "%s FAIL for %s: computed %d != %d expected\n",
                testname, bam_get_qname(aln), value, refvalue);
        nfailures++;
    }
}

int main(int argc, char **argv)
{
    bam_hdr_t *header;
    bam1_t *aln = bam_init1();
	int i;

	for (i = 1; i < argc; i++) {
		samFile *in = sam_open(argv[i], "r");
		if (in == NULL) { perror(argv[1]); return 1; }

		header = sam_hdr_read(in);
		while (sam_read1(in, header, aln) >= 0) {
			check(aln, "cigar2qlen", "XQ",
				  bam_cigar2qlen(aln->core.n_cigar, bam_get_cigar(aln)));
			check(aln, "cigar2rlen", "XR",
				  bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)));
			check(aln, "endpos", "XE", bam_endpos(aln));
		}

		bam_hdr_destroy(header);
		sam_close(in);
	}

    bam_destroy1(aln);

    return (nfailures > 0);
}
