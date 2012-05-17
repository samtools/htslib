#include <stdio.h>
#include <unistd.h>
#include "vcf.h"

int main_vcfview(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8];
	bcf_hdr_t *h;
	htsFile *in;
	bcf1_t *v;

	while ((c = getopt(argc, argv, "l:bSt:o:T:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: vcfview [-bS] [-t ref.fai] [-l level] <in.bcf>\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	h = vcf_hdr_read(in);
	v = bcf_init1();

	{
		htsFile *out;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		out = hts_open(fn_out? fn_out : "-", modew, 0);
		vcf_hdr_write(out, h);
		while (vcf_read1(in, h, v) >= 0) vcf_write1(out, h, v);
		hts_close(out);
	}

	bcf_destroy1(v);
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}

