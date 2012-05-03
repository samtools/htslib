#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "vcf.h"

int main(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8], modew[8];
	vcf_hdr_t *h;
	vcfFile *in, *out;
	vcf1_t *v;

	while ((c = getopt(argc, argv, "l:bSt:o:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: bcf2ls [-bS] [-t ref.fai] [-l level] <in.bcf>\n");
		return 1;
	}
	
	strcpy(moder, "r"); strcpy(modew, "w");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
	if ((flag&1) == 0) strcat(moder, "b");
	if (flag&2) strcat(modew, "b");

	in = vcf_open(argv[optind], moder, fn_ref);
	out = vcf_open(fn_out? fn_out : "-", modew, 0);
	h = vcf_hdr_read(in);
	vcf_hdr_write(out, h);
	v = vcf_init1();
	while (vcf_read1(in, h, v) >= 0) vcf_write1(out, h, v);
	vcf_destroy1(v);
	vcf_hdr_destroy(h);
	vcf_close(out);
	vcf_close(in);
	return 0;
}
