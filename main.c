#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "vcf.h"

int main(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, mode[8];
	vcf_hdr_t *h;
	vcfFile *in, *out;
	vcf1_t *v;

	while ((c = getopt(argc, argv, "l:bSD:o:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 'D': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: bcf2ls [-bS] [-D ref.fai] [-l level] <in.bcf>\n");
		return 1;
	}
	
	strcpy(mode, "r");
	if (clevel >= 0 && clevel <= 9) sprintf(mode + 1, "%d", clevel);
	if ((flag&1) == 0) strcat(mode, "b");

	in = vcf_open(argv[optind], mode, fn_ref);
	out = vcf_open(fn_out? fn_out : "-", (flag&2)? "wb" : "w", 0);
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
