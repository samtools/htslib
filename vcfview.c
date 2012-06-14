#include <stdio.h>
#include <unistd.h>
#include "vcf.h"

int main_vcfview(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8];
	bcf_hdr_t *h;
	htsFile *in;
	bcf1_t *b;

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
		fprintf(stderr, "\nUsage:   vcfview [-bS] [-t ref.fai] [-l level] <in.bcf>|<in.vcf>\n\n");
		fprintf(stderr, "Options: -b        output in BCF\n");
		fprintf(stderr, "         -S        input is VCF\n");
		fprintf(stderr, "         -o FILE   output file name [stdout]\n");
		fprintf(stderr, "         -l INT    compression level [%d]\n", clevel);
		fprintf(stderr, "         -t FILE   list of reference names and lengths [null]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	h = vcf_hdr_read(in);
	b = bcf_init1();

	if ((flag&4) == 0) { // VCF/BCF output
		htsFile *out;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		out = hts_open(fn_out? fn_out : "-", modew, 0);
		vcf_hdr_write(out, h);
		if (optind + 1 < argc && !(flag&1)) { // BAM input and has a region
			int i;
			hts_idx_t *idx;
			if ((idx = bcf_index_load(argv[optind])) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BCF index\n", __func__);
				return 1;
			}
			for (i = optind + 1; i < argc; ++i) {
				hts_iter_t *iter;
				if ((iter = bcf_iter_querys(idx, h, argv[i])) == 0) {
					fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
					continue;
				}
				while (bcf_iter_read((BGZF*)in->fp, iter, b) >= 0) vcf_write1(out, h, b);
				hts_iter_destroy(iter);
			}
			hts_idx_destroy(idx);
		} else while (vcf_read1(in, h, b) >= 0) vcf_write1(out, h, b);
		hts_close(out);
	}

	bcf_destroy1(b);
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}

