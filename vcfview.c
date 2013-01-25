#include <stdio.h>
#include <unistd.h>
#include "vcf.h"

int main_vcfview(int argc, char *argv[])
{
	int i, c, clevel = -1, flag = 0, n_samples = -1, *imap = 0, excl_indel = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8], **samples = 0;
	bcf_hdr_t *h, *hsub = 0;
	htsFile *in;
	bcf1_t *b;

	while ((c = getopt(argc, argv, "l:bSt:o:T:s:GI")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 'G': n_samples = 0; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		case 's': samples = hts_readlines(optarg, &n_samples); break;
		case 'I': excl_indel = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\nUsage:   vcfview [options] <in.bcf>|<in.vcf>|<in.vcf.gz>\n\n");
		fprintf(stderr, "Options: -b           output in BCF\n");
		fprintf(stderr, "         -S           input is VCF\n");
		fprintf(stderr, "         -o FILE      output file name [stdout]\n");
		fprintf(stderr, "         -l INT       compression level [%d]\n", clevel);
		fprintf(stderr, "         -t FILE      list of reference names and lengths [null]\n");
		fprintf(stderr, "         -s FILE/STR  list of samples (STR if started with ':'; FILE otherwise) [null]\n");
		fprintf(stderr, "         -G           drop individual genotype information\n");
		fprintf(stderr, "         -I           exclude INDELs\n");
		fprintf(stderr, "\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0 && !(file_type(argv[optind])&(IS_VCF|IS_VCF_GZ))) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	h = vcf_hdr_read(in);
	if (h == 0) {
		fprintf(stderr, "[E::%s] fail to read the VCF/BCF2 header\n", __func__);
		hts_close(in);
		return 1;
	}
	if (n_samples >= 0) {
		if (n_samples) imap = (int*)malloc(n_samples * sizeof(int));
		hsub = bcf_hdr_subset(h, n_samples, samples, imap);
	}
	b = bcf_init1();

	if ((flag&4) == 0) { // VCF/BCF output
		htsFile *out;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		out = hts_open(fn_out? fn_out : "-", modew, 0);
		vcf_hdr_write(out, hsub? hsub : h);
		if (optind + 1 < argc && !(flag&1)) { // BAM input and has a region
			hts_idx_t *idx;
			if ((idx = bcf_index_load(argv[optind])) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BCF index\n", __func__);
				return 1;
			}
			for (i = optind + 1; i < argc; ++i) {
				hts_itr_t *iter;
				if ((iter = bcf_itr_querys(idx, h, argv[i])) == 0) {
					fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
					continue;
				}
				while (bcf_itr_next((BGZF*)in->fp, iter, b) >= 0) {
					if (excl_indel && !bcf_is_snp(b)) continue;
					if (n_samples >= 0) {
						bcf_subset(h, b, n_samples, imap);
						vcf_write1(out, hsub, b);
					} else vcf_write1(out, h, b);
				}
				hts_itr_destroy(iter);
			}
			hts_idx_destroy(idx);
		} else {
			while (vcf_read1(in, h, b) >= 0) {
				if (excl_indel && !bcf_is_snp(b)) continue;
				if (n_samples >= 0) {
					bcf_subset(h, b, n_samples, imap);
					vcf_write1(out, hsub, b);
				} else vcf_write1(out, h, b);
			}
		}
		hts_close(out);
	}

	bcf_destroy1(b);
	if (n_samples > 0) {
		for (i = 0; i < n_samples; ++i) free(samples[i]);
		free(samples);
		bcf_hdr_destroy(hsub);
		free(imap);
	}
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}

