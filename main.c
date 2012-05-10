#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "vcf.h"

int main(int argc, char *argv[])
{
	int task = 0; // 0 for conversion, 1 for counting and 2 for site frequency
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8];
	vcf_hdr_t *h;
	vcfFile *in;
	vcf1_t *v;

	while ((c = getopt(argc, argv, "l:bSt:o:T:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		case 'T':
			if (strcmp(optarg, "count") == 0) task = 1;
			else if (strcmp(optarg, "freq") == 0) task = 2;
			break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: bcf2ls [-bS] [-t ref.fai] [-l level] [-T count|freq] <in.bcf>\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = vcf_open(argv[optind], moder, fn_ref);
	h = vcf_hdr_read(in);
	v = vcf_init1();

	if (task == 0) {
		vcfFile *out;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		out = vcf_open(fn_out? fn_out : "-", modew, 0);
		vcf_hdr_write(out, h);
		while (vcf_read1(in, h, v) >= 0) vcf_write1(out, h, v);
		vcf_close(out);
	} else if (task == 1) {
		int64_t cnt = 0;
		while (vcf_read1(in, h, v) >= 0) ++cnt;
		printf("%ld\n", (long)cnt);
	} else if (task == 2) { // FIXME: not working for >=10 alleles
		int gt, *n_allele;
		gt = vcf_id2int(h, VCF_DT_ID, "GT");
		n_allele = (int*)alloca(64 * sizeof(int));
		while (vcf_read1(in, h, v) >= 0) {
			int i, j, l;
			vcf_fmt_t *fmt;
			for (i = 0; i < 10; ++i) n_allele[i] = 0;
			fmt = vcf_unpack_fmt(h, v);
			for (i = 0; i < v->n_fmt; ++i)
				if (fmt[i].id == gt) break;
			if (i != v->n_fmt) { // has GT
				int8_t *p = (int8_t*)fmt[i].p;
				for (j = 0; j < v->n_sample; ++j, p += fmt[i].n)
					for (l = 0; l < fmt[i].n; ++l)
						if (p[l]>>1) ++n_allele[(p[l]>>1)-1];
				printf("%s\t%d", h->id[VCF_DT_CTG][v->rid].key, v->pos + 1);
				for (i = 0; i < v->n_allele; ++i) printf("\t%d", n_allele[i]);
				putchar('\n');
			}
			free(fmt);
		}
	}

	vcf_destroy1(v);
	vcf_hdr_destroy(h);
	vcf_close(in);
	return 0;
}
