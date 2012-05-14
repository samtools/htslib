#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"

int main_samview(int argc, char *argv[])
{
	htsFile *in;
	char *fn_ref = 0;
	int flag = 0, c, clevel = -1;
	char moder[8];
	sam_hdr_t *h;
	sam1_t *b;

	while ((c = getopt(argc, argv, "bSl:t:")) >= 0) {
		switch (c) {
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 't': fn_ref = optarg; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samview [-bS] [-t ref.fai] [-l level] <in.bam>|<in.sam> [region]\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	h = sam_hdr_read(in);
	b = sam_init1();

	if ((flag&4) == 0) { // SAM/BAM output
		htsFile *out;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		out = hts_open("-", modew, 0);
		sam_hdr_write(out, h);
		if (optind + 2 >= argc && !(flag&1)) { // BAM input and has a region
			int i;
			sam_idx_t *idx;
			char *tmp;
			tmp = (char*)calloc(1, strlen(argv[optind]) + 5);
			strcat(strcpy(tmp, argv[optind]), ".bai");
			if ((idx = sam_index_load_local(tmp)) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
				return 1;
			}
			for (i = optind + 1; i < argc; ++i) {
				hts_iter_t *iter;
				if ((iter = sam_iter_querys(idx, h, argv[i])) == 0) {
					fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
					continue;
				}
				while (sam_iter_read(in, iter, b) >= 0) sam_write1(out, h, b);
				hts_iter_destroy(iter);
			}
			hts_idx_destroy(idx);
			free(tmp);
		} else {
			while (sam_read1(in, h, b) >= 0) sam_write1(out, h, b);
		}
		hts_close(out);
	}

	sam_destroy1(b);
	sam_hdr_destroy(h);
	hts_close(in);
	return 0;
}
