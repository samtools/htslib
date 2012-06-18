#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "tbx.h"

int main_tabix(int argc, char *argv[])
{
	int c, min_shift = -1, is_force = 0, is_all = 0;
	tbx_conf_t conf = tbx_conf_gff;
	while ((c = getopt(argc, argv, "0fap:s:b:e:S:c:m:")) >= 0)
		if (c == '0') conf.preset |= TBX_UCSC;
		else if (c == 'f') is_force = 1;
		else if (c == 'a') is_all = 1;
		else if (c == 'm') min_shift = atoi(optarg);
		else if (c == 's') conf.sc = atoi(optarg);
		else if (c == 'b') conf.bc = atoi(optarg);
		else if (c == 'e') conf.ec = atoi(optarg);
		else if (c == 'c') conf.meta_char = *optarg;
		else if (c == 'S') conf.line_skip = atoi(optarg);
		else if (c == 'p') {
			if (strcmp(optarg, "gff") == 0) conf = tbx_conf_gff;
			else if (strcmp(optarg, "bed") == 0) conf = tbx_conf_bed;
			else if (strcmp(optarg, "sam") == 0) conf = tbx_conf_sam;
			else if (strcmp(optarg, "vcf") == 0) conf = tbx_conf_vcf;
		}
	if (optind == argc) {
		fprintf(stderr, "Usage: tabix [options] <in.gz>\n");
		return 1;
	}
	if (is_all) { // read without random access
		kstring_t s;
		BGZF *fp;
		s.l = s.m = 0; s.s = 0;
		fp = bgzf_open(argv[optind], "r");
		while (bgzf_getline(fp, '\n', &s) >= 0) puts(s.s);
		bgzf_close(fp);
		free(s.s);
	} else if (optind + 2 > argc) { // create index
		if (!is_force) {
			char *fn;
			FILE *fp;
			fn = (char*)alloca(strlen(argv[optind]) + 5);
			strcat(strcpy(fn, argv[optind]), min_shift <= 0? ".tbi" : ".csi");
			if ((fp = fopen(fn, "rb")) != 0) {
				fclose(fp);
				fprintf(stderr, "[E::%s] the index file exists; use option '-f' to overwrite\n", __func__);
				return 1;
			}
		}
		tbx_index_build(argv[optind], 0, min_shift, &conf);
	} else { // read with random access
	}
	return 0;
}
