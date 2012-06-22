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
		fprintf(stderr, "\nUsage:   tabix [options] <in.gz> [reg1 [...]]\n\n");
		fprintf(stderr, "Options: -p STR    preset: gff, bed, sam or vcf [gff]\n");
		fprintf(stderr, "         -s INT    column number for sequence names (suppressed by -p) [1]\n");
		fprintf(stderr, "         -b INT    column number for region start [4]\n");
		fprintf(stderr, "         -e INT    column number for region end (if no end, set INT to -b) [5]\n");
		fprintf(stderr, "         -0        specify coordinates are zero-based\n");
		fprintf(stderr, "         -S INT    skip first INT lines [0]\n");
		fprintf(stderr, "         -c CHAR   skip lines starting with CHAR [null]\n");
		fprintf(stderr, "         -a        print all records\n");
		fprintf(stderr, "         -f        force to overwrite existing index\n");
		fprintf(stderr, "         -m INT    set the minimal interval size to 1<<INT; 0 for the old tabix index [0]\n");
		fprintf(stderr, "\n");
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
		tbx_index_build(argv[optind], min_shift, &conf);
	} else { // read with random access
		tbx_t *tbx;
		BGZF *fp;
		kstring_t s;
		int i;
		if ((tbx = tbx_index_load(argv[optind])) == 0) return 1;
		if ((fp = bgzf_open(argv[optind], "r")) == 0) return 1;
		s.s = 0; s.l = s.m = 0;
		for (i = optind + 1; i < argc; ++i) {
			hts_itr_t *itr;
			if ((itr = tbx_itr_querys(tbx, argv[i])) == 0) continue;
			while (tbx_itr_next(fp, tbx, itr, &s) >= 0) puts(s.s);
			tbx_itr_destroy(itr);
		}
		free(s.s);
		bgzf_close(fp);
		tbx_destroy(tbx);
	}
	return 0;
}
