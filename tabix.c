#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "tabix.h"

#define PACKAGE_VERSION "0.1.3-3 (r556)"

static int fetch_func(int l, const char *s, void *data)
{
	printf("%s\n", s);
	return 0;
}

int main(int argc, char *argv[])
{
	int c, skip = -1, meta = -1, list_chrms = 0, force = 0;
	ti_conf_t conf = ti_conf_gff;
	while ((c = getopt(argc, argv, "p:s:b:e:0S:c:lf")) >= 0) {
		switch (c) {
		case '0': conf.preset |= TI_FLAG_UCSC; break;
		case 'S': skip = atoi(optarg); break;
		case 'c': meta = optarg[0]; break;
		case 'p':
			if (strcmp(optarg, "gff") == 0) conf = ti_conf_gff;
			else if (strcmp(optarg, "bed") == 0) conf = ti_conf_bed;
			else if (strcmp(optarg, "sam") == 0) conf = ti_conf_sam;
			else if (strcmp(optarg, "vcf") == 0) conf = ti_conf_vcf;
			else if (strcmp(optarg, "psltbl") == 0) conf = ti_conf_psltbl;
			else {
				fprintf(stderr, "[main] unrecognized preset '%s'\n", optarg);
				return 1;
			}
			break;
		case 's': conf.sc = atoi(optarg); break;
		case 'b': conf.bc = atoi(optarg); break;
		case 'e': conf.ec = atoi(optarg); break;
        case 'l': list_chrms = 1; break;
		case 'f': force = 1; break;
		}
	}
	if (skip >= 0) conf.line_skip = skip;
	if (meta >= 0) conf.meta_char = meta;
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: tabix (TAB-delimited file InderXer)\n");
		fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
		fprintf(stderr, "Usage:   tabix <in.tab.bgz> [region1 [region2 [...]]]\n\n");
		fprintf(stderr, "Options: -p STR     preset: gff, bed, sam, vcf, psltbl [gff]\n");
		fprintf(stderr, "         -s INT     sequence name column [1]\n");
		fprintf(stderr, "         -b INT     start column [4]\n");
		fprintf(stderr, "         -e INT     end column [5]\n");
		fprintf(stderr, "         -S INT     skip first INT lines [0]\n");
		fprintf(stderr, "         -c CHAR    symbol for comment/meta lines [#]\n");
		fprintf(stderr, "         -0         zero-based coordinate\n");
		fprintf(stderr, "         -l         list chromosome names\n");
		fprintf(stderr, "         -f         force to overwrite the index\n");
		fprintf(stderr, "\n");
		return 1;
	}
    if (list_chrms)
        return ti_list_chromosomes(argv[optind]);
	if (optind + 1 == argc) {
		if (force == 0) {
			struct stat buf;
			char *fnidx = calloc(strlen(argv[optind]) + 5, 1);
			strcat(strcpy(fnidx, argv[optind]), ".tbi");
			if (stat(fnidx, &buf) == 0) {
				fprintf(stderr, "[tabix] the index file exists. Please use '-f' to overwrite.\n");
				free(fnidx);
				return 1;
			}
			free(fnidx);
		}
		return ti_index_build(argv[optind], &conf);
	}
	{ // retrieve
		BGZF *fp;
		fp = bgzf_open(argv[optind], "r");
		if (fp == 0) {
			fprintf(stderr, "[main] fail to open the data file.\n");
			return 1;
		}
		if (strcmp(argv[optind+1], ".") == 0) { // retrieve all
			kstring_t *str = calloc(1, sizeof(kstring_t));
			while (ti_readline(fp, str) >= 0) { // FIXME: check return code for error
				fputs(str->s, stdout); fputc('\n', stdout);
			}
			free(str->s); free(str);
		} else { // retrieve from specified regions
			ti_index_t *idx;
			int i;
			idx = ti_index_load(argv[optind]);
			if (idx == 0) {
				bgzf_close(fp);
				fprintf(stderr, "[main] fail to load the index.\n");
				return 1;
			}
			for (i = optind + 1; i < argc; ++i) {
				int tid, beg, end;
				if (ti_parse_region(idx, argv[i], &tid, &beg, &end) == 0) {
					ti_fetch(fp, idx, tid, beg, end, 0, fetch_func);
				} else fprintf(stderr, "[main] invalid region: unknown target name or minus interval.\n");
			}
			ti_index_destroy(idx);
		}
		bgzf_close(fp);
	}
	return 0;
}
