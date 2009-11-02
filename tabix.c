#include "bgzf.h"
#include "tabix.h"

#define PACKAGE_VERSION "0.0.0-1 (r500)"

static int fetch_func(int l, const char *s, void *data)
{
	printf("%s\n", s);
	return 0;
}

int main(int argc, char *argv[])
{
	int c;
	ti_conf_t conf = ti_conf_gff;
	while ((c = getopt(argc, argv, "p:s:b:e:")) >= 0) {
		switch (c) {
		case 'p':
			if (strcmp(optarg, "gff") == 0) conf = ti_conf_gff;
			else {
				fprintf(stderr, "[main] unrecognized preset '%s'\n", optarg);
				return 1;
			}
			break;
		case 's': conf.sc = atoi(optarg); break;
		case 'b': conf.bc = atoi(optarg); break;
		case 'e': conf.ec = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: tabix (TAB-delimited file InderXer)\n");
		fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
		fprintf(stderr, "Usage:   tabix <in.tab.bgz> [region1 [region2 [...]]]\n\n");
		fprintf(stderr, "Options: -p STR     preset: gff, bed, sam, vcf [gff]\n");
		fprintf(stderr, "         -s INT     sequence name column [1]\n");
		fprintf(stderr, "         -b INT     start column [4]\n");
		fprintf(stderr, "         -e INT     end column [5]\n\n");
		return 1;
	}
	if (optind + 1 == argc)
		return ti_index_build(argv[optind], &conf);
	{ // retrieve
		BGZF *fp;
		fp = bgzf_open(argv[optind], "r");
		if (fp == 0) {
			fprintf(stderr, "[main] fail to open the data file.\n");
			return 1;
		}
		if (strcmp(argv[optind+1], ".") == 0) { // retrieve all
			kstring_t *str = calloc(1, sizeof(kstring_t));
			while (ti_readline(fp, str) >= 0) // FIXME: check return code for error
				printf("%s\n", str->s);
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
