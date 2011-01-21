#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "tabix.h"

#define PACKAGE_VERSION "0.2.3 (r876)"

int main(int argc, char *argv[])
{
	int c, skip = -1, meta = -1, list_chrms = 0, force = 0, print_header = 0;
	ti_conf_t conf = ti_conf_gff;
	while ((c = getopt(argc, argv, "p:s:b:e:0S:c:lhf")) >= 0) {
		switch (c) {
		case '0': conf.preset |= TI_FLAG_UCSC; break;
		case 'S': skip = atoi(optarg); break;
		case 'c': meta = optarg[0]; break;
		case 'p':
			if (strcmp(optarg, "gff") == 0) conf = ti_conf_gff;
			else if (strcmp(optarg, "bed") == 0) conf = ti_conf_bed;
			else if (strcmp(optarg, "sam") == 0) conf = ti_conf_sam;
			else if (strcmp(optarg, "vcf") == 0 || strcmp(optarg, "vcf4") == 0) conf = ti_conf_vcf;
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
        case 'h': print_header = 1; break;
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
		fprintf(stderr, "         -h         print the VCF header\n");
		fprintf(stderr, "         -l         list chromosome names\n");
		fprintf(stderr, "         -f         force to overwrite the index\n");
		fprintf(stderr, "\n");
		return 1;
	}
    if (list_chrms) {
		ti_index_t *idx;
		int i, n;
		const char **names;
		idx = ti_index_load(argv[optind]);
		if (idx == 0) {
			fprintf(stderr, "[main] fail to load the index file.\n");
			return 1;
		}
		names = ti_seqname(idx, &n);
		for (i = 0; i < n; ++i) printf("%s\n", names[i]);
		free(names);
		ti_index_destroy(idx);
		return 0;
	}

    struct stat stat_tbi,stat_vcf;
    char *fnidx = calloc(strlen(argv[optind]) + 5, 1);
    strcat(strcpy(fnidx, argv[optind]), ".tbi");

	if (optind + 1 == argc) {
		if (force == 0) {
			if (stat(fnidx, &stat_tbi) == 0) 
            {
                // Before complaining, check if the VCF file isn't newer. This is a common source of errors,
                //  people tend not to notice that tabix failed
                stat(argv[optind], &stat_vcf);
                if ( stat_vcf.st_mtime <= stat_tbi.st_mtime )
                {
                    fprintf(stderr, "[tabix] the index file exists. Please use '-f' to overwrite.\n");
                    free(fnidx);
                    return 1;
                }
			}
		}
        if ( is_bgzipped(argv[optind])!=1 )
        {
            fprintf(stderr,"[tabix] was bgzip used to compress this file? %s\n", argv[optind]);
            free(fnidx);
            return 1;
        }
		return ti_index_build(argv[optind], &conf);
	}
	{ // retrieve
        // Common source of errors: new VCF is used with an old index
        stat(fnidx, &stat_tbi);
        stat(argv[optind], &stat_vcf);
        if ( force==0 && stat_vcf.st_mtime >= stat_tbi.st_mtime )
        {
            fprintf(stderr, "[tabix] the index file is older than the vcf file. Please use '-f' to overwrite or reindex.\n");
            free(fnidx);
            return 1;
        }
        free(fnidx);

		tabix_t *t;
		if ((t = ti_open(argv[optind], 0)) == 0) {
			fprintf(stderr, "[main] fail to open the data file.\n");
			return 1;
		}
		if (strcmp(argv[optind+1], ".") == 0) { // retrieve all
			ti_iter_t iter;
			const char *s;
			int len;
			iter = ti_query(t, 0, 0, 0);
			while ((s = ti_read(t, iter, &len)) != 0) {
				fputs(s, stdout); fputc('\n', stdout);
			}
			ti_iter_destroy(iter);
		} else { // retrieve from specified regions
			int i;
			if ( ti_lazy_index_load(t) ) 
            {
                fprintf(stderr,"[tabix] failed to load the index file.\n");
                return 1;
            }

            ti_iter_t iter;
            const char *s;
            int len;
            if ( print_header )
            {
                // If requested, print the header lines here
                iter = ti_query(t, 0, 0, 0);
                while ((s = ti_read(t, iter, &len)) != 0) {
                    if ( *s != '#' ) break;
                    fputs(s, stdout); fputc('\n', stdout);
                }
                ti_iter_destroy(iter);
            }
			for (i = optind + 1; i < argc; ++i) {
				int tid, beg, end;
				if (ti_parse_region(t->idx, argv[i], &tid, &beg, &end) == 0) {
					iter = ti_queryi(t, tid, beg, end);
					while ((s = ti_read(t, iter, &len)) != 0) {
						fputs(s, stdout); fputc('\n', stdout);
					}
					ti_iter_destroy(iter);
				} 
                // else fprintf(stderr, "[main] invalid region: unknown target name or minus interval.\n");
			}
		}
		ti_close(t);
	}
	return 0;
}
