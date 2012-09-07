#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
	readers_t files;
	char **argv;
	int argc;
}
args_t;

static void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(-1);
}

void merge_vcf(args_t *args)
{
	error("To be implemented.");
}

static void usage(void)
{
	fprintf(stderr, "Usage:   vcfmerge [options] <A.vcf.gz> [<B.vcf.gz> [<C.vcf.gz> [...]]]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -c, --collapse <string>           merge sites with differing alleles for <snps|indels|both|any>\n");
	fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    collect statistics in the given region only\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfmerge(int argc, char *argv[])
{
	int c;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"collapse",0,0,'c'},
		{"apply-filters",0,0,'f'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hc:fr:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'c':
				if ( !strcmp(optarg,"snps") ) args->files.collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) args->files.collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) args->files.collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) args->files.collapse |= COLLAPSE_ANY;
				break;
			case 'f': args->files.apply_filters = 1; break;
			case 'r': args->files.region = optarg; break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if (argc == optind) usage();

	if ( argc-optind>2 ) usage();
	while (optind<argc)
	{
		if ( !add_reader(argv[optind], &args->files) ) error("Could not load the index: %s\n", argv[optind]);
		optind++;
	}

	merge_vcf(args);
	destroy_readers(&args->files);
	free(args);
	return 0;
}

