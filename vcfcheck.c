#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
	int n_snps, n_indels;
	int n_ts, n_tv;
	int n_af, *n_ts_af, *n_tv_af;
}
stats_t;

typedef struct
{
	stats_t stats[3];
	readers_t files;
	char **argv;
	int argc;
}
args_t;

void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(-1);
}

void init_stats(args_t *args)
{
	int id,nstats = args->files.nreaders==1 ? 1 : 3;
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		stats->n_af = 20;
		stats->n_ts_af = (int*) calloc(stats->n_af,sizeof(int));
		stats->n_tv_af = (int*) calloc(stats->n_af,sizeof(int));
	}
}
void destroy_stats(args_t *args)
{
	int id,nstats = args->files.nreaders==1 ? 1 : 3;
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		free(stats->n_ts_af);
		free(stats->n_tv_af);
	}
}

void do_indel_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_indels++;
}

inline int acgt2int(char c)
{
	if ( (int)c>96 ) c -= 32;
	if ( c=='A' ) return 0;
	if ( c=='C' ) return 1;
	if ( c=='G' ) return 2;
	if ( c=='T' ) return 3;
	return -1;
}

void do_snp_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_snps++;

	// AF
	bcf1_t *line = reader->line;
	int *ac = (int*) calloc(line->n_allele,sizeof(int));
	calc_ac(reader->header, line, ac, BCF_UN_FMT);
	int i, an=0;
	for (i=0; i<line->n_allele; i++)
		an += ac[i];

	// Ts/Tv
	int ref = acgt2int(*line->d.allele[0]);
	if ( ref>=0 )
	{
		for (i=1; i<line->n_allele; i++)
		{
			if ( !(line->d.var[i].type&VCF_SNP) ) continue;
			int alt = acgt2int(*line->d.allele[i]);
			if ( alt<0 ) continue;

			int iaf = ac[i]*(stats->n_af-1)/an;
			if ( iaf>=stats->n_af ) iaf = stats->n_af-1;

			if ( abs(ref-alt)==2 ) 
			{
				stats->n_ts++;
				if ( iaf>=0 ) stats->n_ts_af[iaf]++;
			}
			else 
			{
				stats->n_tv++;
				if ( iaf>=0 ) stats->n_tv_af[iaf]++;
			}
		}
	}

	free(ac);
}

void check_vcf(args_t *args)
{
	int ret,i;
	readers_t *files = &args->files;
	while ( (ret=next_line(files)) )
	{
		reader_t *reader = NULL;
		bcf1_t *line = NULL;
		for (i=0; i<files->nreaders; i++)
		{
			if ( !(ret&1<<i) ) continue;
			reader = &files->readers[i];
			line = files->readers[i].line;
			break;
		}
		set_variant_types(line);
		if ( line->d.var_type&VCF_SNP ) 
			do_snp_stats(args, &args->stats[ret-1], reader);
		if ( line->d.var_type&VCF_INDEL )
			do_indel_stats(args, &args->stats[ret-1], reader);
	}
}

void print_header(args_t *args)
{
	int i;
	printf("# This file was produced by vcfcheck and can be plotted using plot-vcfcheck.\n");
	printf("# The command line was:\thtscmd %s ", args->argv[0]);
	for (i=1; i<args->argc; i++)
		printf(" %s",args->argv[i]);
	printf("\n#\n");

	printf("# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names\n");
	if ( args->files.nreaders==1 )
		printf("ID\t0\t%s\n", args->files.readers[0].fname);
	else
	{
		printf("ID\t0\t%s\n", args->files.readers[0].fname);
		printf("ID\t1\t%s\n", args->files.readers[1].fname);
		printf("ID\t2\t%s\t%s\n", args->files.readers[0].fname,args->files.readers[1].fname);
	}
}

void print_stats(args_t *args)
{
	int i, id, nstats = args->files.nreaders==1 ? 1 : 3;
	printf("# Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		printf("SN\t%d\tnumber of SNPs:\t%d\n", id, stats->n_snps);
		printf("SN\t%d\tnumber of indels:\t%d\n", id, stats->n_indels);
		printf("SN\t%d\tts/tv:\t%.2f\n", id, stats->n_tv?(float)stats->n_ts/stats->n_tv:(float)0);
	}
	printf("# Ts/Tv by non-reference allele frequency:\n# TsTvAF\t[2]id\t[3]AF\t[4]number of SNPs\t[5]Ts/Tv\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=0; i<stats->n_af; i++)
		{
			int ts = stats->n_ts_af[i];
			int tv = stats->n_tv_af[i];
			if ( !ts && !tv ) continue;
			printf("TsTvAF\t%d\t%.2f\t%d\t%.2f\n", id,(float)i/stats->n_af,ts+tv,tv?(float)ts/tv:0);
		}
	}
}

void usage(void)
{
	fprintf(stderr, "\nAbout:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfcheck.\n");
	fprintf(stderr, "         When two files are given, the program generates separate stats for intersection\n");
	fprintf(stderr, "         and the complements.\n");
	fprintf(stderr, "Usage:   vcfcheck [options] <A.vcf.gz> [<B.vcf.gz>]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -c, --collapse <string>       treat sites with differing alleles as same for <snps|indels|both|any>\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfcheck(int argc, char *argv[])
{
	int c;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"collapse",0,0,'c'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hc:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'c':
				if ( !strcmp(optarg,"snps") ) args->files.collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) args->files.collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) args->files.collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) args->files.collapse |= COLLAPSE_ANY;
				break;
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

	init_stats(args);
	print_header(args);
	check_vcf(args);
	print_stats(args);
	destroy_stats(args);
	destroy_readers(&args->files);
	free(args);
	return 0;
}

