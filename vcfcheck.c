#include <stdio.h>
#include <unistd.h>
#include "vcf.h"
#include "synced_bcf_reader.h"

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
	bcf_unpack(line, BCF_UN_INFO);
	int an_id = bcf_id2int(reader->header, BCF_DT_ID, "AN");
	int ac_id = bcf_id2int(reader->header, BCF_DT_ID, "AC");
	int i, an=0, ac_len=0;
	uint8_t *ac=NULL;
	for (i=0; i<line->n_info; i++)
	{
		bcf_info_t *z = &line->d.info[i];
		if ( z->key == an_id ) an = z->v1.i;
		else if ( z->key == ac_id ) { ac = z->vptr; ac_len = z->len; }
	}

	// Ts/Tv
	int ref = acgt2int(*line->d.allele[0]);
	if ( ref>=0 )
	{
		for (i=1; i<line->n_allele; i++)
		{
			int alt = acgt2int(*line->d.allele[i]);
			if ( alt<0 ) continue;

			int iaf = -1;
			if ( i<=ac_len ) 
			{
				iaf = ac[i-1]*(stats->n_af-1)/an;
				if ( iaf>=stats->n_af ) iaf = stats->n_af-1;
			}
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
		int is_complex = line->rlen-1;
		if ( !is_complex )
		{
			for (i=1; i<line->n_allele; i++)
			{
				if ( strlen(line->d.allele[i])==1 ) continue;
				is_complex = 1; break;
			}
		}
		if ( is_complex ) 
			do_indel_stats(args, &args->stats[ret-1], reader);
		else
			do_snp_stats(args, &args->stats[ret-1], reader);
	}
}

void print_header(args_t *args)
{
	int i;
	printf("# This file was produced by vcfcheck.\n");
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
		printf("ID\t2\t%s\t%s\n", args->files.readers[0].fname,args->files.readers[0].fname);
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
	printf("# Ts/Tv by non-reference allele frequency:\n# TsTvAF\t[2]id\t[3]frequency\t[4]count\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=0; i<stats->n_af; i++)
		{
			if ( !stats->n_ts_af[i] || !stats->n_tv_af[i] ) continue;
			printf("TsTvAF\t%d\t%.2f\t%.2f\n", id,(float)i/stats->n_af, (float)stats->n_ts_af[i]/stats->n_tv_af[i]);
		}
	}
}

void usage(void)
{
	fprintf(stderr, "\nAbout:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfcheck.\n");
	fprintf(stderr, "         When two files are given, the program generates separate stats for intersection\n");
	fprintf(stderr, "         and the complements.\n");
	fprintf(stderr, "Usage:   vcfcheck [options] <A.vcf.gz> [<B.vcf.gz>]\n\n");
	// fprintf(stderr, "Options: -b           output in BCF\n");
	// fprintf(stderr, "         -S           input is VCF\n");
	// fprintf(stderr, "         -o FILE      output file name [stdout]\n");
	// fprintf(stderr, "         -l INT       compression level [%d]\n", clevel);
	// fprintf(stderr, "         -t FILE      list of reference names and lengths [null]\n");
	// fprintf(stderr, "         -s FILE/STR  list of samples (STR if started with ':'; FILE otherwise) [null]\n");
	// fprintf(stderr, "         -G           drop individual genotype information\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfcheck(int argc, char *argv[])
{
	int c, clevel, flag, n_samples=-1;
	char *fn_ref = 0, **samples = 0;
	
	while ((c = getopt(argc, argv, "l:bSt:o:T:s:G")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 'G': n_samples = 0; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 's': samples = hts_readlines(optarg, &n_samples); break;
		}
	}
	if (argc == optind) usage();

	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;
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

