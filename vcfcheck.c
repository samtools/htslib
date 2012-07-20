#include <stdio.h>
#include <unistd.h>
#include "vcf.h"
#include "synced_bcf_reader.h"

typedef struct
{
	int n_snps, n_indels;
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

void do_indel_stats(args_t *args, stats_t *stats, bcf1_t *line)
{
	stats->n_indels++;
}

void do_snp_stats(args_t *args, stats_t *stats, bcf1_t *line)
{
	stats->n_snps++;
}

void check_vcf(args_t *args)
{
	int ret,i;
	readers_t *files = &args->files;
	while ( (ret=next_line(files)) )
	{
		bcf1_t *line = NULL;
		for (i=0; i<files->nreaders; i++)
		{
			if ( !(ret&1<<i) ) continue;
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
			do_indel_stats(args, &args->stats[ret-1], line);
		else
			do_snp_stats(args, &args->stats[ret-1], line);
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
	int id, nstats = args->files.nreaders==1 ? 1 : 3;
	printf("# Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
	for (id=0; id<nstats; id++)
	{
		printf("SN\t%d\tnumber of SNPs:\t%d\n", id, args->stats[id].n_snps);
		printf("SN\t%d\tnumber of indels:\t%d\n", id, args->stats[id].n_indels);
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

	print_header(args);
	check_vcf(args);
	print_stats(args);

	destroy_readers(&args->files);
	free(args);
	return 0;
}

