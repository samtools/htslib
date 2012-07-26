#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
	int n_snps, n_indels, n_mals;
	int *ac_ts, *ac_tv, *ac_snps, *ac_indels, m_ac;
	int *insertions, *deletions, m_indel;	// maximum indel length
	int subst[15];
}
stats_t;

typedef struct
{
	stats_t stats[3];
	int *tmp_ac, ntmp_ac;
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

inline int acgt2int(char c)
{
	if ( (int)c>96 ) c -= 32;
	if ( c=='A' ) return 0;
	if ( c=='C' ) return 1;
	if ( c=='G' ) return 2;
	if ( c=='T' ) return 3;
	return -1;
}
#define int2acgt(i) "ACGT"[i]

void init_stats(args_t *args)
{
	int i,nstats = args->files.nreaders==1 ? 1 : 3;
	for (i=0; i<nstats; i++)
	{
		stats_t *stats = &args->stats[i];
		stats->m_indel    = 60;
		stats->insertions = (int*) calloc(stats->m_indel,sizeof(int));
		stats->deletions  = (int*) calloc(stats->m_indel,sizeof(int));
	}
}
int *realloc_stat(int *array, int ori, int new)
{
	array = (int*) realloc(array, new*sizeof(int));
	int i;
	for (i=ori; i<new; i++) array[i] = 0;
	return array;
}
void realloc_stats(args_t *args, int ac)
{
	int i,nstats = args->files.nreaders==1 ? 1 : 3;
	for (i=0; i<nstats; i++)
	{
		stats_t *stats = &args->stats[i];
		stats->ac_ts     = (int*)realloc_stat(stats->ac_ts,stats->m_ac,ac);
		stats->ac_tv     = (int*)realloc_stat(stats->ac_tv,stats->m_ac,ac);
		stats->ac_snps   = (int*)realloc_stat(stats->ac_snps,stats->m_ac,ac);
		stats->ac_indels = (int*)realloc_stat(stats->ac_indels,stats->m_ac,ac);
		stats->m_ac = ac;
	}
}
void destroy_stats(args_t *args)
{
	int id,nstats = args->files.nreaders==1 ? 1 : 3;
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		if (stats->ac_ts) free(stats->ac_ts);
		if (stats->ac_tv) free(stats->ac_tv);
		if (stats->ac_snps) free(stats->ac_snps);
		if (stats->ac_indels) free(stats->ac_indels);
		free(stats->insertions);
		free(stats->deletions);
	}
	if (args->tmp_ac) free(args->tmp_ac);
}

void init_ac(args_t *args, reader_t *reader)
{
	bcf1_t *line = reader->line;
	if ( args->ntmp_ac < line->n_allele )
	{
		args->tmp_ac = (int*)realloc(args->tmp_ac, line->n_allele*sizeof(int));
		args->ntmp_ac = line->n_allele;
	}
	calc_ac(reader->header, line, args->tmp_ac, BCF_UN_FMT);
	int i, m_ac=0;
	for (i=0; i<line->n_allele; i++)
		if ( m_ac < args->tmp_ac[i] ) m_ac = args->tmp_ac[i]; 
	if ( m_ac >= args->stats[0].m_ac )
		realloc_stats(args, m_ac+1);
}

void do_indel_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_indels++;

	bcf1_t *line = reader->line;
	int i;
	for (i=1; i<line->n_allele; i++)
	{
		if ( line->d.var[i].type!=VCF_INDEL ) continue;
		stats->ac_indels[ args->tmp_ac[i] ]++;
		int len  = line->d.var[i].n;
		int *ptr = stats->insertions;
		if ( len<0 ) 
		{
			len *= -1;
			ptr = stats->deletions;
		}
		if ( --len > stats->m_indel ) len = stats->m_indel-1;
		ptr[len]++;
	}
}
void do_snp_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_snps++;

	bcf1_t *line = reader->line;
	int ref = acgt2int(*line->d.allele[0]);
	if ( ref<0 ) return;

	int i;
	for (i=1; i<line->n_allele; i++)
	{
		if ( !(line->d.var[i].type&VCF_SNP) ) continue;
		int alt = acgt2int(*line->d.allele[i]);
		if ( alt<0 ) continue;
		stats->subst[ref<<2|alt]++;
		int ac = args->tmp_ac[i];
		stats->ac_snps[ac]++;
		if ( abs(ref-alt)==2 ) 
			stats->ac_ts[ac]++;
		else 
			stats->ac_tv[ac]++;
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
		set_variant_types(line);
		init_ac(args, reader);

		if ( line->d.var_type&VCF_SNP ) 
			do_snp_stats(args, &args->stats[ret-1], reader);
		if ( line->d.var_type&VCF_INDEL )
			do_indel_stats(args, &args->stats[ret-1], reader);

		if ( line->n_allele>2 ) args->stats[ret-1].n_mals++;
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
		printf("SN\t%d\tnumber of multiallelic sites:\t%d\n", id, stats->n_mals);

		int ts=0,tv=0;
		for (i=0; i<stats->m_ac; i++) { ts += stats->ac_ts[i]; tv += stats->ac_tv[i];  }
		printf("SN\t%d\tts/tv:\t%.2f\n", id, tv?(float)ts/tv:0);
	}
	printf("# Stats by non-reference allele count:\n# AC\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=0; i<stats->m_ac; i++)
			printf("AC\t%d\t%d\t%d\t%d\t%d\t%d\n", id,i,stats->ac_snps[i],stats->ac_ts[i],stats->ac_tv[i],stats->ac_indels[i]);
	}
	printf("# InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]count\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=stats->m_indel-1; i>=0; i--)
			if ( stats->deletions[i] ) printf("IDD\t%d\t%d\t%d\n", id,-i-1,stats->deletions[i]);
		for (i=0; i<stats->m_indel; i++)
			if ( stats->insertions[i] ) printf("IDD\t%d\t%d\t%d\n", id,i+1,stats->insertions[i]);
	}
	printf("# Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
	for (id=0; id<nstats; id++)
	{
		int t;
		for (t=0; t<15; t++)
		{
			if ( t>>2 == (t&3) ) continue;
			printf("ST\t%d\t%c>%c\t%d\n", id, int2acgt(t>>2),int2acgt(t&3),args->stats[id].subst[t]);
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

