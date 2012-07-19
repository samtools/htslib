#include <stdio.h>
#include <unistd.h>
#include "vcf.h"
#include "synced_bcf_reader.h"

typedef struct
{
	int n;
}
stats_t;

typedef struct
{
	stats_t stats[3];
	int nstats;
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


void check_vcf(args_t *args)
{
	int ret,i;
	readers_t *files = &args->files;
	while ( (ret=next_line(files)) )
	{
		for (i=0; i<files->nreaders; i++)
		{
			if ( !(ret&1<<i) ) continue;
			bcf1_t *line = files->readers[i].line;
			bcf_unpack(line, BCF_UN_STR);
			printf("ret=%d .. %s:%d %s %s\n", ret, files->seqs[files->iseq],line->pos+1, line->d.allele[0],line->d.allele[1]);
			break;
		}
	}

	// bcf1_t *line = bcf_init1();

	// printf("# This file was produced by vcfcheck. The command was:\n");
	// printf("# \thtscmd %s ", args->argv[0]);
	// int i;
	// for (i=1; i<args->argc; i++)
	// 	printf(" %s",args->argv[i]);
	// printf("\n");

	// return;
	// while (vcf_read1(args->file, args->header, line) >= 0)
	// {
	// 	printf("%d\n", line->pos+1);
	// }

	// bcf_destroy1(line);
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
	if (argc == optind) {
		fprintf(stderr, "\nAbout:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfcheck.\n");
		fprintf(stderr, "             When two files are given, gives separate stats for intersection and the complements.\n");
		fprintf(stderr, "Usage:   vcfcheck [options] <A.vcf.gz> [<B.vcf.gz>]\n\n");
		// fprintf(stderr, "Options: -b           output in BCF\n");
		// fprintf(stderr, "         -S           input is VCF\n");
		// fprintf(stderr, "         -o FILE      output file name [stdout]\n");
		// fprintf(stderr, "         -l INT       compression level [%d]\n", clevel);
		// fprintf(stderr, "         -t FILE      list of reference names and lengths [null]\n");
		// fprintf(stderr, "         -s FILE/STR  list of samples (STR if started with ':'; FILE otherwise) [null]\n");
		// fprintf(stderr, "         -G           drop individual genotype information\n");
		fprintf(stderr, "\n");
		return 1;
	}

	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;
	while (optind<argc)
	{
		if ( !add_reader(argv[optind], &args->files) ) error("Could not load the index: %s\n", argv[optind]);
		optind++;
	}
	int i;
	for (i=0; i<args->files.nseqs; i++)
		fprintf(stderr,"chr%s\n", args->files.seqs[i]);

	check_vcf(args);

	destroy_readers(&args->files);
	free(args);
	return 0;
}

