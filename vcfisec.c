#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

#define OP_PLUS 1
#define OP_MINUS 2
#define OP_EQUAL 3

typedef struct
{
    int isec_op, isec_n, *write;
	readers_t *files;
    FILE *fh_log, *fh_sites;
    htsFile **fh_out;
	char **argv, *prefix, **fnames, *write_files, *subset_fname;
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

void isec_vcf(args_t *args)
{
    int ret,i;
    readers_t *files = args->files;
    kstring_t str = {0,0,0};
    htsFile *out_fh = NULL;
    
    // With --subset and only one file present, output the VCF
    int output_vcf = args->subset_fname && files->nreaders==1 ? 1 : 0;
    if ( output_vcf ) 
    {
        out_fh = hts_open("-","w",0);
        vcf_hdr_write(out_fh, files->readers[0].header);
    }

    while ( (ret=bcf_sr_next_line(files)) )
    {
        reader_t *reader = NULL;
        bcf1_t *line = NULL;
        int n = 0;
        for (i=0; i<files->nreaders; i++)
        {
            if ( !(ret&1<<i) ) continue;
            if ( !line ) 
            {
                line = files->readers[i].buffer[0];
                reader = &files->readers[i];
            }
            n++;
        }

        if ( output_vcf )
        {
            vcf_write1(out_fh, files->readers[0].header, line);
            continue;
        }

        if ( args->isec_op==OP_EQUAL )
            if ( n != args->isec_n ) continue;
        else {
            if ( args->isec_op==OP_PLUS )
                if ( n < args->isec_n ) continue;
        }
        else {
            if ( args->isec_op==OP_MINUS )
                if ( n > args->isec_n ) continue;
        }

        str.l = 0;
        kputs(reader->header->id[BCF_DT_CTG][line->rid].key, &str); kputc('\t', &str);
        kputw(line->pos+1, &str); kputc('\t', &str);
        if (line->n_allele > 0) kputs(line->d.allele[0], &str);
        else kputc('.', &str);
        kputc('\t', &str);
        if (line->n_allele > 1) kputs(line->d.allele[1], &str);
        else kputc('.', &str);
        for (i=2; i<line->n_allele; i++)
        {
            kputc(',', &str);
            kputs(line->d.allele[i], &str);
        }
        kputc('\n', &str);
        fwrite(str.s,sizeof(char),str.l,args->fh_sites);

        if ( args->prefix )
        {
            for (i=0; i<files->nreaders; i++)
            {
                if ( !(ret&1<<i) ) continue;
                if ( args->write && !args->write[i] ) continue;
                vcf_write1(args->fh_out[i], files->readers[i].header, files->readers[i].buffer[0]);
            }
        }
    }
    if ( str.s ) free(str.s);
    if ( output_vcf ) hts_close(out_fh);
}

void mkdir_p(const char *path)
{
    char *tmp = strdup(path), *p = tmp+1;
    while (*p)
    {
        while (*p && *p!='/') p++;
        if ( *p )
        {
            *p = 0;
            mkdir(tmp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            *p = '/';
            p++;
        }
        else
            mkdir(tmp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    free(tmp);
}

static void init_data(args_t *args)
{
    if ( args->prefix )
    {
        // Which files to write
        char *p = args->write_files;
        while (p && *p)
        {
            int i;
            if ( !args->write ) args->write = (int*) calloc(args->files->nreaders,sizeof(int));
            if ( sscanf(p,"%d",&i)!=1 ) error("Could not parse --write %s\n", args->write_files);
            if ( i<0 || i>args->files->nreaders ) error("The index is out of range: %d (%s)\n", i, args->write_files);
            args->write[i-1] = 1;
            while (*p && *p!=',') p++;
            if ( *p==',' ) p++;
        }
        
        // Init directories and open files for writing
        mkdir_p(args->prefix);
        int i,ntmp = strlen(args->prefix)+20;
        char *tmp = (char*)malloc(ntmp);
        snprintf(tmp,ntmp,"%s/README.txt",args->prefix);
        args->fh_log = fopen(tmp,"w");
        if ( !args->fh_log ) error("%s: %s\n", tmp,strerror(errno));
        fprintf(args->fh_log,"This file was produced by vcfisec.\n");
        fprintf(args->fh_log,"The command line was:\thtscmd %s ", args->argv[0]);
        for (i=1; i<args->argc; i++)
            fprintf(args->fh_log," %s",args->argv[i]);

        fprintf(args->fh_log,"\n\nUsing the following file names:\n");
        args->fh_out = (htsFile**) malloc(sizeof(htsFile*)*args->files->nreaders);
        args->fnames = (char**) malloc(sizeof(char*)*args->files->nreaders);
        for (i=0; i<args->files->nreaders; i++)
        {   
            if ( args->write && !args->write[i] ) continue;
            snprintf(tmp,ntmp,"%s/%04d.bcf", args->prefix,i);
            fprintf(args->fh_log,"%s\tfor stripped\t%s\n", tmp,args->files->readers[i].fname);
            args->fh_out[i] = hts_open(tmp, "wb", 0);
            if ( !args->fh_out[i] ) error("Could not open %s\n", tmp);
            vcf_hdr_write(args->fh_out[i], args->files->readers[i].header);
            args->fnames[i] = strdup(tmp);
        }
        snprintf(tmp,ntmp,"%s/sites.txt",args->prefix);
        args->fh_sites = fopen(tmp,"w");
        if ( !args->fh_sites ) error("%s: %s\n", tmp,strerror(errno));
        free(tmp);
    }
    else
        args->fh_sites = stdout;
}

static void destroy_data(args_t *args)
{
    if ( args->prefix )
    {
        int i;
        fclose(args->fh_log);
        for (i=0; i<args->files->nreaders; i++)
        {
            if ( args->write && !args->write[i] ) continue;
            hts_close(args->fh_out[i]);
            if ( bcf_index_build(args->fnames[i],14) ) error("Could not index %s\n", args->fnames[i]);
            free(args->fnames[i]);
        }
        free(args->fnames);
        fclose(args->fh_sites);
        if ( args->write ) free(args->write);
    }
}

static void usage(void)
{
	fprintf(stderr, "About:   Create intersections, unions and complements of VCF files\n");
	fprintf(stderr, "Usage:   vcfisec [options] <A.vcf.gz> <B.vcf.gz> ...\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -c, --collapse <string>           treat as identical sites with differing alleles for <snps|indels|both|any>\n");
	fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
	fprintf(stderr, "    -n, --nfiles [+-=]<int>           output positions present in this many (=), this many or more (+), or this many or fewer (-) files\n");
	fprintf(stderr, "    -p, --prefix <dir>                if given, subset each of the input files accordingly, see also -w\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    perform intersection in the given region only\n");
	fprintf(stderr, "    -s, --subset <file>               subset to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
	fprintf(stderr, "    -w, --write <list>                list of files to write with -p given as 1-based indexes. By default, all files are written\n");
	fprintf(stderr, "Examples:\n");
	fprintf(stderr, "   # Extract and write records from A shared by both A and B using exact allele match\n");
	fprintf(stderr, "   vcf isec A.vcf.gz B.vcf.gz -p dir -n =2 -w 1\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "   # Extract records private to A or B comparing by position only\n");
	fprintf(stderr, "   vcf isec A.vcf.gz B.vcf.gz -p dir -n -1 -c any\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfisec(int argc, char *argv[])
{
	int c;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"collapse",1,0,'c'},
		{"apply-filters",0,0,'f'},
		{"nfiles",1,0,'n'},
		{"prefix",1,0,'p'},
		{"write",1,0,'w'},
		{"subset",1,0,'s'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hc:fr:p:n:w:s:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'c':
				if ( !strcmp(optarg,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) args->files->collapse |= COLLAPSE_ANY;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
				break;
			case 'f': args->files->apply_filters = 1; break;
			case 'r': args->files->region = optarg; break;
			case 's': args->subset_fname = optarg; break;
			case 'p': args->prefix = optarg; break;
			case 'w': args->write_files = optarg; break;
			case 'n': 
                {
                    char *p = optarg;
                    if ( *p=='-' ) { args->isec_op = OP_MINUS; p++; }
                    else if ( *p=='+' ) { args->isec_op = OP_PLUS; p++; }
                    else if ( *p=='=' ) { args->isec_op = OP_EQUAL; p++; }
                    else if ( isdigit(*p) ) args->isec_op = OP_EQUAL;
                    else error("Could not parse --nfiles %s\n", optarg);
                    if ( sscanf(p,"%d",&args->isec_n)!=1 ) error("Could not parse --nfiles %s\n", optarg);
                }
                break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if ( argc-optind<1 ) usage();   // no file given
    if ( args->subset_fname )
    {
        if ( !bcf_sr_set_targets(args->files, args->subset_fname) )
            error("Failed to read the targets: %s\n", args->subset_fname);
    }
    else
    {
        if ( argc-optind<2 ) error("Expected multiple files or the --subset option\n");
        if ( !args->isec_op ) error("Missing --nfiles or --subset option\n");
    }
	while (optind<argc)
	{
		if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open: %s\n", argv[optind]);
		optind++;
	}
    init_data(args);
	isec_vcf(args);
    destroy_data(args);
	bcf_sr_destroy(args->files);
	free(args);
	return 0;
}

