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
#define OP_VENN 4
#define OP_COMPLEMENT 5

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

/**
 *  mkdir_p() - create new directory for a file $fname
 *  @fname:   the file name to create the directory for, the part after last "/" is ignored
 */
void mkdir_p(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *path = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(path, n, fmt, ap);
    va_end(ap);

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
    }
    free(tmp);
    free(path);
}

/**
 *  open_file() - open new file creating the file name using vsnprintf
 *  @fname:  if not NULL, on output will point to newly allocated fname string
 *  @mode:   if NULL, only the file name string will be created   
 *  @fmt:    vsnprintf format and args
 *
 *  Returns open file descriptor or NULL if mode is NULL.
 */
FILE *open_file(char **fname, const char *mode, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *str = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(str, n, fmt, ap);
    va_end(ap);

    mkdir_p(str);
    if ( !mode )
    {
        if ( !fname ) error("Uh: expected fname or mode\n");
        *fname = str;
        return NULL;
    }

    FILE *fp = fopen(str,mode);
    if ( fname ) *fname = str;
    else free(str);
    return fp;
}


void isec_vcf(args_t *args)
{
    int ret,i;
    readers_t *files = args->files;
    kstring_t str = {0,0,0};
    htsFile *out_fh = NULL;

    // When only one VCF is output, print VCF to stdout
    int stdout_vcf = args->subset_fname && files->nreaders==1 ? 1 : 0;
    if ( args->isec_op==OP_COMPLEMENT ) 
    {
        stdout_vcf = 0;
        args->isec_n = 1;
    }
    if ( stdout_vcf ) 
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

        switch (args->isec_op) 
        {
            case OP_COMPLEMENT: if ( ret != args->isec_n ) continue; break;
            case OP_EQUAL: if ( n != args->isec_n ) continue; break;
            case OP_PLUS: if ( n < args->isec_n ) continue; break;
            case OP_MINUS: if ( n > args->isec_n ) continue;
        }

        if ( stdout_vcf )
        {
            vcf_write1(out_fh, files->readers[0].header, line);
            continue;
        }
        else if ( args->fh_sites )
        {
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
        }

        if ( args->prefix )
        {
            if ( args->isec_op==OP_VENN )
                vcf_write1(args->fh_out[ret-1], reader->header, line);
            else
            {
                for (i=0; i<files->nreaders; i++)
                {
                    if ( !(ret&1<<i) ) continue;
                    if ( args->write && !args->write[i] ) continue;
                    vcf_write1(args->fh_out[i], files->readers[i].header, files->readers[i].buffer[0]);
                }
            }
        }
    }
    if ( str.s ) free(str.s);
    if ( out_fh ) hts_close(out_fh);
}

static void destroy_data(args_t *args);
static void init_data(args_t *args)
{
    if ( args->prefix )
    {
        // Which files to write: parse the string passed with -w
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
        
        // Init output directory and create the readme file
        args->fh_log = open_file(NULL,"w","%s/README.txt", args->prefix);
        if ( !args->fh_log ) error("%s/README.txt: %s\n", args->prefix, strerror(errno));

        fprintf(args->fh_log,"This file was produced by vcfisec.\n");
        fprintf(args->fh_log,"The command line was:\thtscmd %s ", args->argv[0]);
        int i;
        for (i=1; i<args->argc; i++) fprintf(args->fh_log," %s",args->argv[i]);
        fprintf(args->fh_log,"\n\nUsing the following file names:\n");

        // Open output files and write the legend
        if ( args->isec_op==OP_VENN )
        {
            args->fh_out = (htsFile**) malloc(sizeof(htsFile*)*3);
            args->fnames = (char**) malloc(sizeof(char*)*3);

            #define OPEN_FILE(i,j) { \
                open_file(&args->fnames[i], NULL, "%s/%04d.bcf", args->prefix, i); \
                args->fh_out[i] = hts_open(args->fnames[i], "wb", 0);  \
                if ( !args->fh_out[i] ) error("Could not open %s\n", args->fnames[i]); \
                vcf_hdr_write(args->fh_out[i], args->files->readers[j].header); \
            }
            OPEN_FILE(0,0);
            fprintf(args->fh_log,"%s\tfor records private to\t%s\n", args->fnames[0], args->files->readers[0].fname);
            OPEN_FILE(1,1);
            fprintf(args->fh_log,"%s\tfor records private to\t%s\n", args->fnames[1], args->files->readers[1].fname);
            OPEN_FILE(2,0);
            fprintf(args->fh_log,"%s\tfor records shared by both\t%s %s\n", args->fnames[2], args->files->readers[0].fname, args->files->readers[1].fname);
        }
        else
        {
            // Init one output file for each reader
            args->fh_out = (htsFile**) calloc(args->files->nreaders, sizeof(htsFile*));
            args->fnames = (char**) calloc(args->files->nreaders, sizeof(char*));

            for (i=0; i<args->files->nreaders; i++)
            {   
                if ( args->write && !args->write[i] ) continue;
                if ( args->isec_op==OP_COMPLEMENT && i>0 ) break;
                OPEN_FILE(i,i);
                fprintf(args->fh_log,"%s\tfor stripped\t%s\n", args->fnames[i], args->files->readers[i].fname);
            }
            #undef OPEN_FILE

            args->fh_sites = open_file(NULL, "w", "%s/sites.txt", args->prefix);
            if ( !args->fh_sites ) error("%s/sites.txt: %s\n", args->prefix, strerror(errno));
        }
    }
    else
        args->fh_sites = stdout;
}

static void destroy_data(args_t *args)
{
    if ( args->prefix )
    {
        fclose(args->fh_log);
        int i, n = args->isec_op==OP_VENN ? 3 : args->files->nreaders;
        for (i=0; i<n; i++)
        {
            if ( !args->fnames[i] ) continue;
            hts_close(args->fh_out[i]);
            if ( bcf_index_build(args->fnames[i],14) ) error("Could not index %s\n", args->fnames[i]);
            free(args->fnames[i]);
        }
        free(args->fnames);
        if ( args->fh_sites ) fclose(args->fh_sites);
        if ( args->write ) free(args->write);
    }
}

static void usage(void)
{
	fprintf(stderr, "About:   Create intersections, unions and complements of VCF files\n");
	fprintf(stderr, "Usage:   vcfisec [options] <A.vcf.gz> <B.vcf.gz> ...\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -c, --collapse <string>           treat as identical sites with differing alleles for <snps|indels|both|any>\n");
	fprintf(stderr, "    -C, --complement                  output positions present only in the first file but missing in the others\n");
	fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
	fprintf(stderr, "    -n, --nfiles [+-=]<int>           output positions present in this many (=), this many or more (+), or this many or fewer (-) files\n");
	fprintf(stderr, "    -p, --prefix <dir>                if given, subset each of the input files accordingly, see also -w\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    perform intersection in the given region only\n");
	fprintf(stderr, "    -s, --subset <file>               subset to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
	fprintf(stderr, "    -w, --write <list>                list of files to write with -p given as 1-based indexes. By default, all files are written\n");
	fprintf(stderr, "Examples:\n");
	fprintf(stderr, "   # Create intersection and complements of two sets saving the output in dir/*\n");
	fprintf(stderr, "   vcf isec A.vcf.gz B.vcf.gz -p dir\n");
	fprintf(stderr, "\n");
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
		{"complement",1,0,'C'},
		{"apply-filters",0,0,'f'},
		{"nfiles",1,0,'n'},
		{"prefix",1,0,'p'},
		{"write",1,0,'w'},
		{"subset",1,0,'s'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hc:fr:p:n:w:s:C",loptions,NULL)) >= 0) {
		switch (c) {
			case 'c':
				if ( !strcmp(optarg,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) args->files->collapse |= COLLAPSE_ANY;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
				break;
			case 'f': args->files->apply_filters = 1; break;
			case 'C': args->isec_op = OP_COMPLEMENT; break;
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
    if ( argc-optind==2 && !args->isec_op ) 
    {
        args->isec_op = OP_VENN;
        if ( !args->prefix ) error("Expected the -p option\n");
    }
    if ( !args->subset_fname )
    {
        if ( argc-optind<2  ) error("Expected multiple files or the --subset option\n");
        if ( !args->isec_op ) error("Expected two file names or one of the options --complement, --nfiles or --subset\n");
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

