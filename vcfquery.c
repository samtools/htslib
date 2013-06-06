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

#define T_CHROM   1
#define T_POS     2
#define T_ID      3
#define T_REF     4
#define T_ALT     5
#define T_QUAL    6
#define T_FILTER  7
#define T_INFO    8
#define T_FORMAT  9
#define T_SAMPLE  10
#define T_SEP     11
#define T_IS_TS   12
#define T_TYPE    13
#define T_MASK    14
#define T_GT      15
#define T_TGT     16
#define T_LINE    17

struct _args_t;
typedef struct _args_t args_t;

typedef struct _fmt_t
{
    int type, id, is_gt_field, ready;
    char *key;
    bcf_fmt_t *fmt;
    void (*handler)(args_t *, bcf1_t *, struct _fmt_t *, int, kstring_t *);
}
fmt_t;

struct _args_t
{
    fmt_t *fmt;
    int nfmt, mfmt;
    int nsamples, *samples;
	bcf_srs_t *files;
    bcf_hdr_t *header;
	char **argv, *format, *sample_names, *subset_fname, *vcf_list;
	int argc, list_columns, print_header;
};

char **read_list(char *fname, int *n);
void destroy_list(char **list, int n);

static void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(-1);
}

static void process_chrom(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(args->header->id[BCF_DT_CTG][line->rid].key, str); }
static void process_pos(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+1, str); }
static void process_id(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.id, str); }
static void process_ref(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.allele[0], str); }
static void process_alt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{ 
    int i;
    if ( line->n_allele==1 )
    {
        kputc('.', str);
        return;
    }
    for (i=1; i<line->n_allele; i++)
    {
        if ( i>1 ) kputc(',', str);
        kputs(line->d.allele[i], str); 
    }
}
static void process_qual(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{  
    if ( memcmp(&line->qual, &bcf_float_missing, 4) == 0) kputc('.', str);
    else ksprintf(str, "%g", line->qual);
}
static void process_filter(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{ 
    int i;
    if ( line->d.n_flt ) 
    {
        for (i=0; i<line->d.n_flt; i++) 
        {
            if (i) kputc(';', str);
            kputs(args->header->id[BCF_DT_ID][line->d.flt[i]].key, str);
        }
    } 
    else kputc('.', str);
}
static void process_info(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{
    int i;
    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == fmt->id ) break;

    // output "." if the tag is not present
    if ( i==line->n_info )
    {
        kputc('.', str);
        return;
    }

    bcf_info_t *info = &line->d.info[i];

    // if this is a flag, output 1
    if ( info->len <=0 )
    {
        kputc('1', str);
        return;
    }
        
    if ( info->len == 1 )
    {
        if ( info->type == BCF_BT_FLOAT ) ksprintf(str, "%g", info->v1.f);
        else if ( info->type != BCF_BT_CHAR ) kputw(info->v1.i, str);
        else kputc(info->v1.i, str);
    } 
    else 
        bcf_fmt_array(str, info->len, info->type, info->vptr);
}
static void init_format(args_t *args, bcf1_t *line, fmt_t *fmt)
{
    fmt->id = bcf_id2int(args->header, BCF_DT_ID, fmt->key);
    if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: FORMAT/%s\n", fmt->key);
    fmt->fmt = NULL;
    int i;
    for (i=0; i<(int)line->n_fmt; i++)  
        if ( line->d.fmt[i].id==fmt->id ) { fmt->fmt = &line->d.fmt[i]; break; }
    fmt->ready = 1;
}
static void process_format(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_fmt_array(str, fmt->fmt->n, fmt->fmt->type, fmt->fmt->p + isample*fmt->fmt->size);
}
static void process_gt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_format_gt(fmt->fmt, isample, str);
}
static void process_tgt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }

    assert( fmt->fmt->type==BCF_BT_INT8 );

    int l;
    int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size); // FIXME: does not work with n_alt >= 64
    for (l = 0; l < fmt->fmt->n && x[l] != bcf_int8_vector_end; ++l)
    {
        if (l) kputc("/|"[x[l]&1], str);
        if (x[l]>>1) 
        {
            int ial = (x[l]>>1) - 1;
            kputs(line->d.allele[ial], str);
        }
        else 
            kputc('.', str);
    }
    if (l == 0) kputc('.', str);
}
static void process_sample(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{  
    kputs(args->header->samples[isample], str);
}
static void process_sep(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { if (fmt->key) kputs(fmt->key, str); }
static void process_is_ts(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) 
{ 
    bcf_set_variant_types(line);
    int is_ts = 0;
    if ( line->d.var_type & (VCF_SNP|VCF_MNP) ) 
        is_ts = abs(bcf_acgt2int(*line->d.allele[0])-bcf_acgt2int(*line->d.allele[1])) == 2 ? 1 : 0;
    kputc(is_ts ? '1' : '0', str); 
}
static void process_type(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    bcf_set_variant_types(line);
    int i = 0;
    if ( line->d.var_type == VCF_REF ) { kputs("REF", str); i++; }
    if ( line->d.var_type & VCF_SNP ) { if (i) kputc(',',str); kputs("SNP", str); i++; }
    if ( line->d.var_type & VCF_MNP ) { if (i) kputc(',',str); kputs("MNP", str); i++; }
    if ( line->d.var_type & VCF_INDEL ) { if (i) kputc(',',str); kputs("INDEL", str); i++; }
    if ( line->d.var_type & VCF_OTHER ) { if (i) kputc(',',str); kputs("OTHER", str); i++; }
}
static void process_line(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    vcf_format1(args->header, line, str);
}

static void register_tag(args_t *args, int type, char *key, int is_gtf)
{
    args->nfmt++;
    if ( args->nfmt > args->mfmt )
    {
        args->mfmt += 10;
        args->fmt   = (fmt_t*) realloc(args->fmt, args->mfmt*sizeof(fmt_t));
            
    }
    fmt_t *fmt = &args->fmt[ args->nfmt-1 ];
    fmt->type  = type;
    fmt->key   = key ? strdup(key) : NULL;
    fmt->is_gt_field = is_gtf;
    switch (fmt->type)
    {
        case T_CHROM: fmt->handler = &process_chrom; break;
        case T_POS: fmt->handler = &process_pos; break;
        case T_ID: fmt->handler = &process_id; break;
        case T_REF: fmt->handler = &process_ref; break;
        case T_ALT: fmt->handler = &process_alt; break;
        case T_QUAL: fmt->handler = &process_qual; break;
        case T_FILTER: fmt->handler = &process_filter; args->files->max_unpack |= BCF_UN_FLT; break;
        case T_INFO: fmt->handler = &process_info; args->files->max_unpack |= BCF_UN_INFO; break;
        case T_FORMAT: fmt->handler = &process_format; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_SAMPLE: fmt->handler = &process_sample; break;
        case T_SEP: fmt->handler = &process_sep; break;
        case T_IS_TS: fmt->handler = &process_is_ts; break;
        case T_TYPE: fmt->handler = &process_type; break;
        case T_MASK: fmt->handler = NULL; break;
        case T_GT: fmt->handler = &process_gt; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_TGT: fmt->handler = &process_tgt; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_LINE: fmt->handler = &process_line; break;
        default: error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key )
    {
        if ( fmt->type==T_INFO )
        {
            fmt->id = bcf_id2int(args->header, BCF_DT_ID, key);
            if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: INFO/%s\n", key);
        }
    }
}

static char *parse_tag(args_t *args, char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_') ) q++;
    kstring_t str = {0,0,0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", args->format);
    kputsn(p, q-p, &str);
    if ( is_gtf )
    {
        if ( !strcmp(str.s, "SAMPLE") ) register_tag(args, T_SAMPLE, "SAMPLE", is_gtf);
        else if ( !strcmp(str.s, "GT") ) register_tag(args, T_GT, "GT", is_gtf);
        else if ( !strcmp(str.s, "TGT") ) register_tag(args, T_TGT, "GT", is_gtf);
        else register_tag(args, T_FORMAT, str.s, is_gtf);
    }
    else
    {
        if ( !strcmp(str.s, "CHROM") ) register_tag(args, T_CHROM, str.s, is_gtf);
        else if ( !strcmp(str.s, "POS") ) register_tag(args, T_POS, str.s, is_gtf);
        else if ( !strcmp(str.s, "ID") ) register_tag(args, T_ID, str.s, is_gtf);
        else if ( !strcmp(str.s, "REF") ) register_tag(args, T_REF, str.s, is_gtf);
        else if ( !strcmp(str.s, "ALT") ) register_tag(args, T_ALT, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(args, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "FILTER") ) register_tag(args, T_FILTER, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(args, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "IS_TS") ) register_tag(args, T_IS_TS, str.s, is_gtf);
        else if ( !strcmp(str.s, "TYPE") ) register_tag(args, T_TYPE, str.s, is_gtf);
        else if ( !strcmp(str.s, "MASK") ) register_tag(args, T_MASK, str.s, is_gtf);
        else if ( !strcmp(str.s, "LINE") ) register_tag(args, T_LINE, str.s, is_gtf);
        else if ( !strcmp(str.s, "INFO") ) 
        {
            if ( *q!='/' ) error("Could not parse format string: %s\n", args->format);
            p = ++q;
            str.l = 0;
            while ( *q && isalnum(*q) ) q++;
            if ( q-p==0 ) error("Could not parse format string: %s\n", args->format);
            kputsn(p, q-p, &str);
            register_tag(args, T_INFO, str.s, is_gtf);
        }
        else 
            register_tag(args, T_INFO, str.s, is_gtf);
    }
    free(str.s);
    return q;
}

static char *parse_sep(args_t *args, char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = {0,0,0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' )
    {
        if ( *q=='\\' ) 
        {
            q++;
            if ( *q=='n' ) kputc('\n', &str);
            else if ( *q=='t' ) kputc('\t', &str);
            else kputc(*q, &str);
        }
        else kputc(*q, &str);
        q++;
    }
    if ( !str.l ) error("Could not parse format string: %s\n", args->format);
    register_tag(args, T_SEP, str.s, is_gtf);
    free(str.s);
    return q;
}

static void init_data(args_t *args)
{
    args->header = args->files->readers[0].header;
    int is_gtf = 0;
    char *p = args->format;
    while ( *p )
    {
        //fprintf(stderr,"<%s>\n", p);
        switch (*p) 
        {
            case '[': is_gtf = 1; p++; break;
            case ']': is_gtf = 0; register_tag(args, T_SEP, NULL, 0); p++; break;
            case '%': p = parse_tag(args, p, is_gtf); break;
            default:  p = parse_sep(args, p, is_gtf); break;
        }
    }
    if ( args->sample_names )
    {
        kstring_t str = {0,0,0};
        struct stat sbuf;
        if ( stat(args->sample_names, &sbuf) == 0  )
            error("TODO: reading sample names from a file, please email me\n");
        else
        {
            char *p = args->sample_names, *q = p;
            while (1)
            {
                if ( !*q || *q==',' )
                {
                    str.l = 0;
                    kputsn(p, q-p, &str);
                    int i;
                    if ( (i=bcf_id2int(args->header, BCF_DT_SAMPLE, str.s))==-1 ) error("No such sample in the VCF: [%s]\n", str.s);
                    args->samples = (int*) realloc(args->samples, sizeof(int)*(++args->nsamples));
                    args->samples[ args->nsamples-1 ] = i;
                    p = q+1;
                }
                if ( !*q ) break;
                q++;
            }
        }
        if ( str.m ) free(str.s);
    }
    else
    {
        args->nsamples = args->header->n[BCF_DT_SAMPLE];
        if ( args->nsamples )
        {
            args->samples = (int*)malloc(sizeof(int)*args->nsamples);
            int i;
            for (i=0; i<args->nsamples; i++)
                args->samples[i] = i;
        }
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfmt; i++)
        if ( args->fmt[i].key ) free(args->fmt[i].key);
    if ( args->mfmt ) free(args->fmt);
    if ( args->samples ) free(args->samples);
    args->nfmt = args->mfmt = 0;
    args->fmt = NULL;
}


static void print_header(args_t *args, kstring_t *str)
{
    int i, icol = 0;

    // Supress the header output if LINE is present
    for (i=0; i<args->nfmt; i++)
        if ( args->fmt[i].type == T_LINE ) break;
    if ( i!=args->nfmt )
        return;
    
    kputs("# ", str);
    for (i=0; i<args->nfmt; i++)
    {
        // Genotype fields
        if ( args->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( args->fmt[j].is_gt_field ) j++;
            for (js=0; js<args->nsamples; js++)
            {
                int ks = args->samples[js];
                for (k=i; k<j; k++)
                {
                    if ( args->fmt[k].type == T_SEP )
                    {
                        if ( args->fmt[k].key ) kputs(args->fmt[k].key, str);
                    }
                    else if ( args->fmt[k].type == T_SAMPLE )
                        ksprintf(str, "[%d]%s", ++icol, args->fmt[k].key);
                    else
                        ksprintf(str, "[%d]%s:%s", ++icol, args->header->samples[ks], args->fmt[k].key);
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( args->fmt[i].type == T_SEP )
        {
            if ( args->fmt[i].key ) kputs(args->fmt[i].key, str);
            continue;
        }
        ksprintf(str, "[%d]%s", ++icol, args->fmt[i].key);
    }
    fwrite(str->s, str->l, 1, stdout);
    str->l = 0;
}

static void query_vcf(args_t *args)
{
    kstring_t str = {0,0,0};

    args->files->max_unpack |= BCF_UN_STR;
    if ( args->print_header ) print_header(args, &str);

    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, args->files->max_unpack);

        int i, ir;
        str.l = 0;
        for (i=0; i<args->nfmt; i++)
        {
            // Genotype fields
            if ( args->fmt[i].is_gt_field )
            {
                int j = i, js, k;
                while ( args->fmt[j].is_gt_field ) 
                {
                    args->fmt[j].ready = 0;
                    j++;
                }
                for (js=0; js<args->nsamples; js++)
                {
                    int ks = args->samples[js];
                    for (k=i; k<j; k++)
                    {
                        if ( args->fmt[k].type == T_MASK )
                        {
                            for (ir=0; ir<args->files->nreaders; ir++) 
                                kputc(bcf_sr_has_line(args->files,ir)?'1':'0', &str);
                        }
                        else if ( args->fmt[k].handler )
                            args->fmt[k].handler(args, line, &args->fmt[k], ks, &str);
                    }
                }
                i = j-1;
                continue;
            }
            // Fixed fields
            if ( args->fmt[i].type == T_MASK )
            {
                for (ir=0; ir<args->files->nreaders; ir++) 
                    kputc(bcf_sr_has_line(args->files,ir)?'1':'0', &str);
            }
            else if ( args->fmt[i].handler )
                args->fmt[i].handler(args, line, &args->fmt[i], -1, &str);
        }
        if ( str.l )
            fwrite(str.s, str.l, 1, stdout);
    }
    if ( str.m ) free(str.s);
}

static void list_columns(args_t *args)
{
    int i;
    bcf_sr_t *reader = &args->files->readers[0];
    for (i=0; i<reader->header->n[BCF_DT_SAMPLE]; i++)
        printf("%s\n", reader->header->samples[i]);
}

static char **copy_header(bcf_hdr_t *hdr, int *src, int nsrc)
{
    char **dst = (char**) malloc(sizeof(char*)*nsrc);
    int i;
    for (i=0; i<nsrc; i++) dst[i] = strdup(hdr->samples[src[i]]);
    return dst;
}
static int compare_header(bcf_hdr_t *hdr, int *a, int na, char **b, int nb)
{
    if ( na!=nb ) return na-nb;
    int i;
    for (i=0; i<na; i++)
        if ( strcmp(hdr->samples[a[i]],b[i]) ) return 1;
    return 0;
}


static void usage(void)
{
	fprintf(stderr, "About:   Extracts fields from VCF/BCF file and prints them in user-defined format\n");
	fprintf(stderr, "Usage:   vcfquery [options] <file.vcf.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -a, --annots <list>               alias for -f '%%CHROM\\t%%POS\\t%%MASK\\t%%REF\\t%%ALT\\t%%TYPE\\t' + tab-separated <list> of tags\n");
	fprintf(stderr, "    -c, --collapse <string>           collapse lines with duplicate positions for <snps|indels|both|any>\n");
	fprintf(stderr, "    -f, --format <string>             learn by example, see below\n");
	fprintf(stderr, "    -H, --print-header                print header\n");
	fprintf(stderr, "    -l, --list-columns                list columns\n");
	fprintf(stderr, "    -p, --positions <file>            list positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    output from the given region only\n");
	fprintf(stderr, "    -s, --samples <list|file>         samples to include: comma-separated list or one name per line in a file\n");
	fprintf(stderr, "    -v, --vcf-list <file>             process multiple VCFs listed in the file\n");
	fprintf(stderr, "Expressions:\n");
    fprintf(stderr, "\t%%CHROM          The CHROM column (similarly also other columns, such as POS, ID, QUAL, etc.)\n");
    fprintf(stderr, "\t%%INFO/TAG       Any tag in the INFO column\n");
    fprintf(stderr, "\t%%TYPE           Variant type (REF, SNP, MNP, INDEL, OTHER)\n");
    fprintf(stderr, "\t%%MASK           Indicates presence of the site in other files (with multiple files)\n");
    fprintf(stderr, "\t[]              The brackets loop over all samples\n");
    fprintf(stderr, "\t%%GT             Genotype (e.g. 0/1)\n");
    fprintf(stderr, "\t%%TGT            Translated genotype (e.g. C/A)\n");
    fprintf(stderr, "\t%%LINE           Prints the whole line\n");
    fprintf(stderr, "\t%%SAMPLE         Sample name\n");
    //fprintf(stderr, "\t%*<A><B>        All format fields printed as KEY<A>VALUE<B>\n");
	fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\tvcfquery -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n' file.vcf.gz\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfquery(int argc, char *argv[])
{
	int c, collapse = 0;
    char *region = NULL;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"list-columns",0,0,'l'},
		{"format",1,0,'f'},
		{"region",1,0,'r'},
		{"annots",1,0,'a'},
		{"samples",1,0,'s'},
		{"print-header",0,0,'H'},
		{"positions",1,0,'p'},
		{"collapse",1,0,'c'},
		{"vcf-list",1,0,'v'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hlr:f:a:s:Hp:c:v:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'f': args->format = strdup(optarg); break;
			case 'H': args->print_header = 1; break;
            case 'v': args->vcf_list = optarg; break;
            case 'c':
				if ( !strcmp(optarg,"snps") ) collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) collapse |= COLLAPSE_ANY;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
				break;
			case 'a': 
                {
                    kstring_t str = {0,0,0};
                    kputs("%CHROM\t%POS\t%MASK\t%REF\t%ALT\t%TYPE\t%", &str);
                    char *p = optarg;
                    while ( *p )
                    {
                        if ( *p==',' ) 
                            kputs("\t%", &str);
                        else
                            kputc(*p, &str);
                        p++;
                    }
                    kputc('\n', &str);
                    args->format = str.s;
                    break;
                }
			case 'r': region = optarg; break;
			case 'l': args->list_columns = 1; break;
			case 's': args->sample_names = optarg; break;
			case 'p': args->subset_fname = optarg; break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
    
    if ( args->list_columns )
    {
        if ( optind==argc ) error("Missing the VCF file name\n");
        args->files = bcf_sr_init();
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
        list_columns(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    if ( !args->format ) usage();
    if ( !args->vcf_list )
    {
        args->files = bcf_sr_init();
        args->files->region = region;
        args->files->collapse = collapse;
        if ( optind+1 < argc ) args->files->require_index = 1;
        if ( optind==argc ) usage();
        if ( args->subset_fname )
        {
            args->files->require_index = 1;
            if ( !bcf_sr_set_targets(args->files, args->subset_fname) )
                error("Failed to read the targets: %s\n", args->subset_fname);
        }
        while (optind<argc)
        {
            if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
            optind++;
        }
        init_data(args);
        query_vcf(args);
        free(args->format);
        destroy_data(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    // multiple VCFs
    int i, k, nfiles, prev_nsamples = 0;
    char **fnames, **prev_samples = NULL;
    fnames = read_list(args->vcf_list, &nfiles);
    if ( !nfiles ) error("No files in %s?\n", args->vcf_list);
    for (i=0; i<nfiles; i++)
    {
        args->files = bcf_sr_init();
        args->files->region = region;
        args->files->collapse = collapse;
        if ( optind < argc ) args->files->require_index = 1;
        if ( args->subset_fname )
        {
            args->files->require_index = 1;
            if ( !bcf_sr_set_targets(args->files, args->subset_fname) )
                error("Failed to read the targets: %s\n", args->subset_fname);
        }
        if ( !bcf_sr_add_reader(args->files, fnames[i]) ) error("Failed to open or the file not indexed: %s\n", fnames[i]);
        for (k=optind; k<argc; k++) 
            if ( !bcf_sr_add_reader(args->files, argv[k]) ) error("Failed to open or the file not indexed: %s\n", argv[k]);
        init_data(args);
        if ( i==0 ) 
        {
            prev_samples = copy_header(args->header, args->samples, args->nsamples);
            prev_nsamples = args->nsamples;
        }
        else
        {
            args->print_header = 0;
            if ( compare_header(args->header, args->samples, args->nsamples, prev_samples, prev_nsamples) ) 
                error("Different samples in %s and %s\n", fnames[i-1],fnames[i]);
        }
        query_vcf(args);
        destroy_data(args);
        bcf_sr_destroy(args->files);
    }
    destroy_list(fnames, nfiles);
    destroy_list(prev_samples, prev_nsamples);
    free(args->format);
	free(args);
	return 0;
}


