#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2
typedef struct _token_t token_t;
struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    char *key;          // set only for string constants, otherwise NULL
    double threshold;   // filtering threshold
    int hdr_id;         // BCF header lookup ID
    int idx;            // 1-based index to VCF vector
    void (*setter)(bcf1_t *, struct _token_t *);

    // modified on filter evaluation at each VCF line
    double num_value; 
    char *str_value;
    int pass;           // -1 not applicable, 0 fails, >0 pass
    int missing_value;
};

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

typedef struct _args_t
{
    char *filter_str;
    token_t *filters, **flt_stack;  // input tokens (in RPN) and evaluation stack
    int nfilters, filter_logic;

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_fname;
    int argc, clevel, output_bcf, print_header, update_info, header_only, n_samples, *imap;
    int trim_alts, sites_only, known, novel, multiallelic, biallelic, exclude_ref, private, exclude_uncalled, min_ac, max_ac, calc_ac;
    char *fn_ref, *fn_out, **samples;
    char *include_types, *exclude_types;
    int include, exclude;
    htsFile *out;
}
args_t;


/* general filters - Djikstra's shunting yard algorithm */

#define TOK_VAL  0
#define TOK_LFT  1       // (
#define TOK_RGT  2       // )
#define TOK_LE   3       // less or equal
#define TOK_LT   4       // less than
#define TOK_EQ   5       // equal
#define TOK_BT   6       // bigger than
#define TOK_BE   7       // bigger or equal
#define TOK_OR   8       // |
#define TOK_AND  9       // &
#define TOK_ADD  10      // +
#define TOK_SUB  11      // -
#define TOK_MULT 12      // *
#define TOK_DIV  13      // /

static int op_prec[] = {0,5,5,5,5,5,5,5,1,2,6,6,7,7};

static int filters_next_token(char **str, int *len)
{
    char *tmp = *str;
    while ( *tmp && isspace(*tmp) ) tmp++;
    *str = tmp;
    *len = 0;

    while ( tmp[0] && tmp[1] )
    {
        if ( isspace(tmp[1]) ) break;
        if ( tmp[0]=='<' || tmp[1]=='<' ) break;
        if ( tmp[0]=='>' || tmp[1]=='>' ) break;
        if ( tmp[0]=='=' || tmp[1]=='=' ) break;
        if ( tmp[0]=='&' || tmp[1]=='&' ) break;
        if ( tmp[0]=='|' || tmp[1]=='|' ) break;
        if ( tmp[0]=='(' || tmp[1]=='(' ) break;
        if ( tmp[0]==')' || tmp[1]==')' ) break;
        if ( tmp[0]=='+' || tmp[1]=='-' ) break;
        if ( tmp[0]=='*' || tmp[1]=='/' ) break;
        tmp++;
    }
    if ( tmp > *str )
    {
        *len = tmp - (*str) + 1;
        return TOK_VAL;
    }
    if ( tmp[0]=='<' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_LE; }
        (*str) += 1; return TOK_LT;
    }
    if ( tmp[0]=='>' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_BE; }
        (*str) += 1; return TOK_BT;
    }
    if ( tmp[0]=='=' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_EQ; }
        (*str) += 1; return TOK_EQ;
    }
    if ( tmp[0]=='(' ) { (*str) += 1; return TOK_LFT; }
    if ( tmp[0]==')' ) { (*str) += 1; return TOK_RGT; }
    if ( tmp[0]=='&' && tmp[1]=='&' ) { (*str) += 2; return TOK_AND; }
    if ( tmp[0]=='|' && tmp[1]=='|' ) { (*str) += 2; return TOK_OR; }
    if ( tmp[0]=='&' ) { (*str) += 1; return TOK_AND; }
    if ( tmp[0]=='|' ) { (*str) += 1; return TOK_OR; }
    if ( tmp[0]=='+' ) { (*str) += 1; return TOK_ADD; }
    if ( tmp[0]=='-' ) { (*str) += 1; return TOK_SUB; }
    if ( tmp[0]=='*' ) { (*str) += 1; return TOK_MULT; }
    if ( tmp[0]=='/' ) { (*str) += 1; return TOK_DIV; }

    while ( *tmp && !isspace(*tmp) )
    {
        if ( *tmp=='<' ) break;
        if ( *tmp=='>' ) break;
        if ( *tmp=='=' ) break;
        if ( *tmp=='&' ) break;
        if ( *tmp=='|' ) break;
        if ( *tmp=='(' ) break;
        if ( *tmp==')' ) break;
        if ( *tmp=='+' ) break;
        if ( *tmp=='-' ) break;
        if ( *tmp=='*' ) break;
        if ( *tmp=='/' ) break;
        tmp++;
    }
    *len = tmp - (*str);
    return TOK_VAL;
}

static void filters_set_qual(bcf1_t *line, token_t *tok)
{
    if ( bcf_float_is_missing(line->qual) )     // hmm, how to do this cleanly and avoid the "strict-aliasing rules" warning?
        tok->missing_value = 1;
    else
        tok->num_value = line->qual;
}

/**
 *  bcf_get_info_value() - get single INFO value, int or float
 *  @line:      BCF line
 *  @info_id:   tag ID, as returned by bcf_id2int
 *  @ivec:      0-based index to retrieve
 *  @vptr:      pointer to memory location of sufficient size to accomodate
 *              info_id's type
 *
 *  The returned value is -1 if tag is not present, 0 if present but
 *  values is missing or ivec is out of range, and 1 on success.
 */
int bcf_get_info_value(bcf1_t *line, int info_id, int ivec, void *value)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == info_id ) break;
    if ( j==line->n_info ) return -1;

    bcf_info_t *info = &line->d.info[j];
    if ( info->len == 1 )
    {
        if ( info->type==BCF_HT_INT ) *((int*)value) = info->v1.i;
        else if ( info->type==BCF_HT_REAL ) *((float*)value) = info->v1.f;
        return 1;
    }

    #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
        type_t *p = (type_t *) info->vptr; \
        for (j=0; j<ivec && j<info->len; j++) \
        { \
            if ( is_vector_end ) return 0; \
        } \
        if ( is_missing ) return 0; \
        *((out_type_t*)value) = p[j]; \
        return 1; \
    }
    switch (info->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  int); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, int); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, int); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), float); break;
        default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
    }
    #undef BRANCH
    return -1;  // this shouldn't happen
}

static void filters_set_info_int(bcf1_t *line, token_t *tok)
{
    int value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
        tok->missing_value = 1;
    else
        tok->num_value = value;
}

static void filters_set_info_float(bcf1_t *line, token_t *tok)
{
    float value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
        tok->missing_value = 1;
    else
        tok->num_value = value;
}

static void filters_init1(args_t *args, char *str, int len, token_t *tok)
{
    tok->tok_type = TOK_VAL;
    tok->hdr_id   = -1;
    tok->pass     = -1;

    // is this a string constant?
    if ( str[0]=='"' )
    {
        if ( str[len-1] != '"' ) error("TODO: spaces in [%s]\n", args->filter_str);
        tok->key = (char*) calloc(len-1,sizeof(char));
        memcpy(tok->key,str+1,len-2);
        tok->key[len-2] = 0;
        return;
    }

    if ( !strncmp(str,"QUAL",len) )
    {
        tok->setter = filters_set_qual;
        tok->key = strdup("QUAL");
        return;
    }

    // is this one of the VCF tags? For now do only INFO and QUAL, to be extended...
    kstring_t tmp = {0,0,0};
    kputsn(str, len, &tmp);

    tok->hdr_id = bcf_id2int(args->hdr, BCF_DT_ID, tmp.s);
    if ( tok->hdr_id>=0 ) 
    {
        if ( tmp.s ) free(tmp.s);
        return;
    }

    // is it a substrict VCF vector tag?
    if ( tmp.s[tmp.l-1] == ']' )
    {
        int i;
        for (i=0; i<tmp.l; i++)
            if ( tmp.s[i]=='[' ) { tmp.s[i] = 0; break; }

        tok->hdr_id = bcf_id2int(args->hdr, BCF_DT_ID, tmp.s);
        if ( tok->hdr_id>=0 )
        {
            switch ( bcf_id2type(args->hdr,BCF_HL_INFO,tok->hdr_id) ) 
            {
                case BCF_HT_INT:  tok->setter = &filters_set_info_int; break;
                case BCF_HT_REAL: tok->setter = &filters_set_info_float; break;
                default: error("FIXME: not ready for this, sorry\n");
            }
            tok->idx = atoi(&tmp.s[i+1]);
            if ( tmp.s ) free(tmp.s);
            return;
        }
    }

    // is it a value?
    char *end;
    errno = 0;
    tok->threshold = strtod(tmp.s, &end);
    if ( errno!=0 || end==tmp.s ) error("Error: the tag \"INFO/%s\" is not defined in the VCF header\n", tmp.s);

    if ( tmp.s ) free(tmp.s);
}

void filters_debug(args_t *args, token_t *filters, int nfilters)
{
    if ( !filters ) 
    {
        filters  = args->filters;
        nfilters = args->nfilters;
    }
    fprintf(stderr,"RPN .. ");

    int i;
    for (i=0; i<nfilters; i++)
    {
        if ( filters[i].tok_type == TOK_VAL )
        {
            if ( filters[i].hdr_id >=0 )
                fprintf(stderr," %s", args->hdr->id[BCF_DT_ID][filters[i].hdr_id].key);
            else if ( filters[i].key )
                fprintf(stderr," %s", filters[i].key);
            else
                fprintf(stderr," %e", filters[i].threshold);
        }
        else
            fprintf(stderr," %c", "0()[<=>]|&"[filters[i].tok_type]);
    }
    fprintf(stderr, "\n");
}

void filters_debug_stack(args_t *args, int nstack)
{
    fprintf(stderr,"stack .. ");

    int i;
    for (i=0; i<nstack; i++)
    {
        if ( args->flt_stack[i]->pass < 0 )
        {
            if ( args->flt_stack[i]->str_value )
                fprintf(stderr," %s", args->flt_stack[i]->str_value);
            else
                fprintf(stderr," %e", args->flt_stack[i]->num_value);
        }
        else
            fprintf(stderr," %d", args->flt_stack[i]->pass);
    }
    fprintf(stderr,"\n");
}

// Parse filter expression and convert to reverse polish notation
static void filters_init(args_t *args)
{
    if ( !args->filter_str ) return;

    int nops = 0, mops = 0, *ops = NULL;    // operators stack
    int nout = 0, mout = 0;                 // filter tokens, RPN
    token_t *out = NULL;
    char *tmp = args->filter_str;
    while ( *tmp )
    {
        int len, ret;
        ret = filters_next_token(&tmp, &len);
        // fprintf(stderr,"%c .. [%s] %d\n", "x()[<=>]|&"[ret], tmp, len);
        if ( ret==TOK_LFT )         // left bracket
        {
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( ret==TOK_RGT )    // right bracket
        {
            while ( nops>0 && ops[nops-1]!=TOK_LFT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                out[nout-1].tok_type = ops[nops-1];
                nops--;
            }
            if ( nops<=0 ) error("Could not parse: %s\n", args->filter_str);
            nops--;
        }
        else if ( ret!=TOK_VAL )    // one of the comparison operators
        {
            while ( nops>0 && op_prec[ret] < op_prec[ops[nops-1]] )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                out[nout-1].tok_type = ops[nops-1];
                nops--;
            }
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( !len ) 
        {
            if ( !isspace(*tmp) ) error("Could not parse the expression: %s\n", args->filter_str);
            break;     // all tokens read
        }
        else                        // annotation name or filtering value
        {
            nout++;
            hts_expand0(token_t, nout, mout, out);
            filters_init1(args, tmp, len, &out[nout-1]);
            tmp += len;
        }
    }
    while ( nops>0 )
    {
        if ( ops[nops-1]==TOK_LFT || ops[nops-1]==TOK_RGT ) error("Could not parse the expression: %s\n", args->filter_str);
        nout++;
        hts_expand0(token_t, nout, mout, out);
        out[nout-1].tok_type = ops[nops-1];
        nops--;
    }

    if ( mops ) free(ops);
    args->filters   = out;
    args->nfilters  = nout;
    args->flt_stack = (token_t **)malloc(sizeof(token_t*)*nout);
}

static void filters_destroy(args_t *args)
{
    int i;
    for (i=0; i<args->nfilters; i++)
        if ( args->filters[i].key ) free(args->filters[i].key);
    if (args->filters) free(args->filters);
    if (args->flt_stack) free(args->flt_stack);
}

static int filters_pass(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_INFO);

    // filters_debug(args, NULL, 0);

    int i, j, nstack = 0;
    for (i=0; i<args->nfilters; i++)
    {
        args->filters[i].missing_value = 0;
        args->filters[i].str_value = NULL;
        args->filters[i].pass = -1;
        // filters_debug_stack(args, nstack);

        if ( args->filters[i].tok_type == TOK_VAL )
        {
            if ( args->filters[i].setter ) args->filters[i].setter(line, &args->filters[i]);
            else if ( args->filters[i].hdr_id >=0 ) 
            {
                for (j=0; j<line->n_info; j++)
                    if ( line->d.info[j].key == args->filters[i].hdr_id ) break;

                if ( j==line->n_info ) 
                    args->filters[i].missing_value = 1;
                else if ( line->d.info[j].type==BCF_BT_CHAR )
                    args->filters[i].str_value = (char*)line->d.info[j].vptr;
                else if ( line->d.info[j].type==BCF_BT_FLOAT )
                {
                    args->filters[i].num_value = line->d.info[j].v1.f;
                    args->filters[i].str_value = NULL;
                }
                else
                {
                    args->filters[i].num_value = line->d.info[j].v1.i;
                    args->filters[i].str_value = NULL;
                }
            }
            else if ( args->filters[i].key )
                args->filters[i].str_value = args->filters[i].key;
            else
                args->filters[i].num_value = args->filters[i].threshold;
            args->flt_stack[nstack++] = &args->filters[i];
            continue;
        }
        if ( nstack<2 ) 
            error("Error occurred while processing the filter \"%s\": too few values left on stack (%d)\n", args->filter_str,nstack);

        int is_str  = (args->flt_stack[nstack-1]->str_value ? 1 : 0) + (args->flt_stack[nstack-2]->str_value ? 1 : 0 );

        if ( args->filters[i].tok_type == TOK_OR )
        {
            if ( args->flt_stack[nstack-1]->pass<0 || args->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d OR)\n", args->filter_str,args->flt_stack[nstack-2]->pass,args->flt_stack[nstack-1]->pass);
            args->flt_stack[nstack-2]->pass = args->flt_stack[nstack-1]->pass + args->flt_stack[nstack-2]->pass;
            nstack--;
            continue;
        }
        if ( args->filters[i].tok_type == TOK_AND )
        {
            if ( args->flt_stack[nstack-1]->pass<0 || args->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d AND)\n", args->filter_str,args->flt_stack[nstack-2]->pass,args->flt_stack[nstack-1]->pass);
            args->flt_stack[nstack-2]->pass = args->flt_stack[nstack-1]->pass * args->flt_stack[nstack-2]->pass;
            nstack--;
            continue;
        }

        if ( args->filters[i].tok_type == TOK_ADD )
        {
            args->flt_stack[nstack-2]->num_value += args->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( args->filters[i].tok_type == TOK_SUB )
        {
            args->flt_stack[nstack-2]->num_value -= args->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( args->filters[i].tok_type == TOK_MULT )
        {
            args->flt_stack[nstack-2]->num_value *= args->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( args->filters[i].tok_type == TOK_DIV )
        {
            args->flt_stack[nstack-2]->num_value /= args->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }

        int is_true = 0;
        if ( args->flt_stack[nstack-1]->missing_value || args->flt_stack[nstack-2]->missing_value )
            is_true = 0;
        else if ( args->filters[i].tok_type == TOK_EQ )
        {
            if ( is_str==1 ) error("Comparing string to numeric value: %s\n", args->filter_str);
            if ( is_str==2 ) 
                is_true = strcmp(args->flt_stack[nstack-1]->str_value,args->flt_stack[nstack-2]->str_value) ? 0 : 1;
            else
                is_true = (args->flt_stack[nstack-1]->num_value == args->flt_stack[nstack-2]->num_value) ? 1 : 0;
        }
        else if ( is_str>0 ) error("Wrong operator in string comparison: %s\n", args->filter_str);
        else if ( args->filters[i].tok_type == TOK_LE )
            is_true = ( args->flt_stack[nstack-2]->num_value <= args->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( args->filters[i].tok_type == TOK_LT )
            is_true = ( args->flt_stack[nstack-2]->num_value <  args->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( args->filters[i].tok_type == TOK_EQ )
            is_true = ( args->flt_stack[nstack-2]->num_value == args->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( args->filters[i].tok_type == TOK_BT )
            is_true = ( args->flt_stack[nstack-2]->num_value >  args->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( args->filters[i].tok_type == TOK_BE )
            is_true = ( args->flt_stack[nstack-2]->num_value >= args->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else
            error("FIXME: did not expect this .. tok_type %d = %d\n", i, args->filters[i].tok_type);

        args->flt_stack[nstack-2]->pass = is_true;
        nstack--;
    }
    if ( nstack>1 ) error("Error occurred while processing the filter \"%s\": too many values left on stack (%d)\n", args->filter_str,nstack);
    if ( args->filter_logic==FLT_INCLUDE ) return args->flt_stack[0]->pass;
    return args->flt_stack[0]->pass ? 0 : 1;
}


void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd);

static void init_data(args_t *args)
{
    int i;
    args->hdr = args->files->readers[0].header;
    
    if (args->calc_ac && args->update_info)
    {
        bcf_hdr_append(args->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
        bcf_hdr_append(args->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    }
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "vcfsubset");

    // setup sample data    
    if (args->sample_names)
    {
        struct stat sbuf;
        if ( stat(args->sample_names, &sbuf) == 0  )
            args->samples = hts_readlines(args->sample_names, &args->n_samples);
        else
        {
            int m = 0, n = 0;
            char **s = 0;
            const char *q, *p;
            for (q = p = args->sample_names;; ++p) {
                if (*p == ',' || *p == 0) {
                    if (m == n) {
                        m = m? m<<1 : 16;
                        s = (char**)realloc(s, m * sizeof(void*));
                    }
                    s[n] = (char*)calloc(p - q + 1, 1);
                    strncpy(s[n++], q, p - q);
                    q = p + 1;
                    if (*p == 0) break;
                }
            }
            s = (char**)realloc(s, n * sizeof(void*));
            args->samples = s;
            args->n_samples = n;
        }
    }
    
    if (args->n_samples)
        args->imap = (int*)malloc(args->n_samples * sizeof(int));
    
    // determine variant types to include/exclude
    if (args->include_types || args->exclude_types) {
        if (args->include_types && args->exclude_types) {
            fprintf(stderr, "Error: only supply one of --include_types --exclude-types options\n");
            exit(1);
        }
        char **type_list = 0;
        int m = 0, n = 0;
        const char *q, *p;
        for (q = p = args->include_types ? args->include_types : args->exclude_types;; ++p) {
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    type_list = (char**)realloc(type_list, m * sizeof(void*));
                }
                type_list[n] = (char*)calloc(p - q + 1, 1);
                strncpy(type_list[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
        }
        type_list = (char**)realloc(type_list, n * sizeof(void*));

        if (args->include_types) {
            args->include = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->include |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->include |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->include |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->include |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    exit(1);
                }
            }
        }
        if (args->exclude_types) {
            args->exclude = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->exclude |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->exclude |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->exclude |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->exclude |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    exit(1);
                }
            }
        }
        for (i = 0; i < n; ++i)
            free(type_list[i]);
        free(type_list);
    }

    // setup output
    char modew[8];
    strcpy(modew, "w");
    if (args->clevel >= 0 && args->clevel <= 9) sprintf(modew + 1, "%d", args->clevel);
    if (args->output_bcf) strcat(modew, "b");
    args->out = hts_open(args->fn_out ? args->fn_out : "-", modew, 0);
    
    // headers: hdr=full header, hsub=subset header, hnull=sites only header
    if (args->sites_only)
        args->hnull = bcf_hdr_subset(args->hdr, 0, 0, 0);
    if (args->n_samples > 0)
        args->hsub = bcf_hdr_subset(args->hdr, args->n_samples, args->samples, args->imap);

    filters_init(args);
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->imap ) {
        for (i = 0; i < args->n_samples; ++i)
            free(args->samples[i]);
        free(args->samples);
        free(args->imap);
    }
    if (args->hnull) bcf_hdr_destroy(args->hnull);
    if (args->hsub) bcf_hdr_destroy(args->hsub);
    filters_destroy(args);
}

int subset_vcf(args_t *args, bcf1_t *line)
{
    if ( args->multiallelic && !(line->n_allele>2) ) return 0; // select multiallelic sites
    if ( args->biallelic && !(line->n_allele==2) ) return 0; // skip multiallelic sites
    if (args->novel || args->known)
    {
        if ( args->novel && (line->d.id[0]!='.' || line->d.id[1]!=0) ) return 0; // skip sites which are known, ID != '.'
        if ( args->known && line->d.id[0]=='.' && line->d.id[1]==0 ) return 0;  // skip sites which are novel, ID == '.'
    }
                
    if (args->include || args->exclude)
    {
        if ( args->include && !(line->d.var_type&args->include) ) return 0; // include only given variant types
        if ( args->exclude &&   line->d.var_type&args->exclude  ) return 0; // exclude given variant types
    }

    if ( args->filters && !filters_pass(args,line) ) return 0;
        

    int i, an = 0, n_ac = 0, *ac = (int*) calloc(line->n_allele+1,sizeof(int));
    if (args->calc_ac) {
        bcf_calc_ac(args->hdr, line, ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
        for (i=1; i<=line->n_allele; i++)
            n_ac += ac[i];
        for (i=0; i<line->n_allele; i++)
            an+=ac[i];
    }

    if (args->n_samples)
    {
        int n_ac_sub = 0, *ac_sub = (int*) calloc(line->n_allele+1,sizeof(int));
        bcf_subset(args->hdr, line, args->n_samples, args->imap);
        if (args->calc_ac) {
            bcf_calc_ac(args->hsub, line, ac_sub, BCF_UN_FMT); // recalculate AC and AN
            an = 0;
            for (i=0; i<line->n_allele; i++)
                an+=ac_sub[i];
            for (i=1; i<=line->n_allele; i++)
                n_ac_sub += ac_sub[i];
            if (args->private && !(n_ac_sub > 0 && n_ac == n_ac_sub)) { free(ac); free(ac_sub); return 0; }
            n_ac = n_ac_sub;
            for (i=0; i<=line->n_allele; i++)
                ac[i] = ac_sub[i];
        }
        free(ac_sub);
    }
    if (args->min_ac && args->min_ac>n_ac) { free(ac); return 0; }
    if (args->max_ac && args->max_ac<n_ac) { free(ac); return 0; }
    if (args->exclude_uncalled && n_ac == 0 && ac[0] == 0) { free(ac); return 0; }
    if (args->calc_ac && args->update_info) {
        bcf1_update_info_int32(args->hdr, line, "AC", &ac[1], line->n_allele-1);
        bcf1_update_info_int32(args->hdr, line, "AN", &an, 1);
    }
    free(ac);
    if (args->exclude_ref && n_ac == 0) return 0;
    if (args->trim_alts) bcf_trim_alleles(args->hsub ? args->hsub : args->hdr, line);
    if (args->sites_only) bcf_subset(args->hsub ? args->hsub : args->hdr, line, 0, 0);
bcf_unpack(line,BCF_UN_ALL);
line->d.shared_dirty |= BCF1_DIRTY_INF;
    if (args->output_bcf) bcf1_sync(line);
    return 1;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   View, subset and filter VCF/BCF files.\n");
    fprintf(stderr, "Usage:   vcfsubset [options] <in.bcf>|<in.vcf>|<in.vcf.gz> [region1 [...]]\n");
    fprintf(stderr, "\n");
    // fprintf(stderr, "Input options:\n");
    // fprintf(stderr, "    -S    input is VCF\n");
    // fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o, --out FILE        output file name [stdout]\n");
    // fprintf(stderr, "    -b                    output in BCF\n");
    // fprintf(stderr, "    -l INT                compression level [%d]\n", args->clevel);
    fprintf(stderr, "    -h                    suppress the header in VCF output\n");
    fprintf(stderr, "    -H                    print the header only\n");
    fprintf(stderr, "    -G,                   drop individual genotype information (after subsetting if -s option set)\n");
    fprintf(stderr, "    -t, --targets FILE    restrict to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Subset options:\n");
    fprintf(stderr, "    -a, --trim-alt-alleles      trim alternate alleles not seen in subset\n");
    fprintf(stderr, "    -I, --no-update             do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)\n");
    fprintf(stderr, "    -s, --samples STR/FILE      list of samples (FILE or comma separated list STR) [null]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -R,   --exclude-ref                            exclude sites without a non-reference genotype\n");
    fprintf(stderr, "    -U,   --exclude-uncalled                       exclude sites without a called genotype\n");
    fprintf(stderr, "    -p,   --private                                print sites where only the subset samples carry an non-reference allele\n");
    fprintf(stderr, "    -e,   --exclude-filters <expr>                 include sites for which the expression is true (e.g. 'QUAL>=10 && (DP4[2]+DP4[3] > 2)\n");
    fprintf(stderr, "    -i,   --include-filters <expr>                 same as -e but with the reverse logic\n");
    fprintf(stderr, "    -f,   --apply-filters <list>                   require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -k/n, --known/--novel                          print known/novel sites only (ID is/not '.')\n");
    fprintf(stderr, "    -m/M, --multiallelic/--biallelic               print multiallelic/biallelic sites only\n");
    fprintf(stderr, "    -c/C  --min-ac/--max-ac                        minimum/maximum allele count (INFO/AC) of sites to be printed\n");
    fprintf(stderr, "    -1/2  --singletons/--doubletons                print singleton/doubleton sites only (shortcut for -c1 -C1/-c2 -C2)\n");
    fprintf(stderr, "    -v/V  --include-types/--exclude-types STR      comma-separated list of variant types to include/exclude: snps,indels,mnps,other [null]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfsubset(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;

    static struct option loptions[] = 
    {
        {"exclude-filters",1,0,'e'},
        {"include-filters",1,0,'i'},
        {"trim-alt-alleles",0,0,'a'},
        {"exclude-ref",0,0,'R'},
        {"no-update",0,0,'I'},
        {"private",0,0,'p'},
        {"exclude-uncalled",0,0,'U'},
        {"apply-filters",1,0,'f'},
        {"known",0,0,'k'},
        {"novel",0,0,'n'},
        {"multiallelic",0,0,'m'},
        {"biallelic",0,0,'M'},
        {"samples",1,0,'s'},
        {"out",1,0,'o'},
        {"include-types",1,0,'v'},
        {"exclude-types",1,0,'V'},
        {"targets",1,0,'L'},
        {"min-ac",1,0,'c'},
        {"max-ac",1,0,'C'},
        {"singletons",0,0,'1'},
        {"doubletons",0,0,'2'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "l:bSt:o:L:s:Gf:knv:V:mMaRpUhHc:C:12Ie:i:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'l': args->clevel = atoi(optarg); args->output_bcf = 1; break;
            case 'b': args->output_bcf = 1; break;
            case 'o': args->fn_out = optarg; break;
            case 'h': args->print_header = 0; break;
            case 'H': args->header_only = 1; break;
            
            case 't': args->targets_fname = optarg; break;
            
            case 's': args->sample_names = optarg; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;
            
            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm': args->multiallelic = 1; break;
            case 'M': args->biallelic = 1; break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c': args->min_ac = atoi(optarg); args->calc_ac = 1; break;
            case 'C': args->max_ac = atoi(optarg); args->calc_ac = 1; break;
            case '1': args->min_ac = 1; args->max_ac = 1; args->calc_ac = 1; break;
            case '2': args->min_ac = 2; args->max_ac = 2; args->calc_ac = 1; break;

            case 'R': args->exclude_ref = 1; args->calc_ac = 1; break;
            case 'p': args->private = 1; args->calc_ac = 1; break;
            case 'U': args->exclude_uncalled = 1; args->calc_ac = 1; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Sorry, only one of -i or -e can be given.\n");
    if ( argc<optind+1 ) usage(args);   // none or too many files given
    // read in the regions from the command line
    if ( optind+1 < argc ) {
        args->files->require_index = 1;
        args->files->region = argv[optind+1];
        int i;
        for (i = optind + 1; i < argc; ++i) {
            if ( !args->files->seqs ) {
                args->files->mseqs = args->files->nseqs = 1;
                args->files->seqs = (const char**) malloc(sizeof(const char*));
                args->files->seqs[0] = argv[i];
            }
            else {
                args->files->mseqs += 30;
                args->files->seqs = (const char**) realloc(args->files->seqs, sizeof(const char*)*args->files->mseqs);
                args->files->seqs[args->files->nseqs++] = argv[i];
            }
        }
    }
    if ( args->targets_fname )
    {
        args->files->require_index = 1;
        if ( !bcf_sr_set_targets(args->files, args->targets_fname) )
            error("Failed to read the targets: %s\n", args->targets_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    bcf_hdr_t *out_hdr = args->hnull ? args->hnull : (args->hsub ? args->hsub : args->hdr);
    if (args->print_header)
        vcf_hdr_write(args->out, out_hdr);
    if (!args->header_only) {
        while ( bcf_sr_next_line(args->files) )
        {
            bcf1_t *line = args->files->readers[0].buffer[0];
            if ( subset_vcf(args, line) )
                vcf_write1(args->out, out_hdr, line);
        }
    }
    hts_close(args->out);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
