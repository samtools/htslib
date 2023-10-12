/*
    Copyright (C) 2018-2023 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

/*
    Note: The original code comes from https://github.com/pd3/utils/tree/master/annot-regs.
    In annot-tsv we changed the naming from "destination" file to "target" file, however
    the code still internally uses the original naming convention.
*/

#include <config.h>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include "htslib/hts.h"
#include "htslib/hts_defs.h"
#include "htslib/khash_str2int.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/regidx.h"

#define ANN_NBP     1
#define ANN_FRAC    2
#define ANN_CNT     4

typedef struct
{
    uint32_t n,m;
    char **off, *rmme;
}
cols_t;

typedef struct
{
    void *name2idx;
    cols_t *cols, *annots;
    int dummy;
}
hdr_t;

typedef struct
{
    char *fname;
    hdr_t hdr;
    cols_t *core, *match, *transfer, *annots;
    int *core_idx, *match_idx, *transfer_idx, *annots_idx;
    int *nannots_added; // for --max-annots: the number of annotations added
    int grow_n;
    kstring_t line;     // one buffered line, a byproduct of reading the header
    htsFile *fp;
}
dat_t;

// This is for the special -a annotations, keeps a list of
// source regions that hit the destination region. The start
// coordinates are converted to beg<<1 and end coordinates
// to (end<<1)+1.
#define NBP_SET_BEG(x) ((x)<<1)
#define NBP_SET_END(x) (((x)<<1)+1)
#define NBP_GET(x)     ((x)>>1)
#define NBP_IS_BEG(x)  (((x)&1)==0)
#define NBP_IS_END(x)  (((x)&1)==1)
typedef struct
{
    size_t n,m;          // n is a multiple of two: breakpoints are stored in regs, not regions
    hts_pos_t *regs;
    hts_pos_t beg,end;   // the current destination interval
}
nbp_t;

#define PRINT_MATCHING    1
#define PRINT_NONMATCHING 2
typedef struct
{
    nbp_t *nbp;
    dat_t dst, src;
    char *core_str, *match_str, *transfer_str, *annots_str;
    char *temp_dir, *out_fname;
    BGZF *out_fp;
    int allow_dups, reciprocal, ignore_headers, max_annots, mode;
    double overlap;
    regidx_t *idx;
    regitr_t *itr;
    kstring_t tmp_kstr;
    cols_t *tmp_cols;               // the -t transfer fields to write for each line
    khash_t(str2int) **tmp_hash;    // lookup tables for tmp_cols
}
args_t;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error(const char *format, ...)
{
    va_list ap;
    fflush(stdout);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

static nbp_t *nbp_init(void)
{
    nbp_t *nbp = calloc(1,sizeof(nbp_t));
    if ( !nbp ) error("Out of memory, failed to allocatae %zu bytes\n",sizeof(nbp_t));
    return nbp;
}
static void nbp_destroy(nbp_t *nbp)
{
    free(nbp->regs);
    free(nbp);
}
static inline void nbp_reset(nbp_t *nbp, hts_pos_t beg, hts_pos_t end)
{
    nbp->n   = 0;
    nbp->beg = beg;
    nbp->end = end;
}
static inline void nbp_add(nbp_t *nbp, hts_pos_t beg, hts_pos_t end)
{
    nbp->n += 2;
    if ( nbp->n >= nbp->m )
    {
        nbp->m += 2;
        nbp->regs = realloc(nbp->regs, nbp->m*sizeof(*nbp->regs));
        if ( !nbp->regs ) error("Out of memory, failed to allocate %zu bytes\n",nbp->m*sizeof(*nbp->regs));
    }
    nbp->regs[nbp->n - 2] = NBP_SET_BEG(beg);
    nbp->regs[nbp->n - 1] = NBP_SET_END(end);
}
static int compare_hts_pos(const void *aptr, const void *bptr)
{
    hts_pos_t a = *(const hts_pos_t*) aptr;
    hts_pos_t b = *(const hts_pos_t*) bptr;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}
static hts_pos_t nbp_length(nbp_t *nbp)
{
    qsort(nbp->regs, nbp->n, sizeof(*nbp->regs), compare_hts_pos);
    int i, nopen = 0;
    hts_pos_t beg = 0, length = 0;
    for (i=0; i<nbp->n; i++)
    {
        if ( NBP_IS_BEG(nbp->regs[i]) )
        {
            if ( !nopen ) beg = NBP_GET(nbp->regs[i]);
            nopen++;
        }
        else nopen--;
        assert( nopen>=0 );
        if ( nopen==0 && beg>0 ) length += NBP_GET(nbp->regs[i]) - beg + 1;
    }
    return length;
}

cols_t *cols_split(const char *line, cols_t *cols, char delim)
{
    if ( !cols ) cols = (cols_t*) calloc(1,sizeof(cols_t));
    if ( !cols ) error("Out of memory, failed to allocate %zu bytes\n",sizeof(cols_t));
    if ( cols->rmme ) free(cols->rmme);
    cols->n = 0;
    cols->rmme = strdup(line);
    if ( !cols->rmme ) error("Out of memory\n");
    char *ss = cols->rmme;
    while (1)
    {
        char *se = ss;
        while ( *se && *se!=delim ) se++;
        char tmp = *se;
        *se = 0;
        cols->n++;
        if ( cols->n > cols->m )
        {
            cols->m += 10;
            cols->off = realloc(cols->off, sizeof(*cols->off)*cols->m);
            if ( !cols->off ) error("Out of memory, failed to allocate %zu bytes\n",sizeof(*cols->off)*cols->m);
        }
        cols->off[ cols->n - 1 ] = ss;
        if ( !tmp ) break;
        ss = se + 1;
    }
    return cols;
}
// Can be combined with cols_split() but is much slower.
// The string must exist throughout the life of cols unless initialized with cols_split().
void cols_append(cols_t *cols, char *str)
{
    if ( cols->rmme )
    {
        size_t str_len = strlen(str);
        size_t lst_len = strlen(cols->off[ cols->n - 1 ]);
        size_t tot_len = 2 + str_len + lst_len + (cols->off[ cols->n - 1 ] - cols->rmme);

        cols_t *tmp_cols = (cols_t*)calloc(1,sizeof(cols_t));
        if ( !tmp_cols ) error("Out of memory, failed to allocate %zu bytes\n",sizeof(cols_t));
        tmp_cols->rmme = calloc(tot_len,1);
        tmp_cols->off  = calloc(cols->n+1,sizeof(*tmp_cols->off));
        if ( !tmp_cols->rmme || !tmp_cols->off ) error("Out of memory\n");

        char *ptr = tmp_cols->rmme;
        int i;
        for (i=0; i<cols->n; i++)
        {
            size_t len = strlen(cols->off[i]);
            memcpy(ptr, cols->off[i], len);
            tmp_cols->off[i] = ptr;
            ptr += len + 1;
        }
        memcpy(ptr, str, str_len);
        tmp_cols->off[i] = ptr;

        free(cols->off);
        free(cols->rmme);
        cols->rmme = tmp_cols->rmme;
        cols->off  = tmp_cols->off;
        cols->n    = cols->n+1;
        cols->m    = cols->n;
        free(tmp_cols);
        return;
    }
    cols->n++;
    if ( cols->n > cols->m )
    {
        cols->m++;
        cols->off = realloc(cols->off,sizeof(*cols->off)*cols->m);
        if ( !cols->off ) error("Out of memory, failed to allocate %zu bytes\n",sizeof(*cols->off)*cols->m);
    }
    cols->off[cols->n-1] = str;
}
void cols_clear(cols_t *cols)
{
    if ( !cols ) return;
    free(cols->rmme);
    free(cols->off);
    cols->rmme = NULL;
    cols->off  = NULL;
}
void cols_destroy(cols_t *cols)
{
    if ( !cols ) return;
    cols_clear(cols);
    free(cols);
}

int parse_tab_with_payload(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    static int beg_end_warned = 0;

    if ( line[0]=='#' )
    {
        *((cols_t**)payload) = NULL;
        return -1;
    }

    dat_t *dat = (dat_t*) usr;

    cols_t *cols = cols_split(line, NULL, '\t');
    *((cols_t**)payload) = cols;

    if ( cols->n < dat->core_idx[0] ) error("Expected at least %d columns, found %d: %s\n",dat->core_idx[0]+1,cols->n,line);
    *chr_beg = cols->off[ dat->core_idx[0] ];
    *chr_end = *chr_beg + strlen(*chr_beg) - 1;

    if ( cols->n < dat->core_idx[1] ) error("Expected at least %d columns, found %d: %s\n",dat->core_idx[1]+1,cols->n,line);
    char *tmp, *ptr = cols->off[ dat->core_idx[1] ];
    *beg = strtod(ptr, &tmp);
    if ( tmp==ptr ) error("Expected numeric value, found \"%s\": %s\n",ptr,line);

    if ( cols->n < dat->core_idx[2] ) error("Expected at least %d columns, found %d: %s\n",dat->core_idx[2]+1,cols->n,line);
    ptr = cols->off[ dat->core_idx[2] ];
    *end = strtod(ptr, &tmp);
    if ( tmp==ptr ) error("Expected numeric value, found \"%s\": %s\n",ptr,line);

    if ( *end < *beg )
    {
        if ( !beg_end_warned )
            fprintf(stderr,"Warning: the start coordinate is bigger than the end coordinate:\n\t%s\nThis message is printed only once.\n",line);
        beg_end_warned = 1;
        hts_pos_t tmp = *beg; *beg = *end; *end = tmp;
    }

    return 0;
}
void free_payload(void *payload)
{
    cols_t *cols = *((cols_t**)payload);
    cols_destroy(cols);
}

// Parse header if present (first line has a leading #) or create a dummy header with
// numeric column names. If dummy is set, read first data line (without a leading #)
// and create a dummy header.
void parse_header(dat_t *dat, char *fname, int dummy)
{
    dat->fp = hts_open(fname,"r");
    if ( !dat->fp ) error("Failed to open: %s\n", fname);

    cols_t *cols = NULL;
    while ( hts_getline(dat->fp, KS_SEP_LINE, &dat->line) > 0 )
    {
        if ( dat->line.s[0]=='#' )
        {
            // this is a header or comment line
            if ( dummy ) continue;
            cols = cols_split(dat->line.s, NULL, '\t');
            break;
        }

        // this a data line, we must be in a dummy mode
        cols = cols_split(dat->line.s, NULL, '\t');
        assert(cols && cols->n);
        assert(cols->off[0][0] != '#');

        // create a dummy header with numeric field names
        kstring_t str = {0,0,0};
        int i, n = cols->n;
        for (i=0; i<n; i++)
        {
            if ( i>0 ) kputc('\t', &str);
            kputw(i+1, &str);
        }
        cols_destroy(cols);
        cols = cols_split(str.s, NULL, '\t');
        free(str.s);
        dat->hdr.dummy = 1;

        break;
    }
    if ( !dat->line.l ) error("Failed to read: %s\n", fname);
    assert(cols && cols->n);

    dat->hdr.name2idx = khash_str2int_init();
    int i;
    for (i=0; i<cols->n; i++)
    {
        char *ss = cols->off[i];
        while ( *ss && (*ss=='#' || isspace(*ss)) ) ss++;
        if ( !*ss ) error("Could not parse the header field \"%s\": %s\n", cols->off[i],dat->line.s);
        if ( *ss=='[' )
        {
            char *se = ss+1;
            while ( *se && isdigit(*se) ) se++;
            if ( *se==']' ) ss = se + 1;
        }
        while ( *ss && (*ss=='#' || isspace(*ss)) ) ss++;
        if ( !*ss ) error("Could not parse the header field \"%s\": %s\n", cols->off[i],dat->line.s);
        cols->off[i] = ss;
        khash_str2int_set(dat->hdr.name2idx, cols->off[i], i);
    }
    dat->hdr.cols = cols;
    if ( !dat->hdr.dummy ) dat->line.l = 0;
}
void write_header(args_t *args, dat_t *dat)
{
    if ( dat->hdr.dummy ) return;
    int i;
    kstring_t str = {0,0,0};
    kputc('#', &str);
    for (i=0; i<dat->hdr.cols->n; i++)
    {
        if ( i>0 ) kputc('\t', &str);
        ksprintf(&str,"[%d]", i+1);
        kputs(dat->hdr.cols->off[i], &str);
    }
    if ( dat->hdr.annots )
    {
        for (i=0; i<dat->hdr.annots->n; i++)
        {
            if ( str.l > 1 ) kputc('\t', &str);
            kputs(dat->hdr.annots->off[i], &str);
        }
    }
    kputc('\n',&str);
    if ( bgzf_write(args->out_fp, str.s, str.l) != str.l ) error("Failed to write %zd bytes\n", str.l);
    free(str.s);
}
void destroy_header(dat_t *dat)
{
    if ( dat->hdr.cols ) cols_destroy(dat->hdr.cols);
    khash_str2int_destroy(dat->hdr.name2idx);
}

static int read_next_line(dat_t *dat)
{
    if ( dat->line.l ) return dat->line.l;
    int ret = hts_getline(dat->fp, KS_SEP_LINE, &dat->line);
    if ( ret > 0 ) return dat->line.l;
    if ( ret < -1 ) error("Error encountered while reading %s\n",dat->fname);
    return 0;
}

void sanity_check_columns(char *fname, hdr_t *hdr, cols_t *cols, int **col2idx, int force)
{
    *col2idx = (int*)malloc(sizeof(int)*cols->n);
    if ( !*col2idx ) error("Out of memory, failed to allocate %zu bytes\n",sizeof(int)*cols->n);
    int i, idx;
    for (i=0; i<cols->n; i++)
    {
        if ( khash_str2int_get(hdr->name2idx, cols->off[i], &idx) < 0 )
        {
            if ( !force ) error("The key \"%s\" not found in %s\n", cols->off[i],fname);
            idx = -1;
        }
        (*col2idx)[i] = idx;
    }
}
void init_data(args_t *args)
{
    parse_header(&args->dst, args->dst.fname, args->ignore_headers);
    parse_header(&args->src, args->src.fname, args->ignore_headers);

    // -c, core columns
    if ( !args->core_str ) args->core_str = "chr,beg,end:chr,beg,end";
    cols_t *tmp = cols_split(args->core_str, NULL, ':');
    args->src.core = cols_split(tmp->off[0],NULL,',');
    args->dst.core = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
    sanity_check_columns(args->src.fname, &args->src.hdr, args->src.core, &args->src.core_idx, 0);
    sanity_check_columns(args->dst.fname, &args->dst.hdr, args->dst.core, &args->dst.core_idx, 0);
    if ( args->src.core->n!=3 || args->dst.core->n!=3 ) error("Expected three columns: %s\n", args->core_str);
    cols_destroy(tmp);

    // -m, match columns
    if ( args->match_str )
    {
        tmp = cols_split(args->match_str, NULL, ':');
        args->src.match = cols_split(tmp->off[0],NULL,',');
        args->dst.match = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
        sanity_check_columns(args->src.fname, &args->src.hdr, args->src.match, &args->src.match_idx, 0);
        sanity_check_columns(args->dst.fname, &args->dst.hdr, args->dst.match, &args->dst.match_idx, 0);
        if ( args->src.match->n != args->dst.match->n ) error("Expected equal number of columns: %s\n", args->match_str);
        cols_destroy(tmp);
    }

    // -t, transfer columns
    int i;
    if ( args->transfer_str )
    {
        tmp = cols_split(args->transfer_str, NULL, ':');
        args->src.transfer = cols_split(tmp->off[0],NULL,',');
        args->dst.transfer = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
        sanity_check_columns(args->src.fname, &args->src.hdr, args->src.transfer, &args->src.transfer_idx, 1);
        sanity_check_columns(args->dst.fname, &args->dst.hdr, args->dst.transfer, &args->dst.transfer_idx, 1);
        if ( args->src.transfer->n != args->dst.transfer->n ) error("Expected equal number of columns: %s\n", args->transfer_str);
        for (i=0; i<args->src.transfer->n; i++)
        {
            if ( args->src.transfer_idx[i]==-1 )
            {
                cols_append(args->src.hdr.cols,args->src.transfer->off[i]);
                args->src.transfer_idx[i] = -args->src.hdr.cols->n;    // negative index indicates different ptr location
                args->src.grow_n++;
            }
        }
        for (i=0; i<args->dst.transfer->n; i++)
        {
            if ( args->dst.transfer_idx[i]==-1 )
            {
                cols_append(args->dst.hdr.cols,args->dst.transfer->off[i]);
                args->dst.transfer_idx[i] = args->dst.hdr.cols->n - 1;
                args->dst.grow_n++;
            }
        }
        args->tmp_cols = (cols_t*)calloc(args->src.transfer->n,sizeof(cols_t));
        args->tmp_hash = (khash_t(str2int)**)calloc(args->src.transfer->n,sizeof(khash_t(str2int)*));
        if ( !args->tmp_cols || !args->tmp_hash ) error("Out of memory\n");
        for (i=0; i<args->src.transfer->n; i++)
            args->tmp_hash[i] = khash_str2int_init();
        cols_destroy(tmp);
    }
    else
    {
        args->src.transfer = calloc(1,sizeof(*args->src.transfer));
        if ( !args->src.transfer ) error("Out of memory\n");
    }
    args->src.nannots_added = calloc(args->src.transfer->n,sizeof(*args->src.nannots_added));
    if ( !args->src.nannots_added ) error("Out of memory\n");

    // -a, annotation columns
    if ( args->annots_str )
    {
        tmp = cols_split(args->annots_str, NULL, ':');
        args->src.annots = cols_split(tmp->off[0],NULL,',');
        args->dst.annots = cols_split(tmp->n==2 ? tmp->off[1] : tmp->off[0],NULL,',');
        if ( args->src.annots->n!=args->dst.annots->n ) error("Different number of src and dst columns in %s\n",args->annots_str);
        args->dst.annots_idx = (int*) malloc(sizeof(int)*args->dst.annots->n);
        if ( !args->dst.annots_idx ) error("Out of memory\n");
        for (i=0; i<args->src.annots->n; i++)
        {
            if ( !strcasecmp(args->src.annots->off[i],"nbp") )
            {
                args->dst.annots_idx[i] = ANN_NBP;
                cols_append(args->dst.hdr.cols,"nbp");
            }
            else if ( !strcasecmp(args->src.annots->off[i],"frac") )
            {
                args->dst.annots_idx[i] = ANN_FRAC;
                cols_append(args->dst.hdr.cols,"frac");
            }
            else if ( !strcasecmp(args->src.annots->off[i],"cnt") )
            {
                args->dst.annots_idx[i] = ANN_CNT;
                cols_append(args->dst.hdr.cols,"cnt");
            }
            else error("The annotation \"%s\" is not recognised\n", args->src.annots->off[i]);
        }
        args->nbp = nbp_init();
        cols_destroy(tmp);
    }

    args->idx = regidx_init(NULL, parse_tab_with_payload,free_payload,sizeof(cols_t),&args->src);
    while ( read_next_line(&args->src) )
    {
        if ( regidx_insert(args->idx,args->src.line.s) !=0 ) error("Could not parse the region in %s: %s\n",args->src.fname,args->src.line.s);
        args->src.line.l = 0;
    }
    args->itr = regitr_init(args->idx);
    if ( hts_close(args->src.fp)!=0 ) error("Failed to close: %s\n", args->src.fname);

    int len = args->out_fname ? strlen(args->out_fname) : 0;
    if ( len )
    {
        int compress_output = 0;
        if ( !strcasecmp(".gz",args->out_fname+len-3) || !strcasecmp(".bgz",args->out_fname+len-4) ) compress_output = 1;
        args->out_fp = bgzf_open(args->out_fname, compress_output ? "wg" : "wu");
    }
    else
        args->out_fp = bgzf_open("-","wu");
    if ( !args->out_fp ) error("Could not open file for writing: %s\n",args->out_fname?args->out_fname:"stdout");
}
void destroy_data(args_t *args)
{
    if ( bgzf_close(args->out_fp)!=0 ) error("Failed to close: %s\n", args->out_fname?args->out_fname:"stdout");
    if ( hts_close(args->dst.fp)!=0 ) error("Failed to close: %s\n", args->dst.fname);
    int i;
    for (i=0; i<args->src.transfer->n; i++)
        khash_str2int_destroy(args->tmp_hash[i]);
    free(args->tmp_hash);
    for (i=0; i<args->src.transfer->n; i++) cols_clear(&args->tmp_cols[i]);
    free(args->tmp_cols);
    cols_destroy(args->src.core);
    cols_destroy(args->dst.core);
    cols_destroy(args->src.match);
    cols_destroy(args->dst.match);
    cols_destroy(args->src.transfer);
    cols_destroy(args->dst.transfer);
    if ( args->src.annots ) cols_destroy(args->src.annots);
    if ( args->dst.annots ) cols_destroy(args->dst.annots);
    if ( args->nbp ) nbp_destroy(args->nbp);
    destroy_header(&args->src);
    destroy_header(&args->dst);
    free(args->src.nannots_added);
    free(args->src.core_idx);
    free(args->dst.core_idx);
    free(args->src.match_idx);
    free(args->dst.match_idx);
    free(args->src.transfer_idx);
    free(args->dst.transfer_idx);
    free(args->src.annots_idx);
    free(args->dst.annots_idx);
    free(args->src.line.s);
    free(args->dst.line.s);
    if (args->itr) regitr_destroy(args->itr);
    if (args->idx) regidx_destroy(args->idx);
    free(args->tmp_kstr.s);
}

static inline void write_string(args_t *args, char *str, size_t len)
{
    if ( len==0 ) len = strlen(str);
    if ( len==0 ) str = ".", len = 1;
    if ( bgzf_write(args->out_fp, str, len) != len ) error("Failed to write %zd bytes\n", len);
}
static void write_annots(args_t *args)
{
    if ( !args->dst.annots ) return;

    args->tmp_kstr.l = 0;
    int i;
    hts_pos_t len = nbp_length(args->nbp);
    for (i=0; i<args->dst.annots->n; i++)
    {
        if ( args->dst.annots_idx[i]==ANN_NBP )
        {
            kputc('\t',&args->tmp_kstr);
            kputw(len,&args->tmp_kstr);
        }
        else if ( args->dst.annots_idx[i]==ANN_FRAC )
        {
            kputc('\t',&args->tmp_kstr);
            kputd((double)len/(args->nbp->end - args->nbp->beg + 1),&args->tmp_kstr);
        }
        else if ( args->dst.annots_idx[i]==ANN_CNT )
        {
            kputc('\t',&args->tmp_kstr);
            kputw(args->nbp->n/2,&args->tmp_kstr);
        }
    }
    write_string(args, args->tmp_kstr.s, args->tmp_kstr.l);
}

void process_line(args_t *args, char *line, size_t size)
{
    char *chr_beg, *chr_end;
    hts_pos_t beg, end;
    cols_t *dst_cols = NULL;
    int i,j;
    int ret = parse_tab_with_payload(line, &chr_beg, &chr_end, &beg, &end, &dst_cols, &args->dst);
    if ( ret==-1 )
    {
        cols_destroy(dst_cols);
        return;
    }

    if ( args->nbp ) nbp_reset(args->nbp,beg,end);

    if ( !regidx_overlap(args->idx, chr_beg,beg,end, args->itr) )
    {
        if ( args->mode & PRINT_NONMATCHING )
        {
            write_string(args, line, size);
            write_annots(args);
            write_string(args, "\n", 1);
        }
        cols_destroy(dst_cols);
        return;
    }

    for (i=0; i<args->src.transfer->n; i++)
    {
        args->src.nannots_added[i] = 0;
        args->tmp_cols[i].n = 0;
        kh_clear(str2int, args->tmp_hash[i]);
    }

    int has_match = 0, annot_len = 0;
    while ( regitr_overlap(args->itr) )
    {
        if ( args->overlap )
        {
            double len1 = end - beg + 1;
            double len2 = args->itr->end - args->itr->beg + 1;
            double isec = (args->itr->end < end ? args->itr->end : end) - (args->itr->beg > beg ? args->itr->beg : beg) + 1;
            if ( args->reciprocal )
            {
                if ( isec/len1 < args->overlap || isec/len2 < args->overlap ) continue;
            }
            else
            {
                if ( isec/len1 < args->overlap && isec/len2 < args->overlap ) continue;
            }
        }
        cols_t *src_cols = regitr_payload(args->itr,cols_t*);
        if ( args->dst.match && args->dst.match->n )
        {
            for (i=0; i<args->dst.match->n; i++)
            {
                if ( args->dst.match_idx[i] > dst_cols->n ) error("Expected at least %d columns, found %d: %s\n",args->dst.match_idx[i],dst_cols->n,line);
                char *dst = dst_cols->off[ args->dst.match_idx[i] ];
                char *src = src_cols->off[ args->src.match_idx[i] ];
                if ( strcmp(dst,src) ) break;
            }
            if ( i != args->dst.match->n ) continue;
        }
        has_match = 1;

        if ( args->nbp )
            nbp_add(args->nbp, args->itr->beg >= beg ? args->itr->beg : beg, args->itr->end <= end ? args->itr->end : end);

        int max_annots_reached = 0;
        for (i=0; i<args->src.transfer->n; i++)
        {
            char *str;
            if ( args->src.transfer_idx[i] >= 0 )
                str = src_cols->off[ args->src.transfer_idx[i] ];                   // transfer a value from the src file
            else
                str = args->src.hdr.cols->off[ -args->src.transfer_idx[i] - 1 ];    // non-existent field in src, use a default value

            if ( !str || !*str ) str = ".";     // provide a nill dot value when the field is empty

            if ( !args->allow_dups )
            {
                if ( khash_str2int_has_key(args->tmp_hash[i],str) ) continue;
                khash_str2int_set(args->tmp_hash[i],str,1);
            }
            if ( args->max_annots )
            {
                if ( ++args->src.nannots_added[i] >= args->max_annots ) max_annots_reached = 1;
            }
            cols_append(&args->tmp_cols[i], str);
            annot_len += strlen(str);
        }
        if ( max_annots_reached ) break;
    }

    if ( !has_match )
    {
        if ( args->mode & PRINT_NONMATCHING )
        {
            write_string(args, line, size);
            write_annots(args);
            write_string(args, "\n", 1);
        }
        cols_destroy(dst_cols);
        return;
    }
    if ( !(args->mode & PRINT_MATCHING) )
    {
        cols_destroy(dst_cols);
        return;
    }

    size_t len;
    args->tmp_kstr.l = 0;
    ks_resize(&args->tmp_kstr, annot_len*3 + args->src.transfer->n*2);
    for (i=0; i<args->src.transfer->n; i++)
    {
        char *off = dst_cols->off[ args->dst.transfer_idx[i] ] = args->tmp_kstr.s + args->tmp_kstr.l;
        cols_t *ann = &args->tmp_cols[i];
        if ( !ann->n ) { off[0] = '.'; off[1] = 0; args->tmp_kstr.l += 2; continue; }
        for (j=0; j<ann->n; j++)
        {
            if ( j>0 ) { off[0] = ','; off++; args->tmp_kstr.l++; }
            len = strlen(ann->off[j]);
            memcpy(off, ann->off[j], len);
            off += len;
            args->tmp_kstr.l += len;
        }
        off[0] = 0;
        args->tmp_kstr.l++;
    }
    write_string(args, dst_cols->off[0], 0);
    for (i=1; i<dst_cols->n; i++)
    {
        write_string(args, "\t", 1);
        write_string(args, dst_cols->off[i], 0);
    }
    write_annots(args);
    write_string(args, "\n", 1);
    cols_destroy(dst_cols);
}

static const char *usage_text(void)
{
    return
        "About: Annotate regions of the target file (TGT) with information from overlapping regions\n"
        "       of the source file (SRC). Multiple columns can be transferred (-f) and the transfer\n"
        "       can be conditioned on requiring matching values in one or more columns (-m).\n"
        "       In addition to column transfer (-f) and special annotations (-a), the program can\n"
        "       operate in a simple grep-like mode and print matching lines (when neither -f nor -a\n"
        "       are given) or drop matching lines (-x).\n"
        "       All indexes and coordinates are 1-based and inclusive.\n"
        "\n"
        "Usage: annot-tsv [OPTIONS] -s source.txt -t target.txt > output.txt\n"
        "\n"
        "Common options:\n"
        "   -c, --core SRC:TGT              Core columns in SRC and TGT file [chr,beg,end:chr,beg,end]\n"
        "   -f, --transfer SRC:TGT          Columns to transfer. If SRC column does not exist, interpret\n"
        "                                   as the default value to use. If the TGT column does not exist,\n"
        "                                   a new column is created. If the TGT column does exist, its values\n"
        "                                   are overwritten when overlap is found or left as is otherwise.\n"
        "   -m, --match SRC:TGT             Require match in these columns for annotation transfer\n"
        "   -o, --output FILE               Output file name [STDOUT]\n"
        "   -s, --source-file FILE          Source file to take annotations from\n"
        "   -t, --target-file FILE          Target file to be extend with annotations from -s\n"
        "\n"
        "Other options:\n"
        "       --allow-dups                Add annotations multiple times\n"
        "       --max-annots INT            Adding at most INT annotations per column to save time in big regions\n"
        "       --version                   Print version string and exit\n"
        "   -a, --annotate LIST             Add special annotations, one or more of:\n"
        "                                       cnt  .. number of overlapping regions\n"
        "                                       frac .. fraction of the target region with an overlap\n"
        "                                       nbp  .. number of source base pairs in the overlap\n"
        "   -H, --ignore-headers            Use numeric indexes, ignore the headers completely\n"
        "   -O, --overlap FLOAT             Minimum required overlap (non-reciprocal, unless -r is given)\n"
        "   -r, --reciprocal                Apply the -O requirement to both overlapping intervals\n"
        "   -x, --drop-overlaps             Drop overlapping regions (precludes -f)\n"
        "\n"
        "Examples:\n"
        "   # Header is present, match and transfer by column name\n"
        "   annot-tsv -s src.txt.gz -t tgt.txt.gz -c chr,beg,end:CHR,POS,POS -m type,sample:TYPE,SMPL -f info:INFO\n"
        "\n"
        "   # Header is not present, match and transfer by column index (1-based)\n"
        "   annot-tsv -s src.txt.gz -t tgt.txt.gz -c 1,2,3:1,2,3 -m 4,5:4,5 -f 6:6\n"
        "\n"
        "   # If the TGT part is not given, the program assumes that the SRC:TGT columns are identical\n"
        "   annot-tsv -s src.txt.gz -t tgt.txt.gz -c chr,beg,end -m type,sample -f info\n"
        "\n"
        "   # One of the SRC or TGT file can be streamed from stdin\n"
        "   gunzip -c src.txt.gz | annot-tsv -t tgt.txt.gz -c chr,beg,end -m type,sample -f info\n"
        "   gunzip -c tgt.txt.gz | annot-tsv -s src.txt.gz -c chr,beg,end -m type,sample -f info\n"
        "\n"
        "   # Print matching regions as above but without modifying the records\n"
        "   gunzip -c src.txt.gz | annot-tsv -t tgt.txt.gz -c chr,beg,end -m type,sample\n"
        "\n";
}

int main(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    static struct option loptions[] =
    {
        {"core",required_argument,NULL,'c'},
        {"transfer",required_argument,NULL,'f'},
        {"match",required_argument,NULL,'m'},
        {"output",required_argument,NULL,'o'},
        {"source-file",required_argument,NULL,'s'},
        {"target-file",required_argument,NULL,'t'},
        {"allow-dups",no_argument,NULL,0},
        {"max-annots",required_argument,NULL,2},
        {"version",no_argument,NULL,1},
        {"annotate",required_argument,NULL,'a'},
        {"ignore-headers",no_argument,NULL,'H'},
        {"overlap",required_argument,NULL,'O'},
        {"reciprocal",no_argument,NULL,'r'},
        {"drop-overlaps",no_argument,NULL,'x'},
        {NULL,0,NULL,0}
    };
    char *tmp = NULL;
    int c;
    while ((c = getopt_long(argc, argv, "hc:f:m:o:s:t:a:HO:rx",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  0 : args->allow_dups = 1; break;
            case  1 : printf("%s\n",hts_version()); return 0; break;
            case  2 :
                args->max_annots = strtod(optarg, &tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse --max-annots  %s\n", optarg);
                break;
            case 'H': args->ignore_headers = 1; break;
            case 'r': args->reciprocal = 1; break;
            case 'c': args->core_str  = optarg; break;
            case 't': args->dst.fname = optarg; break;
            case 'm': args->match_str = optarg; break;
            case 'a': args->annots_str = optarg; break;
            case 'o': args->out_fname = optarg; break;
            case 'O':
                args->overlap = strtod(optarg, &tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse --overlap %s\n", optarg);
                if ( args->overlap<0 || args->overlap>1 ) error("Expected value from the interval [0,1]: --overlap %s\n", optarg);
                break;
            case 's': args->src.fname = optarg; break;
            case 'f': args->transfer_str = optarg; break;
            case 'x': args->mode = PRINT_NONMATCHING; break;
            case 'h':
            case '?': printf("\nVersion: %s\n%s\n",hts_version(),usage_text()); exit(EXIT_SUCCESS); break;
            default: error("\nVersion: %s\n%s\n",hts_version(),usage_text()); break;
        }
    }
    if ( argc==1 ) error("\nVersion: %s\n%s\n",hts_version(),usage_text());
    if ( !args->dst.fname && !args->src.fname ) error("Missing the -s and -t options\n");
    if ( !args->dst.fname || !args->src.fname )
    {
        if ( isatty(fileno((FILE *)stdin)) ) error("Missing the %s option\n",args->dst.fname?"-s":"-t");
        // reading from stdin
        if ( !args->dst.fname ) args->dst.fname = "-";
        if ( !args->src.fname ) args->src.fname = "-";
    }
    if ( !args->mode )
    {
        if ( !args->transfer_str && !args->annots_str ) args->mode = PRINT_MATCHING;
        else args->mode = PRINT_MATCHING|PRINT_NONMATCHING;
    }
    if ( (args->transfer_str || args->annots_str) && !(args->mode & PRINT_MATCHING) ) error("The option -x cannot be combined with -f and -a\n");

    init_data(args);
    write_header(args, &args->dst);
    while ( read_next_line(&args->dst) )
    {
        int i;
        for (i=0; i<args->dst.grow_n; i++) kputs("\t.", &args->dst.line);
        process_line(args, args->dst.line.s, args->dst.line.l);
        args->dst.line.l = 0;
    }
    destroy_data(args);
    free(args);

    return 0;
}

