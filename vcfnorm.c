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
#include "faidx.h"

typedef struct
{
    int m,n,f;    // m: allocated size, n: number of elements in the buffer, f: first element
}
rbuf_t;

typedef struct
{
    int32_t dir:4, val:28;
}
cell_t;
typedef struct
{
    int nmat, nref, nseq;
    int ipos, lref, lseq;
    cell_t *mat;
    char *ref, *seq;
    int m_arr, *ipos_arr, *lref_arr, *lseq_arr;
}
aln_aux_t;

typedef struct
{
    aln_aux_t aln;
    char *tseq;
    int nseq;
    bcf1_t **lines;
    rbuf_t rbuf;
    int buf_win;            // maximum distance between two records to consider
    int aln_win;            // the realignment window size (maximum repeat size)
    bcf_srs_t *files;       // using the synced reader only for -r option
    bcf_hdr_t *hdr;
    faidx_t *fai;
	char **argv, *ref_fname, *vcf_fname;
	int argc, rmdup;
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
 *  rbuf_init() - initialize round buffer
 *  @rbuf:  the rbuf_t holder
 *  @size:  the maximum number of elements
 *
 */
inline void rbuf_init(rbuf_t *rbuf, int size)
{
    rbuf->m = size; rbuf->n = rbuf->f = 0;
}
/**
 *  rbuf_last() - get index of the last element of the round buffer
 *  @rbuf:  the rbuf_t holder
 *
 */
inline int rbuf_last(rbuf_t *rbuf)
{
    if ( !rbuf->n ) return -1;
    int i  = rbuf->n + rbuf->f - 1;
    if ( i >= rbuf->m ) i -= rbuf->m;
    return i;
}
/**
 *  rbuf_next() - get index of the next element in the round buffer
 *  @rbuf:  the rbuf_t holder
 *  @i:     pointer to the last rbuf index. Set to -1 before the first call.
 *
 *  Sets i to the next position in the buffer. The return value indicates if
 *  the position points to a valid element (1) or if there are no more elements
 *  after *i (0).
 */
inline int rbuf_next(rbuf_t *rbuf, int *i)
{
    if ( !rbuf->n ) return 0;
    if ( *i==-1 ) { *i = rbuf->f; return 1; }
    int n = (rbuf->f <= *i) ? *i - rbuf->f + 1 : *i + rbuf->m - rbuf->f + 1;
    if ( ++(*i) >= rbuf->m ) *i = 0;
    return n < rbuf->n ? 1 : 0;
}
/**
 *  rbuf_prev() - get index of the previous element in the round buffer
 *  @rbuf:  the rbuf_t holder
 *  @i:     pointer to the last rbuf index. Set to -1 before the first call.
 *
 *  Sets i to the previous position in the buffer. The return value indicates if
 *  the position points to a valid element (1) or if there are no more elements
 *  before *i (0).
 */
inline int rbuf_prev(rbuf_t *rbuf, int *i)
{
    if ( !rbuf->n || *i==rbuf->f ) return 0;
    if ( *i==-1 )
    {
        *i = rbuf_last(rbuf);
        return 1;
    }
    if ( --(*i) < 0 ) *i = rbuf->m - 1;
    return 1;
}
/**
 *  rbuf_add() - register new element in the round buffer
 *  @rbuf:  the rbuf_t holder
 *
 *  Returns index of the newly inserted element.
 */
inline int rbuf_add(rbuf_t *rbuf)
{
    if ( rbuf->n < rbuf->m )
    {
        rbuf->n++;
        int i = rbuf->f + rbuf->n;
        return i <= rbuf->m ? i - 1 : i - rbuf->m - 1;
    }

    rbuf->f++;
    if ( rbuf->f >= rbuf->m ) 
    {
        rbuf->f = 0;
        return rbuf->m - 1;
    }
    return rbuf->f - 1;
}
/**
 *  rbuf_shift() - removes first element of the buffer
 *  @rbuf:  the rbuf_t holder
 *
 *  Returns index of the removed element.
 */
inline int rbuf_shift(rbuf_t *rbuf)
{
    if ( !rbuf->n ) return -1;
    int ret = rbuf->f;
    rbuf->f++;
    if ( rbuf->f >= rbuf->m ) rbuf->f = 0;
    rbuf->n--;
    return ret;
}



static void align(args_t *args, aln_aux_t *aux)
{
    // Needleman-Wunsch global alignment
    char *ref = aux->ref;
    char *seq = aux->seq;
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    int nlen  = nref>nseq ? nref : nseq;
    if ( (nref+1)*(nseq+1) > aux->nmat )
    {
        aux->nmat = (nref+1)*(nseq+1);
        aux->mat  = (cell_t *) realloc(aux->mat, sizeof(cell_t)*aux->nmat);
    }
    cell_t *mat = aux->mat;
    int i, j, k = nref+2, kd = nref+2;
    for (j=0; j<=nref; j++) { mat[j].val = 0; mat[j].dir = 0; }
    for (i=1; i<=nseq; i++)
    {
        mat[k-1].val = 0;
        mat[k-1].dir = 0;
        int jmax = i-1 < nref ? i-1 : nref;
        for (j=1; j<=jmax; j++)
        {
            // prefer insertions to deletions and mismatches
            int max, dir;
            if ( ref[j-1]==seq[i-1] )
            {
                max = mat[k-kd].val + 1; dir = 0;                                           // match
                if ( max < mat[k-kd+1].val - 1 ) { max = mat[k-kd+1].val - 1; dir = 1; }    // insertion
                if ( max < mat[k-1].val - 1 ) { max = mat[k-1].val - 1; dir = -1; }         // deletion

            }
            else
            {
                max = mat[k-kd+1].val - 1; dir = 1;                                         // insertion
                if ( max < mat[k-1].val - 1 ) { max = mat[k-1].val - 1; dir = -1; }         // deletion
                if ( max < mat[k-kd].val - 1 ) { max = mat[k-kd].val - 1; dir = 0; }        // mismatch
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        for (j=jmax+1; j<=nref; j++)
        {
            // prefer deletions to insertions and mismatches
            int max, dir;
            if ( ref[j-1]==seq[i-1] )
            {
                max = mat[k-kd].val + 1; dir = 0;                                           // match
                if ( max < mat[k-1].val - 1 ) { max = mat[k-1].val - 1; dir = -1; }         // deletion
                if ( max < mat[k-kd+1].val - 1 ) { max = mat[k-kd+1].val - 1; dir = 1; }    // insertion

            }
            else
            {
                max = mat[k-1].val - 1; dir = -1;                                           // deletion
                if ( max < mat[k-kd+1].val - 1 ) { max = mat[k-kd+1].val - 1; dir = -1; }   // insertion
                if ( max < mat[k-kd].val - 1 ) { max = mat[k-kd].val - 1; dir = 0; }        // mismatch
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        k++;
    }
    k = (nref+1)*(nseq+1)-1;
    int kmin = nref>nseq ? 2*(nref+1) - nseq : (nseq-nref)*(nref+1);
    while (k>kmin && !mat[k].dir) k -= kd;
    i = k/(nref+1);
    j = k - i*(nref+1);
    assert(i>0 && j>0);
    int l = k, ialn = nlen-1, nout_ref = 0, nout_seq = 0, ipos = 0;
    while ( l>0 )
    {
        if ( j<=0 || mat[l].dir==1 )    // i
        {
            ipos = 0;
            nout_seq++;
            l -= kd - 1;
            i--;
        }
        else if ( i<=0 || mat[l].dir==-1 )  // d
        {
            ipos = 0;
            nout_ref++;
            l--;
            j--;
        }
        else     // m
        {
            ipos = ref[i]==seq[j] ? ipos+1 : 1;
            nout_seq++;
            nout_ref++;
            l -= kd;
            i--;
            j--;
        }
        ialn--;
    }
    assert(ipos>0);
    aux->ipos = ( nout_ref>ipos && nout_seq>ipos && ref[ipos-1]==seq[ipos-1] ) ?  ipos : ipos-1;
    aux->lref = nout_ref;
    aux->lseq = nout_seq;
}

void realign(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_STR);

    int i, ref_len = strlen(line->d.allele[0]), len = ref_len;
    for (i=1; i<line->n_allele; i++)
    {
        int l = strlen(line->d.allele[i]);
        if ( len < l ) len = l;
    }
    if ( len==1 ) return;    // SNP

    int win = line->pos < args->aln_win ? line->pos - 1 : args->aln_win;
    len += win + 2;
    if ( args->nseq < len ) 
    {
        args->nseq = len;
        args->tseq = (char*) realloc(args->tseq, sizeof(char)*len);
    }
    char *ref  = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win, line->pos+ref_len, &len);
    assert(ref[win]==line->d.allele[0][0]);

    if ( args->aln.m_arr < line->n_allele )
    {
        args->aln.m_arr = line->n_allele;
        args->aln.ipos_arr = (int*) realloc(args->aln.ipos_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lref_arr = (int*) realloc(args->aln.lref_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lseq_arr = (int*) realloc(args->aln.lseq_arr, sizeof(int)*args->aln.m_arr);
    }
    int *ipos = args->aln.ipos_arr;
    int *lref = args->aln.lref_arr;
    int *lseq = args->aln.lseq_arr;
    int min_ipos = INT_MAX, max_lref = 0, max_lseq = 0;

    for (i=0; i<win; i++)
        args->tseq[i] = ref[i];
    int j, k;
    for (j=1; j<line->n_allele; j++)
    {
        char *t = line->d.allele[j];
        k = i;
        while (*t) { args->tseq[k++] = *t; t++; }
        args->tseq[k++]   = ref[len-1];
        args->tseq[k]     = 0;

        args->aln.ref  = ref;
        args->aln.seq  = args->tseq;
        args->aln.nref = len;
        args->aln.nseq = k;

        align(args, &args->aln);

        ipos[j] = args->aln.ipos;
        lref[j] = args->aln.lref;
        lseq[j] = args->aln.lseq - args->aln.ipos;
        if ( min_ipos > ipos[j] ) min_ipos = ipos[j];
        if ( max_lref < lref[j] ) max_lref = lref[j];
        if ( max_lseq < lseq[j] ) max_lseq = lseq[j];
    }

    line->pos += min_ipos - win;

    char *rmme = line->d.als;
    kstring_t str; str.l = 0; str.m = max_lref-args->aln.ipos;
    if ( str.m<max_lseq ) str.m = max_lseq;
    str.m = line->n_allele*(str.m+1);
    str.s = (char*) malloc(sizeof(char)*str.m);
    line->d.m_als = str.m;
    line->d.allele[0] = line->d.als = str.s;
    char *t = str.s;
    for (k=min_ipos; k<max_lref; k++) { *t = ref[k]; t++; } 
    *t = 0; t++;
    for (k=1; k<line->n_allele; k++)
    {
        char *p = t;
        int nprefix = ipos[k] - min_ipos;
        int nshift  = win - min_ipos;   
        int nsuffix = max_lref - lref[k];
        for (j=0; j<nprefix; j++) { *t = ref[j+min_ipos]; t++; }
        for (j=0; j<lseq[k] && j<nshift; j++) { *t = ref[j+min_ipos+nprefix]; t++; }
        nshift = lseq[k] - nshift;
        for (j=0; j<nshift; j++) { *t = line->d.allele[k][j+nprefix]; t++; }    // warning: this may be buggy
        for (j=0; j<nsuffix; j++) { *t = ref[lref[k]+j]; t++; }
        *t = 0; t++;
        line->d.allele[k] = p;
    }
    free(rmme);
    free(ref);
}

void flush_buffer(args_t *args, htsFile *file, int n)
{
    int i, k, prev_rid = -1, prev_pos = 0, prev_type = 0;
    for (i=0; i<n; i++)
    {
        k = rbuf_shift(&args->rbuf);
        // todo: merge with next record if POS and the type are same. For now, just discard if asked to do so.
        if ( args->rmdup )
        {
            bcf_set_variant_types(args->lines[k]);
            if ( prev_rid>=0 && prev_rid==args->lines[k]->rid && prev_pos==args->lines[k]->pos && prev_type==args->lines[k]->d.var_type )
                continue;
            prev_rid  = args->lines[k]->rid;
            prev_pos  = args->lines[k]->pos;
            prev_type = args->lines[k]->d.var_type;
        }
        vcf_write1(file, args->hdr, args->lines[k]);
    }
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    rbuf_init(&args->rbuf, 100);
    args->lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
    args->fai = fai_load(args->ref_fname);
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->rbuf.m; i++)
        if ( args->lines[i] ) bcf_destroy1(args->lines[i]);
    free(args->lines);
    fai_destroy(args->fai);
    if ( args->nseq ) free(args->tseq);
    if ( args->aln.nmat ) free(args->aln.mat);
    if ( args->aln.m_arr ) { free(args->aln.ipos_arr); free(args->aln.lref_arr); free(args->aln.lseq_arr); }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void normalize_vcf(args_t *args)
{
    htsFile *out = hts_open("-","w",0);
    vcf_hdr_write(out, args->hdr);

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];

        // still on the same chromosome?
        int i, j, ilast = rbuf_last(&args->rbuf); 
        if ( ilast>=0 && line->rid != args->lines[ilast]->rid ) flush_buffer(args, out, args->rbuf.n); // new chromosome

        realign(args, line);

        // insert into sorted buffer
        i = j = ilast = rbuf_add(&args->rbuf);
        if ( !args->lines[i] ) args->lines[i] = bcf_init1();
        SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->lines[i]);
        while ( rbuf_prev(&args->rbuf,&i) )
        {
            if ( args->lines[i]->pos > args->lines[j]->pos ) SWAP(bcf1_t*, args->lines[i], args->lines[j]);
            j = i;
        }

        // find out how many sites to flush
        j = 0;
        for (i=-1; rbuf_next(&args->rbuf,&i); )
        {
            if ( args->lines[ilast]->pos - args->lines[i]->pos < args->buf_win ) break;
            j++;
        }
        if ( j==args->rbuf.m ) j = 1;
        if ( j>0 ) flush_buffer(args, out, j);
    }
    flush_buffer(args, out, args->rbuf.n);
    hts_close(out);
}

static void usage(void)
{
	fprintf(stderr, "About:   Left-align and normalize indels.\n");
	fprintf(stderr, "Usage:   vcfnorm [options] -f ref.fa <file.vcf.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -D, --remove-duplicates           remove duplicate lines of the same type. [Todo: merge genotypes, don't just throw away.]\n");
	fprintf(stderr, "    -f, --fasta-ref <file>            reference sequence\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    perform intersection in the given region only\n");
	fprintf(stderr, "    -w, --win <int,int>               alignment window and buffer window [50,1000]\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfnorm(int argc, char *argv[])
{
	int c;
	args_t *args  = (args_t*) calloc(1,sizeof(args_t));
	args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->aln_win = 50;
    args->buf_win = 1000;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"fasta-ref",1,0,'f'},
		{"region",1,0,'r'},
		{"win",1,0,'w'},
		{"remove-duplicates",0,0,'D'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hr:f:w:D",loptions,NULL)) >= 0) {
		switch (c) {
			case 'D': args->rmdup = 1; break;
			case 'f': args->ref_fname = optarg; break;
			case 'r': args->files->region = optarg; break;
            case 'w': { if (sscanf(optarg,"%d,%d",&args->aln_win,&args->buf_win)!=2) error("Could not parse --win %s\n", optarg); break; }
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if ( argc!=optind+1 || !args->ref_fname ) usage();   // none or too many files given
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    init_data(args);
    normalize_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
	free(args);
	return 0;
}

