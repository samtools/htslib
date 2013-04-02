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
    char *tseq, *seq;
    int mseq;
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


void _vcfnorm_debug_print(aln_aux_t *aux)
{
    cell_t *mat = aux->mat;
    char *ref = aux->ref;
    char *seq = aux->seq; 
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    int nlen  = nref>nseq ? nref : nseq;
    int k     = (nref+1)*(nseq+1)-1;
    int kd    = nref+2;
    int i = k/(nref+1);
    int j = k - i*(nref+1);
    assert(i>0 && j>0);
    int l = k, ialn = 0, nout_ref = 0, nout_seq = 0, ipos = 0;
    char *aln_ref = (char*) malloc(sizeof(char)*(nlen+1));
    char *aln_seq = (char*) malloc(sizeof(char)*(nlen+1));
    while ( l>0 )
    {
        if ( j<=0 || mat[l].dir==1 )    // i
        {
            aln_ref[ialn] = '-';
            aln_seq[ialn] = seq[nseq-i]; 
            ipos = 0;
            nout_seq++;
            l -= kd - 1;
            i--;
        }
        else if ( i<=0 || mat[l].dir==-1 )  // d
        {
            aln_ref[ialn] = ref[nref-j]; 
            aln_seq[ialn] = '-'; 
            ipos = 0;
            nout_ref++;
            l--;
            j--;
        }
        else     // m
        {
            aln_ref[ialn] = ref[nref-j];
            aln_seq[ialn] = seq[nseq-i]; 
            ipos = ref[nref-j+1]==seq[nseq-i+1] ? ipos+1 : 1;
            nout_seq++;
            nout_ref++;
            l -= kd;
            i--;
            j--;
        }
        ialn++;
    }
    aln_ref[ialn] = aln_seq[ialn] = 0;
    fprintf(stderr, "ref: %s\n", ref);
    fprintf(stderr, "seq: %s\n", seq); 
    fprintf(stderr, "-> %s\n", aln_ref);
    fprintf(stderr, "-> %s\n", aln_seq); 
    free(aln_ref);
    free(aln_seq);

    fprintf(stderr, "      ");
    for (j=0; j<nref; j++) fprintf(stderr, "   %c ", ref[nref-j-1]); fprintf(stderr, "\n"); 
    for (i=0; i<=nseq; i++)
    {
        fprintf(stderr, "%c", i==0 ? ' ' : seq[nseq-i]);
        for (j=0; j<=nref; j++)
        {
            char dir = ' ';
            if ( mat[i*(nref+1)+j].dir==1 ) dir = 'i';
            else if ( mat[i*(nref+1)+j].dir==-1 ) dir = 'd';
            fprintf(stderr, " %3d%c", (int)mat[i*(nref+1)+j].val, dir);
        }
        fprintf(stderr,"\n");
    } 
}

static int align(args_t *args, aln_aux_t *aux)
{
    // Needleman-Wunsch global alignment. Note that the sequences are aligned from
    //  the end where matches are preferred, gaps are pushed to the front (left-aligned)
    char *ref = aux->ref;
    char *seq = aux->seq;
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    if ( (nref+1)*(nseq+1) > aux->nmat )
    {
        aux->nmat = (nref+1)*(nseq+1);
        aux->mat  = (cell_t *) realloc(aux->mat, sizeof(cell_t)*aux->nmat);
        if ( !aux->mat ) 
            error("Could not allocate %ld bytes of memory at %d\n", sizeof(cell_t)*aux->nmat, args->files->readers[0].buffer[0]->pos+1);
    }
    const int GAP_OPEN = -1, GAP_CLOSE = -1, GAP_EXT = 0, MATCH = 1, MISM = -1, DI = 1, DD = -1, DM = 0;
    cell_t *mat = aux->mat;
    int i, j, k = nref+2, kd = nref+2;
    mat[0].val = 20; mat[0].dir = DM;   // the last ref and alt bases match
    for (j=1; j<=nref; j++) { mat[j].val = 0; mat[j].dir = DM; }
    for (i=1; i<=nseq; i++)
    {
        mat[k-1].val = 0;
        mat[k-1].dir = DM;
        int jmax = i-1 < nref ? i-1 : nref;
        for (j=1; j<=jmax; j++)
        {
            // prefer insertions to deletions and mismatches
            int max, dir, score;
            if ( ref[nref-j]==seq[nseq-i] )
            {
                // match
                max = mat[k-kd].val + MATCH;
                if ( mat[k-kd].dir!=DM ) max += GAP_CLOSE;
                dir = DM;

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT: mat[k-kd+1].val + GAP_OPEN;
                if ( max < score )  { max = score; dir = DI; }

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN; 
                if ( max < score ) { max = score; dir = DD; }
            }
            else
            {
                // insertion
                max = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN; 
                dir = DI; 

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DD; }

                // mismatch
                score = mat[k-kd].val + MISM;
                if ( mat[k-kd].dir!=DM ) score += GAP_CLOSE;
                if ( max < score ) { max = score; dir = DM; }
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        for (j=jmax+1; j<=nref; j++)
        {
            // prefer deletions to insertions and mismatches
            int max, dir, score;
            if ( ref[nref-j]==seq[nseq-i] )
            {
                // match
                max = mat[k-kd].val + MATCH;
                if ( mat[k-kd].dir!=DM ) max += GAP_CLOSE;
                dir = DM; 

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT: mat[k-1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DD; }

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN;
                if ( max < score )  { max = score; dir = DI; }
            }
            else
            {
                // deletion
                max = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN; 
                dir = DD;

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DI; }

                // mismatch
                score = mat[k-kd].val + MISM;
                if ( mat[k-kd].dir!=DM ) score += GAP_CLOSE;
                if ( max < score ) { max = score; dir = DM; }
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        k++;
    }

    // _vcfnorm_debug_print(aux);

    // skip as much of the matching sequence at the beggining as possible
    k = (nref+1)*(nseq+1)-1;
    int kmin = nref>nseq ? 2*(nref+1) - nseq : (nseq-nref)*(nref+1);
    int ipos = 0;
    while (k>kmin && !mat[k].dir) { k -= kd; ipos++; }

    i = k/(nref+1);
    j = k - i*(nref+1);
    assert(i>0 && j>0);
    int l = k, nout_ref = ipos, nout_seq = ipos, nsuffix = 0;
    while ( l>0 )
    {
        if ( j<=0 || mat[l].dir==DI )    // insertion
        {
            nsuffix = 0;
            nout_seq++;
            l -= kd - 1;
            i--;
        }
        else if ( i<=0 || mat[l].dir==DD )  // deletion
        {
            nsuffix = 0;
            nout_ref++;
            l--;
            j--;
        }
        else     // match/mismatch
        {
            nsuffix = ref[nref-i]==seq[nseq-j] ? nsuffix + 1 : 0;
            nout_seq++;
            nout_ref++;
            l -= kd;
            i--;
            j--;
        }
    }
    if ( !ipos ) return -1; // the window is too small

    aux->ipos = ipos - 1;
    aux->lref = nout_ref - nsuffix;
    aux->lseq = nout_seq - nsuffix;

    // The indels and complex events do not have to be padded
    if ( aux->lref - aux->ipos > 1 && aux->lseq - aux->ipos > 1 && ref[aux->ipos]==seq[aux->ipos] ) aux->ipos++;
    return 0;
}

int realign(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_STR);

    int i, ref_len = strlen(line->d.allele[0]), len = ref_len;
    for (i=1; i<line->n_allele; i++)
    {
        int l = strlen(line->d.allele[i]);
        if ( len < l ) len = l;
    }
    if ( len==1 ) return 0;    // SNP

    // Sanity check: exclude broken VCFs with long stretches of N's
    if ( len>1000 )
    {
        for (i=0; i<ref_len; i++)
            if ( line->d.allele[0][i]=='N' ) return -1;
    }

    int win = line->pos < args->aln_win ? line->pos - 1 : args->aln_win;
    len += win + 2;
    if ( args->mseq < len*(line->n_allele-1) ) 
    {
        args->mseq = len*(line->n_allele-1);
        args->seq  = (char*) realloc(args->seq, sizeof(char)*args->mseq);
    }
    int ref_winlen;
    char *ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win, line->pos+ref_len, &ref_winlen);
    assert( ref_winlen==ref_len+win+1 );

    // Sanity check: the reference sequence must match the REF allele
    if ( strncasecmp(&ref[win],line->d.allele[0],ref_len) ) 
        error("\nSanity check failed, the reference sequence differs at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);

    if ( args->aln.m_arr < line->n_allele )
    {
        args->aln.m_arr = line->n_allele;
        args->aln.ipos_arr = (int*) realloc(args->aln.ipos_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lref_arr = (int*) realloc(args->aln.lref_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lseq_arr = (int*) realloc(args->aln.lseq_arr, sizeof(int)*args->aln.m_arr);
    }
    int *ipos = args->aln.ipos_arr;
    int *iref = args->aln.lref_arr;
    int *iseq = args->aln.lseq_arr;
    int min_pos = INT_MAX, max_ref = 0;

    int j, k;
    for (j=1; j<line->n_allele; j++)
    {
        // get the ALT ready for alignment
        args->tseq = args->seq + (j-1)*len;
        for (i=0; i<win; i++) args->tseq[i] = ref[i];
        char *t = line->d.allele[j];
        while (*t) { args->tseq[i++] = *t; t++; }
        args->tseq[i++] = ref[ref_winlen-1];
        args->tseq[i]   = 0;

        args->aln.ref  = ref;
        args->aln.seq  = args->tseq;
        args->aln.nref = ref_winlen;
        args->aln.nseq = i;

        if ( align(args, &args->aln)<0 ) 
        {
            // something went wrong - output the original line
            free(ref);
            return 0;
        }

        // fprintf(stderr, "%s  \t nref=%d\n", ref, ref_winlen);
        // fprintf(stderr, "%s  \t nseq=%d\n", args->tseq, i);
        // fprintf(stderr, "pos=%d win=%d  ipos=%d lref=%d lseq=%d\n", line->pos+1, win, args->aln.ipos, args->aln.lref, args->aln.lseq);
        // fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lref; k++) fprintf(stderr, "%c", ref[k]); fprintf(stderr, "\n");
        // fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lseq; k++) fprintf(stderr, "%c", args->tseq[k]); fprintf(stderr, "\n");
        // fprintf(stderr, "\n"); 

        ipos[j] = args->aln.ipos;   // position before the first difference (w.r.t. the window)
        iref[j] = args->aln.lref;   // length of the REF alignment (or the index after the last aligned position)
        iseq[j] = args->aln.lseq;
        if ( max_ref < iref[j] ) max_ref = iref[j];
        if ( min_pos > ipos[j] ) min_pos = ipos[j];
        assert( iseq[j]<=len );
    }

    // Check if the record's position must be changed
    int nmv = win - min_pos;
    assert( nmv>=0 );   // assuming that it will never align more to the right
    line->pos -= nmv;

    // todo: 
    //      - modify the alleles only if needed. For now redoing always to catch errors
    //      - in some cases the realignment does not really improve things, see the case at 2:114 in test/norm.vcf

    // REF
    kstring_t str = {0,0,0};
    kputsn_(&ref[min_pos], max_ref-min_pos, &str); 
    kputc_(0, &str);
    // ALTs
    for (k=1; k<line->n_allele; k++)
    {
        // prefix the sequence with REF bases if the other alleles were aligned more to the left
        int nprefix = ipos[k] - min_pos;
        if ( nprefix ) kputsn_(&ref[min_pos], nprefix, &str);

        // the ALT sequence 
        int nseq = iseq[k] - ipos[k];
        if ( nseq )
        {
            char *alt = args->seq + (k-1)*len + ipos[k];
            kputsn_(alt, nseq, &str);
        }

        // suffix invoked by other deletions which must be added to match the REF
        int nsuffix = max_ref - iref[k];
        if ( nsuffix ) kputsn_(&ref[iref[k]], nsuffix, &str);
        kputc_(0, &str);
    }
    // create new block of alleles 
    char *rmme = line->d.als;
    line->d.allele[0] = line->d.als = str.s;
    line->d.m_als = str.m;
    char *t = str.s;
    for (k=1; k<line->n_allele; k++)
    {
        while (*t) t++;
        line->d.allele[k] = ++t;
    }
    free(rmme);
    free(ref);
    return 1;
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
    if ( args->mseq ) free(args->seq);
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
        if ( realign(args, line)<0 ) continue;   // exclude broken VCF lines

        // still on the same chromosome?
        int i, j, ilast = rbuf_last(&args->rbuf); 
        if ( ilast>=0 && line->rid != args->lines[ilast]->rid ) flush_buffer(args, out, args->rbuf.n); // new chromosome

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
        if ( args->rbuf.n==args->rbuf.m ) j = 1;
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

