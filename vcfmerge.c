/*
    Known issues:
        - shared block needs to be updated on BCF output (BCF output currently not supported)
        - Number=A,G tags not treated
 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

#include "khash.h"
KHASH_MAP_INIT_STR(strdict, int)
typedef khash_t(strdict) strdict_t;

#define SKIP_DONE 1
#define SKIP_DIFF 2

#define IS_VL_G(hdr,id) (bcf_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_G)
#define IS_VL_A(hdr,id) (bcf_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_A)

// Auxiliary merge data for selecting the right combination
//  of buffered records across multiple readers. maux1_t 
//  corresponds to one buffered line.
typedef struct
{
    int skip;
    int *map;   // mapping from input alleles to the output array
    int mmap;   // size of map array (only buffer[i].n_allele is actually used)
    int als_differ;
}
maux1_t;
typedef struct
{
    int n;          // number of readers
    char **als;     // merged output alleles
    int nals, mals; // size of the output array
    int *cnt, ncnt; // number of records that refer to the alleles
    int *nbuf;      // readers have buffers of varying lengths
    int *smpl_ploidy, *smpl_nGsize; // ploidy and derived number of values in Number=G tags, updated for each line (todo: cache for missing cases)
    int *flt, mflt, minf;
    bcf_info_t *inf;// out_line's INFO fields
    bcf_fmt_t **fmt_map; // i-th output FORMAT field corresponds in j-th reader to i*nreader+j, first row is reserved for GT
    int nfmt_map;        // number of rows in the fmt_map array
    void *tmp_arr;
    int ntmp_arr;
    maux1_t **d;    // d[i][j] i-th reader, j-th buffer line
    bcf_srs_t *files;
    int *has_line;  // which files are being merged
}
maux_t;

typedef struct
{
    maux_t *maux;
    int header_only, collapse;
    char *header_fname;
    strdict_t *tmph;
    kstring_t tmps;
    bcf_srs_t *files;
    bcf1_t *out_line;
    htsFile *out_fh;
    bcf_hdr_t *out_hdr;
    char **argv;
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

int bcf_hdr_sync(bcf_hdr_t *h);
void bcf_hdr_add_sample(bcf_hdr_t *h, char *s);

void bcf_hdr_merge(bcf_hdr_t *hw, const bcf_hdr_t *_hr, const char *clash_prefix)
{
    bcf_hdr_t *hr = (bcf_hdr_t*)_hr;

    // header lines
    int i, nw_ori = hw->nhrec;
    for (i=0; i<hr->nhrec; i++)
    {
        if ( hr->hrec[i]->type==BCF_HL_GEN && hr->hrec[i]->value )
        {
            int j;
            for (j=0; j<nw_ori; j++)
            {
                if ( hw->hrec[j]->type!=BCF_HL_GEN ) continue;
                if ( !strcmp(hr->hrec[i]->key,hw->hrec[j]->key) && !strcmp(hr->hrec[i]->value,hw->hrec[j]->value) ) break;
            }
            if ( j>=nw_ori )
                bcf_hdr_add_hrec(hw, bcf_hrec_dup(hr->hrec[i]));
        }
        else
        {
            bcf_hrec_t *rec = bcf_hdr_get_hrec(hw, hr->hrec[i]->type, hr->hrec[i]->vals[0]);
            if ( !rec )
                bcf_hdr_add_hrec(hw, bcf_hrec_dup(hr->hrec[i]));
        }
    }

    // samples
    for (i=0; i<hr->n[BCF_DT_SAMPLE]; i++)
    {
        char *name = strdup(hr->samples[i]);
        if ( bcf_id2int(hw, BCF_DT_SAMPLE, name)!=-1 )
        {
            // there is a sample with the same name
            free(name);
            int len = strlen(hr->samples[i]) + strlen(clash_prefix) + 1;
            name = (char*) malloc(sizeof(char)*(len+1));
            sprintf(name,"%s:%s",clash_prefix,hr->samples[i]);
        }
        bcf_hdr_add_sample(hw,name);
    }
}

void debug_als(char **als, int nals)
{
    int k; for (k=0; k<nals; k++) printf("%s ", als[k]); 
    printf("\n");
}

/**
 * normalize_alleles() - create smallest possible representation of the alleles
 * @als:    alleles to be merged, first is REF (rw)
 * @nals:   number of $a alleles
 *
 * Best explained on an example:
 *      In:  REF=GTTT  ALT=GTT
 *      Out: REF=GT    ALT=G
 *
 * Note: the als array will be modified
 */
void normalize_alleles(char **als, int nals)
{
    int j, i = 1, done = 0, rlen = strlen(als[0]);
    while ( i<rlen )
    {
        for (j=1; j<nals; j++)
        {
            int len = strlen(als[j]);
            if ( i>=len ) done = 1;
            if ( als[j][len-i] != als[0][rlen-i] ) { done = 1; break; }
        }
        if ( done ) break;
        i++;
    }
    if ( i>1 )
    {
        i--;
        als[0][rlen-i] = 0;
        for (j=1; j<nals; j++) als[j][strlen(als[j])-i] = 0;
    }
}

 /**
 * merge_alleles() - merge two REF,ALT records, $a and $b into $b.
 * @a:      alleles to be merged, first is REF
 * @na:     number of $a alleles
 * @map:    map from the original $a indexes to new $b indexes (0-based)
 * @b:      alleles to be merged, the array will be expanded as required
 * @nb:     number of $b alleles
 * @mb:     size of $b
 *
 * Returns $b expanded to incorporate $a alleles and sets $map. Best explained 
 * on an example:
 *      In:     REF   ALT
 *           a: ACG,  AC,A    (1bp and 2bp deletion)
 *           b: ACGT, A       (3bp deletion)
 *      Out:
 *           b: ACGT, A,ACT,AT (3bp, 1bp and 2bp deletion)
 *           map: 0,2,3
 * Here the mapping from the original $a alleles to the new $b alleles is 0->0, 
 * 1->2, and 2->3.
 */
char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb)
{
    // reference allele never changes
    map[0] = 0;

    int rla = !a[0][1] ? 1 : strlen(a[0]);
    int rlb = !b[0][1] ? 1 : strlen(b[0]);

    // the most common case: same SNPs
    if ( na==2 && *nb==2 && rla==1 && rlb==1 && a[1][0]==b[1][0] && !a[1][1] && !b[1][1] )
    {
        map[1] = 1;
        return b;
    }

    // Sanity check: reference prefixes must be identical
    if ( strncmp(a[0],b[0],rla<rlb?rla:rlb) ) error("The REF prefixes differ: %s vs %s (%d,%d)\n", a[0],b[0],rla,rlb);

    int n = *nb + na;
    hts_expand0(char*,n,*mb,b);

    // $b alleles need expanding
    int i,j;
    if ( rla>rlb )
    {
        for (i=0; i<*nb; i++)
        {
            int l = strlen(b[i]);
            b[i] = (char*) realloc(b[i],l+rla-rlb+1);
            memcpy(b[i]+l,a[0]+rlb,rla-rlb+1);
        }
    }

    // now check if the $a alleles are present and if not add them
    for (i=1; i<na; i++)
    {
        char *ai;
        if ( rlb>rla )  // $a alleles need expanding
        {
            int l = strlen(a[i]);
            ai = (char*) malloc(l+rlb-rla+1);
            memcpy(ai,a[i],l);
            memcpy(ai+l,b[0]+rla,rlb-rla+1);
        }
        else
            ai = a[i];

        for (j=1; j<*nb; j++)
            if ( !strcmp(ai,b[j]) ) break;

        if ( j<*nb ) // $b already has the same allele
        {
            map[i] = j;
            if ( rlb>rla ) free(ai);
            continue; 
        }
        // new allele
        map[i] = *nb;
        b[*nb] = rlb>rla ? ai : strdup(ai);
        (*nb)++;
    }
    return b;
}

maux_t *maux_init(bcf_srs_t *files)
{
    maux_t *ma = (maux_t*) calloc(1,sizeof(maux_t));
    ma->n      = files->nreaders;
    ma->nbuf   = (int *) calloc(ma->n,sizeof(int));
    ma->d      = (maux1_t**) calloc(ma->n,sizeof(maux1_t*));
    ma->files  = files;
    int i, n_smpl = 0;
    for (i=0; i<ma->n; i++)
        n_smpl += files->readers[i].header->n[BCF_DT_SAMPLE];
    ma->smpl_ploidy = (int*) calloc(n_smpl,sizeof(int));
    ma->smpl_nGsize = (int*) malloc(n_smpl*sizeof(int));
    ma->has_line = (int*) malloc(ma->n*sizeof(int));
    return ma;
}
void maux_destroy(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) // for each reader
    {
        if ( !ma->d[i] ) continue;
        int j;
        for (j=0; j<ma->nbuf[i]; j++)  // for each buffered line
            if ( ma->d[i][j].map ) free(ma->d[i][j].map);
        free(ma->d[i]);
    }
    if (ma->ntmp_arr) free(ma->tmp_arr);
    if (ma->nfmt_map) free(ma->fmt_map);
    // ma->inf freed in bcf_destroy1
    free(ma->d);
    free(ma->nbuf);
    for (i=0; i<ma->mals; i++) free(ma->als[i]);
    free(ma->als);
    free(ma->cnt);
    free(ma->smpl_ploidy);
    free(ma->smpl_nGsize);
    free(ma->has_line);
    free(ma);
}
void maux_expand1(maux_t *ma, int i)
{
    if ( ma->nbuf[i] <= ma->files->readers[i].nbuffer )
    {
        int n = ma->files->readers[i].nbuffer + 1;
        ma->d[i] = (maux1_t*) realloc(ma->d[i], sizeof(maux1_t)*n);
        memset(ma->d[i]+ma->nbuf[i],0,sizeof(maux1_t)*(n-ma->nbuf[i]));
        ma->nbuf[i] = n;
    }
}
void maux_reset(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) maux_expand1(ma, i);
    for (i=1; i<ma->ncnt; i++) ma->cnt[i] = 0;
}
void maux_debug(maux_t *ma, int ir, int ib)
{
    printf("[%d,%d]\t", ir,ib);
    int i;
    for (i=0; i<ma->nals; i++)
    {
        printf(" %s [%d]", ma->als[i], ma->cnt[i]);
    }
    printf("\n");
}

void merge_chrom2qual(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    kstring_t *tmps = &args->tmps;
    tmps->l = 0;

    maux_t *ma = args->maux;
    int *al_idxs = (int*) calloc(ma->nals,sizeof(int));
    out->qual = 0;

    // CHROM, POS, ID, QUAL
    out->pos = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;

        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;

        // alleles
        int j;
        for (j=1; j<line->n_allele; j++) 
            al_idxs[ ma->d[i][0].map[j] ] = 1;

        // position
        if ( out->pos==-1 )
        {
            const char *chr = hdr->id[BCF_DT_CTG][line->rid].key;
            out->rid = bcf_name2id(out_hdr, chr);
            if ( strcmp(chr,out_hdr->id[BCF_DT_CTG][out->rid].key) ) error("Uh\n"); 
            out->pos = line->pos;
        }

        // ID
        if ( line->d.id[0]!='.' || line->d.id[1] )
        {
            kitr = kh_get(strdict, tmph, line->d.id);
            if ( kitr == kh_end(tmph) )
            {
                if ( tmps->l ) kputc(';', tmps);
                kputs(line->d.id, tmps);
                kh_put(strdict, tmph, line->d.id, &ret);
            }
        }

        // set QUAL to the max qual value. Not exactly correct, but good enough for now
        if ( out->qual < files->readers[i].buffer[0]->qual ) out->qual = files->readers[i].buffer[0]->qual;
    }

    // set ID
    if ( !tmps->l ) kputs(".", tmps);
    if ( out->d.id ) free(out->d.id);
    out->d.id = strdup(tmps->s);

    // set alleles
    out->n_allele = 0;
    for (i=1; i<ma->nals; i++) 
    {
        if ( !al_idxs[i] ) continue;
        out->n_allele++;

        // Adjust the indexes, the allele map could be created for multiple collapsed records, 
        //  some of which might be unused for this output line
        int ir, j;
        for (ir=0; ir<files->nreaders; ir++)
        {
            if ( !ma->has_line[ir] ) continue;
            bcf1_t *line = files->readers[ir].buffer[0];
            for (j=1; j<line->n_allele; j++)
                if ( ma->d[ir][0].map[j]==i ) ma->d[ir][0].map[j] = out->n_allele;
        }
    }
    // Expand the arrays and realloc the alleles string. Note that all alleles are in a single allocated block.
    out->n_allele++;
    hts_expand(char*, out->n_allele, out->d.m_allele, out->d.allele);
    int k = 0;  // new size
    for (i=0; i<ma->nals; i++)
        if ( i==0 || al_idxs[i] ) out->d.allele[k++] = ma->als[i];
    normalize_alleles(out->d.allele, out->n_allele);
    assert( k==out->n_allele );

    k = 0;
    for (i=0; i<out->n_allele; i++) k += strlen(out->d.allele[i]) + 1;
    hts_expand(char, k, out->d.m_als, out->d.als);

    // Copy the alleles 
    char *dst = out->d.als;
    for (i=0; i<out->n_allele; i++)
    {
        char *src = out->d.allele[i];
        out->d.allele[i] = dst;
        while ( *src ) { *dst = *src; dst++; src++; } 
        *dst = 0;
        dst++;
    }
    free(al_idxs);
}

void merge_filter(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->d.n_flt = 0;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i]) continue;

        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        bcf_unpack(line, BCF_UN_ALL);

        int k;
        for (k=0; k<line->d.n_flt; k++)
        {
            const char *flt = hdr->id[BCF_DT_ID][line->d.flt[k]].key;
            kitr = kh_get(strdict, tmph, flt);
            if ( kitr == kh_end(tmph) )
            {
                int id = bcf_id2int(out_hdr, BCF_DT_ID, flt);
                if ( id==-1 ) error("The filter not defined: %s\n", flt);
                hts_expand(int,out->d.n_flt+1,ma->mflt,ma->flt);
                ma->flt[out->d.n_flt] = id;
                out->d.n_flt++;
                kh_put(strdict, tmph, flt, &ret);
            }
        }
    }
    // Check if PASS is not mixed with other filters
    if ( out->d.n_flt>1 )
    {
        int id = bcf_id2int(out_hdr, BCF_DT_ID, "PASS");
        for (i=0; i<out->d.n_flt; i++)
            if ( ma->flt[i]==id ) break;
        if ( i<out->d.n_flt )
        {
            out->d.n_flt--;
            for (; i<out->d.n_flt; i++) ma->flt[i] = ma->flt[i+1];
        }
    }
    out->d.flt = ma->flt;
}

void merge_info(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, j, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->n_info = 0;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_info; j++) 
        {
            bcf_info_t *inf = &line->d.info[j];

            const char *key = hdr->id[BCF_DT_ID][inf->key].key;
            kitr = kh_get(strdict, tmph, key);
            if ( kitr == kh_end(tmph) )
            {
                int id = bcf_id2int(out_hdr, BCF_DT_ID, key);
                if ( id==-1 ) error("Error: The INFO field not defined: %s\n", key);
                hts_expand(bcf_info_t,out->n_info+1,ma->minf,ma->inf);
                ma->inf[out->n_info].key  = id;
                ma->inf[out->n_info].type = inf->type;
                ma->inf[out->n_info].len  = inf->len;
                ma->inf[out->n_info].vptr = inf->vptr;
                ma->inf[out->n_info].v1.i = inf->v1.i;
                ma->inf[out->n_info].v1.f = inf->v1.f;
                ma->inf[out->n_info].vptr_off  = inf->vptr_off;
                ma->inf[out->n_info].vptr_len  = inf->vptr_len;
                ma->inf[out->n_info].vptr_free = inf->vptr_free;
                out->n_info++;
                kh_put(strdict, tmph, key, &ret);
            }
            // todo: G-tags, A-tags
        }
    }
    out->d.info = ma->inf;
    out->d.m_info = ma->minf;
    for (i=out->n_info; i<out->d.m_info; i++) out->d.info[i].vptr_free = 0;
}

// Only existing AN, AC will be modified. If not present, the line stays unchanged
void update_AN_AC(bcf_hdr_t *hdr, bcf1_t *line)
{
    int i;
    int AN_id = bcf_id2int(hdr, BCF_DT_ID, "AN");
    int AC_id = bcf_id2int(hdr, BCF_DT_ID, "AC");
    if ( AN_id<0 && AC_id<0 ) return;

    bcf_info_t *AN_ptr = NULL, *AC_ptr = NULL;
    if ( AN_id>=0 )
    {
        for (i=0; i<line->n_info; i++)
            if ( AN_id==line->d.info[i].key ) 
            {
                AN_ptr = &line->d.info[i];
                break;
            }
    }
    if ( AC_id>=0 )
    {
        for (i=0; i<line->n_info; i++)
            if ( AC_id==line->d.info[i].key ) 
            {
                AC_ptr = &line->d.info[i];
                break;
            }
    }
    if ( !AN_ptr && !AC_ptr ) return;

    int32_t an = 0, *tmp = (int32_t*) malloc(sizeof(int)*line->n_allele);
    int ret = bcf_calc_ac(hdr, line, tmp, BCF_UN_FMT);
    if ( ret>0 )
    {
        for (i=0; i<line->n_allele; i++) an += tmp[i];
        if ( AN_ptr ) bcf1_update_info_int32(hdr, line, "AN", &an, 1);
        if ( AC_ptr ) bcf1_update_info_int32(hdr, line, "AC", tmp+1, line->n_allele-1);
    }
    free(tmp);
}

void merge_GT(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = out_hdr->n[BCF_DT_SAMPLE];
    
    int nsize = 0, msize = sizeof(int32_t);
    for (i=0; i<files->nreaders; i++)
    {
        if ( !fmt_map[i] ) continue;
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }
    memset(ma->smpl_ploidy,0,nsamples*sizeof(int));

    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        int32_t *tmp  = (int32_t *) ma->tmp_arr + ismpl*nsize;

        int j, k;
        if ( !fmt_ori )
        {
            // missing values: assume maximum ploidy
            for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++)
            {
                for (k=0; k<nsize; k++) { tmp[k] = 0; ma->smpl_ploidy[ismpl+j]++; }
                tmp += nsize;
            }
            ismpl += hdr->n[BCF_DT_SAMPLE];
            continue;
        }

        #define BRANCH(type_t, missing, vector_end) { \
            type_t *p_ori  = (type_t*) fmt_ori->p; \
            if ( !ma->d[i][0].als_differ ) \
            { \
                /* the allele numbering is unchanged */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (k=0; k<fmt_ori->n; k++) \
                    { \
                        if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                        ma->smpl_ploidy[ismpl+j]++; \
                        if ( p_ori[k]==missing ) tmp[k] = 0; /* missing allele */ \
                        else tmp[k] = p_ori[k]; \
                    } \
                    for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                    tmp += nsize; \
                    p_ori += fmt_ori->n; \
                } \
                ismpl += hdr->n[BCF_DT_SAMPLE]; \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
            { \
                for (k=0; k<fmt_ori->n; k++) \
                { \
                    if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                    ma->smpl_ploidy[ismpl+j]++; \
                    if ( !(p_ori[k]>>1) || p_ori[k]==missing ) tmp[k] = 0; /* missing allele */ \
                    else \
                    { \
                        int al = (p_ori[k]>>1) - 1; \
                        al = al<=0 ? al + 1 : ma->d[i][0].map[al] + 1; \
                        tmp[k] = (al << 1) | ((p_ori[k])&1); \
                    } \
                } \
                for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                tmp += nsize; \
                p_ori += fmt_ori->n; \
            } \
            ismpl += hdr->n[BCF_DT_SAMPLE]; \
        }
        switch (fmt_ori->type)
        {
            case BCF_BT_INT8: BRANCH(int8_t,   bcf_int8_missing,  bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: error("Unexpected case: %d\n", fmt_ori->type);
        }
        #undef BRANCH
    }
    bcf1_update_format_int32(out_hdr, out, "GT", (int32_t*)ma->tmp_arr, nsamples*nsize);
    for (i=0; i<nsamples; i++)
    {
        assert( ma->smpl_ploidy[i]>0 && ma->smpl_ploidy[i]<=2 );
        ma->smpl_nGsize[i] = ma->smpl_ploidy[i]==1 ?  out->n_allele : out->n_allele*(out->n_allele + 1)/2;
    }
}

void merge_format_field(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = out_hdr->n[BCF_DT_SAMPLE];

    const char *key = NULL;
    int nsize = 0, length = BCF_VL_FIXED, type = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        if ( !fmt_map[i] ) continue;
        if ( !key ) key = files->readers[i].header->id[BCF_DT_ID][fmt_map[i]->id].key;
        type = fmt_map[i]->type;
        if ( IS_VL_G(files->readers[i].header, fmt_map[i]->id) ) { length = BCF_VL_G; nsize = out->n_allele*(out->n_allele + 1)/2; break; }
        if ( IS_VL_A(files->readers[i].header, fmt_map[i]->id) ) { length = BCF_VL_A; nsize = out->n_allele - 1; break; }
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    int msize = sizeof(float)>sizeof(int32_t) ? sizeof(float) : sizeof(int32_t);
    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }

    // Fill the temp array for all samples by collecting values from all files
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        if ( fmt_ori ) type = fmt_ori->type;

        // set the values
        #define BRANCH(tgt_type_t, src_type_t, src_is_missing, src_is_vector_end, tgt_set_missing, tgt_set_vector_end) { \
            int j, l, k; \
            tgt_type_t *tgt = (tgt_type_t *) ma->tmp_arr + ismpl*nsize; \
            if ( !fmt_ori ) \
            { \
                /* the field is not present in this file, set missing values */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    tgt_set_missing; tgt++; for (l=1; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                } \
                ismpl += hdr->n[BCF_DT_SAMPLE]; \
                continue; \
            } \
            assert( ma->has_line[i] ); \
            /* if ( !ma->has_line[i] ) continue; */ \
            bcf1_t *line    = reader->buffer[0]; \
            src_type_t *src = (src_type_t*) fmt_ori->p; \
            if ( (length!=BCF_VL_G && length!=BCF_VL_A) || (line->n_allele==out->n_allele && !ma->d[i][0].als_differ) ) \
            { \
                /* alleles unchanged, copy over */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (l=0; l<fmt_ori->n; l++) \
                    { \
                        if ( src_is_vector_end ) break; \
                        else if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        tgt++; src++; \
                    } \
                    for (k=l; k<nsize; k++) { tgt_set_vector_end; tgt++; } \
                    for (k=l; k<fmt_ori->n; k++) { src++; } \
                } \
                ismpl += hdr->n[BCF_DT_SAMPLE]; \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            if ( length==BCF_VL_G ) \
            { \
                /* Number=G tags */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    for (l=0; l<ma->smpl_nGsize[ismpl+j]; l++) { tgt_set_missing; tgt++; } \
                    for (; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                    int iori,jori, inew,jnew; \
                    for (iori=0; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori]; \
                        for (jori=0; jori<=iori; jori++) \
                        { \
                            jnew = ma->d[i][0].map[jori]; \
                            int kori = iori*(iori+1)/2 + jori; \
                            int knew = inew*(inew+1)/2 + jnew; \
                            src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + kori; \
                            tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + knew; \
                            if ( src_is_vector_end ) \
                            { \
                                iori = line->n_allele; \
                                break; \
                            } \
                            if ( src_is_missing ) tgt_set_missing; \
                            else *tgt = *src; \
                        } \
                    } \
                } \
            } \
            else \
            { \
                /* Number=A tags */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    for (l=0; l<nsize; l++) { tgt_set_missing; tgt++; } \
                    int iori,inew; \
                    for (iori=1; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori] - 1; \
                        tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                        if ( src_is_vector_end ) break; \
                        if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        src++; \
                    } \
                } \
            } \
            ismpl += hdr->n[BCF_DT_SAMPLE]; \
        }
        switch (type)
        {
            case BCF_BT_INT8:  BRANCH(int32_t,  int8_t, *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT16: BRANCH(int32_t, int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_FLOAT: BRANCH(float, float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*tgt), bcf_float_set_vector_end(*tgt)); break;
            default: error("Unexpected case: %d, %s\n", type, key);
        }
        #undef BRANCH
    }
    if ( type==BCF_BT_FLOAT )
        bcf1_update_format_float(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else
        bcf1_update_format_int32(out_hdr, out, key, (int32_t*)ma->tmp_arr, nsamples*nsize);
}

void merge_format(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    if ( !ma->nfmt_map ) 
    {
        ma->nfmt_map = 2;
        ma->fmt_map  = (bcf_fmt_t**) calloc(ma->nfmt_map*files->nreaders, sizeof(bcf_fmt_t*));
    }
    else
        memset(ma->fmt_map, 0, ma->nfmt_map*files->nreaders*sizeof(bcf_fmt_t**));

    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    int i, j, ret, has_GT = 0, max_ifmt = 0; // max fmt index
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_fmt; j++) 
        {
            // Wat this tag already seen?
            bcf_fmt_t *fmt = &line->d.fmt[j];
            const char *key = hdr->id[BCF_DT_ID][fmt->id].key;
            kitr = kh_get(strdict, tmph, key);

            int ifmt;
            if ( kitr != kh_end(tmph) )
                ifmt = kh_value(tmph, kitr);    // seen
            else
            {
                // new FORMAT tag
                if ( key[0]=='G' && key[1]=='T' && key[2]==0 ) { has_GT = 1; ifmt = 0; }
                else 
                {
                    ifmt = ++max_ifmt;
                    if ( max_ifmt >= ma->nfmt_map )
                    {
                        ma->fmt_map = (bcf_fmt_t**) realloc(ma->fmt_map, sizeof(bcf_fmt_t*)*(max_ifmt+1)*files->nreaders);
                        memset(ma->fmt_map+ma->nfmt_map*files->nreaders, 0, (max_ifmt-ma->nfmt_map+1)*files->nreaders*sizeof(bcf_fmt_t*));
                        ma->nfmt_map = max_ifmt+1;
                    }
                }
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_value(tmph, kitr) = ifmt;
            }
            ma->fmt_map[ifmt*files->nreaders+i] = fmt;
        }
        // Check if the allele numbering must be changed
        for (j=1; j<reader->buffer[0]->n_allele; j++)
            if ( ma->d[i][0].map[j]!=j ) break;
        ma->d[i][0].als_differ = j==reader->buffer[0]->n_allele ? 0 : 1;
    }

    out->n_sample = out_hdr->n[BCF_DT_SAMPLE];
    if ( has_GT )
        merge_GT(args, ma->fmt_map, out);
    update_AN_AC(out_hdr, out);

    for (i=1; i<=max_ifmt; i++)
        merge_format_field(args, &ma->fmt_map[i*files->nreaders], out);
    out->d.indiv_dirty = 1;
}

// The core merging function, one or none line from each reader
void merge_line(args_t *args)
{
    bcf1_t *out = args->out_line;
    bcf_clear1(out);
    out->unpacked = BCF_UN_ALL;

    merge_chrom2qual(args, out);
    merge_filter(args, out);
    merge_info(args, out);
    merge_format(args, out);

    vcf_write1(args->out_fh, args->out_hdr, out);
}


void debug_buffers(FILE *fp, bcf_srs_t *files);
void debug_buffer(FILE *fp, bcf_sr_t *reader);

#define SWAP(type_t,a,b) { type_t tmp = (a); (a) = (b); (b) = tmp; }

// Clean the reader's buffer to and make it ready for the next next_line() call.
// Moves finished records (SKIP_DONE flag set) at the end of the buffer and put
// the rest to the beggining. Then shorten the buffer so that the last element
// points to the last unfinished record. There are two special cases: the last
// line of the buffer typically has a different position and must stay at the
// end; next, the first record of the buffer must be one of those already
// printed, as it will be discarded by next_line().
//
void shake_buffer(maux_t *maux, int ir, int pos)
{
    bcf_sr_t *reader = &maux->files->readers[ir];
    maux1_t *m = maux->d[ir];

    if ( !reader->buffer ) return;

    int i;
    // FILE *fp = stdout;
    // fprintf(fp,"<going to shake> nbuf=%d\t", reader->nbuffer); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]); fprintf(fp,"\n");
    // debug_buffer(fp,reader);
    // fprintf(fp,"--\n");

    int a = 1, b = reader->nbuffer;
    if ( reader->buffer[b]->pos != pos ) b--;   // move the last line separately afterwards

    while ( a<b )
    {
        if ( !(m[a].skip&SKIP_DONE) ) { a++; continue; }
        if ( m[b].skip&SKIP_DONE ) { b--; continue; }
        SWAP(bcf1_t*, reader->buffer[a], reader->buffer[b]);
        SWAP(maux1_t, m[a], m[b]);
        a++;
        b--;
    }

    // position $a to the after the first unfinished record 
    while ( a<=reader->nbuffer && !(m[a].skip&SKIP_DONE) ) a++;

    if ( a<reader->nbuffer )
    {
        // there is a gap between the unfinished lines at the beggining and the
        // last line. The last line must be brought forward to fill the gap
        if ( reader->buffer[reader->nbuffer]->pos != pos ) 
        {
            SWAP(bcf1_t*, reader->buffer[a], reader->buffer[reader->nbuffer]);
            SWAP(maux1_t, m[a], m[reader->nbuffer]);
            reader->nbuffer = a;
        }
    }

    if ( !(m[0].skip&SKIP_DONE) && reader->buffer[0]->pos==pos )
    {
        // the first record is unfinished, replace it with an empty line
        // from the end of the buffer or else next_line will remove it
        if ( reader->nbuffer + 1 >= maux->nbuf[ir] ) 
        {
            reader->nbuffer++;
            maux_expand1(maux, ir);
            reader->nbuffer--;
            m = maux->d[ir];
        }
        if ( reader->nbuffer+1 >= reader->mbuffer ) 
            error("Uh, did not expect this: %d vs %d\n", reader->nbuffer,reader->mbuffer);

        if ( reader->buffer[reader->nbuffer]->pos!=pos ) 
        {
            // 4way swap
            bcf1_t *tmp = reader->buffer[0];
            reader->buffer[0] = reader->buffer[reader->nbuffer+1];
            reader->buffer[reader->nbuffer+1] = reader->buffer[reader->nbuffer];
            reader->buffer[reader->nbuffer] = tmp;
            m[reader->nbuffer].skip   = m[0].skip;
            m[reader->nbuffer+1].skip = SKIP_DIFF;
            reader->nbuffer++;
        }
        else
        {
            SWAP(bcf1_t*, reader->buffer[0], reader->buffer[reader->nbuffer+1]);
            SWAP(maux1_t, m[0], m[reader->nbuffer+1]);
        }
    }

    // debug_buffer(fp,reader);
    // fprintf(fp,"<shaken>\t"); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]);
    // fprintf(fp,"\n\n");

    // set position of finished buffer[0] line to -1, otherwise swapping may
    // bring it back after next_line()
    reader->buffer[0]->pos = -1;

    // trim the buffer, remove finished lines from the end
    i = reader->nbuffer;
    while ( i>=1 && m[i--].skip&SKIP_DONE )
        reader->nbuffer--;
}

// Determine which line should be merged from which reader: go through all
// readers and all buffered lines, expand REF,ALT and try to match lines with
// the same ALTs. A step towards output independent on input ordering of the
// lines.
void merge_buffer(args_t *args)
{
    bcf_srs_t *files = args->files;
    int i, pos = -1, var_type = 0;
    maux_t *maux = args->maux;
    maux_reset(maux);

    // set the current position
    for (i=0; i<files->nreaders; i++)
    {
        if ( bcf_sr_has_line(files,i) )
        {
            pos = files->readers[i].buffer[0]->pos;
            bcf_set_variant_types(files->readers[i].buffer[0]);
            var_type = files->readers[i].buffer[0]->d.var_type;
            break;
        }
    }

    // go through all files and all lines at this position and normalize
    // relevant alleles
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        if ( !reader->buffer ) continue;
        int j;
        for (j=0; j<=reader->nbuffer; j++)
        {
            bcf1_t *line = reader->buffer[j];
            bcf_set_variant_types(line);

            // select relevant lines
            maux->d[i][j].skip = SKIP_DIFF;
            if ( pos!=line->pos ) 
            {
                if ( j==0 ) maux->d[i][j].skip |= SKIP_DONE; // left from previous run, force to ignore
                continue; 
            }
            if ( args->collapse==COLLAPSE_NONE && var_type!=line->d.var_type ) continue;
            if ( var_type&VCF_SNP && !(line->d.var_type&VCF_SNP) && !(args->collapse&COLLAPSE_ANY) ) continue;
            if ( var_type&VCF_INDEL && !(line->d.var_type&VCF_INDEL) && !(args->collapse&COLLAPSE_ANY) ) continue;
            maux->d[i][j].skip = 0;

            hts_expand(int, line->n_allele, maux->d[i][j].mmap, maux->d[i][j].map);
            int k;
            if ( !maux->nals )    // first record, copy the alleles to the output
            {
                maux->nals = line->n_allele;
                hts_expand0(char*, maux->nals, maux->mals, maux->als);
                hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
                for (k=0; k<maux->nals; k++)
                {
                    maux->als[k] = strdup(line->d.allele[k]);
                    maux->d[i][j].map[k] = k;
                    maux->cnt[k] = 1;
                }
                pos = line->pos;
                continue;
            }

            // normalize alleles
            maux->als = merge_alleles(line->d.allele, line->n_allele, maux->d[i][j].map, maux->als, &maux->nals, &maux->mals);
            hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
            for (k=1; k<line->n_allele; k++)
                maux->cnt[ maux->d[i][j].map[k] ]++;
        }
    }

    // Select records that have the same alleles; the input ordering of indels
    // must not matter. Multiple VCF lines can be emitted from this loop.
    // We expect only very few alleles and not many records with the same
    // position in the buffers, therefore the nested loops should not slow us
    // much.
    while (1)
    {
        // take the most frequent allele present in multiple files
        int icnt = 1;
        for (i=1; i<maux->nals; i++) 
            if ( maux->cnt[i] > maux->cnt[icnt] ) icnt = i;
        if ( maux->cnt[icnt]<0 ) break;

        int nmask = 0;
        for (i=0; i<files->nreaders; i++)
        {
            maux->has_line[i] = 0;

            // first pass: try to find lines with the same allele
            bcf_sr_t *reader = &files->readers[i];
            if ( !reader->buffer ) continue;
            int j;
            for (j=0; j<=reader->nbuffer; j++)
            {
                if ( maux->d[i][j].skip ) continue;
                int k;
                for (k=1; k<reader->buffer[j]->n_allele; k++)
                    if ( icnt==maux->d[i][j].map[k] ) break;
                if ( k<reader->buffer[j]->n_allele ) break;
            }
            if ( j>reader->nbuffer )
            {
                // no matching allele found in this file, any allele must do
                for (j=0; j<=reader->nbuffer; j++)
                    if ( !maux->d[i][j].skip ) break;
            }
            if ( j<=reader->nbuffer ) 
            {
                // found a suitable line for merging, place it at the beggining
                if ( j>0 ) 
                {
                    SWAP(bcf1_t*, reader->buffer[0], reader->buffer[j]); 
                    SWAP(maux1_t, maux->d[i][0], maux->d[i][j]); 
                }
                // mark as finished so that it's ignored next time
                maux->d[i][0].skip |= SKIP_DONE;
                maux->has_line[i] = 1;
                nmask++;
            }
        }
        if ( !nmask ) break;    // done, no more lines suitable for merging found 
        merge_line(args);       // merge and output the line
        maux->cnt[icnt] = -1;   // do not pick this allele again, mark it as finished
    }

    // clean the alleles
    for (i=0; i<maux->nals; i++)
    {
        free(maux->als[i]);
        maux->als[i] = 0;
    }
    maux->nals = 0;

    // get the buffers ready for the next next_line() call
    for (i=0; i<files->nreaders; i++)
        shake_buffer(maux, i, pos);
}

void merge_vcf(args_t *args)
{
    args->out_fh = hts_open("-","w",0);
    args->out_hdr = bcf_hdr_init();

    if ( args->header_fname )
    {
        if ( bcf_hdr_set(args->out_hdr,args->header_fname) ) error("Could not read/parse the header: %s\n", args->header_fname);
    }
    else
    {
        int i;
        for (i=0; i<args->files->nreaders; i++)
        {
            char buf[10]; snprintf(buf,10,"%d",i+1);
            bcf_hdr_merge(args->out_hdr, args->files->readers[i].header,buf);
        }
        kstring_t str = {0,0,0};
        ksprintf(&str,"##vcfmergeVersion=%s\n", HTS_VERSION);
        bcf_hdr_append(args->out_hdr,str.s);

        str.l = 0;
        ksprintf(&str,"##vcfmergeCommand=%s", args->argv[0]);
        for (i=1; i<args->argc; i++) ksprintf(&str, " %s", args->argv[i]);
        kputc('\n', &str);
        bcf_hdr_append(args->out_hdr,str.s);
        free(str.s);

        bcf_hdr_sync(args->out_hdr);
        bcf_hdr_fmt_text(args->out_hdr);
    }

    vcf_hdr_write(args->out_fh, args->out_hdr);
    if ( args->header_only )
    {
        bcf_hdr_destroy(args->out_hdr);
        hts_close(args->out_fh);
        return;
    }

    args->maux = maux_init(args->files);
    args->out_line = bcf_init1();
    args->tmph = kh_init(strdict);
    int ret;
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        merge_buffer(args);
    }
    maux_destroy(args->maux);
    bcf_hdr_destroy(args->out_hdr);
    hts_close(args->out_fh);
    bcf_destroy1(args->out_line);
    kh_destroy(strdict, args->tmph);
    if ( args->tmps.m ) free(args->tmps.s);
}

static void usage(void)
{
    fprintf(stderr, "Usage:   vcfmerge [options] <A.vcf.gz> <B.vcf.gz> ...\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        --use-header <file>           use the provided header\n");
    fprintf(stderr, "        --print-header <file>         print only header of the output file and exit\n");
    fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
    fprintf(stderr, "    -m, --merge <string>              merge sites with differing alleles for <snps|indels|both|any>\n");
    fprintf(stderr, "    -r, --region <chr|chr:from-to>    merge in the given region only\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfmerge(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"merge",1,0,'m'},
        {"apply-filters",0,0,'f'},
        {"use-header",1,0,1},
        {"print-header",0,0,2},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hm:fr:1:2",loptions,NULL)) >= 0) {
        switch (c) {
            case 'm':
                if ( !strcmp(optarg,"snps") ) args->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) args->collapse |= COLLAPSE_ANY;
                break;
            case 'f': args->files->apply_filters = 1; break;
            case 'r': args->files->region = optarg; break;
            case  1 : args->header_fname = optarg; break;
            case  2 : args->header_only = 1; break;
            case 'h': 
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if (argc == optind) usage();
    if ( argc-optind<2 ) usage();

    args->files->require_index = 1;

    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open: %s\n", argv[optind]);
        optind++;
    }
    merge_vcf(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

