#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

#include "khash.h"
KHASH_MAP_INIT_STR(strdict, const char *)
typedef khash_t(strdict) strdict_t;

extern uint32_t bcf_missing_float;
#define SET_FLOAT_MISSING(var) *((uint32_t*)&(var)) = (uint32_t) bcf_missing_float;

#define SKIP_DONE 1
#define SKIP_DIFF 2

#define IS_VL_G(hdr,_id) (((hdr)->id[BCF_DT_ID][_id].val->info[BCF_HL_FMT]>>8 & 0xf) == BCF_VL_G)
#define IS_VL_A(hdr,_id) (((hdr)->id[BCF_DT_ID][_id].val->info[BCF_HL_FMT]>>8 & 0xf) == BCF_VL_A)

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
    float *weight;  // weighting factor for recalculating qualities
    int *flt, mflt, minf;
    bcf_info_t *inf;// out_line's INFO fields
    bcf_fmt_t *fmt; // out_line's indiv fields
    int mfmt, *fmt_map, mfmt_map;
    maux1_t **d;    // d[i][j] i-th reader, j-th buffer line
    readers_t *files;
}
maux_t;

typedef struct
{
    maux_t *maux;
    int header_only, collapse;
    char *header_fname;
    strdict_t *tmph;
    kstring_t tmps;
    readers_t *files;
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
        if ( hr->hrec[i]->type==BCF_HL_GEN )
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
    if ( strncmp(a[0],b[0],rla<rlb?rla:rlb) ) error("The REF prefixes differ: %s vs %s (%d,%d)\n", a[0],b[0]);

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

maux_t *maux_init(readers_t *files)
{
    maux_t *ma = (maux_t*) calloc(1,sizeof(maux_t));
    ma->n      = files->nreaders;
    ma->nbuf   = (int *) calloc(ma->n,sizeof(int));
    ma->d      = (maux1_t**) calloc(ma->n,sizeof(maux1_t*));
    ma->files  = files;
    ma->weight = (float*) calloc(ma->n,sizeof(float));
    int i, n = 0;
    for (i=0; i<ma->n; i++)
        n += files->readers[i].header->n[BCF_DT_SAMPLE];
    for (i=0; i<ma->n; i++)
        ma->weight[i] = n ? (float)files->readers[i].header->n[BCF_DT_SAMPLE]/n : (float)1/ma->n;
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
    if (ma->mfmt_map) free(ma->fmt_map);
    for (i=0; i<ma->mfmt; i++)
        if ( ma->fmt[i].p ) free(ma->fmt[i].p);
    // ma->fmt freed in bcf_destroy1
    // ma->inf freed in bcf_destroy1
    free(ma->d);
    free(ma->nbuf);
    for (i=0; i<ma->mals; i++) free(ma->als[i]);
    free(ma->als);
    free(ma->cnt);
    free(ma->weight);
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

void merge_chrom2qual(args_t *args, int mask, bcf1_t *out)
{
    readers_t *files = args->files;
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
        if ( !(mask&1<<i) ) continue;

        reader_t *reader = &files->readers[i];
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
    int k = 0;
    for (i=1; i<ma->nals; i++) 
    {
        if ( !al_idxs[i] ) continue;
        out->n_allele++;

        // Adjust the indexes, the allele map could be created for multiple collapsed records, 
        //  some of which might be unused for this output line
        int ir, j;
        for (ir=0; ir<files->nreaders; ir++)
        {
            bcf1_t *line = files->readers[ir].buffer[0];
            for (j=1; j<line->n_allele; j++)
                if ( ma->d[ir][0].map[j]==i ) ma->d[ir][0].map[j] = out->n_allele;
        }
    }
    out->n_allele++;
    out->d.allele = (char **) malloc(sizeof(char*)*out->n_allele);
    for (i=0; i<ma->nals; i++)
        if ( i==0 || al_idxs[i] ) out->d.allele[k++] = strdup(ma->als[i]);
    normalize_alleles(out->d.allele, out->n_allele);
    free(al_idxs);
}

void merge_filter(args_t *args, int mask, bcf1_t *out)
{
    readers_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->d.n_flt = 0;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !(mask&1<<i) ) continue;

        reader_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;

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

void merge_info(args_t *args, int mask, bcf1_t *out)
{
    readers_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, j, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->n_info = 0;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !(mask&1<<i) ) continue;
        reader_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        bcf_unpack(line, BCF_UN_ALL);
        for (j=0; j<line->n_info; j++) 
        {
            bcf_info_t *inf = &line->d.info[j];

            const char *key = hdr->id[BCF_DT_ID][inf->key].key;
            kitr = kh_get(strdict, tmph, key);
            if ( kitr == kh_end(tmph) )
            {
                int id = bcf_id2int(out_hdr, BCF_DT_ID, key);
                if ( id==-1 ) error("The INFO not defined: %s\n", key);
                hts_expand(bcf_info_t,out->n_info+1,ma->minf,ma->inf);
                ma->inf[out->n_info].key  = id;
                ma->inf[out->n_info].type = inf->type;
                ma->inf[out->n_info].len  = inf->len;
                ma->inf[out->n_info].vptr = inf->vptr;
                ma->inf[out->n_info].v1.i = inf->v1.i;
                ma->inf[out->n_info].v1.f = inf->v1.f;
                out->n_info++;
                kh_put(strdict, tmph, key, &ret);
            }
            // todo: AN,AC, G-tags, A-tags
        }
    }
    out->d.info = ma->inf;
}

void merge_GT(args_t *args, int mask, bcf_fmt_t *fmt, int *fmt_map, bcf1_t *out)
{
    readers_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = out_hdr->n[BCF_DT_SAMPLE];
    if ( !fmt->p ) 
        fmt->p = (uint8_t*) malloc(sizeof(uint8_t)*nsamples*fmt->size);
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        bcf1_t *line   = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;

        int j, k = -1;
        if ( mask&1<<i ) k = fmt_map[i] - 1;
        if ( k<0 )
        {
            // missing values
            for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++)
                fmt->p[fmt->size*(ismpl+j)] = (uint8_t)INT8_MIN;
            ismpl += hdr->n[BCF_DT_SAMPLE];
            continue;
        }
        bcf_fmt_t *fmt_ori = &line->d.fmt[k];
        uint8_t *p_out = &fmt->p[fmt->size*ismpl], *p_ori = fmt_ori->p;

        if ( !ma->d[i][0].als_differ )
        {
            // the allele numbering is unchanged
            for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++)
            {
                for (k=0; k<fmt_ori->size; k++, p_out++, p_ori++)
                    *p_out = *p_ori;
                for (; k<fmt->size; k++, p_out++)
                    *p_out =  (uint8_t)INT8_MIN;
            }
            ismpl += hdr->n[BCF_DT_SAMPLE];
            continue;
        }
        // allele numbering needs to be changed
        for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++)
        {
            for (k=0; k<fmt_ori->size; k++, p_out++, p_ori++)
            {
                if ( *p_ori == (uint8_t)INT8_MIN ) break;
                int al = (*p_ori>>1) - 1;
                al = al<=0 ? al + 1 : ma->d[i][0].map[al] + 1;
                *p_out = (al << 1) | ((*p_ori)&1);
            }
            for (; k<fmt->size; k++, p_out++)
                *p_out =  (uint8_t)INT8_MIN;
        }
        ismpl += hdr->n[BCF_DT_SAMPLE];
    }
}

void merge_format_field(args_t *args, int mask, bcf_fmt_t *fmt, int *fmt_map, bcf1_t *out)
{
    readers_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = out_hdr->n[BCF_DT_SAMPLE];
    if ( !fmt->p ) 
        fmt->p = (uint8_t*) malloc(sizeof(uint8_t)*nsamples*fmt->size);
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        bcf1_t *line   = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        int j, k = -1, l;
        if ( mask&1<<i ) k = fmt_map[i] - 1;
        #define BRANCH(type_t, missing) { \
            type_t *p_out = (type_t*) &fmt->p[fmt->size*ismpl]; \
            if ( k<0 ) \
            { \
                /* the field is not present in this file */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (l=0; l<fmt->n; l++) missing; \
                    p_out += fmt->n; \
                } \
                ismpl += hdr->n[BCF_DT_SAMPLE]; \
                continue; \
            } \
            bcf_fmt_t *fmt_ori = &line->d.fmt[k]; \
            type_t *p_ori = (type_t*) fmt_ori->p; \
            if ( (!IS_VL_G(out_hdr,fmt->id) && !IS_VL_A(out_hdr,fmt->id)) || (line->n_allele==out->n_allele && !ma->d[i][0].als_differ) ) \
            { \
                /* alleles unchanged, copy over */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (l=0; l<fmt->n; l++) p_out[l] = p_ori[l]; \
                    p_out = (type_t *)(fmt->size + (void*)p_out); \
                    p_ori = (type_t *)(fmt_ori->size + (void*)p_ori); \
                } \
                ismpl += hdr->n[BCF_DT_SAMPLE]; \
                continue; \
            } \
            /* allele numbering needs to be changed: Number=G tags */ \
            if ( IS_VL_G(out_hdr,fmt->id) ) \
            { \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (l=0; l<fmt->n; l++) missing; \
                    int iori,jori, inew,jnew; \
                    for (iori=0; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori]; \
                        for (jori=0; jori<=iori; jori++) \
                        { \
                            jnew = ma->d[i][0].map[jori]; \
                            int kori = iori*(iori+1)/2 + jori; \
                            int knew = inew*(inew+1)/2 + jnew; \
                            p_out[knew] = p_ori[kori]; \
                        } \
                    } \
                    p_out = (type_t *)(fmt->size + (void*)p_out); \
                    p_ori = (type_t *)(fmt_ori->size + (void*)p_ori); \
                } \
            } \
            else \
            { \
                /* Number=A tags */ \
                for (j=0; j<hdr->n[BCF_DT_SAMPLE]; j++) \
                { \
                    for (l=0; l<fmt->n; l++) missing; \
                    int iori,inew; \
                    for (iori=1; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori] - 1; \
                        p_out[inew] = p_ori[iori - 1]; \
                    } \
                    p_out = (type_t *)(fmt->size + (void*)p_out); \
                    p_ori = (type_t *)(fmt_ori->size + (void*)p_ori); \
                } \
            } \
            ismpl += hdr->n[BCF_DT_SAMPLE]; \
        }
        switch (fmt->type)
        {
            case BCF_BT_INT8: BRANCH(uint8_t, p_out[l] = INT8_MIN); break;
            case BCF_BT_INT16: BRANCH(uint16_t, p_out[l] = INT16_MIN); break;
            case BCF_BT_INT32: BRANCH(uint32_t, p_out[l] = INT32_MIN); break;
            case BCF_BT_FLOAT: BRANCH(float, SET_FLOAT_MISSING(p_out[l]) ); break;
            default: error("Unexpected case: %d\n", fmt->type);
        }
        #undef BRANCH
    }

    #if 0
    // debug
    for (ismpl=0; ismpl<nsamples; ismpl++) 
    {
        uint8_t *p = &fmt->p[ismpl*fmt->size];
        int x=0;
        for (i=0; i<fmt->size; i++)
            if ( !p[i] ) x++; 
    }
    #endif
}

void merge_format(args_t *args, int mask, bcf1_t *out)
{
    readers_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    memset(ma->fmt_map,0,ma->mfmt_map*sizeof(int));
    int nsamples = out_hdr->n[BCF_DT_SAMPLE];
    int i, j, n_fmt = 0;
    int nG_new = out->n_allele*(out->n_allele + 1)/2;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !(mask&1<<i) ) continue;
        reader_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_fmt; j++) 
        {
            bcf_fmt_t *fmt = &line->d.fmt[j];
            int id = bcf_id2int(out_hdr, BCF_DT_ID, hdr->id[BCF_DT_ID][fmt->id].key);
            int k;
            for (k=0; k<n_fmt; k++)
                if ( ma->fmt[k].id ==id ) break;
            int new_n = IS_VL_G(hdr,fmt->id) ? nG_new : (IS_VL_A(hdr,fmt->id) ? out->n_allele - 1 : fmt->n);
            int new_size = new_n * fmt->size / fmt->n;
            if ( k<n_fmt ) // existing format field
            {
                if ( new_size < fmt->size ) 
                {
                    ma->fmt[k].size = new_size;
                    ma->fmt[k].type = fmt->type;
                    ma->fmt[k].p = (uint8_t*) realloc(ma->fmt[k].p, sizeof(uint8_t)*nsamples);
                }
                //printf("already there: %d %s  id=%d  n=%d  size=%d  type=%d\n", k,hdr->id[BCF_DT_ID][fmt->id].key, id,fmt->n,fmt->size,fmt->type);
            }
            else    // new format fields
            {
                n_fmt++;
                hts_expand0(bcf_fmt_t,(n_fmt+1),ma->mfmt,ma->fmt);
                hts_expand0(int,(n_fmt+1)*files->nreaders,ma->mfmt_map,ma->fmt_map);
                ma->fmt[k].id   = id;
                ma->fmt[k].n    = new_n;
                ma->fmt[k].size = new_size;
                ma->fmt[k].type = fmt->type;
                if ( ma->fmt[k].p ) { free(ma->fmt[k].p); ma->fmt[k].p = NULL; }
                // printf("add: %d %s  id=%d  n=%d  size=%d  type=%d info=%d\n", k,hdr->id[BCF_DT_ID][fmt->id].key, id,fmt->n,fmt->size,fmt->type, (hdr->id[BCF_DT_ID][fmt->id].val->info[BCF_HL_FMT]>>8)&0xf);
            }
            ma->fmt_map[k*files->nreaders+i] = j+1;
        }
        // Check if the allele numbering must be changed
        for (j=1; j<reader->buffer[0]->n_allele; j++)
            if ( ma->d[i][0].map[j]!=j ) break;
        ma->d[i][0].als_differ = j==reader->buffer[0]->n_allele ? 0 : 1;
    }
    int gt_id = bcf_id2int(out_hdr,BCF_DT_ID,"GT");
    out->n_fmt = n_fmt;
    for (j=0; j<n_fmt; j++)
    {
        int id = ma->fmt[j].id;
        if ( id==gt_id )
            merge_GT(args, mask, &ma->fmt[j], &ma->fmt_map[j*files->nreaders], out);
        else 
            merge_format_field(args, mask, &ma->fmt[j], &ma->fmt_map[j*files->nreaders], out);
    }
    out->n_sample = nsamples;
    out->d.fmt = ma->fmt;
}

// The core merging function, one or none line from each reader
void merge_line(args_t *args, int mask)
{
    bcf1_t *out = args->out_line;

    out->shared.l = out->indiv.l = 0;
    out->unpacked = BCF_UN_ALL;
    out->unpack_ptr = NULL;

    merge_chrom2qual(args, mask, out);
    merge_filter(args, mask, out);
    merge_info(args, mask, out);
    merge_format(args, mask, out);

    vcf_write1(args->out_fh, args->out_hdr, out);

    int i;
    for (i=0; i<out->n_allele; i++) free(out->d.allele[i]);
    free(out->d.allele);
    out->d.allele = 0;
    out->n_allele = 0;
}


void debug_buffers(FILE *fp, readers_t *files);
void debug_buffer(FILE *fp, reader_t *reader);

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
    reader_t *reader = &maux->files->readers[ir];
    maux1_t *m = maux->d[ir];

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
void merge_buffer(args_t *args, int mask)
{
    readers_t *files = args->files;
    int i, pos = -1, var_type = 0;
    maux_t *maux = args->maux;
    maux_reset(maux);

    // set the current position
    for (i=0; i<files->nreaders; i++)
    {
        if ( mask&1<<i)
        {
            pos = files->readers[i].buffer[0]->pos;
            set_variant_types(files->readers[i].buffer[0]);
            var_type = files->readers[i].buffer[0]->d.var_type;
            break;
        }
    }

    // go through all files and all lines at this position and normalize
    // relevant alleles
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        int j;
        for (j=0; j<=reader->nbuffer; j++)
        {
            bcf1_t *line = reader->buffer[j];
            set_variant_types(line);

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

        mask = 0;
        for (i=0; i<files->nreaders; i++)
        {
            // first pass: try to find lines with the same allele
            reader_t *reader = &files->readers[i];
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
                mask |= 1<<i;
            }
        }
        if ( !mask ) break;     // done, no more lines suitable for merging found 
        merge_line(args,mask);  // merge and output the line
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
        merge_buffer(args, ret);
        // printf("<merge done>\n");
        // debug_buffers(stdout, &args->files);
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

#if 0
int j;
uint32_t i = bcf_missing_float;
float f;
uint8_t *p8;

fprintf(stderr,"sizeof(float) = %ld, sizeof(uint32_t) = %ld\n", sizeof(float),sizeof(uint32_t));
fprintf(stderr,"%d %x %e\n", i,i,(float)i);
p8 = (uint8_t*)&i; for (j=0; j<4; j++) fprintf(stderr," %x", p8[j]); fprintf(stderr,"\n");

*((uint32_t*)&f) = bcf_missing_float;
fprintf(stderr,"%d %x %e\n", (int)f,(int)f,f);
p8 = (uint8_t*)&f; for (j=0; j<4; j++) fprintf(stderr," %x", p8[j]); fprintf(stderr,"\n");

f = bcf_missing_float;
fprintf(stderr,"%d %x %e\n", (int)f,(int)f,f);
p8 = (uint8_t*)&f; for (j=0; j<4; j++) fprintf(stderr," %x", p8[j]); fprintf(stderr,"\n");

return 1;
#endif

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
    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Could not load the index: %s\n", argv[optind]);
        optind++;
    }
    merge_vcf(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

