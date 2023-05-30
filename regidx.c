/*
    Copyright (C) 2014-2019 Genome Research Ltd.

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
#include <strings.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash_str2int.h"
#include "htslib/regidx.h"
#include "hts_internal.h"

#define MAX_COOR_0 REGIDX_MAX   // CSI and hts_itr_query limit, 0-based

#define iBIN(x) ((x)>>13)

typedef struct
{
    hts_pos_t beg, end;
}
reg_t;

typedef struct
{
    hts_pos_t pos;  // position
    uint32_t  ireg; // index to reglist.reg and reglist.dat
}
pos_t;

typedef struct reglist_t reglist_t;

typedef struct
{
    hts_pos_t beg, end; // query region
    uint32_t ireg;      // index of active region
    regidx_t *ridx;
    reglist_t *list;
    int active;
}
itr_t_;

// List of regions for one chromosome.
struct reglist_t
{
    uint32_t *idx, nidx;    // index to list.reg+1
    uint32_t nreg, mreg;    // n:used, m:allocated
    reg_t *reg;             // regions
    uint8_t *dat;           // payload data
    char *seq;              // sequence name
    int unsorted;
};

// Container of all sequences
struct regidx_t
{
    int nseq, mseq;         // n:used, m:alloced
    reglist_t *seq;         // regions for each sequence
    void *seq2regs;         // hash for fast lookup from chr name to regions
    char **seq_names;
    regidx_free_f free;     // function to free any data allocated by regidx_parse_f
    regidx_parse_f parse;   // parse one input line
    void *usr;              // user data to pass to regidx_parse_f
    int payload_size;
    void *payload;          // temporary payload data set by regidx_parse_f (sequence is not known beforehand)
    kstring_t str;
};

int regidx_seq_nregs(regidx_t *idx, const char *seq)
{
    int iseq;
    if ( khash_str2int_get(idx->seq2regs, seq, &iseq)!=0 ) return 0; // no such sequence
    return idx->seq[iseq].nreg;
}

int regidx_nregs(regidx_t *idx)
{
    int i, nreg = 0;
    for (i=0; i<idx->nseq; i++) nreg += idx->seq[i].nreg;
    return nreg;
}

char **regidx_seq_names(regidx_t *idx, int *n)
{
    *n = idx->nseq;
    return idx->seq_names;
}

int regidx_insert_list(regidx_t *idx, char *line, char delim)
{
    kstring_t tmp = KS_INITIALIZE;
    char *ss = line;
    while ( *ss )
    {
        char *se = ss;
        while ( *se && *se!=delim ) se++;
        kputsn(ss, se-ss, ks_clear(&tmp));
        if ( regidx_insert(idx,tmp.s) < 0 )
        {
            ks_free(&tmp);
            return -1;
        }
        if ( !*se ) break;
        ss = se+1;
    }
    ks_free(&tmp);
    return 0;
}

static inline int cmp_regs(reg_t *a, reg_t *b)
{
    if ( a->beg < b->beg ) return -1;
    if ( a->beg > b->beg ) return 1;
    if ( a->end < b->end ) return 1;    // longer intervals come first
    if ( a->end > b->end ) return -1;
    if ( a < b ) return -1; // this is are just for qsort reproducibility across platforms
    if ( a > b ) return 1;
    return 0;
}
static int cmp_reg_ptrs(const void *a, const void *b)
{
    return cmp_regs((reg_t*)a,(reg_t*)b);
}
static int cmp_reg_ptrs2(const void *a, const void *b)
{
    return cmp_regs(*((reg_t**)a),*((reg_t**)b));
}

int regidx_push(regidx_t *idx, char *chr_beg, char *chr_end, hts_pos_t beg, hts_pos_t end, void *payload)
{
    if (beg < 0) beg = 0;
    if (end < 0) end = 0;
    if ( beg > MAX_COOR_0 ) beg = MAX_COOR_0;
    if ( end > MAX_COOR_0 ) end = MAX_COOR_0;

    int rid;
    if (kputsn(chr_beg, chr_end-chr_beg+1, ks_clear(&idx->str)) < 0) return -1;
    if ( khash_str2int_get(idx->seq2regs, idx->str.s, &rid)!=0 )
    {
        // new chromosome
        int m_tmp = idx->mseq;
        if (hts_resize(char*, idx->nseq + 1, &m_tmp,
                       &idx->seq_names, HTS_RESIZE_CLEAR) < 0) {
            return -1;
        }
        if (hts_resize(reglist_t, idx->nseq + 1, &idx->mseq,
                       &idx->seq, HTS_RESIZE_CLEAR) < 0) {
            return -1;
        }
        assert(m_tmp == idx->mseq);
        idx->seq_names[idx->nseq] = strdup(idx->str.s);
        rid = khash_str2int_inc(idx->seq2regs, idx->seq_names[idx->nseq]);
        idx->nseq++;
    }

    reglist_t *list = &idx->seq[rid];
    list->seq = idx->seq_names[rid];
    int mreg = list->mreg;
    if (hts_resize(reg_t, list->nreg + 1, &list->mreg, &list->reg, 0) < 0)
        return -1;
    list->reg[list->nreg].beg = beg;
    list->reg[list->nreg].end = end;
    if ( idx->payload_size ) {
        if ( mreg != list->mreg ) {
            uint8_t *new_dat = realloc(list->dat, idx->payload_size*list->mreg);
            if (!new_dat) return -1;
            list->dat = new_dat;
        }
        memcpy(list->dat + idx->payload_size*list->nreg, payload, idx->payload_size);
    }
    list->nreg++;
    if ( !list->unsorted && list->nreg>1 && cmp_regs(&list->reg[list->nreg-2],&list->reg[list->nreg-1])>0 ) list->unsorted = 1;
    return 0;
}

int regidx_insert(regidx_t *idx, char *line)
{
    if ( !line ) return 0;
    char *chr_from, *chr_to;
    hts_pos_t beg,end;
    int ret = idx->parse(line,&chr_from,&chr_to,&beg,&end,idx->payload,idx->usr);
    if ( ret==-2 ) return -1;   // error
    if ( ret==-1 ) return 0;    // skip the line
    return regidx_push(idx, chr_from,chr_to,beg,end,idx->payload);
}

regidx_t *regidx_init_string(const char *str, regidx_parse_f parser, regidx_free_f free_f, size_t payload_size, void *usr_dat)
{
    kstring_t tmp = KS_INITIALIZE;
    regidx_t *idx = (regidx_t*) calloc(1,sizeof(regidx_t));
    if ( !idx ) return NULL;

    idx->free  = free_f;
    idx->parse = parser ? parser : regidx_parse_tab;
    idx->usr   = usr_dat;
    idx->seq2regs = khash_str2int_init();
    if (!idx->seq2regs) goto fail;
    idx->payload_size = payload_size;
    if ( payload_size ) {
        idx->payload = malloc(payload_size);
        if (!idx->payload) goto fail;
    }

    const char *ss = str;
    while ( *ss )
    {
        while ( *ss && isspace_c(*ss) ) ss++;
        const char *se = ss;
        while ( *se && *se!='\r' && *se!='\n' ) se++;
        if (kputsn(ss, se-ss, ks_clear(&tmp)) < 0) goto fail;
        if (regidx_insert(idx, tmp.s) < 0) goto fail;
        while ( *se && isspace_c(*se) ) se++;
        ss = se;
    }
    ks_free(&tmp);
    return idx;

 fail:
    regidx_destroy(idx);
    ks_free(&tmp);
    return NULL;
}

regidx_t *regidx_init(const char *fname, regidx_parse_f parser, regidx_free_f free_f, size_t payload_size, void *usr_dat)
{
    if ( !parser )
    {
        if ( !fname ) parser = regidx_parse_tab;
        else
        {
            int len = strlen(fname);
            if ( len>=7 && !strcasecmp(".bed.gz",fname+len-7) )
                parser = regidx_parse_bed;
            else if ( len>=8 && !strcasecmp(".bed.bgz",fname+len-8) )
                parser = regidx_parse_bed;
            else if ( len>=4 && !strcasecmp(".bed",fname+len-4) )
                parser = regidx_parse_bed;
            else if ( len>=4 && !strcasecmp(".vcf",fname+len-4) )
                parser = regidx_parse_vcf;
            else if ( len>=7 && !strcasecmp(".vcf.gz",fname+len-7) )
                parser = regidx_parse_vcf;
            else
                parser = regidx_parse_tab;
        }
    }

    kstring_t str = KS_INITIALIZE;
    htsFile *fp = NULL;
    int ret;
    regidx_t *idx = (regidx_t*) calloc(1,sizeof(regidx_t));
    if (!idx) return NULL;
    idx->free  = free_f;
    idx->parse = parser;
    idx->usr   = usr_dat;
    idx->seq2regs = khash_str2int_init();
    if (!idx->seq2regs) goto error;
    idx->payload_size = payload_size;
    if ( payload_size ) {
        idx->payload = malloc(payload_size);
        if (!idx->payload) goto error;
    }

    if ( !fname ) return idx;

    fp = hts_open(fname,"r");
    if ( !fp ) goto error;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &str)) > 0 ) {
        if ( regidx_insert(idx, str.s) ) goto error;
    }
    if (ret < -1) goto error;

    ret = hts_close(fp);
    fp = NULL;
    if ( ret != 0 ) {
        hts_log_error("Close failed .. %s", fname);
        goto error;
    }
    ks_free(&str);
    return idx;

error:
    ks_free(&str);
    if ( fp ) hts_close(fp);
    regidx_destroy(idx);
    return NULL;
}

void regidx_destroy(regidx_t *idx)
{
    int i, j;
    if (!idx) return;
    for (i=0; i<idx->nseq; i++)
    {
        reglist_t *list = &idx->seq[i];
        if ( idx->free )
        {
            for (j=0; j<list->nreg; j++)
                idx->free((char *)list->dat + idx->payload_size*j);
        }
        free(list->dat);
        free(list->reg);
        free(list->idx);
    }
    free(idx->seq_names);
    free(idx->seq);
    free(idx->str.s);
    free(idx->payload);
    khash_str2int_destroy_free(idx->seq2regs);
    free(idx);
}

static int reglist_build_index_(regidx_t *regidx, reglist_t *list)
{
    int i;
    if ( list->unsorted ) {
        if ( !regidx->payload_size ) {
            qsort(list->reg,list->nreg,sizeof(reg_t),cmp_reg_ptrs);
        } else {
            reg_t **ptr = malloc(sizeof(*ptr)*list->nreg);
            if (!ptr) return -1;
            for (i=0; i<list->nreg; i++) ptr[i] = list->reg + i;
            qsort(ptr,list->nreg,sizeof(*ptr),cmp_reg_ptrs2);

            uint8_t *tmp_dat = malloc(regidx->payload_size*list->nreg);
            if (!tmp_dat) { free(ptr); return -1; }
            for (i=0; i<list->nreg; i++) {
                size_t iori = ptr[i] - list->reg;
                memcpy(tmp_dat+i*regidx->payload_size,
                       list->dat+iori*regidx->payload_size,
                       regidx->payload_size);
            }
            free(list->dat);
            list->dat = tmp_dat;

            reg_t *tmp_reg = (reg_t*) malloc(sizeof(reg_t)*list->nreg);
            if (!tmp_reg) { free(ptr); return -1; }
            for (i=0; i<list->nreg; i++) {
                size_t iori = ptr[i] - list->reg;
                tmp_reg[i] = list->reg[iori];
            }
            free(ptr);
            free(list->reg);
            list->reg  = tmp_reg;
            list->mreg = list->nreg;
        }
        list->unsorted = 0;
    }

    list->nidx = 0;
    uint32_t j,k, midx = 0;
    // Find highest index bin.  It's possible that we could just look at
    // the last region, but go through the list in case some entries overlap.
    for (j=0; j<list->nreg; j++) {
        int iend = iBIN(list->reg[j].end);
        if (midx <= iend) midx = iend;
    }
    midx++;
    uint32_t *new_idx = calloc(midx, sizeof(uint32_t));
    if (!new_idx) return -1;
    free(list->idx); // Should be NULL on entry, but just in case...
    list->idx = new_idx;
    list->nidx = midx;

    for (j=0; j<list->nreg; j++) {
        int ibeg = iBIN(list->reg[j].beg);
        int iend = iBIN(list->reg[j].end);
        if ( ibeg==iend ) {
            if ( !list->idx[ibeg] ) list->idx[ibeg] = j + 1;
        } else {
            for (k=ibeg; k<=iend; k++)
                if ( !list->idx[k] ) list->idx[k] = j + 1;
        }
    }

    return 0;
}

int regidx_overlap(regidx_t *regidx, const char *chr, hts_pos_t beg, hts_pos_t end, regitr_t *regitr)
{
    if ( regitr ) regitr->seq = NULL;

    int iseq, ireg;
    if ( khash_str2int_get(regidx->seq2regs, chr, &iseq)!=0 ) return 0;    // no such sequence

    reglist_t *list = &regidx->seq[iseq];
    if ( !list->nreg ) return 0;

    if ( list->nreg==1 )
    {
        if ( beg > list->reg[0].end ) return 0;
        if ( end < list->reg[0].beg ) return 0;
        ireg = 0;
    }
    else
    {
        if ( !list->idx ) {
            if (reglist_build_index_(regidx,list) < 0) return -1;
        }

        int ibeg = iBIN(beg);
        if ( ibeg >= list->nidx ) return 0;     // beg is too big

        // find a matching region
        uint32_t i = list->idx[ibeg];
        if ( !i )
        {
            int iend = iBIN(end);
            if ( iend > list->nidx ) iend = list->nidx;
            for (i=ibeg; i<=iend; i++)
                if ( list->idx[i] ) break;
            if ( i>iend ) return 0;
            i = list->idx[i];
        }
        for (ireg=i-1; ireg<list->nreg; ireg++)
        {
            if ( list->reg[ireg].beg > end ) return 0;   // no match, past the query region
            if ( list->reg[ireg].end >= beg && list->reg[ireg].beg <= end ) break; // found
        }

        if ( ireg >= list->nreg ) return 0;   // no match
    }

    if ( !regitr ) return 1;    // match, but no more info to save

    // may need to iterate over the matching regions later
    itr_t_ *itr = (itr_t_*)regitr->itr;
    itr->ridx = regidx;
    itr->list = list;
    itr->beg  = beg;
    itr->end  = end;
    itr->ireg = ireg;
    itr->active = 0;

    regitr->seq = list->seq;
    regitr->beg = list->reg[ireg].beg;
    regitr->end = list->reg[ireg].end;
    if ( regidx->payload_size )
        regitr->payload = list->dat + regidx->payload_size*ireg;

    return 1;
}

int regidx_parse_bed(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace_c(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace_c(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se )
    {
        // just the chromosome name
        *beg = 0;
        *end = MAX_COOR_0;
        return 0;
    }

    ss = se+1;
    *beg = hts_parse_decimal(ss, &se, 0);
    if ( ss==se ) { hts_log_error("Could not parse bed line: %s", line); return -2; }

    ss = se+1;
    *end = hts_parse_decimal(ss, &se, 0) - 1;
    if ( ss==se ) { hts_log_error("Could not parse bed line: %s", line); return -2; }

    return 0;
}

int regidx_parse_tab(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace_c(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace_c(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se )
    {
        // just the chromosome name
        *beg = 0;
        *end = MAX_COOR_0;
        return 0;
    }

    ss = se+1;
    *beg = hts_parse_decimal(ss, &se, 0);
    if ( ss==se ) { hts_log_error("Could not parse tab line: %s", line); return -2; }
    if ( *beg==0 ) { hts_log_error("Could not parse tab line, expected 1-based coordinate: %s", line); return -2; }
    (*beg)--;

    if ( !se[0] || !se[1] )
        *end = *beg;
    else
    {
        ss = se+1;
        *end = hts_parse_decimal(ss, &se, 0);
        if ( ss==se || (*se && !isspace_c(*se)) ) *end = *beg;
        else if ( *end==0 ) { hts_log_error("Could not parse tab line, expected 1-based coordinate: %s", line); return -2; }
        else (*end)--;
    }
    return 0;
}

int regidx_parse_vcf(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    int ret = regidx_parse_tab(line, chr_beg, chr_end, beg, end, payload, usr);
    if ( !ret ) *end = *beg;
    return ret;
}

int regidx_parse_reg(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace_c(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && *se!=':' ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se )
    {
        *beg = 0;
        *end = MAX_COOR_0;
        return 0;
    }

    ss = se+1;
    *beg = hts_parse_decimal(ss, &se, 0);
    if ( ss==se ) { hts_log_error("Could not parse reg line: %s", line); return -2; }
    if ( *beg==0 ) { hts_log_error("Could not parse reg line, expected 1-based coordinate: %s", line); return -2; }
    (*beg)--;

    if ( !se[0] || !se[1] )
        *end = se[0]=='-' ? MAX_COOR_0 : *beg;
    else
    {
        ss = se+1;
        *end = hts_parse_decimal(ss, &se, 0);
        if ( ss==se ) *end = *beg;
        else if ( *end==0 ) { hts_log_error("Could not parse reg line, expected 1-based coordinate: %s", line); return -2; }
        else (*end)--;
    }
    return 0;
}

regitr_t *regitr_init(regidx_t *regidx)
{
    regitr_t *regitr = (regitr_t*) calloc(1,sizeof(regitr_t));
    if (!regitr) return NULL;
    regitr->itr  = (itr_t_*) calloc(1,sizeof(itr_t_));
    if (!regitr->itr) {
        free(regitr);
        return NULL;
    }
    itr_t_ *itr = (itr_t_*) regitr->itr;
    itr->ridx = regidx;
    itr->list = NULL;
    return regitr;
}

void regitr_reset(regidx_t *regidx, regitr_t *regitr)
{
    itr_t_ *itr = (itr_t_*) regitr->itr;
    memset(itr,0,sizeof(itr_t_));
    itr->ridx = regidx;
}

void regitr_destroy(regitr_t *regitr)
{
    free(regitr->itr);
    free(regitr);
}

int regitr_overlap(regitr_t *regitr)
{
    if ( !regitr || !regitr->seq || !regitr->itr ) return 0;

    itr_t_ *itr = (itr_t_*) regitr->itr;
    if ( !itr->active )
    {
        // is this the first call after regidx_overlap?
        itr->active = 1;
        itr->ireg++;
        return 1;
    }

    reglist_t *list = itr->list;

    int i;
    for (i=itr->ireg; i<list->nreg; i++)
    {
        if ( list->reg[i].beg > itr->end ) return 0;   // no match, past the query region
        if ( list->reg[i].end >= itr->beg && list->reg[i].beg <= itr->end ) break; // found
    }

    if ( i >= list->nreg ) return 0;   // no match

    itr->ireg = i + 1;
    regitr->seq = list->seq;
    regitr->beg = list->reg[i].beg;
    regitr->end = list->reg[i].end;
    if ( itr->ridx->payload_size )
        regitr->payload = (char *)list->dat + itr->ridx->payload_size*i;

    return 1;
}

int regitr_loop(regitr_t *regitr)
{
    if ( !regitr || !regitr->itr ) return 0;

    itr_t_ *itr = (itr_t_*) regitr->itr;
    regidx_t *regidx = itr->ridx;

    if ( !itr->list )    // first time here
    {
        itr->list = regidx->seq;
        itr->ireg = 0;
    }

    size_t iseq = itr->list - regidx->seq;
    if ( iseq >= regidx->nseq ) return 0;

    if ( itr->ireg >= itr->list->nreg )
    {
        iseq++;
        if ( iseq >= regidx->nseq ) return 0; // no more sequences, done
        itr->ireg = 0;
        itr->list = &regidx->seq[iseq];
    }

    regitr->seq = itr->list->seq;
    regitr->beg = itr->list->reg[itr->ireg].beg;
    regitr->end = itr->list->reg[itr->ireg].end;
    if ( regidx->payload_size )
        regitr->payload = (char *)itr->list->dat + regidx->payload_size*itr->ireg;
    itr->ireg++;

    return 1;
}


void regitr_copy(regitr_t *dst, regitr_t *src)
{
    itr_t_ *dst_itr = (itr_t_*) dst->itr;
    itr_t_ *src_itr = (itr_t_*) src->itr;
    *dst_itr = *src_itr;
    *dst = *src;
    dst->itr = dst_itr;
}
