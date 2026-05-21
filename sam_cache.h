/*  sam_cache.h -- Functions to create a cache of reads for depth handling

    Copyright (C) 2026 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE.  */
#ifndef HTSLIB_SAMCACHE_H
#define HTSLIB_SAMCACHE_H

#include "htslib/sam.h"
#include "htslib/khash.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ce_t {//cache element
    uint64_t ord;   //ordinal
    bam1_t *r;
    struct ce_t *next, *prev;
    uint64_t len;
    kstring_t log;
} ce_t;

typedef struct cache_t {//cache
    int f, m, n;    //free, max, chunks(of 1024?)
    ce_t **p;
    struct ce_t *head, *tail;
} cache_t;

typedef struct pair_exp {
    int mtid, tid;
    hts_pos_t mpos, pos;
} pair_exp;
KHASH_MAP_INIT_STR(kh_pair, pair_exp)

typedef struct rc_t {
    cache_t cache;  //cache of mem space
    ce_t *head, *tail;  //alignments
    ce_t *head_sel, *tail_sel;  //selected alignments
    ce_t *head_nsel, *tail_nsel;  //non-selected alignments
    ce_t *head_ins, *tail_ins;  //inserted alignments
    uint64_t selcnt, nselcnt, inscnt,rcnt;
    uint64_t ord;   //last ordinal
    int trgr;  //sts: 0 not ready 1 caching 2 wnd full 3 ready 4 end
    int wndsz, maxdpth, itr;
    hts_pos_t w_st, w_en, dp_en;
    khash_t(kh_pair) *selpair;
    int dp_sz, /*dp_st,*/ tid;
    int *inc, inc_sz;
    int *dpth;
} rc_t;


void destroycache(htsFile *fp);
int setupcache(htsFile *fp, int wndsz, int maxdpth);

void retcache(rc_t *c, ce_t* elem);
ce_t* getcache(htsFile *fp);
int addtoreadcache(rc_t *c, ce_t *e, int *sts);
int getfromreadcache(rc_t *c, bam1_t *b, hts_pos_t *end);
int processcache(rc_t *c);

//wrapper / for iterators
void* getsamcache(hts_itr_t *itr, void *data);
int getfromreadcache_iter(void *c, void *s, int *tid, hts_pos_t *beg, hts_pos_t* end);
void *getcache_iter(void *data);
void *getreadbuffer_iter(void *e);
void notifyend_iter(void *c, void *e);
int addtoreadcache_iter(void *c, void *s, int *sts);
int processcache_iter(void *c);
void resetcache_itr(rc_t *c);

#ifdef CACHE_DBG_LOG
extern FILE *clogfp;
//this is closed by system on exit!
#define LG(...) {if (!clogfp) clogfp = fopen("/tmp/op","w"); if (clogfp) { fprintf(clogfp, __VA_ARGS__);}}
#define LGlog(s,...) ksprintf(s,__VA_ARGS__)
#else
#define LG(...) ;
#define LGlog(s,...) ;
#endif //CACHE_DBG_LOG

#ifdef __cplusplus
}
#endif

#endif //HTSLIB_SAMCACHE_H

