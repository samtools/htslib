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
    hts_pos_t len;
#ifdef CACHE_DBG_LOG
    kstring_t log;
#endif //CACHE_DBG_LOG
} ce_t;

typedef struct cache_t {//cache
    int f, m, n;    //free elements, max elements, count of elem chunks (of 1024)
    ce_t **p;       //array holding chunks of elements
    struct ce_t *head, *tail;
} cache_t;

typedef struct pair_exp {
    int mtid, tid;
    hts_pos_t mpos, pos;
} pair_exp;
KHASH_MAP_INIT_STR(pair, pair_exp)

/// cache status
typedef enum cs {NOTREADY = 0, CACHING, WNDFULL, READY, END} cs;

typedef struct rc_t {//read cache
    cache_t cache;  //cache of mem space
    ce_t *head, *tail;  //alignments
    ce_t *head_sel, *tail_sel;  //selected alignments
    ce_t *head_nsel, *tail_nsel;  //non-selected alignments
    ce_t *head_ins, *tail_ins;  //inserted alignments
    uint64_t ord;   //last ordinal
    cs sts;
    int wndsz, maxdpth, itr;    //size of cache window, depth limit, thr' iterator or not
    hts_pos_t w_st, w_en, dp_en, inc_sz;    //wnd start, end, dpth buffer end, size of inc. buffer
    khash_t(pair) *selpair;     //hash holding name of selected reads for pair selection
    int dp_sz, tid;
    int *inc;   //buffer holding inc val (1), attempt to force intrinsics
    int *dpth;  //depth buffer
} rc_t;

/// @brief setup cache
/// @param fp file pointer to which cache is assigned and used
/// @param wndsz size of cache window
/// @param maxdpth depth limit
/// @return 0 on success others on failure
int setup_readcache(htsFile *fp, int wndsz, int maxdpth);
void destroy_readcache(htsFile *fp);
// return an element to cache
void ret_cache(rc_t *c, ce_t* elem);
// get a cached storage from cache
ce_t* get_cache(htsFile *fp);
//notify end of read
void notify_end(void *c, void *e);
//add a read to cache
int addto_readcache(rc_t *c, ce_t *e, cs *sts);
//retrieve a selected read from cached ones
int getfrom_readcache(rc_t *c, bam1_t *b, hts_pos_t *end);
//process cached reads and select required ones
int process_readcache(rc_t *c);
//get status of cache
cs get_readcache_status(rc_t *c);
//wrapper / for iterators
//get/set access status - thr' iterator or not; if thr' iterator, cache handling is done in itr_nxt
void set_iter_access(htsFile *fp);
int get_iter_access(htsFile *fp);
//get cache pointer, for use in iterator
void* get_sam_readcache(hts_itr_t *itr, void *data);
//retrieve a selected read from cached ones, wrapper for iterator
int getfrom_readcache_iter(void *c, void *s, int *tid, hts_pos_t *beg, hts_pos_t* end);
// get a cached storage from cache, wrapper for iterator
void *get_cache_iter(void *data);
//retrieve bam storage from cache
void *get_readbuffer_iter(void *e);
//notify end of read, wrapper for iterator
void notify_end_iter(void *c, void *e);
//add a read to cache, wrapper for iterator
int addto_readcache_iter(void *c, void *s, cs *sts);
//process cached reads and select required ones, wrapper for iterator
int process_readcache_iter(void *c);
//resets cache status and depth buffer, end of tid/region
void reset_readcache_iter(rc_t *c);

#ifdef CACHE_DBG_LOG
extern FILE *cachelog;
//this is closed by system on exit!
#define LG(...) {if (!cachelog) cachelog = fopen("/tmp/op","w"); if (cachelog) { fprintf(cachelog, __VA_ARGS__);}}
#define LGlog(s,...) ksprintf(s,__VA_ARGS__)
#else
#define LG(...) ;
#define LGlog(s,...) ;
#endif //CACHE_DBG_LOG

#ifdef __cplusplus
}
#endif

#endif //HTSLIB_SAMCACHE_H

