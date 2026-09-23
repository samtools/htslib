/*  sam_cache.c -- Functions to create a cache of reads for depth handling

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

#include <assert.h>

#include "sam_cache.h"
#include "htslib/hts_alloc.h"


#ifdef CACHE_DBG_LOG
FILE *cachelog = NULL;
#endif //CACHE_DBG_LOG

#define CHUNK 1024  //no of items in a cache allocation

/// @brief ensures depth buffer is big enough
/// @param c pointer to read cache
/// @param sz required size of depth buffer
/// @return 0 on success -ve on failure
static int ensure_depthbuffer(rc_t *c, hts_pos_t sz)
{
    if (c->dp_sz && sz <= c->dp_sz)
        return 0;   //big enough
    int *dpth = hts_realloc_p(c->dpth, sizeof(int), sz);
    if (!dpth)
        return -1;
    c->dpth = dpth;
    c->dp_sz = sz;
    return 0;
}
/// @brief setup read cache
/// @param fp htsFile pointer to which cache to be attached
/// @param wndsz size of cache window
/// @param maxdpth depth limit
/// @return 0 on success and non-zero on failure
int setup_readcache(htsFile *fp, int wndsz, int maxdpth)
{
    int i, j, first = 0;
    rc_t *c = (rc_t*)fp->c;
    ce_t *elem = NULL, *tail = NULL, **p = NULL;
    const int def_wndsz = 3500, def_dpth = 1000;
    /*create cache if it doesn't exists. set window and depth size during
    initialisation. when free cache slots are 1, allocate next chunk.*/
    if (!c) { //create cache
        first = 1;
        wndsz = wndsz <= 0 ? def_wndsz : wndsz;
        maxdpth = maxdpth <= 0 ? def_dpth : maxdpth;
        if (!(c = hts_calloc(sizeof(rc_t), 1)))
            goto fail;
        fp->c = c;
    } else if (c->cache.f > 1) {    //re-init or retrieval
        //update depth buffer / dpth settings if needed
        if (wndsz || maxdpth) {
            if (wndsz > 0) {        //succeeding window size arg
                //avoid retaining unnecessary buffer if wndsz is less than set
                if (wndsz < c->wndsz) {
                    free(c->dpth);
                    c->dpth = NULL;
                    c->dp_sz = 0;
                }
                if (ensure_depthbuffer(c, wndsz+1))
                    goto fail;
                memset(c->dpth, 0, (wndsz+1) * sizeof(int));
                c->wndsz = wndsz;
            }
            if (maxdpth > 0) {      //succeeding depth arg
                c->maxdpth = maxdpth;
            }
#ifdef CACHE_DBG_LOG
            LG("setup_readcache wnd:%d dpth:%d\n", c->wndsz, c->maxdpth);
#endif //CACHE_DBG_LOG
        }
        return 0;
    }
    //make cache storage
    if (!c->cache.m) {  //initial
        p = hts_calloc_ps(sizeof(ce_t*), c->cache.n, 1);
    } else {            //growing
        p = hts_realloc_ps(c->cache.p, sizeof(ce_t*), c->cache.n, 1);
    }
    if (!p)
        goto fail;
    c->cache.p = p;     //array holding chunks of storage
    //allocate cache elements
    if((elem = hts_calloc(sizeof(ce_t), CHUNK))) {
        c->cache.p[c->cache.n++] = elem;
        if (!c->cache.m) {  //initial
            if (!(elem->r = bam_init1()))
                goto fail;
            c->cache.head = tail = elem;
            i = 1;
#ifdef CACHE_DBG_LOG
            ks_initialize(&elem->log);
#endif //CACHE_DBG_LOG
        } else {            //growing
            tail = c->cache.tail;
            i = 0;
        }
        //initialize and add to tail
        for (; i < CHUNK; ++i) {
            if (!((elem + i )->r = bam_init1()))
                goto fail;
            tail->next = elem + i;
            tail = tail->next;
#ifdef CACHE_DBG_LOG
            ks_initialize(&elem->log);
#endif //CACHE_DBG_LOG
        }
        c->cache.m += CHUNK;
        c->cache.f += CHUNK;
        c->cache.tail = tail;
    } else
        goto fail;

    if (first) {    //setup starting params
        c->w_st = c->w_en = -1;
        c->tid = -2;
        if (!c->selpair && !(c->selpair = kh_init(pair)))
            goto fail;
        //todo check the +1 allocations
        if (wndsz) {    //window size arg in use
            if (ensure_depthbuffer(c, wndsz+1))
                goto fail;
            memset(c->dpth, 0, (wndsz+1) * sizeof(int));
            c->wndsz = wndsz;
        }
        if (maxdpth) {  //depth arg in use
            c->maxdpth = maxdpth;
        }
#ifdef CACHE_DBG_LOG
        LG("setup_readcache wnd:%d dpth:%d\n", c->wndsz, c->maxdpth);
#endif //CACHE_DBG_LOG
    }

    return 0;

fail:
    if (c) {
        for (i = 0; i < c->cache.n; ++i) {
            elem = c->cache.p[i];
            for (j = 0; j < CHUNK; ++j) {
                bam_destroy1(elem[j].r);
#ifdef CACHE_DBG_LOG
                ks_free(&(elem[j].log));
#endif //CACHE_DBG_LOG
            }
            free(elem);
            c->cache.p[i] = NULL;
        }
        free(c->cache.p);
        c->cache.p = NULL;
        free(c->dpth);
        c->dpth = NULL;
        free(c);
        fp->c = NULL;
    }
    return 1;
}
/// @brief destroys the read cache
/// @param fp htsFile pointer
void destroy_readcache(htsFile *fp)
{
    int i, j;
    ce_t *elem = NULL;
    rc_t *c = (rc_t*) fp->c;
    khint_t iter;

    if(!c)      //cache not in use
        return;

    for (iter = kh_begin(c->selpair); iter != kh_end(c->selpair); ++iter) {
        if (kh_exist(c->selpair, iter)) {
            kh_del(pair, c->selpair, iter);
        }
    }

    while (c->head) {   //clear remaining reads
        elem = c->head->next;
        LGlog(&c->head->log, "%s", "cleanup");
        ret_cache(c, c->head);
        c->head = elem;
    }
    while (c->head_nsel) {  //clear non-selected reads
        elem = c->head_nsel->next;
        LGlog(&c->head_nsel->log, "%s", "cleanup");
        ret_cache(c, c->head_nsel);
        c->head_nsel = elem;
    }
    //cleanup cache
    for (i = 0; i < c->cache.n; ++i) {
        elem = c->cache.p[i];
        for (j = 0; j < CHUNK; ++j) {
            bam_destroy1(elem[j].r);
#ifdef CACHE_DBG_LOG
            ks_free(&elem[j].log);
#endif //CACHE_DBG_LOG
        }
        free(elem);
    }
    free(c->cache.p);
    free(c->dpth);
    free(c->inc);
    kh_destroy(pair, c->selpair);
    free(c);
    fp->c = NULL;
}

//implementation / internals
/// @brief marks the access as thr' iterator
/// @param fp htsFile pointer
void set_iter_access(htsFile *fp) {
    /*when cache is used thr' iterator, cached read handling is done in itr_nxt
      when it is used on whole file, it is done in sam_read1. this flag
      helps to identify these scenarios and use cache appropriately*/
    if (fp->c) {    //cache is in use
        rc_t *c = (rc_t*)fp->c;
        c->itr = 1;
    }
}
/// @brief retrieves the how the cache is accessed, thr' iterator/for whole file
/// or not in use at all
/// @param fp htsFile pointer for cache access
/// @return 1 if thr' iterator and 0 if not in use / not thr' iterator
int get_iter_access(htsFile *fp) {
    if (fp->c) {    //cache in use
        rc_t *c = (rc_t*) fp->c;
        return c->itr;
    }
    return 0;       //cache not in use
}
/// @brief cache's status
/// @param c read cache
/// @return cache status enum showing status of cache
cs get_readcache_status(rc_t *c) {
    return c ? c->sts : NOTREADY;
}

//todo htsopt3 to try
//-1 on failure and 1 when required and 0 on skip
/// @brief update the depth buffer based on reads length
/// @param c read cache
/// @param e cache element holding the read under processing
/// @param chk 1 checks whether the read is required or not; 0 to update depth
/// @return -ve - error, 0 - read not required, 1 - read required
static int update_depth(rc_t *c, ce_t *e, int chk)
{
    int *dpth = NULL;
    uint32_t *cgr = bam_get_cigar(e->r), i, j;
    int clen, off, req = 0;
    hts_pos_t st, en, len = 0, h;
    if (!c->dpth) {     //setup depth buffer
        if (ensure_depthbuffer(c, c->dp_sz + 1))
            goto fail;
        c->dp_en = c->w_st + c->dp_sz;
        off = 0;
    }
    st = c->w_st; en = c->dp_en;
    if (st > e->r->core.pos)
        goto fail;  //not sorted?
    //inc buffer is an attempt to force compiler to use intrinsics
    if (e->len > c->inc_sz) {   //grow increment buffer as required
        int *inc = hts_realloc_p(c->inc, sizeof(int), e->len);
        if (!inc) goto fail;
        c->inc = inc;
        c->inc_sz = e->len;
        for(h = 0; h < c->inc_sz; ++h)
            *(c->inc + h) = 1;
    }
    if (st < e->r->core.pos) {
        off = e->r->core.pos - st;
        len = off + e->len;
    } else {
        off = st - e->r->core.pos;
        len = e->len - off;
        off *= -1;
    }
    if (en < (c->w_st+len)) {   //goes over the end of buffer, increase it
        len = c->w_st + len - en;
        hts_pos_t bkpsz = c->dp_sz;
        if (ensure_depthbuffer(c, len + c->dp_sz))
            goto fail;
        memset(c->dpth + bkpsz, 0, len * sizeof(int));
        c->dp_en = en = c->w_st + c->dp_sz;
    }
    len = 0;
    dpth = c->dpth;
    if (chk) {  //check depth
        for (i = 0; i < e->r->core.n_cigar; ++i) {
            if (!(bam_cigar_type(bam_cigar_op(cgr[i])) & 2)) {
                continue;       //not consuming ref
            }
            //deletion is counted!
            //check depth for each position is above required limit or not
            //read required if depth is <= the limit
            clen = bam_cigar_oplen(cgr[i]);
            for (j = 0; j < clen; ++j) {
                if (off + j >= 0)
                    req |= dpth[off + j] + 1 <= c->maxdpth;
            }
            off += clen;
            if (req) {  //required, no need to check further
                break;
            }
        }
    } else {    //update depth
        req = 1;
        for (i = 0; i < e->r->core.n_cigar; ++i) {
            if (!(bam_cigar_type(bam_cigar_op(cgr[i])) & 2)) {   //not consuming ref
                continue;
            }
            clen = bam_cigar_oplen(cgr[i]);
            //inc buffer is an attempt to force compiler to use intrinsics
            for (j = 0; j < clen; ++j) {
                dpth[off + len + j] += c->inc[j];
            }
            len += clen;
        }
    }
    if (req)
        return 1;
    return 0;

fail:
    return -1;
}
/// @brief return the storage back to cache
/// @param c read cache
/// @param elem element/space in cache
void ret_cache(rc_t *c, ce_t* elem)
{
    if ((elem->r->core.flag & BAM_FPAIRED) && !(elem->r->core.flag & BAM_FUNMAP)) {
        //paired and mate mapped, remove from expected pair
        khiter_t it = kh_get(pair, c->selpair, bam_get_qname(elem->r));
        if (it != kh_end(c->selpair) && kh_exist(c->selpair, it)) {
            kh_del(pair, c->selpair, it);
        }
    }
    //add as head in cache
    elem->prev = NULL;
    elem->ord = 0;
    elem->len = 0;
    elem->next = c->cache.head;
    c->cache.head->prev = elem;
    c->cache.head = elem;
    ++c->cache.f;
#ifdef CACHE_DBG_LOG
    ks_clear(&elem->log);
#endif //CACHE_DBG_LOG
}

/// @brief get cache element / storage from preallocated cache
/// @param fp htsFile pointer to setup / retrieve cache
/// @return ce_t* on success or NULL on failure
ce_t* get_cache(htsFile *fp)
{
    rc_t *c = (rc_t*)fp->c;
    ce_t *ret = NULL;

    //create enough space
    if (setup_readcache(fp, 0, 0))   //passing 0 to avoid re-initialization
        goto fail;

    ret = c->cache.head;
    c->cache.head = c->cache.head->next;
    c->cache.head->prev = NULL;
    ret->prev = NULL;
    ret->next = NULL;
    --c->cache.f;
    return ret;
fail:
    return NULL;
}
/// @brief mark end of input
/// @param p read cache pointer
/// @param e allocated and unused space
void notify_end(void *p, void *e)
{
    rc_t *c = (rc_t*)p;
    c->sts = END;   //end
    c->tid = -3;    //reset that it doesn't match to any/initial vals
    //real end of input, unlike iterator where it could be just end of a region
    ret_cache(c, (ce_t*)e);
}

/// @brief add a read to cache
/// @param c pointer to read cache
/// @param e cache element containing the read to be cached
/// @param sts to return status of cache post caching
/// @return -1 on failure 0 on success
int addto_readcache(rc_t *c, ce_t *e, cs *sts)
{
    int unmap = 0;
    if (!(e->r->core.flag & BAM_FUNMAP)) {
        if (c->w_st == -1) {
            //starting, use pos of 1st or one being added as start of window
            c->w_st = c->head ? c->head->r->core.pos : e->r->core.pos;
            c->w_en = c->w_st + c->wndsz;   //end of wnd
            c->dp_en = c->w_st + c->dp_sz;  //end of depth buffer
        }
    } else {
        unmap = 1;  //unmapped, add w/o depth check
    }
    e->ord = ++(c->ord);
    e->len = bam_cigar2rlen(e->r->core.n_cigar, bam_get_cigar(e->r));
    LG("+ %s %"PRIu64"\t\t%"PRIhts_pos" %"PRIhts_pos" %"PRIu64" %"PRIhts_pos"\n", bam_get_qname(e->r), e->ord, c->w_st, e->r->core.pos, e->len, c->w_en);
    if (!c->head) {
        c->head = c->tail = e;
    } else {
        ce_t *p = c->tail;
        ce_t *tmpn = NULL;
        //add to the tail
        if (!unmap) {
            if (p->r->core.tid == e->r->core.tid && p->r->core.pos > e->r->core.pos) {
                hts_log_error("Unsorted data");
                return -1;   //not sorted!
            }
        }
        if (p) {    //useful if read is sorted based on len and being inserted
            tmpn = p->next;
            p->next = e;
            e->prev = p;
            e->next = tmpn;
            if (tmpn)
                tmpn->prev = e;
            if (p == c->tail)
                c->tail = e;
        } else { //either last or 1st
            if (c->head == c->tail && !c->tail) {   //none in list
                c->tail = c->head = e;
                e->prev = e->next = NULL;
            } else {    //insert 1st
                tmpn = c->head;
                c->head = e;
                e->prev = NULL;
                e->next = tmpn;
                if(tmpn)
                    tmpn->prev = e;
            }
        }
    }
    //todo do we need a limit on max no of items that are cached? like the whole file is for same pos, probably cant be loaded!
    if (c->tid == e->r->core.tid) {
        if (c->w_en < e->r->core.pos) {  //post window, process and advance
            LG("wnd full\n");
            c->sts = WNDFULL;   //wnd full, go for processing
        }
        else
            c->sts = CACHING;   //caching
    } else if (c->tid != -2) {
        LG("tid change\n");
        c->sts = READY;         //ready for processing
    }
    else
        c->sts = CACHING;       //caching

    c->tid = e->r->core.tid;
    if (sts)
        *sts = c->sts;
    LGlog(&e->log, "%s,%"PRIu64",added,%d,%d,%"PRIhts_pos",%"PRIhts_pos",%"PRIu64",%"PRIhts_pos",", bam_get_qname(e->r), e->ord,e->r->core.tid, e->r->core.flag,e->r->core.pos, e->r->core.mpos,e->len, e->len+e->r->core.pos);
    return 0;
}

/// @brief get read from processed cache
/// @param c pointer to read cache
/// @param b pointer to bam data, to which read data is copied
/// @param end end of read, for iterators
/// @return -1 on failure, 0 when nothing to retrieve and 1 with read retrieved
int getfrom_readcache(rc_t *c, bam1_t *b, hts_pos_t *end)
{
    if (!c || c->sts < READY) {    //not ready!
        return 0;
    }
    //todo at some point, removal from selpair need to be done based on pos as well
    uint64_t sel = UINT64_MAX, ins = UINT64_MAX;
    ce_t *e = c->head_sel, *f = c->head_ins, *p = NULL;

    //get from selected or inserted list, based on ordinal
    if (e)
        sel = e->ord;
    if (f)
        ins = f->ord;
    if(sel < ins)
        p = e;
    else
        p = f;

    if (p && (c->sts == READY || c->sts == END)) {
        //send only upto start of wnd to maintain the order, except when it is end
        if (!bam_copy1(b, p->r))
            return -1;
        if (p == e) {   //remove from sel list
            c->head_sel = p->next;
            if (!c->head_sel) {
                c->tail_sel = NULL;
            }
        } else {        //remove from ins list
            c->head_ins = p->next;
            if (!c->head_ins) {
                c->tail_ins = NULL;
            }
        }
        if (!c->head_sel && !c->head_ins && c->sts != END) {
            c->sts = NOTREADY;    //not ready
        } else {
            if (c->head_sel) hts_prefetch(c->head_sel);
            if (c->head_ins) hts_prefetch(c->head_ins);
        }
        LG("- %s %"PRIu64"\n", bam_get_qname(b), p->ord);
        LGlog(&p->log, "%s", ",retrieved");
        if (end) *end = p->len + p->r->core.pos;
        ret_cache(c, p); //return storage to cache
        return 1;
    }
    return 0;
}
/// @brief find a read matching to given one from non-selected list
/// @param c pointer to read cache
/// @param e pointer to read for which pair need to be found
/// @param ep pointer to previous one of the pair, to maintain list
/// @return NULL when not found and cache element pointer when found
static inline ce_t* find_nsel(rc_t *c, ce_t *e, ce_t **ep)
{
    ce_t *s = c->head_nsel;
    *ep = NULL;
    while (s) {
        if (s->next)
            hts_prefetch(s->next);
        if (s->ord > e->ord)
            break;  //not found
        if (s->r->core.pos == e->r->core.mpos &&
            s->r->core.mpos == e->r->core.pos &&
            s->r->core.tid == e->r->core.mtid &&
            s->r->core.mtid == e->r->core.tid &&
            !strcmp(bam_get_qname(s->r), bam_get_qname(e->r)))
            return s;   //found
        *ep = s;
        s = s->next;
    }
    return NULL;
}
/// @brief move the read/cache element from main list to selected/unselected/insert list
/// @param c pointer to cache
/// @param ep pointer to previous element to maintain the list
/// @param e element being moved
/// @param en next element
/// @param sel 1 to move read to selected list 0 to move to nselected list
/// @param ins 1 to move to insert list, relevant with sel = 1
static inline void move_read(rc_t *c, ce_t *ep, ce_t *e, ce_t* en, int sel, int ins)
{
    int paired = (e->r->core.flag & BAM_FPAIRED) &&
        !(e->r->core.flag & BAM_FUNMAP) && !(e->r->core.flag & BAM_FMUNMAP) &&
        (e->r->core.mtid != -1) && (e->r->core.mpos != -1);
    if (!ins) {     //remove from cache
        if (en)
            en->prev = ep;
        if (c->head == e)
            c->head = en;
        if (c->tail == e)
            c->tail = ep;
        if (!c->head)
            c->tail = c->head;
        if(ep) {
            ep->next = en;
        }
    } else {        //remove from nsel
        if (ep) {
            ep->next = en;
        } else {
            c->head_nsel = en;
        }
        if (en)
            en->prev = ep;
        else
            c->tail_nsel = ep ? ep : NULL;

    }
    e->next = NULL;
    e->prev = NULL;

    if (sel) {      //moving to sel/ins list
        //insert in required pos, starting from tail, in order of ordinal
        ce_t *s = ins? c->tail_ins : c->tail_sel, *p = NULL;
        if (s && s->ord < e->ord) { //shortcut
            s->next = e;
            e->next = NULL;
            e->prev = s;
            if(ins)
                c->tail_ins = e;
            else
                c->tail_sel = e;
            return;
        }
        while (s) {
            if (s->ord < e->ord) {  //add in ascending order
                break;
            }
            s = s->prev;
        }
        if (!s) {   //as head
            if (ins) {
                p = c->head_ins;
                c->head_ins = e;
            }
            else {
                p = c->head_sel;
                c->head_sel = e;
            }
            e->prev = NULL;
            e->next = p;
            if(p)
                p->prev = e;
            else {  //update tail
                if (ins)
                    c->tail_ins = e;
                else
                    c->tail_sel = e;
            }
            return;
        } else {
            p = s->next;
            s->next = e;
            e->prev = s;
            e->next = p;
            if (p)
                p->prev = e;
            return;
        }
        return;
    } else if (paired) {
        //move to nsel if paired, otherwise discard and return cache
        //add to non-selected list, for pair lookup
        ce_t *s = c->tail_nsel, *p = NULL;
        if (s && s->r->core.pos < e->r->core.pos &&
            s->r->core.tid == e->r->core.tid) {   //shortcut
            s->next = e;
            e->next = NULL;
            e->prev = s;
            c->tail_nsel = e;
            return;
        }
        //find pos and fit, in order of increasing pos, that it is easy to remove
        while (s && (s->r->core.tid == e->r->core.tid)) {
            if (s->r->core.pos < e->r->core.pos) {
                break;
            }
            s = s->prev;
        }
        if (!s) {   //add as head
            p = c->head_nsel;
            c->head_nsel = e;
            e->prev = NULL;
            e->next = p;
            if(p)
                p->prev = e;
            if (!p)
                c->tail_nsel = e;
            return;
        } else {
            p = s->next;
            s->next = e;
            e->prev = s;
            e->next = p;
            if (p)
                p->prev = e;
            return;
        }
        return;
    } else { //non selected, non paired reads, release them
        LGlog(&e->log, "%s", "npair,disc");
        if(ep)
            ep->next = en;
        if (en)
            en->prev = ep;
        ret_cache(c, e);
        return;
    }
}
/// @brief reset cache status, for next tid/iterator...
/// @param c pointer to read cache
static inline void reset_depth(rc_t* c)
{
    c->w_st = -1;
    if (c->dp_sz <= 0 || !c->dpth)
        return;
    memset(c->dpth, 0, c->dp_sz * sizeof(int));

    ce_t *en = NULL;
    //clear all from previous tid
    while (c->head && c->head->r->core.tid != c->tid) {
        en = c->head->next;
        LGlog(&c->head->log, "%s", "h-reset");
        ret_cache(c, c->head);
        c->head = en;
    }
    if (!c->head) c->tail = NULL;
    else c->head->prev = NULL;

    //clear whole non selected ones
    while (c->head_nsel) {
        en = c->head_nsel->next;
        LGlog(&c->head_nsel->log, "%s", "n-reset");
        ret_cache(c, c->head_nsel);
        c->head_nsel = en;
    }

    LG("reset: t %"PRIu64" s %"PRIu64" i %"PRIu64" n %"PRIu64"; nxt %d\n", c->rcnt, c->selcnt,c->inscnt, c->nselcnt, c->tid);
    c->tail_nsel = NULL;
}
/// @brief process the cached reads and find required ones
/// @param c read cache
/// @return 0 on success and -ve on error
int process_readcache(rc_t *c)
{
    ce_t *e = NULL, *ep = NULL, *en = NULL;
    hts_pos_t pos, off;
    khiter_t pairitr;
    if (!c->head)
        return 0;

    hts_pos_t endpos = c->w_en < c->tail->r->core.pos ? c->tail->r->core.pos - 1 : c->w_en;
    hts_pos_t lastpos = 0;

    if(c->sts == END) {
        //fine tune endpos for last iteration, by looking for valid len which may not be the tail one!
        e = c->head;
        lastpos = c->head->r->core.pos + c->head->len;
        while (e) {
            pos = e->r->core.pos + e->len;
            if (lastpos < pos)
                lastpos = pos;
            e = e->next;
        }
        if (endpos < lastpos) {
            endpos =  lastpos - 1;
        }
    }
    pos = c->w_st;
    while (pos <= endpos) {
        if (!(e = c->head))
            break;
        ep = NULL;
        //discard any irrelevant ones
        while (e && (e->r->core.pos <= pos) && ((c->tid != e->r->core.tid) || c->sts == WNDFULL)) {
            en = e->next;
            if (e->r->core.flag & BAM_FUNMAP) { //unmapped, nothing further to check
                LGlog(&e->log, "%s%"PRIhts_pos, "sel-unmapped,", pos);
                LG("* x %"PRIhts_pos" selunmapped\n", e->ord);
                move_read(c, ep, e, en, 1, 0);
            } else if (e->r->core.pos + e->len - 1 < pos) {
                //not relevant for this pos or succeeding ones
                //check whether pair is selected before discarding
                pairitr = kh_get(pair, c->selpair, bam_get_qname(e->r));
                if (pairitr != kh_end(c->selpair)) {
                    pair_exp *p = &kh_val(c->selpair, pairitr);
                    if (p->mpos == e->r->core.pos &&
                        p->mtid == e->r->core.tid &&
                        p->pos == e->r->core.mpos &&
                        p->tid == e->r->core.mtid) {   //pair already selected
                        kh_del(pair, c->selpair, pairitr);    //remove from expected pairs
                        LGlog(&e->log, "%s%"PRIhts_pos, "sel as paired,", pos);
                        move_read(c, ep, e, en, 1, 0);   //select
                        LG("* x %"PRIhts_pos" selpaired\n", e->ord);
                        //no depth update!
                    } else {
                        LGlog(&e->log, "%s%"PRIhts_pos, "nsel,", pos);
                        move_read(c, ep, e, en, 0, 0);   //move to unselected
                        LG("* x %"PRIhts_pos" nsel\n", e->ord);
                    }
                } else {
                    LGlog(&e->log, "%s%"PRIhts_pos",", "nsel,", pos);
                    move_read(c, ep, e, en, 0, 0);   //move to unselected
                    LG("* x %"PRIhts_pos" nsel2\n", e->ord);
                }
            } else {
                ep = e;
            }
            e = en;
        }
        e = c->head;
        ep = NULL;
        off = pos - c->w_st;
        if (c->head == c->tail && c->sts == WNDFULL) {
            //lastone --> all from wnd are done and last one to be considered in nxt iteration
            endpos = pos;
            break;
        } else {
            if (c->dp_sz <= off) {
                hts_pos_t bkp = c->dp_sz, ln = 100;
                if (ensure_depthbuffer(c, c->dp_sz + ln)) {
                    goto fail;
                }
                memset(c->dpth + bkp, 0, ln * sizeof(int));
                c->dp_en += ln;
            }
        }
        if (c->dpth[off] >= c->maxdpth) { //have enough depth
            ++pos;
            continue;
        }
        //find the last one covering the pos
        while (e && (e->r->core.pos <= pos) && ((c->tid != e->r->core.tid) || c->sts == WNDFULL)) {
            en = e->next;
            if (e->r->core.flag & BAM_FUNMAP) {
                LGlog(&e->log, "%s%"PRIhts_pos, "sel-unmapped,", pos);
                LG("* x %"PRIhts_pos" selunmapped\n", e->ord);
                move_read(c, ep, e, en, 1, 0);
                e = en;
                continue;
            }
            if (e->r->core.pos + e->len - 1 >= pos) {
                ep = e;
            }
            e = en;
        }
        if (!ep) {  //nothing!
            ++pos;
            continue;
        }
        LG("* x %"PRIhts_pos" sel @ %"PRIhts_pos"\n", ep->ord, pos);
        LGlog(&ep->log, "%s%"PRIhts_pos, "sel,", pos);
        move_read(c, ep->prev, ep, ep->next, 1, 0);
        if (update_depth(c, ep, 0) < 0) {
            goto fail;
        }
        //get/set for pair
        if (ep->r->core.flag & BAM_FPAIRED && !(ep->r->core.flag & BAM_FMUNMAP)) {
            if (ep->r->core.mpos >= pos) {  //upcoming mate, add to hash
                int r = -1;
                pairitr = kh_put(pair, c->selpair, bam_get_qname(ep->r), &r);
                if (r == -1)
                    goto fail;
                pair_exp *p = &kh_val(c->selpair, pairitr);
                p->pos = ep->r->core.pos; p->tid = ep->r->core.tid;
                p->mpos = ep->r->core.mpos; p->mtid = ep->r->core.mtid;
            } else {    //mate already passed, selected or not?
                int sel = 0;
                pairitr = kh_get(pair, c->selpair, bam_get_qname(ep->r));
                if (pairitr != kh_end(c->selpair)) {
                    pair_exp *p = &kh_val(c->selpair, pairitr);
                    if (p->mpos == ep->r->core.pos &&
                        p->mtid == ep->r->core.tid &&
                        p->pos == ep->r->core.mpos &&
                        p->tid == ep->r->core.mtid) {   //pair already selected
                        sel = 1;
                    }
                }
                if (!sel) { //find and select from unselected ones
                    ce_t *o = NULL, *op = NULL;
                    if ((o = find_nsel(c, ep, &op))) {
                        LGlog(&o->log, "%s%"PRIhts_pos",%"PRIhts_pos",", "sel from nsel as paired ", ep->ord, pos);
                        move_read(c, op, o, o->next, 1, 1);   //select
                        //no depth update!
                    } else {
                        LG("searching for pair failed, %"PRIu64"\n", ep->ord);
                    }
                }
                if (pairitr != kh_end(c->selpair)) {    //remove from hash
                    kh_del(pair, c->selpair, pairitr);  //remove from expected pairs
                }
            }
        }
    }
    assert((c->head && c->tail) || (c->sts != WNDFULL));
    //check for any pairs discarded due to sufficient depth and select
    LG("pos at lpexit %"PRIhts_pos", head %llu\n", pos, c->head?c->head->ord:0);
    e = c->head;
    ep = NULL;
    if (c->head)
        --pos;  //out of loop --> ++pos
    while (e && (e->r->core.pos <= pos) && ((c->tid != e->r->core.tid) || c->sts == WNDFULL)) {
        en = e->next;
        pairitr = kh_get(pair, c->selpair, bam_get_qname(e->r));
        if (pairitr != kh_end(c->selpair)) {
            pair_exp *p = &kh_val(c->selpair, pairitr);
            if (p->mpos == e->r->core.pos &&
                p->mtid == e->r->core.tid &&
                p->pos == e->r->core.mpos &&
                p->tid == e->r->core.mtid) {   //pair already selected
                kh_del(pair, c->selpair, pairitr);    //remove from expected pairs
                LGlog(&e->log, "%s%"PRIhts_pos, "sel as pair on wnd move,", pos);
                move_read(c, ep, e, en, 1, 0);   //select
                //no depth update!
            } else {
                ep = e;
            }
        } else {
            ep = e;
        }
        e = en;
    }
    if (c->sts != WNDFULL) { //when it is not wnd full, tid change or end
        reset_depth(c);
    } else {    //update window
        en = NULL;
        hts_pos_t adj = c->tail ? c->tail->r->core.pos - c->w_en : 0; //last one, out of window - current end
        hts_pos_t new_st = c->w_st + adj;
        hts_pos_t bkp_st = c->w_st;
        uint64_t last = c->tail_sel ? c->tail_sel->ord : 0;
        uint64_t lasti = c->tail_ins ? c->tail_ins->ord : 0;
        if (lasti > last)
            last = lasti;
        int rem = 0;
        if (c->head_nsel) {
            while ( c->head_nsel->r->core.pos + c->head_nsel->len - 1 < new_st || c->head_nsel->ord < last) {
                rem = 1;
                en = c->head_nsel->next;
                LG("* nsel discarded %s %"PRIu64"\n", bam_get_qname(c->head_nsel->r), c->head_nsel->ord);
                LGlog(&c->head_nsel->log,"%s","nsel-disc");
                ret_cache(c, c->head_nsel);
                if(!(c->head_nsel = en)) {
                    c->tail_nsel = NULL;
                    break;
                } else {
                    c->head_nsel->prev = NULL;
                }
            }
        }
        if (c->head) {//////
            while ( c->head->r->core.pos + c->head->len - 1 < new_st || c->head->ord < last) {
                rem = 1;
                en = c->head->next;
                LG("* discarded %s %"PRIu64" last %"PRIu64" new_st %"PRIhts_pos"\n", bam_get_qname(c->head->r), c->head->ord, last, new_st);
                LGlog(&c->head->log,"%s,%"PRIhts_pos",%s", "wndchange", pos,"disc");
                ret_cache(c, c->head);
                if(!(c->head = en)) {
                    c->tail = NULL;
                    break;
                } else {
                    c->head->prev = NULL;
                }
            }
        }
        if (rem) {
          LG("* wnd full, removed items from head_nsel\n");
        }
        else {
           LG("* wnd full, 0 removed items from head_nsel, [%"PRIhts_pos"-%"PRIhts_pos"] %"PRIhts_pos"\n", c->w_st, c->w_en, c->head_nsel?c->head_nsel->r->core.pos : 0);
        }
        c->w_st = c->head ? c->head->r->core.pos : new_st;    //move wnd
        c->w_en = c->w_st + c->wndsz;
        adj = c->w_st - bkp_st;
        if (adj >= c->dp_sz) {
            memset(c->dpth, 0, c->dp_sz * sizeof(int));
            c->dp_en = c->w_st + c->dp_sz;
        } else {
            LG("adj %"PRIhts_pos", mv %"PRIhts_pos"-%"PRIhts_pos",", adj, c->w_st+adj, c->w_st+c->dp_sz);
            LG("0 set %"PRIhts_pos" - %"PRIhts_pos"\n", c->w_st+c->dp_sz-adj,c->w_st+c->dp_sz);
            memmove(c->dpth, c->dpth + adj, (c->dp_sz - adj) * sizeof(int));
            memset(c->dpth + c->dp_sz - adj, 0, adj * sizeof(int));
            c->dp_en += adj;
        }
        LG("* wnd moved, %"PRIhts_pos" - %"PRIhts_pos", dpth %"PRIhts_pos" - %"PRIhts_pos"; s %"PRIu64" i %"PRIu64" ns %"PRIu64"\n", c->w_st, c->w_en, c->w_st, c->dp_en, c->selcnt, c->inscnt, c->nselcnt);
        c->sts = READY;    //reset full status n get already processed
    }

    return 0;
fail:
    return -1;
}

//wrappers for iterators
//these wrappers help to avoid complexity in hts / iterator code by moving the
//cache structure access to this file

/// @brief get read cache
/// @param itr iterator which is in use
/// @param data custom data with iterator (htsFile/kstring based on invoker)
/// @return read cache pointer
void* get_sam_readcache(hts_itr_t *itr, void *data)
{
    htsFile *fp = NULL;
    /*invoked from itr_nxt, which is used by utilities like tabix as well.
    cache is in use only for sam data and to identify the invocation usecache
    flag is used. this flag is set when iterator is used in sam context.*/
    if (itr && itr->usecache)
        fp = (htsFile*)data;
    return fp ? fp->c : NULL;
}
/// @brief wrapper to get read from cache
/// @param p read cache pointer
/// @param s bam storage for output
/// @param tid tid of output read
/// @param beg beg of output read
/// @param end end of output read
/// @return -1 on failure, 0 when nothing to retrieve and 1 with read retrieved
int getfrom_readcache_iter(void *p, void *s, int *tid, hts_pos_t *beg, hts_pos_t* end)
{
    rc_t *c = (rc_t*)p;
    bam1_t *b = (bam1_t*)s;
    int ret = getfrom_readcache(c, b, end);
    if (ret > 0) {
        *tid = b->core.tid;
        *beg = b->core.pos;
    }
    return ret;
}
/// @brief wrapper to get cache storage
/// @param data custom data for iterator, htsfile pointer for sam data
/// @return storage released from cache or NULL
void *get_cache_iter(void *data)
{
    htsFile *fp = (htsFile*)data;
    void *p = get_cache(fp);
    return p;
}
/// @brief retrieves bam pointer from cache storage
/// @param p cache storage retrieved
/// @return bam pointer as void *
void *get_readbuffer_iter(void *p)
{
    ce_t* e = (ce_t*)p;
    return e->r;
}
/// @brief mark end of input, end of region or file
/// @param p read cache pointer
/// @param e cache storage in use
void notify_end_iter(void *p, void *e)
{
    rc_t *c = (rc_t*)p;
    c->sts = END;   //end
    c->tid = -2;    //reset as in start
    ret_cache(c, (ce_t*)e);
}
/// @brief wrapper to add read to cache
/// @param c cache
/// @param s read to be added
/// @param sts status of cache, output
/// @return -1 on failure 0 on success
int addto_readcache_iter(void *c, void *s, cs *sts)
{
    return addto_readcache(c, s, sts);
}
/// @brief wrapper to process cached data
/// @param c read cache
/// @return 0 on success and -ve on failure
int process_readcache_iter(void *c)
{
    return process_readcache((rc_t*)c);
}
/// @brief wrapper to reset cache status
/// @param c read cache
void reset_readcache_iter(rc_t *c)
{
    reset_depth(c);
    c->sts = NOTREADY;
    c->tid = -2;
    c->w_st = c->w_en = -1;
}
