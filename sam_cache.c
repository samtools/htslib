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

#ifdef CACHE_DBG_LOG
FILE *clogfp = NULL;
#endif //CACHE_DBG_LOG

//setup/destroy
int setupcache(htsFile *fp, int wndsz, int maxdpth)
{
    int i, j;
    rc_t *c = (rc_t*)fp->c;
    ce_t *elem = NULL, *tail = NULL;
    const int inc = 1024;
    if (!c) { //create cache
        if (!(c = calloc(1, sizeof(rc_t))))
            goto fail;
        fp->c = c;
    }

    if (!(c->cache.p = malloc(sizeof(ce_t*))))
        goto fail;
    if((elem = calloc(inc, sizeof(ce_t)))) {
        c->cache.p[c->cache.n++] = elem;
        c->cache.head = tail = elem;
        if (!(elem->r = bam_init1()))
            goto fail;
        ks_initialize(&elem->log);
        for (i = 1; i < inc; ++i) {
            if (!((elem + i )->r = bam_init1()))
                goto fail;
            ks_initialize(&(elem+i)->log);
            tail->next = elem + i;
            tail = tail->next;
        }
        c->cache.m += inc;
        c->cache.f += inc;
        c->cache.tail = tail;
    } else
        goto fail;

    c->wndsz = wndsz;
    c->maxdpth = maxdpth;
    c->w_st = c->w_en = -1;
    c->tid = -2;    //start
    c->selpair = kh_init(kh_pair);
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    return 0;

fail:
    if (c) {
        for (i = 0; i < c->cache.n; ++i) {
            elem = c->cache.p[i];
            for (j = 0; j < inc; ++j) {
                bam_destroy1(elem[j].r);
                ks_free(&(elem[j].log));
            }
            free(elem);
        }
        free(c->cache.p);
        free(c);
        fp->c = NULL;
    }
    return 1;
}

void destroycache(htsFile *fp)
{
    int i, j;
    const int inc = 1024;
    ce_t *elem = NULL;
    rc_t *c = (rc_t*) fp->c;
    khint_t iter;
    if(!c)
        return;
    for (iter = kh_begin(c->selpair); iter != kh_end(c->selpair); ++iter) {
        if (kh_exist(c->selpair, iter)) {
            kh_del(kh_pair, c->selpair, iter);
        }
    }
    for (i = 0; i < c->cache.n; ++i) {
        elem = c->cache.p[i];
        for (j = 0; j < inc; ++j) {
            bam_destroy1(elem[j].r);
            ks_free(&elem[j].log);
        }
        free(elem);
    }
    free(c->cache.p);
    free(c->dpth);
    free(c->inc);
    kh_destroy(kh_pair, c->selpair);
    free(c);
    fp->c = NULL;
}

//implementation / internals
//todo htsopt3 to try
//-1 on failure and 1 when required and 0 on skip
static int updatedepth(rc_t *c, ce_t *e, int chk)
{
    int *dpth = NULL;
    uint32_t *cgr = bam_get_cigar(e->r), i, j;
    int clen, off, req = 0;
    hts_pos_t st, en, len = 0;
    if (!c->dpth) { //setup depth buffer
        c->dp_sz = c->wndsz;
        if (!(c->dpth = calloc(c->dp_sz, sizeof(int))))
            goto fail;
        c->dp_en = c->w_st + c->dp_sz; off = 0;
    }
    dpth = c->dpth;
    st = c->w_st; en = c->dp_en;
    if (st > e->r->core.pos)
        goto fail;  //not sorted?
    if (e->len > c->inc_sz) {   //grow increment buffer as required
        c->inc = realloc(c->inc, sizeof(int) * e->len);
        c->inc_sz = e->len;
        for(i = 0; i < c->inc_sz; ++i)
            *(c->inc + i) = 1;
    }
    off = e->r->core.pos - st;
    {
        len = off + e->len;
        if (en < (c->w_st+len)) {    //extra space possible, every cigars may not contribute!
            len = c->w_st + len - en;
            if (!(dpth = realloc(c->dpth, (len+c->dp_sz) * sizeof(int)))) {
                goto fail;
            }
            memset(dpth + c->dp_sz, 0, len * sizeof(int));
            c->dp_sz += len;
            c->dpth = dpth;
            c->dp_en = en = c->w_st + c->dp_sz;
        }
        len = 0;
    }
    if (chk) {  //check depth
        for (i = 0; i < e->r->core.n_cigar; ++i) {
            if (!(bam_cigar_type(bam_cigar_op(cgr[i])) & 2)) {   //not consuming ref
                continue;
            }
            //deletion is counted!
            clen = bam_cigar_oplen(cgr[i]);
            for (j = 0; j < clen; ++j) {
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

//1 if required 0 if not -ve on failure
static inline int readrequired(rc_t *c, ce_t *e)
{
    return updatedepth(c, e, 1);
}

//return cache element to cache
void retcache(rc_t *c, ce_t* elem)
{
    if ((elem->r->core.flag & BAM_FPAIRED) && !(elem->r->core.flag & BAM_FUNMAP)) { //paired and mate mapped, remove from expected pair
        khiter_t it = kh_get(kh_pair, c->selpair, bam_get_qname(elem->r));
        if (it != kh_end(c->selpair) && kh_exist(c->selpair, it)) {
            //clear from expected pair hash
            kh_del(kh_pair, c->selpair, it);
        }
    }

#ifdef CACHE_DBG_LOG
    if (elem->ord){//todo dbg
        LG("LG ret %"PRIu64" %s\n", elem->ord, elem->log.l?elem->log.s:"");
    }
#endif //CACHE_DBG_LOG

    elem->prev = NULL;
    elem->ord = 0;
    elem->len = 0;
    elem->next = c->cache.head;
    c->cache.head->prev = elem;
    c->cache.head = elem;
    ++c->cache.f;
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    ks_clear(&elem->log);

}

/// @brief get cache elemnt / storage from preallocated cache
/// @param fp htsFile pointer to setup / retrieve cache
/// @return ce_t* on success or NULL on failure
ce_t* getcache(htsFile *fp)
{
    rc_t *c = (rc_t*)fp->c;
    ce_t *elem = NULL, *tail = NULL, *ret = NULL, **p = NULL;
    int i = 0;
    const int inc = 1024;
    if (!c) { //create cache
        if (!(c = calloc(1, sizeof(rc_t))))
            return NULL;
        fp->c = c;
    }
    if (c->cache.f <= 1) { //grow cache
        assert((c->cache.f == 1 && c->cache.head == c->cache.tail) || (c->cache.f != 1 && c->cache.head != c->cache.tail));
        if (!(p = realloc(c->cache.p, (c->cache.n + 1) * sizeof(ce_t*))))
            return NULL;
        c->cache.p = p;
        if((elem = calloc(inc, sizeof(ce_t)))) {
            c->cache.p[c->cache.n++] = elem;
            if (!c->cache.m){
                c->cache.head = tail = elem;
                if (!(elem->r = bam_init1()))
                    return NULL;
                i = 1;
                ks_initialize(&elem->log);
            } else {
                tail = c->cache.tail;
                i = 0;
            }
            for (; i < inc; ++i) {
                if (!((elem + i )->r = bam_init1()))
                    return NULL;
                tail->next = elem + i;
                tail = tail->next;
                ks_initialize(&elem->log);
            }
            c->cache.m += inc;
            c->cache.f += inc;
            c->cache.tail = tail;
            assert(!c->cache.tail->next);
            assert(c->cache.tail == elem+(i?i-1:0));
        } else
            return NULL;
    }
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    ret = c->cache.head;
    c->cache.head = c->cache.head->next;
    c->cache.head->prev = NULL;
    ret->prev = NULL;
    ret->next = NULL;
    --c->cache.f;
    assert(c->cache.m == c->cache.f+1+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
   return ret;
}

/// @brief add a read to cache
/// @param c pointer to read cache
/// @param e cache element containing the read to be cached
/// @param sts to return status of cache post caching
/// @return -1 on failure 0 on success
int addtoreadcache(rc_t *c, ce_t *e, int *sts)
{
    int unmap = 0;
    if (!(e->r->core.flag & BAM_FUNMAP)) {
        if (c->w_st == -1) {
            c->w_st = c->head ? c->head->r->core.pos : e->r->core.pos;    //starting
            c->w_en = c->w_st + c->wndsz;
            c->dp_en = c->w_st + c->dp_sz;
        }
    } else {
        unmap = 1;
    }
    assert(c->cache.m == c->cache.f+1+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    e->ord = ++(c->ord);
    e->len = bam_cigar2rlen(e->r->core.n_cigar, bam_get_cigar(e->r));
    LG("+ %s %"PRIu64"\t\t%"PRIhts_pos" %"PRIhts_pos" %"PRIu64" %"PRIhts_pos"\n", bam_get_qname(e->r), e->ord, c->w_st, e->r->core.pos, e->len, c->w_en);
    if (!c->head) {
        c->head = c->tail = e;
    } else {
        ce_t *p = c->tail;
        ce_t *tmpn = NULL;
        if (!unmap) {
            if (p->r->core.tid == e->r->core.tid && p->r->core.pos > e->r->core.pos) {
                hts_log_error("Unsorted data");
                return -1;   //not sorted!
            }
            while (p && p->r->core.tid == e->r->core.tid && p->r->core.pos == e->r->core.pos && p->len <= e->len) {
                //consider flags to make sure dup/fail/sec/supp etc. won't mask others
                if (p->r->core.flag>>8 >= e->r->core.flag>>8) {
                    LG("\t%"PRIu64" %d - %"PRIu64" %d, continuing\n", p->ord, p->r->core.flag, e->ord, e->r->core.flag);
                    p = p->prev;
                }
                else {
                    LG("\t%"PRIu64" %d %"PRIu64" %d, break\n", p->ord, p->r->core.flag, e->ord, e->r->core.flag);
                    break;
                }
            }
        }
        if (p) {
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
    ++c->rcnt;
    //todo do we need a limit on max no of items that are cached? like the whole file is for same pos, probably cant be loaded!
    if (c->tid == e->r->core.tid) {
        if (c->w_en < e->r->core.pos) {  //post window, process and advance
            LG("wnd full\n");
            c->trgr = 2; //wnd full, go for processing
        }
        else
            c->trgr = 1;    //caching
    } else if (c->tid != -2) {
        LG("tid change\n");
        c->trgr = 3; //ready for processing
    }
    else
        c->trgr = 1;    //caching

    c->tid = e->r->core.tid;
    if (sts) *sts = c->trgr;
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    LGlog(&e->log, "%s,%"PRIu64",added,%d,%"PRIhts_pos",%"PRIhts_pos",%"PRIu64",%"PRIhts_pos",", bam_get_qname(e->r), e->ord,e->r->core.flag,e->r->core.pos+1, e->r->core.mpos+1,e->len, e->len+e->r->core.pos+1);
    return 0;
}

/// @brief get read from procesed cache
/// @param c pointer to read cache
/// @param b pointer to bam data, to which read data is copied
/// @param end end of read, for iterators
/// @return -1 on failure, 0 when nothing to retrieve and 1 with read retrieved
int getfromreadcache(rc_t *c, bam1_t *b, hts_pos_t *end)
{
    if (!c || c->trgr < 3) {    //either ready or end
        return 0;
    }
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
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

    if (p && (c->trgr == 3 || c->trgr == 4)) { //send only upto start of wnd to maintain the order, except when it is end
        if (!bam_copy1(b, p->r))
            return -1;
        if (p == e) {
            c->selcnt--;
            c->head_sel = p->next;
            if (!c->head_sel) {
                c->tail_sel = NULL;
            }
        }
        else {
            c->head_ins = p->next;
            c->inscnt--;
            if (!c->head_ins) {
                c->tail_ins = NULL;
            }
        }
        if (!c->head_sel && !c->head_ins && c->trgr != 4)
            c->trgr = 0;    //not ready
        else {
            if (c->head_sel) hts_prefetch(c->head_sel);
            if (c->head_ins) hts_prefetch(c->head_ins);
        }
        LG("- %s %"PRIu64"\n", bam_get_qname(b), p->ord);
        LGlog(&p->log, "%s", ",retrieved");
        if (end) *end = p->len + p->r->core.pos;
        retcache(c, p);
        assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);

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
            hts_prefetch(s->next->r);
        if (s->ord > e->ord)
            break;  //not found
        if (s->r->core.pos == e->r->core.mpos &&
            s->r->core.mpos == e->r->core.pos &&
            s->r->core.tid == e->r->core.mtid &&
            s->r->core.mtid == e->r->core.mtid && !strcmp(bam_get_qname(s->r), bam_get_qname(e->r)))
            return s;
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
static inline void moveread(rc_t *c, ce_t *ep, ce_t *e, ce_t* en, int sel, int ins)
{
    int paired = (e->r->core.flag & BAM_FPAIRED) &&
        !(e->r->core.flag & BAM_FUNMAP) && !(e->r->core.flag & BAM_FMUNMAP) &&
        e->r->core.mtid != -1 && e->r->core.mpos != -1;
    if (!ins) { //remove from cache
        assert(c->head == e);   //always be head as they are moved to either selected or nsel
        if (en)
            en->prev = e->prev;
        c->head = en;
        if (!c->head)
            c->tail = c->head;
        c->rcnt--;
    } else { //remove from nsel
        if (ep) {
            ep->next = en;
            if (en)
                en->prev = ep;
        }
        else {
            c->head_nsel = en;
            if (en)
                en->prev = NULL;
        }
        if (!en)
            c->tail_nsel = ep ? ep : NULL;

        c->nselcnt--;
    }
    e->next = NULL;
    e->prev = NULL;

    if (sel) {
        LG("mv %"PRIu64" sel\n", e->ord);
        ins ? c->inscnt++ : c->selcnt++;
        //insert in required pos, starting from tail
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
    } else if (paired) {    //move to nsel if paired, otherwise discard and return cache
        c->nselcnt++;
        LG("mv %"PRIu64" nsel\n", e->ord);
        //add to non-selected list, for pair lookup
        ce_t *s = c->tail_nsel, *p = NULL;
        assert (!s || s->r->core.tid == e->r->core.tid);
        if (s && s->r->core.pos < e->r->core.pos) {   //shortcut
            s->next = e;
            e->next = NULL;
            e->prev = s;
            c->tail_nsel = e;
            return;
        }
        //find pos and fit, in order of increasing pos, that it is easy to remove
        while (s) {
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
        LGlog(&e->log, "%s", "npair,,,disc");
        retcache(c, e);
        return;
    }
}
/// @brief reset cache status, for next tid/iterator...
/// @param c pointer to read cache
static inline void resetdepth(rc_t* c)
{
    c->w_st = -1;
    if (c->dp_sz <= 0 || !c->dpth)
        return;
    LG("reset\n");
    memset(c->dpth, 0, c->dp_sz * sizeof(int));
    ce_t *en = NULL;
    //clear all from previous tid
    while (c->head && c->head->r->core.tid != c->tid) {
        en = c->head->next;
        c->rcnt--;
        retcache(c, c->head);
        c->head = en;
    }
    if (!c->head) c->tail = NULL;
    else c->head->prev = NULL;

    //clear whole non selected ones
    while (c->head_nsel) {
        en = c->head_nsel->next;
        c->nselcnt--;
        retcache(c, c->head_nsel);
        c->head_nsel = en;
    }
    assert(!c->head_nsel && !c->nselcnt);
    if(c->head_sel) assert(c->head_sel->r->core.tid == c->tail_sel->r->core.tid);
    if(c->head_nsel) assert(c->head_nsel->r->core.tid == c->tail_nsel->r->core.tid);
    if(c->head_ins) assert(c->head_ins->r->core.tid == c->tail_ins->r->core.tid);

    LG("reset: t %"PRIu64" s %"PRIu64" i %"PRIu64" n %"PRIu64"; nxt %d\n", c->rcnt, c->selcnt,c->inscnt, c->nselcnt, c->tid);
    c->tail_nsel = NULL;
}

/// @brief process the cache and find reads relevant based on depth
/// @param c pointer to read cache
/// @return -ve on error, 0 on success
int processcache(rc_t *c)
{
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    assert(!c->inscnt && !c->selcnt);
    assert(!c->head_sel && !c->tail_sel);
    assert(!c->head_ins && !c->tail_ins);

    /* there is a window which starts at pos of 1st read and adds reads upto
    a read past the end pos to cache. a read past the end, a tid change or end
    of file could be a trigger to start processing the cached reads. they are
    expected as sorted by pos and kept in descending order of reference
    consumption (given by bam_cigar2rlen) if have same pos.

    iterate through cached reads, check if it is expected as pair of an
    earlier selected read otherwise use cigar to calc depth and decide whether
    read is needed or not. move required reads to sel list and unwanted paired
    ones to nsel list. check whether selected read's pair is awaited or already
    processed based on pos values. add to selected pair hash for awaiting ones
    and search in nsel list for already processed ones. if a pair is found, move
    to ins list. sel and ins list are in order of reading / as in source.

    processing will stop at tid change or once whole cache is processed. window
    will be adjusted to include the read found outside the end pos. start will
    move by same amount. if there is gap b/w start and next read, start moves to
    that pos and end gets extended. any read which falls outside the new start
    in nsel list is discarded.

    the 1st one from sel / ins list is removed and passed until both are empty.
    the read and processig continues with new window extremities.
    */
    ce_t *e = c->head, *en = NULL;
    pair_exp *p = NULL;
    khiter_t it;
    int sel = 0, ret = -1;
    int chkpair, foundpair;
    LG("pc: t %"PRIu64" s %"PRIu64" i %"PRIu64" n %"PRIu64"\n", c->rcnt, c->selcnt,c->inscnt, c->nselcnt);
    while (e) {
        sel = foundpair = 0;
        chkpair = (e->r->core.flag & BAM_FPAIRED) && !(e->r->core.flag & BAM_FMUNMAP);
        en = e->next;
        if (c->tid == e->r->core.tid) {
            if (c->trgr != 2) { //if not wnd full, it is either end or tid change
                resetdepth(c);
                ret = 0;
                break;   //last one / one that triggered the processing; on next iteration
            } else {    //wnd full, process and move wnd
                if (e->r->core.pos >= c->w_en) { //done enough
                    LG("* wnd full,[%"PRIhts_pos" - %"PRIhts_pos"] processed upto %"PRIhts_pos"\n", c->w_st, c->w_en, e->r->core.pos);
                    break;
                }
            }
        }
        LG("* checking %s %"PRIu64"\n", bam_get_qname(e->r), e->ord);
        if (e->r->core.flag & BAM_FUNMAP) {//unmapped, select anyway
            //selectread(c, ep, e, en, 0);   //add to selected list
            moveread(c, NULL, e, en, 1, 0);   //add to selected list
            LG("* s unmap %s %"PRIu64"\n", bam_get_qname(e->r), e->ord);
            LGlog(&e->log, "%s", "sel,umap,,");
            e = en;
            continue;
        }
        //LGlog(&e->log,",");
        if (chkpair) {    //paired and mate mapped
            //have to remove from map as ce_t are freed; also they can't be modified while in cache
            if ((it = kh_get(kh_pair, c->selpair, bam_get_qname(e->r))) != kh_end(c->selpair)) {
                if (kh_exist(c->selpair, it)) { //iterate and find pair to this
                    p = &kh_val(c->selpair, it);
                    if (p->mpos == e->r->core.pos &&
                        p->mtid == e->r->core.tid &&
                        p->pos == e->r->core.mpos &&
                        p->tid == e->r->core.mtid) {   //pair already selected
                        kh_del(kh_pair, c->selpair, it);    //remove from expected pairs
                        foundpair = 1;
                        moveread(c, NULL, e, en, 1, 0);   //select this
                        sel = 1;
                    } else {
                        kh_del(kh_pair, c->selpair, it);    //remove from expected pairs
                        //not possible to have duplicate on qname, chk n confirm
                        it = kh_end(c->selpair);
                    }
                }
            }
        }
        if (!sel) {
            //check depth
            int r = 0;
            if ((r = readrequired(c, e)) > 0) {  //read required
                moveread(c, NULL, e, en, 1, 0);   //select this
                sel = 1;
            } else if (r < 0) {
                goto fail;
            }
        }
        if (sel) {
            if (updatedepth(c, e, 0) == -1)
                goto fail;
            LG("* s %s %"PRIu64" wnd:%"PRIhts_pos"-%"PRIhts_pos"", bam_get_qname(e->r), e->ord, c->w_st, c->w_en);
            LGlog(&e->log, "%s", "sel,");
            if (chkpair && !foundpair) { //1st one or pair not selected
                if (e->r->core.pos <= e->r->core.mpos) {    //add only if it is yet to be processed, sorted data!
                    int r = 0;
                    it = kh_put(kh_pair, c->selpair, bam_get_qname(e->r), &r);
                    if (r == -1)
                        goto fail;
                    pair_exp *p = &kh_val(c->selpair, it);
                    p->pos = e->r->core.pos; p->tid = e->r->core.tid;
                    p->mpos = e->r->core.mpos; p->mtid = e->r->core.mtid;
                    LG(" PAIR expected");
                    LGlog(&e->log, "%s", "paired,,");
                } else {
                    //do it after finishing the loop, to avoid issues with ep/epp...
                    //have to insert them based on ord., if not found, discard. if eq. limit there if done here.
                    ce_t *o = NULL, *op = NULL;
                    if ((o = find_nsel(c, e, &op))) {
                        if (o->r->core.pos >= c->w_st) {    //only if order can be maintained
                            moveread(c, op, o, o->next, 1, 1);
                            if (updatedepth(c, o, 0) == -1)
                                goto fail;
                            LG(" inserted PAIR\n* s %s %"PRIu64" (inspair)", bam_get_qname(o->r), o->ord);
                            LGlog(&o->log, "%s", "paired,inserted,");
                        }
                        LGlog(&e->log, "%s", "paired,nsel,");
                    } else {
                        LG(" no PAIR");
                        LGlog(&e->log, "%s", "paired,notfound,");
                    }
                }
            } else if (foundpair) {
                LG(" found PAIR");
                LGlog(&e->log, "%s", "paired,found,");
            } else {
                LG(" no PAIR");
                LGlog(&e->log, "%s", "notpaired,NA,");
            }
            LG("\n");
        }
        else {
            LG("* d %s %"PRIu64"\n", bam_get_qname(e->r), e->ord);
            LGlog(&e->log, "%s", "nsel,");
            moveread(c, NULL, e, en, 0, 0);   //remove as non-selected
        }
        e = en; //chk with next one
    }
    if (c->trgr == 2) {    //2 --> wnd full, processed, move wnd
        en = NULL;
        hts_pos_t adj = c->tail ? c->tail->r->core.pos - c->w_en : 0; //last one, out of window - current end
        hts_pos_t new_st = c->w_st + adj;
        hts_pos_t bkp_st = c->w_st;
        int rem = 0;
        while (c->head_nsel && c->head_nsel->r->core.pos < new_st) { //holding until wnd passes mate pos, but anything after this which has already passed out is held until this is cleared!
            rem = 1;
            en = c->head_nsel->next;
            LG("* nsel discarded %s %"PRIu64"\n", bam_get_qname(c->head_nsel->r), c->head_nsel->ord);
            LGlog(&c->head_nsel->log,"%s",",,,nsel-disc,");
            c->nselcnt--;
            retcache(c, c->head_nsel);
            if(!(c->head_nsel = en)) c->tail_nsel = NULL;
        }
        if (rem) {
          LG("* wnd full, removed items from head_nsel\n")//fprintf(fp1, "* wnd full, removed items from head_nsel\n");
        }
        else {
           LG("* wnd full, 0 removed items from head_nsel, [%"PRIhts_pos"-%"PRIhts_pos"] %"PRIhts_pos"\n", c->w_st, c->w_en, c->head_nsel?c->head_nsel->r->core.pos : 0)//fprintf(fp1, "* wnd full, removed items from head_nsel\n")
        }
        assert(e == c->head);
        c->w_st = c->head ? c->head->r->core.pos : new_st;    //move wnd
        c->w_en = c->w_st + c->wndsz;
        adj = c->w_st - bkp_st;
        if (adj >= c->dp_sz) {
            memset(c->dpth, 0, c->dp_sz * sizeof(int));
            c->dp_en = c->w_st + c->dp_sz;
            LG("0 dpth buffer\n");
        } else {
            LG("adj %"PRIhts_pos", mv %"PRIhts_pos"-%"PRIhts_pos",", adj, c->w_st+adj, c->w_st+c->dp_sz);
            LG("0 set %"PRIhts_pos" - %"PRIhts_pos"\n", c->w_st+c->dp_sz-adj,c->w_st+c->dp_sz);
            memmove(c->dpth, c->dpth + adj, (c->dp_sz - adj) * sizeof(int));
            memset(c->dpth + c->dp_sz - adj, 0, adj * sizeof(int));
            c->dp_en += adj;
            assert(c->dp_en == (c->w_st+c->dp_sz));
        }
        LG("* wnd moved, %"PRIhts_pos" - %"PRIhts_pos", dpth %"PRIhts_pos" - %"PRIhts_pos"; s %"PRIu64" i %"PRIu64" ns %"PRIu64"\n", c->w_st, c->w_en, c->w_st, c->dp_en, c->selcnt, c->inscnt, c->nselcnt);
        c->trgr = 3;    //reset full status n get already processedn
    }
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    LG("pc2: t %"PRIu64" s %"PRIu64" i %"PRIu64" n %"PRIu64"\n", c->rcnt, c->selcnt,c->inscnt, c->nselcnt);

    if (c->head_sel)
        hts_prefetch(c->head_sel);
    if (c->head_ins)
        hts_prefetch(c->head_ins);
    return ret;
fail:
    LG(" FAIL\n");
    return -1;
}

//wrappers / for iterators
void* getsamcache(hts_itr_t *itr, void *data)
{
    htsFile *fp = NULL;
    if (itr && itr->usecache)
        fp = (htsFile*)data;
    return fp ? fp->c : NULL;
}

int getfromreadcache_iter(void *p, void *s, int *tid, hts_pos_t *beg, hts_pos_t* end) //?
{
    rc_t *c = (rc_t*)p;
    bam1_t *b = (bam1_t*)s;
    int ret = getfromreadcache(c, b, end);
    if (ret > 0) {
        *tid = b->core.tid;
        *beg = b->core.pos;
    }

    if (!ret) {
        assert(!c->inscnt && !c->selcnt);
        assert(!c->head_sel && !c->tail_sel);
        assert(!c->head_ins && !c->tail_ins);
    }

    return ret;
}

void *getcache_iter(void *data) //?
{
    htsFile *fp = (htsFile*)data;
    rc_t *c = (rc_t*)fp->c;
    assert(c->cache.m == c->cache.f+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    void *p = getcache(fp);
    assert(c->cache.m == c->cache.f+1+c->rcnt+c->inscnt+c->selcnt+c->nselcnt);
    return p;
}

void *getreadbuffer_iter(void *p)
{
    ce_t* e = (ce_t*)p;
    return e->r;
}

void notifyend_iter(void *p, void *e)
{
    rc_t *c = (rc_t*)p;
    c->trgr = 4;//end
    c->tid = -2;//reset that it doesn't match to any
    retcache(c, (ce_t*)e);
}

int addtoreadcache_iter(void *c, void *s, int *sts)
{
    return addtoreadcache(c, s, sts);
}

int processcache_iter(void *c)
{
    return processcache((rc_t*)c);
}

void resetcache_itr(rc_t *c)
{
    resetdepth(c);
    c->trgr = 0;
    c->tid = -2;
    c->w_st = c->w_en = -1;
}
