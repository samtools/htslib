/*
Copyright (c) 2013-2020, 2023-2024 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * The index is a gzipped tab-delimited text file with one line per slice.
 * The columns are:
 * 1: reference number (0 to N-1, as per BAM ref_id)
 * 2: reference position of 1st read in slice (1..?)
 * 3: number of reads in slice
 * 4: offset of container start (relative to end of SAM header, so 1st
 *    container is offset 0).
 * 5: slice number within container (ie which landmark).
 *
 * In memory, we hold this in a nested containment list. Each list element is
 * a cram_index struct. Each element in turn can contain its own list of
 * cram_index structs.
 *
 * Any start..end range which is entirely contained within another (and
 * earlier as it is sorted) range will be held within it. This ensures that
 * the outer list will never have containments and we can safely do a
 * binary search to find the first range which overlaps any given coordinate.
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "../htslib/bgzf.h"
#include "../htslib/hfile.h"
#include "../hts_internal.h"
#include "cram.h"
#include "os.h"

#if 0
static void dump_index_(cram_index *e, int level) {
    int i, n;
    n = printf("%*s%d / %d .. %d, ", level*4, "", e->refid, e->start, e->end);
    printf("%*soffset %"PRId64" %p %p\n", MAX(0,50-n), "", e->offset, e, e->e_next);
    for (i = 0; i < e->nslice; i++) {
        dump_index_(&e->e[i], level+1);
    }
}

static void dump_index(cram_fd *fd) {
    int i;
    for (i = 0; i < fd->index_sz; i++) {
        dump_index_(&fd->index[i], 0);
    }
}
#endif

// Thread a linked list through the nested containment list.
// This makes navigating it and finding the "next" index entry
// trivial.
static cram_index *link_index_(cram_index *e, cram_index *e_last) {
    int i;
    if (e_last)
        e_last->e_next = e;

    // We don't want to link in the top-level cram_index with
    // offset=0 and start/end = INT_MIN/INT_MAX.
    if (e->offset)
        e_last = e;

    for (i = 0; i < e->nslice; i++)
        e_last = link_index_(&e->e[i], e_last);

    return e_last;
}

static void link_index(cram_fd *fd) {
    int i;
    cram_index *e_last = NULL;

    for (i = 0; i < fd->index_sz; i++) {
        e_last = link_index_(&fd->index[i], e_last);
    }

    if (e_last)
        e_last->e_next = NULL;
}

static int kget_int32(kstring_t *k, size_t *pos, int32_t *val_p) {
    int sign = 1;
    int32_t val = 0;
    size_t p = *pos;

    while (p < k->l && (k->s[p] == ' ' || k->s[p] == '\t'))
        p++;

    if (p < k->l && k->s[p] == '-')
        sign = -1, p++;

    if (p >= k->l || !(k->s[p] >= '0' && k->s[p] <= '9'))
        return -1;

    while (p < k->l && k->s[p] >= '0' && k->s[p] <= '9') {
        int digit = k->s[p++]-'0';
        val = val*10 + digit;
    }

    *pos = p;
    *val_p = sign*val;

    return 0;
}

static int kget_int64(kstring_t *k, size_t *pos, int64_t *val_p) {
    int sign = 1;
    int64_t val = 0;
    size_t p = *pos;

    while (p < k->l && (k->s[p] == ' ' || k->s[p] == '\t'))
        p++;

    if (p < k->l && k->s[p] == '-')
        sign = -1, p++;

    if (p >= k->l || !(k->s[p] >= '0' && k->s[p] <= '9'))
        return -1;

    while (p < k->l && k->s[p] >= '0' && k->s[p] <= '9') {
        int digit = k->s[p++]-'0';
        val = val*10 + digit;
    }

    *pos = p;
    *val_p = sign*val;

    return 0;
}

/*
 * Loads a CRAM .crai index into memory.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, const char *fn, const char *fn_idx) {

    char *tfn_idx = NULL;
    char buf[65536];
    ssize_t len;
    kstring_t kstr = {0};
    hFILE *fp;
    cram_index *idx;
    cram_index **idx_stack = NULL, *ep, e;
    int idx_stack_alloc = 0, idx_stack_ptr = 0;
    size_t pos = 0;

    /* Check if already loaded */
    if (fd->index)
        return 0;

    fd->index = calloc((fd->index_sz = 1), sizeof(*fd->index));
    if (!fd->index)
        return -1;

    idx = &fd->index[0];
    idx->refid = -1;
    idx->start = INT_MIN;
    idx->end   = INT_MAX;

    idx_stack = calloc(++idx_stack_alloc, sizeof(*idx_stack));
    if (!idx_stack)
        goto fail;

    idx_stack[idx_stack_ptr] = idx;

    // Support pathX.cram##idx##pathY.crai
    const char *fn_delim = strstr(fn, HTS_IDX_DELIM);
    if (fn_delim && !fn_idx)
        fn_idx = fn_delim + strlen(HTS_IDX_DELIM);

    if (!fn_idx) {
        if (hts_idx_check_local(fn, HTS_FMT_CRAI, &tfn_idx) == 0 && hisremote(fn))
            tfn_idx = hts_idx_getfn(fn, ".crai");

        if (!tfn_idx) {
            hts_log_error("Could not retrieve index file for '%s'", fn);
            goto fail;
        }
        fn_idx = tfn_idx;
    }

    if (!(fp = hopen(fn_idx, "r"))) {
        hts_log_error("Could not open index file '%s'", fn_idx);
        goto fail;
    }

    // Load the file into memory
    while ((len = hread(fp, buf, sizeof(buf))) > 0) {
        if (kputsn(buf, len, &kstr) < 0)
            goto fail;
    }

    if (len < 0 || kstr.l < 2)
        goto fail;

    if (hclose(fp) < 0)
        goto fail;

    // Uncompress if required
    if (kstr.s[0] == 31 && (uc)kstr.s[1] == 139) {
        size_t l = 0;
        char *s = zlib_mem_inflate(kstr.s, kstr.l, &l);
        if (!s)
            goto fail;

        free(kstr.s);
        kstr.s = s;
        kstr.l = l;
        kstr.m = l; // conservative estimate of the size allocated
        if (kputsn("", 0, &kstr) < 0) // ensure kstr.s is NUL-terminated
            goto fail;
    }


    // Parse it line at a time
    while (pos < kstr.l) {
        /* 1.1 layout */
        if (kget_int32(&kstr, &pos, &e.refid) == -1)
            goto fail;

        if (kget_int32(&kstr, &pos, &e.start) == -1)
            goto fail;

        if (kget_int32(&kstr, &pos, &e.end) == -1)
            goto fail;

        if (kget_int64(&kstr, &pos, &e.offset) == -1)
            goto fail;

        if (kget_int32(&kstr, &pos, &e.slice) == -1)
            goto fail;

        if (kget_int32(&kstr, &pos, &e.len) == -1)
            goto fail;

        e.end += e.start-1;
        //printf("%d/%d..%d-offset=%" PRIu64 ",len=%d,slice=%d\n", e.refid, e.start, e.end, e.offset, e.len, e.slice);

        if (e.refid < -1) {
            hts_log_error("Malformed index file, refid %d", e.refid);
            goto fail;
        }

        if (e.refid != idx->refid) {
            if (fd->index_sz < e.refid+2) {
                cram_index *new_idx;
                int new_sz = e.refid+2;
                size_t index_end = fd->index_sz * sizeof(*fd->index);
                new_idx = realloc(fd->index,
                                  new_sz * sizeof(*fd->index));
                if (!new_idx)
                    goto fail;

                fd->index = new_idx;
                fd->index_sz = new_sz;
                memset(((char *)fd->index) + index_end, 0,
                       fd->index_sz * sizeof(*fd->index) - index_end);
            }
            idx = &fd->index[e.refid+1];
            idx->refid = e.refid;
            idx->start = INT_MIN;
            idx->end   = INT_MAX;
            idx->nslice = idx->nalloc = 0;
            idx->e = NULL;
            idx_stack[(idx_stack_ptr = 0)] = idx;
        }

        while (!(e.start >= idx->start && e.end <= idx->end) ||
               (idx->start == 0 && idx->refid == -1)) {
            idx = idx_stack[--idx_stack_ptr];
        }

        // Now contains, so append
        if (idx->nslice+1 >= idx->nalloc) {
            cram_index *new_e;
            idx->nalloc = idx->nalloc ? idx->nalloc*2 : 16;
            new_e = realloc(idx->e, idx->nalloc * sizeof(*idx->e));
            if (!new_e)
                goto fail;

            idx->e = new_e;
        }

        e.nalloc = e.nslice = 0; e.e = NULL;
        *(ep = &idx->e[idx->nslice++]) = e;
        idx = ep;

        if (++idx_stack_ptr >= idx_stack_alloc) {
            cram_index **new_stack;
            idx_stack_alloc *= 2;
            new_stack = realloc(idx_stack, idx_stack_alloc*sizeof(*idx_stack));
            if (!new_stack)
                goto fail;
            idx_stack = new_stack;
        }
        idx_stack[idx_stack_ptr] = idx;

        while (pos < kstr.l && kstr.s[pos] != '\n')
            pos++;
        pos++;
    }

    free(idx_stack);
    free(kstr.s);
    free(tfn_idx);

    // Convert NCList to linear linked list
    link_index(fd);

    //dump_index(fd);

    return 0;

 fail:
    free(kstr.s);
    free(idx_stack);
    free(tfn_idx);
    cram_index_free(fd); // Also sets fd->index = NULL
    return -1;
}

static void cram_index_free_recurse(cram_index *e) {
    if (e->e) {
        int i;
        for (i = 0; i < e->nslice; i++) {
            cram_index_free_recurse(&e->e[i]);
        }
        free(e->e);
    }
}

void cram_index_free(cram_fd *fd) {
    int i;

    if (!fd->index)
        return;

    for (i = 0; i < fd->index_sz; i++) {
        cram_index_free_recurse(&fd->index[i]);
    }
    free(fd->index);

    fd->index = NULL;
}

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position, or one immediately preceding it if none is found in
 * the index to overlap this position. (Our index may have missing
 * entries, but we require at least one per reference.)
 *
 * If the index finds multiple slices overlapping this position we
 * return the first one only. Subsequent calls should specify
 * "from" as the last slice we checked to find the next one. Otherwise
 * set "from" to be NULL to find the first one.
 *
 * Refid can also be any of the special HTS_IDX_ values.
 * For backwards compatibility, refid -1 is equivalent to HTS_IDX_NOCOOR.
 *
 * Returns the cram_index pointer on success
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, hts_pos_t pos,
                             cram_index *from) {
    int i, j, k;
    cram_index *e;

    if (from) {
        // Continue from a previous search.
        // We switch to just scanning the linked list, as the nested
        // lists are typically short.
        if (refid == HTS_IDX_NOCOOR)
            refid = -1;

        e = from->e_next;
        if (e && e->refid == refid && e->start <= pos)
            return e;
        else
            return NULL;
    }

    switch(refid) {
    case HTS_IDX_NONE:
    case HTS_IDX_REST:
        // fail, or already there, dealt with elsewhere.
        return NULL;

    case -1:
    case HTS_IDX_NOCOOR:
        refid = -1;
        pos = 0;
        break;

    case HTS_IDX_START: {
        int64_t min_idx = INT64_MAX;
        for (i = 0, j = -1; i < fd->index_sz; i++) {
            if (fd->index[i].e && fd->index[i].e[0].offset < min_idx) {
                min_idx = fd->index[i].e[0].offset;
                j = i;
            }
        }
        if (j < 0)
            return NULL;
        return fd->index[j].e;
    }

    default:
        if (refid < HTS_IDX_NONE || refid+1 >= fd->index_sz)
            return NULL;
    }

    from = &fd->index[refid+1];

    // Ref with nothing aligned against it.
    if (!from->e)
        return NULL;

    // This sequence is covered by the index, so binary search to find
    // the optimal starting block.
    i = 0, j = fd->index[refid+1].nslice-1;
    for (k = j/2; k != i; k = (j-i)/2 + i) {
        if (from->e[k].refid > refid) {
            j = k;
            continue;
        }

        if (from->e[k].refid < refid) {
            i = k;
            continue;
        }

        if (from->e[k].start >= pos) {
            j = k;
            continue;
        }

        if (from->e[k].start < pos) {
            i = k;
            continue;
        }
    }
    // i==j or i==j-1. Check if j is better.
    if (j >= 0 && from->e[j].start < pos && from->e[j].refid == refid)
        i = j;

    /* The above found *a* bin overlapping, but not necessarily the first */
    while (i > 0 && from->e[i-1].end >= pos)
        i--;

    /* We may be one bin before the optimum, so check */
    while (i+1 < from->nslice &&
           (from->e[i].refid < refid ||
            from->e[i].end < pos))
        i++;

    e = &from->e[i];

    return e;
}

// Return the index entry for last slice on a specific reference.
cram_index *cram_index_last(cram_fd *fd, int refid, cram_index *from) {
    int slice;

    if (refid+1 < 0 || refid+1 >= fd->index_sz)
        return NULL;

    if (!from)
        from = &fd->index[refid+1];

    // Ref with nothing aligned against it.
    if (!from->e)
        return NULL;

    slice = fd->index[refid+1].nslice - 1;

    // e is the last entry in the nested containment list, but it may
    // contain further slices within it.
    cram_index *e = &from->e[slice];
    while (e->e_next)
        e = e->e_next;

    return e;
}

/*
 * Find the last container overlapping pos 'end', and the file offset of
 * its end (equivalent to the start offset of the container following it).
 */
cram_index *cram_index_query_last(cram_fd *fd, int refid, hts_pos_t end) {
    cram_index *e = NULL, *prev_e;
    do {
        prev_e = e;
        e = cram_index_query(fd, refid, end, prev_e);
    } while (e);

    if (!prev_e)
        return NULL;
    e = prev_e;

    // Note: offset of e and e->e_next may be the same if we're using a
    // multi-ref container where a single container generates multiple
    // index entries.
    //
    // We need to keep iterating until offset differs in order to find
    // the genuine file offset for the end of container.
    do {
        prev_e = e;
        e = e->e_next;
    } while (e && e->offset == prev_e->offset);

    return prev_e;
}

/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * In theory we call cram_index_query multiple times, once per slice
 * overlapping the range. However slices may be absent from the index
 * which makes this problematic. Instead we find the left-most slice
 * and then read from then on, skipping decoding of slices and/or
 * whole containers when they don't overlap the specified cram_range.
 *
 * This function also updates the cram_fd range field.
 *
 * Returns 0 on success
 *        -1 on general failure
 *        -2 on no-data (empty chromosome)
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r) {
    int ret = 0;
    cram_index *e;

    if (r->refid == HTS_IDX_NONE) {
        ret = -2; goto err;
    }

    // Ideally use an index, so see if we have one.
    if ((e = cram_index_query(fd, r->refid, r->start, NULL))) {
        if (0 != cram_seek(fd, e->offset, SEEK_SET)) {
            if (0 != cram_seek(fd, e->offset - fd->first_container, SEEK_CUR)) {
                ret = -1; goto err;
            }
        }
    } else {
        // Absent from index, but this most likely means it simply has no data.
        ret = -2; goto err;
    }

    pthread_mutex_lock(&fd->range_lock);
    fd->range = *r;
    if (r->refid == HTS_IDX_NOCOOR) {
        fd->range.refid = -1;
        fd->range.start = 0;
    } else if (r->refid == HTS_IDX_START || r->refid == HTS_IDX_REST) {
        fd->range.refid = -2; // special case in cram_next_slice
    }
    pthread_mutex_unlock(&fd->range_lock);

    if (fd->ctr) {
        cram_free_container(fd->ctr);
        if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
            cram_free_container(fd->ctr_mt);
        fd->ctr = NULL;
        fd->ctr_mt = NULL;
        fd->ooc = 0;
        fd->eof = 0;
    }

    return 0;

 err:
    // It's unlikely fd->range will be accessed after EOF or error,
    // but this maintains identical behaviour to the previous code.
    pthread_mutex_lock(&fd->range_lock);
    fd->range = *r;
    pthread_mutex_unlock(&fd->range_lock);
    return ret;
}


/*
 * A specialised form of cram_index_build (below) that deals with slices
 * having multiple references in this (ref_id -2). In this scenario we
 * decode the slice to look at the RI data series instead.
 *
 * Returns 0 on success
 *        -1 on read failure
 *        -2 on wrong sort order
 *        -4 on write failure
 */
static int cram_index_build_multiref(cram_fd *fd,
                                     cram_container *c,
                                     cram_slice *s,
                                     BGZF *fp,
                                     off_t cpos,
                                     int32_t landmark,
                                     int sz) {
    int i, ref = -2;
    int64_t ref_start = 0, ref_end;
    char buf[1024];

    if (fd->mode != 'w') {
        if (0 != cram_decode_slice(fd, c, s, fd->header))
            return -1;
    }

    ref_end = INT_MIN;

    int32_t last_ref = -9;
    int32_t last_pos = -9;
    for (i = 0; i < s->hdr->num_records; i++) {
        if (s->crecs[i].ref_id == last_ref && s->crecs[i].apos < last_pos) {
            hts_log_error("CRAM file is not sorted by chromosome / position");
            return -2;
        }
        last_ref = s->crecs[i].ref_id;
        last_pos = s->crecs[i].apos;

        if (s->crecs[i].ref_id == ref) {
            if (ref_end < s->crecs[i].aend)
                ref_end = s->crecs[i].aend;
            continue;
        }

        if (ref != -2) {
            snprintf(buf, sizeof(buf),
                     "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
                     ref, ref_start, ref_end - ref_start + 1,
                     (int64_t)cpos, landmark, sz);
            if (bgzf_write(fp, buf, strlen(buf)) < 0)
                return -4;
        }

        ref = s->crecs[i].ref_id;
        ref_start = s->crecs[i].apos;
        ref_end   = s->crecs[i].aend;
    }

    if (ref != -2) {
        snprintf(buf, sizeof(buf),
                 "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
                 ref, ref_start, ref_end - ref_start + 1,
                 (int64_t)cpos, landmark, sz);
        if (bgzf_write(fp, buf, strlen(buf)) < 0)
            return -4;
    }

    return 0;
}

/*
 * Adds a single slice to the index.
 */
int cram_index_slice(cram_fd *fd,
                     cram_container *c,
                     cram_slice *s,
                     BGZF *fp,
                     off_t cpos,
                     off_t spos, // relative to cpos
                     off_t sz) {
    int ret;
    char buf[1024];

    if (sz > INT_MAX) {
        hts_log_error("CRAM slice is too big (%"PRId64" bytes)",
                      (int64_t) sz);
        return -1;
    }

    if (s->hdr->ref_seq_id == -2) {
        ret = cram_index_build_multiref(fd, c, s, fp, cpos, spos, sz);
    } else {
        snprintf(buf, sizeof(buf),
                 "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
                 s->hdr->ref_seq_id, s->hdr->ref_seq_start,
                 s->hdr->ref_seq_span, (int64_t)cpos, (int)spos, (int)sz);
        ret = (bgzf_write(fp, buf, strlen(buf)) >= 0)? 0 : -4;
    }

    return ret;
}

/*
 * Adds a single container to the index.
 */
static
int cram_index_container(cram_fd *fd,
                         cram_container *c,
                         BGZF *fp,
                         off_t cpos) {
    int j;
    off_t spos;

    // 2.0 format
    for (j = 0; j < c->num_landmarks; j++) {
        cram_slice *s;
        off_t sz;
        int ret;

        spos = htell(fd->fp);
        if (spos - cpos - (off_t) c->offset != c->landmark[j]) {
            hts_log_error("CRAM slice offset %"PRId64" does not match"
                          " landmark %d in container header (%"PRId32")",
                          (int64_t) (spos - cpos - (off_t) c->offset),
                          j, c->landmark[j]);
            return -1;
        }

        if (!(s = cram_read_slice(fd))) {
            return -1;
        }

        sz = htell(fd->fp) - spos;
        ret = cram_index_slice(fd, c, s, fp, cpos, c->landmark[j], sz);

        cram_free_slice(s);

        if (ret < 0) {
            return ret;
        }
    }

    return 0;
}


/*
 * Builds an index file.
 *
 * fd is a newly opened cram file that we wish to index.
 * fn_base is the filename of the associated CRAM file.
 * fn_idx is the filename of the index file to be written;
 * if NULL, we add ".crai" to fn_base to get the index filename.
 *
 * Returns 0 on success,
 *         negative on failure (-1 for read failure, -4 for write failure)
 */
int cram_index_build(cram_fd *fd, const char *fn_base, const char *fn_idx) {
    cram_container *c;
    off_t cpos, hpos;
    BGZF *fp;
    kstring_t fn_idx_str = {0};
    int64_t last_ref = -9, last_start = -9;

    // Useful for cram_index_build_multiref
    cram_set_option(fd, CRAM_OPT_REQUIRED_FIELDS, SAM_RNAME | SAM_POS | SAM_CIGAR);

    if (! fn_idx) {
        kputs(fn_base, &fn_idx_str);
        kputs(".crai", &fn_idx_str);
        fn_idx = fn_idx_str.s;
    }

    if (!(fp = bgzf_open(fn_idx, "wg"))) {
        perror(fn_idx);
        free(fn_idx_str.s);
        return -4;
    }

    free(fn_idx_str.s);

    cpos = htell(fd->fp);
    while ((c = cram_read_container(fd))) {
        if (fd->err) {
            perror("Cram container read");
            return -1;
        }

        hpos = htell(fd->fp);

        if (!(c->comp_hdr_block = cram_read_block(fd)))
            return -1;
        assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

        c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
        if (!c->comp_hdr)
            return -1;

        if (c->ref_seq_id == last_ref && c->ref_seq_start < last_start) {
            hts_log_error("CRAM file is not sorted by chromosome / position");
            return -2;
        }
        last_ref = c->ref_seq_id;
        last_start = c->ref_seq_start;

        if (cram_index_container(fd, c, fp, cpos) < 0) {
            bgzf_close(fp);
            return -1;
        }

        off_t next_cpos = htell(fd->fp);
        if (next_cpos != hpos + c->length) {
            hts_log_error("Length %"PRId32" in container header at offset %lld does not match block lengths (%lld)",
                          c->length, (long long) cpos, (long long) next_cpos - hpos);
            return -1;
        }
        cpos = next_cpos;

        cram_free_container(c);
    }
    if (fd->err) {
        bgzf_close(fp);
        return -1;
    }

    return (bgzf_close(fp) >= 0)? 0 : -4;
}

// internal recursive step
static int64_t cram_num_containers_between_(cram_index *e, int64_t *last_pos,
                                            int64_t nct,
                                            off_t cstart, off_t cend,
                                            int64_t *first, int64_t *last) {
    int64_t nc = 0, i;

    if (e->offset) {
        if (e->offset != *last_pos) {
            if (e->offset >= cstart && (!cend || e->offset <= cend)) {
                if (first && *first < 0)
                    *first = nct;
                if (last)
                    *last = nct;
            }
            nc++;
        }
        // else a new multi-ref in same container
        *last_pos = e->offset;
    }

    for (i = 0; i < e->nslice; i++)
        nc += cram_num_containers_between_(&e->e[i], last_pos, nc + nct,
                                           cstart, cend, first, last);

    return nc;
}

/*! Returns the number of containers in the CRAM file within given offsets.
 *
 * The cstart and cend offsets are the locations of the start of containers
 * as returned by index_container_offset.
 *
 * If non-NULL, first and last will hold the inclusive range of container
 * numbers, counting from zero.
 *
 * @return
 * Returns the number of containers, equivalent to *last-*first+1.
 */
int64_t cram_num_containers_between(cram_fd *fd,
                                    off_t cstart, off_t cend,
                                    int64_t *first, int64_t *last) {
    int64_t nc = 0, i;
    int64_t last_pos = -99;
    int64_t l_first = -1, l_last = -1;

    for (i = 0; i < fd->index_sz; i++) {
        int j = i+1 == fd->index_sz ? 0 : i+1; // maps "*" to end
        nc += cram_num_containers_between_(&fd->index[j], &last_pos, nc,
                                           cstart, cend, &l_first, &l_last);
    }

    if (first)
        *first = l_first;
    if (last)
        *last = l_last;

    return l_last - l_first + 1;
}

/*
 * Queries the total number of distinct containers in the index.
 * Note there may be more containers in the file than in the index, as we
 * are not required to have an index entry for every one.
 */
int64_t cram_num_containers(cram_fd *fd) {
    return cram_num_containers_between(fd, 0, 0, NULL, NULL);
}


/*! Returns the byte offset for the start of the n^th container.
 *
 * The index must have previously been loaded, otherwise <0 is returned.
 */
static cram_index *cram_container_num2offset_(cram_index *e, int num,
                                              int64_t *last_pos, int *nc) {
    if (e->offset) {
        if (e->offset != *last_pos) {
            if (*nc == num)
                return e;
            (*nc)++;
        }
        // else a new multi-ref in same container
        *last_pos = e->offset;
    }

    int i;
    for (i = 0; i < e->nslice; i++) {
        cram_index *tmp = cram_container_num2offset_(&e->e[i], num,
                                                     last_pos, nc);
        if (tmp)
            return tmp;
    }


    return NULL;
}

off_t cram_container_num2offset(cram_fd *fd, int64_t num) {
    int nc = 0, i;
    int64_t last_pos = -9;
    cram_index *e = NULL;

    for (i = 0; i < fd->index_sz; i++) {
        int j = i+1 == fd->index_sz ? 0 : i+1; // maps "*" to end
        if (!fd->index[j].nslice)
            continue;
        if ((e = cram_container_num2offset_(&fd->index[j], num,
                                            &last_pos, &nc)))
            break;
    }

    return e ? e->offset : -1;
}


/*! Returns the container number for the first container at offset >= pos.
 *
 * The index must have previously been loaded, otherwise <0 is returned.
 */
static cram_index *cram_container_offset2num_(cram_index *e, off_t pos,
                                              int64_t *last_pos, int *nc) {
    if (e->offset) {
        if (e->offset != *last_pos) {
            if (e->offset >= pos)
                return e;
            (*nc)++;
        }
        // else a new multi-ref in same container
        *last_pos = e->offset;
    }

    int i;
    for (i = 0; i < e->nslice; i++) {
        cram_index *tmp = cram_container_offset2num_(&e->e[i], pos,
                                                     last_pos, nc);
        if (tmp)
            return tmp;
    }


    return NULL;
}

int64_t cram_container_offset2num(cram_fd *fd, off_t pos) {
    int nc = 0, i;
    int64_t last_pos = -9;
    cram_index *e = NULL;

    for (i = 0; i < fd->index_sz; i++) {
        int j = i+1 == fd->index_sz ? 0 : i+1; // maps "*" to end
        if (!fd->index[j].nslice)
            continue;
        if ((e = cram_container_offset2num_(&fd->index[j], pos,
                                            &last_pos, &nc)))
            break;
    }

    return e ? nc : -1;
}

/*!
 * Returns the file offsets of CRAM containers covering a specific region
 * query.  Note both offsets are the START of the container.
 *
 * first will point to the start of the first overlapping container
 * last will point to the start of the last overlapping container
 *
 * Returns 0 on success
 *        <0 on failure
 */
int cram_index_extents(cram_fd *fd, int refid, hts_pos_t start, hts_pos_t end,
                       off_t *first, off_t *last) {
    cram_index *ci;

    if (first) {
        if (!(ci = cram_index_query(fd, refid, start, NULL)))
            return -1;
        *first = ci->offset;
    }

    if (last) {
        if (!(ci = cram_index_query_last(fd, refid, end)))
            return -1;
        *last = ci->offset;
    }

    return 0;
}
