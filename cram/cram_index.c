/*
Copyright (c) 2013-2020 Genome Research Ltd.
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
    printf("%*soffset %"PRId64"\n", MAX(0,50-n), "", e->offset);
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

        while (!(e.start >= idx->start && e.end <= idx->end) || idx->end == 0) {
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

    // dump_index(fd);

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
 * return the first one only. Subsequent calls should specifying
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

    switch(refid) {
    case HTS_IDX_NONE:
    case HTS_IDX_REST:
        // fail, or already there, dealt with elsewhere.
        return NULL;

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

    if (!from)
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

    return &from->e[slice];
}

cram_index *cram_index_query_last(cram_fd *fd, int refid, hts_pos_t end) {
    cram_index *first = cram_index_query(fd, refid, end, NULL);
    cram_index *last =  cram_index_last(fd, refid, NULL);
    if (!first || !last)
        return NULL;

    while (first < last && (first+1)->start <= end)
        first++;

    while (first->e) {
        int count = 0;
        int nslices = first->nslice;
        first = first->e;
        while (++count < nslices && (first+1)->start <= end)
            first++;
    }

    // Compute the start location of next container.
    //
    // This is useful for stitching containers together in the multi-region
    // iterator.  Sadly we can't compute this from the single index line.
    //
    // Note we can have neighbouring index entries at the same location
    // for when we have multi-reference mode and/or multiple slices per
    // container.
    cram_index *next = first;
    do {
        if (next >= last) {
            // Next non-empty reference
            while (++refid+1 < fd->index_sz)
                if (fd->index[refid+1].nslice)
                    break;
            if (refid+1 >= fd->index_sz) {
                next = NULL;
            } else {
                next = fd->index[refid+1].e;
                last = fd->index[refid+1].e + fd->index[refid+1].nslice;
            }
        } else {
            next++;
        }
    } while (next && next->offset == first->offset);

    first->next = next ? next->offset : 0;

    return first;
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
            sprintf(buf, "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
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
        sprintf(buf, "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
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
        sprintf(buf, "%d\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%d\t%d\n",
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
        if (spos - cpos - c->offset != c->landmark[j]) {
            hts_log_error("CRAM slice offset %"PRId64" does not match"
                          " landmark %d in container header (%d)",
                          spos - cpos - c->offset, j, c->landmark[j]);
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

        cpos = htell(fd->fp);
        assert(cpos == hpos + c->length);

        cram_free_container(c);
    }
    if (fd->err) {
        bgzf_close(fp);
        return -1;
    }

    return (bgzf_close(fp) >= 0)? 0 : -4;
}
