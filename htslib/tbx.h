/// @file htslib/tbx.h
/// Tabix API functions.
/*
    Copyright (C) 2009, 2012-2015, 2019 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#ifndef HTSLIB_TBX_H
#define HTSLIB_TBX_H

#include "hts.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TBX_MAX_SHIFT 31

#define TBX_GENERIC 0
#define TBX_SAM     1
#define TBX_VCF     2
#define TBX_GAF     3
#define TBX_UCSC    0x10000

typedef struct tbx_conf_t {
    int32_t preset;
    int32_t sc, bc, ec; // seq col., beg col. and end col.
    int32_t meta_char, line_skip;
} tbx_conf_t;

typedef struct tbx_t {
    tbx_conf_t conf;
    hts_idx_t *idx;
    void *dict;
} tbx_t;

HTSLIB_EXPORT
extern const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf, tbx_conf_gaf;

    #define tbx_itr_destroy(iter) hts_itr_destroy(iter)
    #define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
    #define tbx_itr_querys(tbx, s) tbx_itr_querys1((tbx), (s))
    #define tbx_itr_next(htsfp, tbx, itr, r) hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))
    #define tbx_bgzf_itr_next(bgzfp, tbx, itr, r) hts_itr_next((bgzfp), (itr), (r), (tbx))

    HTSLIB_EXPORT
    hts_itr_t *tbx_itr_querys1(tbx_t *tbx, const char *region);

    HTSLIB_EXPORT
    int tbx_name2id(tbx_t *tbx, const char *ss);

    /* Internal helper function used by tbx_itr_next() */
    HTSLIB_EXPORT
    BGZF *hts_get_bgzfp(htsFile *fp);

    HTSLIB_EXPORT
    int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end);

    /// Construct a multi-region iterator over a tabix indexed text file.
    /** Returns an hts_itr_t that, when driven by hts_itr_multi_next, yields
        records from all the regions in reglist in genome order. Adjacent
        and nearby regions are coalesced into combined tabix index lookups
        by hts_itr_multi_next, so this is meaningfully faster than calling
        tbx_itr_queryi once per region when the region list is large or
        dense. See test/bench_tbx_regions.c for a reproducible measurement
        of the speedup curve.

        Reglist entries may specify their target contig either way:
          - reglist[i].tid set to a tid (use tbx_name2id() to resolve)
            with reglist[i].reg == NULL, or
          - reglist[i].reg set to a reference name string such as "chr1"
            (or "." for all-with-coords records, "*" for unmapped),
            which this function resolves internally via tbx_name2id().
        Each reglist[i].intervals must hold reglist[i].count
        hts_pair_pos_t entries.

        Semantics differ from running tbx_itr_queryi() multiple times in
        one respect: each underlying file record is yielded at most once,
        even when multiple intervals in reglist cover it. Duplicate or
        overlapping intervals produce a single emission per matching
        record, not one per matching interval. Callers that need
        per-interval multiplicity must call tbx_itr_queryi() per interval.

        Ownership: the caller's reglist is deep-copied internally and is
        never mutated or freed by this function regardless of outcome.
        The caller retains ownership of reglist (and is responsible for
        freeing it) as well as of tbx and fp. tbx must remain valid for
        the lifetime of the iterator. Destroy the iterator with
        hts_itr_destroy.

        Concurrency: only one multi-region tabix iterator may be active
        on a given htsFile at a time. The iterator uses the BGZF
        private_data slot on fp to thread tbx_t through hts_itr_multi_next;
        constructing a second multi-region iterator on the same fp before
        destroying the first would clobber that slot and yield wrong
        records from the first. Single-region iterators (tbx_itr_queryi)
        on the same fp are unaffected. This matches the existing
        constraint on binary BCF iteration, which uses the same slot.

        Lifetime: destroy the iterator (hts_itr_destroy) before destroying
        the tbx (tbx_destroy) or closing the fp (hts_close). The iterator
        retains a pointer to tbx via the BGZF private_data slot; destroying
        tbx first leaves a dangling pointer in fp's BGZF cache that would
        be returned by a later bgzf_get_private_data on the same fp.

        Returns NULL if fp, tbx, or reglist is NULL; if count <= 0; if
        fp is not BGZF; if memory cannot be allocated; or if the iterator
        cannot otherwise be constructed. */
    HTSLIB_EXPORT
    hts_itr_t *tbx_itr_regions(htsFile *fp, tbx_t *tbx,
                                hts_reglist_t *reglist, int count);

/// Build an index of the lines in a BGZF-compressed file
/** The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index(BGZF *fp, int min_shift, const tbx_conf_t *conf);
/*
 * All tbx_index_build* methods return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
    HTSLIB_EXPORT
    int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);

    HTSLIB_EXPORT
    int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf);

    HTSLIB_EXPORT
    int tbx_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads, const tbx_conf_t *conf);


/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index

    Equivalent to tbx_index_load3(fn, NULL, HTS_IDX_SAVE_REMOTE);
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load(const char *fn);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    Equivalent to tbx_index_load3(fn, fnidx, HTS_IDX_SAVE_REMOTE);
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load2(const char *fn, const char *fnidx);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @param flags  Flags to alter behaviour (see description)
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    The @p flags parameter can be set to a combination of the following
    values:

        HTS_IDX_SAVE_REMOTE   Save a local copy of any remote indexes
        HTS_IDX_SILENT_FAIL   Fail silently if the index is not present

    The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    HTSLIB_EXPORT
    tbx_t *tbx_index_load3(const char *fn, const char *fnidx, int flags);

    HTSLIB_EXPORT
    const char **tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values

    HTSLIB_EXPORT
    void tbx_destroy(tbx_t *tbx);

#ifdef __cplusplus
}
#endif

#endif
