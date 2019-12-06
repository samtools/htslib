/*  region.c -- Functions to create and free region lists

    Copyright (C) 2019 Genome Research Ltd.

    Author: Valeriu Ohan <vo2@sanger.ac.uk>

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include "htslib/hts.h"
#include "htslib/khash.h"

typedef struct reglist
{
    uint32_t n, m;
    hts_pair_pos_t *a;
    int tid;
} reglist_t;

KHASH_MAP_INIT_INT(reg, reglist_t)
typedef kh_reg_t reghash_t;

static int compare_hts_pair_pos_t (const void *av, const void *bv)
{
    hts_pair_pos_t *a = (hts_pair_pos_t *) av;
    hts_pair_pos_t *b = (hts_pair_pos_t *) bv;
    if (a->beg < b->beg) return -1;
    if (a->beg > b->beg) return  1;
    if (a->end < b->end) return -1;
    if (a->end > b->end) return  1;

    return 0;
}

#if 0
/**
 * Good to have around for debugging
 */
static void reg_print(reghash_t *h) {
    reglist_t *p;
    khint_t k;
    uint32_t i;
    khint32_t key;

    if (!h) {
        fprintf(stderr, "Hash table is empty!\n");
        return;
    }
    for (k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            key = kh_key(h,k);
            fprintf(stderr, "Region: key %u tid %d\n", key, p->tid);
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                for (i=0; i<p->n; i++) {
                    fprintf(stderr, "\tinterval[%d]: %"PRIhts_pos"-%"PRIhts_pos"\n", i,
                            p->a[i].beg, p->a[i].end);
                }
            } else {
                fprintf(stderr, "Region key %u has no intervals!\n", key);
            }
        }
    }
}
#endif

/**
 * Sort and merge overlapping or adjacent intervals.
 */
static int reg_compact(reghash_t *h) {
    khint_t i;
    uint32_t j, new_n;
    reglist_t *p;
    int count = 0;

    if (!h)
        return 0;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || !(p->n))
            continue;

        qsort(p->a, p->n, sizeof(p->a[0]), compare_hts_pair_pos_t);
        for (new_n = 0, j = 1; j < p->n; j++) {
            if (p->a[new_n].end < p->a[j].beg) {
                p->a[++new_n].beg = p->a[j].beg;
                p->a[new_n].end = p->a[j].end;
            } else {
                if (p->a[new_n].end < p->a[j].end)
                    p->a[new_n].end = p->a[j].end;
            }
        }
        ++new_n;
        if (p->n > new_n) {
            // Shrink array to required size.
            hts_pair_pos_t *new_a = realloc(p->a, new_n * sizeof(p->a[0]));
            if (new_a) p->a = new_a;
        }
        p->n = new_n;
        count++;
    }

    return count;
}

static int reg_insert(reghash_t *h, int tid, hts_pos_t beg, hts_pos_t end) {

    khint_t k;
    reglist_t *p;

    if (!h)
        return -1;

    // Put reg in the hash table if not already there
    k = kh_get(reg, h, tid);
    if (k == kh_end(h)) { // absent from the hash table
        int ret;
        k = kh_put(reg, h, tid, &ret);
        if (-1 == ret) {
            return -1;
        }
        memset(&kh_val(h, k), 0, sizeof(reglist_t));
        kh_val(h, k).tid = tid;
    }
    p = &kh_val(h, k);

    // Add beg and end to the list
    if (p->n == p->m) {
        uint32_t new_m = p->m ? p->m<<1 : 4;
        if (new_m == 0) return -1;
        hts_pair_pos_t *new_a = realloc(p->a, new_m * sizeof(p->a[0]));
        if (new_a == NULL) return -1;
        p->m = new_m;
        p->a = new_a;
    }
    p->a[p->n].beg = beg;
    p->a[p->n++].end = end;

    return 0;
}

static void reg_destroy(reghash_t *h) {

    khint_t k;

    if (!h)
        return;

    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k).a);
        }
    }
    kh_destroy(reg, h);
}

/**
 * Take a char array of reg:interval elements and produce a hts_reglis_t with r_count elements.
 */
hts_reglist_t *hts_reglist_create(char **argv, int argc, int *r_count, void *hdr,  hts_name2id_f getid) {

    if (!argv || argc < 1)
        return NULL;

    reghash_t *h = NULL;
    reglist_t *p;
    hts_reglist_t *h_reglist = NULL;

    khint_t k;
    int i, l_count = 0, tid;
    const char *q;
    hts_pos_t beg, end;

    /* First, transform the char array into a hash table */
    h = kh_init(reg);
    if (!h) {
        hts_log_error("Error when creating the region hash table");
        return NULL;
    }

    for (i=0; i<argc; i++) {
        if (!strcmp(argv[i], ".")) {
            q = argv[i] + 1;
            tid = HTS_IDX_START; beg = 0; end = INT64_MAX;
        } else if (!strcmp(argv[i], "*")) {
            q = argv[i] + 1;
            tid = HTS_IDX_NOCOOR; beg = 0; end = INT64_MAX;
        } else {
            q = hts_parse_region(argv[i], &tid, &beg, &end, getid, hdr,
                                 HTS_PARSE_THOUSANDS_SEP);
        }
        if (!q) {
            if (tid < -1) {
                hts_log_error("Failed to parse header");
                goto fail;
            } else {
                // not parsable as a region
                hts_log_warning("Region '%s' specifies an unknown reference name. Continue anyway", argv[i]);
                continue;
            }
        }

        if (reg_insert(h, tid, beg, end) != 0) {
            hts_log_error("Error when inserting region='%s' in the bed hash table at address=%p", argv[i], (void *) h);
            goto fail;
        }
    }

    *r_count = reg_compact(h);
    if (!*r_count)
        goto fail;

    /* Transform the hash table into a list */
    h_reglist = (hts_reglist_t *)calloc(*r_count, sizeof(hts_reglist_t));
    if (!h_reglist)
        goto fail;

    for (k = kh_begin(h); k < kh_end(h) && l_count < *r_count; k++) {
        if (!kh_exist(h,k) || !(p = &kh_val(h,k)))
            continue;

        h_reglist[l_count].tid = p->tid;
        h_reglist[l_count].intervals = p->a;
        h_reglist[l_count].count = p->n;
        p->a = NULL; // As we stole it.

        // After reg_compact(), list is ordered and non-overlapping, so...
        if (p->n > 0) {
            h_reglist[l_count].min_beg = h_reglist[l_count].intervals[0].beg;
            h_reglist[l_count].max_end = h_reglist[l_count].intervals[p->n - 1].end;
        } else {
            h_reglist[l_count].min_beg = 0;
            h_reglist[l_count].max_end = 0;
        }

        l_count++;
    }
    reg_destroy(h);

    return h_reglist;

fail:
    reg_destroy(h);
    if(h_reglist) hts_reglist_free(h_reglist, l_count);

    return NULL;
}

void hts_reglist_free(hts_reglist_t *reglist, int count) {

    int i;
    if(reglist) {
        for (i = 0; i < count; i++) {
            if (reglist[i].intervals)
                free(reglist[i].intervals);
        }
        free(reglist);
    }
}
