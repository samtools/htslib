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

#include <config.h>

#include "htslib/hts.h"
#include "htslib/khash.h"

typedef struct reglist
{
    uint32_t n, m;
    uint64_t *a;
} reglist_t;

KHASH_MAP_INIT_STR(reg, reglist_t)
typedef kh_reg_t reghash_t;

static int compare_uint64 (const void * a, const void * b)
{
    if (*(uint64_t *)a < *(uint64_t *)b) return -1;
    if (*(uint64_t *)a > *(uint64_t *)b) return 1;

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
    const char *reg;
    uint32_t beg, end;

    if (!h) {
        fprintf(stderr, "Hash table is empty!\n");
        return;
    }
    for (k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            reg = kh_key(h,k);
            fprintf(stderr, "Region: '%s'\n", reg);
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                for (i=0; i<p->n; i++) {
                    beg = (uint32_t)(p->a[i]>>32);
                    end = (uint32_t)(p->a[i]);
                    fprintf(stderr, "\tinterval[%d]: %d-%d\n", i, beg, end);
                }
            } else {
                fprintf(stderr, "Region '%s' has no intervals!\n", reg);
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

        qsort(p->a, p->n, sizeof(uint64_t), compare_uint64);
        for (new_n = 0, j = 1; j < p->n; j++) {
            if ((uint32_t)p->a[new_n] < (uint32_t)(p->a[j]>>32)) {
                p->a[++new_n] = p->a[j];
            } else {
                if ((uint32_t)p->a[new_n] < (uint32_t)p->a[j])
                    p->a[new_n] = (p->a[new_n] & 0xFFFFFFFF00000000) | (uint32_t)(p->a[j]);
            }
        }
        p->n = ++new_n;
        count++;
    }

    return count;
}

static int reg_insert(reghash_t *h, char *reg, unsigned int beg, unsigned int end) {

    khint_t k;
    reglist_t *p;

    if (!h)
        return -1;

    // Put reg in the hash table if not already there
    k = kh_get(reg, h, reg); //looks strange, but only the second reg is the actual region name.
    if (k == kh_end(h)) { // absent from the hash table
        int ret;
        char *s = strdup(reg);
        if (NULL == s) return -1;
        k = kh_put(reg, h, s, &ret);
        if (-1 == ret) {
            free(s);
            return -1;
        }
        memset(&kh_val(h, k), 0, sizeof(reglist_t));
    }
    p = &kh_val(h, k);

    // Add beg and end to the list
    if (p->n == p->m) {
        uint32_t new_m = p->m ? p->m<<1 : 4;
        if (new_m == 0) return -1;
        uint64_t *new_a = realloc(p->a, new_m * sizeof(uint64_t));
        if (new_a == NULL) return -1;
        p->m = new_m;
        p->a = new_a;
    }
    p->a[p->n++] = (uint64_t)beg<<32 | end;

    return 0;
}

static void reg_destroy(reghash_t *h) {

    khint_t k;

    if (!h)
        return;

    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k).a);
            free((char*)kh_key(h, k));
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
    int i, l_count = 0;
    uint32_t j;
    char reg[1024];
    const char *q;
    int beg, end;

    /* First, transform the char array into a hash table */
    h = kh_init(reg);
    if (!h) {
        hts_log_error("Error when creating the region hash table");
        return NULL;
    }

    for (i=0; i<argc; i++) {
        q = hts_parse_reg(argv[i], &beg, &end);
        if (q) {
            if (q - argv[i] > sizeof(reg) - 1) {
                hts_log_error("Region name '%s' is too long (bigger than %d)", argv[i], (int) sizeof(reg) - 1);
                continue;
            }
            memcpy(reg, argv[i], q - argv[i]);
            reg[q - argv[i]] = 0;
        } else {
            // not parsable as a region, but possibly a sequence named "foo:a"
            if (strlen(argv[i]) > sizeof(reg) - 1) {
                hts_log_error("Region name '%s' is too long (bigger than %d)", argv[i], (int) sizeof(reg) - 1);
                continue;
            }
            strcpy(reg, argv[i]);
            beg = 0; end = INT_MAX;
        }

        if (reg_insert(h, reg, beg, end) != 0) {
            hts_log_error("Error when inserting region='%s' in the bed hash table at address=%p", argv[i], (void *) h);
            goto fail;
        }
    }

    *r_count = reg_compact(h);
    if (!*r_count)
        return NULL;

    /* Transform the hash table into a list */
    h_reglist = (hts_reglist_t *)calloc(*r_count, sizeof(hts_reglist_t));
    if (!h_reglist)
        return NULL;

    for (k = kh_begin(h); k < kh_end(h) && l_count < *r_count; k++) {
        if (!kh_exist(h,k) || !(p = &kh_val(h,k)))
            continue;

        char *reg_name = (char *)kh_key(h,k);
        if (!strcmp(reg_name, ".")) {
            h_reglist[l_count].tid = HTS_IDX_START;
        } else if (!strcmp(reg_name, "*")) {
            h_reglist[l_count].tid = HTS_IDX_NOCOOR;
        } else {
            h_reglist[l_count].tid = getid(hdr, reg_name);
            if (h_reglist[l_count].tid < 0)
                hts_log_warning("Region '%s' specifies an unknown reference name. Continue anyway", reg_name);
        }

        h_reglist[l_count].intervals = (hts_pair32_t *)calloc(p->n, sizeof(hts_pair32_t));
        if(!(h_reglist[l_count].intervals)) {
            hts_log_error("Could not allocate memory for intervals for region='%s'", kh_key(h,k));
            goto fail;
        }
        h_reglist[l_count].count = p->n;
        h_reglist[l_count].max_end = 0;

        for (j = 0; j < p->n; j++) {
            h_reglist[l_count].intervals[j].beg = (uint32_t)(p->a[j]>>32);
            h_reglist[l_count].intervals[j].end = (uint32_t)(p->a[j]);

            if (h_reglist[l_count].intervals[j].end > h_reglist[l_count].max_end)
                h_reglist[l_count].max_end = h_reglist[l_count].intervals[j].end;
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
