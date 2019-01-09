/*
Copyright (c) 2018 Genome Research Ltd.
Authors: James Bonfield <jkb@sanger.ac.uk>, Valeriu Ohan <vo2@sanger.ac.uk>

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

#include <config.h>

#include <string.h>
#include <assert.h>
#include "textutils_internal.h"
#include "header.h"

// Hash table for removing multiple lines from the header
KHASH_SET_INIT_STR(rm)
typedef khash_t(rm) rmhash_t;

#define MAX_ERROR_QUOTE 320 // Prevent over-long error messages
static void bam_hrecs_error(const char *msg, const char *line, size_t len, size_t lno) {
    int j;

    if (len > MAX_ERROR_QUOTE)
        len = MAX_ERROR_QUOTE;
    for (j = 0; j < len && line[j] != '\n'; j++)
        ;
    hts_log_error("%s at line %zd: \"%.*s\"", msg, lno, j, line);
}

/* ==== Static methods ==== */

static int bam_hrecs_init_type_order(bam_hrecs_t *hrecs, char *type_list) {
    if (!hrecs)
        return -1;

    if (!type_list) {
        hrecs->type_count = 5;
        hrecs->type_order = calloc(hrecs->type_count, 3);
        if (!hrecs->type_order)
            return -1;
        memcpy(hrecs->type_order[0], "HD", 2);
        memcpy(hrecs->type_order[1], "SQ", 2);
        memcpy(hrecs->type_order[2], "RG", 2);
        memcpy(hrecs->type_order[3], "PG", 2);
        memcpy(hrecs->type_order[4], "CO", 2);
    }

    return 0;
}

/* Updates the hash tables in the bam_hrecs_t structure.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int bam_hrecs_update_hashes(bam_hrecs_t *hrecs,
                                   int type,
                                   bam_hrec_type_t *h_type) {
    /* Add to reference hash? */
    if ((type>>8) == 'S' && (type&0xff) == 'Q') {
        bam_hrec_tag_t *tag;
        bam_hrec_sq_t *new_ref;
        int nref = hrecs->nref;

        new_ref = realloc(hrecs->ref, (hrecs->nref+1)*sizeof(*hrecs->ref));
        if (!new_ref)
            return -1;
        hrecs->ref = new_ref;

        tag = h_type->tag;
        hrecs->ref[nref].name = NULL;
        hrecs->ref[nref].len  = 0;
        hrecs->ref[nref].ty = h_type;
        hrecs->ref[nref].tag  = tag;

        while (tag) {
            if (tag->str[0] == 'S' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                hrecs->ref[nref].name = tag->str+3;
            } else if (tag->str[0] == 'L' && tag->str[1] == 'N') {
                hrecs->ref[nref].len = atoi(tag->str+3);
            }
            tag = tag->next;
        }

        if (hrecs->ref[nref].name) {
            khint_t k;
            int r;
            k = kh_put(m_s2i, hrecs->ref_hash, hrecs->ref[nref].name, &r);
            if (-1 == r) return -1;
            kh_val(hrecs->ref_hash, k) = nref;
        } else {
            return -1; // SN should be present, according to spec.
        }

        if (hrecs->refs_changed < 0 || hrecs->refs_changed > hrecs->nref)
            hrecs->refs_changed = hrecs->nref;
        hrecs->nref++;
    }

    /* Add to read-group hash? */
    if ((type>>8) == 'R' && (type&0xff) == 'G') {
        bam_hrec_tag_t *tag;
        bam_hrec_rg_t *new_rg;
        int nrg = hrecs->nrg;

        new_rg = realloc(hrecs->rg, (hrecs->nrg+1)*sizeof(*hrecs->rg));
        if (!new_rg)
            return -1;
        hrecs->rg = new_rg;

        tag = h_type->tag;
        hrecs->rg[nrg].name = NULL;
        hrecs->rg[nrg].name_len = 0;
        hrecs->rg[nrg].ty   = h_type;
        hrecs->rg[nrg].tag  = tag;
        hrecs->rg[nrg].id   = nrg;

        while (tag) {
            if (tag->str[0] == 'I' && tag->str[1] == 'D') {
                assert(tag->len >= 3);
                hrecs->rg[nrg].name = tag->str + 3;
                hrecs->rg[nrg].name_len = tag->len - 3;
            }
            tag = tag->next;
        }

        if (hrecs->rg[nrg].name) {
            khint_t k;
            int r;
            k = kh_put(m_s2i, hrecs->rg_hash, hrecs->rg[nrg].name, &r);
            if (-1 == r) return -1;
            kh_val(hrecs->rg_hash, k) = nrg;
        } else {
            return -1; // ID should be present, according to spec.
        }

        hrecs->nrg++;
    }

    /* Add to program hash? */
    if ((type>>8) == 'P' && (type&0xff) == 'G') {
        bam_hrec_tag_t *tag;
        bam_hrec_pg_t *new_pg;
        int npg = hrecs->npg;

        new_pg = realloc(hrecs->pg, (hrecs->npg+1)*sizeof(*hrecs->pg));
        if (!new_pg)
            return -1;
        hrecs->pg = new_pg;

        tag = h_type->tag;
        hrecs->pg[npg].name = NULL;
        hrecs->pg[npg].name_len = 0;
        hrecs->pg[npg].ty  = h_type;
        hrecs->pg[npg].tag  = tag;
        hrecs->pg[npg].id   = npg;
        hrecs->pg[npg].prev_id = -1;

        while (tag) {
            if (tag->str[0] == 'I' && tag->str[1] == 'D') {
                assert(tag->len >= 3);
                hrecs->pg[npg].name = tag->str + 3;
                hrecs->pg[npg].name_len = tag->len - 3;
            } else if (tag->str[0] == 'P' && tag->str[1] == 'P') {
                // Resolve later if needed
                khint_t k;
                char tmp = tag->str[tag->len]; tag->str[tag->len] = 0;
                k = kh_get(m_s2i, hrecs->pg_hash, tag->str+3);
                tag->str[tag->len] = tmp;

                if (k != kh_end(hrecs->pg_hash)) {
                    int p_id = kh_val(hrecs->pg_hash, k);
                    hrecs->pg[npg].prev_id = hrecs->pg[p_id].id;

                    /* Unmark previous entry as a PG termination */
                    if (hrecs->npg_end > 0 &&
                        hrecs->pg_end[hrecs->npg_end-1] == p_id) {
                        hrecs->npg_end--;
                    } else {
                        int i;
                        for (i = 0; i < hrecs->npg_end; i++) {
                            if (hrecs->pg_end[i] == p_id) {
                                memmove(&hrecs->pg_end[i], &hrecs->pg_end[i+1],
                                        (hrecs->npg_end-i-1)*sizeof(*hrecs->pg_end));
                                hrecs->npg_end--;
                            }
                        }
                    }
                } else {
                    hrecs->pg[npg].prev_id = -1;
                }
            }
            tag = tag->next;
        }

        if (hrecs->pg[npg].name) {
            khint_t k;
            int r;
            k = kh_put(m_s2i, hrecs->pg_hash, hrecs->pg[npg].name, &r);
            if (-1 == r) return -1;
            kh_val(hrecs->pg_hash, k) = npg;
        } else {
            return -1; // ID should be present, according to spec.
        }

        /* Add to npg_end[] array. Remove later if we find a PP line */
        if (hrecs->npg_end >= hrecs->npg_end_alloc) {
            int *new_pg_end;
            int  new_alloc = hrecs->npg_end_alloc ? hrecs->npg_end_alloc*2 : 4;

            new_pg_end = realloc(hrecs->pg_end, new_alloc * sizeof(int));
            if (!new_pg_end)
                return -1;
            hrecs->npg_end_alloc = new_alloc;
            hrecs->pg_end = new_pg_end;
        }
        hrecs->pg_end[hrecs->npg_end++] = npg;

        hrecs->npg++;
    }

    return 0;
}

static int bam_hrecs_remove_hash_entry(bam_hrecs_t *hrecs, int type, bam_hrec_type_t *h_type) {
    if (!hrecs || !h_type)
        return -1;

    bam_hrec_tag_t *tag;
    char *key = NULL;
    khint_t k;

    /* Remove from reference hash */
    if ((type>>8) == 'S' && (type&0xff) == 'Q') {
        tag = h_type->tag;

        while (tag) {
            if (tag->str[0] == 'S' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                key = tag->str + 3;
                k = kh_get(m_s2i, hrecs->ref_hash, key);
                if (k != kh_end(hrecs->ref_hash)) {
                    int idx = kh_val(hrecs->ref_hash, k);
                    if (idx < hrecs->nref-1)
                        memmove(&hrecs->ref[idx], &hrecs->ref[idx+1],
                                sizeof(bam_hrec_sq_t)*(hrecs->nref - idx - 1));
                    kh_del(m_s2i, hrecs->ref_hash, k);
                    hrecs->nref--;
                    if (hrecs->refs_changed < 0 || hrecs->refs_changed > idx)
                        hrecs->refs_changed = idx;
                }
                break;
            }
            tag = tag->next;
        }
    }

    /* Remove from read-group hash */
    if ((type>>8) == 'R' && (type&0xff) == 'G') {
        tag = h_type->tag;

        while (tag) {
            if (tag->str[0] == 'I' && tag->str[1] == 'D') {
                assert(tag->len >= 3);
                key = tag->str + 3;
                k = kh_get(m_s2i, hrecs->rg_hash, key);
                if (k != kh_end(hrecs->rg_hash)) {
                    if (kh_val(hrecs->rg_hash, k) < hrecs->nrg-1)
                        memmove(&hrecs->rg[kh_val(hrecs->rg_hash, k)], &hrecs->rg[kh_val(hrecs->rg_hash, k)+1], sizeof(bam_hrec_rg_t)*(hrecs->nrg - kh_val(hrecs->rg_hash, k) - 1));
                    kh_del(m_s2i, hrecs->rg_hash, k);
                    hrecs->nrg--;
                }
                break;
            }
            tag = tag->next;
        }
    }

    return 0;
}

/*! Add header record with a va_list interface.
 *
 * Adds a single record to a SAM header.
 *
 * This takes a header record type, a va_list argument and one or more
 * key,value pairs, ending with the NULL key.
 *
 * Eg. sam_hdr_vadd(h, "SQ", args, "ID", "foo", "LN", "100", NULL).
 *
 * The purpose of the additional va_list parameter is to permit other
 * varargs functions to call this while including their own additional
 * parameters; an example is in sam_hdr_add_PG().
 *
 * Note: this function invokes va_arg at least once, making the value
 * of ap indeterminate after the return. The caller should call
 * va_start/va_end before/after calling this function or use va_copy.
 *
 * @return
 * Returns >= 0 on success;
 *        -1 on failure
 */
static int bam_hrecs_vadd(bam_hrecs_t *hrecs, const char *type, va_list ap, ...) {
    va_list args;
    bam_hrec_type_t *h_type;
    bam_hrec_tag_t *h_tag, *last=NULL;
    int new;
    khint32_t type_i = (type[0]<<8) | type[1], k;

    if (!(h_type = pool_alloc(hrecs->type_pool)))
        return -1;
    if (-1 == (k = kh_put(bam_hrecs_t, hrecs->h, type_i, &new)))
        return -1;

    // Form the ring, either with self or other lines of this type
    if (!new) {
        bam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
        p = t->prev;

        assert(p->next == t);
        p->next = h_type;
        h_type->prev = p;

        t->prev = h_type;
        h_type->next = t;
        h_type->order = p->order + 1;
    } else {
        kh_val(hrecs->h, k) = h_type;
        h_type->prev = h_type->next = h_type;
        h_type->order = 0;
    }
    h_type->skip = 0;
    h_type->tag = NULL;
    h_type->comm = NULL;

    // Any ... varargs
    va_start(args, ap);
    for (;;) {
        char *k, *v = NULL;

        if (!(k = (char *)va_arg(args, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(v = (char *)va_arg(args, char *)))
            break;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(v);
            h_tag->str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!h_tag->str || snprintf(h_tag->str, h_tag->len+1, "%2.2s:%s", k, v) < 0)
                return -1;
        } else {
            h_tag->len = strlen(k);
            h_tag->str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!h_tag->str || snprintf(h_tag->str, h_tag->len+1, "%s", k) < 0)
                return -1;
        }

        h_tag->next = NULL;
        if (last)
            last->next = h_tag;
        else
            h_type->tag = h_tag;

        last = h_tag;
    }
    va_end(args);

    // Plus the specified va_list params
    for (;;) {
        char *k, *v = NULL;

        if (!(k = (char *)va_arg(ap, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(v = (char *)va_arg(ap, char *)))
            break;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(v);
            h_tag->str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!h_tag->str || snprintf(h_tag->str, h_tag->len+1, "%2.2s:%s", k, v) < 0)
                return -1;
        } else {
            h_tag->len = strlen(k);
            h_tag->str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!h_tag->str || snprintf(h_tag->str, h_tag->len+1, "%s", k) < 0)
                return -1;
        }

        h_tag->next = NULL;
        if (last)
            last->next = h_tag;
        else
            h_type->tag = h_tag;

        last = h_tag;
    }

    int itype = (type[0]<<8) | type[1];
    if (-1 == bam_hrecs_update_hashes(hrecs, itype, h_type))
        return -1;

    hrecs->dirty = 1;

    return h_type->order;
}

/*
 * Function for deallocating a list of tags
 */

static void bam_hrecs_free_tags(bam_hrecs_t *hrecs, bam_hrec_tag_t *tag) {
    if (!hrecs || !tag)
        return;
    if (tag->next)
        bam_hrecs_free_tags(hrecs, tag->next);

    pool_free(hrecs->tag_pool, tag);
}

static int bam_hrecs_remove_line(bam_hrecs_t *hrecs, const char *type_name, bam_hrec_type_t *type_found) {
    if (!hrecs || !type_name || !type_found)
        return -1;

    int itype = (type_name[0]<<8) | type_name[1];
    khint_t k = kh_get(bam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h))
        return -1;

    bam_hrecs_free_tags(hrecs, type_found->tag);
    /* single element in the list */
    if (type_found->prev == type_found || type_found->next == type_found) {
        kh_del(bam_hrecs_t, hrecs->h, k);
    } else {
        type_found->prev->next = type_found->next;
        type_found->next->prev = type_found->prev;
        if (!type_found->order) { //first element
            kh_val(hrecs->h, k) = type_found->next;
            type_found->next->order = 0; //mark the new element as the first
        }
    }

    if (!strncmp(type_name, "SQ", 2) || !strncmp(type_name, "RG", 2))
        bam_hrecs_remove_hash_entry(hrecs, itype, type_found);

    pool_free(hrecs->type_pool, type_found);

    hrecs->dirty = 1;

    return 0;
}

static int bam_hrec_rebuild_comments(const bam_hrec_type_t *t, kstring_t *ks) {
    if (!t->comm) {
        if (EOF == kputsn_("@CO\t", 4, ks))
            return -1;
        if (t->tag)
            if (EOF == kputsn_(t->tag->str, t->tag->len, ks))
                return -1;
        if (EOF == kputc('\n', ks))
            return -1;
    } else {
        if (bam_hrec_rebuild_comments(t->comm, ks))
            return -1;
        if (EOF == kputsn_("@CO\t", 4, ks))
            return -1;
        if (t->tag)
            if (EOF == kputsn_(t->tag->str, t->tag->len, ks))
                return -1;
        if (EOF == kputc('\n', ks))
            return -1;
    }

    return 0;
}

static int bam_hrecs_rebuild_lines(const bam_hrecs_t *hrecs, khint_t k, kstring_t *ks) {
    const bam_hrec_type_t *t1, *t2;

    t1 = t2 = kh_val(hrecs->h, k);
    do {
        bam_hrec_tag_t *tag;
        char c[2];

        if (t1->skip) {
            t1 = t1->next;
            continue;
        }

        /* if there is a comment line attached, print it first */
        if (t1->comm) {
            if (bam_hrec_rebuild_comments(t1->comm, ks))
                return -1;
        }

        if (EOF == kputc_('@', ks))
            return -1;
        c[0] = kh_key(hrecs->h, k)>>8;
        c[1] = kh_key(hrecs->h, k)&0xff;
        if (EOF == kputsn_(c, 2, ks))
            return -1;
        for (tag = t1->tag; tag; tag=tag->next) {
            if (EOF == kputc_('\t', ks))
                return -1;
            if (EOF == kputsn_(tag->str, tag->len, ks))
                return -1;
        }
        if (EOF == kputc('\n', ks))
            return -1;
        t1 = t1->next;
    } while (t1 != t2);

    return 0;
}

static int bam_hrecs_parse_lines(bam_hrecs_t *hrecs, const char *hdr, size_t len) {
    size_t i, lno;
    bam_hrec_type_t *last_comm = NULL;

    if (!hrecs || len > SSIZE_MAX)
        return -1;

    if (!len)
        len = strlen(hdr);

    if (len < 3) {
        if (len == 0 || *hdr == '\0') return 0;
        bam_hrecs_error("Header line too short", hdr, len, 1);
        return -1;
    }

    for (i = 0, lno = 1; i < len - 3 && hdr[i] != '\0'; i++, lno++) {
        khint32_t type;
        khint_t k;

        int l_start = i, new;
        bam_hrec_type_t *h_type;
        bam_hrec_tag_t *h_tag, *last;

        if (hdr[i] != '@') {
            bam_hrecs_error("Header line does not start with '@'",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }

        type = (hdr[i+1]<<8) | hdr[i+2];
        if (!isalpha_c(hdr[i+1]) || !isalpha_c(hdr[i+2])) {
            bam_hrecs_error("Header line does not have a two character key",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }

        i += 3;
        if (i == len || hdr[i] == '\n')
            continue;

        // Add the header line type
        if (!(h_type = pool_alloc(hrecs->type_pool)))
            return -1;
        if (-1 == (k = kh_put(bam_hrecs_t, hrecs->h, type, &new)))
            return -1;

        // Form the ring, either with self or other lines of this type
        if (!new) {
            bam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
            p = t->prev;

            assert(p->next == t);
            p->next = h_type;
            h_type->prev = p;

            t->prev = h_type;
            h_type->next = t;
            h_type->order = p->order+1;
        } else {
            kh_val(hrecs->h, k) = h_type;
            h_type->prev = h_type->next = h_type;
            h_type->order = 0;
        }
        h_type->skip = 0;

        // Parse the tags on this line
        last = NULL;
        if ((type>>8) == 'C' && (type&0xff) == 'O') {
            size_t j;

            h_type->comm = last_comm;
            if (last_comm)
                last_comm->skip = 1;
            last_comm = h_type;

            if (i == len || hdr[i] != '\t') {
                bam_hrecs_error("Missing tab",
                              &hdr[l_start], len - l_start, lno);
                return -1;
            }

            for (j = ++i; j < len && hdr[j] != '\0' && hdr[j] != '\n'; j++)
                ;

            if (!(h_type->tag = h_tag = pool_alloc(hrecs->tag_pool)))
                return -1;
            h_tag->str = string_ndup(hrecs->str_pool, &hdr[i], j-i);
            h_tag->len = j-i;
            h_tag->next = NULL;
            if (!h_tag->str)
                return -1;

            i = j;

        } else {
            h_type->comm = last_comm;
            if (last_comm) {
                last_comm->skip = 1;
                last_comm = NULL;
            }

            do {
                size_t j;

                if (i == len || hdr[i] != '\t') {
                    bam_hrecs_error("Missing tab",
                                  &hdr[l_start], len - l_start, lno);
                    return -1;
                }

                for (j = ++i; j < len && hdr[j] != '\0' && hdr[j] != '\n' && hdr[j] != '\t'; j++)
                    ;

                if (j - i < 3 || hdr[i + 2] != ':') {
                    bam_hrecs_error("Malformed key:value pair",
                                   &hdr[l_start], len - l_start, lno);
                    return -1;
                }

                if (!(h_tag = pool_alloc(hrecs->tag_pool)))
                    return -1;
                h_tag->str = string_ndup(hrecs->str_pool, &hdr[i], j-i);
                h_tag->len = j-i;
                h_tag->next = NULL;
                if (!h_tag->str)
                    return -1;

                if (last)
                    last->next = h_tag;
                else
                    h_type->tag = h_tag;

                last = h_tag;
                i = j;
            } while (i < len && hdr[i] != '\0' && hdr[i] != '\n');
        }

        /* Update RG/SQ hashes */
        if (-1 == bam_hrecs_update_hashes(hrecs, type, h_type))
            return -1;
    }

    return 0;
}

/*! Update bam_hdr_t target_name and target_len arrays
 *
 *  @return 0 on success; -1 on failure
 */
int update_target_arrays(bam_hdr_t *bh, const bam_hrecs_t *hrecs,
                         int refs_changed) {
    if (!bh || !hrecs)
        return -1;

    if (refs_changed < 0)
        return 0;

    // Grow arrays if necessary
    if (bh->n_targets < hrecs->nref) {
        char **new_names = realloc(bh->target_name,
                                   hrecs->nref * sizeof(*new_names));
        if (!new_names)
            return -1;
        bh->target_name = new_names;
        uint32_t *new_lens = realloc(bh->target_len,
                                     hrecs->nref * sizeof(*new_lens));
        if (!new_lens)
            return -1;
        bh->target_len = new_lens;
    }

    // Update names and lengths where changed
    // hrecs->refs_changed is the first ref that has been updated, so ones
    // before that can be skipped.
    int i;
    for (i = refs_changed; i < hrecs->nref; i++) {
        if (i >= bh->n_targets
            || strcmp(bh->target_name[i], hrecs->ref[i].name) != 0) {
            if (i < bh->n_targets)
                free(bh->target_name[i]);
            bh->target_name[i] = strdup(hrecs->ref[i].name);
            if (!bh->target_name[i])
                return -1;
        }
        bh->target_len[i] = hrecs->ref[i].len;
    }

    // Free up any names that have been removed
    for (; i < bh->n_targets; i++) {
        free(bh->target_name[i]);
    }

    bh->n_targets = hrecs->nref;
    return 0;
}

static int rebuild_target_arrays(bam_hdr_t *bh) {
    if (!bh || !bh->hrecs)
        return -1;

    bam_hrecs_t *hrecs = bh->hrecs;
    if (hrecs->refs_changed < 0)
        return 0;

    if (update_target_arrays(bh, hrecs, hrecs->refs_changed) != 0)
        return -1;

    hrecs->refs_changed = -1;
    return 0;
}

int bam_hdr_parse(bam_hdr_t *bh) {
    bam_hrecs_t *hrecs = bam_hrecs_new();

    if (!hrecs)
        return -1;

    // Parse existing header text
    if (bh->text && bh->l_text > 0) {
        if (bam_hrecs_parse_lines(hrecs, bh->text, bh->l_text) != 0) {
            bam_hrecs_free(hrecs);
            return -1;
        }
    }

    bh->hrecs = hrecs;

    if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
        return -1;

    bam_hdr_link_pg(bh);

    return 0;
}

/** Remove outdated header text

    @param bh     BAM header

    This is called when API functions have changed the header so that the
    text version is no longer valid.
 */
static void redact_header_text(bam_hdr_t *bh) {
    assert(bh->hrecs && bh->hrecs->dirty);
    bh->l_text = 0;
    free(bh->text);
    bh->text = NULL;
}

/* ==== Public methods ==== */

int bam_hdr_length(bam_hdr_t *bh) {
    if (-1 == bam_hdr_rebuild(bh))
        return -1;

    return bh->l_text;
}

const char *bam_hdr_str(bam_hdr_t *bh) {
    if (-1 == bam_hdr_rebuild(bh))
        return NULL;

    return bh->text;
}

/*
 * Reconstructs the kstring from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hdr_rebuild(bam_hdr_t *bh) {
    bam_hrecs_t *hrecs;
    if (!bh || !(hrecs = bh->hrecs))
        return -1;

    if (hrecs->refs_changed >= 0) {
        if (rebuild_target_arrays(bh) < 0)
            return -1;
    }

    /* If header text wasn't changed or header is empty, don't rebuild it. */
    if (!hrecs->dirty)
        return 0;

    kstring_t ks = KS_INITIALIZER;
    if (bam_hrecs_rebuild_text(hrecs, &ks) != 0) {
        KS_FREE(&ks);
        return -1;
    }

    hrecs->dirty = 0;

    /* Sync */
    free(bh->text);
    bh->l_text = ks_len(&ks);
    bh->text = ks_release(&ks);

    return 0;
}

/*
 * Appends a formatted line to an existing SAM header.
 * Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 * optional new-line. If it contains more than 1 line then multiple lines
 * will be added in order.
 *
 * Input text is of maximum length len or as terminated earlier by a NUL.
 * len may be 0 if unknown, in which case lines must be NUL-terminated.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hdr_add_lines(bam_hdr_t *bh, const char *lines, int len) {
    bam_hrecs_t *hrecs;

    if (!bh || !lines)
        return -1;

    if (len == 0 && *lines == '\0')
        return 0;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (bam_hrecs_parse_lines(hrecs, lines, len) != 0)
        return -1;

    if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
        return -1;

    hrecs->dirty = 1;
    redact_header_text(bh);

    return 0;
}

/*
 * Adds a single line to a SAM header.
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. bam_hdr_add(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * Returns index for specific entry on success (eg 2nd SQ, 4th RG)
 *        -1 on failure
 */
int bam_hdr_add_line(bam_hdr_t *bh, const char *type, ...) {
    va_list args;
    bam_hrecs_t *hrecs;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    va_start(args, type);
    int ret = bam_hrecs_vadd(hrecs, type, args, NULL);
    va_end(args);

    if (ret >= 0) {
        if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
            return -1;

        if (hrecs->dirty)
            redact_header_text(bh);
    }

    return ret;
}


/*
 * Returns a complete line of formatted text for a specific head type/ID
 * combination. If ID is NULL then it returns the first line of the specified
 * type.
 *
 * The returned string is malloced and should be freed by the calling
 * function with free().
 *
 * Returns NULL if no type/ID is found.
 */
char *bam_hdr_find_line(bam_hdr_t *bh, const char *type, ...) {
    //                    ) {
    va_list args;
    bam_hrecs_t *hrecs;
    if (!bh)
        return NULL;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return NULL;
        hrecs = bh->hrecs;
    }

    const char *ID_key, *ID_value;
    va_start(args, type);
    ID_key = (char *)va_arg(args, char *);
    if (ID_key)
        ID_value = (char *)va_arg(args, char *);
    else
        ID_value = NULL;
    va_end(args);

    bam_hrec_type_t *ty = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!ty) {
        hts_log_warning("Could not find type '%s'", type);
        return NULL;
    }

    kstring_t ks = KS_INITIALIZER;
    bam_hrec_tag_t *tag;
    int r = 0;
    // Paste together the line from the hashed copy
    r |= (kputc_('@', &ks) == EOF);
    r |= (kputs(type, &ks) == EOF);
    for (tag = ty->tag; tag; tag = tag->next) {
        r |= (kputc_('\t', &ks) == EOF);
        r |= (kputsn(tag->str, tag->len, &ks) == EOF);
    }

    if (r) {
        KS_FREE(&ks);
        return NULL;
    }

    return ks_str(&ks);
}

/*
 * Remove a line from the header by specifying a tag:value that uniquely
 * identifies a line, i.e. the @SQ line containing "SN:ref1".
 * @SQ line is uniquely identified by SN tag.
 * @RG line is uniquely identified by ID tag.
 * @PG line is uniquely identified by ID tag.
 *
 * Returns 0 on success and -1 on error
 */

int bam_hdr_remove_line_key(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value) {
    bam_hrecs_t *hrecs;
    if (!bh || !type)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!strncmp(type, "PG", 2)) {
        hts_log_warning("Removing PG lines is not supported!");
        return -1;
    }

    bam_hrec_type_t *type_found = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!type_found)
        return 0;

    int ret = bam_hrecs_remove_line(hrecs, type, type_found);
    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

/*
 * Remove a line from the header by specifying the position in the type
 * group, i.e. 3rd @SQ line.
 *
 * Returns 0 on success and -1 on error
 */

int bam_hdr_remove_line_pos(bam_hdr_t *bh, const char *type, int position) {
    bam_hrecs_t *hrecs;
    if (!bh || !type || position < 0)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!strncmp(type, "PG", 2)) {
        hts_log_warning("Removing PG lines is not supported!");
        return -1;
    }

    bam_hrec_type_t *type_beg, *type_end, *type_found = NULL;
    int itype = (type[0]<<8) | type[1];

    khint_t k = kh_get(bam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h))
        return -1;

    type_beg = type_end = kh_val(hrecs->h, k);
    do {
        if (!type_found) {
            if(type_beg->order == position-1) {
                type_found = type_beg;
                break;
            }
        }
        type_beg = type_beg->next;
    } while (type_beg != type_end);

    if (!type_found) // nothing to remove
        return 0;
    int ret = bam_hrecs_remove_line(hrecs, type, type_found);
    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

int bam_hdr_update_line(bam_hdr_t *bh, const char *type,
        const char *ID_key, const char *ID_value, ...) {
    bam_hrecs_t *hrecs;
    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    int ret;
    bam_hrec_type_t *ty = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!ty)
        return -1;

    va_list args;
    va_start(args, ID_value);
    ret = bam_hrecs_update(hrecs, ty, args);
    va_end(args);

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

int bam_hdr_keep_line(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value) {
    bam_hrecs_t *hrecs;
    if (!bh || !type)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    bam_hrec_type_t *step;
    int ret = 1, remove_all = 0;

    if (!strncmp(type, "PG", 2) || !strncmp(type, "CO", 2)) {
        hts_log_warning("Removing PG or CO lines is not supported!");
        return -1;
    }

    bam_hrec_type_t *type_found = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!type_found) { // remove all line of this type
        int itype = (type[0]<<8)|(type[1]);
        khint_t k = kh_get(bam_hrecs_t, hrecs->h, itype);
        if (k == kh_end(hrecs->h))
            return 0;
        type_found =  kh_val(hrecs->h, k);
        if (!type_found)
            return 0;
        remove_all = 1;
    }

    step = type_found->next;
    while (step != type_found) {
        bam_hrec_type_t *to_remove = step;
        step = step->next;
        ret &= bam_hrecs_remove_line(hrecs, type, to_remove);
    }

    if (remove_all)
        ret &= bam_hrecs_remove_line(hrecs, type, type_found);

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return 0;
}

int bam_hdr_remove_lines(bam_hdr_t *bh, const char *type, const char *id, void *h) {
    bam_hrecs_t *hrecs;
    rmhash_t *rh;
    if (!bh || !type || !id || !(rh = (rmhash_t *)h))
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    int itype = (type[0]<<8)|(type[1]);
    khint_t k = kh_get(bam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h)) { // nothing to remove from
        hts_log_warning("Type '%s' does not exist in the header", type);
        return 0;
    }

    bam_hrec_type_t *head = kh_val(hrecs->h, k);
    if (!head) {
        hts_log_error("Header inconsistency");
        return -1;
    }

    int ret = 0;
    bam_hrec_type_t *step = head->next;
    while (step != head) {
        bam_hrec_tag_t *tag = bam_hrecs_find_key(step, id, NULL);
        if (tag && tag->str && tag->len >= 3) {
           k = kh_get(rm, rh, tag->str+3);
           if (k == kh_end(rh)) { // value is not in the hash table, so remove
               bam_hrec_type_t *to_remove = step;
               step = step->next;
               ret |= bam_hrecs_remove_line(hrecs, type, to_remove);
           } else {
               step = step->next;
           }
        } else { // tag is not on the line, so skip to next line
            step = step->next;
        }
    }

    // process the first line
    bam_hrec_tag_t * tag = bam_hrecs_find_key(head, id, NULL);
    if (tag && tag->str && tag->len >= 3) {
       k = kh_get(rm, rh, tag->str+3);
       if (k == kh_end(rh)) { // value is not in the hash table, so remove
           bam_hrec_type_t *to_remove = head;
           head = head->next;
           ret |= bam_hrecs_remove_line(hrecs, type, to_remove);
       }
    }

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

/* ==== Key:val level methods ==== */

const char *bam_hdr_find_tag(bam_hdr_t *bh,
        const char *type,
        const char *ID_key,
        const char *ID_value,
        const char *key) {
    bam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return NULL;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return NULL;
        hrecs = bh->hrecs;
    }

    bam_hrec_type_t *ty = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!ty)
        return NULL;

    bam_hrec_tag_t *tag = bam_hrecs_find_key(ty, key, NULL);
    if (!tag || !tag->str || tag->len < 4)
        return NULL;

    kstring_t ks = KS_INITIALIZER;
    if (kputsn(tag->str+3, tag->len-3, &ks) == EOF) {
        KS_FREE(&ks);
        return NULL;
    }

    return ks_str(&ks);
}

int bam_hdr_remove_tag(bam_hdr_t *bh,
        const char *type,
        const char *ID_key,
        const char *ID_value,
        const char *key) {
    bam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    bam_hrec_type_t *ty = bam_hrecs_find_type(hrecs, type, ID_key, ID_value);
    if (!ty)
        return -1;

    int ret = bam_hrecs_remove_key(hrecs, ty, key);
    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

/*
 * Reconstructs a kstring from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hrecs_rebuild_text(const bam_hrecs_t *hrecs, kstring_t *ks) {
    khint_t k;
    int i;

    ks->l = 0;

    if (!hrecs->h || !hrecs->h->size) {
        return kputsn("", 0, ks) >= 0 ? 0 : -1;
    }

    /* process the array keys first */
    for (i = 0; i < hrecs->type_count; i++) {

        k = kh_get(bam_hrecs_t, hrecs->h, K(hrecs->type_order[i]));
        if (!kh_exist(hrecs->h, k))
            continue;

        if (bam_hrecs_rebuild_lines(hrecs, k, ks))
            return -1;
    }

    /* process the other keys from the hash table */
    for (k = kh_begin(hrecs->h); k != kh_end(hrecs->h); k++) {

        if (!kh_exist(hrecs->h, k))
            continue;

        for (i = 0; i < hrecs->type_count; i++)
            if (kh_key(hrecs->h, k) == K(hrecs->type_order[i]))
                break;
        if (i < hrecs->type_count)
            continue;

        if (bam_hrecs_rebuild_lines(hrecs, k, ks))
            return -1;
    }

    return 0;
}

/*
 * Looks up a reference sequence by name and returns the numerical ID.
 * Returns -1 if unknown reference.
 */
int bam_hdr_name2ref(bam_hdr_t *bh, const char *ref) {
    bam_hrecs_t *hrecs;
    khint_t k;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!hrecs->ref_hash)
        return -1;

    k = kh_get(m_s2i, hrecs->ref_hash, ref);
    return k == kh_end(hrecs->ref_hash) ? bam_name2id(bh, ref) : kh_val(hrecs->ref_hash, k);
}

/*
 * Fixes any PP links in @PG headers.
 * If the entries are in order then this doesn't need doing, but incase
 * our header is out of order this goes through the hrecs->pg[] array
 * setting the prev_id field.
 *
 * Note we can have multiple complete chains. This code should identify the
 * tails of these chains as these are the entries we have to link to in
 * subsequent PP records.
 *
 * Returns 0 on success
 *        -1 on failure (indicating broken PG/PP records)
 */
int bam_hdr_link_pg(bam_hdr_t *bh) {
    bam_hrecs_t *hrecs;
    int i, j, ret = 0, *new_pg_end;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    hrecs->npg_end_alloc = hrecs->npg;
    new_pg_end = realloc(hrecs->pg_end, hrecs->npg * sizeof(*new_pg_end));
    if (!new_pg_end)
        return -1;
    hrecs->pg_end = new_pg_end;

    for (i = 0; i < hrecs->npg; i++)
        hrecs->pg_end[i] = i;

    for (i = 0; i < hrecs->npg; i++) {
        khint_t k;
        bam_hrec_tag_t *tag;
        char tmp;

        for (tag = hrecs->pg[i].tag; tag; tag = tag->next) {
            if (tag->str[0] == 'P' && tag->str[1] == 'P')
                break;
        }
        if (!tag) {
            /* Chain start points */
            continue;
        }

        tmp = tag->str[tag->len]; tag->str[tag->len] = 0;
        k = kh_get(m_s2i, hrecs->pg_hash, tag->str+3);
        tag->str[tag->len] = tmp;

        if (k == kh_end(hrecs->pg_hash)) {
            ret = -1;
            continue;
        }

        hrecs->pg[i].prev_id = hrecs->pg[kh_val(hrecs->pg_hash, k)].id;
        hrecs->pg_end[kh_val(hrecs->pg_hash, k)] = -1;
    }

    for (i = j = 0; i < hrecs->npg; i++) {
        if (hrecs->pg_end[i] != -1)
            hrecs->pg_end[j++] = hrecs->pg_end[i];
    }
    hrecs->npg_end = j;

    /* mark as dirty or empty for rebuild */
    hrecs->dirty = 1;
    redact_header_text(bh);

    return ret;
}

/*
 * Returns a unique ID from a base name.
 *
 * The value returned is valid until the next call to
 * this function.
 */
const char *bam_hdr_pg_id(bam_hdr_t *bh, const char *name) {
    bam_hrecs_t *hrecs;
    if (!bh)
        return NULL;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return NULL;
        hrecs = bh->hrecs;
    }

    khint_t k = kh_get(m_s2i, hrecs->pg_hash, name);
    if (k == kh_end(hrecs->pg_hash))
        return name;

    do {
        sprintf(hrecs->ID_buf, "%.1000s.%d", name, hrecs->ID_cnt++);
        k = kh_get(m_s2i, hrecs->pg_hash, hrecs->ID_buf);
    } while (k != kh_end(hrecs->pg_hash));

    return hrecs->ID_buf;
}

/*
 * Add an @PG line.
 *
 * If we wish complete control over this use bam_hdr_add_line() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple @PG records if we have multiple PG chains.
 *
 * Call it as per bam_hdr_add_line() with a series of key,value pairs ending
 * in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hdr_add_pg(bam_hdr_t *bh, const char *name, ...) {
    bam_hrecs_t *hrecs;
    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (bam_hdr_parse(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    va_list args;

    if (hrecs->npg_end) {
        /* Copy ends array to avoid us looping while modifying it */
        int *end = malloc(hrecs->npg_end * sizeof(int));
        int i, nends = hrecs->npg_end;

        if (!end)
            return -1;

        memcpy(end, hrecs->pg_end, nends * sizeof(*end));

        for (i = 0; i < nends; i++) {
            va_start(args, name);
            if (-1 == bam_hrecs_vadd(hrecs, "PG", args,
                                     "ID", bam_hdr_pg_id(bh, name),
                                     "PN", name,
                                     "PP", hrecs->pg[end[i]].name,
                                     NULL)) {
                free(end);
                return  -1;
            }
            va_end(args);
        }

        free(end);
    } else {
        va_start(args, name);
        if (-1 == bam_hrecs_vadd(hrecs, "PG", args,
                                 "ID", bam_hdr_pg_id(bh, name),
                                 "PN", name,
                                 NULL))
            return -1;
        va_end(args);
    }

    hrecs->dirty = 1;
    redact_header_text(bh);

    return 0;
}

/*! Increments a reference count on bh.
 *
 * This permits multiple files to share the same header, all calling
 * bam_hdr_destroy when done, without causing errors for other open files.
 */
void bam_hdr_incr_ref(bam_hdr_t *bh) {
    if (!bh)
        return;
    bh->ref_count++;
}

/* ==== Internal methods ==== */

/*
 * Creates an empty SAM header.  Allocates space for the SAM header
 * structures (hash tables) ready to be populated.
 *
 * Returns a bam_hrecs_t struct on success (free with bam_hrecs_free())
 *         NULL on failure
 */
bam_hrecs_t *bam_hrecs_new() {
    bam_hrecs_t *hrecs = calloc(1, sizeof(*hrecs));

    if (!hrecs)
        return NULL;

    hrecs->h = kh_init(bam_hrecs_t);
    if (!hrecs->h)
        goto err;

    hrecs->ID_cnt = 1;

    hrecs->nref = 0;
    hrecs->ref  = NULL;
    if (!(hrecs->ref_hash = kh_init(m_s2i)))
        goto err;
    hrecs->refs_changed = -1;

    hrecs->nrg = 0;
    hrecs->rg  = NULL;
    if (!(hrecs->rg_hash = kh_init(m_s2i)))
        goto err;

    hrecs->npg = 0;
    hrecs->pg  = NULL;
    hrecs->npg_end = hrecs->npg_end_alloc = 0;
    hrecs->pg_end = NULL;
    if (!(hrecs->pg_hash = kh_init(m_s2i)))
        goto err;

    if (!(hrecs->tag_pool = pool_create(sizeof(bam_hrec_tag_t))))
        goto err;

    if (!(hrecs->type_pool = pool_create(sizeof(bam_hrec_type_t))))
        goto err;

    if (!(hrecs->str_pool = string_pool_create(8192)))
        goto err;

    if (bam_hrecs_init_type_order(hrecs, NULL))
        goto err;

    return hrecs;

err:
    if (hrecs->h)
        kh_destroy(bam_hrecs_t, hrecs->h);

    if (hrecs->tag_pool)
        pool_destroy(hrecs->tag_pool);

    if (hrecs->type_pool)
        pool_destroy(hrecs->type_pool);

    if (hrecs->str_pool)
        string_pool_destroy(hrecs->str_pool);

    free(hrecs);

    return NULL;
}
#if 0
/*
 * Produces a duplicate copy of source and returns it.
 * Returns NULL on failure
 */
bam_hrecs_t *bam_hrecs_dup(bam_hrecs_t *source) {
        return NULL;
}
#endif
/*! Deallocates all storage used by a bam_hrecs_t struct.
 *
 * This also decrements the header reference count. If after decrementing
 * it is still non-zero then the header is assumed to be in use by another
 * caller and the free is not done.
 *
 */
void bam_hrecs_free(bam_hrecs_t *hrecs) {
    if (!hrecs)
        return;

    if (hrecs->h)
        kh_destroy(bam_hrecs_t, hrecs->h);

    if (hrecs->ref_hash)
        kh_destroy(m_s2i, hrecs->ref_hash);

    if (hrecs->ref)
        free(hrecs->ref);

    if (hrecs->rg_hash)
        kh_destroy(m_s2i, hrecs->rg_hash);

    if (hrecs->rg)
        free(hrecs->rg);

    if (hrecs->pg_hash)
        kh_destroy(m_s2i, hrecs->pg_hash);

    if (hrecs->pg)
        free(hrecs->pg);

    if (hrecs->pg_end)
        free(hrecs->pg_end);

    if (hrecs->type_pool)
        pool_destroy(hrecs->type_pool);

    if (hrecs->tag_pool)
        pool_destroy(hrecs->tag_pool);

    if (hrecs->str_pool)
        string_pool_destroy(hrecs->str_pool);

    if (hrecs->type_order)
        free(hrecs->type_order);

    free(hrecs);
}

/*
 * Internal method already used by the CRAM code
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
bam_hrec_type_t *bam_hrecs_find_type(bam_hrecs_t *hrecs, const char *type,
                                     const char *ID_key, const char *ID_value) {
    if (!hrecs || !type)
        return NULL;
    bam_hrec_type_t *t1, *t2;
    int itype = (type[0]<<8)|(type[1]);
    khint_t k;

    /* Special case for types we have prebuilt hashes on */
    if (ID_key) {
        if (!ID_value)
            return NULL;

        if (type[0]   == 'S' && type[1]   == 'Q' &&
            ID_key[0] == 'S' && ID_key[1] == 'N') {
            k = kh_get(m_s2i, hrecs->ref_hash, ID_value);
            return k != kh_end(hrecs->ref_hash)
                ? hrecs->ref[kh_val(hrecs->ref_hash, k)].ty
                : NULL;
        }

        if (type[0]   == 'R' && type[1]   == 'G' &&
            ID_key[0] == 'I' && ID_key[1] == 'D') {
            k = kh_get(m_s2i, hrecs->rg_hash, ID_value);
            return k != kh_end(hrecs->rg_hash)
                ? hrecs->rg[kh_val(hrecs->rg_hash, k)].ty
                : NULL;
        }

        if (type[0]   == 'P' && type[1]   == 'G' &&
            ID_key[0] == 'I' && ID_key[1] == 'D') {
            k = kh_get(m_s2i, hrecs->pg_hash, ID_value);
            return k != kh_end(hrecs->pg_hash)
                ? hrecs->pg[kh_val(hrecs->pg_hash, k)].ty
                : NULL;
        }
    }

    k = kh_get(bam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h))
        return NULL;

    if (!ID_key)
        return kh_val(hrecs->h, k);

    t1 = t2 = kh_val(hrecs->h, k);
    do {
        bam_hrec_tag_t *tag;
        for (tag = t1->tag; tag; tag = tag->next) {
            if (tag->str[0] == ID_key[0] && tag->str[1] == ID_key[1]) {
                char *cp1 = tag->str+3;
                const char *cp2 = ID_value;
                while (*cp1 && *cp1 == *cp2)
                    cp1++, cp2++;
                if (*cp2 || *cp1)
                    continue;
                return t1;
            }
        }
        t1 = t1->next;
    } while (t1 != t2);

    return NULL;
}

/*
 * Adds or updates tag key,value pairs in a header line.
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line.
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hrecs_update(bam_hrecs_t *hrecs, bam_hrec_type_t *type, va_list ap) {
    if (!hrecs)
        return -1;

    for (;;) {
        char *k, *v;
        bam_hrec_tag_t *tag, *prev = NULL;

        if (!(k = (char *)va_arg(ap, char *)))
            break;
        if (!(v = va_arg(ap, char *)))
            v = "";

        tag = bam_hrecs_find_key(type, k, &prev);
        if (!tag) {
            if (!(tag = pool_alloc(hrecs->tag_pool)))
                return -1;
            if (prev)
                prev->next = tag;
            else
                type->tag = tag;

            tag->next = NULL;
        }

        tag->len = 3 + strlen(v);
        tag->str = string_alloc(hrecs->str_pool, tag->len+1);
        if (!tag->str)
            return -1;

        if (snprintf(tag->str, tag->len+1, "%2.2s:%s", k, v) < 0)
            return -1;
    }

    hrecs->dirty = 1; //mark text as dirty and force a rebuild

    return 0;
}

/*
 * Looks for a specific key in a single sam header line identified by *type.
 * If prev is non-NULL it also fills this out with the previous tag, to
 * permit use in key removal. *prev is set to NULL when the tag is the first
 * key in the list. When a tag isn't found, prev (if non NULL) will be the last
 * tag in the existing list.
 *
 * Returns the tag pointer on success
 *         NULL on failure
 */
bam_hrec_tag_t *bam_hrecs_find_key(bam_hrec_type_t *type,
                                   const char *key,
                                   bam_hrec_tag_t **prev) {
    bam_hrec_tag_t *tag, *p = NULL;
    if (!type)
        return NULL;

    for (tag = type->tag; tag; p = tag, tag = tag->next) {
        if (tag->str[0] == key[0] && tag->str[1] == key[1]) {
            if (prev)
                *prev = p;
            return tag;
        }
    }

    if (prev)
        *prev = p;

    return NULL;
}

int bam_hrecs_remove_key(bam_hrecs_t *hrecs,
                         bam_hrec_type_t *type,
                         const char *key) {
    bam_hrec_tag_t *tag, *prev;
    if (!hrecs)
        return -1;
    tag = bam_hrecs_find_key(type, key, &prev);
    if (!tag)
        return 0; // Not there anyway

    if (!prev) { //first tag
        type->tag = tag->next;
    } else {
        prev->next = tag->next;
    }
    pool_free(hrecs->tag_pool, tag);
    hrecs->dirty = 1; //mark text as dirty and force a rebuild

    return 1;
}

/*
 * Looks up a read-group by name and returns a pointer to the start of the
 * associated tag list.
 *
 * Returns NULL on failure
 */
bam_hrec_rg_t *bam_hrecs_find_rg(bam_hrecs_t *hrecs, const char *rg) {
    khint_t k = kh_get(m_s2i, hrecs->rg_hash, rg);
    return k == kh_end(hrecs->rg_hash)
        ? NULL
        : &hrecs->rg[kh_val(hrecs->rg_hash, k)];
}

void bam_hrecs_dump(bam_hrecs_t *hrecs) {
    khint_t k;
    int i;

    printf("===DUMP===\n");
    for (k = kh_begin(hrecs->h); k != kh_end(hrecs->h); k++) {
        bam_hrec_type_t *t1, *t2;
        char c[2];

        if (!kh_exist(hrecs->h, k))
            continue;

        t1 = t2 = kh_val(hrecs->h, k);
        c[0] = kh_key(hrecs->h, k)>>8;
        c[1] = kh_key(hrecs->h, k)&0xff;
        printf("Type %.2s, count %d\n", c, t1->prev->order+1);

        do {
            bam_hrec_tag_t *tag;
            printf(">>>%d ", t1->order);
            for (tag = t1->tag; tag; tag=tag->next) {
                if (strncmp(c, "CO", 2))
                    printf("\"%.2s\":\"%.*s\"\t", tag->str, tag->len-3, tag->str+3);
                else
                    printf("%s", tag->str);
            }
            putchar('\n');
            t1 = t1->next;
        } while (t1 != t2);
    }

    /* Dump out PG chains */
    printf("\n@PG chains:\n");
    for (i = 0; i < hrecs->npg_end; i++) {
        int j;
        printf("  %d:", i);
        for (j = hrecs->pg_end[i]; j != -1; j = hrecs->pg[j].prev_id) {
            printf("%s%d(%.*s)",
                   j == hrecs->pg_end[i] ? " " : "->",
                   j, hrecs->pg[j].name_len, hrecs->pg[j].name);
        }
        printf("\n");
    }

    puts("===END DUMP===");
}

/*
 * Returns the sort order:
 */
enum sam_sort_order bam_hrecs_sort_order(bam_hrecs_t *hrecs) {
    khint_t k;
    enum sam_sort_order so;

    so = ORDER_UNKNOWN;
    k = kh_get(bam_hrecs_t, hrecs->h, K("HD"));
    if (k != kh_end(hrecs->h)) {
        bam_hrec_type_t *ty = kh_val(hrecs->h, k);
        bam_hrec_tag_t *tag;
        for (tag = ty->tag; tag; tag = tag->next) {
            if (tag->str[0] == 'S' && tag->str[1] == 'O') {
                if (strcmp(tag->str+3, "unsorted") == 0)
                    so = ORDER_UNSORTED;
                else if (strcmp(tag->str+3, "queryname") == 0)
                    so = ORDER_NAME;
                else if (strcmp(tag->str+3, "coordinate") == 0)
                    so = ORDER_COORD;
                else if (strcmp(tag->str+3, "unknown") != 0)
                    hts_log_error("Unknown sort order field: %s", tag->str+3);
            }
        }
    }

    return so;
}

enum sam_group_order bam_hrecs_group_order(bam_hrecs_t *hrecs) {
    khint_t k;
    enum sam_group_order go;

    go = ORDER_NONE;
    k = kh_get(bam_hrecs_t, hrecs->h, K("HD"));
    if (k != kh_end(hrecs->h)) {
        bam_hrec_type_t *ty = kh_val(hrecs->h, k);
        bam_hrec_tag_t *tag;
        for (tag = ty->tag; tag; tag = tag->next) {
            if (tag->str[0] == 'G' && tag->str[1] == 'O') {
                if (strcmp(tag->str+3, "query") == 0)
                    go = ORDER_QUERY;
                else if (strcmp(tag->str+3, "reference") == 0)
                    go = ORDER_REFERENCE;
            }
        }
    }

    return go;
}

// Legacy functions from htslb/cram.h, included here for API compatibility.
typedef bam_hdr_t SAM_hdr;

SAM_hdr *sam_hdr_parse_(const char *hdr, int len) {
    bam_hdr_t *bh = bam_hdr_init();
    if (!bh) return NULL;

    if (bam_hdr_add_lines(bh, hdr, len) != 0) {
        bam_hdr_destroy(bh);
        return NULL;
    }

    return bh;
}

void sam_hdr_free(SAM_hdr *hdr) {
    bam_hdr_destroy(hdr);
}

int sam_hdr_length(SAM_hdr *hdr) {
    return bam_hdr_length(hdr);
}

char *sam_hdr_str(SAM_hdr *hdr) {
    return (char *) bam_hdr_str(hdr);
}
