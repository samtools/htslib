/*
Copyright (c) 2018-2019 Genome Research Ltd.
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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <string.h>
#include <assert.h>
#include <errno.h>
#include "textutils_internal.h"
#include "header.h"

// Hash table for removing multiple lines from the header
KHASH_SET_INIT_STR(rm)
// Used for long refs in SAM files
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

typedef khash_t(rm) rmhash_t;

static int sam_hdr_link_pg(sam_hdr_t *bh);

static int sam_hrecs_vupdate(sam_hrecs_t *hrecs, sam_hrec_type_t *type, va_list ap);
static int sam_hrecs_update(sam_hrecs_t *hrecs, sam_hrec_type_t *type, ...);


#define MAX_ERROR_QUOTE 320 // Prevent over-long error messages
static void sam_hrecs_error(const char *msg, const char *line, size_t len, size_t lno) {
    int j;

    if (len > MAX_ERROR_QUOTE)
        len = MAX_ERROR_QUOTE;
    for (j = 0; j < len && line[j] != '\n'; j++)
        ;
    hts_log_error("%s at line %zd: \"%.*s\"", msg, lno, j, line);
}

/* ==== Static methods ==== */

static int sam_hrecs_init_type_order(sam_hrecs_t *hrecs, char *type_list) {
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

static int sam_hrecs_add_ref_altnames(sam_hrecs_t *hrecs, int nref, const char *list) {
    const char *token;
    ks_tokaux_t aux;

    if (!list)
        return 0;

    for (token = kstrtok(list, ",", &aux); token; token = kstrtok(NULL, NULL, &aux)) {
        if (aux.p == token)
            continue;

        char *name = string_ndup(hrecs->str_pool, token, aux.p - token);
        if (!name)
            return -1;
        int r;
        khint_t k = kh_put(m_s2i, hrecs->ref_hash, name, &r);
        if (r < 0) return -1;

        if (r > 0)
            kh_val(hrecs->ref_hash, k) = nref;
        else if (kh_val(hrecs->ref_hash, k) != nref)
            hts_log_warning("Duplicate entry AN:\"%s\" in sam header", name);
    }

    return 0;
}

static void sam_hrecs_remove_ref_altnames(sam_hrecs_t *hrecs, int expected, const char *list) {
    const char *token, *sn;
    ks_tokaux_t aux;
    kstring_t str = KS_INITIALIZE;

    if (expected < 0 || expected >= hrecs->nref)
        return;
    sn = hrecs->ref[expected].name;

    for (token = kstrtok(list, ",", &aux); token; token = kstrtok(NULL, NULL, &aux)) {
        kputsn(token, aux.p - token, ks_clear(&str));
        khint_t k = kh_get(m_s2i, hrecs->ref_hash, str.s);
        if (k != kh_end(hrecs->ref_hash)
            && kh_val(hrecs->ref_hash, k) == expected
            && strcmp(sn, str.s) != 0)
            kh_del(m_s2i, hrecs->ref_hash, k);
    }

    free(str.s);
}

/* Updates the hash tables in the sam_hrecs_t structure.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int sam_hrecs_update_hashes(sam_hrecs_t *hrecs,
                                   int type,
                                   sam_hrec_type_t *h_type) {
    /* Add to reference hash? */
    if ((type>>8) == 'S' && (type&0xff) == 'Q') {
        sam_hrec_tag_t *tag = h_type->tag;
        int nref = hrecs->nref;
        const char *name = NULL;
        const char *altnames = NULL;
        hts_pos_t len = -1;
        int r;
        khint_t k;

        while (tag) {
            if (tag->str[0] == 'S' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                name = tag->str+3;
            } else if (tag->str[0] == 'L' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                len = strtoll(tag->str+3, NULL, 10);
            } else if (tag->str[0] == 'A' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                altnames = tag->str+3;
            }
            tag = tag->next;
        }

        if (!name) {
            hts_log_error("Header includes @SQ line with no SN: tag");
            return -1; // SN should be present, according to spec.
        }

        if (len == -1) {
            hts_log_error("Header includes @SQ line \"%s\" with no LN: tag",
                          name);
            return -1; // LN should be present, according to spec.
        }

        // Seen already?
        k = kh_get(m_s2i, hrecs->ref_hash, name);
        if (k < kh_end(hrecs->ref_hash)) {
            nref = kh_val(hrecs->ref_hash, k);

            // Check for hash entry added by sam_hrecs_refs_from_targets_array()
            if (hrecs->ref[nref].ty == NULL) {
                // Attach header line to existing stub entry.
                hrecs->ref[nref].ty = h_type;
                // Check lengths match; correct if not.
                if (len != hrecs->ref[nref].len) {
                    char tmp[32];
                    snprintf(tmp, sizeof(tmp), "%" PRIhts_pos,
                             hrecs->ref[nref].len);
                    if (sam_hrecs_update(hrecs, h_type, "LN", tmp, NULL) < 0)
                        return -1;
                }
                if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
                    return -1;

                if (hrecs->refs_changed < 0 || hrecs->refs_changed > nref)
                    hrecs->refs_changed = nref;
                return 0;
            }

            // Check to see if an existing entry is being updated
            if (hrecs->ref[nref].ty == h_type) {
                hrecs->ref[nref].len = len;
                hrecs->ref[nref].name = name;
                if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
                    return -1;

                if (hrecs->refs_changed < 0 || hrecs->refs_changed > nref)
                    hrecs->refs_changed = nref;
                return 0;
            }

            // If here, the name is a duplicate.
            // Check to see if it matches the SN: tag from the earlier record.
            if (strcmp(hrecs->ref[nref].name, name) == 0) {
                hts_log_error("Duplicate entry \"%s\" in sam header",
                                name);
                return -1;
            }

            // Clash with an already-seen altname
            // As SN: should be preferred to AN: add this as a new
            // record and update the hash entry to point to it.
            hts_log_warning("Ref name SN:\"%s\" is a duplicate of an existing AN key", name);
            nref = hrecs->nref;
        }

        if (nref == hrecs->ref_sz) {
            size_t new_sz = hrecs->ref_sz >= 4 ? hrecs->ref_sz + (hrecs->ref_sz / 4) : 32;
            sam_hrec_sq_t *new_ref = realloc(hrecs->ref, sizeof(*hrecs->ref) * new_sz);
            if (!new_ref)
                return -1;
            hrecs->ref = new_ref;
            hrecs->ref_sz = new_sz;
        }

        hrecs->ref[nref].name = name;
        hrecs->ref[nref].len  = len;
        hrecs->ref[nref].ty = h_type;

        k = kh_put(m_s2i, hrecs->ref_hash, hrecs->ref[nref].name, &r);
        if (-1 == r) return -1;
        kh_val(hrecs->ref_hash, k) = nref;

        if (sam_hrecs_add_ref_altnames(hrecs, nref, altnames) < 0)
            return -1;

        if (hrecs->refs_changed < 0 || hrecs->refs_changed > hrecs->nref)
            hrecs->refs_changed = hrecs->nref;
        hrecs->nref++;
    }

    /* Add to read-group hash? */
    if ((type>>8) == 'R' && (type&0xff) == 'G') {
        sam_hrec_tag_t *tag = sam_hrecs_find_key(h_type, "ID", NULL);
        int nrg = hrecs->nrg, r;
        khint_t k;

        if (!tag) {
            hts_log_error("Header includes @RG line with no ID: tag");
            return -1;  // ID should be present, according to spec.
        }
        assert(tag->str && tag->len >= 3);

        // Seen already?
        k = kh_get(m_s2i, hrecs->rg_hash, tag->str + 3);
        if (k < kh_end(hrecs->rg_hash)) {
            nrg = kh_val(hrecs->rg_hash, k);
            assert(hrecs->rg[nrg].ty != NULL);
            if (hrecs->rg[nrg].ty != h_type) {
                hts_log_warning("Duplicate entry \"%s\" in sam header",
                                tag->str + 3);
            } else {
                hrecs->rg[nrg].name = tag->str + 3;
                hrecs->rg[nrg].name_len = tag->len - 3;
            }
            return 0;
        }

        if (nrg == hrecs->rg_sz) {
            size_t new_sz = hrecs->rg_sz >= 4 ? hrecs->rg_sz + hrecs->rg_sz / 4 : 4;
            sam_hrec_rg_t *new_rg = realloc(hrecs->rg, sizeof(*hrecs->rg) * new_sz);
            if (!new_rg)
                return -1;
            hrecs->rg = new_rg;
            hrecs->rg_sz = new_sz;
        }

        hrecs->rg[nrg].name = tag->str + 3;
        hrecs->rg[nrg].name_len = tag->len - 3;
        hrecs->rg[nrg].ty   = h_type;
        hrecs->rg[nrg].id   = nrg;

        k = kh_put(m_s2i, hrecs->rg_hash, hrecs->rg[nrg].name, &r);
        if (-1 == r) return -1;
        kh_val(hrecs->rg_hash, k) = nrg;

        hrecs->nrg++;
    }

    /* Add to program hash? */
    if ((type>>8) == 'P' && (type&0xff) == 'G') {
        sam_hrec_tag_t *tag;
        sam_hrec_pg_t *new_pg;
        int npg = hrecs->npg;

        if (npg == hrecs->pg_sz) {
            size_t new_sz = hrecs->pg_sz >= 4 ? hrecs->pg_sz + hrecs->pg_sz / 4 : 4;
            new_pg = realloc(hrecs->pg, sizeof(*hrecs->pg) * new_sz);
            if (!new_pg)
                return -1;
            hrecs->pg = new_pg;
            hrecs->pg_sz = new_sz;
        }

        tag = h_type->tag;
        hrecs->pg[npg].name = NULL;
        hrecs->pg[npg].name_len = 0;
        hrecs->pg[npg].ty  = h_type;
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
                k = kh_get(m_s2i, hrecs->pg_hash, tag->str+3);

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

static int sam_hrecs_remove_hash_entry(sam_hrecs_t *hrecs, int type, sam_hrec_type_t *h_type) {
    if (!hrecs || !h_type)
        return -1;

    sam_hrec_tag_t *tag;
    const char *key = NULL;
    khint_t k;

    /* Remove name and any alternative names from reference hash */
    if ((type>>8) == 'S' && (type&0xff) == 'Q') {
        const char *altnames = NULL;

        tag = h_type->tag;

        while (tag) {
            if (tag->str[0] == 'S' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                key = tag->str + 3;
            } else if (tag->str[0] == 'A' && tag->str[1] == 'N') {
                assert(tag->len >= 3);
                altnames = tag->str + 3;
            }
            tag = tag->next;
        }

        if (key) {
            k = kh_get(m_s2i, hrecs->ref_hash, key);
            if (k != kh_end(hrecs->ref_hash)) {
                int idx = kh_val(hrecs->ref_hash, k);
                if (idx + 1 < hrecs->nref)
                    memmove(&hrecs->ref[idx], &hrecs->ref[idx+1],
                            sizeof(sam_hrec_sq_t)*(hrecs->nref - idx - 1));
                if (altnames)
                    sam_hrecs_remove_ref_altnames(hrecs, idx, altnames);
                kh_del(m_s2i, hrecs->ref_hash, k);
                hrecs->nref--;
                if (hrecs->refs_changed < 0 || hrecs->refs_changed > idx)
                    hrecs->refs_changed = idx;
                for (k = 0; k < kh_end(hrecs->ref_hash); k++) {
                    if (kh_exist(hrecs->ref_hash, k)
                        && kh_value(hrecs->ref_hash, k) > idx) {
                        kh_value(hrecs->ref_hash, k)--;
                    }
                }
            }
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
                    int idx = kh_val(hrecs->rg_hash, k);
                    if (idx + 1 < hrecs->nrg)
                        memmove(&hrecs->rg[idx], &hrecs->rg[idx+1], sizeof(sam_hrec_rg_t)*(hrecs->nrg - idx - 1));
                    kh_del(m_s2i, hrecs->rg_hash, k);
                    hrecs->nrg--;
                    for (k = 0; k < kh_end(hrecs->rg_hash); k++) {
                        if (kh_exist(hrecs->rg_hash, k)
                            && kh_value(hrecs->rg_hash, k) > idx) {
                            kh_value(hrecs->rg_hash, k)--;
                        }
                    }
                }
                break;
            }
            tag = tag->next;
        }
    }

    return 0;
}

/** Add a header record to the global line ordering
 *
 * If @p after is not NULL, the new record will be inserted after this one,
 * otherwise it will go at the end.
 *
 * An exception is an HD record, which will always be put first unless
 * one is already present.
 */
static void sam_hrecs_global_list_add(sam_hrecs_t *hrecs,
                                      sam_hrec_type_t *h_type,
                                      sam_hrec_type_t *after) {
    const khint32_t hd_type = 'H' << 8 | 'D';
    int update_first_line = 0;

    // First line seen
    if (!hrecs->first_line) {
        hrecs->first_line = h_type->global_next = h_type->global_prev = h_type;
        return;
    }

    // @HD goes at the top (unless there's one already)
    if (h_type->type == hd_type && hrecs->first_line->type != hd_type) {
        after = hrecs->first_line->global_prev;
        update_first_line = 1;
    }

    // If no instructions given, put it at the end
    if (!after)
        after = hrecs->first_line->global_prev;

    h_type->global_prev = after;
    h_type->global_next = after->global_next;
    h_type->global_prev->global_next = h_type;
    h_type->global_next->global_prev = h_type;

    if (update_first_line)
        hrecs->first_line = h_type;
}

/*! Add header record with a va_list interface.
 *
 * Adds a single record to a SAM header.
 *
 * This takes a header record type, a va_list argument and one or more
 * key,value pairs, ending with the NULL key.
 *
 * Eg. sam_hrecs_vadd(h, "SQ", args, "ID", "foo", "LN", "100", NULL).
 *
 * The purpose of the additional va_list parameter is to permit other
 * varargs functions to call this while including their own additional
 * parameters; an example is in sam_hdr_add_pg().
 *
 * Note: this function invokes va_arg at least once, making the value
 * of ap indeterminate after the return. The caller should call
 * va_start/va_end before/after calling this function or use va_copy.
 *
 * @return
 * Returns >= 0 on success;
 *        -1 on failure
 */
static int sam_hrecs_vadd(sam_hrecs_t *hrecs, const char *type, va_list ap, ...) {
    va_list args;
    sam_hrec_type_t *h_type;
    sam_hrec_tag_t *h_tag, *last=NULL;
    int new;
    khint32_t type_i = (type[0]<<8) | type[1], k;

    if (!strncmp(type, "HD", 2) && (h_type = sam_hrecs_find_type_id(hrecs, "HD", NULL, NULL)))
        return sam_hrecs_vupdate(hrecs, h_type, ap);

    if (!(h_type = pool_alloc(hrecs->type_pool)))
        return -1;
    k = kh_put(sam_hrecs_t, hrecs->h, type_i, &new);
    if (new < 0)
        return -1;

    h_type->type = type_i;

    // Form the ring, either with self or other lines of this type
    if (!new) {
        sam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
        p = t->prev;

        assert(p->next == t);
        p->next = h_type;
        h_type->prev = p;

        t->prev = h_type;
        h_type->next = t;
    } else {
        kh_val(hrecs->h, k) = h_type;
        h_type->prev = h_type->next = h_type;
    }
    h_type->tag = NULL;

    // Add to global line ordering after any existing line of the same type,
    // or at the end if no line of this type exists yet.
    sam_hrecs_global_list_add(hrecs, h_type, !new ? h_type->prev : NULL);

    // Check linked-list invariants
    assert(h_type->prev->next == h_type);
    assert(h_type->next->prev == h_type);
    assert(h_type->global_prev->global_next == h_type);
    assert(h_type->global_next->global_prev == h_type);

    // Any ... varargs
    va_start(args, ap);
    for (;;) {
        char *key, *val = NULL, *str;

        if (!(key = (char *)va_arg(args, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(val = (char *)va_arg(args, char *)))
            break;
        if (*val == '\0')
            continue;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(val);
            str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!str || snprintf(str, h_tag->len+1, "%2.2s:%s", key, val) < 0)
                return -1;
            h_tag->str = str;
        } else {
            h_tag->len = strlen(key);
            h_tag->str = string_ndup(hrecs->str_pool, key, h_tag->len);
            if (!h_tag->str)
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
        char *key, *val = NULL, *str;

        if (!(key = (char *)va_arg(ap, char *)))
            break;
        if (strncmp(type, "CO", 2) && !(val = (char *)va_arg(ap, char *)))
            break;

        if (!(h_tag = pool_alloc(hrecs->tag_pool)))
            return -1;

        if (strncmp(type, "CO", 2)) {
            h_tag->len = 3 + strlen(val);
            str = string_alloc(hrecs->str_pool, h_tag->len+1);
            if (!str || snprintf(str, h_tag->len+1, "%2.2s:%s", key, val) < 0)
                return -1;
            h_tag->str = str;
        } else {
            h_tag->len = strlen(key);
            h_tag->str = string_ndup(hrecs->str_pool, key, h_tag->len);
            if (!h_tag->str)
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
    if (-1 == sam_hrecs_update_hashes(hrecs, itype, h_type))
        return -1;

    if (!strncmp(type, "PG", 2))
        hrecs->pgs_changed = 1;

    hrecs->dirty = 1;

    return 0;
}

// As sam_hrecs_vadd(), but without the extra va_list parameter
static int sam_hrecs_add(sam_hrecs_t *hrecs, const char *type, ...) {
    va_list args;
    int res;
    va_start(args, type);
    res = sam_hrecs_vadd(hrecs, type, args, NULL);
    va_end(args);
    return res;
}

/*
 * Function for deallocating a list of tags
 */

static void sam_hrecs_free_tags(sam_hrecs_t *hrecs, sam_hrec_tag_t *tag) {
    if (!hrecs || !tag)
        return;
    if (tag->next)
        sam_hrecs_free_tags(hrecs, tag->next);

    pool_free(hrecs->tag_pool, tag);
}

static int sam_hrecs_remove_line(sam_hrecs_t *hrecs, const char *type_name, sam_hrec_type_t *type_found) {
    if (!hrecs || !type_name || !type_found)
        return -1;

    int itype = (type_name[0]<<8) | type_name[1];
    khint_t k = kh_get(sam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h))
        return -1;

    // Remove from global list (remembering it could be the only line)
    if (hrecs->first_line == type_found) {
        hrecs->first_line = (type_found->global_next != type_found
                             ? type_found->global_next : NULL);
    }
    type_found->global_next->global_prev = type_found->global_prev;
    type_found->global_prev->global_next = type_found->global_next;

    /* single element in the list */
    if (type_found->prev == type_found || type_found->next == type_found) {
        kh_del(sam_hrecs_t, hrecs->h, k);
    } else {
        type_found->prev->next = type_found->next;
        type_found->next->prev = type_found->prev;
        if (kh_val(hrecs->h, k) == type_found) { //first element
            kh_val(hrecs->h, k) = type_found->next;
        }
    }

    if (!strncmp(type_name, "SQ", 2) || !strncmp(type_name, "RG", 2))
        sam_hrecs_remove_hash_entry(hrecs, itype, type_found);

    sam_hrecs_free_tags(hrecs, type_found->tag);
    pool_free(hrecs->type_pool, type_found);

    hrecs->dirty = 1;

    return 0;
}

// Paste together a line from the parsed data structures
static int build_header_line(const sam_hrec_type_t *ty, kstring_t *ks) {
    sam_hrec_tag_t *tag;
    int r = 0;
    char c[2]= { ty->type >> 8, ty->type & 0xff };

    r |= (kputc_('@', ks) == EOF);
    r |= (kputsn(c, 2, ks) == EOF);
    for (tag = ty->tag; tag; tag = tag->next) {
        r |= (kputc_('\t', ks) == EOF);
        r |= (kputsn(tag->str, tag->len, ks) == EOF);
    }

    return r;
}

static int sam_hrecs_rebuild_lines(const sam_hrecs_t *hrecs, kstring_t *ks) {
    const sam_hrec_type_t *t1, *t2;

    if (!hrecs->first_line)
        return kputsn("", 0, ks) >= 0 ? 0 : -1;

    t1 = t2 = hrecs->first_line;
    do {
        if (build_header_line(t1, ks) != 0)
            return -1;
        if (kputc('\n', ks) < 0)
            return -1;

        t1 = t1->global_next;
    } while (t1 != t2);

    return 0;
}

static int sam_hrecs_parse_lines(sam_hrecs_t *hrecs, const char *hdr, size_t len) {
    size_t i, lno;

    if (!hrecs || len > SSIZE_MAX)
        return -1;

    if (!len)
        len = strlen(hdr);

    if (len < 3) {
        if (len == 0 || *hdr == '\0') return 0;
        sam_hrecs_error("Header line too short", hdr, len, 1);
        return -1;
    }

    for (i = 0, lno = 1; i < len - 3 && hdr[i] != '\0'; i++, lno++) {
        khint32_t type;
        khint_t k;

        int l_start = i, new;
        sam_hrec_type_t *h_type;
        sam_hrec_tag_t *h_tag, *last;

        if (hdr[i] != '@') {
            sam_hrecs_error("Header line does not start with '@'",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }

        type = (((uint8_t) hdr[i+1])<<8) | (uint8_t) hdr[i+2];
        if (!isalpha_c(hdr[i+1]) || !isalpha_c(hdr[i+2])) {
            sam_hrecs_error("Header line does not have a two character key",
                          &hdr[l_start], len - l_start, lno);
            return -1;
        }

        i += 3;
        if (i == len || hdr[i] == '\n')
            continue;

        // Add the header line type
        if (!(h_type = pool_alloc(hrecs->type_pool)))
            return -1;
        k = kh_put(sam_hrecs_t, hrecs->h, type, &new);
        if (new < 0)
            return -1;

        h_type->type = type;

        // Add to end of global list
        sam_hrecs_global_list_add(hrecs, h_type, NULL);

        // Form the ring, either with self or other lines of this type
        if (!new) {
            sam_hrec_type_t *t = kh_val(hrecs->h, k), *p;
            p = t->prev;

            assert(p->next == t);
            p->next = h_type;
            h_type->prev = p;

            t->prev = h_type;
            h_type->next = t;
        } else {
            kh_val(hrecs->h, k) = h_type;
            h_type->prev = h_type->next = h_type;
        }

        // Parse the tags on this line
        last = NULL;
        if ((type>>8) == 'C' && (type&0xff) == 'O') {
            size_t j;

            if (i == len || hdr[i] != '\t') {
                sam_hrecs_error("Missing tab",
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
            do {
                size_t j;

                if (i == len || hdr[i] != '\t') {
                    sam_hrecs_error("Missing tab",
                                  &hdr[l_start], len - l_start, lno);
                    return -1;
                }

                for (j = ++i; j < len && hdr[j] != '\0' && hdr[j] != '\n' && hdr[j] != '\t'; j++)
                    ;

                if (j - i < 3 || hdr[i + 2] != ':') {
                    sam_hrecs_error("Malformed key:value pair",
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
        if (-1 == sam_hrecs_update_hashes(hrecs, type, h_type))
            return -1;
    }

    return 0;
}

/*! Update sam_hdr_t target_name and target_len arrays
 *
 *  @return 0 on success; -1 on failure
 */
int sam_hdr_update_target_arrays(sam_hdr_t *bh, const sam_hrecs_t *hrecs,
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
    khint_t k;
    khash_t(s2i) *long_refs = (khash_t(s2i) *) bh->sdict;
    for (i = refs_changed; i < hrecs->nref; i++) {
        if (i >= bh->n_targets
            || strcmp(bh->target_name[i], hrecs->ref[i].name) != 0) {
            if (i < bh->n_targets)
                free(bh->target_name[i]);
            bh->target_name[i] = strdup(hrecs->ref[i].name);
            if (!bh->target_name[i])
                return -1;
        }
        if (hrecs->ref[i].len < UINT32_MAX) {
            bh->target_len[i] = hrecs->ref[i].len;

            if (!long_refs)
                continue;

            // Check if we have an old length, if so remove it.
            k = kh_get(s2i, long_refs, bh->target_name[i]);
            if (k < kh_end(long_refs))
                kh_del(s2i, long_refs, k);
        } else {
            bh->target_len[i] = UINT32_MAX;
            if (bh->hrecs != hrecs) {
                // Called from sam_hdr_dup; need to add sdict entries
                if (!long_refs) {
                    if (!(bh->sdict = long_refs = kh_init(s2i)))
                        return -1;
                }

                // Add / update length
                int absent;
                k = kh_put(s2i, long_refs, bh->target_name[i], &absent);
                if (absent < 0)
                    return -1;
                kh_val(long_refs, k) = hrecs->ref[i].len;
            }
        }
    }

    // Free up any names that have been removed
    for (; i < bh->n_targets; i++) {
        if (long_refs) {
            k = kh_get(s2i, long_refs, bh->target_name[i]);
            if (k < kh_end(long_refs))
                kh_del(s2i, long_refs, k);
        }
        free(bh->target_name[i]);
    }

    bh->n_targets = hrecs->nref;
    return 0;
}

static int rebuild_target_arrays(sam_hdr_t *bh) {
    if (!bh || !bh->hrecs)
        return -1;

    sam_hrecs_t *hrecs = bh->hrecs;
    if (hrecs->refs_changed < 0)
        return 0;

    if (sam_hdr_update_target_arrays(bh, hrecs, hrecs->refs_changed) != 0)
        return -1;

    hrecs->refs_changed = -1;
    return 0;
}

/// Populate hrecs refs array from header target_name, target_len arrays
/**
 * @return 0 on success; -1 on failure
 *
 * Pre-fills the refs hash from the target arrays.  For BAM files this
 * will ensure that they are in the correct order as the target arrays
 * are the canonical source for converting target ids to names and lengths.
 *
 * The added entries do not link to a header line. sam_hrecs_update_hashes()
 * will add the links later for lines found in the text header.
 *
 * This should be called before the text header is parsed.
 */
static int sam_hrecs_refs_from_targets_array(sam_hrecs_t *hrecs,
                                             const sam_hdr_t *bh) {
    int32_t tid = 0;

    if (!hrecs || !bh)
        return -1;

    // This should always be called before parsing the text header
    // so the ref array should start off empty, and we don't have to try
    // to reconcile any existing data.
    if (hrecs->nref > 0) {
        hts_log_error("Called with non-empty ref array");
        return -1;
    }

    if (hrecs->ref_sz < bh->n_targets) {
        sam_hrec_sq_t *new_ref = realloc(hrecs->ref,
                                         bh->n_targets * sizeof(*new_ref));
        if (!new_ref)
            return -1;

        hrecs->ref = new_ref;
        hrecs->ref_sz = bh->n_targets;
    }

    for (tid = 0; tid < bh->n_targets; tid++) {
        khint_t k;
        int r;
        hrecs->ref[tid].name = string_dup(hrecs->str_pool, bh->target_name[tid]);
        if (!hrecs->ref[tid].name) goto fail;
        if (bh->target_len[tid] < UINT32_MAX || !bh->sdict) {
            hrecs->ref[tid].len  = bh->target_len[tid];
        } else {
            khash_t(s2i) *long_refs = (khash_t(s2i) *) bh->sdict;
            k = kh_get(s2i, long_refs, hrecs->ref[tid].name);
            if (k < kh_end(long_refs)) {
                hrecs->ref[tid].len = kh_val(long_refs, k);
            } else {
                hrecs->ref[tid].len = UINT32_MAX;
            }
        }
        hrecs->ref[tid].ty   = NULL;
        k = kh_put(m_s2i, hrecs->ref_hash, hrecs->ref[tid].name, &r);
        if (r < 0) goto fail;
        if (r == 0) {
            hts_log_error("Duplicate entry \"%s\" in target list",
                            hrecs->ref[tid].name);
            return -1;
        } else {
            kh_val(hrecs->ref_hash, k) = tid;
        }
    }
    hrecs->nref = bh->n_targets;
    return 0;

 fail: {
        int32_t i;
        hts_log_error("%s", strerror(errno));
        for (i = 0; i < tid; i++) {
            khint_t k;
            if (!hrecs->ref[i].name) continue;
            k = kh_get(m_s2i, hrecs->ref_hash, hrecs->ref[tid].name);
            if (k < kh_end(hrecs->ref_hash)) kh_del(m_s2i, hrecs->ref_hash, k);
        }
        hrecs->nref = 0;
        return -1;
    }
}

/*
 * Add SQ header records for any references in the hrecs->ref array that
 * were added by sam_hrecs_refs_from_targets_array() but have not
 * been linked to an @SQ line by sam_hrecs_update_hashes() yet.
 *
 * This may be needed either because:
 *
 *   - A bam file was read that had entries in its refs list with no
 *     corresponding @SQ line.
 *
 *   - A program constructed a sam_hdr_t which has target_name and target_len
 *     array entries with no corresponding @SQ line in text.
 */
static int add_stub_ref_sq_lines(sam_hrecs_t *hrecs) {
    int tid;
    char len[32];

    for (tid = 0; tid < hrecs->nref; tid++) {
        if (hrecs->ref[tid].ty == NULL) {
            snprintf(len, sizeof(len), "%"PRIhts_pos, hrecs->ref[tid].len);
            if (sam_hrecs_add(hrecs, "SQ",
                              "SN", hrecs->ref[tid].name,
                              "LN", len, NULL) != 0)
                return -1;

            // Check that the stub has actually been filled
            if(hrecs->ref[tid].ty == NULL) {
                hts_log_error("Reference stub with tid=%d, name=\"%s\", len=%"PRIhts_pos" could not be filled",
                        tid, hrecs->ref[tid].name, hrecs->ref[tid].len);
                return -1;
            }
        }
    }
    return 0;
}

int sam_hdr_fill_hrecs(sam_hdr_t *bh) {
    sam_hrecs_t *hrecs = sam_hrecs_new();

    if (!hrecs)
        return -1;

    if (bh->target_name && bh->target_len && bh->n_targets > 0) {
        if (sam_hrecs_refs_from_targets_array(hrecs, bh) != 0) {
            sam_hrecs_free(hrecs);
            return -1;
        }
    }

    // Parse existing header text
    if (bh->text && bh->l_text > 0) {
        if (sam_hrecs_parse_lines(hrecs, bh->text, bh->l_text) != 0) {
            sam_hrecs_free(hrecs);
            return -1;
        }
    }

    if (add_stub_ref_sq_lines(hrecs) < 0) {
        sam_hrecs_free(hrecs);
        return -1;
    }

    bh->hrecs = hrecs;

    if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
        return -1;

    return 0;
}

/** Remove outdated header text

    @param bh     BAM header

    This is called when API functions have changed the header so that the
    text version is no longer valid.
 */
static void redact_header_text(sam_hdr_t *bh) {
    assert(bh->hrecs && bh->hrecs->dirty);
    bh->l_text = 0;
    free(bh->text);
    bh->text = NULL;
}

/** Find nth header record of a given type

    @param type   Header type (SQ, RG etc.)
    @param idx    0-based index

    @return sam_hrec_type_t pointer to the record on success
            NULL if no record exists with the given type and index
 */

static sam_hrec_type_t *sam_hrecs_find_type_pos(sam_hrecs_t *hrecs,
                                                const char *type, int idx) {
    sam_hrec_type_t *first, *itr;

    if (idx < 0)
        return NULL;

    if (type[0] == 'S' && type[1] == 'Q')
        return idx < hrecs->nref ? hrecs->ref[idx].ty : NULL;

    if (type[0] == 'R' && type[1] == 'G')
        return idx < hrecs->nrg ? hrecs->rg[idx].ty : NULL;

    if (type[0] == 'P' && type[1] == 'G')
        return idx < hrecs->npg ? hrecs->pg[idx].ty : NULL;

    first = itr = sam_hrecs_find_type_id(hrecs, type, NULL, NULL);
    if (!first)
        return NULL;

    while (idx > 0) {
        itr = itr->next;
        if (itr == first)
            break;
        --idx;
    }

    return idx == 0 ? itr : NULL;
}

/* ==== Public methods ==== */

size_t sam_hdr_length(sam_hdr_t *bh) {
    if (!bh || -1 == sam_hdr_rebuild(bh))
        return SIZE_MAX;

    return bh->l_text;
}

const char *sam_hdr_str(sam_hdr_t *bh) {
    if (!bh || -1 == sam_hdr_rebuild(bh))
        return NULL;

    return bh->text;
}

int sam_hdr_nref(const sam_hdr_t *bh) {
    if (!bh)
        return -1;

    return bh->hrecs ? bh->hrecs->nref : bh->n_targets;
}

/*
 * Reconstructs the text representation from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_rebuild(sam_hdr_t *bh) {
    sam_hrecs_t *hrecs;
    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs))
        return bh->text ? 0 : -1;

    if (hrecs->refs_changed >= 0) {
        if (rebuild_target_arrays(bh) < 0) {
            hts_log_error("Header target array rebuild has failed");
            return -1;
        }
    }

    /* If header text wasn't changed or header is empty, don't rebuild it. */
    if (!hrecs->dirty)
        return 0;

    if (hrecs->pgs_changed)
        sam_hdr_link_pg(bh);

    kstring_t ks = KS_INITIALIZE;
    if (sam_hrecs_rebuild_text(hrecs, &ks) != 0) {
        ks_free(&ks);
        hts_log_error("Header text rebuild has failed");
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
int sam_hdr_add_lines(sam_hdr_t *bh, const char *lines, size_t len) {
    sam_hrecs_t *hrecs;

    if (!bh || !lines)
        return -1;

    if (len == 0 && *lines == '\0')
        return 0;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (sam_hrecs_parse_lines(hrecs, lines, len) != 0)
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
 * Eg. sam_hdr_add_line(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_add_line(sam_hdr_t *bh, const char *type, ...) {
    va_list args;
    sam_hrecs_t *hrecs;

    if (!bh || !type)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    va_start(args, type);
    int ret = sam_hrecs_vadd(hrecs, type, args, NULL);
    va_end(args);

    if (ret == 0) {
        if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
            return -1;

        if (hrecs->dirty)
            redact_header_text(bh);
    }

    return ret;
}

/*
 * Returns a complete line of formatted text for a specific head type/ID
 * combination. If ID_key is NULL then it returns the first line of the specified
 * type.
 */
int sam_hdr_find_line_id(sam_hdr_t *bh, const char *type,
                      const char *ID_key, const char *ID_val, kstring_t *ks) {
    sam_hrecs_t *hrecs;
    if (!bh || !type)
        return -2;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *ty = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_val);
    if (!ty)
        return -1;

    ks->l = 0;
    if (build_header_line(ty, ks) < 0) {
        return -2;
    }

    return 0;
}

int sam_hdr_find_line_pos(sam_hdr_t *bh, const char *type,
                          int pos, kstring_t *ks) {
    sam_hrecs_t *hrecs;
    if (!bh || !type)
        return -2;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *ty = sam_hrecs_find_type_pos(hrecs, type, pos);
    if (!ty)
        return -1;

    ks->l = 0;
    if (build_header_line(ty, ks) < 0) {
        return -2;
    }

    return 0;
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

int sam_hdr_remove_line_id(sam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value) {
    sam_hrecs_t *hrecs;
    if (!bh || !type)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!strncmp(type, "PG", 2)) {
        hts_log_warning("Removing PG lines is not supported!");
        return -1;
    }

    sam_hrec_type_t *type_found = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_value);
    if (!type_found)
        return 0;

    int ret = sam_hrecs_remove_line(hrecs, type, type_found);
    if (ret == 0) {
        if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
            return -1;

        if (hrecs->dirty)
            redact_header_text(bh);
    }

    return ret;
}

/*
 * Remove a line from the header by specifying the position in the type
 * group, i.e. 3rd @SQ line.
 *
 * Returns 0 on success and -1 on error
 */

int sam_hdr_remove_line_pos(sam_hdr_t *bh, const char *type, int position) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || position <= 0)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!strncmp(type, "PG", 2)) {
        hts_log_warning("Removing PG lines is not supported!");
        return -1;
    }

    sam_hrec_type_t *type_found = sam_hrecs_find_type_pos(hrecs, type,
                                                          position);
    if (!type_found)
        return -1;

    int ret = sam_hrecs_remove_line(hrecs, type, type_found);
    if (ret == 0) {
        if (hrecs->refs_changed >= 0 && rebuild_target_arrays(bh) != 0)
            return -1;

        if (hrecs->dirty)
            redact_header_text(bh);
    }

    return ret;
}

/*
 * Check if sam_hdr_update_line() is being used to change the name of
 * a record, and if the new name is going to clash with an existing one.
 *
 * If ap includes repeated keys, we go with the last one as sam_hrecs_vupdate()
 * will go through them all and leave the final one in place.
 *
 * Returns 0 if the name does not change
 *         1 if the name changes but does not clash
 *        -1 if the name changes and the new one is already in use
 */
static int check_for_name_update(sam_hrecs_t *hrecs, sam_hrec_type_t *rec,
                                 va_list ap, const char **old_name,
                                 const char **new_name,
                                 char id_tag_out[3],
                                 khash_t(m_s2i) **hash_out) {
    char *key, *val;
    const char *id_tag;
    sam_hrec_tag_t *tag, *prev;
    khash_t(m_s2i) *hash;
    khint_t k;
    int ret = 0;

    if        (rec->type == TYPEKEY("SQ")) {
        id_tag = "SN"; hash = hrecs->ref_hash;
    } else if (rec->type == TYPEKEY("RG")) {
        id_tag = "ID"; hash = hrecs->rg_hash;
    } else if (rec->type == TYPEKEY("PG")) {
        id_tag = "ID"; hash = hrecs->pg_hash;
    } else {
        return 0;
    }

    memcpy(id_tag_out, id_tag, 3);
    *hash_out = hash;

    tag = sam_hrecs_find_key(rec, id_tag, &prev);
    if (!tag)
        return 0;
    assert(tag->len >= 3);
    *old_name = tag->str + 3;

    while ((key = va_arg(ap, char *)) != NULL) {
        val = va_arg(ap, char *);
        if (!val) val = "";
        if (strcmp(key, id_tag) != 0) continue;
        if (strcmp(val, tag->str + 3) == 0) { ret = 0; continue; }
        k = kh_get(m_s2i, hash, val);
        ret = k < kh_end(hash) ? -1 : 1;
        *new_name = val;
    }
    return ret;
}

int sam_hdr_update_line(sam_hdr_t *bh, const char *type,
        const char *ID_key, const char *ID_value, ...) {
    sam_hrecs_t *hrecs;
    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    int ret, rename;
    sam_hrec_type_t *ty = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_value);
    if (!ty)
        return -1;

    va_list args;
    const char *old_name = "?", *new_name = "?";
    char id_tag[3];
    khash_t(m_s2i) *hash = NULL;
    va_start(args, ID_value);
    rename = check_for_name_update(hrecs, ty, args,
                                   &old_name, &new_name, id_tag, &hash);
    va_end(args);
    if (rename < 0) {
        hts_log_error("Cannot rename @%s \"%s\" to \"%s\" : already exists",
                      type, old_name, new_name);
        return -1;
    }
    if (rename > 0 && TYPEKEY(type) == TYPEKEY("PG")) {
        // This is just too complicated
        hts_log_error("Renaming @PG records is not supported");
        return -1;
    }
    va_start(args, ID_value);
    ret = sam_hrecs_vupdate(hrecs, ty, args);
    va_end(args);

    if (ret)
        return ret;

    // TODO Account for @SQ-AN altnames

    if (rename) {
        // Adjust the hash table to point to the new name
        // sam_hrecs_update_hashes() should sort out everything else
        khint_t k = kh_get(m_s2i, hash, old_name);
        sam_hrec_tag_t *new_tag = sam_hrecs_find_key(ty, id_tag, NULL);
        int r, pos;
        assert(k < kh_end(hash));        // Or we wouldn't have found it earlier
        assert(new_tag && new_tag->str); // id_tag should exist
        assert(new_tag->len > 3);
        pos = kh_val(hash, k);
        kh_del(m_s2i, hash, k);
        k = kh_put(m_s2i, hash, new_tag->str + 3, &r);
        if (r < 1) {
            hts_log_error("Failed to rename item in hash table");
            return -1;
        }
        kh_val(hash, k) = pos;
    }

    ret = sam_hrecs_update_hashes(hrecs, TYPEKEY(type), ty);

    if (!ret && hrecs->refs_changed >= 0)
        ret = rebuild_target_arrays(bh);

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

int sam_hdr_remove_except(sam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value) {
    sam_hrecs_t *hrecs;
    if (!bh || !type)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *step;
    int ret = 1, remove_all = (ID_key == NULL);

    if (!strncmp(type, "PG", 2) || !strncmp(type, "CO", 2)) {
        hts_log_warning("Removing PG or CO lines is not supported!");
        return -1;
    }

    sam_hrec_type_t *type_found = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_value);
    if (!type_found) { // remove all line of this type
        int itype = (type[0]<<8)|(type[1]);
        khint_t k = kh_get(sam_hrecs_t, hrecs->h, itype);
        if (k == kh_end(hrecs->h))
            return 0;
        type_found =  kh_val(hrecs->h, k);
        if (!type_found)
            return 0;
        remove_all = 1;
    }

    step = type_found->next;
    while (step != type_found) {
        sam_hrec_type_t *to_remove = step;
        step = step->next;
        ret &= sam_hrecs_remove_line(hrecs, type, to_remove);
    }

    if (remove_all)
        ret &= sam_hrecs_remove_line(hrecs, type, type_found);

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return 0;
}

int sam_hdr_remove_lines(sam_hdr_t *bh, const char *type, const char *id, void *h) {
    sam_hrecs_t *hrecs;
    rmhash_t *rh = (rmhash_t *)h;

    if (!bh || !type)
        return -1;
    if (!rh) // remove all lines
        return sam_hdr_remove_except(bh, type, NULL, NULL);
    if (!id)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    int itype = (type[0]<<8)|(type[1]);
    khint_t k = kh_get(sam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h)) // nothing to remove from
        return 0;

    sam_hrec_type_t *head = kh_val(hrecs->h, k);
    if (!head) {
        hts_log_error("Header inconsistency");
        return -1;
    }

    int ret = 0;
    sam_hrec_type_t *step = head->next;
    while (step != head) {
        sam_hrec_tag_t *tag = sam_hrecs_find_key(step, id, NULL);
        if (tag && tag->str && tag->len >= 3) {
           k = kh_get(rm, rh, tag->str+3);
           if (k == kh_end(rh)) { // value is not in the hash table, so remove
               sam_hrec_type_t *to_remove = step;
               step = step->next;
               ret |= sam_hrecs_remove_line(hrecs, type, to_remove);
           } else {
               step = step->next;
           }
        } else { // tag is not on the line, so skip to next line
            step = step->next;
        }
    }

    // process the first line
    sam_hrec_tag_t * tag = sam_hrecs_find_key(head, id, NULL);
    if (tag && tag->str && tag->len >= 3) {
       k = kh_get(rm, rh, tag->str+3);
       if (k == kh_end(rh)) { // value is not in the hash table, so remove
           sam_hrec_type_t *to_remove = head;
           head = head->next;
           ret |= sam_hrecs_remove_line(hrecs, type, to_remove);
       }
    }

    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

int sam_hdr_count_lines(sam_hdr_t *bh, const char *type) {
    int count;
    sam_hrec_type_t *first_ty, *itr_ty;

    if (!bh || !type)
        return -1;

    if (!bh->hrecs) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
    }

    // Deal with types that have counts
    switch (type[0]) {
    case 'S':
        if (type[1] == 'Q')
            return bh->hrecs->nref;
        break;
    case 'R':
        if (type[1] == 'G')
            return bh->hrecs->nrg;
        break;
    case 'P':
        if (type[1] == 'G')
            return bh->hrecs->npg;
        break;
    default:
        break;
    }

    first_ty = sam_hrecs_find_type_id(bh->hrecs, type, NULL, NULL);
    if (!first_ty)
        return 0;

    count = 1;
    for (itr_ty = first_ty->next;
         itr_ty && itr_ty != first_ty; itr_ty = itr_ty->next) {
        count++;
    }

    return count;
}

int sam_hdr_line_index(sam_hdr_t *bh,
                       const char *type,
                       const char *key) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return -2;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    khint_t k;
    int idx = -1;
    switch (type[0]) {
    case 'S':
        if (type[1] == 'Q') {
            k = kh_get(m_s2i, hrecs->ref_hash, key);
            if (k != kh_end(hrecs->ref_hash))
                idx = kh_val(hrecs->ref_hash, k);
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    case 'R':
        if (type[1] == 'G') {
            k = kh_get(m_s2i, hrecs->rg_hash, key);
            if (k != kh_end(hrecs->rg_hash))
                idx = kh_val(hrecs->rg_hash, k);
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    case 'P':
        if (type[1] == 'G') {
            k = kh_get(m_s2i, hrecs->pg_hash, key);
            if (k != kh_end(hrecs->pg_hash))
                idx = kh_val(hrecs->pg_hash, k);
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    default:
        hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
    }

    return idx;
}

const char *sam_hdr_line_name(sam_hdr_t *bh,
                              const char *type,
                              int pos) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || pos < 0)
        return NULL;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return NULL;
        hrecs = bh->hrecs;
    }

    switch (type[0]) {
    case 'S':
        if (type[1] == 'Q') {
            if (pos < hrecs->nref)
                return hrecs->ref[pos].name;
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    case 'R':
        if (type[1] == 'G') {
            if (pos < hrecs->nrg)
                return hrecs->rg[pos].name;
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    case 'P':
        if (type[1] == 'G') {
            if (pos < hrecs->npg)
                return hrecs->pg[pos].name;
        } else {
            hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
        }
        break;
    default:
        hts_log_warning("Type '%s' not supported. Only @SQ, @RG and @PG lines are indexed", type);
    }

    return NULL;
}

/* ==== Key:val level methods ==== */

int sam_hdr_find_tag_id(sam_hdr_t *bh,
                     const char *type,
                     const char *ID_key,
                     const char *ID_value,
                     const char *key,
                     kstring_t *ks) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return -2;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *ty = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_value);
    if (!ty)
        return -1;

    sam_hrec_tag_t *tag = sam_hrecs_find_key(ty, key, NULL);
    if (!tag || !tag->str || tag->len < 4)
        return -1;

    ks->l = 0;
    if (kputsn(tag->str+3, tag->len-3, ks) == EOF) {
        return -2;
    }

    return 0;
}

int sam_hdr_find_tag_pos(sam_hdr_t *bh,
                     const char *type,
                     int pos,
                     const char *key,
                     kstring_t *ks) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return -2;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *ty = sam_hrecs_find_type_pos(hrecs, type, pos);
    if (!ty)
        return -1;

    sam_hrec_tag_t *tag = sam_hrecs_find_key(ty, key, NULL);
    if (!tag || !tag->str || tag->len < 4)
        return -1;

    ks->l = 0;
    if (kputsn(tag->str+3, tag->len-3, ks) == EOF) {
        return -2;
    }

    return 0;
}

int sam_hdr_remove_tag_id(sam_hdr_t *bh,
        const char *type,
        const char *ID_key,
        const char *ID_value,
        const char *key) {
    sam_hrecs_t *hrecs;
    if (!bh || !type || !key)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    sam_hrec_type_t *ty = sam_hrecs_find_type_id(hrecs, type, ID_key, ID_value);
    if (!ty)
        return -1;

    int ret = sam_hrecs_remove_key(hrecs, ty, key);
    if (!ret && hrecs->dirty)
        redact_header_text(bh);

    return ret;
}

/*
 * Reconstructs a kstring from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hrecs_rebuild_text(const sam_hrecs_t *hrecs, kstring_t *ks) {
    ks->l = 0;

    if (!hrecs->h || !hrecs->h->size) {
        return kputsn("", 0, ks) >= 0 ? 0 : -1;
    }
    if (sam_hrecs_rebuild_lines(hrecs, ks) != 0)
        return -1;

    return 0;
}

/*
 * Looks up a reference sequence by name and returns the numerical ID.
 * Returns -1 if unknown reference; -2 if header could not be parsed.
 */
int sam_hdr_name2tid(sam_hdr_t *bh, const char *ref) {
    sam_hrecs_t *hrecs;
    khint_t k;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -2;
        hrecs = bh->hrecs;
    }

    if (!hrecs->ref_hash)
        return -1;

    k = kh_get(m_s2i, hrecs->ref_hash, ref);
    return k == kh_end(hrecs->ref_hash) ? -1 : kh_val(hrecs->ref_hash, k);
}

const char *sam_hdr_tid2name(const sam_hdr_t *h, int tid) {
    sam_hrecs_t *hrecs;

    if (!h || tid < 0)
        return NULL;

    if ((hrecs = h->hrecs) != NULL && tid < hrecs->nref) {
        return hrecs->ref[tid].name;
    } else {
        if (tid < h->n_targets)
            return h->target_name[tid];
    }

    return NULL;
}

hts_pos_t sam_hdr_tid2len(const sam_hdr_t *h, int tid) {
    sam_hrecs_t *hrecs;

    if (!h || tid < 0)
        return 0;

    if ((hrecs = h->hrecs) != NULL && tid < hrecs->nref) {
        return hrecs->ref[tid].len;
    } else {
        if (tid < h->n_targets) {
            if (h->target_len[tid] < UINT32_MAX || !h->sdict) {
                return h->target_len[tid];
            } else {
                khash_t(s2i) *long_refs = (khash_t(s2i) *) h->sdict;
                khint_t k = kh_get(s2i, long_refs, h->target_name[tid]);
                if (k < kh_end(long_refs)) {
                    return kh_val(long_refs, k);
                } else {
                    return UINT32_MAX;
                }
            }
        }
    }

    return 0;
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
static int sam_hdr_link_pg(sam_hdr_t *bh) {
    sam_hrecs_t *hrecs;
    int i, j, ret = 0, *new_pg_end;

    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    if (!hrecs->pgs_changed)
        return 0;

    hrecs->npg_end_alloc = hrecs->npg;
    new_pg_end = realloc(hrecs->pg_end, hrecs->npg * sizeof(*new_pg_end));
    if (!new_pg_end)
        return -1;
    hrecs->pg_end = new_pg_end;

    for (i = 0; i < hrecs->npg; i++)
        hrecs->pg_end[i] = i;

    for (i = 0; i < hrecs->npg; i++) {
        khint_t k;
        sam_hrec_tag_t *tag;

        assert(hrecs->pg[i].ty != NULL);
        for (tag = hrecs->pg[i].ty->tag; tag; tag = tag->next) {
            if (tag->str[0] == 'P' && tag->str[1] == 'P')
                break;
        }
        if (!tag) {
            /* Chain start points */
            continue;
        }

        k = kh_get(m_s2i, hrecs->pg_hash, tag->str+3);

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
    hrecs->pgs_changed = 0;

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
const char *sam_hdr_pg_id(sam_hdr_t *bh, const char *name) {
    sam_hrecs_t *hrecs;
    size_t name_len;
    const size_t name_extra = 17;
    if (!bh || !name)
        return NULL;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return NULL;
        hrecs = bh->hrecs;
    }

    khint_t k = kh_get(m_s2i, hrecs->pg_hash, name);
    if (k == kh_end(hrecs->pg_hash))
        return name;

    name_len = strlen(name);
    if (name_len > 1000) name_len = 1000;
    if (hrecs->ID_buf_sz < name_len + name_extra) {
        char *new_ID_buf = realloc(hrecs->ID_buf, name_len + name_extra);
        if (new_ID_buf == NULL)
            return NULL;
        hrecs->ID_buf = new_ID_buf;
        hrecs->ID_buf_sz = name_len + name_extra;
    }

    do {
        snprintf(hrecs->ID_buf, hrecs->ID_buf_sz, "%.1000s.%d", name, hrecs->ID_cnt++);
        k = kh_get(m_s2i, hrecs->pg_hash, hrecs->ID_buf);
    } while (k != kh_end(hrecs->pg_hash));

    return hrecs->ID_buf;
}

/*
 * Add an @PG line.
 *
 * If we wish complete control over this use sam_hdr_add_line() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple @PG records if we have multiple PG chains.
 *
 * Call it as per sam_hdr_add_line() with a series of key,value pairs ending
 * in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_add_pg(sam_hdr_t *bh, const char *name, ...) {
    sam_hrecs_t *hrecs;
    const char *specified_id = NULL, *specified_pn = NULL, *specified_pp = NULL;
    const char *key, *val;
    if (!bh)
        return -1;

    if (!(hrecs = bh->hrecs)) {
        if (sam_hdr_fill_hrecs(bh) != 0)
            return -1;
        hrecs = bh->hrecs;
    }

    va_list args;
    // Check for ID / PN / PP tags in varargs list
    va_start(args, name);
    while ((key = va_arg(args, const char *)) != NULL) {
        val = va_arg(args, const char *);
        if (!val) break;
        if (strcmp(key, "PN") == 0 && *val != '\0')
            specified_pn = val;
        else if (strcmp(key, "PP") == 0 && *val != '\0')
            specified_pp = val;
        else if (strcmp(key, "ID") == 0 && *val != '\0')
            specified_id = val;
    }
    va_end(args);

    if (specified_id && hrecs->pg_hash) {
        khint_t k = kh_get(m_s2i, hrecs->pg_hash, specified_id);
        if (k != kh_end(hrecs->pg_hash)) {
            hts_log_error("Header @PG ID:%s already present", specified_id);
            return -1;
        }
    }

    if (specified_pp && hrecs->pg_hash) {
        khint_t k = kh_get(m_s2i, hrecs->pg_hash, specified_pp);
        if (k == kh_end(hrecs->pg_hash)) {
            hts_log_error("Header @PG ID:%s referred to by PP tag not present",
                          specified_pp);
            return -1;
        }
    }

    if (!specified_pp && hrecs->npg_end) {
        /* Copy ends array to avoid us looping while modifying it */
        int *end = malloc(hrecs->npg_end * sizeof(int));
        int i, nends = hrecs->npg_end;

        if (!end)
            return -1;

        memcpy(end, hrecs->pg_end, nends * sizeof(*end));

        for (i = 0; i < nends; i++) {
            const char *id = !specified_id ? sam_hdr_pg_id(bh, name) : "";
            if (!id) {
                free(end);
                return -1;
            }
            va_start(args, name);
            if (-1 == sam_hrecs_vadd(hrecs, "PG", args,
                                     "ID", id,
                                     "PN", !specified_pn ? name : "",
                                     "PP", hrecs->pg[end[i]].name,
                                     NULL)) {
                free(end);
                return  -1;
            }
            va_end(args);
        }

        free(end);
    } else {
        const char *id = !specified_id ? sam_hdr_pg_id(bh, name) : "";
        if (!id)
            return -1;
        va_start(args, name);
        if (-1 == sam_hrecs_vadd(hrecs, "PG", args,
                                 "ID", id,
                                 "PN", !specified_pn ? name : "",
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
 * sam_hdr_destroy when done, without causing errors for other open files.
 */
void sam_hdr_incr_ref(sam_hdr_t *bh) {
    if (!bh)
        return;
    bh->ref_count++;
}

/* ==== Internal methods ==== */

/*
 * Creates an empty SAM header.  Allocates space for the SAM header
 * structures (hash tables) ready to be populated.
 *
 * Returns a sam_hrecs_t struct on success (free with sam_hrecs_free())
 *         NULL on failure
 */
sam_hrecs_t *sam_hrecs_new() {
    sam_hrecs_t *hrecs = calloc(1, sizeof(*hrecs));

    if (!hrecs)
        return NULL;

    hrecs->h = kh_init(sam_hrecs_t);
    if (!hrecs->h)
        goto err;

    hrecs->ID_cnt = 1;

    hrecs->nref = 0;
    hrecs->ref_sz = 0;
    hrecs->ref  = NULL;
    if (!(hrecs->ref_hash = kh_init(m_s2i)))
        goto err;
    hrecs->refs_changed = -1;

    hrecs->nrg = 0;
    hrecs->rg_sz = 0;
    hrecs->rg  = NULL;
    if (!(hrecs->rg_hash = kh_init(m_s2i)))
        goto err;

    hrecs->npg = 0;
    hrecs->pg_sz = 0;
    hrecs->pg  = NULL;
    hrecs->npg_end = hrecs->npg_end_alloc = 0;
    hrecs->pg_end = NULL;
    if (!(hrecs->pg_hash = kh_init(m_s2i)))
        goto err;

    if (!(hrecs->tag_pool = pool_create(sizeof(sam_hrec_tag_t))))
        goto err;

    if (!(hrecs->type_pool = pool_create(sizeof(sam_hrec_type_t))))
        goto err;

    if (!(hrecs->str_pool = string_pool_create(65536)))
        goto err;

    if (sam_hrecs_init_type_order(hrecs, NULL))
        goto err;

    return hrecs;

err:
    if (hrecs->h)
        kh_destroy(sam_hrecs_t, hrecs->h);

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
sam_hrecs_t *sam_hrecs_dup(sam_hrecs_t *source) {
        return NULL;
}
#endif
/*! Deallocates all storage used by a sam_hrecs_t struct.
 *
 * This also decrements the header reference count. If after decrementing
 * it is still non-zero then the header is assumed to be in use by another
 * caller and the free is not done.
 *
 */
void sam_hrecs_free(sam_hrecs_t *hrecs) {
    if (!hrecs)
        return;

    if (hrecs->h)
        kh_destroy(sam_hrecs_t, hrecs->h);

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

    if (hrecs->ID_buf)
        free(hrecs->ID_buf);

    free(hrecs);
}

/*
 * Internal method already used by the CRAM code
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
sam_hrec_type_t *sam_hrecs_find_type_id(sam_hrecs_t *hrecs, const char *type,
                                     const char *ID_key, const char *ID_value) {
    if (!hrecs || !type)
        return NULL;
    sam_hrec_type_t *t1, *t2;
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

    k = kh_get(sam_hrecs_t, hrecs->h, itype);
    if (k == kh_end(hrecs->h))
        return NULL;

    if (!ID_key)
        return kh_val(hrecs->h, k);

    t1 = t2 = kh_val(hrecs->h, k);
    do {
        sam_hrec_tag_t *tag;
        for (tag = t1->tag; tag; tag = tag->next) {
            if (tag->str[0] == ID_key[0] && tag->str[1] == ID_key[1]) {
                const char *cp1 = tag->str+3;
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
 * va_list contains multiple key,value pairs ending in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hrecs_vupdate(sam_hrecs_t *hrecs, sam_hrec_type_t *type, va_list ap) {
    if (!hrecs)
        return -1;

    for (;;) {
        char *k, *v, *str;
        sam_hrec_tag_t *tag, *prev = NULL;

        if (!(k = (char *)va_arg(ap, char *)))
            break;
        if (!(v = va_arg(ap, char *)))
            v = "";

        tag = sam_hrecs_find_key(type, k, &prev);
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
        str = string_alloc(hrecs->str_pool, tag->len+1);
        if (!str)
            return -1;

        if (snprintf(str, tag->len+1, "%2.2s:%s", k, v) < 0)
            return -1;

        tag->str = str;
    }

    hrecs->dirty = 1; //mark text as dirty and force a rebuild

    return 0;
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
static int sam_hrecs_update(sam_hrecs_t *hrecs, sam_hrec_type_t *type, ...) {
    va_list args;
    int res;
    va_start(args, type);
    res = sam_hrecs_vupdate(hrecs, type, args);
    va_end(args);
    return res;
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
sam_hrec_tag_t *sam_hrecs_find_key(sam_hrec_type_t *type,
                                   const char *key,
                                   sam_hrec_tag_t **prev) {
    sam_hrec_tag_t *tag, *p = NULL;
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

int sam_hrecs_remove_key(sam_hrecs_t *hrecs,
                         sam_hrec_type_t *type,
                         const char *key) {
    sam_hrec_tag_t *tag, *prev;
    if (!hrecs)
        return -1;
    tag = sam_hrecs_find_key(type, key, &prev);
    if (!tag)
        return 0; // Not there anyway

    if (type->type == TYPEKEY("SQ") && tag->str[0] == 'A' && tag->str[1] == 'N') {
        assert(tag->len >= 3);
        sam_hrec_tag_t *sn_tag = sam_hrecs_find_key(type, "SN", NULL);
        if (sn_tag) {
            assert(sn_tag->len >= 3);
            khint_t k = kh_get(m_s2i, hrecs->ref_hash, sn_tag->str + 3);
            if (k != kh_end(hrecs->ref_hash))
                sam_hrecs_remove_ref_altnames(hrecs, kh_val(hrecs->ref_hash, k), tag->str + 3);
        }
    }

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
sam_hrec_rg_t *sam_hrecs_find_rg(sam_hrecs_t *hrecs, const char *rg) {
    khint_t k = kh_get(m_s2i, hrecs->rg_hash, rg);
    return k == kh_end(hrecs->rg_hash)
        ? NULL
        : &hrecs->rg[kh_val(hrecs->rg_hash, k)];
}

#if DEBUG_HEADER
void sam_hrecs_dump(sam_hrecs_t *hrecs) {
    khint_t k;
    int i;

    printf("===DUMP===\n");
    for (k = kh_begin(hrecs->h); k != kh_end(hrecs->h); k++) {
        sam_hrec_type_t *t1, *t2;
        char c[2];
        int idx = 0;

        if (!kh_exist(hrecs->h, k))
            continue;

        t1 = t2 = kh_val(hrecs->h, k);
        c[0] = kh_key(hrecs->h, k)>>8;
        c[1] = kh_key(hrecs->h, k)&0xff;
        printf("Type %.2s\n", c);

        do {
            sam_hrec_tag_t *tag;
            printf(">>>%d ", idx++);
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
#endif

/*
 * Returns the sort order:
 */
enum sam_sort_order sam_hrecs_sort_order(sam_hrecs_t *hrecs) {
    khint_t k;
    enum sam_sort_order so;

    so = ORDER_UNKNOWN;
    k = kh_get(sam_hrecs_t, hrecs->h, TYPEKEY("HD"));
    if (k != kh_end(hrecs->h)) {
        sam_hrec_type_t *ty = kh_val(hrecs->h, k);
        sam_hrec_tag_t *tag;
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

enum sam_group_order sam_hrecs_group_order(sam_hrecs_t *hrecs) {
    khint_t k;
    enum sam_group_order go;

    go = ORDER_NONE;
    k = kh_get(sam_hrecs_t, hrecs->h, TYPEKEY("HD"));
    if (k != kh_end(hrecs->h)) {
        sam_hrec_type_t *ty = kh_val(hrecs->h, k);
        sam_hrec_tag_t *tag;
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
