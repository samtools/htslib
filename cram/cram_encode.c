/*
Copyright (c) 2012-2020 Genome Research Ltd.
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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <inttypes.h>

#include "cram.h"
#include "os.h"
#include "../sam_internal.h" // for nibble2base
#include "../htslib/hts.h"
#include "../htslib/hts_endian.h"

KHASH_MAP_INIT_STR(m_s2u64, uint64_t)

#define Z_CRAM_STRAT Z_FILTERED
//#define Z_CRAM_STRAT Z_RLE
//#define Z_CRAM_STRAT Z_HUFFMAN_ONLY
//#define Z_CRAM_STRAT Z_DEFAULT_STRATEGY

static int process_one_read(cram_fd *fd, cram_container *c,
                            cram_slice *s, cram_record *cr,
                            bam_seq_t *b, int rnum, kstring_t *MD);

/*
 * Returns index of val into key.
 * Basically strchr(key, val)-key;
 */
static int sub_idx(char *key, char val) {
    int i;

    for (i = 0; *key && *key++ != val; i++);
    return i;
}

/*
 * Encodes a compression header block into a generic cram_block structure.
 *
 * Returns cram_block ptr on success
 *         NULL on failure
 */
cram_block *cram_encode_compression_header(cram_fd *fd, cram_container *c,
                                           cram_block_compression_hdr *h) {
    cram_block *cb  = cram_new_block(COMPRESSION_HEADER, 0);
    cram_block *map = cram_new_block(COMPRESSION_HEADER, 0);
    int i, mc, r = 0;

    if (!cb || !map)
        return NULL;

    /*
     * This is a concatenation of several blocks of data:
     * header + landmarks, preservation map, read encoding map, and the tag
     * encoding map.
     * All 4 are variable sized and we need to know how large these are
     * before creating the compression header itself as this starts with
     * the total size (stored as a variable length string).
     */

/*
 * LARGE_POS used in this code is purely a debugging mechanism for testing
 * whether the htslib API can cope with 64-bit quantities.  These are
 * possible in SAM, but not *yet* in BAM or CRAM.
 *
 * DO NOT ENABLE LARGE_POS for anything other than debugging / testing.
 *
 * At some point it is expected these ifdefs will become a version check
 * instead.
 */

    // Duplicated from container itself, and removed in 1.1
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        r |= itf8_put_blk(cb, h->ref_seq_id);
#ifdef LARGE_POS
        r |= ltf8_put_blk(cb, h->ref_seq_start);
        r |= ltf8_put_blk(cb, h->ref_seq_span);
#else
        r |= itf8_put_blk(cb, h->ref_seq_start);
        r |= itf8_put_blk(cb, h->ref_seq_span);
#endif
        r |= itf8_put_blk(cb, h->num_records);
        r |= itf8_put_blk(cb, h->num_landmarks);
        for (i = 0; i < h->num_landmarks; i++) {
            r |= itf8_put_blk(cb, h->landmark[i]);
        }
    }

    if (h->preservation_map)
        kh_destroy(map, h->preservation_map);

    /* Create in-memory preservation map */
    /* FIXME: should create this when we create the container */
    {
        khint_t k;
        int r;

        if (!(h->preservation_map = kh_init(map)))
            return NULL;

        k = kh_put(map, h->preservation_map, "RN", &r);
        if (-1 == r) return NULL;
        kh_val(h->preservation_map, k).i = !fd->lossy_read_names;

        if (CRAM_MAJOR_VERS(fd->version) == 1) {
            k = kh_put(map, h->preservation_map, "PI", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = 0;

            k = kh_put(map, h->preservation_map, "UI", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = 1;

            k = kh_put(map, h->preservation_map, "MI", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = 1;

        } else {
            // Technically SM was in 1.0, but wasn't in Java impl.
            k = kh_put(map, h->preservation_map, "SM", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = 0;

            k = kh_put(map, h->preservation_map, "TD", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = 0;

            k = kh_put(map, h->preservation_map, "AP", &r);
            if (-1 == r) return NULL;
            kh_val(h->preservation_map, k).i = h->AP_delta;

            if (fd->no_ref || fd->embed_ref) {
                // Reference Required == No
                k = kh_put(map, h->preservation_map, "RR", &r);
                if (-1 == r) return NULL;
                kh_val(h->preservation_map, k).i = 0;
            }
        }
    }

    /* Encode preservation map; could collapse this and above into one */
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (h->preservation_map) {
        khint_t k;

        for (k = kh_begin(h->preservation_map);
             k != kh_end(h->preservation_map);
             k++) {
            const char *key;
            khash_t(map) *pmap = h->preservation_map;


            if (!kh_exist(pmap, k))
                continue;

            key = kh_key(pmap, k);
            BLOCK_APPEND(map, key, 2);

            switch(CRAM_KEY(key[0], key[1])) {
            case CRAM_KEY('M','I'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('U','I'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('P','I'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('A','P'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('R','N'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('R','R'):
                BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
                break;

            case CRAM_KEY('S','M'): {
                char smat[5], *mp = smat;
                *mp++ =
                    (sub_idx("CGTN", h->substitution_matrix[0][0]) << 6) |
                    (sub_idx("CGTN", h->substitution_matrix[0][1]) << 4) |
                    (sub_idx("CGTN", h->substitution_matrix[0][2]) << 2) |
                    (sub_idx("CGTN", h->substitution_matrix[0][3]) << 0);
                *mp++ =
                    (sub_idx("AGTN", h->substitution_matrix[1][0]) << 6) |
                    (sub_idx("AGTN", h->substitution_matrix[1][1]) << 4) |
                    (sub_idx("AGTN", h->substitution_matrix[1][2]) << 2) |
                    (sub_idx("AGTN", h->substitution_matrix[1][3]) << 0);
                *mp++ =
                    (sub_idx("ACTN", h->substitution_matrix[2][0]) << 6) |
                    (sub_idx("ACTN", h->substitution_matrix[2][1]) << 4) |
                    (sub_idx("ACTN", h->substitution_matrix[2][2]) << 2) |
                    (sub_idx("ACTN", h->substitution_matrix[2][3]) << 0);
                *mp++ =
                    (sub_idx("ACGN", h->substitution_matrix[3][0]) << 6) |
                    (sub_idx("ACGN", h->substitution_matrix[3][1]) << 4) |
                    (sub_idx("ACGN", h->substitution_matrix[3][2]) << 2) |
                    (sub_idx("ACGN", h->substitution_matrix[3][3]) << 0);
                *mp++ =
                    (sub_idx("ACGT", h->substitution_matrix[4][0]) << 6) |
                    (sub_idx("ACGT", h->substitution_matrix[4][1]) << 4) |
                    (sub_idx("ACGT", h->substitution_matrix[4][2]) << 2) |
                    (sub_idx("ACGT", h->substitution_matrix[4][3]) << 0);
                BLOCK_APPEND(map, smat, 5);
                break;
            }

            case CRAM_KEY('T','D'): {
                r |= itf8_put_blk(map, BLOCK_SIZE(h->TD_blk));
                BLOCK_APPEND(map,
                             BLOCK_DATA(h->TD_blk),
                             BLOCK_SIZE(h->TD_blk));
                break;
            }

            default:
                hts_log_warning("Unknown preservation key '%.2s'", key);
                break;
            }

            mc++;
        }
    }
    r |= itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    r |= itf8_put_blk(cb, mc);
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));

    /* rec encoding map */
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (h->codecs[DS_BF]) {
        if (-1 == h->codecs[DS_BF]->store(h->codecs[DS_BF], map, "BF",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_CF]) {
        if (-1 == h->codecs[DS_CF]->store(h->codecs[DS_CF], map, "CF",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_RL]) {
        if (-1 == h->codecs[DS_RL]->store(h->codecs[DS_RL], map, "RL",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_AP]) {
        if (-1 == h->codecs[DS_AP]->store(h->codecs[DS_AP], map, "AP",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_RG]) {
        if (-1 == h->codecs[DS_RG]->store(h->codecs[DS_RG], map, "RG",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_MF]) {
        if (-1 == h->codecs[DS_MF]->store(h->codecs[DS_MF], map, "MF",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_NS]) {
        if (-1 == h->codecs[DS_NS]->store(h->codecs[DS_NS], map, "NS",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_NP]) {
        if (-1 == h->codecs[DS_NP]->store(h->codecs[DS_NP], map, "NP",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_TS]) {
        if (-1 == h->codecs[DS_TS]->store(h->codecs[DS_TS], map, "TS",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_NF]) {
        if (-1 == h->codecs[DS_NF]->store(h->codecs[DS_NF], map, "NF",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_TC]) {
        if (-1 == h->codecs[DS_TC]->store(h->codecs[DS_TC], map, "TC",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_TN]) {
        if (-1 == h->codecs[DS_TN]->store(h->codecs[DS_TN], map, "TN",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_TL]) {
        if (-1 == h->codecs[DS_TL]->store(h->codecs[DS_TL], map, "TL",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_FN]) {
        if (-1 == h->codecs[DS_FN]->store(h->codecs[DS_FN], map, "FN",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_FC]) {
        if (-1 == h->codecs[DS_FC]->store(h->codecs[DS_FC], map, "FC",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_FP]) {
        if (-1 == h->codecs[DS_FP]->store(h->codecs[DS_FP], map, "FP",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_BS]) {
        if (-1 == h->codecs[DS_BS]->store(h->codecs[DS_BS], map, "BS",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_IN]) {
        if (-1 == h->codecs[DS_IN]->store(h->codecs[DS_IN], map, "IN",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_DL]) {
        if (-1 == h->codecs[DS_DL]->store(h->codecs[DS_DL], map, "DL",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_BA]) {
        if (-1 == h->codecs[DS_BA]->store(h->codecs[DS_BA], map, "BA",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_BB]) {
        if (-1 == h->codecs[DS_BB]->store(h->codecs[DS_BB], map, "BB",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_MQ]) {
        if (-1 == h->codecs[DS_MQ]->store(h->codecs[DS_MQ], map, "MQ",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_RN]) {
        if (-1 == h->codecs[DS_RN]->store(h->codecs[DS_RN], map, "RN",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_QS]) {
        if (-1 == h->codecs[DS_QS]->store(h->codecs[DS_QS], map, "QS",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_QQ]) {
        if (-1 == h->codecs[DS_QQ]->store(h->codecs[DS_QQ], map, "QQ",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_RI]) {
        if (-1 == h->codecs[DS_RI]->store(h->codecs[DS_RI], map, "RI",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (CRAM_MAJOR_VERS(fd->version) != 1) {
        if (h->codecs[DS_SC]) {
            if (-1 == h->codecs[DS_SC]->store(h->codecs[DS_SC], map, "SC",
                                              fd->version))
                return NULL;
            mc++;
        }
        if (h->codecs[DS_RS]) {
            if (-1 == h->codecs[DS_RS]->store(h->codecs[DS_RS], map, "RS",
                                              fd->version))
                return NULL;
            mc++;
        }
        if (h->codecs[DS_PD]) {
            if (-1 == h->codecs[DS_PD]->store(h->codecs[DS_PD], map, "PD",
                                              fd->version))
                return NULL;
            mc++;
        }
        if (h->codecs[DS_HC]) {
            if (-1 == h->codecs[DS_HC]->store(h->codecs[DS_HC], map, "HC",
                                              fd->version))
                return NULL;
            mc++;
        }
    }
    if (h->codecs[DS_TM]) {
        if (-1 == h->codecs[DS_TM]->store(h->codecs[DS_TM], map, "TM",
                                          fd->version))
            return NULL;
        mc++;
    }
    if (h->codecs[DS_TV]) {
        if (-1 == h->codecs[DS_TV]->store(h->codecs[DS_TV], map, "TV",
                                          fd->version))
            return NULL;
        mc++;
    }
    r |= itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    r |= itf8_put_blk(cb, mc);
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));

    /* tag encoding map */
#if 0
    mp = map; mc = 0;
    if (h->tag_encoding_map) {
        HashItem *hi;
        HashIter *iter = HashTableIterCreate();
        if (!iter)
            return NULL;

        while ((hi = HashTableIterNext(h->tag_encoding_map, iter))) {
            cram_map *m = hi->data.p;
            int sz;

            mp += itf8_put(mp, (hi->key[0]<<16)|(hi->key[1]<<8)|hi->key[2]);
            if (-1 == (sz = m->codec->store(m->codec, mp, NULL, fd->version)))
                return NULL;
            mp += sz;
            mc++;
        }

        HashTableIterDestroy(iter);
    }
#else
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (c->tags_used) {
        khint_t k;

        for (k = kh_begin(c->tags_used); k != kh_end(c->tags_used); k++) {
            int key;
            if (!kh_exist(c->tags_used, k))
                continue;

            key = kh_key(c->tags_used, k);
            cram_codec *cd = kh_val(c->tags_used, k)->codec;

            r |= itf8_put_blk(map, key);
            if (-1 == cd->store(cd, map, NULL, fd->version))
                return NULL;

            mc++;
        }
    }
#endif
    r |= itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    r |= itf8_put_blk(cb, mc);
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));

    hts_log_info("Wrote compression block header in %d bytes", (int)BLOCK_SIZE(cb));

    BLOCK_UPLEN(cb);

    cram_free_block(map);

    if (r >= 0)
        return cb;

 block_err:
    return NULL;
}


/*
 * Encodes a slice compression header.
 *
 * Returns cram_block on success
 *         NULL on failure
 */
cram_block *cram_encode_slice_header(cram_fd *fd, cram_slice *s) {
    char *buf;
    char *cp;
    cram_block *b = cram_new_block(MAPPED_SLICE, 0);
    int j;

    if (!b)
        return NULL;

    cp = buf = malloc(16+5*(8+s->hdr->num_blocks));
    if (NULL == buf) {
        cram_free_block(b);
        return NULL;
    }

    cp += itf8_put(cp, s->hdr->ref_seq_id);
#ifdef LARGE_POS
    cp += ltf8_put(cp, s->hdr->ref_seq_start);
    cp += ltf8_put(cp, s->hdr->ref_seq_span);
#else
    cp += itf8_put(cp, s->hdr->ref_seq_start);
    cp += itf8_put(cp, s->hdr->ref_seq_span);
#endif
    cp += itf8_put(cp, s->hdr->num_records);
    if (CRAM_MAJOR_VERS(fd->version) == 2)
        cp += itf8_put(cp, s->hdr->record_counter);
    else if (CRAM_MAJOR_VERS(fd->version) >= 3)
        cp += ltf8_put(cp, s->hdr->record_counter);
    cp += itf8_put(cp, s->hdr->num_blocks);
    cp += itf8_put(cp, s->hdr->num_content_ids);
    for (j = 0; j < s->hdr->num_content_ids; j++) {
        cp += itf8_put(cp, s->hdr->block_content_ids[j]);
    }
    if (s->hdr->content_type == MAPPED_SLICE)
        cp += itf8_put(cp, s->hdr->ref_base_id);

    if (CRAM_MAJOR_VERS(fd->version) != 1) {
        memcpy(cp, s->hdr->md5, 16); cp += 16;
    }

    assert(cp-buf <= 16+5*(8+s->hdr->num_blocks));

    b->data = (unsigned char *)buf;
    b->comp_size = b->uncomp_size = cp-buf;

    return b;
}


/*
 * Encodes a single read.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_encode_slice_read(cram_fd *fd,
                                  cram_container *c,
                                  cram_block_compression_hdr *h,
                                  cram_slice *s,
                                  cram_record *cr,
                                  int64_t *last_pos) {
    int r = 0;
    int32_t i32;
    unsigned char uc;

    //fprintf(stderr, "Encode seq %d, %d/%d FN=%d, %s\n", rec, core->byte, core->bit, cr->nfeature, s->name_ds->str + cr->name);

    //printf("BF=0x%x\n", cr->flags);
    //      bf = cram_flag_swap[cr->flags];
    i32 = fd->cram_flag_swap[cr->flags & 0xfff];
    r |= h->codecs[DS_BF]->encode(s, h->codecs[DS_BF], (char *)&i32, 1);

    i32 = cr->cram_flags & CRAM_FLAG_MASK;
    r |= h->codecs[DS_CF]->encode(s, h->codecs[DS_CF], (char *)&i32, 1);

    if (CRAM_MAJOR_VERS(fd->version) != 1 && s->hdr->ref_seq_id == -2)
        r |= h->codecs[DS_RI]->encode(s, h->codecs[DS_RI], (char *)&cr->ref_id, 1);

    r |= h->codecs[DS_RL]->encode(s, h->codecs[DS_RL], (char *)&cr->len, 1);

    if (c->pos_sorted) {
#ifdef LARGE_POS
        int64_t i64;
        i64 = cr->apos - *last_pos;
        r |= h->codecs[DS_AP]->encode(s, h->codecs[DS_AP], (char *)&i64, 1);
#else
        i32 = cr->apos - *last_pos;
        r |= h->codecs[DS_AP]->encode(s, h->codecs[DS_AP], (char *)&i32, 1);
#endif
        *last_pos = cr->apos;
    } else {
#ifdef LARGE_POS
        int64_t i64;
        i64 = cr->apos;
        r |= h->codecs[DS_AP]->encode(s, h->codecs[DS_AP], (char *)&i64, 1);
#else
        i32 = cr->apos;
        r |= h->codecs[DS_AP]->encode(s, h->codecs[DS_AP], (char *)&i32, 1);
#endif
    }

    r |= h->codecs[DS_RG]->encode(s, h->codecs[DS_RG], (char *)&cr->rg, 1);

    if (cr->cram_flags & CRAM_FLAG_DETACHED) {
        i32 = cr->mate_flags;
        r |= h->codecs[DS_MF]->encode(s, h->codecs[DS_MF], (char *)&i32, 1);

        r |= h->codecs[DS_NS]->encode(s, h->codecs[DS_NS],
                                      (char *)&cr->mate_ref_id, 1);

#ifdef LARGE_POS
        r |= h->codecs[DS_NP]->encode(s, h->codecs[DS_NP],
                                      (char *)&cr->mate_pos, 1);

        r |= h->codecs[DS_TS]->encode(s, h->codecs[DS_TS],
                                      (char *)&cr->tlen, 1);
#else
        i32 = cr->mate_pos;
        r |= h->codecs[DS_NP]->encode(s, h->codecs[DS_NP],
                                      (char *)&i32, 1);
        i32 = cr->tlen;
        r |= h->codecs[DS_TS]->encode(s, h->codecs[DS_TS],
                                      (char *)&i32, 1);
#endif
    } else if (cr->cram_flags & CRAM_FLAG_MATE_DOWNSTREAM) {
        r |= h->codecs[DS_NF]->encode(s, h->codecs[DS_NF],
                                      (char *)&cr->mate_line, 1);
    }

    /* Aux tags */
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        int j;
        uc = cr->ntags;
        r |= h->codecs[DS_TC]->encode(s, h->codecs[DS_TC], (char *)&uc, 1);

        for (j = 0; j < cr->ntags; j++) {
            uint32_t i32 = s->TN[cr->TN_idx + j]; // id
            r |= h->codecs[DS_TN]->encode(s, h->codecs[DS_TN], (char *)&i32, 1);
        }
    } else {
        r |= h->codecs[DS_TL]->encode(s, h->codecs[DS_TL], (char *)&cr->TL, 1);
    }

    // qual
    // QS codec : Already stored in block[2].

    // features (diffs)
    if (!(cr->flags & BAM_FUNMAP)) {
        int prev_pos = 0, j;

        r |= h->codecs[DS_FN]->encode(s, h->codecs[DS_FN],
                                      (char *)&cr->nfeature, 1);
        for (j = 0; j < cr->nfeature; j++) {
            cram_feature *f = &s->features[cr->feature + j];

            uc = f->X.code;
            r |= h->codecs[DS_FC]->encode(s, h->codecs[DS_FC], (char *)&uc, 1);
            i32 = f->X.pos - prev_pos;
            r |= h->codecs[DS_FP]->encode(s, h->codecs[DS_FP], (char *)&i32, 1);
            prev_pos = f->X.pos;

            switch(f->X.code) {
                //char *seq;

            case 'X':
                //fprintf(stderr, "    FC=%c FP=%d base=%d\n", f->X.code, i32, f->X.base);

                uc = f->X.base;
                r |= h->codecs[DS_BS]->encode(s, h->codecs[DS_BS],
                                              (char *)&uc, 1);
                break;
            case 'S':
                // Already done
                //r |= h->codecs[DS_SC]->encode(s, h->codecs[DS_SC],
                //                              BLOCK_DATA(s->soft_blk) + f->S.seq_idx,
                //                              f->S.len);

                //if (IS_CRAM_3_VERS(fd)) {
                //    r |= h->codecs[DS_BB]->encode(s, h->codecs[DS_BB],
                //                                  BLOCK_DATA(s->seqs_blk) + f->S.seq_idx,
                //                                  f->S.len);
                //}
                break;
            case 'I':
                //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
                //r |= h->codecs[DS_IN]->encode(s, h->codecs[DS_IN],
                //                           seq, f->S.len);
                //if (IS_CRAM_3_VERS(fd)) {
                //    r |= h->codecs[DS_BB]->encode(s, h->codecs[DS_BB],
                //                                  BLOCK_DATA(s->seqs_blk) + f->I.seq_idx,
                //                                  f->I.len);
                //}
                break;
            case 'i':
                uc = f->i.base;
                r |= h->codecs[DS_BA]->encode(s, h->codecs[DS_BA],
                                              (char *)&uc, 1);
                //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
                //r |= h->codecs[DS_IN]->encode(s, h->codecs[DS_IN],
                //                           seq, 1);
                break;
            case 'D':
                i32 = f->D.len;
                r |= h->codecs[DS_DL]->encode(s, h->codecs[DS_DL],
                                              (char *)&i32, 1);
                break;

            case 'B':
                //                  // Used when we try to store a non ACGTN base or an N
                //                  // that aligns against a non ACGTN reference

                uc  = f->B.base;
                r |= h->codecs[DS_BA]->encode(s, h->codecs[DS_BA],
                                              (char *)&uc, 1);

                //                  Already added
                //                  uc  = f->B.qual;
                //                  r |= h->codecs[DS_QS]->encode(s, h->codecs[DS_QS],
                //                                           (char *)&uc, 1);
                break;

            case 'b':
                // string of bases
                r |= h->codecs[DS_BB]->encode(s, h->codecs[DS_BB],
                                              (char *)BLOCK_DATA(s->seqs_blk)
                                                      + f->b.seq_idx,
                                              f->b.len);
                break;

            case 'Q':
                //                  Already added
                //                  uc  = f->B.qual;
                //                  r |= h->codecs[DS_QS]->encode(s, h->codecs[DS_QS],
                //                                           (char *)&uc, 1);
                break;

            case 'N':
                i32 = f->N.len;
                r |= h->codecs[DS_RS]->encode(s, h->codecs[DS_RS],
                                              (char *)&i32, 1);
                break;

            case 'P':
                i32 = f->P.len;
                r |= h->codecs[DS_PD]->encode(s, h->codecs[DS_PD],
                                              (char *)&i32, 1);
                break;

            case 'H':
                i32 = f->H.len;
                r |= h->codecs[DS_HC]->encode(s, h->codecs[DS_HC],
                                              (char *)&i32, 1);
                break;


            default:
                hts_log_error("Unhandled feature code %c", f->X.code);
                return -1;
            }
        }

        r |= h->codecs[DS_MQ]->encode(s, h->codecs[DS_MQ],
                                      (char *)&cr->mqual, 1);
    } else {
        char *seq = (char *)BLOCK_DATA(s->seqs_blk) + cr->seq;
        if (cr->len)
            r |= h->codecs[DS_BA]->encode(s, h->codecs[DS_BA], seq, cr->len);
    }

    return r ? -1 : 0;
}


/*
 * Applies various compression methods to specific blocks, depending on
 * known observations of how data series compress.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_compress_slice(cram_fd *fd, cram_container *c, cram_slice *s) {
    int level = fd->level, i;
    int method = 1<<GZIP | 1<<GZIP_RLE, methodF = method;

    /* Compress the CORE Block too, with minimal zlib level */
    if (level > 5 && s->block[0]->uncomp_size > 500)
        cram_compress_block(fd, s->block[0], NULL, 1<<GZIP, 1);

    if (fd->use_bz2)
        method |= 1<<BZIP2;

    if (fd->use_rans)
        method |= (1<<RANS0) | (1<<RANS1);

    if (fd->use_lzma)
        method |= (1<<LZMA);

    /* Faster method for data series we only need entropy encoding on */
    methodF = method & ~(1<<GZIP | 1<<BZIP2 | 1<<LZMA);
    if (level >= 6)
        methodF = method;


    /* Specific compression methods for certain block types */
    if (cram_compress_block(fd, s->block[DS_IN], fd->m[DS_IN], //IN (seq)
                            method, level))
        return -1;

    if (fd->level == 0) {
        /* Do nothing */
    } else if (fd->level == 1) {
        if (cram_compress_block(fd, s->block[DS_QS], fd->m[DS_QS],
                                methodF, 1))
            return -1;
        for (i = DS_aux; i <= DS_aux_oz; i++) {
            if (s->block[i])
                if (cram_compress_block(fd, s->block[i], fd->m[i],
                                        method, 1))
                    return -1;
        }
    } else if (fd->level < 3) {
        if (cram_compress_block(fd, s->block[DS_QS], fd->m[DS_QS],
                                method, 1))
            return -1;
        if (cram_compress_block(fd, s->block[DS_BA], fd->m[DS_BA],
                                method, 1))
            return -1;
        if (s->block[DS_BB])
            if (cram_compress_block(fd, s->block[DS_BB], fd->m[DS_BB],
                                    method, 1))
                return -1;
        for (i = DS_aux; i <= DS_aux_oz; i++) {
            if (s->block[i])
                if (cram_compress_block(fd, s->block[i], fd->m[i],
                                        method, level))
                    return -1;
        }
    } else {
        if (cram_compress_block(fd, s->block[DS_QS], fd->m[DS_QS],
                                method, level))
            return -1;
        if (cram_compress_block(fd, s->block[DS_BA], fd->m[DS_BA],
                                method, level))
            return -1;
        if (s->block[DS_BB])
            if (cram_compress_block(fd, s->block[DS_BB], fd->m[DS_BB],
                                    method, level))
                return -1;
        for (i = DS_aux; i <= DS_aux_oz; i++) {
            if (s->block[i])
                if (cram_compress_block(fd, s->block[i], fd->m[i],
                                        method, level))
                    return -1;
        }
    }

    // NAME: best is generally xz, bzip2, zlib then rans1
    if (cram_compress_block(fd, s->block[DS_RN], fd->m[DS_RN],
                            method & ~(1<<RANS0 | 1<<GZIP_RLE),
                            level))
        return -1;

    // NS shows strong local correlation as rearrangements are localised
    if (s->block[DS_NS] != s->block[0])
        if (cram_compress_block(fd, s->block[DS_NS], fd->m[DS_NS],
                                method, level))
            return -1;


    /*
     * Compress any auxiliary tags with their own per-tag metrics
     */
    {
        int i;
        for (i = 0; i < s->naux_block; i++) {
            if (!s->aux_block[i] || s->aux_block[i] == s->block[0])
                continue;

            if (s->aux_block[i]->method != RAW)
                continue;

            if (cram_compress_block(fd, s->aux_block[i], s->aux_block[i]->m,
                                    method, level))
                return -1;
        }
    }

    /*
     * Minimal compression of any block still uncompressed, bar CORE
     */
    {
        int i;
        for (i = 1; i < s->hdr->num_blocks && i < DS_END; i++) {
            if (!s->block[i] || s->block[i] == s->block[0])
                continue;

            if (s->block[i]->method != RAW)
                continue;

            if (cram_compress_block(fd, s->block[i], fd->m[i],
                                    methodF, level))
                return -1;
        }
    }

    return 0;
}

/*
 * Encodes a single slice from a container
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_encode_slice(cram_fd *fd, cram_container *c,
                             cram_block_compression_hdr *h, cram_slice *s) {
    int rec, r = 0;
    int64_t last_pos;
    int embed_ref;
    enum cram_DS_ID id;

    embed_ref = fd->embed_ref && s->hdr->ref_seq_id != -1 ? 1 : 0;

    /*
     * Slice external blocks:
     * ID 0 => base calls (insertions, soft-clip)
     * ID 1 => qualities
     * ID 2 => names
     * ID 3 => TS (insert size), NP (next frag)
     * ID 4 => tag values
     * ID 6 => tag IDs (TN), if CRAM_V1.0
     * ID 7 => TD tag dictionary, if !CRAM_V1.0
     */

    /* Create cram slice header */
    s->hdr->ref_base_id = embed_ref ? DS_ref : -1;
    s->hdr->record_counter = c->num_records + c->record_counter;
    c->num_records += s->hdr->num_records;

    int ntags = c->tags_used ? c->tags_used->n_occupied : 0;
    s->block = calloc(DS_END + ntags, sizeof(s->block[0]));
    s->hdr->block_content_ids = malloc(DS_END * sizeof(int32_t));
    if (!s->block || !s->hdr->block_content_ids)
        return -1;

    // Create first fixed blocks, always external.
    // CORE
    if (!(s->block[0] = cram_new_block(CORE, 0)))
        return -1;

    // TN block for CRAM v1
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        if (h->codecs[DS_TN]->codec == E_EXTERNAL) {
            if (!(s->block[DS_TN] = cram_new_block(EXTERNAL,DS_TN))) return -1;
            h->codecs[DS_TN]->u.external.content_id = DS_TN;
        } else {
            s->block[DS_TN] = s->block[0];
        }
    }

    // Embedded reference
    if (embed_ref) {
        if (!(s->block[DS_ref] = cram_new_block(EXTERNAL, DS_ref)))
            return -1;
        s->ref_id = DS_ref; // needed?
        BLOCK_APPEND(s->block[DS_ref],
                     c->ref + s->hdr->ref_seq_start - c->ref_start,
                     s->hdr->ref_seq_span);
    }

    /*
     * All the data-series blocks if appropriate.
     */
    for (id = DS_BF; id < DS_TN; id++) {
        if (h->codecs[id] && (h->codecs[id]->codec == E_EXTERNAL ||
                              h->codecs[id]->codec == E_BYTE_ARRAY_STOP ||
                              h->codecs[id]->codec == E_BYTE_ARRAY_LEN)) {
            switch (h->codecs[id]->codec) {
            case E_EXTERNAL:
                if (!(s->block[id] = cram_new_block(EXTERNAL, id)))
                    return -1;
                h->codecs[id]->u.external.content_id = id;
                break;

            case E_BYTE_ARRAY_STOP:
                if (!(s->block[id] = cram_new_block(EXTERNAL, id)))
                    return -1;
                h->codecs[id]->u.byte_array_stop.content_id = id;
                break;

            case E_BYTE_ARRAY_LEN: {
                cram_codec *cc;

                cc = h->codecs[id]->u.e_byte_array_len.len_codec;
                if (cc->codec == E_EXTERNAL) {
                    int eid = cc->u.external.content_id;
                    if (!(s->block[eid] = cram_new_block(EXTERNAL, eid)))
                        return -1;
                    cc->u.external.content_id = eid;
                    cc->out = s->block[eid];
                }

                cc = h->codecs[id]->u.e_byte_array_len.val_codec;
                if (cc->codec == E_EXTERNAL) {
                    int eid = cc->u.external.content_id;
                    if (!s->block[eid])
                        if (!(s->block[eid] = cram_new_block(EXTERNAL, eid)))
                            return -1;
                    cc->u.external.content_id = eid;
                    cc->out = s->block[eid];
                }
                break;
            }
            default:
                break;
            }
        } else {
            if (!(id == DS_BB && !h->codecs[DS_BB]))
                s->block[id] = s->block[0];
        }
        if (h->codecs[id])
            h->codecs[id]->out = s->block[id];
    }

    /*
     * Add in the external tag blocks too.
     */
    if (c->tags_used) {
        int n;
        s->hdr->num_blocks = DS_END;
        for (n = 0; n < s->naux_block; n++)
            s->block[s->hdr->num_blocks++] = s->aux_block[n];
    }

    /* Encode reads */
    last_pos = s->hdr->ref_seq_start;
    for (rec = 0; rec < s->hdr->num_records; rec++) {
        cram_record *cr = &s->crecs[rec];
        if (cram_encode_slice_read(fd, c, h, s, cr, &last_pos) == -1)
            return -1;
    }

    s->block[0]->uncomp_size = s->block[0]->byte + (s->block[0]->bit < 7);
    s->block[0]->comp_size = s->block[0]->uncomp_size;

    // Make sure the fixed blocks point to the correct sources
    s->block[DS_IN] = s->base_blk; s->base_blk = NULL;
    s->block[DS_QS] = s->qual_blk; s->qual_blk = NULL;
    s->block[DS_RN] = s->name_blk; s->name_blk = NULL;
    s->block[DS_SC] = s->soft_blk; s->soft_blk = NULL;

    // Ensure block sizes are up to date.
    for (id = 1; id < s->hdr->num_blocks; id++) {
        if (!s->block[id] || s->block[id] == s->block[0])
            continue;

        if (s->block[id]->uncomp_size == 0)
            BLOCK_UPLEN(s->block[id]);
    }

    // Compress it all
    if (cram_compress_slice(fd, c, s) == -1)
        return -1;

    // Collapse empty blocks and create hdr_block
    {
        int i, j;

        s->hdr->block_content_ids = realloc(s->hdr->block_content_ids,
                                            s->hdr->num_blocks * sizeof(int32_t));
        if (!s->hdr->block_content_ids)
            return -1;

        for (i = j = 1; i < s->hdr->num_blocks; i++) {
            if (!s->block[i] || s->block[i] == s->block[0])
                continue;
            if (s->block[i]->uncomp_size == 0) {
                cram_free_block(s->block[i]);
                s->block[i] = NULL;
                continue;
            }
            s->block[j] = s->block[i];
            s->hdr->block_content_ids[j-1] = s->block[i]->content_id;
            j++;
        }
        s->hdr->num_content_ids = j-1;
        s->hdr->num_blocks = j;

        if (!(s->hdr_block = cram_encode_slice_header(fd, s)))
            return -1;
    }

    return r ? -1 : 0;

 block_err:
    return -1;
}

/*
 * Returns the number of expected read names for this record.
 */
static int expected_template_count(bam_seq_t *b) {
    int expected = bam_flag(b) & BAM_FPAIRED ? 2 : 1;

    uint8_t *TC = (uint8_t *)bam_aux_get(b, "TC");
    if (TC) {
        int n = bam_aux2i(TC);
        if (expected < n)
            expected = n;
    }

    if (!TC && bam_aux_get(b, "SA")) {
        // We could count the semicolons, but we'd have to do this for
        // read1, read2 and read(not-1-or-2) combining the results
        // together.  This is a cheap and safe alternative for now.
        expected = INT_MAX;
    }

    return expected;
}

/*
 * Lossily reject read names.
 *
 * The rule here is that if all reads for this template reside in the
 * same slice then we can lose the name.  Otherwise we keep them as we
 * do not know when (or if) the other reads will turn up.
 *
 * Note there may be only 1 read (non-paired library) or more than 2
 * reads (paired library with supplementary reads), or other weird
 * setups.  We need to know how many are expected.  Ways to guess:
 *
 * - Flags (0x1 - has > 1 read)
 * - TC aux field (not mandatory)
 * - SA tags (count semicolons, NB per fragment so sum - hard)
 * - RNEXT/PNEXT uniqueness count. (not implemented, tricky)
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int lossy_read_names(cram_fd *fd, cram_container *c, cram_slice *s,
                            int bam_start) {
    int r1, r2, ret = -1;

    // Initialise cram_flags
    for (r2 = 0; r2 < s->hdr->num_records; r2++)
        s->crecs[r2].cram_flags = 0;

    if (!fd->lossy_read_names)
        return 0;

    khash_t(m_s2u64) *names = kh_init(m_s2u64);
    if (!names)
        goto fail;

    // 1: Iterate through names to count frequency
    for (r1 = bam_start, r2 = 0; r2 < s->hdr->num_records; r1++, r2++) {
        //cram_record *cr = &s->crecs[r2];
        bam_seq_t *b = c->bams[r1];
        khint_t k;
        int n;
        uint64_t e;
        union {
            uint64_t i64;
            struct {
                int32_t e,c; // expected & observed counts.
            } counts;
        } u;

        e = expected_template_count(b);
        u.counts.e = e; u.counts.c = 1;

        k = kh_put(m_s2u64, names, bam_name(b), &n);
        if (n == -1)
            goto fail;

        if (n == 0) {
            // not a new name
            u.i64 = kh_val(names, k);
            if (u.counts.e != e) {
                // different expectation or already hit the max
                //fprintf(stderr, "Err computing no. %s recs\n", bam_name(b));
                kh_val(names, k) = 0;
            } else {
                u.counts.c++;
                if (u.counts.e == u.counts.c) {
                    // Reached expected count.
                    kh_val(names, k) = -1;
                } else {
                    kh_val(names, k) = u.i64;
                }
            }
        } else {
            // new name
            kh_val(names, k) = u.i64;
        }
    }

    // 2: Remove names if all present (hd.i == -1)
    for (r1 = bam_start, r2 = 0; r2 < s->hdr->num_records; r1++, r2++) {
        cram_record *cr = &s->crecs[r2];
        bam_seq_t *b = c->bams[r1];
        khint_t k;

        k = kh_get(m_s2u64, names, bam_name(b));

        if (k == kh_end(names))
            goto fail;

        if (kh_val(names, k) == -1)
            cr->cram_flags = CRAM_FLAG_DISCARD_NAME;
    }

    ret = 0;
 fail: // ret==-1

    if (names)
        kh_destroy(m_s2u64, names);

    return ret;
}

/*
 * Adds the reading names.  We do this here as a separate pass rather
 * than per record in the process_one_read calls as that function can
 * go back and change the CRAM_FLAG_DETACHED status of a previously
 * processed read if it subsequently determines the TLEN field is
 * incorrect.  Given DETACHED reads always try to decode read names,
 * we need to know their status before generating the read-name block.
 *
 * Output is an update s->name_blk, and cr->name / cr->name_len
 * fields.
 */
static int add_read_names(cram_fd *fd, cram_container *c, cram_slice *s,
                          int bam_start) {
    int r1, r2;
    int keep_names = !fd->lossy_read_names;

    for (r1 = bam_start, r2 = 0;
         r1 < c->curr_c_rec && r2 < s->hdr->num_records;
         r1++, r2++) {
        cram_record *cr = &s->crecs[r2];
        bam_seq_t *b = c->bams[r1];

        cr->name        = BLOCK_SIZE(s->name_blk);
        if ((cr->cram_flags & CRAM_FLAG_DETACHED) || keep_names) {
            BLOCK_APPEND(s->name_blk, bam_name(b), bam_name_len(b));
            cr->name_len    = bam_name_len(b);
        } else {
            // Can only discard duplicate names if not detached
            cr->name_len = 0;
        }
        if (cram_stats_add(c->stats[DS_RN], cr->name_len) < 0)
            goto block_err;
    }

    return 0;

 block_err:
    return -1;
}

/*
 * Encodes all slices in a container into blocks.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_encode_container(cram_fd *fd, cram_container *c) {
    int i, j, slice_offset;
    cram_block_compression_hdr *h = c->comp_hdr;
    cram_block *c_hdr;
    int multi_ref = 0;
    int r1, r2, sn, nref;
    spare_bams *spares;

//#define goto_err {fprintf(stderr, "ERR at %s:%d\n", __FILE__, __LINE__);goto err;}
#define goto_err goto err

    /* Cache references up-front if we have unsorted access patterns */
    pthread_mutex_lock(&fd->ref_lock);
    nref = fd->refs->nref;
    pthread_mutex_unlock(&fd->ref_lock);

    if (!fd->no_ref && c->refs_used) {
        for (i = 0; i < nref; i++) {
            if (c->refs_used[i])
                cram_get_ref(fd, i, 1, 0);
        }
    }

    /* To create M5 strings */
    /* Fetch reference sequence */
    if (!fd->no_ref) {
        if (!c->bams || !c->bams[0])
            goto_err;
        bam_seq_t *b = c->bams[0];

        char *ref = cram_get_ref(fd, bam_ref(b), 1, 0);
        if (!ref && bam_ref(b) >= 0) {
            hts_log_error("Failed to load reference #%d", bam_ref(b));
            return -1;
        }
        if ((c->ref_id = bam_ref(b)) >= 0) {
            c->ref_seq_id = c->ref_id;
            c->ref       = fd->refs->ref_id[c->ref_seq_id]->seq;
            c->ref_start = 1;
            c->ref_end   = fd->refs->ref_id[c->ref_seq_id]->length;
        } else {
            c->ref_seq_id = c->ref_id; // FIXME remove one var!
        }
    } else {
        c->ref_id = bam_ref(c->bams[0]);
        cram_ref_incr(fd->refs, c->ref_id);
        c->ref_seq_id = c->ref_id;
    }

    /* Turn bams into cram_records and gather basic stats */
    for (r1 = sn = 0; r1 < c->curr_c_rec; sn++) {
        cram_slice *s = c->slices[sn];
        int64_t first_base = INT64_MAX, last_base = INT64_MIN;

        int r1_start = r1;

        assert(sn < c->curr_slice);

        // Discover which read names *may* be safely removed.
        // Ie which ones have all their records in this slice.
        if (lossy_read_names(fd, c, s, r1_start) != 0)
            return -1;

        // Tracking of MD tags so we can spot when the auto-generated values
        // will differ from the current stored ones.  The kstring here is
        // simply to avoid excessive malloc and free calls.  All initialisation
        // is done within process_one_read().
        kstring_t MD = {0};

        // Iterate through records creating the cram blocks for some
        // fields and just gathering stats for others.
        for (r2 = 0; r1 < c->curr_c_rec && r2 < s->hdr->num_records; r1++, r2++) {
            cram_record *cr = &s->crecs[r2];
            bam_seq_t *b = c->bams[r1];

            /* If multi-ref we need to cope with changing reference per seq */
            if (c->multi_seq && !fd->no_ref) {
                if (bam_ref(b) != c->ref_seq_id && bam_ref(b) >= 0) {
                    if (c->ref_seq_id >= 0)
                        cram_ref_decr(fd->refs, c->ref_seq_id);

                    if (!cram_get_ref(fd, bam_ref(b), 1, 0)) {
                        hts_log_error("Failed to load reference #%d", bam_ref(b));
                        free(MD.s);
                        return -1;
                    }

                    c->ref_seq_id = bam_ref(b); // overwritten later by -2
                    if (!fd->refs->ref_id[c->ref_seq_id]->seq)
                        return -1;
                    c->ref       = fd->refs->ref_id[c->ref_seq_id]->seq;
                    c->ref_start = 1;
                    c->ref_end   = fd->refs->ref_id[c->ref_seq_id]->length;
                }
            }

            if (process_one_read(fd, c, s, cr, b, r2, &MD) != 0) {
                free(MD.s);
                return -1;
            }

            if (first_base > cr->apos)
                first_base = cr->apos;

            if (last_base < cr->aend)
                last_base = cr->aend;
        }
        free(MD.s);

        // Process_one_read doesn't add read names as it can change
        // its mind during the loop on the CRAM_FLAG_DETACHED setting
        // of earlier records (if it detects the auto-generation of
        // TLEN is incorrect).  This affects which read-names can be
        // lossily compressed, so we do these in another pass.
        if (add_read_names(fd, c, s, r1_start) < 0)
            return -1;

        if (c->multi_seq) {
            s->hdr->ref_seq_id    = -2;
            s->hdr->ref_seq_start = 0;
            s->hdr->ref_seq_span  = 0;
        } else {
            s->hdr->ref_seq_id    = c->ref_id;
            s->hdr->ref_seq_start = first_base;
            s->hdr->ref_seq_span  = MAX(0, last_base - first_base + 1);
        }
        s->hdr->num_records = r2;

        // Processed a slice, now stash the aux blocks so the next
        // slice can start aggregating them from the start again.
        if (c->tags_used->n_occupied) {
            int ntags = c->tags_used->n_occupied;
            s->aux_block = calloc(ntags, sizeof(*s->aux_block));
            if (!s->aux_block)
                return -1;

            khint_t k;

            s->naux_block = 0;
            for (k = kh_begin(c->tags_used); k != kh_end(c->tags_used); k++) {
                if (!kh_exist(c->tags_used, k))
                    continue;

                cram_tag_map *tm = kh_val(c->tags_used, k);
                if (!tm) goto_err;
                if (!tm->blk) continue;
                s->aux_block[s->naux_block++] = tm->blk;
                tm->blk = NULL;
            }
            assert(s->naux_block <= c->tags_used->n_occupied);
        }
    }

    if (c->multi_seq && !fd->no_ref) {
        if (c->ref_seq_id >= 0)
            cram_ref_decr(fd->refs, c->ref_seq_id);
    }

    /* Link our bams[] array onto the spare bam list for reuse */
    spares = malloc(sizeof(*spares));
    if (!spares) goto_err;
    pthread_mutex_lock(&fd->bam_list_lock);
    spares->bams = c->bams;
    spares->next = fd->bl;
    fd->bl = spares;
    pthread_mutex_unlock(&fd->bam_list_lock);
    c->bams = NULL;

    /* Detect if a multi-seq container */
    cram_stats_encoding(fd, c->stats[DS_RI]);
    multi_ref = c->stats[DS_RI]->nvals > 1;
    pthread_mutex_lock(&fd->metrics_lock);
    fd->last_RI_count = c->stats[DS_RI]->nvals;
    pthread_mutex_unlock(&fd->metrics_lock);


    if (multi_ref) {
        hts_log_info("Multi-ref container");
        c->ref_seq_id = -2;
        c->ref_seq_start = 0;
        c->ref_seq_span = 0;
    }


    /* Compute MD5s */
    for (i = 0; i < c->curr_slice; i++) {
        cram_slice *s = c->slices[i];

        if (CRAM_MAJOR_VERS(fd->version) != 1) {
            if (s->hdr->ref_seq_id >= 0 && c->multi_seq == 0 && !fd->no_ref) {
                hts_md5_context *md5 = hts_md5_init();
                if (!md5)
                    return -1;
                hts_md5_update(md5,
                               c->ref + s->hdr->ref_seq_start - c->ref_start,
                               s->hdr->ref_seq_span);
                hts_md5_final(s->hdr->md5, md5);
                hts_md5_destroy(md5);
            } else {
                memset(s->hdr->md5, 0, 16);
            }
        }
    }

    c->num_records = 0;
    c->num_blocks = 1; // cram_block_compression_hdr
    c->length = 0;

    //fprintf(stderr, "=== BF ===\n");
    h->codecs[DS_BF] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_BF]),
                                         c->stats[DS_BF], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_BF]->nvals && !h->codecs[DS_BF]) goto_err;

    //fprintf(stderr, "=== CF ===\n");
    h->codecs[DS_CF] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_CF]),
                                         c->stats[DS_CF], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_CF]->nvals && !h->codecs[DS_CF]) goto_err;

    //fprintf(stderr, "=== RN ===\n");
    //h->codecs[DS_RN] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_RN]),
    //                                     c->stats[DS_RN], E_BYTE_ARRAY, NULL,
    //                                     fd->version);

    //fprintf(stderr, "=== AP ===\n");
    if (c->pos_sorted) {
        h->codecs[DS_AP] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_AP]),
                                             c->stats[DS_AP],
#ifdef LARGE_POS
                                             E_LONG,
#else
                                             E_INT,
#endif
                                             NULL, fd->version);
    } else {
        int p[2] = {0, c->max_apos};
        h->codecs[DS_AP] = cram_encoder_init(E_BETA, NULL, E_INT, p,
                                             fd->version);
    }
    if (!h->codecs[DS_AP]) goto_err;

    //fprintf(stderr, "=== RG ===\n");
    h->codecs[DS_RG] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_RG]),
                                         c->stats[DS_RG], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_RG]->nvals && !h->codecs[DS_RG]) goto_err;

    //fprintf(stderr, "=== MQ ===\n");
    h->codecs[DS_MQ] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_MQ]),
                                         c->stats[DS_MQ], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_MQ]->nvals && !h->codecs[DS_MQ]) goto_err;

    //fprintf(stderr, "=== NS ===\n");
    h->codecs[DS_NS] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_NS]),
                                         c->stats[DS_NS], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_NS]->nvals && !h->codecs[DS_NS]) goto_err;

    //fprintf(stderr, "=== MF ===\n");
    h->codecs[DS_MF] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_MF]),
                                         c->stats[DS_MF], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_MF]->nvals && !h->codecs[DS_MF]) goto_err;

    //fprintf(stderr, "=== TS ===\n");
    h->codecs[DS_TS] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_TS]),
                                         c->stats[DS_TS],
#ifdef LARGE_POS
                                         E_LONG,
#else
                                         E_INT,
#endif
                                         NULL, fd->version);
    if (c->stats[DS_TS]->nvals && !h->codecs[DS_TS]) goto_err;

    //fprintf(stderr, "=== NP ===\n");
    h->codecs[DS_NP] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_NP]),
                                         c->stats[DS_NP],
#ifdef LARGE_POS
                                         E_LONG,
#else
                                         E_INT,
#endif
                                         NULL, fd->version);
    if (c->stats[DS_NP]->nvals && !h->codecs[DS_NP]) goto_err;

    //fprintf(stderr, "=== NF ===\n");
    h->codecs[DS_NF] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_NF]),
                                         c->stats[DS_NF], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_NF]->nvals && !h->codecs[DS_NF]) goto_err;

    //fprintf(stderr, "=== RL ===\n");
    h->codecs[DS_RL] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_RL]),
                                         c->stats[DS_RL], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_RL]->nvals && !h->codecs[DS_RL]) goto_err;

    //fprintf(stderr, "=== FN ===\n");
    h->codecs[DS_FN] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_FN]),
                                         c->stats[DS_FN], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_FN]->nvals && !h->codecs[DS_FN]) goto_err;

    //fprintf(stderr, "=== FC ===\n");
    h->codecs[DS_FC] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_FC]),
                                         c->stats[DS_FC], E_BYTE, NULL,
                                         fd->version);
    if (c->stats[DS_FC]->nvals && !h->codecs[DS_FC]) goto_err;

    //fprintf(stderr, "=== FP ===\n");
    h->codecs[DS_FP] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_FP]),
                                         c->stats[DS_FP], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_FP]->nvals && !h->codecs[DS_FP]) goto_err;

    //fprintf(stderr, "=== DL ===\n");
    h->codecs[DS_DL] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_DL]),
                                         c->stats[DS_DL], E_INT, NULL,
                                         fd->version);
    if (c->stats[DS_DL]->nvals && !h->codecs[DS_DL]) goto_err;

    //fprintf(stderr, "=== BA ===\n");
    h->codecs[DS_BA] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_BA]),
                                         c->stats[DS_BA], E_BYTE, NULL,
                                         fd->version);
    if (c->stats[DS_BA]->nvals && !h->codecs[DS_BA]) goto_err;

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        cram_byte_array_len_encoder e;

        e.len_encoding = E_EXTERNAL;
        e.len_dat = (void *)DS_BB_len;
        //e.len_dat = (void *)DS_BB;

        e.val_encoding = E_EXTERNAL;
        e.val_dat = (void *)DS_BB;

        h->codecs[DS_BB] = cram_encoder_init(E_BYTE_ARRAY_LEN, NULL,
                                             E_BYTE_ARRAY, (void *)&e,
                                             fd->version);
        if (!h->codecs[DS_BB]) goto_err;
    } else {
        h->codecs[DS_BB] = NULL;
    }

    //fprintf(stderr, "=== BS ===\n");
    h->codecs[DS_BS] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_BS]),
                                         c->stats[DS_BS], E_BYTE, NULL,
                                         fd->version);
    if (c->stats[DS_BS]->nvals && !h->codecs[DS_BS]) goto_err;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        h->codecs[DS_TL] = NULL;
        h->codecs[DS_RI] = NULL;
        h->codecs[DS_RS] = NULL;
        h->codecs[DS_PD] = NULL;
        h->codecs[DS_HC] = NULL;
        h->codecs[DS_SC] = NULL;

        //fprintf(stderr, "=== TC ===\n");
        h->codecs[DS_TC] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_TC]),
                                             c->stats[DS_TC], E_BYTE, NULL,
                                             fd->version);
        if (c->stats[DS_TC]->nvals && !h->codecs[DS_TC]) goto_err;

        //fprintf(stderr, "=== TN ===\n");
        h->codecs[DS_TN] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_TN]),
                                             c->stats[DS_TN], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_TN]->nvals && !h->codecs[DS_TN]) goto_err;
    } else {
        h->codecs[DS_TC] = NULL;
        h->codecs[DS_TN] = NULL;

        //fprintf(stderr, "=== TL ===\n");
        h->codecs[DS_TL] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_TL]),
                                             c->stats[DS_TL], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_TL]->nvals && !h->codecs[DS_TL]) goto_err;


        //fprintf(stderr, "=== RI ===\n");
        h->codecs[DS_RI] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_RI]),
                                             c->stats[DS_RI], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_RI]->nvals && !h->codecs[DS_RI]) goto_err;

        //fprintf(stderr, "=== RS ===\n");
        h->codecs[DS_RS] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_RS]),
                                             c->stats[DS_RS], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_RS]->nvals && !h->codecs[DS_RS]) goto_err;

        //fprintf(stderr, "=== PD ===\n");
        h->codecs[DS_PD] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_PD]),
                                             c->stats[DS_PD], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_PD]->nvals && !h->codecs[DS_PD]) goto_err;

        //fprintf(stderr, "=== HC ===\n");
        h->codecs[DS_HC] = cram_encoder_init(cram_stats_encoding(fd, c->stats[DS_HC]),
                                             c->stats[DS_HC], E_INT, NULL,
                                             fd->version);
        if (c->stats[DS_HC]->nvals && !h->codecs[DS_HC]) goto_err;

        //fprintf(stderr, "=== SC ===\n");
        if (1) {
            int i2[2] = {0, DS_SC};

            h->codecs[DS_SC] = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
                                                 E_BYTE_ARRAY, (void *)i2,
                                                 fd->version);
        } else {
            // Appears to be no practical benefit to using this method,
            // but it may work better if we start mixing SC, IN and BB
            // elements into the same external block.
            cram_byte_array_len_encoder e;

            e.len_encoding = E_EXTERNAL;
            e.len_dat = (void *)DS_SC_len;

            e.val_encoding = E_EXTERNAL;
            e.val_dat = (void *)DS_SC;

            h->codecs[DS_SC] = cram_encoder_init(E_BYTE_ARRAY_LEN, NULL,
                                                 E_BYTE_ARRAY, (void *)&e,
                                                 fd->version);
        }
        if (!h->codecs[DS_SC]) goto_err;
    }

    //fprintf(stderr, "=== IN ===\n");
    {
        int i2[2] = {0, DS_IN};
        h->codecs[DS_IN] = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
                                             E_BYTE_ARRAY, (void *)i2,
                                             fd->version);
        if (!h->codecs[DS_IN]) goto_err;
    }

    h->codecs[DS_QS] = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE,
                                         (void *)DS_QS,
                                         fd->version);
    if (!h->codecs[DS_QS]) goto_err;
    {
        int i2[2] = {0, DS_RN};
        h->codecs[DS_RN] = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
                                             E_BYTE_ARRAY, (void *)i2,
                                             fd->version);
        if (!h->codecs[DS_RN]) goto_err;
    }


    /* Encode slices */
    for (i = 0; i < c->curr_slice; i++) {
        hts_log_info("Encode slice %d", i);

        if (cram_encode_slice(fd, c, h, c->slices[i]) != 0)
            return -1;
    }

    /* Create compression header */
    {
        h->ref_seq_id    = c->ref_seq_id;
        h->ref_seq_start = c->ref_seq_start;
        h->ref_seq_span  = c->ref_seq_span;
        h->num_records   = c->num_records;
        h->AP_delta = c->pos_sorted;
        memcpy(h->substitution_matrix, CRAM_SUBST_MATRIX, 20);

        if (!(c_hdr = cram_encode_compression_header(fd, c, h)))
            return -1;
    }

    /* Compute landmarks */
    /* Fill out slice landmarks */
    c->num_landmarks = c->curr_slice;
    c->landmark = malloc(c->num_landmarks * sizeof(*c->landmark));
    if (!c->landmark)
        return -1;

    /*
     * Slice offset starts after the first block, so we need to simulate
     * writing it to work out the correct offset
     */
    {
        slice_offset = c_hdr->method == RAW
            ? c_hdr->uncomp_size
            : c_hdr->comp_size;
        slice_offset += 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
            itf8_size(c_hdr->content_id) +
            itf8_size(c_hdr->comp_size) +
            itf8_size(c_hdr->uncomp_size);
    }

    c->ref_seq_id    = c->slices[0]->hdr->ref_seq_id;
    c->ref_seq_start = c->slices[0]->hdr->ref_seq_start;
    c->ref_seq_span  = c->slices[0]->hdr->ref_seq_span;
    for (i = 0; i < c->curr_slice; i++) {
        cram_slice *s = c->slices[i];

        c->num_blocks += s->hdr->num_blocks + 1; // slice header
        c->landmark[i] = slice_offset;

        if (s->hdr->ref_seq_start + s->hdr->ref_seq_span >
            c->ref_seq_start + c->ref_seq_span) {
            c->ref_seq_span = s->hdr->ref_seq_start + s->hdr->ref_seq_span
                - c->ref_seq_start;
        }

        slice_offset += s->hdr_block->method == RAW
            ? s->hdr_block->uncomp_size
            : s->hdr_block->comp_size;

        slice_offset += 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
            itf8_size(s->hdr_block->content_id) +
            itf8_size(s->hdr_block->comp_size) +
            itf8_size(s->hdr_block->uncomp_size);

        for (j = 0; j < s->hdr->num_blocks; j++) {
            slice_offset += 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
                itf8_size(s->block[j]->content_id) +
                itf8_size(s->block[j]->comp_size) +
                itf8_size(s->block[j]->uncomp_size);

            slice_offset += s->block[j]->method == RAW
                ? s->block[j]->uncomp_size
                : s->block[j]->comp_size;
        }
    }
    c->length += slice_offset; // just past the final slice

    c->comp_hdr_block = c_hdr;

    if (c->ref_seq_id >= 0) {
        cram_ref_decr(fd->refs, c->ref_seq_id);
    }

    /* Cache references up-front if we have unsorted access patterns */
    if (!fd->no_ref && c->refs_used) {
        for (i = 0; i < fd->refs->nref; i++) {
            if (c->refs_used[i])
                cram_ref_decr(fd->refs, i);
        }
    }

    return 0;

 err:
    return -1;
}


/*
 * Adds a feature code to a read within a slice. For purposes of minimising
 * memory allocations and fragmentation we have one array of features for all
 * reads within the slice. We return the index into this array for this new
 * feature.
 *
 * Returns feature index on success
 *         -1 on failure.
 */
static int cram_add_feature(cram_container *c, cram_slice *s,
                            cram_record *r, cram_feature *f) {
    if (s->nfeatures >= s->afeatures) {
        s->afeatures = s->afeatures ? s->afeatures*2 : 1024;
        s->features = realloc(s->features, s->afeatures*sizeof(*s->features));
        if (!s->features)
            return -1;
    }

    if (!r->nfeature++) {
        r->feature = s->nfeatures;
        if (cram_stats_add(c->stats[DS_FP], f->X.pos) < 0)
            return -1;
    } else {
        if (cram_stats_add(c->stats[DS_FP],
                           f->X.pos - s->features[r->feature + r->nfeature-2].X.pos) < 0)
            return -1;

    }
    if (cram_stats_add(c->stats[DS_FC], f->X.code) < 0)
        return -1;

    s->features[s->nfeatures++] = *f;

    return 0;
}

static int cram_add_substitution(cram_fd *fd, cram_container *c,
                                 cram_slice *s, cram_record *r,
                                 int pos, char base, char qual, char ref) {
    cram_feature f;

    // seq=ACGTN vs ref=ACGT or seq=ACGT vs ref=ACGTN
    if (fd->L2[(uc)base]<4 || (fd->L2[(uc)base]<5 && fd->L2[(uc)ref]<4)) {
        f.X.pos = pos+1;
        f.X.code = 'X';
        f.X.base = fd->cram_sub_matrix[ref&0x1f][base&0x1f];
        if (cram_stats_add(c->stats[DS_BS], f.X.base) < 0)
            return -1;
    } else {
        f.B.pos = pos+1;
        f.B.code = 'B';
        f.B.base = base;
        f.B.qual = qual;
        if (cram_stats_add(c->stats[DS_BA], f.B.base) < 0) return -1;
        if (cram_stats_add(c->stats[DS_QS], f.B.qual) < 0) return -1;
        BLOCK_APPEND_CHAR(s->qual_blk, qual);
    }
    return cram_add_feature(c, s, r, &f);

 block_err:
    return -1;
}

static int cram_add_bases(cram_fd *fd, cram_container *c,
                          cram_slice *s, cram_record *r,
                          int pos, int len, char *base) {
    cram_feature f;

    f.b.pos = pos+1;
    f.b.code = 'b';
    f.b.seq_idx = base - (char *)BLOCK_DATA(s->seqs_blk);
    f.b.len = len;

    return cram_add_feature(c, s, r, &f);
}

static int cram_add_base(cram_fd *fd, cram_container *c,
                         cram_slice *s, cram_record *r,
                         int pos, char base, char qual) {
    cram_feature f;
    f.B.pos = pos+1;
    f.B.code = 'B';
    f.B.base = base;
    f.B.qual = qual;
    if (cram_stats_add(c->stats[DS_BA], base) < 0) return -1;
    if (cram_stats_add(c->stats[DS_QS], qual) < 0) return -1;
    BLOCK_APPEND_CHAR(s->qual_blk, qual);
    return cram_add_feature(c, s, r, &f);

 block_err:
    return -1;
}

static int cram_add_quality(cram_fd *fd, cram_container *c,
                            cram_slice *s, cram_record *r,
                            int pos, char qual) {
    cram_feature f;
    f.Q.pos = pos+1;
    f.Q.code = 'Q';
    f.Q.qual = qual;
    if (cram_stats_add(c->stats[DS_QS], qual) < 0) return -1;
    BLOCK_APPEND_CHAR(s->qual_blk, qual);
    return cram_add_feature(c, s, r, &f);

 block_err:
    return -1;
}

static int cram_add_deletion(cram_container *c, cram_slice *s, cram_record *r,
                             int pos, int len, char *base) {
    cram_feature f;
    f.D.pos = pos+1;
    f.D.code = 'D';
    f.D.len = len;
    if (cram_stats_add(c->stats[DS_DL], len) < 0) return -1;
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_softclip(cram_container *c, cram_slice *s, cram_record *r,
                             int pos, int len, char *base, int version) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'S';
    f.S.len = len;
    switch (CRAM_MAJOR_VERS(version)) {
    case 1:
        f.S.seq_idx = BLOCK_SIZE(s->base_blk);
        BLOCK_APPEND(s->base_blk, base, len);
        BLOCK_APPEND_CHAR(s->base_blk, '\0');
        break;

    case 2:
    default:
        f.S.seq_idx = BLOCK_SIZE(s->soft_blk);
        if (base) {
            BLOCK_APPEND(s->soft_blk, base, len);
        } else {
            int i;
            for (i = 0; i < len; i++)
                BLOCK_APPEND_CHAR(s->soft_blk, 'N');
        }
        BLOCK_APPEND_CHAR(s->soft_blk, '\0');
        break;

        //default:
        //    // v3.0 onwards uses BB data-series
        //    f.S.seq_idx = BLOCK_SIZE(s->soft_blk);
    }
    return cram_add_feature(c, s, r, &f);

 block_err:
    return -1;
}

static int cram_add_hardclip(cram_container *c, cram_slice *s, cram_record *r,
                             int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'H';
    f.S.len = len;
    if (cram_stats_add(c->stats[DS_HC], len) < 0) return -1;
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_skip(cram_container *c, cram_slice *s, cram_record *r,
                         int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'N';
    f.S.len = len;
    if (cram_stats_add(c->stats[DS_RS], len) < 0) return -1;
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_pad(cram_container *c, cram_slice *s, cram_record *r,
                        int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'P';
    f.S.len = len;
    if (cram_stats_add(c->stats[DS_PD], len) < 0) return -1;
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_insertion(cram_container *c, cram_slice *s, cram_record *r,
                              int pos, int len, char *base) {
    cram_feature f;
    f.I.pos = pos+1;
    if (len == 1) {
        char b = base ? *base : 'N';
        f.i.code = 'i';
        f.i.base = b;
        if (cram_stats_add(c->stats[DS_BA], b) < 0) return -1;
    } else {
        f.I.code = 'I';
        f.I.len = len;
        f.S.seq_idx = BLOCK_SIZE(s->base_blk);
        if (base) {
            BLOCK_APPEND(s->base_blk, base, len);
        } else {
            int i;
            for (i = 0; i < len; i++)
                BLOCK_APPEND_CHAR(s->base_blk, 'N');
        }
        BLOCK_APPEND_CHAR(s->base_blk, '\0');
    }
    return cram_add_feature(c, s, r, &f);

 block_err:
    return -1;
}

/*
 * Encodes auxiliary data, CRAM 1.0 format.
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present (FIXME)
 */
static char *cram_encode_aux_1_0(cram_fd *fd, bam_seq_t *b, cram_container *c,
                                 cram_slice *s, cram_record *cr) {
    char *aux, *tmp, *rg = NULL;
    int aux_size = bam_blk_size(b) -
        ((char *)bam_aux(b) - (char *)&bam_ref(b));

    /* Worst case is 1 nul char on every ??:Z: string, so +33% */
    BLOCK_GROW(s->aux_blk, aux_size*1.34+1);
    tmp = (char *)BLOCK_END(s->aux_blk);

    aux = (char *)bam_aux(b);
    cr->TN_idx = s->nTN;

    while (aux[0] != 0) {
        int32_t i32;
        int r;

        if (aux[0] == 'R' && aux[1] == 'G' && aux[2] == 'Z') {
            rg = &aux[3];
            while (*aux++);
            continue;
        }
        if (aux[0] == 'M' && aux[1] == 'D' && aux[2] == 'Z') {
            while (*aux++);
            continue;
        }
        if (aux[0] == 'N' && aux[1] == 'M') {
            switch(aux[2]) {
            case 'A': case 'C': case 'c': aux+=4; break;
            case 'I': case 'i': case 'f': aux+=7; break;
            default:
                hts_log_error("Unhandled type code for NM tag");
                return NULL;
            }
            continue;
        }

        cr->ntags++;

        i32 = (aux[0]<<16) | (aux[1]<<8) | aux[2];
        kh_put(m_tagmap, c->tags_used, i32, &r);
        if (-1 == r)
            return NULL;

        if (s->nTN >= s->aTN) {
            s->aTN = s->aTN ? s->aTN*2 : 1024;
            if (!(s->TN = realloc(s->TN, s->aTN * sizeof(*s->TN))))
                return NULL;
        }
        s->TN[s->nTN++] = i32;
        if (cram_stats_add(c->stats[DS_TN], i32) < 0)
            goto block_err;

        switch(aux[2]) {
        case 'A': case 'C': case 'c':
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            *tmp++=*aux++;
            break;

        case 'S': case 's':
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            *tmp++=*aux++; *tmp++=*aux++;
            break;

        case 'I': case 'i': case 'f':
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            break;

        case 'd':
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            break;

        case 'Z': case 'H':
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            while ((*tmp++=*aux++));
            *tmp++ = '\t'; // stop byte
            break;

        case 'B': {
            int type = aux[3], blen;
            uint32_t count = (uint32_t)((((unsigned char *)aux)[4]<< 0) +
                                        (((unsigned char *)aux)[5]<< 8) +
                                        (((unsigned char *)aux)[6]<<16) +
                                        (((unsigned char *)aux)[7]<<24));
            // skip TN field
            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

            // We use BYTE_ARRAY_LEN with external length, so store that first
            switch (type) {
            case 'c': case 'C':
                blen = count;
                break;
            case 's': case 'S':
                blen = 2*count;
                break;
            case 'i': case 'I': case 'f':
                blen = 4*count;
                break;
            default:
                hts_log_error("Unknown sub-type '%c' for aux type 'B'", type);
                return NULL;
            }

            tmp += itf8_put(tmp, blen+5);

            *tmp++=*aux++; // sub-type & length
            *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

            // The tag data itself
            memcpy(tmp, aux, blen); tmp += blen; aux += blen;

            //cram_stats_add(c->aux_B_stats, blen);
            break;
        }
        default:
            hts_log_error("Unknown aux type '%c'", aux[2]);
            return NULL;
        }
    }
    if (cram_stats_add(c->stats[DS_TC], cr->ntags) < 0)
        goto block_err;

    cr->aux = BLOCK_SIZE(s->aux_blk);
    cr->aux_size = (uc *)tmp - (BLOCK_DATA(s->aux_blk) + cr->aux);
    BLOCK_SIZE(s->aux_blk) = (uc *)tmp - BLOCK_DATA(s->aux_blk);
    assert(s->aux_blk->byte <= s->aux_blk->alloc);

    return rg;

 block_err:
    return NULL;
}

/*
 * Encodes auxiliary data. Largely duplicated from above, but done so to
 * keep it simple and avoid a myriad of version ifs.
 *
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present, also sets "*err" to non-zero
 */
static char *cram_encode_aux(cram_fd *fd, bam_seq_t *b, cram_container *c,
                             cram_slice *s, cram_record *cr,
                             int verbatim_NM, int verbatim_MD,
                             int NM, kstring_t *MD,
                             int *err) {
    char *aux, *orig, *rg = NULL;
    int aux_size = bam_get_l_aux(b);
    cram_block *td_b = c->comp_hdr->TD_blk;
    int TD_blk_size = BLOCK_SIZE(td_b), new;
    char *key;
    khint_t k;

    if (err) *err = 1;

    orig = aux = (char *)bam_aux(b);

    // Copy aux keys to td_b and aux values to slice aux blocks
    while (aux - orig < aux_size && aux[0] != 0) {
        int r;

        // RG:Z
        if (aux[0] == 'R' && aux[1] == 'G' && aux[2] == 'Z') {
            rg = &aux[3];
            while (*aux++);
            continue;
        }

        // MD:Z
        if (aux[0] == 'M' && aux[1] == 'D' && aux[2] == 'Z') {
            if (cr->len && !fd->no_ref && !(cr->flags & BAM_FUNMAP) && !verbatim_MD) {
                if (MD && MD->s && strncasecmp(MD->s, aux+3, orig + aux_size - (aux+3)) == 0) {
                    while (*aux++);
                    continue;
                }
            }
        }

        // NM:i
        if (aux[0] == 'N' && aux[1] == 'M') {
            if (cr->len && !fd->no_ref && !(cr->flags & BAM_FUNMAP) && !verbatim_NM) {
                int NM_ = bam_aux2i((uint8_t *)aux+2);
                if (NM_ == NM) {
                    switch(aux[2]) {
                    case 'A': case 'C': case 'c': aux+=4; break;
                    case 'S': case 's':           aux+=5; break;
                    case 'I': case 'i': case 'f': aux+=7; break;
                    default:
                        hts_log_error("Unhandled type code for NM tag");
                        return NULL;
                    }
                    continue;
                }
            }
        }

        BLOCK_APPEND(td_b, aux, 3);

        // Container level tags_used, for TD series
        // Maps integer key ('X0i') to cram_tag_map struct.
        int key = (aux[0]<<16)|(aux[1]<<8)|aux[2];
        k = kh_put(m_tagmap, c->tags_used, key, &r);
        if (-1 == r)
            return NULL;
        else if (r != 0)
            kh_val(c->tags_used, k) = NULL;

        if (r == 1) {
            khint_t k_global;

            // Global tags_used for cram_metrics support
            pthread_mutex_lock(&fd->metrics_lock);
            k_global = kh_put(m_metrics, fd->tags_used, key, &r);
            if (-1 == r) {
                pthread_mutex_unlock(&fd->metrics_lock);
                return NULL;
            }
            if (r >= 1) {
                kh_val(fd->tags_used, k_global) = cram_new_metrics();
                if (!kh_val(fd->tags_used, k_global)) {
                    kh_del(m_metrics, fd->tags_used, k_global);
                    pthread_mutex_unlock(&fd->metrics_lock);
                    goto err;
                }
            }

            pthread_mutex_unlock(&fd->metrics_lock);

            int i2[2] = {'\t',key};
            size_t sk = key;
            cram_tag_map *m = calloc(1, sizeof(*m));
            if (!m)
                goto_err;
            kh_val(c->tags_used, k) = m;

            cram_codec *c;

            // Use a block content id based on the tag id.
            // Codec type depends on tag data type.
            switch(aux[2]) {
            case 'Z': case 'H':
                // string as byte_array_stop
                c = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
                                      E_BYTE_ARRAY, (void *)i2,
                                      fd->version);
                break;

            case 'A': case 'c': case 'C': {
                // byte array len, 1 byte
                cram_byte_array_len_encoder e;
                cram_stats st;

                e.len_encoding = E_HUFFMAN;
                e.len_dat = NULL;
                memset(&st, 0, sizeof(st));
                if (cram_stats_add(&st, 1) < 0) goto block_err;
                cram_stats_encoding(fd, &st);

                e.val_encoding = E_EXTERNAL;
                e.val_dat = (void *)sk;

                c = cram_encoder_init(E_BYTE_ARRAY_LEN, &st,
                                      E_BYTE_ARRAY, (void *)&e,
                                      fd->version);
                break;
            }

            case 's': case 'S': {
                // byte array len, 2 byte
                cram_byte_array_len_encoder e;
                cram_stats st;

                e.len_encoding = E_HUFFMAN;
                e.len_dat = NULL;
                memset(&st, 0, sizeof(st));
                if (cram_stats_add(&st, 2) < 0) goto block_err;
                cram_stats_encoding(fd, &st);

                e.val_encoding = E_EXTERNAL;
                e.val_dat = (void *)sk;

                c = cram_encoder_init(E_BYTE_ARRAY_LEN, &st,
                                      E_BYTE_ARRAY, (void *)&e,
                                      fd->version);
                break;
            }
            case 'i': case 'I': case 'f': {
                // byte array len, 4 byte
                cram_byte_array_len_encoder e;
                cram_stats st;

                e.len_encoding = E_HUFFMAN;
                e.len_dat = NULL;
                memset(&st, 0, sizeof(st));
                if (cram_stats_add(&st, 4) < 0) goto block_err;
                cram_stats_encoding(fd, &st);

                e.val_encoding = E_EXTERNAL;
                e.val_dat = (void *)sk;

                c = cram_encoder_init(E_BYTE_ARRAY_LEN, &st,
                                      E_BYTE_ARRAY, (void *)&e,
                                      fd->version);
                break;
            }

            case 'B': {
                // Byte array of variable size, but we generate our tag
                // byte stream at the wrong stage (during reading and not
                // after slice header construction). So we use
                // BYTE_ARRAY_LEN with the length codec being external
                // too.
                cram_byte_array_len_encoder e;

                e.len_encoding = E_EXTERNAL;
                e.len_dat = (void *)sk; // or key+128 for len?

                e.val_encoding = E_EXTERNAL;
                e.val_dat = (void *)sk;

                c = cram_encoder_init(E_BYTE_ARRAY_LEN, NULL,
                                      E_BYTE_ARRAY, (void *)&e,
                                      fd->version);
                break;
            }

            default:
                hts_log_error("Unsupported SAM aux type '%c'", aux[2]);
                c = NULL;
            }

            if (!c)
                goto_err;

            m->codec = c;

            // Link to fd-global tag metrics
            pthread_mutex_lock(&fd->metrics_lock);
            m->m = k_global ? (cram_metrics *)kh_val(fd->tags_used, k_global) : NULL;
            pthread_mutex_unlock(&fd->metrics_lock);
        }

        cram_tag_map *tm = (cram_tag_map *)kh_val(c->tags_used, k);
        if (!tm) goto_err;
        cram_codec *codec = tm->codec;
        if (!tm->codec) goto_err;

        switch(aux[2]) {
        case 'A': case 'C': case 'c':
            if (!tm->blk) {
                if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                    return NULL;
                codec->u.e_byte_array_len.val_codec->out = tm->blk;
            }

            aux+=3;
            //codec->encode(s, codec, aux, 1);
            // Functionally equivalent, but less code.
            BLOCK_APPEND_CHAR(tm->blk, *aux);
            aux++;
            break;

        case 'S': case 's':
            if (!tm->blk) {
                if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                    return NULL;
                codec->u.e_byte_array_len.val_codec->out = tm->blk;
            }

            aux+=3;
            //codec->encode(s, codec, aux, 2);
            BLOCK_APPEND(tm->blk, aux, 2);
            aux+=2;
            break;

        case 'I': case 'i': case 'f':
            if (!tm->blk) {
                if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                    return NULL;
                codec->u.e_byte_array_len.val_codec->out = tm->blk;
            }

            aux+=3;
            //codec->encode(s, codec, aux, 4);
            BLOCK_APPEND(tm->blk, aux, 4);
            aux+=4;
            break;

        case 'd':
            if (!tm->blk) {
                if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                    return NULL;
                codec->u.e_byte_array_len.val_codec->out = tm->blk;
            }

            aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
            //codec->encode(s, codec, aux, 8);
            BLOCK_APPEND(tm->blk, aux, 8);
            aux+=8;
            break;

        case 'Z': case 'H':
            {
                if (!tm->blk) {
                    if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                        return NULL;
                    codec->out = tm->blk;
                }

                char *aux_s;
                aux += 3;
                aux_s = aux;
                while (*aux++);
                if (codec->encode(s, codec, aux_s, aux - aux_s) < 0)
                    return NULL;
            }
            break;

        case 'B': {
            int type = aux[3], blen;
            uint32_t count = (uint32_t)((((unsigned char *)aux)[4]<< 0) +
                                        (((unsigned char *)aux)[5]<< 8) +
                                        (((unsigned char *)aux)[6]<<16) +
                                        (((unsigned char *)aux)[7]<<24));
            if (!tm->blk) {
                if (!(tm->blk = cram_new_block(EXTERNAL, key)))
                    return NULL;
                codec->u.e_byte_array_len.len_codec->out = tm->blk;
                codec->u.e_byte_array_len.val_codec->out = tm->blk;
            }

            // skip TN field
            aux+=3;

            // We use BYTE_ARRAY_LEN with external length, so store that first
            switch (type) {
            case 'c': case 'C':
                blen = count;
                break;
            case 's': case 'S':
                blen = 2*count;
                break;
            case 'i': case 'I': case 'f':
                blen = 4*count;
                break;
            default:
                hts_log_error("Unknown sub-type '%c' for aux type 'B'", type);
                return NULL;
            }

            blen += 5; // sub-type & length

            if (codec->encode(s, codec, aux, blen) < 0)
                return NULL;
            aux += blen;
            break;
        }
        default:
            hts_log_error("Unknown aux type '%c'", aux[2]);
            return NULL;
        }
        tm->blk->m = tm->m;
    }

    // FIXME: sort BLOCK_DATA(td_b) by char[3] triples

    // And and increment TD hash entry
    BLOCK_APPEND_CHAR(td_b, 0);

    // Duplicate key as BLOCK_DATA() can be realloced to a new pointer.
    key = string_ndup(c->comp_hdr->TD_keys,
                      (char *)BLOCK_DATA(td_b) + TD_blk_size,
                      BLOCK_SIZE(td_b) - TD_blk_size);
    if (!key)
        goto block_err;
    k = kh_put(m_s2i, c->comp_hdr->TD_hash, key, &new);
    if (new < 0) {
        return NULL;
    } else if (new == 0) {
        BLOCK_SIZE(td_b) = TD_blk_size;
    } else {
        kh_val(c->comp_hdr->TD_hash, k) = c->comp_hdr->nTL;
        c->comp_hdr->nTL++;
    }

    cr->TL = kh_val(c->comp_hdr->TD_hash, k);
    if (cram_stats_add(c->stats[DS_TL], cr->TL) < 0)
        goto block_err;

    if (err) *err = 0;
    return rg;

 err:
 block_err:
    return NULL;
}

/*
 * During cram_next_container or before the final flush at end of
 * file, we update the current slice headers and increment the slice
 * number to the next slice.
 *
 * See cram_next_container() and cram_close().
 */
void cram_update_curr_slice(cram_container *c) {
    cram_slice *s = c->slice;
    if (c->multi_seq) {
        s->hdr->ref_seq_id    = -2;
        s->hdr->ref_seq_start = 0;
        s->hdr->ref_seq_span  = 0;
    } else {
        s->hdr->ref_seq_id    = c->curr_ref;
        s->hdr->ref_seq_start = c->first_base;
        s->hdr->ref_seq_span  = MAX(0, c->last_base - c->first_base + 1);
    }
    s->hdr->num_records   = c->curr_rec;

    if (c->curr_slice == 0) {
        if (c->ref_seq_id != s->hdr->ref_seq_id)
            c->ref_seq_id  = s->hdr->ref_seq_id;
        c->ref_seq_start = c->first_base;
    }

    c->curr_slice++;
}

/*
 * Handles creation of a new container or new slice, flushing any
 * existing containers when appropriate.
 *
 * Really this is next slice, which may or may not lead to a new container.
 *
 * Returns cram_container pointer on success
 *         NULL on failure.
 */
static cram_container *cram_next_container(cram_fd *fd, bam_seq_t *b) {
    cram_container *c = fd->ctr;
    int i;

    /* First occurrence */
    if (c->curr_ref == -2)
        c->curr_ref = bam_ref(b);

    if (c->slice)
        cram_update_curr_slice(c);

    /* Flush container */
    if (c->curr_slice == c->max_slice ||
        (bam_ref(b) != c->curr_ref && !c->multi_seq)) {
        c->ref_seq_span = fd->last_base - c->ref_seq_start + 1;
        hts_log_info("Flush container %d/%"PRId64"..%"PRId64,
                     c->ref_seq_id, c->ref_seq_start,
                     c->ref_seq_start + c->ref_seq_span -1);

        /* Encode slices */
        if (-1 == cram_flush_container_mt(fd, c))
            return NULL;
        if (!fd->pool) {
            // Move to sep func, as we need cram_flush_container for
            // the closing phase to flush the partial container.
            for (i = 0; i < c->max_slice; i++) {
                cram_free_slice(c->slices[i]);
                c->slices[i] = NULL;
            }

            c->slice = NULL;
            c->curr_slice = 0;

            /* Easy approach for purposes of freeing stats */
            cram_free_container(c);
        }

        c = fd->ctr = cram_new_container(fd->seqs_per_slice,
                                         fd->slices_per_container);
        if (!c)
            return NULL;
        c->record_counter = fd->record_counter;
        c->curr_ref = bam_ref(b);
    }

    c->last_pos = c->first_base = c->last_base = bam_pos(b)+1;

    /* New slice */
    c->slice = c->slices[c->curr_slice] =
        cram_new_slice(MAPPED_SLICE, c->max_rec);
    if (!c->slice)
        return NULL;

    if (c->multi_seq) {
        c->slice->hdr->ref_seq_id = -2;
        c->slice->hdr->ref_seq_start = 0;
        c->slice->last_apos = 1;
    } else {
        c->slice->hdr->ref_seq_id = bam_ref(b);
        // wrong for unsorted data, will fix during encoding.
        c->slice->hdr->ref_seq_start = bam_pos(b)+1;
        c->slice->last_apos = bam_pos(b)+1;
    }

    c->curr_rec = 0;
    c->s_num_bases = 0;
    c->n_mapped = 0;

    return c;
}

/*
 * Converts a single bam record into a cram record.
 * Possibly used within a thread.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int process_one_read(cram_fd *fd, cram_container *c,
                            cram_slice *s, cram_record *cr,
                            bam_seq_t *b, int rnum, kstring_t *MD) {
    int i, fake_qual = -1, NM = 0;
    char *cp, *rg;
    char *ref, *seq, *qual;

    // Any places with N in seq and/or reference can lead to ambiguous
    // interpretation of the SAM NM:i tag.  So we store these verbatim
    // to ensure valid data round-trips the same regardless of who
    // defines it as valid.
    // Similarly when alignments go beyond end of the reference.
    int verbatim_NM = fd->store_nm;
    int verbatim_MD = fd->store_md;

    // FIXME: multi-ref containers

    ref = c->ref;
    cr->flags       = bam_flag(b);
    cr->len         = bam_seq_len(b);
    if (!bam_aux_get(b, "MD"))
        MD = NULL;
    else
        MD->l = 0;

    //fprintf(stderr, "%s => %d\n", rg ? rg : "\"\"", cr->rg);

    cr->ref_id      = bam_ref(b);
    if (cram_stats_add(c->stats[DS_RI], cr->ref_id) < 0)
        goto block_err;
    if (cram_stats_add(c->stats[DS_BF], fd->cram_flag_swap[cr->flags & 0xfff]) < 0)
        goto block_err;

    // Non reference based encoding means storing the bases verbatim as features, which in
    // turn means every base also has a quality already stored.
    if (!fd->no_ref || CRAM_MAJOR_VERS(fd->version) >= 3)
        cr->cram_flags |= CRAM_FLAG_PRESERVE_QUAL_SCORES;

    if (cr->len <= 0 && CRAM_MAJOR_VERS(fd->version) >= 3)
        cr->cram_flags |= CRAM_FLAG_NO_SEQ;
    //cram_stats_add(c->stats[DS_CF], cr->cram_flags & CRAM_FLAG_MASK);

    c->num_bases   += cr->len;
    cr->apos        = bam_pos(b)+1;
    if (c->pos_sorted) {
        if (cr->apos < s->last_apos) {
            c->pos_sorted = 0;
        } else {
            if (cram_stats_add(c->stats[DS_AP], cr->apos - s->last_apos) < 0)
                goto block_err;
            s->last_apos = cr->apos;
        }
    } else {
        //cram_stats_add(c->stats[DS_AP], cr->apos);
    }
    c->max_apos += (cr->apos > c->max_apos) * (cr->apos - c->max_apos);

    /*
     * This seqs_ds is largely pointless and it could reuse the same memory
     * over and over.
     * s->base_blk is what we need for encoding.
     */
    cr->seq         = BLOCK_SIZE(s->seqs_blk);
    cr->qual        = BLOCK_SIZE(s->qual_blk);
    BLOCK_GROW(s->seqs_blk, cr->len+1);
    BLOCK_GROW(s->qual_blk, cr->len);

    // Convert BAM nibble encoded sequence to string of base pairs
    seq = cp = (char *)BLOCK_END(s->seqs_blk);
    *seq = 0;
    nibble2base(bam_seq(b), cp, cr->len);
    BLOCK_SIZE(s->seqs_blk) += cr->len;

    qual = cp = (char *)bam_qual(b);

    /* Copy and parse */
    if (!(cr->flags & BAM_FUNMAP)) {
        uint32_t *cig_to, *cig_from;
        int64_t apos = cr->apos-1, spos = 0;
        int64_t MD_last = apos; // last position of edit in MD tag

        cr->cigar       = s->ncigar;
        cr->ncigar      = bam_cigar_len(b);
        while (cr->cigar + cr->ncigar >= s->cigar_alloc) {
            s->cigar_alloc = s->cigar_alloc ? s->cigar_alloc*2 : 1024;
            s->cigar = realloc(s->cigar, s->cigar_alloc * sizeof(*s->cigar));
            if (!s->cigar)
                return -1;
        }

        cig_to = (uint32_t *)s->cigar;
        cig_from = (uint32_t *)bam_cigar(b);

        cr->feature = 0;
        cr->nfeature = 0;
        for (i = 0; i < cr->ncigar; i++) {
            enum cigar_op cig_op = cig_from[i] & BAM_CIGAR_MASK;
            uint32_t cig_len = cig_from[i] >> BAM_CIGAR_SHIFT;
            cig_to[i] = cig_from[i];

            /* Can also generate events from here for CRAM diffs */

            switch (cig_op) {
                int l;

                // Don't trust = and X ops to be correct.
            case BAM_CMATCH:
            case BAM_CBASE_MATCH:
            case BAM_CBASE_MISMATCH:
                //fprintf(stderr, "\nBAM_CMATCH\nR: %.*s\nS: %.*s\n",
                //      cig_len, &ref[apos], cig_len, &seq[spos]);
                l = 0;
                if (!fd->no_ref && cr->len) {
                    int end = cig_len+apos < c->ref_end
                        ? cig_len : c->ref_end - apos;
                    char *sp = &seq[spos];
                    char *rp = &ref[apos];
                    char *qp = &qual[spos];
                    if (end > cr->len) {
                        hts_log_error("CIGAR and query sequence are of different length");
                        return -1;
                    }
                    for (l = 0; l < end; l++) {
                        // This case is just too disputed and different tools
                        // interpret these in different ways.  We give up and
                        // store verbatim.
                        if (rp[l] == 'N' && sp[l] == 'N')
                            verbatim_NM = verbatim_MD = 1;
                        if (rp[l] != sp[l]) {
                            // Build our own MD tag if one is on the sequence, so
                            // we can ensure it matches and thus can be discarded.
                            if (MD && ref) {
                                if (kputuw(apos+l - MD_last, MD) < 0) goto err;
                                if (kputc(rp[l], MD) < 0) goto err;
                                MD_last = apos+l+1;
                            }
                            NM++;
                            if (!sp[l])
                                break;
                            if (0 && CRAM_MAJOR_VERS(fd->version) >= 3) {
                                // Disabled for the time being as it doesn't
                                // seem to gain us much.
                                int ol=l;
                                while (l<end && rp[l] != sp[l])
                                    l++;
                                if (l-ol > 1) {
                                    if (cram_add_bases(fd, c, s, cr, spos+ol,
                                                       l-ol, &seq[spos+ol]))
                                        return -1;
                                    l--;
                                } else {
                                    l = ol;
                                    if (cram_add_substitution(fd, c, s, cr,
                                                              spos+l, sp[l],
                                                              qp[l], rp[l]))
                                        return -1;
                                }
                            } else {
                                if (cram_add_substitution(fd, c, s, cr, spos+l,
                                                          sp[l], qp[l], rp[l]))
                                    return -1;
                            }
                        }
                    }
                    spos += l;
                    apos += l;
                }

                if (l < cig_len && cr->len) {
                    if (fd->no_ref) {
                        if (CRAM_MAJOR_VERS(fd->version) == 3) {
                            if (cram_add_bases(fd, c, s, cr, spos,
                                               cig_len-l, &seq[spos]))
                                return -1;
                            spos += cig_len-l;
                        } else {
                            for (; l < cig_len && seq[spos]; l++, spos++) {
                                if (cram_add_base(fd, c, s, cr, spos,
                                                  seq[spos], qual[spos]))
                                    return -1;
                            }
                        }
                    } else {
                        /* off end of sequence or non-ref based output */
                        verbatim_NM = verbatim_MD = 1;
                        for (; l < cig_len && seq[spos]; l++, spos++) {
                            if (cram_add_base(fd, c, s, cr, spos,
                                              seq[spos], qual[spos]))
                                return -1;
                        }
                    }
                    apos += cig_len;
                } else if (!cr->len) {
                    /* Seq "*" */
                    verbatim_NM = verbatim_MD = 1;
                    apos += cig_len;
                    spos += cig_len;
                }
                break;

            case BAM_CDEL:
                if (MD && ref) {
                    if (kputuw(apos - MD_last, MD) < 0) goto err;
                    if (apos < c->ref_end) {
                        if (kputc_('^', MD) < 0) goto err;
                        if (kputsn(&ref[apos], MIN(c->ref_end - apos, cig_len), MD) < 0)
                            goto err;
                    }
                }
                NM += cig_len;

                if (cram_add_deletion(c, s, cr, spos, cig_len, &seq[spos]))
                    return -1;
                apos += cig_len;
                MD_last = apos;
                break;

            case BAM_CREF_SKIP:
                if (cram_add_skip(c, s, cr, spos, cig_len, &seq[spos]))
                    return -1;
                apos += cig_len;
                MD_last += cig_len;
                break;

            case BAM_CINS:
                if (cram_add_insertion(c, s, cr, spos, cig_len,
                                       cr->len ? &seq[spos] : NULL))
                    return -1;
                if (fd->no_ref && cr->len) {
                    for (l = 0; l < cig_len; l++, spos++) {
                        cram_add_quality(fd, c, s, cr, spos, qual[spos]);
                    }
                } else {
                    spos += cig_len;
                }
                NM += cig_len;
                break;

            case BAM_CSOFT_CLIP:
                if (cram_add_softclip(c, s, cr, spos, cig_len,
                                      cr->len ? &seq[spos] : NULL,
                                      fd->version))
                    return -1;

                if (fd->no_ref &&
                    !(cr->cram_flags & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
                    if (cr->len) {
                        for (l = 0; l < cig_len; l++, spos++) {
                            cram_add_quality(fd, c, s, cr, spos, qual[spos]);
                        }
                    } else {
                        for (l = 0; l < cig_len; l++, spos++) {
                            cram_add_quality(fd, c, s, cr, spos, -1);
                        }
                    }
                } else {
                    spos += cig_len;
                }
                break;

            case BAM_CHARD_CLIP:
                if (cram_add_hardclip(c, s, cr, spos, cig_len, &seq[spos]))
                    return -1;
                break;

            case BAM_CPAD:
                if (cram_add_pad(c, s, cr, spos, cig_len, &seq[spos]))
                    return -1;
                break;

            default:
                hts_log_error("Unknown CIGAR op code %d", cig_op);
                return -1;
            }
        }
        if (cr->len && spos != cr->len) {
            hts_log_error("CIGAR and query sequence are of different length");
            return -1;
        }
        fake_qual = spos;
        cr->aend = fd->no_ref ? apos : MIN(apos, c->ref_end);
        if (cram_stats_add(c->stats[DS_FN], cr->nfeature) < 0)
            goto block_err;

        if (MD && ref)
            if (kputuw(apos - MD_last, MD) < 0) goto err;
    } else {
        // Unmapped
        cr->cram_flags |= CRAM_FLAG_PRESERVE_QUAL_SCORES;
        cr->cigar  = 0;
        cr->ncigar = 0;
        cr->nfeature = 0;
        cr->aend = cr->apos;
        for (i = 0; i < cr->len; i++)
            if (cram_stats_add(c->stats[DS_BA], seq[i]) < 0)
                goto block_err;
        fake_qual = 0;
    }

    cr->ntags      = 0; //cram_stats_add(c->stats[DS_TC], cr->ntags);
    int err = 0;
    if (CRAM_MAJOR_VERS(fd->version) == 1)
        rg = cram_encode_aux_1_0(fd, b, c, s, cr);
    else
        rg = cram_encode_aux(fd, b, c, s, cr, verbatim_NM, verbatim_MD, NM, MD, &err);
    if (err)
        goto block_err;

    /* Read group, identified earlier */
    if (rg) {
        sam_hrec_rg_t *brg = sam_hrecs_find_rg(fd->header->hrecs, rg);
        cr->rg = brg ? brg->id : -1;
    } else if (CRAM_MAJOR_VERS(fd->version) == 1) {
        sam_hrec_rg_t *brg = sam_hrecs_find_rg(fd->header->hrecs, "UNKNOWN");
        assert(brg);
    } else {
        cr->rg = -1;
    }
    if (cram_stats_add(c->stats[DS_RG], cr->rg) < 0)
        goto block_err;

    /*
     * Append to the qual block now. We do this here as
     * cram_add_substitution() can generate BA/QS events which need to
     * be in the qual block before we append the rest of the data.
     */
    if (cr->cram_flags & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
        /* Special case of seq "*" */
        if (cr->len == 0) {
            cr->len = fake_qual;
            BLOCK_GROW(s->qual_blk, cr->len);
            cp = (char *)BLOCK_END(s->qual_blk);
            memset(cp, 255, cr->len);
        } else {
            BLOCK_GROW(s->qual_blk, cr->len);
            cp = (char *)BLOCK_END(s->qual_blk);
            char *from = (char *)&bam_qual(b)[0];
            char *to = &cp[0];
            memcpy(to, from, cr->len);
            //for (i = 0; i < cr->len; i++) cp[i] = from[i];
        }
        BLOCK_SIZE(s->qual_blk) += cr->len;
    } else {
        if (cr->len == 0)
            cr->len = fake_qual >= 0 ? fake_qual : cr->aend - cr->apos + 1;
    }

    if (cram_stats_add(c->stats[DS_RL], cr->len) < 0)
        goto block_err;

    /* Now we know apos and aend both, update mate-pair information */
    {
        int new;
        khint_t k;
        int sec = (cr->flags & BAM_FSECONDARY) ? 1 : 0;

        //fprintf(stderr, "Checking %"PRId64"/%.*s\t", rnum,
        //      cr->name_len, DSTRING_STR(s->name_ds)+cr->name);
        if (cr->flags & BAM_FPAIRED) {
            char *key = string_ndup(s->pair_keys, bam_name(b), bam_name_len(b));
            if (!key)
                return -1;

            k = kh_put(m_s2i, s->pair[sec], key, &new);
            if (-1 == new)
                return -1;
            else if (new > 0)
                kh_val(s->pair[sec], k) = rnum;
        } else {
            new = 1;
            k = 0; // Prevents false-positive warning from gcc -Og
        }

        if (new == 0) {
            cram_record *p = &s->crecs[kh_val(s->pair[sec], k)];
            int64_t aleft, aright;
            int sign;

            aleft = MIN(cr->apos, p->apos);
            aright = MAX(cr->aend, p->aend);
            if (cr->apos < p->apos) {
                sign = 1;
            } else if (cr->apos > p->apos) {
                sign = -1;
            } else if (cr->flags & BAM_FREAD1) {
                sign = 1;
            } else {
                sign = -1;
            }

            // This vs p: tlen, matepos, flags. Permit TLEN 0 and/or TLEN +/-
            // a small amount, if appropriate options set.
            if ((bam_ins_size(b) &&
                 llabs(bam_ins_size(b) - sign*(aright-aleft+1)) > fd->tlen_approx) ||
                (!bam_ins_size(b) && !fd->tlen_zero))
                goto detached;

            if ((!fd->tlen_zero && MAX(bam_mate_pos(b)+1, 0) != p->apos) &&
                !(fd->tlen_zero && bam_mate_pos(b) == 0))
                goto detached;

            if (((bam_flag(b) & BAM_FMUNMAP) != 0) !=
                ((p->flags & BAM_FUNMAP) != 0))
                goto detached;

            if (((bam_flag(b) & BAM_FMREVERSE) != 0) !=
                ((p->flags & BAM_FREVERSE) != 0))
                goto detached;


            // p vs this: tlen, matepos, flags
            if (p->ref_id != cr->ref_id &&
                !(fd->tlen_zero && p->ref_id == -1))
                goto detached;

            if ((p->tlen && llabs(p->tlen - -sign*(aright-aleft+1)) > fd->tlen_approx) ||
                (!p->tlen && !fd->tlen_zero))
                goto detached;

            if (p->mate_pos != cr->apos &&
                !(fd->tlen_zero && p->mate_pos == 0))
                goto detached;

            if (((p->flags & BAM_FMUNMAP) != 0) !=
                ((p->mate_flags & CRAM_M_UNMAP) != 0))
                goto detached;

            if (((p->flags & BAM_FMREVERSE) != 0) !=
                ((p->mate_flags & CRAM_M_REVERSE) != 0))
                goto detached;

            // Supplementary reads are just too ill defined
            if ((cr->flags & BAM_FSUPPLEMENTARY) ||
                (p->flags & BAM_FSUPPLEMENTARY))
                goto detached;

            // When in lossy name mode, if a read isn't detached we
            // cannot store the name.  The corollary is that when we
            // must store the name, it must be detached (inefficient).
            if (fd->lossy_read_names &&
                (!(cr->cram_flags & CRAM_FLAG_DISCARD_NAME) ||
                 !((p->cram_flags & CRAM_FLAG_DISCARD_NAME))))
                goto detached;

            /*
             * The fields below are unused when encoding this read as it is
             * no longer detached.  In theory they may get referred to when
             * processing a 3rd or 4th read in this template?, so we set them
             * here just to be sure.
             *
             * They do not need cram_stats_add() calls those as they are
             * not emitted.
             */
            cr->mate_pos = p->apos;
            cr->tlen = sign*(aright-aleft+1);
            cr->mate_flags =
                ((p->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)   * CRAM_M_UNMAP +
                ((p->flags & BAM_FMREVERSE) == BAM_FMREVERSE) * CRAM_M_REVERSE;

            // Decrement statistics aggregated earlier
            if (p->cram_flags & CRAM_FLAG_STATS_ADDED) {
                cram_stats_del(c->stats[DS_NP], p->mate_pos);
                cram_stats_del(c->stats[DS_MF], p->mate_flags);
                cram_stats_del(c->stats[DS_TS], p->tlen);
                cram_stats_del(c->stats[DS_NS], p->mate_ref_id);
            }

            /* Similarly we could correct the p-> values too, but these will no
             * longer have any code that refers back to them as the new 'p'
             * for this template is our current 'cr'.
             */
            //p->mate_pos = cr->apos;
            //p->mate_flags =
            //  ((cr->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)  * CRAM_M_UNMAP +
            //  ((cr->flags & BAM_FMREVERSE) == BAM_FMREVERSE)* CRAM_M_REVERSE;
            //p->tlen = p->apos - cr->aend;

            // Clear detached from cr flags
            cr->cram_flags &= ~CRAM_FLAG_DETACHED;
            if (cram_stats_add(c->stats[DS_CF], cr->cram_flags & CRAM_FLAG_MASK) < 0)
                goto block_err;

            // Clear detached from p flags and set downstream
            if (p->cram_flags & CRAM_FLAG_STATS_ADDED) {
                cram_stats_del(c->stats[DS_CF], p->cram_flags & CRAM_FLAG_MASK);
                p->cram_flags &= ~CRAM_FLAG_STATS_ADDED;
            }

            p->cram_flags  &= ~CRAM_FLAG_DETACHED;
            p->cram_flags  |=  CRAM_FLAG_MATE_DOWNSTREAM;
            if (cram_stats_add(c->stats[DS_CF], p->cram_flags & CRAM_FLAG_MASK) < 0)
                goto block_err;

            p->mate_line = rnum - (kh_val(s->pair[sec], k) + 1);
            if (cram_stats_add(c->stats[DS_NF], p->mate_line) < 0)
                goto block_err;

            kh_val(s->pair[sec], k) = rnum;
        } else {
        detached:
            //fprintf(stderr, "unpaired\n");

            /* Derive mate flags from this flag */
            cr->mate_flags = 0;
            if (bam_flag(b) & BAM_FMUNMAP)
                cr->mate_flags |= CRAM_M_UNMAP;
            if (bam_flag(b) & BAM_FMREVERSE)
                cr->mate_flags |= CRAM_M_REVERSE;

            if (cram_stats_add(c->stats[DS_MF], cr->mate_flags) < 0)
                goto block_err;

            cr->mate_pos    = MAX(bam_mate_pos(b)+1, 0);
            if (cram_stats_add(c->stats[DS_NP], cr->mate_pos) < 0)
                goto block_err;

            cr->tlen        = bam_ins_size(b);
            if (cram_stats_add(c->stats[DS_TS], cr->tlen) < 0)
                goto block_err;

            cr->cram_flags |= CRAM_FLAG_DETACHED;
            if (cram_stats_add(c->stats[DS_CF], cr->cram_flags & CRAM_FLAG_MASK) < 0)
                goto block_err;
            if (cram_stats_add(c->stats[DS_NS], bam_mate_ref(b)) < 0)
                goto block_err;

            cr->cram_flags |= CRAM_FLAG_STATS_ADDED;
        }
    }

    cr->mqual       = bam_map_qual(b);
    if (cram_stats_add(c->stats[DS_MQ], cr->mqual) < 0)
        goto block_err;

    cr->mate_ref_id = bam_mate_ref(b);

    if (!(bam_flag(b) & BAM_FUNMAP)) {
        if (c->first_base > cr->apos)
            c->first_base = cr->apos;

        if (c->last_base < cr->aend)
            c->last_base = cr->aend;
    }

    return 0;

 block_err:
 err:
    return -1;
}

/*
 * Write iterator: put BAM format sequences into a CRAM file.
 * We buffer up a containers worth of data at a time.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b) {
    cram_container *c;

    if (!fd->ctr) {
        fd->ctr = cram_new_container(fd->seqs_per_slice,
                                     fd->slices_per_container);
        if (!fd->ctr)
            return -1;
        fd->ctr->record_counter = fd->record_counter;
    }
    c = fd->ctr;

    if (!c->slice || c->curr_rec == c->max_rec ||
        (bam_ref(b) != c->curr_ref && c->curr_ref >= -1) ||
        (c->s_num_bases >= fd->bases_per_slice)) {
        int slice_rec, curr_rec, multi_seq = fd->multi_seq == 1;
        int curr_ref = c->slice ? c->curr_ref : bam_ref(b);

        /*
         * Start packing slices when we routinely have under 1/4tr full.
         *
         * This option isn't available if we choose to embed references
         * since we can only have one per slice.
         *
         * The multi_seq var here refers to our intention for the next slice.
         * This slice has already been encoded so we output as-is.
         */
        if (fd->multi_seq == -1 && c->curr_rec < c->max_rec/4+10 &&
            fd->last_slice && fd->last_slice < c->max_rec/4+10 &&
            !fd->embed_ref) {
            if (!c->multi_seq)
                hts_log_info("Multi-ref enabled for next container");
            multi_seq = 1;
        } else if (fd->multi_seq == 1) {
            pthread_mutex_lock(&fd->metrics_lock);
            if (fd->last_RI_count <= c->max_slice && fd->multi_seq_user != 1) {
                multi_seq = 0;
                hts_log_info("Multi-ref disabled for next container");
            }
            pthread_mutex_unlock(&fd->metrics_lock);
        }

        slice_rec = c->slice_rec;
        curr_rec  = c->curr_rec;

        if (CRAM_MAJOR_VERS(fd->version) == 1 ||
            c->curr_rec == c->max_rec || fd->multi_seq != 1 || !c->slice ||
            c->s_num_bases >= fd->bases_per_slice) {
            if (NULL == (c = cram_next_container(fd, b))) {
                if (fd->ctr) {
                    // prevent cram_close attempting to flush
                    cram_free_container(fd->ctr);
                    fd->ctr = NULL;
                }
                return -1;
            }
        }

        /*
         * Due to our processing order, some things we've already done we
         * cannot easily undo. So when we first notice we should be packing
         * multiple sequences per container we emit the small partial
         * container as-is and then start a fresh one in a different mode.
         */
        if (multi_seq == 0 && fd->multi_seq == 1 && fd->multi_seq_user == -1) {
            // User selected auto-mode, we're currently using multi-seq, but
            // have detected we don't need to.  Switch back to auto.
            fd->multi_seq = -1;
        } else if (multi_seq) {
            // We detected we need multi-seq
            fd->multi_seq = 1;
            c->multi_seq = 1;
            c->pos_sorted = 0; // required atm for multi_seq slices

            if (!c->refs_used) {
                pthread_mutex_lock(&fd->ref_lock);
                c->refs_used = calloc(fd->refs->nref, sizeof(int));
                pthread_mutex_unlock(&fd->ref_lock);
                if (!c->refs_used)
                    return -1;
            }
        }

        fd->last_slice = curr_rec - slice_rec;
        c->slice_rec = c->curr_rec;

        // Have we seen this reference before?
        if (bam_ref(b) >= 0 && curr_ref >= 0 && bam_ref(b) != curr_ref && !fd->embed_ref &&
            !fd->unsorted && multi_seq) {

            if (!c->refs_used) {
                pthread_mutex_lock(&fd->ref_lock);
                c->refs_used = calloc(fd->refs->nref, sizeof(int));
                pthread_mutex_unlock(&fd->ref_lock);
                if (!c->refs_used)
                    return -1;
            } else if (c->refs_used && c->refs_used[bam_ref(b)]) {
                pthread_mutex_lock(&fd->ref_lock);
                fd->unsorted = 1;
                pthread_mutex_unlock(&fd->ref_lock);
                fd->multi_seq = 1;
            }
        }

        c->curr_ref = bam_ref(b);
        if (c->refs_used && c->curr_ref >= 0) c->refs_used[c->curr_ref]++;
    }

    if (!c->bams) {
        /* First time through, allocate a set of bam pointers */
        pthread_mutex_lock(&fd->bam_list_lock);
        if (fd->bl) {
            spare_bams *spare = fd->bl;
            c->bams = spare->bams;
            fd->bl = spare->next;
            free(spare);
        } else {
            c->bams = calloc(c->max_c_rec, sizeof(bam_seq_t *));
            if (!c->bams) {
                pthread_mutex_unlock(&fd->bam_list_lock);
                return -1;
            }
        }
        pthread_mutex_unlock(&fd->bam_list_lock);
    }

    /* Copy or alloc+copy the bam record, for later encoding */
    if (c->bams[c->curr_c_rec]) {
        if (bam_copy1(c->bams[c->curr_c_rec], b) == NULL)
            return -1;
    } else {
        c->bams[c->curr_c_rec] = bam_dup1(b);
        if (c->bams[c->curr_c_rec] == NULL)
            return -1;
    }
    c->curr_rec++;
    c->curr_c_rec++;
    c->s_num_bases += bam_seq_len(b);
    c->n_mapped += (bam_flag(b) & BAM_FUNMAP) ? 0 : 1;
    fd->record_counter++;

    return 0;
}
