/*
Copyright (c) 2015, 2018-2020, 2022-2023 Genome Research Ltd.
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

/*! \file
 * External CRAM interface.
 *
 * Internally we're happy to use macros and to grub around in the cram
 * structures.  This isn't very sustainable for an externally usable
 * ABI though, so we have anonymous structs and accessor functions too
 * to permit software such as samtools reheader to manipulate cram
 * containers and blocks in a robust manner.
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
#include <stdint.h>

#if defined(HAVE_EXTERNAL_LIBHTSCODECS)
#include <htscodecs/rANS_static4x16.h.h>
#else
#include "../htscodecs/htscodecs/rANS_static4x16.h"
#endif

#include "../htslib/hfile.h"
#include "cram.h"

/*
 *-----------------------------------------------------------------------------
 * cram_fd
 */
sam_hdr_t *cram_fd_get_header(cram_fd *fd) { return fd->header; }
void cram_fd_set_header(cram_fd *fd, sam_hdr_t *hdr) { fd->header = hdr; }

int cram_fd_get_version(cram_fd *fd) { return fd->version; }
void cram_fd_set_version(cram_fd *fd, int vers) { fd->version = vers; }

int cram_major_vers(cram_fd *fd) { return CRAM_MAJOR_VERS(fd->version); }
int cram_minor_vers(cram_fd *fd) { return CRAM_MINOR_VERS(fd->version); }

hFILE *cram_fd_get_fp(cram_fd *fd) { return fd->fp; }
void cram_fd_set_fp(cram_fd *fd, hFILE *fp) { fd->fp = fp; }


/*
 *-----------------------------------------------------------------------------
 * cram_container
 */
int32_t cram_container_get_length(cram_container *c) {
    return c->length;
}

void cram_container_set_length(cram_container *c, int32_t length) {
    c->length = length;
}


int32_t cram_container_get_num_blocks(cram_container *c) {
    return c->num_blocks;
}

void cram_container_set_num_blocks(cram_container *c, int32_t num_blocks) {
    c->num_blocks = num_blocks;
}

int32_t cram_container_get_num_records(cram_container *c) {
    return c->num_records;
}

int64_t cram_container_get_num_bases(cram_container *c) {
    return c->num_bases;
}


/* Returns the landmarks[] array and the number of elements
 * in num_landmarks.
 */
int32_t *cram_container_get_landmarks(cram_container *c, int32_t *num_landmarks) {
    *num_landmarks = c->num_landmarks;
    return c->landmark;
}

/* Sets the landmarks[] array (pointer copy, not a memory dup) and
 * num_landmarks value.
 */
void cram_container_set_landmarks(cram_container *c, int32_t num_landmarks,
                                  int32_t *landmarks) {
    c->num_landmarks = num_landmarks;
    c->landmark = landmarks;
}


/* Returns true if the container is empty (EOF marker) */
int cram_container_is_empty(cram_fd *fd) {
    return fd->empty_container;
}


/*
 *-----------------------------------------------------------------------------
 * cram_block_compression_hdr
 */

/*
 * Utility function to edit an RG id.
 * This is only possible if there is one single RG value used and it
 * is in the container compression header using HUFFMAN or BETA
 * codec.  In this case it is essentially hard coded and needs no
 * editing of external (or worse, CORE) blocks.
 *
 * Returns 0 on success
 *        -1 on failure
 */
// Or arbitrary set compression header constant?

static int cram_block_compression_hdr_set_DS(cram_block_compression_hdr *ch,
                                             int ds, int new_rg) {
    if (!ch || !ch->codecs[ds])
        return -1;

    switch (ch->codecs[ds]->codec) {
    case E_HUFFMAN:
        if (ch->codecs[ds]->u.huffman.ncodes != 1)
            return -1;
        ch->codecs[ds]->u.huffman.codes[0].symbol = new_rg;
        return 0;

    case E_BETA:
        if (ch->codecs[ds]->u.beta.nbits != 0)
            return -1;
        ch->codecs[ds]->u.beta.offset = -new_rg;
        return 0;

    default:
        break;
    }

    return -1;
}

int cram_block_compression_hdr_set_rg(cram_block_compression_hdr *ch, int new_rg) {
    return cram_block_compression_hdr_set_DS(ch, DS_RG, new_rg);
}

/*
 * Converts a cram_block_compression_hdr struct used for decoding to
 * one used for encoding.  Maybe this should be a transparent
 * operation applied on-demand.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_block_compression_hdr_decoder2encoder(cram_fd *fd,
                                               cram_block_compression_hdr *ch) {
    int i;

    if (!ch)
        return -1;

    for (i = 0; i < DS_END; i++) {
        cram_codec *co = ch->codecs[i];
        if (!co)
            continue;

        if (-1 == cram_codec_decoder2encoder(fd, co))
            return -1;
    }

    return 0;
}

typedef struct {
    cram_block_compression_hdr *hdr;
    cram_map *curr_map;
    int idx;
    int is_tag; // phase 2 using tag_encoding_map
} cram_codec_iter;

static void cram_codec_iter_init(cram_block_compression_hdr *hdr,
                                 cram_codec_iter *iter) {
    iter->hdr = hdr;
    iter->curr_map = NULL;
    iter->idx = 0;
    iter->is_tag = 0;
}

// See enum cram_DS_ID in cram/cram_structs
static int cram_ds_to_key(enum cram_DS_ID ds) {
    switch(ds) {
    case DS_RN: return 256*'R'+'N';
    case DS_QS: return 256*'Q'+'S';
    case DS_IN: return 256*'I'+'N';
    case DS_SC: return 256*'S'+'C';
    case DS_BF: return 256*'B'+'F';
    case DS_CF: return 256*'C'+'F';
    case DS_AP: return 256*'A'+'P';
    case DS_RG: return 256*'R'+'G';
    case DS_MQ: return 256*'M'+'Q';
    case DS_NS: return 256*'N'+'S';
    case DS_MF: return 256*'M'+'F';
    case DS_TS: return 256*'T'+'S';
    case DS_NP: return 256*'N'+'P';
    case DS_NF: return 256*'N'+'F';
    case DS_RL: return 256*'R'+'L';
    case DS_FN: return 256*'F'+'N';
    case DS_FC: return 256*'F'+'C';
    case DS_FP: return 256*'F'+'P';
    case DS_DL: return 256*'D'+'L';
    case DS_BA: return 256*'B'+'A';
    case DS_BS: return 256*'B'+'S';
    case DS_TL: return 256*'T'+'L';
    case DS_RI: return 256*'R'+'I';
    case DS_RS: return 256*'R'+'S';
    case DS_PD: return 256*'P'+'D';
    case DS_HC: return 256*'H'+'C';
    case DS_BB: return 256*'B'+'B';
    case DS_QQ: return 256*'Q'+'Q';
    case DS_TN: return 256*'T'+'N';
    case DS_TC: return 256*'T'+'C';
    case DS_TM: return 256*'T'+'M';
    case DS_TV: return 256*'T'+'V';
    default: break;
    }

    return -1; // unknown
}

static cram_codec *cram_codec_iter_next(cram_codec_iter *iter,
                                        int *key) {
    cram_codec *cc = NULL;
    cram_block_compression_hdr *hdr = iter->hdr;

    if (!iter->is_tag) {
        // 1: Iterating through main data-series
        do {
            cc = hdr->codecs[iter->idx++];
        } while(!cc && iter->idx < DS_END);
        if (cc) {
            *key = cram_ds_to_key(iter->idx-1);
            return cc;
        }

        // Reset index for phase 2
        iter->idx = 0;
        iter->is_tag = 1;
    }

    do {
        if (!iter->curr_map)
            iter->curr_map = hdr->tag_encoding_map[iter->idx++];

        cc = iter->curr_map ? iter->curr_map->codec : NULL;
        if (cc) {
            *key = iter->curr_map->key;
            iter->curr_map = iter->curr_map->next;
            return cc;
        }
    } while (iter->idx <= CRAM_MAP_HASH);

    // End of codecs
    return NULL;
}

/*
 * A list of data-series, used to create a linked list threaded through
 * a single array.
 */
typedef struct ds_list {
    int data_series;
    int next;
} ds_list;

KHASH_MAP_INIT_INT(cid, int64_t)

// Opaque struct for the CRAM block content-id -> data-series map.
struct cram_cid2ds_t {
    ds_list *ds;        // array of data-series with linked lists threading through it
    int ds_size;
    int ds_idx;
    khash_t(cid) *hash; // key=content_id,  value=index to ds array
    int *ds_a;          // serialised array of data-series returned by queries.
};

void cram_cid2ds_free(cram_cid2ds_t *cid2ds) {
    if (cid2ds) {
        if (cid2ds->hash)
            kh_destroy(cid, cid2ds->hash);
        free(cid2ds->ds);
        free(cid2ds->ds_a);
        free(cid2ds);
    }
}

/*
 * Map cram block numbers to data-series.  It's normally a 1:1 mapping,
 * but in rare cases it can be 1:many (or even many:many).
 * The key is the block number and the value is an index into the data-series
 * array, which we iterate over until reaching a negative value.
 *
 * Provide cid2ds as NULL to allocate a new map or pass in an existing one
 * to append to this map.  The new (or existing) map is returned.
 *
 * Returns the cid2ds (newly allocated or as provided) on success,
 *         NULL on failure.
 */
cram_cid2ds_t *cram_update_cid2ds_map(cram_block_compression_hdr *hdr,
                                      cram_cid2ds_t *cid2ds) {
    cram_cid2ds_t *c2d = cid2ds;
    if (!c2d) {
        c2d = calloc(1, sizeof(*c2d));
        if (!c2d)
            return NULL;

        c2d->hash = kh_init(cid);
        if (!c2d->hash)
            goto err;
    }

    // Iterate through codecs.  Initially primary two-left ones in
    // rec_encoding_map, and then the three letter in tag_encoding_map.
    cram_codec_iter citer;
    cram_codec_iter_init(hdr, &citer);
    cram_codec *codec;
    int key;

    while ((codec = cram_codec_iter_next(&citer, &key))) {
        // Having got a codec, we can then use cram_codec_to_id to get
        // the block IDs utilised by that codec.  This is then our
        // map for allocating data blocks to data series, but for shared
        // blocks we can't separate out how much is used by each DS.
        int bnum[2];
        cram_codec_get_content_ids(codec, bnum);

        khiter_t k;
        int ret, i;
        for (i = 0; i < 2; i++) {
            if (bnum[i] > -2) {
                k = kh_put(cid, c2d->hash, bnum[i], &ret);
                if (ret < 0)
                    goto err;

                if (c2d->ds_idx >= c2d->ds_size) {
                    c2d->ds_size += 100;
                    c2d->ds_size *= 2;
                    ds_list *ds_new = realloc(c2d->ds,
                                              c2d->ds_size * sizeof(*ds_new));
                    if (!ds_new)
                        goto err;
                    c2d->ds = ds_new;
                }

                if (ret == 0) {
                    // Shared content_id, so add to list of DS

                    // Maybe data-series should be part of the hash key?
                    //
                    // So top-32 bit is content-id, bot-32 bit is key.
                    // Sort hash by key and then can group all the data-series
                    // known together. ??
                    //
                    // Brute force for now, scan to see if recorded.
                    // Typically this is minimal effort as we almost always
                    // have 1 data-series per block content-id, so the list to
                    // search is of size 1.
                    int dsi = kh_value(c2d->hash, k);
                    while (dsi >= 0) {
                        if (c2d->ds[dsi].data_series == key)
                            break;
                        dsi = c2d->ds[dsi].next;
                    }

                    if (dsi == -1) {
                        // Block content_id seen before, but not with this DS
                        c2d->ds[c2d->ds_idx].data_series = key;
                        c2d->ds[c2d->ds_idx].next = kh_value(c2d->hash, k);
                        kh_value(c2d->hash, k) = c2d->ds_idx;
                        c2d->ds_idx++;
                    }
                } else {
                    // First time this content id has been used
                    c2d->ds[c2d->ds_idx].data_series = key;
                    c2d->ds[c2d->ds_idx].next = -1;
                    kh_value(c2d->hash, k) = c2d->ds_idx;
                    c2d->ds_idx++;
                }
            }
        }
    }

    return c2d;

 err:
    if (c2d != cid2ds)
        cram_cid2ds_free(c2d);
    return NULL;
}

/*
 * Return a list of data series observed as belonging to a block with
 * the specified content_id.  *n is the number of data series
 * returned, or 0 if block is unused.
 * Block content_id of -1 is used to indicate the CORE block.
 *
 * The pointer returned is owned by the cram_cid2ds state and should
 * not be freed by the caller.
 */
int *cram_cid2ds_query(cram_cid2ds_t *c2d, int content_id, int *n) {
    *n = 0;
    if (!c2d || !c2d->hash)
        return NULL;

    khiter_t k = kh_get(cid, c2d->hash, content_id);
    if (k == kh_end(c2d->hash))
        return NULL;

    if (!c2d->ds_a) {
        c2d->ds_a = malloc(c2d->ds_idx * sizeof(int));
        if (!c2d->ds_a)
            return NULL;
    }

    int dsi = kh_value(c2d->hash, k); // initial ds array index from hash
    int idx = 0;
    while (dsi >= 0) {
        c2d->ds_a[idx++] = c2d->ds[dsi].data_series;
        dsi = c2d->ds[dsi].next;      // iterate over list within ds array
    }

    *n = idx;
    return c2d->ds_a;
}

/*
 * Produces a description of the record and tag encodings held within
 * a compression header and appends to 'ks'.
 *
 * Returns 0 on success,
 *        <0 on failure.
 */
int cram_describe_encodings(cram_block_compression_hdr *hdr, kstring_t *ks) {
    cram_codec_iter citer;
    cram_codec_iter_init(hdr, &citer);
    cram_codec *codec;
    int key, r = 0;

    while ((codec = cram_codec_iter_next(&citer, &key))) {
        char key_s[4] = {0};
        int key_i = 0;
        if (key>>16) key_s[key_i++] = key>>16;
        key_s[key_i++] = (key>>8)&0xff;
        key_s[key_i++] = key&0xff;
        r |= ksprintf(ks, "\t%s\t", key_s) < 0;
        r |= cram_codec_describe(codec, ks) < 0;
        r |= kputc('\n', ks) < 0;
    }

    return r ? -1 : 0;
}

/*
 *-----------------------------------------------------------------------------
 * cram_slice
 */
int32_t cram_slice_hdr_get_num_blocks(cram_block_slice_hdr *hdr) {
    return hdr->num_blocks;
}

int cram_slice_hdr_get_embed_ref_id(cram_block_slice_hdr *h) {
    return h->ref_base_id;
}

void cram_slice_hdr_get_coords(cram_block_slice_hdr *h,
                               int *refid, hts_pos_t *start, hts_pos_t *span) {
    if (refid)
        *refid = h->ref_seq_id;
    if (start)
        *start = h->ref_seq_start;
    if (span)
        *span  = h->ref_seq_span;
}

/*
 *-----------------------------------------------------------------------------
 * cram_block
 */
int32_t cram_block_get_content_id(cram_block *b)  {
    return b->content_type == CORE ? -1 : b->content_id;
}
int32_t cram_block_get_comp_size(cram_block *b)   { return b->comp_size; }
int32_t cram_block_get_uncomp_size(cram_block *b) { return b->uncomp_size; }
int32_t cram_block_get_crc32(cram_block *b)       { return b->crc32; }
void *  cram_block_get_data(cram_block *b)        { return BLOCK_DATA(b); }
int32_t cram_block_get_size(cram_block *b)        { return BLOCK_SIZE(b); }
enum cram_block_method cram_block_get_method(cram_block *b) {
    return (enum cram_block_method)b->orig_method;
}
enum cram_content_type cram_block_get_content_type(cram_block *b) {
    return b->content_type;
}

void cram_block_set_content_id(cram_block *b, int32_t id) { b->content_id = id; }
void cram_block_set_comp_size(cram_block *b, int32_t size) { b->comp_size = size; }
void cram_block_set_uncomp_size(cram_block *b, int32_t size) { b->uncomp_size = size; }
void cram_block_set_crc32(cram_block *b, int32_t crc) { b->crc32 = crc; }
void cram_block_set_data(cram_block *b, void *data) { BLOCK_DATA(b) = data; }
void cram_block_set_size(cram_block *b, int32_t size) { BLOCK_SIZE(b) = size; }

int cram_block_append(cram_block *b, const void *data, int size) {
    BLOCK_APPEND(b, data, size);
    return 0;

 block_err:
    return -1;
}
void cram_block_update_size(cram_block *b) { BLOCK_UPLEN(b); }

// Offset is known as "size" internally, but it can be confusing.
size_t cram_block_get_offset(cram_block *b) { return BLOCK_SIZE(b); }
void cram_block_set_offset(cram_block *b, size_t offset) { BLOCK_SIZE(b) = offset; }

/*
 * Given a compressed block of data in a specified compression method,
 * fill out the 'cm' field with meta-data gleaned from the compressed
 * block.
 *
 * If comp is CRAM_COMP_UNKNOWN, we attempt to auto-detect the compression
 * format, but this doesn't work for all methods.
 *
 * Retuns the detected or specified comp method, and fills out *cm
 * if non-NULL.
 */
cram_method_details *cram_expand_method(uint8_t *data, int32_t size,
                                        enum cram_block_method comp) {
    cram_method_details *cm = calloc(1, sizeof(*cm));
    if (!cm)
        return NULL;

    const char *xz_header = "\xFD""7zXZ"; // including nul

    if (comp == CRAM_COMP_UNKNOWN) {
        // Auto-detect
        if (size > 1 && data[0] == 0x1f && data[1] == 0x8b)
            comp = CRAM_COMP_GZIP;
        else if (size > 3 && data[1] == 'B' && data[2] == 'Z'
                 && data[3] == 'h')
            comp = CRAM_COMP_BZIP2;
        else if (size > 6 && memcmp(xz_header, data, 6) == 0)
            comp = CRAM_COMP_LZMA;
        else
            comp = CRAM_COMP_UNKNOWN;
    }
    cm->method = comp;

    // Interrogate the compressed data stream to fill out additional fields.
    switch (comp) {
    case CRAM_COMP_GZIP:
        if (size > 8) {
            if (data[8] == 4)
                cm->level = 1;
            else if (data[8] == 2)
                cm->level = 9;
            else
                cm->level = 5;
        }
        break;

    case CRAM_COMP_BZIP2:
        if (size > 3 && data[3] >= '1' && data[3] <= '9')
            cm->level = data[3]-'0';
        break;

    case CRAM_COMP_RANS4x8:
        cm->Nway = 4;
        if (size > 0 && data[0] == 1)
            cm->order = 1;
        else
            cm->order = 0;
        break;

    case CRAM_COMP_RANSNx16:
        if (size > 0) {
            cm->order  = data[0] & 1;
            cm->Nway   = data[0] & RANS_ORDER_X32    ? 32 : 4;
            cm->rle    = data[0] & RANS_ORDER_RLE    ?  1 : 0;
            cm->pack   = data[0] & RANS_ORDER_PACK   ?  1 : 0;
            cm->cat    = data[0] & RANS_ORDER_CAT    ?  1 : 0;
            cm->stripe = data[0] & RANS_ORDER_STRIPE ?  1 : 0;
            cm->nosz   = data[0] & RANS_ORDER_NOSZ   ?  1 : 0;
        }
        break;

    case CRAM_COMP_ARITH:
        if (size > 0) {
            // Not in a public header, but the same transforms as rANSNx16
            cm->order  = data[0] & 3;
            cm->rle    = data[0] & RANS_ORDER_RLE    ?  1 : 0;
            cm->pack   = data[0] & RANS_ORDER_PACK   ?  1 : 0;
            cm->cat    = data[0] & RANS_ORDER_CAT    ?  1 : 0;
            cm->stripe = data[0] & RANS_ORDER_STRIPE ?  1 : 0;
            cm->nosz   = data[0] & RANS_ORDER_NOSZ   ?  1 : 0;
            cm->ext    = data[0] & 4 /*external*/    ?  1 : 0;
        }
        break;

    case CRAM_COMP_TOK3:
        if (size > 8) {
            if (data[8] == 1)
                cm->level = 11;
            else if (data[8] == 0)
                cm->level = 1;
        }
        break;

    default:
        break;
    }

    return cm;
}

/*
 *-----------------------------------------------------------------------------
 * cram_codecs
 */

// -2 is unused.
// -1 is CORE
// >= 0 is the block with that Content ID
void cram_codec_get_content_ids(cram_codec *c, int ids[2]) {
    ids[0] = cram_codec_to_id(c, &ids[1]);
}

/*
 *-----------------------------------------------------------------------------
 * Utility functions
 */

/*
 * Copies the blocks representing the next num_slice slices from a
 * container from 'in' to 'out'.  It is expected that the file pointer
 * is just after the read of the cram_container and cram compression
 * header.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_copy_slice(cram_fd *in, cram_fd *out, int32_t num_slice) {
    int32_t i, j;

    for (i = 0; i < num_slice; i++) {
        cram_block *blk;
        cram_block_slice_hdr *hdr;

        if (!(blk = cram_read_block(in)))
            return -1;
        if (!(hdr = cram_decode_slice_header(in, blk))) {
            cram_free_block(blk);
            return -1;
        }
        if (cram_write_block(out, blk) != 0) {
            cram_free_block(blk);
            return -1;
        }
        cram_free_block(blk);

        int num_blocks = cram_slice_hdr_get_num_blocks(hdr);
        for (j = 0; j < num_blocks; j++) {
            blk = cram_read_block(in);
            if (!blk || cram_write_block(out, blk) != 0) {
                if (blk) cram_free_block(blk);
                return -1;
            }
            cram_free_block(blk);
        }
        cram_free_slice_header(hdr);
    }

    return 0;
}

/*
 * Renumbers RG numbers in a cram compression header.
 *
 * CRAM stores RG as the Nth number in the header, rather than a
 * string holding the ID: tag.  This is smaller in space, but means
 * "samtools cat" to join files together that contain single but
 * different RG lines needs a way of renumbering them.
 *
 * The file descriptor is expected to be immediately after the
 * cram_container structure (ie before the cram compression header).
 * Due to the nature of the CRAM format, this needs to read and write
 * the blocks itself.  Note that there may be multiple slices within
 * the container, meaning multiple compression headers to manipulate.
 * Changing RG may change the size of the compression header and
 * therefore the length field in the container.  Hence we rewrite all
 * blocks just in case and also emit the adjusted container.
 *
 * The current implementation can only cope with renumbering a single
 * RG (and only then if it is using HUFFMAN or BETA codecs).  In
 * theory it *may* be possible to renumber multiple RGs if they use
 * HUFFMAN to the CORE block or use an external block unshared by any
 * other data series.  So we have an API that can be upgraded to
 * support this, but do not implement it for now.  An example
 * implementation of RG as an EXTERNAL block would be to find that
 * block and rewrite it, returning the number of blocks consumed.
 *
 * Returns 0 on success;
 *        -1 if unable to edit;
 *        -2 on other errors (eg I/O).
 */
int cram_transcode_rg(cram_fd *in, cram_fd *out,
                      cram_container *c,
                      int nrg, int *in_rg, int *out_rg) {
    int new_rg = *out_rg, old_size, new_size;
    cram_block *o_blk, *n_blk;
    cram_block_compression_hdr *ch;

    if (nrg != 1) {
        hts_log_error("CRAM transcode supports only a single RG");
        return -2;
    }

    // Produce a new block holding the updated compression header,
    // with RG transcoded to a new value. (Single only supported.)
    o_blk = cram_read_block(in);
    old_size = cram_block_size(o_blk);
    ch = cram_decode_compression_header(in, o_blk);
    if (cram_block_compression_hdr_set_rg(ch, new_rg) != 0)
        return -1;
    if (cram_block_compression_hdr_decoder2encoder(in, ch) != 0)
        return -1;
    n_blk = cram_encode_compression_header(in, c, ch, in->embed_ref);
    cram_free_compression_header(ch);

    /*
     * Warning: this has internal knowledge of the cram compression
     * header format.
     *
     * The decoder doesn't set c->tags_used, so the encoder puts a two
     * byte blank segment.  This means n_blk is too short.  We skip
     * through the decoded old block (o_blk) and copy from there.
     */
    char *cp = cram_block_get_data(o_blk);
    char *op = cp;
    char *endp = cp + cram_block_get_uncomp_size(o_blk);
    //fprintf(stderr, "sz = %d\n", (int)(endp-cp));
    int32_t i32, err = 0;

    i32 = in->vv.varint_get32(&cp, endp, &err);
    cp += i32;
    i32 = in->vv.varint_get32(&cp, endp, &err);
    cp += i32;
    op = cp;
    i32 = in->vv.varint_get32(&cp, endp, &err);
    i32 += (cp-op);
    if (err)
        return -2;

    //fprintf(stderr, "remaining %d bytes\n", i32);
    cram_block_set_size(n_blk, cram_block_get_size(n_blk)-2);
    cram_block_append(n_blk, op, i32);
    cram_block_update_size(n_blk);

    new_size = cram_block_size(n_blk);

    //fprintf(stderr, "size %d -> %d\n", old_size, new_size);

    // Now we've constructedthe updated compression header,
    // amend the container too (it may have changed size).
    int32_t *landmarks, num_landmarks;
    landmarks = cram_container_get_landmarks(c, &num_landmarks);

    if (old_size != new_size) {
        int diff = new_size - old_size, j;

        for (j = 0; j < num_landmarks; j++)
            landmarks[j] += diff;
        //cram_container_set_landmarks(c, num_landmarks, landmarks);
        cram_container_set_length(c, cram_container_get_length(c) + diff);
    }

    // Finally write it all out; container, compression header,
    // and then all the remaining slice blocks.
    if (cram_write_container(out, c) != 0)
        return -2;

    cram_write_block(out, n_blk);
    cram_free_block(o_blk);
    cram_free_block(n_blk);

    // Container num_blocks can be invalid, due to a bug.
    // Instead we iterate in slice context instead.
    return cram_copy_slice(in, out, num_landmarks);
}


/*!
 * Returns the refs_t structure used by a cram file handle.
 *
 * This may be used in conjunction with option CRAM_OPT_SHARED_REF to
 * share reference memory between multiple file handles.
 *
 * @return
 * Returns NULL if none exists or the file handle is not a CRAM file.
 */
refs_t *cram_get_refs(htsFile *fd) {
    return fd->format.format == cram
        ? fd->fp.cram->refs
        : NULL;
}
