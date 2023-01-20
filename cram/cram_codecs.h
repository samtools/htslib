/*
Copyright (c) 2012-2015, 2018, 2020, 2023 Genome Research Ltd.
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

#ifndef CRAM_CODECS_H
#define CRAM_CODECS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct cram_codec;

/*
 * Slow but simple huffman decoder to start with.
 * Read a bit at a time, keeping track of {length, value}
 * eg. 1 1 0 1 => {1,1},  {2,3}, {3,6}, {4,13}
 *
 * Keep track of this through the huffman code table.
 * For fast scanning we have an index of where the first code of length X
 * appears.
 */
typedef struct {
    int64_t symbol;
    int32_t p; // next code start value, minus index to codes[]
    int32_t code;
    int32_t len;
} cram_huffman_code;

typedef struct {
    int ncodes;
    cram_huffman_code *codes;
    int option;
} cram_huffman_decoder;

#define MAX_HUFF 128
typedef struct {
    cram_huffman_code *codes;
    int nvals;
    int val2code[MAX_HUFF+1]; // value to code lookup for small values
    int option;
} cram_huffman_encoder;

typedef struct {
    int32_t offset;
    int32_t nbits;
} cram_beta_decoder;

// A PACK transform, packing multiple values into a single byte
typedef struct {
    int32_t nbits;
    enum cram_encoding sub_encoding;
    void *sub_codec_dat;
    struct cram_codec *sub_codec;
    int nval;  // number of items in maps
    uint32_t rmap[256]; // 0,1,2,3 -> P,A,C,K
    int map[256];       // P,A,C,K -> 0,1,2,3 // NB: max input is uint8_tb? Or use hash?
} cram_xpack_decoder;
typedef cram_xpack_decoder cram_xpack_encoder;

// Transforms symbols X,Y,Z to bytes 0,1,2.
typedef struct {
    enum cram_encoding len_encoding;
    enum cram_encoding lit_encoding;
    void *len_dat;
    void *lit_dat;
    struct cram_codec *len_codec;
    struct cram_codec *lit_codec;
    int cur_len;
    int cur_lit;
    int rep_score[256];
    char *to_flush;
    size_t to_flush_size;
} cram_xrle_decoder;
typedef cram_xrle_decoder cram_xrle_encoder;

// DELTA + zigzag + varint encoding
typedef struct {
    // FIXME: define endian here too.  Require little endian?
    int64_t last;
    uint8_t word_size; // 1, 2, 4, 8
    //uint8_t sign;      // true if input data is already signed
    enum cram_encoding sub_encoding;
    void *sub_codec_dat;
    struct cram_codec *sub_codec;
} cram_xdelta_decoder;
typedef cram_xdelta_decoder cram_xdelta_encoder;

typedef struct {
    int32_t offset;
} cram_gamma_decoder;

typedef struct {
    int32_t offset;
    int32_t k;
} cram_subexp_decoder;

typedef struct {
    int32_t content_id;
    enum cram_external_type type;
} cram_external_decoder;

typedef struct {
    int32_t content_id;
    int64_t offset;
    enum cram_external_type type;
} cram_varint_decoder;

typedef struct {
    struct cram_codec *len_codec;
    struct cram_codec *val_codec;
} cram_byte_array_len_decoder;

typedef struct {
    unsigned char stop;
    int32_t content_id;
} cram_byte_array_stop_decoder;

typedef struct {
    enum cram_encoding len_encoding;
    enum cram_encoding val_encoding;
    void *len_dat;
    void *val_dat;
    struct cram_codec *len_codec;
    struct cram_codec *val_codec;
} cram_byte_array_len_encoder;

typedef struct {
    int64_t val;
} cram_const_codec;

/*
 * A generic codec structure.
 */
struct cram_codec {
    enum cram_encoding codec;
    cram_block *out;
    varint_vec *vv;
    int codec_id;
    void (*free)(struct cram_codec *codec);
    int (*decode)(cram_slice *slice, struct cram_codec *codec,
                  cram_block *in, char *out, int *out_size);
    int (*encode)(cram_slice *slice, struct cram_codec *codec,
                  char *in, int in_size);
    int (*store)(struct cram_codec *codec, cram_block *b, char *prefix,
                 int version);
    int (*size)(cram_slice *slice, struct cram_codec *codec);
    int (*flush)(struct cram_codec *codec);
    cram_block *(*get_block)(cram_slice *slice, struct cram_codec *codec);
    int (*describe)(struct cram_codec *codec, kstring_t *ks);

    union {
        cram_huffman_decoder         huffman;
        cram_external_decoder        external;
        cram_beta_decoder            beta;
        cram_gamma_decoder           gamma;
        cram_subexp_decoder          subexp;
        cram_byte_array_len_decoder  byte_array_len;
        cram_byte_array_stop_decoder byte_array_stop;
        cram_xpack_decoder           xpack;
        cram_xrle_decoder            xrle;
        cram_xdelta_decoder          xdelta;
        cram_const_codec             xconst;
        cram_varint_decoder          varint;

        cram_huffman_encoder         e_huffman;
        cram_external_decoder        e_external;
        cram_byte_array_stop_decoder e_byte_array_stop;
        cram_byte_array_len_encoder  e_byte_array_len;
        cram_beta_decoder            e_beta;
        cram_xpack_decoder           e_xpack;
        cram_xrle_decoder            e_xrle;
        cram_xdelta_decoder          e_xdelta;
        cram_const_codec             e_xconst;
        cram_varint_decoder          e_varint;
    } u;
};

const char *cram_encoding2str(enum cram_encoding t);

cram_codec *cram_decoder_init(cram_block_compression_hdr *hdr,
                              enum cram_encoding codec, char *data, int size,
                              enum cram_external_type option,
                              int version, varint_vec *vv);
cram_codec *cram_encoder_init(enum cram_encoding codec, cram_stats *st,
                              enum cram_external_type option, void *dat,
                              int version, varint_vec *vv);

//int cram_decode(void *codes, char *in, int in_size, char *out, int *out_size);
//void cram_decoder_free(void *codes);

//#define GET_BIT_MSB(b,v) (void)(v<<=1, v|=(b->data[b->byte] >> b->bit)&1, (--b->bit == -1) && (b->bit = 7, b->byte++))

#define GET_BIT_MSB(b,v) (void)(v<<=1, v|=(b->data[b->byte] >> b->bit)&1, b->byte += (--b->bit<0), b->bit&=7)

/*
 * Check that enough bits are left in a block to satisy a bit-based decoder.
 * Return  0 if there are enough
 *         1 if not.
 */

static inline int cram_not_enough_bits(cram_block *blk, int nbits) {
    if (nbits < 0 ||
        (blk->byte >= blk->uncomp_size && nbits > 0) ||
        (blk->uncomp_size - blk->byte <= INT32_MAX / 8 + 1 &&
         (blk->uncomp_size - blk->byte) * 8 + blk->bit - 7 < nbits)) {
        return 1;
    }
    return 0;
}

/*
 * Returns the content_id used by this codec, also in id2 if byte_array_len.
 * Returns -1 for the CORE block and -2 for unneeded.
 * id2 is only filled out for BYTE_ARRAY_LEN which uses 2 codecs.
 */
int cram_codec_to_id(cram_codec *c, int *id2);

/*
 * cram_codec structures are specialised for decoding or encoding.
 * Unfortunately this makes turning a decoder into an encoder (such as
 * when transcoding files) problematic.
 *
 * This function converts a cram decoder codec into an encoder version
 * in-place (ie it modifiers the codec itself).
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int cram_codec_decoder2encoder(cram_fd *fd, cram_codec *c);

#ifdef __cplusplus
}
#endif

#endif /* CRAM_CODECS_H */
