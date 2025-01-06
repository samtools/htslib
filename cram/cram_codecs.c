/*
Copyright (c) 2012-2021,2023 Genome Research Ltd.
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
 * FIXME: add checking of cram_external_type to return NULL on unsupported
 * {codec,type} tuples.
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <errno.h>
#include <stddef.h>

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
#include "../fuzz_settings.h"
#endif

#include "../htslib/hts_endian.h"

#if defined(HAVE_EXTERNAL_LIBHTSCODECS)
#include <htscodecs/varint.h>
#include <htscodecs/pack.h>
#include <htscodecs/rle.h>
#else
#include "../htscodecs/htscodecs/varint.h"
#include "../htscodecs/htscodecs/pack.h"
#include "../htscodecs/htscodecs/rle.h"
#endif

#include "cram.h"

/*
 * ---------------------------------------------------------------------------
 * Block bit-level I/O functions.
 * All defined static here to promote easy inlining by the compiler.
 */

#if 0
/* Get a single bit, MSB first */
static signed int get_bit_MSB(cram_block *block) {
    unsigned int val;

    if (block->byte > block->alloc)
        return -1;

    val = block->data[block->byte] >> block->bit;
    if (--block->bit == -1) {
        block->bit = 7;
        block->byte++;
        //printf("(%02X)", block->data[block->byte]);
    }

    //printf("-B%d-", val&1);

    return val & 1;
}
#endif

/*
 * Count number of successive 0 and 1 bits
 */
static int get_one_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
        b = block->data[block->byte] >> block->bit;
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            if (block->byte == block->uncomp_size && (b&1))
                return -1;
        }
        n++;
    } while (b&1);

    return n-1;
}

static int get_zero_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
        b = block->data[block->byte] >> block->bit;
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            if (block->byte == block->uncomp_size && !(b&1))
                return -1;
        }
        n++;
    } while (!(b&1));

    return n-1;
}

#if 0
/* Stores a single bit */
static void store_bit_MSB(cram_block *block, unsigned int bit) {
    if (block->byte >= block->alloc) {
        block->alloc = block->alloc ? block->alloc*2 : 1024;
        block->data = realloc(block->data, block->alloc);
    }

    if (bit)
        block->data[block->byte] |= (1 << block->bit);

    if (--block->bit == -1) {
        block->bit = 7;
        block->byte++;
        block->data[block->byte] = 0;
    }
}
#endif

#if 0
/* Rounds to the next whole byte boundary first */
static void store_bytes_MSB(cram_block *block, char *bytes, int len) {
    if (block->bit != 7) {
        block->bit = 7;
        block->byte++;
    }

    while (block->byte + len >= block->alloc) {
        block->alloc = block->alloc ? block->alloc*2 : 1024;
        block->data = realloc(block->data, block->alloc);
    }

    memcpy(&block->data[block->byte], bytes, len);
    block->byte += len;
}
#endif

/* Local optimised copy for inlining */
static inline int64_t get_bits_MSB(cram_block *block, int nbits) {
    uint64_t val = 0;
    int i;

#if 0
    // Fits within the current byte */
    if (nbits <= block->bit+1) {
        val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
        if ((block->bit -= nbits) == -1) {
            block->bit = 7;
            block->byte++;
        }
        return val;
    }

    // partial first byte
    val = block->data[block->byte] & ((1<<(block->bit+1))-1);
    nbits -= block->bit+1;
    block->bit = 7;
    block->byte++;

    // whole middle bytes
    while (nbits >= 8) {
        val = (val << 8) | block->data[block->byte++];
        nbits -= 8;
    }

    val <<= nbits;
    val |= (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
    block->bit -= nbits;
    return val;
#endif

#if 0
    /* Inefficient implementation! */
    //printf("{");
    for (i = 0; i < nbits; i++)
        //val = (val << 1) | get_bit_MSB(block);
        GET_BIT_MSB(block, val);
#endif

#if 1
    /* Combination of 1st two methods */
    if (nbits <= block->bit+1) {
        val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
        if ((block->bit -= nbits) == -1) {
            block->bit = 7;
            block->byte++;
        }
        return val;
    }

    switch(nbits) {
//  case 15: GET_BIT_MSB(block, val); // fall through
//  case 14: GET_BIT_MSB(block, val); // fall through
//  case 13: GET_BIT_MSB(block, val); // fall through
//  case 12: GET_BIT_MSB(block, val); // fall through
//  case 11: GET_BIT_MSB(block, val); // fall through
//  case 10: GET_BIT_MSB(block, val); // fall through
//  case  9: GET_BIT_MSB(block, val); // fall through
    case  8: GET_BIT_MSB(block, val); // fall through
    case  7: GET_BIT_MSB(block, val); // fall through
    case  6: GET_BIT_MSB(block, val); // fall through
    case  5: GET_BIT_MSB(block, val); // fall through
    case  4: GET_BIT_MSB(block, val); // fall through
    case  3: GET_BIT_MSB(block, val); // fall through
    case  2: GET_BIT_MSB(block, val); // fall through
    case  1: GET_BIT_MSB(block, val);
        break;

    default:
        for (i = 0; i < nbits; i++)
            //val = (val << 1) | get_bit_MSB(block);
            GET_BIT_MSB(block, val);
    }
#endif

    //printf("=0x%x}", val);

    return val;
}

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static int store_bits_MSB(cram_block *block, uint64_t val, int nbits) {
    //fprintf(stderr, " store_bits: %02x %d\n", val, nbits);

    /*
     * Use slow mode until we tweak the huffman generator to never generate
     * codes longer than 24-bits.
     */
    unsigned int mask;

    if (block->byte+8 >= block->alloc) {
        if (block->byte) {
            block->alloc *= 2;
            block->data = realloc(block->data, block->alloc + 8);
            if (!block->data)
                return -1;
        } else {
            block->alloc = 1024;
            block->data = realloc(block->data, block->alloc + 8);
            if (!block->data)
                return -1;
            block->data[0] = 0; // initialise first byte of buffer
        }
    }

    /* fits in current bit-field */
    if (nbits <= block->bit+1) {
        block->data[block->byte] |= (val << (block->bit+1-nbits));
        if ((block->bit-=nbits) == -1) {
            block->bit = 7;
            block->byte++;
            block->data[block->byte] = 0;
        }
        return 0;
    }

    block->data[block->byte] |= (val >> (nbits -= block->bit+1));
    block->bit = 7;
    block->byte++;
    block->data[block->byte] = 0;

    mask = 1<<(nbits-1);
    do {
        if (val & mask)
            block->data[block->byte] |= (1 << block->bit);
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            block->data[block->byte] = 0;
        }
        mask >>= 1;
    } while(--nbits);

    return 0;
}

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left.This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
static char *cram_extract_block(cram_block *b, int size) {
    char *cp = (char *)b->data + b->idx;
    b->idx += size;
    if (b->idx > b->uncomp_size)
        return NULL;

    return cp;
}

/*
 * ---------------------------------------------------------------------------
 * EXTERNAL
 *
 * In CRAM 3.0 and earlier, E_EXTERNAL use the data type to determine the
 * size of the object being returned.  This type is hard coded in the
 * spec document (changing from uint32 to uint64 requires a spec change)
 * and there is no data format introspection so implementations have
 * to determine which size to use based on version numbers.   It also
 * doesn't support signed data.
 *
 * With CRAM 4.0 onwards the size and sign of the data is no longer stated
 * explicitly in the specification.  Instead EXTERNAL is replaced by three
 * new encodings, for bytes and signed / unsigned integers which used a
 * variable sized encoding.
 *
 * For simplicity we use the same encode and decode functions for
 * bytes (CRAM4) and external (CRAM3). Given we already had code to
 * replace codec + type into a function pointer it makes little
 * difference how we ended up at that function.  However we disallow
 * this codec to operate on integer data for CRAM4 onwards.
 */
int cram_external_decode_int(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->u.external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int32_t *)out = c->vv->varint_get32(&cp, (char *)b->data + b->uncomp_size, &err);
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_external_decode_long(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->u.external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int64_t *)out = c->vv->varint_get64(&cp, (char *)b->data + b->uncomp_size, &err);
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_external_decode_char(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out,
                              int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->u.external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
        return -1;

    if (out)
        memcpy(out, cp, *out_size);
    return 0;
}

static int cram_external_decode_block(cram_slice *slice, cram_codec *c,
                                      cram_block *in, char *out_,
                                      int *out_size) {
    char *cp;
    cram_block *out = (cram_block *)out_;
    cram_block *b = NULL;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->u.external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
        return -1;

    BLOCK_APPEND(out, cp, *out_size);
    return 0;

 block_err:
    return -1;
}

void cram_external_decode_free(cram_codec *c) {
    if (c)
        free(c);
}


int cram_external_decode_size(cram_slice *slice, cram_codec *c) {
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->u.external.content_id);
    if (!b)
        return -1;

    return b->uncomp_size;
}

cram_block *cram_external_get_block(cram_slice *slice, cram_codec *c) {
    return cram_get_block_by_id(slice, c->u.external.content_id);
}

int cram_external_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "EXTERNAL(id=%d)",
                    c->u.external.content_id) < 0 ? -1 : 0;
}

cram_codec *cram_external_decode_init(cram_block_compression_hdr *hdr,
                                      char *data, int size,
                                      enum cram_encoding codec,
                                      enum cram_external_type option,
                                      int version, varint_vec *vv) {
    cram_codec *c = NULL;
    char *cp = data;

    if (size < 1)
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_EXTERNAL;
    if (CRAM_MAJOR_VERS(version) >= 4) {
        // Version 4 does not permit integer data to be encoded as a
        // series of bytes.  This is used purely for bytes, either
        // singular or declared as arrays
        switch (codec) {
        case E_EXTERNAL:
            if (option == E_BYTE_ARRAY_BLOCK)
                c->decode = cram_external_decode_block;
            else if (option == E_BYTE || option == E_BYTE_ARRAY)
                c->decode = cram_external_decode_char;
            else
                goto malformed;
            break;
        default:
            goto malformed;
        }
    } else {
        // CRAM 3 and earlier encodes integers as EXTERNAL.  We need
        // use the option field to indicate the input data format so
        // we know which serialisation format to use.
        if (option == E_INT)
            c->decode = cram_external_decode_int;
        else if (option == E_LONG)
            c->decode = cram_external_decode_long;
        else if (option == E_BYTE_ARRAY || option == E_BYTE)
            c->decode = cram_external_decode_char;
        else
            c->decode = cram_external_decode_block;
    }
    c->free   = cram_external_decode_free;
    c->size   = cram_external_decode_size;
    c->get_block = cram_external_get_block;
    c->describe = cram_external_describe;

    c->u.external.content_id = vv->varint_get32(&cp, data+size, NULL);

    if (cp - data != size)
        goto malformed;

    c->u.external.type = option;

    return c;

 malformed:
    hts_log_error("Malformed external header stream");
    free(c);
    return NULL;
}

int cram_external_encode_int(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    uint32_t *i32 = (uint32_t *)in;
    return c->vv->varint_put32_blk(c->out, *i32) >= 0 ? 0 : -1;
}

int cram_external_encode_sint(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    int32_t *i32 = (int32_t *)in;
    return c->vv->varint_put32s_blk(c->out, *i32) >= 0 ? 0 : -1;
}

int cram_external_encode_long(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    uint64_t *i64 = (uint64_t *)in;
    return c->vv->varint_put64_blk(c->out, *i64) >= 0 ? 0 : -1;
}

int cram_external_encode_slong(cram_slice *slice, cram_codec *c,
                               char *in, int in_size) {
    int64_t *i64 = (int64_t *)in;
    return c->vv->varint_put64s_blk(c->out, *i64) >= 0 ? 0 : -1;
}

int cram_external_encode_char(cram_slice *slice, cram_codec *c,
                              char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    return 0;

 block_err:
    return -1;
}

void cram_external_encode_free(cram_codec *c) {
    if (!c)
        return;
    free(c);
}

int cram_external_encode_store(cram_codec *c, cram_block *b, char *prefix,
                               int version) {
    char tmp[99], *tp = tmp, *tpend = tmp+99;
    int len = 0, r = 0, n;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += c->vv->varint_put32(tp, tpend, c->u.e_external.content_id);
    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;
    len += (n = c->vv->varint_put32_blk(b, tp-tmp));   r |= n;
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    if (r > 0)
        return len;

 block_err:
    return -1;
}

cram_codec *cram_external_encode_init(cram_stats *st,
                                      enum cram_encoding codec,
                                      enum cram_external_type option,
                                      void *dat,
                                      int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_EXTERNAL;
    c->free = cram_external_encode_free;
    if (CRAM_MAJOR_VERS(version) >= 4) {
        // Version 4 does not permit integer data to be encoded as a
        // series of bytes.  This is used purely for bytes, either
        // singular or declared as arrays
        switch (codec) {
        case E_EXTERNAL:
            if (option != E_BYTE && option != E_BYTE_ARRAY)
                return NULL;
            c->encode = cram_external_encode_char;
            break;
        default:
            return NULL;
        }
    } else {
        // CRAM 3 and earlier encodes integers as EXTERNAL.  We need
        // use the option field to indicate the input data format so
        // we know which serialisation format to use.
        if (option == E_INT)
            c->encode = cram_external_encode_int;
        else if (option == E_LONG)
            c->encode = cram_external_encode_long;
        else if (option == E_BYTE_ARRAY || option == E_BYTE)
            c->encode = cram_external_encode_char;
        else
            abort();
    }
    c->store = cram_external_encode_store;
    c->flush = NULL;

    c->u.e_external.content_id = (size_t)dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * VARINT
 *
 * In CRAM 3.0 and earlier, E_EXTERNAL stored both integers in ITF8
 * format as well as bytes.  In CRAM 4 EXTERNAL is only for bytes and
 * byte arrays, with two dedicated encodings for integers:
 * VARINT_SIGNED and VARINT_UNSIGNED.  These also differ a little to
 * EXTERNAL with the addition of an offset field, meaning we can store
 * values in, say, the range -2 to 1 million without needing to use
 * a signed zig-zag transformation.
 */
int cram_varint_decode_int(cram_slice *slice, cram_codec *c,
                           cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the data block */
    b = cram_get_block_by_id(slice, c->u.varint.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int32_t *)out = c->vv->varint_get32(&cp,
                                          (char *)b->data + b->uncomp_size,
                                          &err) + c->u.varint.offset;
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_varint_decode_sint(cram_slice *slice, cram_codec *c,
                            cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the data block */
    b = cram_get_block_by_id(slice, c->u.varint.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int32_t *)out = c->vv->varint_get32s(&cp,
                                           (char *)b->data + b->uncomp_size,
                                           &err) + c->u.varint.offset;
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_varint_decode_long(cram_slice *slice, cram_codec *c,
                            cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the data block */
    b = cram_get_block_by_id(slice, c->u.varint.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int64_t *)out = c->vv->varint_get64(&cp,
                                          (char *)b->data + b->uncomp_size,
                                          &err) + c->u.varint.offset;
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_varint_decode_slong(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the data block */
    b = cram_get_block_by_id(slice, c->u.varint.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int64_t *)out = c->vv->varint_get64s(&cp,
                                           (char *)b->data + b->uncomp_size,
                                           &err) + c->u.varint.offset;
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

void cram_varint_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

int cram_varint_decode_size(cram_slice *slice, cram_codec *c) {
    cram_block *b;

    /* Find the data block */
    b = cram_get_block_by_id(slice, c->u.varint.content_id);
    if (!b)
        return -1;

    return b->uncomp_size;
}

cram_block *cram_varint_get_block(cram_slice *slice, cram_codec *c) {
    return cram_get_block_by_id(slice, c->u.varint.content_id);
}

int cram_varint_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "VARINT(id=%d,offset=%"PRId64",type=%d)",
                    c->u.varint.content_id,
                    c->u.varint.offset,
                    c->u.varint.type)
        < 0 ? -1 : 0;
}

cram_codec *cram_varint_decode_init(cram_block_compression_hdr *hdr,
                                    char *data, int size,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data, *cp_end = data+size;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = codec;

    // Function pointer choice is theoretically by codec type.
    // Given we have some vars as int32 and some as int64 we
    // use option too for sizing, although on disk format
    // does not change.
    switch(codec) {
    case E_VARINT_UNSIGNED:
        c->decode = (option == E_INT)
            ? cram_varint_decode_int
            : cram_varint_decode_long;
        break;
    case E_VARINT_SIGNED:
        c->decode = (option == E_INT)
            ? cram_varint_decode_sint
            : cram_varint_decode_slong;
        break;
    default:
        return NULL;
    }

    c->free   = cram_varint_decode_free;
    c->size   = cram_varint_decode_size;
    c->get_block = cram_varint_get_block;
    c->describe = cram_varint_describe;

    c->u.varint.content_id = vv->varint_get32 (&cp, cp_end, NULL);
    c->u.varint.offset     = vv->varint_get64s(&cp, cp_end, NULL);

    if (cp - data != size) {
        fprintf(stderr, "Malformed varint header stream\n");
        free(c);
        return NULL;
    }

    c->u.varint.type = option;

    return c;
}

int cram_varint_encode_int(cram_slice *slice, cram_codec *c,
                           char *in, int in_size) {
    uint32_t *i32 = (uint32_t *)in;
    return c->vv->varint_put32_blk(c->out, *i32 - c->u.varint.offset) >= 0
        ? 0 : -1;
}

int cram_varint_encode_sint(cram_slice *slice, cram_codec *c,
                            char *in, int in_size) {
    int32_t *i32 = (int32_t *)in;
    return c->vv->varint_put32s_blk(c->out, *i32 - c->u.varint.offset) >= 0
        ? 0 : -1;
}

int cram_varint_encode_long(cram_slice *slice, cram_codec *c,
                            char *in, int in_size) {
    uint64_t *i64 = (uint64_t *)in;
    return c->vv->varint_put64_blk(c->out, *i64 - c->u.varint.offset) >= 0
        ? 0 : -1;
}

int cram_varint_encode_slong(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    int64_t *i64 = (int64_t *)in;
    return c->vv->varint_put64s_blk(c->out, *i64 - c->u.varint.offset) >= 0
        ? 0 : -1;
}

void cram_varint_encode_free(cram_codec *c) {
    if (!c)
        return;
    free(c);
}

int cram_varint_encode_store(cram_codec *c, cram_block *b, char *prefix,
                             int version) {
    char tmp[99], *tp = tmp;
    int len = 0;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += c->vv->varint_put32 (tp, NULL, c->u.e_varint.content_id);
    tp += c->vv->varint_put64s(tp, NULL, c->u.e_varint.offset);
    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    return len;

 block_err:
    return -1;
}

cram_codec *cram_varint_encode_init(cram_stats *st,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    void *dat,
                                    int version, varint_vec *vv) {
    cram_codec *c;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->u.e_varint.offset = 0;
    if (st) {
        // Marginal difference so far! Not worth the hassle?
        if (st->min_val < 0 && st->min_val >= -127
            && st->max_val / -st->min_val > 100) {
            c->u.e_varint.offset = -st->min_val;
            codec = E_VARINT_UNSIGNED;
        } else if (st->min_val > 0) {
            c->u.e_varint.offset = -st->min_val;
        }
    }

    c->codec = codec;
    c->free = cram_varint_encode_free;

    // Function pointer choice is theoretically by codec type.
    // Given we have some vars as int32 and some as int64 we
    // use option too for sizing, although on disk format
    // does not change.
    switch (codec) {
    case E_VARINT_UNSIGNED:
        c->encode = (option == E_INT)
            ? cram_varint_encode_int
            : cram_varint_encode_long;
        break;
    case E_VARINT_SIGNED:
        c->encode = (option == E_INT)
            ? cram_varint_encode_sint
            : cram_varint_encode_slong;
        break;
    default:
        return NULL;
    }
    c->store = cram_varint_encode_store;
    c->flush = NULL;

    c->u.e_varint.content_id = (size_t)dat;

    return c;
}
/*
 * ---------------------------------------------------------------------------
 * CONST_BYTE and CONST_INT
 */
int cram_const_decode_byte(cram_slice *slice, cram_codec *c,
                           cram_block *in, char *out, int *out_size) {
    int i, n;

    for (i = 0, n = *out_size; i < n; i++)
        out[i] = c->u.xconst.val;

    return 0;
}

int cram_const_decode_int(cram_slice *slice, cram_codec *c,
                          cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++)
        out_i[i] = c->u.xconst.val;

    return 0;
}

int cram_const_decode_long(cram_slice *slice, cram_codec *c,
                           cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++)
        out_i[i] = c->u.xconst.val;

    return 0;
}

void cram_const_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

int cram_const_decode_size(cram_slice *slice, cram_codec *c) {
    return 0;
}

int cram_const_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "CONST(val=%"PRId64")",
                    c->u.xconst.val) < 0 ? -1 : 0;
}

cram_codec *cram_const_decode_init(cram_block_compression_hdr *hdr,
                                   char *data, int size,
                                   enum cram_encoding codec,
                                   enum cram_external_type option,
                                   int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = codec;
    if (codec == E_CONST_BYTE)
        c->decode = cram_const_decode_byte;
    else if (option == E_INT)
        c->decode = cram_const_decode_int;
    else
        c->decode = cram_const_decode_long;
    c->free   = cram_const_decode_free;
    c->size   = cram_const_decode_size;
    c->get_block = NULL;
    c->describe = cram_const_describe;

    c->u.xconst.val = vv->varint_get64s(&cp, data+size, NULL);

    if (cp - data != size) {
        fprintf(stderr, "Malformed const header stream\n");
        free(c);
        return NULL;
    }

    return c;
}

int cram_const_encode(cram_slice *slice, cram_codec *c,
                      char *in, int in_size) {
    return 0;
}

int cram_const_encode_store(cram_codec *c, cram_block *b, char *prefix,
                            int version) {
    char tmp[99], *tp = tmp;
    int len = 0;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += c->vv->varint_put64s(tp, NULL, c->u.xconst.val);
    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    return len;

 block_err:
    return -1;
}

cram_codec *cram_const_encode_init(cram_stats *st,
                                   enum cram_encoding codec,
                                   enum cram_external_type option,
                                   void *dat,
                                   int version, varint_vec *vv) {
    cram_codec *c;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec = codec;
    c->free = cram_const_decode_free; // as as decode
    c->encode = cram_const_encode; // a nop
    c->store = cram_const_encode_store;
    c->flush = NULL;
    c->u.e_xconst.val = st->min_val;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
int cram_beta_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n = *out_size;

    if (c->u.beta.nbits) {
        if (cram_not_enough_bits(in, c->u.beta.nbits * n))
            return -1;

        for (i = 0; i < n; i++)
            out_i[i] = get_bits_MSB(in, c->u.beta.nbits) - c->u.beta.offset;
    } else {
        for (i = 0; i < n; i++)
            out_i[i] = -c->u.beta.offset;
    }

    return 0;
}

int cram_beta_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n = *out_size;

    if (c->u.beta.nbits) {
        if (cram_not_enough_bits(in, c->u.beta.nbits * n))
            return -1;

        for (i = 0; i < n; i++)
            out_i[i] = get_bits_MSB(in, c->u.beta.nbits) - c->u.beta.offset;
    } else {
        for (i = 0; i < n; i++)
            out_i[i] = -c->u.beta.offset;
    }

    return 0;
}

int cram_beta_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n = *out_size;


    if (c->u.beta.nbits) {
        if (cram_not_enough_bits(in, c->u.beta.nbits * n))
            return -1;

        if (out)
            for (i = 0; i < n; i++)
                out[i] = get_bits_MSB(in, c->u.beta.nbits) - c->u.beta.offset;
        else
            for (i = 0; i < n; i++)
                get_bits_MSB(in, c->u.beta.nbits);
    } else {
        if (out)
            for (i = 0; i < n; i++)
                out[i] = -c->u.beta.offset;
    }

    return 0;
}

void cram_beta_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

int cram_beta_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "BETA(offset=%d, nbits=%d)",
                    c->u.beta.offset, c->u.beta.nbits)
        < 0 ? -1 : 0;
}

cram_codec *cram_beta_decode_init(cram_block_compression_hdr *hdr,
                                  char *data, int size,
                                  enum cram_encoding codec,
                                  enum cram_external_type option,
                                  int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BETA;
    if (option == E_INT || option == E_SINT)
        c->decode = cram_beta_decode_int;
    else if (option == E_LONG || option == E_SLONG)
        c->decode = cram_beta_decode_long;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_beta_decode_char;
    else {
        hts_log_error("BYTE_ARRAYs not supported by this codec");
        free(c);
        return NULL;
    }
    c->free   = cram_beta_decode_free;
    c->describe = cram_beta_describe;

    c->u.beta.nbits = -1;
    c->u.beta.offset = vv->varint_get32(&cp, data + size, NULL);
    if (cp < data + size) // Ensure test below works
        c->u.beta.nbits  = vv->varint_get32(&cp, data + size, NULL);

    if (cp - data != size
        || c->u.beta.nbits < 0 || c->u.beta.nbits > 8 * sizeof(int)) {
        hts_log_error("Malformed beta header stream");
        free(c);
        return NULL;
    }

    return c;
}

int cram_beta_encode_store(cram_codec *c, cram_block *b,
                           char *prefix, int version) {
    int len = 0, r = 0, n;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;
    // codec length
    len += (n = c->vv->varint_put32_blk(b, c->vv->varint_size(c->u.e_beta.offset)
                                         + c->vv->varint_size(c->u.e_beta.nbits)));
    r |= n;
    len += (n = c->vv->varint_put32_blk(b, c->u.e_beta.offset)); r |= n;
    len += (n = c->vv->varint_put32_blk(b, c->u.e_beta.nbits));  r |= n;

    if (r > 0) return len;

 block_err:
    return -1;
}

int cram_beta_encode_long(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    int64_t *syms = (int64_t *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, syms[i] + c->u.e_beta.offset,
                            c->u.e_beta.nbits);

    return r;
}

int cram_beta_encode_int(cram_slice *slice, cram_codec *c,
                         char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, syms[i] + c->u.e_beta.offset,
                            c->u.e_beta.nbits);

    return r;
}

int cram_beta_encode_char(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, syms[i] + c->u.e_beta.offset,
                            c->u.e_beta.nbits);

    return r;
}

void cram_beta_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_beta_encode_init(cram_stats *st,
                                  enum cram_encoding codec,
                                  enum cram_external_type option,
                                  void *dat,
                                  int version, varint_vec *vv) {
    cram_codec *c;
    hts_pos_t min_val, max_val;
    int len = 0;
    int64_t range;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec  = E_BETA;
    c->free   = cram_beta_encode_free;
    if (option == E_INT || option == E_SINT)
        c->encode = cram_beta_encode_int;
    else if (option == E_LONG || option == E_SLONG)
        c->encode = cram_beta_encode_long;
    else
        c->encode = cram_beta_encode_char;
    c->store  = cram_beta_encode_store;
    c->flush = NULL;

    if (dat) {
        min_val = ((hts_pos_t *)dat)[0];
        max_val = ((hts_pos_t *)dat)[1];
    } else {
        min_val = INT_MAX;
        max_val = INT_MIN;
        int i;
        for (i = 0; i < MAX_STAT_VAL; i++) {
            if (!st->freqs[i])
                continue;
            if (min_val > i)
                min_val = i;
            max_val = i;
        }
        if (st->h) {
            khint_t k;

            for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
                if (!kh_exist(st->h, k))
                    continue;

                i = kh_key(st->h, k);
                if (min_val > i)
                    min_val = i;
                if (max_val < i)
                    max_val = i;
            }
        }
    }

    if (max_val < min_val)
        goto err;

    range = (int64_t) max_val - min_val;
    switch (option) {
    case E_SINT:
        if (min_val < INT_MIN || range > INT_MAX)
            goto err;
        break;

    case E_INT:
        if (max_val > UINT_MAX || range > UINT_MAX)
            goto err;
        break;

    default:
        break;
    }

    c->u.e_beta.offset = -min_val;
    while (range) {
        len++;
        range >>= 1;
    }
    c->u.e_beta.nbits = len;

    return c;

 err:
    free(c);
    return NULL;
}

/*
 * ---------------------------------------------------------------------------
 * XPACK: Packing multiple values into a single byte.  A fast transform that
 * reduces time taken by entropy encoder and may also improve compression.
 *
 * This also has the additional requirement that the data series is not
 * interleaved with another, permitting efficient encoding and decoding
 * of all elements enmasse instead of needing to only extract the bits
 * necessary per item.
 */
int cram_xpack_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n = *out_size;

    if (c->u.xpack.nbits) {
        for (i = 0; i < n; i++)
            out_i[i] = c->u.xpack.rmap[get_bits_MSB(in, c->u.xpack.nbits)];
    } else {
        for (i = 0; i < n; i++)
            out_i[i] = c->u.xpack.rmap[0];
    }

    return 0;
}

int cram_xpack_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n = *out_size;

    if (c->u.xpack.nbits) {
        if (cram_not_enough_bits(in, c->u.xpack.nbits * n))
            return -1;

        for (i = 0; i < n; i++)
            out_i[i] = c->u.xpack.rmap[get_bits_MSB(in, c->u.xpack.nbits)];
    } else {
        for (i = 0; i < n; i++)
            out_i[i] = c->u.xpack.rmap[0];
    }

    return 0;
}

static int cram_xpack_decode_expand_char(cram_slice *slice, cram_codec *c) {
    cram_block *b = slice->block_by_id[512 + c->codec_id];
    if (b)
        return 0;

    // get sub-codec data.
    cram_block *sub_b = c->u.xpack.sub_codec->get_block(slice, c->u.xpack.sub_codec);
    if (!sub_b)
        return -1;

    // Allocate local block to expand into
    b = slice->block_by_id[512 + c->codec_id] = cram_new_block(0, 0);
    if (!b)
        return -1;
    int n = sub_b->uncomp_size * 8/c->u.xpack.nbits;
    BLOCK_GROW(b, n);
    b->uncomp_size = n;

    uint8_t p[256];
    int z;
    for (z = 0; z < 256; z++)
        p[z] = c->u.xpack.rmap[z];
    hts_unpack(sub_b->data, sub_b->uncomp_size, b->data, b->uncomp_size,
               8 / c->u.xpack.nbits, p);

    return 0;

 block_err:
    return -1;
}

int cram_xpack_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // FIXME: we need to ban data-series interleaving in the spec for this to work.

    // Remember this may be called when threaded and multi-slice per container.
    // Hence one cram_codec instance, multiple slices, multiple blocks.
    // We therefore have to cache appropriate block info in slice and not codec.
    //    b = cram_get_block_by_id(slice, c->external.content_id);
    if (c->u.xpack.nval > 1) {
        cram_xpack_decode_expand_char(slice, c);
        cram_block *b = slice->block_by_id[512 + c->codec_id];
        if (!b)
            return -1;

        if (out)
            memcpy(out, b->data + b->byte, *out_size);
        b->byte += *out_size;
    } else {
        memset(out, c->u.xpack.rmap[0], *out_size);
    }

    return 0;
}

void cram_xpack_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.xpack.sub_codec)
        c->u.xpack.sub_codec->free(c->u.xpack.sub_codec);

    //free(slice->block_by_id[512 + c->codec_id]);
    //slice->block_by_id[512 + c->codec_id] = 0;

    free(c);
}

int cram_xpack_decode_size(cram_slice *slice, cram_codec *c) {
    cram_xpack_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id]->uncomp_size;
}

cram_block *cram_xpack_get_block(cram_slice *slice, cram_codec *c) {
    cram_xpack_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id];
}

cram_codec *cram_xpack_decode_init(cram_block_compression_hdr *hdr,
                                   char *data, int size,
                                   enum cram_encoding codec,
                                   enum cram_external_type option,
                                   int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;
    char *endp = data+size;

    if (!(c = calloc(1, sizeof(*c))))
        return NULL;

    c->codec  = E_XPACK;
    if (option == E_LONG)
        c->decode = cram_xpack_decode_long;
    else if (option == E_INT)
        c->decode = cram_xpack_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_xpack_decode_char;
    else {
        fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
        goto malformed;
    }
    c->free = cram_xpack_decode_free;
    c->size = cram_xpack_decode_size;
    c->get_block = cram_xpack_get_block;
    c->describe = NULL;

    c->u.xpack.nbits = vv->varint_get32(&cp, endp, NULL);
    c->u.xpack.nval  = vv->varint_get32(&cp, endp, NULL);
    if (c->u.xpack.nbits >= 8  || c->u.xpack.nbits < 0 ||
        c->u.xpack.nval  > 256 || c->u.xpack.nval < 0)
        goto malformed;
    int i;
    for (i = 0; i < c->u.xpack.nval; i++) {
        uint32_t v = vv->varint_get32(&cp, endp, NULL);
        if (v >= 256)
            goto malformed;
        c->u.xpack.rmap[i] = v; // reverse map: e.g 0-3 to P,A,C,K
    }

    int encoding = vv->varint_get32(&cp, endp, NULL);
    int sub_size = vv->varint_get32(&cp, endp, NULL);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.xpack.sub_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
                                             option, version, vv);
    if (c->u.xpack.sub_codec == NULL)
        goto malformed;
    cp += sub_size;

    if (cp - data != size
        || c->u.xpack.nbits < 0 || c->u.xpack.nbits > 8 * sizeof(int64_t)) {
    malformed:
        fprintf(stderr, "Malformed xpack header stream\n");
        cram_xpack_decode_free(c);
        return NULL;
    }

    return c;
}

int cram_xpack_encode_flush(cram_codec *c) {
    // Pack the buffered up data
    int meta_len;
    uint64_t out_len;
    uint8_t out_meta[1024];
    uint8_t *out = hts_pack(BLOCK_DATA(c->out), BLOCK_SIZE(c->out),
                            out_meta, &meta_len, &out_len);

    // We now need to pass this through the next layer of transform
    if (c->u.e_xpack.sub_codec->encode(NULL, // also indicates flush incoming
                                     c->u.e_xpack.sub_codec,
                                     (char *)out, out_len))
        return -1;

    int r = 0;
    if (c->u.e_xpack.sub_codec->flush)
        r = c->u.e_xpack.sub_codec->flush(c->u.e_xpack.sub_codec);

    free(out);
    return r;
}

int cram_xpack_encode_store(cram_codec *c, cram_block *b,
                            char *prefix, int version) {
    int len = 0, r = 0, n;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    // Store sub-codec
    cram_codec *tc = c->u.e_xpack.sub_codec;
    cram_block *tb = cram_new_block(0, 0);
    if (!tb)
        return -1;
    int len2 = tc->store(tc, tb, NULL, version);

    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;

    // codec length
    int len1 = 0, i;
    for (i = 0; i < c->u.e_xpack.nval; i++)
        len1 += (n = c->vv->varint_size(c->u.e_xpack.rmap[i])), r |= n;
    len += (n = c->vv->varint_put32_blk(b, c->vv->varint_size(c->u.e_xpack.nbits)
                                        +  c->vv->varint_size(c->u.e_xpack.nval)
                                        + len1 + len2)); r |= n;

    // The map and sub-codec
    len += (n = c->vv->varint_put32_blk(b, c->u.e_xpack.nbits)); r |= n;
    len += (n = c->vv->varint_put32_blk(b, c->u.e_xpack.nval));  r |= n;
    for (i = 0; i < c->u.e_xpack.nval; i++)
        len += (n = c->vv->varint_put32_blk(b, c->u.e_xpack.rmap[i])), r |= n;

    BLOCK_APPEND(b, BLOCK_DATA(tb), BLOCK_SIZE(tb));

    cram_free_block(tb);

    return r > 0 ? len + len2 : -1;

 block_err:
    return -1;
}

// Same as cram_beta_encode_long
int cram_xpack_encode_long(cram_slice *slice, cram_codec *c,
                           char *in, int in_size) {
    int64_t *syms = (int64_t *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, c->u.e_xpack.map[syms[i]], c->u.e_xpack.nbits);

    return r;
}

int cram_xpack_encode_int(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, c->u.e_xpack.map[syms[i]], c->u.e_xpack.nbits);

    return r;
}

int cram_xpack_encode_char(cram_slice *slice, cram_codec *c,
                           char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    return 0;

 block_err:
    return -1;
}

void cram_xpack_encode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.e_xpack.sub_codec)
        c->u.e_xpack.sub_codec->free(c->u.e_xpack.sub_codec);

    cram_free_block(c->out);

    free(c);
}

cram_codec *cram_xpack_encode_init(cram_stats *st,
                                   enum cram_encoding codec,
                                   enum cram_external_type option,
                                   void *dat,
                                   int version, varint_vec *vv) {
    cram_codec *c;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_XPACK;
    c->free   = cram_xpack_encode_free;
    if (option == E_LONG)
        c->encode = cram_xpack_encode_long;
    else if (option == E_INT)
        c->encode = cram_xpack_encode_int;
    else
        c->encode = cram_xpack_encode_char;
    c->store  = cram_xpack_encode_store;
    c->flush  = cram_xpack_encode_flush;

    cram_xpack_encoder *e = (cram_xpack_encoder *)dat;
    c->u.e_xpack.nbits = e->nbits;
    c->u.e_xpack.nval = e->nval;
    c->u.e_xpack.sub_codec = cram_encoder_init(e->sub_encoding, NULL,
                                               E_BYTE_ARRAY, e->sub_codec_dat,
                                               version, vv);

    // Initialise fwd and rev maps
    memcpy(c->u.e_xpack.map, e->map, sizeof(e->map)); // P,A,C,K to 0,1,2,3
    int i, n;
    for (i = n = 0; i < 256; i++)
        if (e->map[i] != -1)
            c->u.e_xpack.rmap[n++] = i;               // 0,1,2,3 to P,A,C,K
    if (n != e->nval) {
        fprintf(stderr, "Incorrectly specified number of map items in PACK\n");
        return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * XDELTA: subtract successive values, zig-zag to turn +/- to + only,
 * and then var-int encode the result.
 *
 * This also has the additional requirement that the data series is not
 * interleaved with another, permitting efficient encoding and decoding
 * of all elements enmasse instead of needing to only extract the bits
 * necessary per item.
 */

static uint8_t  zigzag8 (int8_t  x) { return (x << 1) ^ (x >>  7); }
static uint16_t zigzag16(int16_t x) { return (x << 1) ^ (x >> 15); }
static uint32_t zigzag32(int32_t x) { return (x << 1) ^ (x >> 31); }

//static int8_t  unzigzag8 (uint8_t  x) { return (x >> 1) ^ -(x & 1); }
static int16_t unzigzag16(uint16_t x) { return (x >> 1) ^ -(x & 1); }
static int32_t unzigzag32(uint32_t x) { return (x >> 1) ^ -(x & 1); }

int cram_xdelta_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    return -1;
}

int cram_xdelta_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // Slow value-by-value method for now
    uint32_t *out32 = (uint32_t *)out;
    int i;
    for (i = 0; i < *out_size; i++) {
        uint32_t v;
        int one = 1;
        if (c->u.e_xdelta.sub_codec->decode(slice, c->u.e_xdelta.sub_codec, in,
                                          (char *)&v, &one) < 0)
            return -1;
        uint32_t d = unzigzag32(v);
        c->u.xdelta.last = out32[i] = d + c->u.xdelta.last;
    }

    return 0;
}

static int cram_xdelta_decode_expand_char(cram_slice *slice, cram_codec *c) {
    return -1;
}

int cram_xdelta_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    return -1;
}

static inline int16_t le_int2(int16_t i) {
    int16_t s;
    i16_to_le(i, (uint8_t *)&s);
    return s;
}

int cram_xdelta_decode_block(cram_slice *slice, cram_codec *c, cram_block *in,
                             char *out_, int *out_size) {
    cram_block *out = (cram_block *)out_;
    cram_block *b = c->u.e_xdelta.sub_codec->get_block(slice, c->u.e_xdelta.sub_codec);
    int i = 0;

    const int w = c->u.xdelta.word_size;
    uint32_t npad = (w - *out_size%w)%w;
    uint32_t out_sz = *out_size + npad;
    c->u.xdelta.last = 0;  // reset for each new array

    for (i = 0; i < out_sz; i += w) {
        uint16_t v;
        // Need better interface
        char *cp = (char *)b->data + b->byte;
        char *cp_end = (char *)b->data + b->uncomp_size;
        int err = 0;
        v = c->vv->varint_get32(&cp, cp_end, &err);
        if (err)
            return -1;
        b->byte = cp - (char *)b->data;

        switch(w) {
        case 2: {
            int16_t d = unzigzag16(v), z;
            c->u.xdelta.last = d + c->u.xdelta.last;
            z = le_int2(c->u.xdelta.last);
            BLOCK_APPEND(out, &z, 2-npad);
            npad = 0;
            break;
        }
        default:
            fprintf(stderr, "Unsupported word size by XDELTA\n");
            return -1;
        }
    }

    return 0;

 block_err:
    return -1;
}

void cram_xdelta_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.xdelta.sub_codec)
        c->u.xdelta.sub_codec->free(c->u.xdelta.sub_codec);

    free(c);
}

int cram_xdelta_decode_size(cram_slice *slice, cram_codec *c) {
    cram_xdelta_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id]->uncomp_size;
}

cram_block *cram_xdelta_get_block(cram_slice *slice, cram_codec *c) {
    cram_xdelta_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id];
}

cram_codec *cram_xdelta_decode_init(cram_block_compression_hdr *hdr,
                                    char *data, int size,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;
    char *endp = data+size;

    if (!(c = calloc(1, sizeof(*c))))
        return NULL;

    c->codec  = E_XDELTA;
    if (option == E_LONG)
        c->decode = cram_xdelta_decode_long;
    else if (option == E_INT)
        c->decode = cram_xdelta_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_xdelta_decode_char;
    else if (option == E_BYTE_ARRAY_BLOCK) {
        option = E_BYTE_ARRAY;
        c->decode = cram_xdelta_decode_block;
    } else {
        free(c);
        return NULL;
    }
    c->free = cram_xdelta_decode_free;
    c->size = cram_xdelta_decode_size;
    c->get_block = cram_xdelta_get_block;
    c->describe = NULL;

    c->u.xdelta.word_size = vv->varint_get32(&cp, endp, NULL);
    c->u.xdelta.last = 0;

    int encoding = vv->varint_get32(&cp, endp, NULL);
    int sub_size = vv->varint_get32(&cp, endp, NULL);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.xdelta.sub_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
                                              option, version, vv);
    if (c->u.xdelta.sub_codec == NULL)
        goto malformed;
    cp += sub_size;

    if (cp - data != size) {
    malformed:
        fprintf(stderr, "Malformed xdelta header stream\n");
        cram_xdelta_decode_free(c);
        return NULL;
    }

    return c;
}

int cram_xdelta_encode_flush(cram_codec *c) {
    int r = -1;
    cram_block *b = cram_new_block(0, 0);
    if (!b)
        return -1;

    switch (c->u.e_xdelta.word_size) {
    case 2: {
        // Delta + zigzag transform.
        // Subtracting two 8-bit values has a 9-bit result (-255 to 255).
        // However think of it as turning a wheel clockwise or anti-clockwise.
        // If it has 256 gradations then a -ve rotation followed by a +ve
        // rotation of the same amount reverses it regardless.
        //
        // Similarly the zig-zag transformation doesn't invent any extra bits,
        // so the entire thing can be done in-situ.  This may permit faster
        // SIMD loops if we break apart the steps.

        // uint16_t last = 0, d;
        // for (i = 0; i < n; i++) {
        //     d = io[i] - last;
        //     last = io[i];
        //     io[i] = zigzag16(vd);
        // }

        // --- vs ---

        // for (i = n-1; i >= 1; i--)
        //     io[i] -= io[i-1];
        // for (i = 0; i < n; i++)
        //     io[i] = zigzag16(io[i]);

        // varint: need array variant for speed here.
        // With zig-zag
        int i, n = BLOCK_SIZE(c->out)/2;;
        uint16_t *dat = (uint16_t *)BLOCK_DATA(c->out), last = 0;

        if (n*2 < BLOCK_SIZE(c->out)) {
            // half word
            last = *(uint8_t *)dat;
            c->vv->varint_put32_blk(b, zigzag16(last));
            dat = (uint16_t *)(((uint8_t *)dat)+1);
        }

        for (i = 0; i < n; i++) {
            uint16_t d = dat[i] - last; // possibly unaligned
            last = dat[i];
            c->vv->varint_put32_blk(b, zigzag16(d));
        }

        break;
    }

    case 4: {
        int i, n = BLOCK_SIZE(c->out)/4;;
        uint32_t *dat = (uint32_t *)BLOCK_DATA(c->out), last = 0;

        for (i = 0; i < n; i++) {
            uint32_t d = dat[i] - last;
            last = dat[i];
            c->vv->varint_put32_blk(b, zigzag32(d));
        }

        break;
    }

    case 1: {
        int i, n = BLOCK_SIZE(c->out);;
        uint8_t *dat = (uint8_t *)BLOCK_DATA(c->out), last = 0;

        for (i = 0; i < n; i++) {
            uint32_t d = dat[i] - last;
            last = dat[i];
            c->vv->varint_put32_blk(b, zigzag8(d));
        }

        break;
    }

    default:
        goto err;
    }

    if (c->u.e_xdelta.sub_codec->encode(NULL, c->u.e_xdelta.sub_codec,
                                      (char *)b->data, b->byte))
        goto err;

    r = 0;

 err:
    cram_free_block(b);
    return r;

}

int cram_xdelta_encode_store(cram_codec *c, cram_block *b,
                            char *prefix, int version) {
    int len = 0, r = 0, n;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    // Store sub-codec
    cram_codec *tc = c->u.e_xdelta.sub_codec;
    cram_block *tb = cram_new_block(0, 0);
    if (!tb)
        return -1;
    int len2 = tc->store(tc, tb, NULL, version);

    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;

    // codec length
    len += (n = c->vv->varint_put32_blk(b, c->vv->varint_size(c->u.e_xdelta.word_size)
                                        + len2)); r |= n;

    // This and sub-codec
    len += (n = c->vv->varint_put32_blk(b, c->u.e_xdelta.word_size)); r |= n;
    BLOCK_APPEND(b, BLOCK_DATA(tb), BLOCK_SIZE(tb));

    cram_free_block(tb);

    return r > 0 ? len + len2 : -1;

 block_err:
    return -1;
}

// Same as cram_beta_encode_long
int cram_xdelta_encode_long(cram_slice *slice, cram_codec *c,
                           char *in, int in_size) {
    return -1;
}

int cram_xdelta_encode_int(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    return -1;
}

int cram_xdelta_encode_char(cram_slice *slice, cram_codec *c,
                            char *in, int in_size) {
    char *dat = malloc(in_size*5);
    if (!dat)
        return -1;
    char *cp = dat, *cp_end = dat + in_size*5;

    c->u.e_xdelta.last = 0; // reset for each new array
    if (c->u.e_xdelta.word_size == 2) {
        int i, part;

        part = in_size%2;
        if (part) {
            uint16_t z = in[0];
            c->u.e_xdelta.last = le_int2(z);
            cp += c->vv->varint_put32(cp, cp_end, zigzag16(c->u.e_xdelta.last));
        }

        uint16_t *in16 = (uint16_t *)(in+part);
        for (i = 0; i < in_size/2; i++) {
            uint16_t d = le_int2(in16[i]) - c->u.e_xdelta.last;
            c->u.e_xdelta.last = le_int2(in16[i]);
            cp += c->vv->varint_put32(cp, cp_end, zigzag16(d));
        }
    }
    if (c->u.e_xdelta.sub_codec->encode(slice, c->u.e_xdelta.sub_codec,
                                      (char *)dat, cp-dat)) {
        free(dat);
        return -1;
    }

    free(dat);
    return 0;
}

void cram_xdelta_encode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.e_xdelta.sub_codec)
        c->u.e_xdelta.sub_codec->free(c->u.e_xdelta.sub_codec);

    cram_free_block(c->out);

    free(c);
}

cram_codec *cram_xdelta_encode_init(cram_stats *st,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    void *dat,
                                    int version, varint_vec *vv) {
    cram_codec *c;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_XDELTA;
    c->free   = cram_xdelta_encode_free;
    if (option == E_LONG)
        c->encode = cram_xdelta_encode_long;
    else if (option == E_INT)
        c->encode = cram_xdelta_encode_int;
    else
        c->encode = cram_xdelta_encode_char;
    c->store  = cram_xdelta_encode_store;
    c->flush  = cram_xdelta_encode_flush;

    cram_xdelta_encoder *e = (cram_xdelta_encoder *)dat;
    c->u.e_xdelta.word_size = e->word_size;
    c->u.e_xdelta.last = 0;
    c->u.e_xdelta.sub_codec = cram_encoder_init(e->sub_encoding, NULL,
                                                E_BYTE_ARRAY,
                                                e->sub_codec_dat,
                                                version, vv);

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * XRLE
 *
 * This also has the additional requirement that the data series is not
 * interleaved with another, permitting efficient encoding and decoding
 * of all elements enmasse instead of needing to only extract the bits
 * necessary per item.
 */
int cram_xrle_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // TODO if and when needed
    return -1;
}

int cram_xrle_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // TODO if and when needed
    return -1;
}

// Expands an XRLE transform and caches result in slice->block_by_id[]
static int cram_xrle_decode_expand_char(cram_slice *slice, cram_codec *c) {
    cram_block *b = slice->block_by_id[512 + c->codec_id];
    if (b)
        return 0;

    b = slice->block_by_id[512 + c->codec_id] = cram_new_block(0, 0);
    if (!b)
        return -1;
    cram_block *lit_b = c->u.xrle.lit_codec->get_block(slice, c->u.xrle.lit_codec);
    if (!lit_b)
        return -1;
    unsigned char *lit_dat = lit_b->data;
    unsigned int lit_sz = lit_b->uncomp_size;
    unsigned int len_sz = c->u.xrle.len_codec->size(slice, c->u.xrle.len_codec);

    cram_block *len_b = c->u.xrle.len_codec->get_block(slice, c->u.xrle.len_codec);
    if (!len_b)
        return -1;
    unsigned char *len_dat = len_b->data;

    uint8_t rle_syms[256];
    int rle_nsyms = 0;
    int i;
    for (i = 0; i < 256; i++) {
        if (c->u.xrle.rep_score[i] > 0)
            rle_syms[rle_nsyms++] = i;
    }

    uint64_t out_sz;
    int nb = var_get_u64(len_dat, len_dat+len_sz, &out_sz);
    if (!(b->data = malloc(out_sz)))
        return -1;
    hts_rle_decode(lit_dat, lit_sz,
                   len_dat+nb, len_sz-nb,
                   rle_syms, rle_nsyms,
                   b->data, &out_sz);
    b->uncomp_size = out_sz;

    return 0;
}

int cram_xrle_decode_size(cram_slice *slice, cram_codec *c) {
    cram_xrle_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id]->uncomp_size;
}

cram_block *cram_xrle_get_block(cram_slice *slice, cram_codec *c) {
    cram_xrle_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id];
}

int cram_xrle_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int n = *out_size;

    cram_xrle_decode_expand_char(slice, c);
    cram_block *b = slice->block_by_id[512 + c->codec_id];

    memcpy(out, b->data + b->idx, n);
    b->idx += n;
    return 0;

    // Old code when not cached
    while (n > 0) {
        if (c->u.xrle.cur_len == 0) {
            unsigned char lit;
            int one = 1;
            if (c->u.xrle.lit_codec->decode(slice, c->u.xrle.lit_codec, in,
                                          (char *)&lit, &one) < 0)
                return -1;
            c->u.xrle.cur_lit = lit;

            if (c->u.xrle.rep_score[lit] > 0) {
                if (c->u.xrle.len_codec->decode(slice, c->u.xrle.len_codec, in,
                                              (char *)&c->u.xrle.cur_len, &one) < 0)
                    return -1;
            } // else cur_len still zero
            //else fprintf(stderr, "%d\n", lit);

            c->u.xrle.cur_len++;
        }

        if (n >= c->u.xrle.cur_len) {
            memset(out, c->u.xrle.cur_lit, c->u.xrle.cur_len);
            out += c->u.xrle.cur_len;
            n -= c->u.xrle.cur_len;
            c->u.xrle.cur_len = 0;
        } else {
            memset(out, c->u.xrle.cur_lit, n);
            out += n;
            c->u.xrle.cur_len -= n;
            n = 0;
        }
    }

    return 0;
}

void cram_xrle_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.xrle.len_codec)
        c->u.xrle.len_codec->free(c->u.xrle.len_codec);

    if (c->u.xrle.lit_codec)
        c->u.xrle.lit_codec->free(c->u.xrle.lit_codec);

    free(c);
}

cram_codec *cram_xrle_decode_init(cram_block_compression_hdr *hdr,
                                  char *data, int size,
                                  enum cram_encoding codec,
                                  enum cram_external_type option,
                                  int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;
    char *endp = data+size;
    int err = 0;

    if (!(c = calloc(1, sizeof(*c))))
        return NULL;

    c->codec  = E_XRLE;
    if (option == E_LONG)
        c->decode = cram_xrle_decode_long;
    else if (option == E_INT)
        c->decode = cram_xrle_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_xrle_decode_char;
    else {
        fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
        free(c);
        return NULL;
    }
    c->free   = cram_xrle_decode_free;
    c->size   = cram_xrle_decode_size;
    c->get_block = cram_xrle_get_block;
    c->describe = NULL;
    c->u.xrle.cur_len = 0;
    c->u.xrle.cur_lit = -1;

    // RLE map
    int i, j, nrle = vv->varint_get32(&cp, endp, &err);
    memset(c->u.xrle.rep_score, 0, 256*sizeof(*c->u.xrle.rep_score));
    for (i = 0; i < nrle && i < 256; i++) {
        j = vv->varint_get32(&cp, endp, &err);
        if (j >= 0 && j < 256)
            c->u.xrle.rep_score[j] = 1;
    }

    // Length and literal sub encodings
    c->u.xrle.len_encoding = vv->varint_get32(&cp, endp, &err);
    int sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.xrle.len_codec = cram_decoder_init(hdr, c->u.xrle.len_encoding,
                                            cp, sub_size, E_INT, version, vv);
    if (c->u.xrle.len_codec == NULL)
        goto malformed;
    cp += sub_size;

    c->u.xrle.lit_encoding = vv->varint_get32(&cp, endp, &err);
    sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.xrle.lit_codec = cram_decoder_init(hdr, c->u.xrle.lit_encoding,
                                            cp, sub_size, option, version, vv);
    if (c->u.xrle.lit_codec == NULL)
        goto malformed;
    cp += sub_size;

    if (err)
        goto malformed;

    return c;

 malformed:
    fprintf(stderr, "Malformed xrle header stream\n");
    cram_xrle_decode_free(c);
    return NULL;
}

int cram_xrle_encode_flush(cram_codec *c) {
    uint8_t *out_lit, *out_len;
    uint64_t out_lit_size, out_len_size;
    uint8_t rle_syms[256];
    int rle_nsyms = 0, i;

    for (i = 0; i < 256; i++)
        if (c->u.e_xrle.rep_score[i] > 0)
            rle_syms[rle_nsyms++] = i;

    if (!c->u.e_xrle.to_flush) {
        c->u.e_xrle.to_flush = (char *)BLOCK_DATA(c->out);
        c->u.e_xrle.to_flush_size = BLOCK_SIZE(c->out);
    }

    out_len = malloc(c->u.e_xrle.to_flush_size+8);
    if (!out_len)
        return -1;

    int nb = var_put_u64(out_len, NULL, c->u.e_xrle.to_flush_size);

    out_lit = hts_rle_encode((uint8_t *)c->u.e_xrle.to_flush, c->u.e_xrle.to_flush_size,
                             out_len+nb, &out_len_size,
                             rle_syms, &rle_nsyms,
                             NULL, &out_lit_size);
    out_len_size += nb;


    // TODO: can maybe "gift" the sub codec the data block, to remove
    // one level of memcpy.
    if (c->u.e_xrle.len_codec->encode(NULL,
                                      c->u.e_xrle.len_codec,
                                      (char *)out_len, out_len_size))
        return -1;

    if (c->u.e_xrle.lit_codec->encode(NULL,
                                      c->u.e_xrle.lit_codec,
                                      (char *)out_lit, out_lit_size))
        return -1;

    free(out_len);
    free(out_lit);

    return 0;
}

int cram_xrle_encode_store(cram_codec *c, cram_block *b,
                            char *prefix, int version) {
    int len = 0, r = 0, n;
    cram_codec *tc;
    cram_block *b_rle, *b_len, *b_lit;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    // List of symbols to RLE
    b_rle = cram_new_block(0, 0);
    if (!b_rle)
        return -1;
    int i, nrle = 0, len1 = 0;
    for (i = 0; i < 256; i++) {
        if (c->u.e_xrle.rep_score[i] > 0) {
            nrle++;
            len1 += (n = c->vv->varint_put32_blk(b_rle,i)); r |= n;
        }
    }

    // Store length and literal sub-codecs to get encoded length
    tc = c->u.e_xrle.len_codec;
    b_len = cram_new_block(0, 0);
    if (!b_len)
        return -1;
    int len2 = tc->store(tc, b_len, NULL, version);

    tc = c->u.e_xrle.lit_codec;
    b_lit = cram_new_block(0, 0);
    if (!b_lit)
        return -1;
    int len3 = tc->store(tc, b_lit, NULL, version);

    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;
    len += (n = c->vv->varint_put32_blk(b, len1 + len2 + len3
                                        + c->vv->varint_size(nrle))); r |= n;
    len += (n = c->vv->varint_put32_blk(b, nrle)); r |= n;
    BLOCK_APPEND(b, BLOCK_DATA(b_rle), BLOCK_SIZE(b_rle));
    BLOCK_APPEND(b, BLOCK_DATA(b_len), BLOCK_SIZE(b_len));
    BLOCK_APPEND(b, BLOCK_DATA(b_lit), BLOCK_SIZE(b_lit));

    cram_free_block(b_rle);
    cram_free_block(b_len);
    cram_free_block(b_lit);

    if (r > 0)
        return len + len1 + len2 + len3;

 block_err:
    return -1;
}

int cram_xrle_encode_long(cram_slice *slice, cram_codec *c,
                           char *in, int in_size) {
    // TODO if and when needed
    return -1;
}

int cram_xrle_encode_int(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    // TODO if and when needed
    return -1;
}

int cram_xrle_encode_char(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    if (c->u.e_xrle.to_flush) {
        if (!c->out && !(c->out = cram_new_block(0, 0)))
            return -1;
        BLOCK_APPEND(c->out, c->u.e_xrle.to_flush, c->u.e_xrle.to_flush_size);
        c->u.e_xrle.to_flush = NULL;
        c->u.e_xrle.to_flush_size = 0;
    }

    if (c->out && BLOCK_SIZE(c->out) > 0) {
        // Gathering data
        BLOCK_APPEND(c->out, in, in_size);
        return 0;
    }

    // else cache copy of the data we're about to send to flush instead.
    c->u.e_xrle.to_flush = in;
    c->u.e_xrle.to_flush_size = in_size;
    return 0;

 block_err:
    return -1;
}

void cram_xrle_encode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.e_xrle.len_codec)
        c->u.e_xrle.len_codec->free(c->u.e_xrle.len_codec);
    if (c->u.e_xrle.lit_codec)
        c->u.e_xrle.lit_codec->free(c->u.e_xrle.lit_codec);

    cram_free_block(c->out);

    free(c);
}

cram_codec *cram_xrle_encode_init(cram_stats *st,
                                  enum cram_encoding codec,
                                  enum cram_external_type option,
                                  void *dat,
                                  int version, varint_vec *vv) {
    cram_codec *c;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_XRLE;
    c->free   = cram_xrle_encode_free;
    if (option == E_LONG)
        c->encode = cram_xrle_encode_long;
    else if (option == E_INT)
        c->encode = cram_xrle_encode_int;
    else
        c->encode = cram_xrle_encode_char;
    c->store  = cram_xrle_encode_store;
    c->flush  = cram_xrle_encode_flush;

    cram_xrle_encoder *e = (cram_xrle_encoder *)dat;

    c->u.e_xrle.len_codec = cram_encoder_init(e->len_encoding, NULL,
                                              E_BYTE, e->len_dat,
                                              version, vv);
    c->u.e_xrle.lit_codec = cram_encoder_init(e->lit_encoding, NULL,
                                              E_BYTE, e->lit_dat,
                                              version, vv);
    c->u.e_xrle.cur_lit = -1;
    c->u.e_xrle.cur_len = -1;
    c->u.e_xrle.to_flush = NULL;
    c->u.e_xrle.to_flush_size = 0;

    memcpy(c->u.e_xrle.rep_score, e->rep_score, 256*sizeof(*c->u.e_xrle.rep_score));

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * SUBEXP
 */
int cram_subexp_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int n, count;
    int k = c->u.subexp.k;

    for (count = 0, n = *out_size; count < n; count++) {
        int i = 0, tail;
        int val;

        /* Get number of 1s */
        //while (get_bit_MSB(in) == 1) i++;
        i = get_one_bits_MSB(in);
        if (i < 0 || cram_not_enough_bits(in, i > 0 ? i + k - 1 : k))
            return -1;
        /*
         * Val is
         * i > 0:  2^(k+i-1) + k+i-1 bits
         * i = 0:  k bits
         */
        if (i) {
            tail = i + k-1;
            val = 0;
            while (tail) {
                //val = val<<1; val |= get_bit_MSB(in);
                GET_BIT_MSB(in, val);
                tail--;
            }
            val += 1 << (i + k-1);
        } else {
            tail = k;
            val = 0;
            while (tail) {
                //val = val<<1; val |= get_bit_MSB(in);
                GET_BIT_MSB(in, val);
                tail--;
            }
        }

        out_i[count] = val - c->u.subexp.offset;
    }

    return 0;
}

void cram_subexp_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

int cram_subexp_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "SUBEXP(offset=%d,k=%d)",
                    c->u.subexp.offset,
                    c->u.subexp.k)
        < 0 ? -1 : 0;
}

cram_codec *cram_subexp_decode_init(cram_block_compression_hdr *hdr,
                                    char *data, int size,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (option != E_INT) {
        hts_log_error("This codec only supports INT encodings");
        return NULL;
    }

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_SUBEXP;
    c->decode = cram_subexp_decode;
    c->free   = cram_subexp_decode_free;
    c->describe = cram_subexp_describe;
    c->u.subexp.k = -1;

    c->u.subexp.offset = vv->varint_get32(&cp, data + size, NULL);
    c->u.subexp.k      = vv->varint_get32(&cp, data + size, NULL);

    if (cp - data != size || c->u.subexp.k < 0) {
        hts_log_error("Malformed subexp header stream");
        free(c);
        return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * GAMMA
 */
int cram_gamma_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++) {
        int nz = 0;
        int val;
        //while (get_bit_MSB(in) == 0) nz++;
        nz = get_zero_bits_MSB(in);
        if (cram_not_enough_bits(in, nz))
            return -1;
        val = 1;
        while (nz > 0) {
            //val <<= 1; val |= get_bit_MSB(in);
            GET_BIT_MSB(in, val);
            nz--;
        }

        out_i[i] = val - c->u.gamma.offset;
    }

    return 0;
}

void cram_gamma_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

int cram_gamma_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "GAMMA(offset=%d)", c->u.subexp.offset)
        < 0 ? -1 : 0;
}

cram_codec *cram_gamma_decode_init(cram_block_compression_hdr *hdr,
                                   char *data, int size,
                                   enum cram_encoding codec,
                                   enum cram_external_type option,
                                   int version, varint_vec *vv) {
    cram_codec *c = NULL;
    char *cp = data;

    if (option != E_INT) {
        hts_log_error("This codec only supports INT encodings");
        return NULL;
    }

    if (size < 1)
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_GAMMA;
    c->decode = cram_gamma_decode;
    c->free   = cram_gamma_decode_free;
    c->describe = cram_gamma_describe;

    c->u.gamma.offset = vv->varint_get32(&cp, data+size, NULL);

    if (cp - data != size)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed gamma header stream");
    free(c);
    return NULL;
}

/*
 * ---------------------------------------------------------------------------
 * HUFFMAN
 */

static int code_sort(const void *vp1, const void *vp2) {
    const cram_huffman_code *c1 = (const cram_huffman_code *)vp1;
    const cram_huffman_code *c2 = (const cram_huffman_code *)vp2;

    if (c1->len != c2->len)
        return c1->len - c2->len;
    else
        return c1->symbol < c2->symbol ? -1 : (c1->symbol > c2->symbol ? 1 : 0);
}

void cram_huffman_decode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->u.huffman.codes)
        free(c->u.huffman.codes);
    free(c);
}

int cram_huffman_decode_null(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    return -1;
}

int cram_huffman_decode_char0(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out, int *out_size) {
    int i, n;

    if (!out)
        return 0;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
        out[i] = c->u.huffman.codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_char(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int i, n, ncodes = c->u.huffman.ncodes;
    const cram_huffman_code * const codes = c->u.huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
        int idx = 0;
        int val = 0, len = 0, last_len = 0;

        for (;;) {
            int dlen = codes[idx].len - last_len;
            if (cram_not_enough_bits(in, dlen))
                return -1;

            //val <<= dlen;
            //val  |= get_bits_MSB(in, dlen);
            //last_len = (len += dlen);

            last_len = (len += dlen);
            for (; dlen; dlen--) GET_BIT_MSB(in, val);

            idx = val - codes[idx].p;
            if (idx >= ncodes || idx < 0)
                return -1;

            if (codes[idx].code == val && codes[idx].len == len) {
                if (out) out[i] = codes[idx].symbol;
                break;
            }
        }
    }

    return 0;
}

int cram_huffman_decode_int0(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->u.huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
        out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_int(cram_slice *slice, cram_codec *c,
                            cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n, ncodes = c->u.huffman.ncodes;
    const cram_huffman_code * const codes = c->u.huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
        int idx = 0;
        int val = 0, len = 0, last_len = 0;

        // Now one bit at a time for remaining checks
        for (;;) {
            int dlen = codes[idx].len - last_len;
            if (cram_not_enough_bits(in, dlen))
                return -1;

            //val <<= dlen;
            //val  |= get_bits_MSB(in, dlen);
            //last_len = (len += dlen);

            last_len = (len += dlen);
            for (; dlen; dlen--) GET_BIT_MSB(in, val);

            idx = val - codes[idx].p;
            if (idx >= ncodes || idx < 0)
                return -1;

            if (codes[idx].code == val && codes[idx].len == len) {
                out_i[i] = codes[idx].symbol;
                break;
            }
        }
    }

    return 0;
}

int cram_huffman_decode_long0(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->u.huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
        out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_long(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n, ncodes = c->u.huffman.ncodes;
    const cram_huffman_code * const codes = c->u.huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
        int idx = 0;
        int val = 0, len = 0, last_len = 0;

        // Now one bit at a time for remaining checks
        for (;;) {
            int dlen = codes[idx].len - last_len;
            if (cram_not_enough_bits(in, dlen))
                return -1;

            //val <<= dlen;
            //val  |= get_bits_MSB(in, dlen);
            //last_len = (len += dlen);

            last_len = (len += dlen);
            for (; dlen; dlen--) GET_BIT_MSB(in, val);

            idx = val - codes[idx].p;
            if (idx >= ncodes || idx < 0)
                return -1;

            if (codes[idx].code == val && codes[idx].len == len) {
                out_i[i] = codes[idx].symbol;
                break;
            }
        }
    }

    return 0;
}

int cram_huffman_describe(cram_codec *c, kstring_t *ks) {
    int r = 0, n;
    r |= ksprintf(ks, "HUFFMAN(codes={") < 0;
    for (n = 0; n < c->u.huffman.ncodes; n++) {
        r |= ksprintf(ks, "%s%"PRId64, n?",":"",
                      c->u.huffman.codes[n].symbol);
    }
    r |= ksprintf(ks, "},lengths={") < 0;
    for (n = 0; n < c->u.huffman.ncodes; n++) {
        r |= ksprintf(ks, "%s%d", n?",":"",
                      c->u.huffman.codes[n].len);
    }
    r |= ksprintf(ks, "})") < 0;
    return r;
}

/*
 * Initialises a huffman decoder from an encoding data stream.
 */
cram_codec *cram_huffman_decode_init(cram_block_compression_hdr *hdr,
                                     char *data, int size,
                                     enum cram_encoding codec,
                                     enum cram_external_type option,
                                     int version, varint_vec *vv) {
    int32_t ncodes = 0, i, j;
    char *cp = data, *data_end = &data[size];
    cram_codec *h;
    cram_huffman_code *codes = NULL;
    int32_t val, last_len, max_len = 0;
    uint32_t max_val; // needs one more bit than val
    const int max_code_bits = sizeof(val) * 8 - 1;
    int err = 0;

    if (option == E_BYTE_ARRAY_BLOCK) {
        hts_log_error("BYTE_ARRAYs not supported by this codec");
        return NULL;
    }

    ncodes = vv->varint_get32(&cp, data_end, &err);
    if (ncodes < 0) {
        hts_log_error("Invalid number of symbols in huffman stream");
        return NULL;
    }
    if (ncodes >= SIZE_MAX / sizeof(*codes)) {
        errno = ENOMEM;
        return NULL;
    }
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (ncodes > FUZZ_ALLOC_LIMIT / sizeof(*codes)) {
        errno = ENOMEM;
        return NULL;
    }
#endif
    h = calloc(1, sizeof(*h));
    if (!h)
        return NULL;

    h->codec  = E_HUFFMAN;
    h->free   = cram_huffman_decode_free;

    h->u.huffman.ncodes = ncodes;
    h->u.huffman.option = option;
    if (ncodes) {
        codes = h->u.huffman.codes = malloc(ncodes * sizeof(*codes));
        if (!codes) {
            free(h);
            return NULL;
        }
    } else {
        codes = h->u.huffman.codes = NULL;
    }

    /* Read symbols and bit-lengths */
    if (option == E_LONG) {
        for (i = 0; i < ncodes; i++)
            codes[i].symbol = vv->varint_get64(&cp, data_end, &err);
    } else if (option == E_INT || option == E_BYTE) {
        for (i = 0; i < ncodes; i++)
            codes[i].symbol = vv->varint_get32(&cp, data_end, &err);
    } else {
        goto malformed;
    }

    if (err)
        goto malformed;

    i = vv->varint_get32(&cp, data_end, &err);
    if (i != ncodes)
        goto malformed;

    if (ncodes == 0) {
        /* NULL huffman stream.  Ensure it returns an error if
           anything tries to use it. */
        h->decode = cram_huffman_decode_null;
        return h;
    }

    for (i = 0; i < ncodes; i++) {
        codes[i].len = vv->varint_get32(&cp, data_end, &err);
        if (err)
            break;
        if (codes[i].len < 0) {
            hts_log_error("Huffman code length (%d) is negative", codes[i].len);
            goto malformed;
        }
        if (max_len < codes[i].len)
            max_len = codes[i].len;
    }
    if (err || cp - data != size || max_len >= ncodes)
        goto malformed;

    /* 31 is max. bits available in val */
    if (max_len > max_code_bits) {
        hts_log_error("Huffman code length (%d) is greater "
                      "than maximum supported (%d)", max_len, max_code_bits);
        goto malformed;
    }

    /* Sort by bit length and then by symbol value */
    qsort(codes, ncodes, sizeof(*codes), code_sort);

    /* Assign canonical codes */
    val = -1, last_len = 0, max_val = 0;
    for (i = 0; i < ncodes; i++) {
        val++;
        if (val > max_val)
            goto malformed;

        if (codes[i].len > last_len) {
            val <<= (codes[i].len - last_len);
            last_len = codes[i].len;
            max_val = (1U << codes[i].len) - 1;
        }
        codes[i].code = val;
    }

    /*
     * Compute the next starting point, offset by the i'th value.
     * For example if codes 10, 11, 12, 13 are 30, 31, 32, 33 then
     * codes[10..13].p = 30 - 10.
     */
    last_len = 0;
    for (i = j = 0; i < ncodes; i++) {
        if (codes[i].len > last_len) {
            j = codes[i].code - i;
            last_len = codes[i].len;
        }
        codes[i].p = j;
    }

    // puts("==HUFF LEN==");
    // for (i = 0; i <= last_len+1; i++) {
    //     printf("len %d=%d prefix %d\n", i, h->u.huffman.lengths[i], h->u.huffman.prefix[i]);
    // }
    // puts("===HUFFMAN CODES===");
    // for (i = 0; i < ncodes; i++) {
    //     int j;
    //     printf("%d: %d %d %d ", i, codes[i].symbol, codes[i].len, codes[i].code);
    //     j = codes[i].len;
    //     while (j) {
    //         putchar(codes[i].code & (1 << --j) ? '1' : '0');
    //     }
    //     printf(" %d\n", codes[i].code);
    // }

    if (option == E_BYTE || option == E_BYTE_ARRAY) {
        if (h->u.huffman.codes[0].len == 0)
            h->decode = cram_huffman_decode_char0;
        else
            h->decode = cram_huffman_decode_char;
    } else if (option == E_LONG || option == E_SLONG) {
        if (h->u.huffman.codes[0].len == 0)
            h->decode = cram_huffman_decode_long0;
        else
            h->decode = cram_huffman_decode_long;
    } else if (option == E_INT || option == E_SINT || option == E_BYTE) {
        if (h->u.huffman.codes[0].len == 0)
            h->decode = cram_huffman_decode_int0;
        else
            h->decode = cram_huffman_decode_int;
    } else {
        return NULL;
    }
    h->describe = cram_huffman_describe;

    return (cram_codec *)h;

 malformed:
    hts_log_error("Malformed huffman header stream");
    free(codes);
    free(h);
    return NULL;
}

int cram_huffman_encode_char0(cram_slice *slice, cram_codec *c,
                              char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_char(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    int i, code, len, r = 0;
    unsigned char *syms = (unsigned char *)in;

    while (in_size--) {
        int sym = *syms++;
        if (sym >= -1 && sym < MAX_HUFF) {
            i = c->u.e_huffman.val2code[sym+1];
            assert(c->u.e_huffman.codes[i].symbol == sym);
            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        } else {
            /* Slow - use a lookup table for when sym < MAX_HUFF? */
            for (i = 0; i < c->u.e_huffman.nvals; i++) {
                if (c->u.e_huffman.codes[i].symbol == sym)
                    break;
            }
            if (i == c->u.e_huffman.nvals)
                return -1;

            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        }

        r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

int cram_huffman_encode_int0(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_int(cram_slice *slice, cram_codec *c,
                            char *in, int in_size) {
    int i, code, len, r = 0;
    int *syms = (int *)in;

    while (in_size--) {
        int sym = *syms++;

        if (sym >= -1 && sym < MAX_HUFF) {
            i = c->u.e_huffman.val2code[sym+1];
            assert(c->u.e_huffman.codes[i].symbol == sym);
            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        } else {
            /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
            for (i = 0; i < c->u.e_huffman.nvals; i++) {
                if (c->u.e_huffman.codes[i].symbol == sym)
                    break;
            }
            if (i == c->u.e_huffman.nvals)
                return -1;

            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        }

        r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

int cram_huffman_encode_long0(cram_slice *slice, cram_codec *c,
                              char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_long(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    int i, code, len, r = 0;
    int64_t *syms = (int64_t *)in;

    while (in_size--) {
        int sym = *syms++;

        if (sym >= -1 && sym < MAX_HUFF) {
            i = c->u.e_huffman.val2code[sym+1];
            assert(c->u.e_huffman.codes[i].symbol == sym);
            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        } else {
            /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
            for (i = 0; i < c->u.e_huffman.nvals; i++) {
                if (c->u.e_huffman.codes[i].symbol == sym)
                    break;
            }
            if (i == c->u.e_huffman.nvals)
                return -1;

            code = c->u.e_huffman.codes[i].code;
            len  = c->u.e_huffman.codes[i].len;
        }

        r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

void cram_huffman_encode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->u.e_huffman.codes)
        free(c->u.e_huffman.codes);
    free(c);
}

/*
 * Encodes a huffman tree.
 * Returns number of bytes written.
 */
int cram_huffman_encode_store(cram_codec *c, cram_block *b, char *prefix,
                              int version) {
    int i, len = 0, r = 0, n;
    cram_huffman_code *codes = c->u.e_huffman.codes;
    /*
     * Up to code length 127 means 2.5e+26 bytes of data required (worst
     * case huffman tree needs symbols with freqs matching the Fibonacci
     * series). So guaranteed 1 byte per code.
     *
     * Symbols themselves could be 5 bytes (eg -1 is 5 bytes in itf8).
     *
     * Therefore 6*ncodes + 5 + 5 + 1 + 5 is max memory
     */
    char *tmp = malloc(6*c->u.e_huffman.nvals+16);
    char *tp = tmp, *tpend = tmp+6*c->u.e_huffman.nvals+16;

    if (!tmp)
        return -1;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += c->vv->varint_put32(tp, tpend, c->u.e_huffman.nvals);
    if (c->u.e_huffman.option == E_LONG) {
        for (i = 0; i < c->u.e_huffman.nvals; i++) {
            tp += c->vv->varint_put64(tp, tpend, codes[i].symbol);
        }
    } else if (c->u.e_huffman.option == E_SLONG) {
        for (i = 0; i < c->u.e_huffman.nvals; i++) {
            tp += c->vv->varint_put64s(tp, tpend, codes[i].symbol);
        }
    } else if (c->u.e_huffman.option == E_INT || c->u.e_huffman.option == E_BYTE) {
        for (i = 0; i < c->u.e_huffman.nvals; i++) {
            tp += c->vv->varint_put32(tp, tpend, codes[i].symbol);
        }
    } else if (c->u.e_huffman.option == E_SINT) {
        for (i = 0; i < c->u.e_huffman.nvals; i++) {
            tp += c->vv->varint_put32s(tp, tpend, codes[i].symbol);
        }
    } else {
        return -1;
    }

    tp += c->vv->varint_put32(tp, tpend, c->u.e_huffman.nvals);
    for (i = 0; i < c->u.e_huffman.nvals; i++)
        tp += c->vv->varint_put32(tp, tpend, codes[i].len);

    len += (n = c->vv->varint_put32_blk(b, c->codec)); r |= n;
    len += (n = c->vv->varint_put32_blk(b, tp-tmp));   r |= n;
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    free(tmp);

    if (r > 0)
        return len;

 block_err:
    return -1;
}

cram_codec *cram_huffman_encode_init(cram_stats *st,
                                     enum cram_encoding codec,
                                     enum cram_external_type option,
                                     void *dat,
                                     int version, varint_vec *vv) {
    int *vals = NULL, *freqs = NULL, *lens = NULL, code, len;
    int *new_vals, *new_freqs;
    int i, max_val = 0, min_val = INT_MAX, k;
    size_t nvals, vals_alloc = 0;
    cram_codec *c;
    cram_huffman_code *codes;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_HUFFMAN;

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
        if (!st->freqs[i])
            continue;
        if (nvals >= vals_alloc) {
            vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
            new_vals  = realloc(vals,  vals_alloc * sizeof(int));
            if (!new_vals) goto nomem;
            vals = new_vals;
            new_freqs = realloc(freqs, vals_alloc * sizeof(int));
            if (!new_freqs) goto nomem;
            freqs = new_freqs;
        }
        vals[nvals] = i;
        freqs[nvals] = st->freqs[i];
        assert(st->freqs[i] > 0);
        if (max_val < i) max_val = i;
        if (min_val > i) min_val = i;
        nvals++;
    }
    if (st->h) {
        khint_t k;

        for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
            if (!kh_exist(st->h, k))
                continue;
            if (nvals >= vals_alloc) {
                vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
                new_vals  = realloc(vals,  vals_alloc * sizeof(int));
                if (!new_vals) goto nomem;
                vals = new_vals;
                new_freqs = realloc(freqs, vals_alloc * sizeof(int));
                if (!new_freqs) goto nomem;
                freqs = new_freqs;
            }
            vals[nvals]= kh_key(st->h, k);
            freqs[nvals] = kh_val(st->h, k);
            assert(freqs[nvals] > 0);
            if (max_val < i) max_val = i;
            if (min_val > i) min_val = i;
            nvals++;
        }
    }

    assert(nvals > 0);

    new_freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    if (!new_freqs) goto nomem;
    freqs = new_freqs;
    lens = calloc(2*nvals, sizeof(*lens));
    if (!lens) goto nomem;

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
        int low1 = INT_MAX, low2 = INT_MAX;
        int ind1 = 0, ind2 = 0;
        for (i = 0; i < nvals; i++) {
            if (freqs[i] < 0)
                continue;
            if (low1 > freqs[i])
                low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
            else if (low2 > freqs[i])
                low2 = freqs[i], ind2 = i;
        }
        if (low2 == INT_MAX)
            break;

        freqs[nvals] = low1 + low2;
        lens[ind1] = nvals;
        lens[ind2] = nvals;
        freqs[ind1] *= -1;
        freqs[ind2] *= -1;
        nvals++;
    }
    nvals = nvals/2+1;

    /* Assign lengths */
    for (i = 0; i < nvals; i++) {
        int code_len = 0;
        for (k = lens[i]; k; k = lens[k])
            code_len++;
        lens[i] = code_len;
        freqs[i] *= -1;
        //fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], lens[i]);
    }


    /* Sort, need in a struct */
    if (!(codes = malloc(nvals * sizeof(*codes))))
        goto nomem;
    for (i = 0; i < nvals; i++) {
        codes[i].symbol = vals[i];
        codes[i].len = lens[i];
    }
    qsort(codes, nvals, sizeof(*codes), code_sort);

    /*
     * Generate canonical codes from lengths.
     * Sort by length.
     * Start with 0.
     * Every new code of same length is +1.
     * Every new code of new length is +1 then <<1 per extra length.
     *
     * /\
     * a/\
     * /\/\
     * bcd/\
     *    ef
     *
     * a 1  0
     * b 3  4 (0+1)<<2
     * c 3  5
     * d 3  6
     * e 4  14  (6+1)<<1
     * f 5  15
     */
    code = 0; len = codes[0].len;
    for (i = 0; i < nvals; i++) {
        while (len != codes[i].len) {
            code<<=1;
            len++;
        }
        codes[i].code = code++;

        if (codes[i].symbol >= -1 && codes[i].symbol < MAX_HUFF)
            c->u.e_huffman.val2code[codes[i].symbol+1] = i;

        //fprintf(stderr, "sym %d, code %d, len %d\n",
        //      codes[i].symbol, codes[i].code, codes[i].len);
    }

    free(lens);
    free(vals);
    free(freqs);

    c->u.e_huffman.codes = codes;
    c->u.e_huffman.nvals = nvals;
    c->u.e_huffman.option = option;

    c->free = cram_huffman_encode_free;
    if (option == E_BYTE || option == E_BYTE_ARRAY) {
        if (c->u.e_huffman.codes[0].len == 0)
            c->encode = cram_huffman_encode_char0;
        else
            c->encode = cram_huffman_encode_char;
    } else if (option == E_INT || option == E_SINT) {
        if (c->u.e_huffman.codes[0].len == 0)
            c->encode = cram_huffman_encode_int0;
        else
            c->encode = cram_huffman_encode_int;
    } else if (option == E_LONG || option == E_SLONG) {
        if (c->u.e_huffman.codes[0].len == 0)
            c->encode = cram_huffman_encode_long0;
        else
            c->encode = cram_huffman_encode_long;
    } else {
        return NULL;
    }
    c->store = cram_huffman_encode_store;
    c->flush = NULL;

    return c;

 nomem:
    hts_log_error("Out of memory");
    free(vals);
    free(freqs);
    free(lens);
    free(c);
    return NULL;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_LEN
 */
int cram_byte_array_len_decode(cram_slice *slice, cram_codec *c,
                               cram_block *in, char *out,
                               int *out_size) {
    /* Fetch length */
    int32_t len = 0, one = 1;
    int r;

    r = c->u.byte_array_len.len_codec->decode(slice, c->u.byte_array_len.len_codec,
                                              in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (!r && c->u.byte_array_len.val_codec && len >= 0) {
        r = c->u.byte_array_len.val_codec->decode(slice,
                                                  c->u.byte_array_len.val_codec,
                                                  in, out, &len);
    } else {
        return -1;
    }

    *out_size = len;

    return r;
}

void cram_byte_array_len_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->u.byte_array_len.len_codec)
        c->u.byte_array_len.len_codec->free(c->u.byte_array_len.len_codec);

    if (c->u.byte_array_len.val_codec)
        c->u.byte_array_len.val_codec->free(c->u.byte_array_len.val_codec);

    free(c);
}

int cram_byte_array_len_describe(cram_codec *c, kstring_t *ks) {
    int r = 0;
    r |= ksprintf(ks, "BYTE_ARRAY_LEN(len_codec={") < 0;
    cram_byte_array_len_decoder *l = &c->u.byte_array_len;
    r |=  l->len_codec->describe
        ? l->len_codec->describe(l->len_codec, ks)
        : (ksprintf(ks, "?")<0);
    r |= ksprintf(ks, "},val_codec={") < 0;
    r |=  l->val_codec->describe
        ? l->val_codec->describe(l->val_codec, ks)
        : (ksprintf(ks, "?")<0);
    r |= ksprintf(ks, "}") < 0;

    return r;
}

cram_codec *cram_byte_array_len_decode_init(cram_block_compression_hdr *hdr,
                                            char *data, int size,
                                            enum cram_encoding codec,
                                            enum cram_external_type option,
                                            int version, varint_vec *vv) {
    cram_codec *c;
    char *cp   = data;
    char *endp = data + size;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BYTE_ARRAY_LEN;
    c->decode = cram_byte_array_len_decode;
    c->free   = cram_byte_array_len_decode_free;
    c->describe = cram_byte_array_len_describe;
    c->u.byte_array_len.len_codec = NULL;
    c->u.byte_array_len.val_codec = NULL;

    int encoding = vv->varint_get32(&cp, endp, NULL);
    int sub_size = vv->varint_get32(&cp, endp, NULL);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.byte_array_len.len_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
                                                      E_INT, version, vv);
    if (c->u.byte_array_len.len_codec == NULL)
        goto no_codec;
    cp += sub_size;

    encoding = vv->varint_get32(&cp, endp, NULL);
    sub_size = vv->varint_get32(&cp, endp, NULL);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->u.byte_array_len.val_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
                                                      option, version, vv);
    if (c->u.byte_array_len.val_codec == NULL)
        goto no_codec;
    cp += sub_size;

    if (cp - data != size)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed byte_array_len header stream");
 no_codec:
    cram_byte_array_len_decode_free(c);
    return NULL;
}

int cram_byte_array_len_encode(cram_slice *slice, cram_codec *c,
                               char *in, int in_size) {
    int32_t i32 = in_size;
    int r = 0;

    r |= c->u.e_byte_array_len.len_codec->encode(slice,
                                                 c->u.e_byte_array_len.len_codec,
                                                 (char *)&i32, 1);
    r |= c->u.e_byte_array_len.val_codec->encode(slice,
                                                 c->u.e_byte_array_len.val_codec,
                                                 in, in_size);
    return r;
}

void cram_byte_array_len_encode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->u.e_byte_array_len.len_codec)
        c->u.e_byte_array_len.len_codec->free(c->u.e_byte_array_len.len_codec);

    if (c->u.e_byte_array_len.val_codec)
        c->u.e_byte_array_len.val_codec->free(c->u.e_byte_array_len.val_codec);

    free(c);
}

int cram_byte_array_len_encode_store(cram_codec *c, cram_block *b,
                                     char *prefix, int version) {
    int len = 0, len2, len3, r = 0, n;
    cram_codec *tc;
    cram_block *b_len = NULL, *b_val = NULL;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tc = c->u.e_byte_array_len.len_codec;
    b_len = cram_new_block(0, 0);
    if (!b_len) goto block_err;
    len2 = tc->store(tc, b_len, NULL, version);
    if (len2 < 0) goto block_err;

    tc = c->u.e_byte_array_len.val_codec;
    b_val = cram_new_block(0, 0);
    if (!b_val) goto block_err;
    len3 = tc->store(tc, b_val, NULL, version);
    if (len3 < 0) goto block_err;

    len += (n = c->vv->varint_put32_blk(b, c->codec));  r |= n;
    len += (n = c->vv->varint_put32_blk(b, len2+len3)); r |= n;
    BLOCK_APPEND(b, BLOCK_DATA(b_len), BLOCK_SIZE(b_len));
    BLOCK_APPEND(b, BLOCK_DATA(b_val), BLOCK_SIZE(b_val));

    cram_free_block(b_len);
    cram_free_block(b_val);

    if (r > 0)
        return len + len2 + len3;

 block_err:
    if (b_len) cram_free_block(b_len);
    if (b_val) cram_free_block(b_val);
    return -1;
}

cram_codec *cram_byte_array_len_encode_init(cram_stats *st,
                                            enum cram_encoding codec,
                                            enum cram_external_type option,
                                            void *dat,
                                            int version, varint_vec *vv) {
    cram_codec *c;
    cram_byte_array_len_encoder *e = (cram_byte_array_len_encoder *)dat;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_BYTE_ARRAY_LEN;
    c->free = cram_byte_array_len_encode_free;
    c->encode = cram_byte_array_len_encode;
    c->store = cram_byte_array_len_encode_store;
    c->flush = NULL;

    c->u.e_byte_array_len.len_codec = cram_encoder_init(e->len_encoding,
                                                        st, E_INT,
                                                        e->len_dat,
                                                        version, vv);
    c->u.e_byte_array_len.val_codec = cram_encoder_init(e->val_encoding,
                                                        NULL, E_BYTE_ARRAY,
                                                        e->val_dat,
                                                        version, vv);

    if (!c->u.e_byte_array_len.len_codec ||
        !c->u.e_byte_array_len.val_codec) {
        cram_byte_array_len_encode_free(c);
        return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
static int cram_byte_array_stop_decode_char(cram_slice *slice, cram_codec *c,
                                            cram_block *in, char *out,
                                            int *out_size) {
    char *cp, ch;
    cram_block *b = NULL;

    b = cram_get_block_by_id(slice, c->u.byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
        return -1;

    cp = (char *)b->data + b->idx;
    if (out) {
        while ((ch = *cp) != (char)c->u.byte_array_stop.stop) {
            if (cp - (char *)b->data >= b->uncomp_size)
                return -1;
            *out++ = ch;
            cp++;
        }
    } else {
        // Consume input, but produce no output
        while ((ch = *cp) != (char)c->u.byte_array_stop.stop) {
            if (cp - (char *)b->data >= b->uncomp_size)
                return -1;
            cp++;
        }
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

int cram_byte_array_stop_decode_block(cram_slice *slice, cram_codec *c,
                                      cram_block *in, char *out_,
                                      int *out_size) {
    cram_block *b;
    cram_block *out = (cram_block *)out_;
    unsigned char *cp, *cp_end;
    unsigned char stop;

    b = cram_get_block_by_id(slice, c->u.byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
        return -1;
    cp = b->data + b->idx;
    cp_end = b->data + b->uncomp_size;

    stop = c->u.byte_array_stop.stop;
    if (cp_end - cp < out->alloc - out->byte) {
        unsigned char *out_cp = BLOCK_END(out);
        while (cp != cp_end && *cp != stop)
            *out_cp++ = *cp++;
        BLOCK_SIZE(out) = out_cp - BLOCK_DATA(out);
    } else {
        unsigned char *cp_start;
        for (cp_start = cp; cp != cp_end && *cp != stop; cp++)
            ;
        BLOCK_APPEND(out, cp_start, cp - cp_start);
        BLOCK_GROW(out, cp - cp_start);
    }

    *out_size = cp - (b->data + b->idx);
    b->idx = cp - b->data + 1;

    return 0;

 block_err:
    return -1;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

int cram_byte_array_stop_describe(cram_codec *c, kstring_t *ks) {
    return ksprintf(ks, "BYTE_ARRAY_STOP(stop=%d,id=%d)",
                    c->u.byte_array_stop.stop,
                    c->u.byte_array_stop.content_id)
        < 0 ? -1 : 0;
}

cram_codec *cram_byte_array_stop_decode_init(cram_block_compression_hdr *hdr,
                                             char *data, int size,
                                             enum cram_encoding codec,
                                             enum cram_external_type option,
                                             int version, varint_vec *vv) {
    cram_codec *c = NULL;
    unsigned char *cp = (unsigned char *)data;
    int err = 0;

    if (size < (CRAM_MAJOR_VERS(version) == 1 ? 5 : 2))
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    switch (option) {
    case E_BYTE_ARRAY_BLOCK:
        c->decode = cram_byte_array_stop_decode_block;
        break;
    case E_BYTE_ARRAY:
        c->decode = cram_byte_array_stop_decode_char;
        break;
    default:
        hts_log_error("The byte_array_stop codec only supports BYTE_ARRAYs");
        free(c);
        return NULL;
    }
    c->free   = cram_byte_array_stop_decode_free;
    c->describe = cram_byte_array_stop_describe;

    c->u.byte_array_stop.stop = *cp++;
    if (CRAM_MAJOR_VERS(version) == 1) {
        c->u.byte_array_stop.content_id = cp[0] + (cp[1]<<8) + (cp[2]<<16)
            + ((unsigned int) cp[3]<<24);
        cp += 4;
    } else {
        c->u.byte_array_stop.content_id = vv->varint_get32((char **)&cp, data+size, &err);
    }

    if ((char *)cp - data != size || err)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed byte_array_stop header stream");
    free(c);
    return NULL;
}

int cram_byte_array_stop_encode(cram_slice *slice, cram_codec *c,
                                char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    BLOCK_APPEND_CHAR(c->out, c->u.e_byte_array_stop.stop);
    return 0;

 block_err:
    return -1;
}

void cram_byte_array_stop_encode_free(cram_codec *c) {
    if (!c)
        return;
    free(c);
}

int cram_byte_array_stop_encode_store(cram_codec *c, cram_block *b,
                                      char *prefix, int version) {
    int len = 0;
    char buf[20], *cp = buf;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    cp += c->vv->varint_put32(cp, buf+20, c->codec);

    if (CRAM_MAJOR_VERS(version) == 1) {
        cp += c->vv->varint_put32(cp, buf+20, 5);
        *cp++ = c->u.e_byte_array_stop.stop;
        *cp++ = (c->u.e_byte_array_stop.content_id >>  0) & 0xff;
        *cp++ = (c->u.e_byte_array_stop.content_id >>  8) & 0xff;
        *cp++ = (c->u.e_byte_array_stop.content_id >> 16) & 0xff;
        *cp++ = (c->u.e_byte_array_stop.content_id >> 24) & 0xff;
    } else {
        cp += c->vv->varint_put32(cp, buf+20, 1 +
                                  c->vv->varint_size(c->u.e_byte_array_stop.content_id));
        *cp++ = c->u.e_byte_array_stop.stop;
        cp += c->vv->varint_put32(cp, buf+20, c->u.e_byte_array_stop.content_id);
    }

    BLOCK_APPEND(b, buf, cp-buf);
    len += cp-buf;

    return len;

 block_err:
    return -1;
}

cram_codec *cram_byte_array_stop_encode_init(cram_stats *st,
                                             enum cram_encoding codec,
                                             enum cram_external_type option,
                                             void *dat,
                                             int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_BYTE_ARRAY_STOP;
    c->free = cram_byte_array_stop_encode_free;
    c->encode = cram_byte_array_stop_encode;
    c->store = cram_byte_array_stop_encode_store;
    c->flush = NULL;

    c->u.e_byte_array_stop.stop = ((int *)dat)[0];
    c->u.e_byte_array_stop.content_id = ((int *)dat)[1];

    return c;
}

/*
 * ---------------------------------------------------------------------------
 */

const char *cram_encoding2str(enum cram_encoding t) {
    switch (t) {
    case E_NULL:            return "NULL";
    case E_EXTERNAL:        return "EXTERNAL";
    case E_GOLOMB:          return "GOLOMB";
    case E_HUFFMAN:         return "HUFFMAN";
    case E_BYTE_ARRAY_LEN:  return "BYTE_ARRAY_LEN";
    case E_BYTE_ARRAY_STOP: return "BYTE_ARRAY_STOP";
    case E_BETA:            return "BETA";
    case E_SUBEXP:          return "SUBEXP";
    case E_GOLOMB_RICE:     return "GOLOMB_RICE";
    case E_GAMMA:           return "GAMMA";

    case E_VARINT_UNSIGNED: return "VARINT_UNSIGNED";
    case E_VARINT_SIGNED:   return "VARINT_SIGNED";
    case E_CONST_BYTE:      return "CONST_BYTE";
    case E_CONST_INT:       return "CONST_INT";

    case E_NUM_CODECS:
    default:                return "?";
    }
}

static cram_codec *(*decode_init[])(cram_block_compression_hdr *hdr,
                                    char *data,
                                    int size,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    int version, varint_vec *vv) = {
    // CRAM 3.0 valid codecs
    NULL, // null codec
    cram_external_decode_init,
    NULL, // golomb
    cram_huffman_decode_init,
    cram_byte_array_len_decode_init,
    cram_byte_array_stop_decode_init,
    cram_beta_decode_init,
    cram_subexp_decode_init,
    NULL, // golomb rice
    cram_gamma_decode_init,

    // Gap between CRAM 3 and CRAM 4; 9 to 39 inclusive
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,

    NULL,                      // was xbyte
    cram_varint_decode_init,   // varint unsigned
    cram_varint_decode_init,   // varint signed
    cram_const_decode_init,    // const byte
    cram_const_decode_init,    // const int

    // Gap to CRAM 4 transfomrations; 45 to 49 inclusive
    NULL, NULL, NULL, NULL, NULL,

    NULL, // xhuffman
    cram_xpack_decode_init,
    cram_xrle_decode_init,
    cram_xdelta_decode_init,
};

cram_codec *cram_decoder_init(cram_block_compression_hdr *hdr,
                              enum cram_encoding codec,
                              char *data, int size,
                              enum cram_external_type option,
                              int version, varint_vec *vv) {
    if (codec >= E_NULL && codec < E_NUM_CODECS && decode_init[codec]) {
        cram_codec *r = decode_init[codec](hdr, data, size, codec,
                                           option, version, vv);
        if (r) {
            r->vv = vv;
            r->codec_id = hdr->ncodecs++;
        }
        return r;
    } else {
        hts_log_error("Unimplemented codec of type %s", cram_encoding2str(codec));
        return NULL;
    }
}

static cram_codec *(*encode_init[])(cram_stats *stx,
                                    enum cram_encoding codec,
                                    enum cram_external_type option,
                                    void *opt,
                                    int version, varint_vec *vv) = {
    // CRAM 3.0 valid codecs
    NULL, // null codec
    cram_external_encode_init, // int/bytes in cram 3, byte only in cram 4
    NULL, // golomb
    cram_huffman_encode_init,
    cram_byte_array_len_encode_init,
    cram_byte_array_stop_encode_init,
    cram_beta_encode_init,
    NULL, // subexponential (we support decode only)
    NULL, // golomb rice
    NULL, // gamma (we support decode only)

    // Gap between CRAM 3 and CRAM 4; 9 to 39 inclusive
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,

    NULL, // was xbyte
    cram_varint_encode_init, // varint unsigned
    cram_varint_encode_init, // varint signed
    cram_const_encode_init,  // const byte
    cram_const_encode_init,  // const int

    // Gap to CRAM 4 transfomrations; 45 to 49 inclusive
    NULL, NULL, NULL, NULL, NULL,

    NULL, // xhuffman
    cram_xpack_encode_init,
    cram_xrle_encode_init,
    cram_xdelta_encode_init,
};

cram_codec *cram_encoder_init(enum cram_encoding codec,
                              cram_stats *st,
                              enum cram_external_type option,
                              void *dat,
                              int version, varint_vec *vv) {
    if (st && !st->nvals)
        return NULL;

    // cram_stats_encoding assumes integer data, but if option
    // is E_BYTE then tweak the requested encoding.  This ought
    // to be fixed in cram_stats_encoding instead.
    if (option == E_BYTE || option == E_BYTE_ARRAY ||
       option == E_BYTE_ARRAY_BLOCK) {
       if (codec == E_VARINT_SIGNED || codec == E_VARINT_UNSIGNED)
           codec = E_EXTERNAL;
       else if (codec == E_CONST_INT)
           codec = E_CONST_BYTE;
    }

    if (encode_init[codec]) {
        cram_codec *r;
        if ((r = encode_init[codec](st, codec, option, dat, version, vv)))
            r->out = NULL;
        if (!r) {
            hts_log_error("Unable to initialise codec of type %s", cram_encoding2str(codec));
            return NULL;
        }
        r->vv = vv;
        return r;
    } else {
        hts_log_error("Unimplemented codec of type %s", cram_encoding2str(codec));
        abort();
    }
}

/*
 * Returns the content_id used by this codec, also in id2 if byte_array_len.
 * Returns -1 for the CORE block and -2 for unneeded.
 * id2 is only filled out for BYTE_ARRAY_LEN which uses 2 codecs.
 */
int cram_codec_to_id(cram_codec *c, int *id2) {
    int bnum1, bnum2 = -2;

    switch (c->codec) {
    case E_CONST_INT:
    case E_CONST_BYTE:
        bnum1 = -2; // no blocks used
        break;

    case E_HUFFMAN:
        bnum1 = c->u.huffman.ncodes == 1 ? -2 : -1;
        break;

    case E_GOLOMB:
    case E_BETA:
    case E_SUBEXP:
    case E_GOLOMB_RICE:
    case E_GAMMA:
        // CORE block
        bnum1 = -1;
        break;

    case E_EXTERNAL:
    case E_VARINT_UNSIGNED:
    case E_VARINT_SIGNED:
        bnum1 = c->u.external.content_id;
        break;

    case E_BYTE_ARRAY_LEN:
        bnum1 = cram_codec_to_id(c->u.byte_array_len.len_codec, NULL);
        bnum2 = cram_codec_to_id(c->u.byte_array_len.val_codec, NULL);
        break;

    case E_BYTE_ARRAY_STOP:
        bnum1 = c->u.byte_array_stop.content_id;
        break;

    case E_NULL:
        bnum1 = -2;
        break;

    default:
        hts_log_error("Unknown codec type %d", c->codec);
        bnum1 = -1;
    }

    if (id2)
        *id2 = bnum2;
    return bnum1;
}


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
int cram_codec_decoder2encoder(cram_fd *fd, cram_codec *c) {
    int j;

    switch (c->codec) {
    case E_CONST_INT:
    case E_CONST_BYTE:
        // shares struct with decode
        c->store = cram_const_encode_store;
        break;

    case E_EXTERNAL:
        // shares struct with decode
        c->free = cram_external_encode_free;
        c->store = cram_external_encode_store;
        if (c->decode == cram_external_decode_int)
            c->encode = cram_external_encode_int;
        else if (c->decode == cram_external_decode_long)
            c->encode = cram_external_encode_long;
        else if (c->decode == cram_external_decode_char)
            c->encode = cram_external_encode_char;
        else if (c->decode == cram_external_decode_block)
            c->encode = cram_external_encode_char;
        else
            return -1;
        break;

    case E_VARINT_SIGNED:
    case E_VARINT_UNSIGNED:
        // shares struct with decode
        c->free = cram_varint_encode_free;
        c->store = cram_varint_encode_store;
        if (c->decode == cram_varint_decode_int)
            c->encode = cram_varint_encode_int;
        else if (c->decode == cram_varint_decode_sint)
            c->encode = cram_varint_encode_sint;
        else if (c->decode == cram_varint_decode_long)
            c->encode = cram_varint_encode_long;
        else if (c->decode == cram_varint_decode_slong)
            c->encode = cram_varint_encode_slong;
        else
            return -1;
        break;

    case E_HUFFMAN: {
        // New structure, so switch.
        // FIXME: we huffman and e_huffman structs amended, we could
        // unify this.
        cram_codec *t = malloc(sizeof(*t));
        if (!t) return -1;
        t->vv     = c->vv;
        t->codec = E_HUFFMAN;
        t->free = cram_huffman_encode_free;
        t->store = cram_huffman_encode_store;
        t->u.e_huffman.codes = c->u.huffman.codes;
        t->u.e_huffman.nvals = c->u.huffman.ncodes;
        t->u.e_huffman.option = c->u.huffman.option;
        for (j = 0; j < t->u.e_huffman.nvals; j++) {
            int32_t sym = t->u.e_huffman.codes[j].symbol;
            if (sym >= -1 && sym < MAX_HUFF)
                t->u.e_huffman.val2code[sym+1] = j;
        }

        if (c->decode == cram_huffman_decode_char0)
            t->encode = cram_huffman_encode_char0;
        else if (c->decode == cram_huffman_decode_char)
            t->encode = cram_huffman_encode_char;
        else if (c->decode == cram_huffman_decode_int0)
            t->encode = cram_huffman_encode_int0;
        else if (c->decode == cram_huffman_decode_int)
            t->encode = cram_huffman_encode_int;
        else if (c->decode == cram_huffman_decode_long0)
            t->encode = cram_huffman_encode_long0;
        else if (c->decode == cram_huffman_decode_long)
            t->encode = cram_huffman_encode_long;
        else {
            free(t);
            return -1;
        }
        *c = *t;
        free(t);
        break;
    }

    case E_BETA:
        // shares struct with decode
        c->free = cram_beta_encode_free;
        c->store = cram_beta_encode_store;
        if (c->decode == cram_beta_decode_int)
            c->encode = cram_beta_encode_int;
        else if (c->decode == cram_beta_decode_long)
            c->encode = cram_beta_encode_long;
        else if (c->decode == cram_beta_decode_char)
            c->encode = cram_beta_encode_char;
        else
            return -1;
        break;

    case E_XPACK: {
        // shares struct with decode
        cram_codec t = *c;
        t.free = cram_xpack_encode_free;
        t.store = cram_xpack_encode_store;
        if (t.decode == cram_xpack_decode_long)
            t.encode = cram_xpack_encode_long;
        else if (t.decode == cram_xpack_decode_int)
            t.encode = cram_xpack_encode_int;
        else if (t.decode == cram_xpack_decode_char)
            t.encode = cram_xpack_encode_char;
        else
            return -1;
        t.u.e_xpack.sub_codec = t.u.xpack.sub_codec;
        if (cram_codec_decoder2encoder(fd, t.u.e_xpack.sub_codec) == -1)
            return -1;
        *c = t;
        break;
    }

    case E_BYTE_ARRAY_LEN: {
        cram_codec *t = malloc(sizeof(*t));
        if (!t) return -1;
        t->vv     = c->vv;
        t->codec  = E_BYTE_ARRAY_LEN;
        t->free   = cram_byte_array_len_encode_free;
        t->store  = cram_byte_array_len_encode_store;
        t->encode = cram_byte_array_len_encode;
        t->u.e_byte_array_len.len_codec = c->u.byte_array_len.len_codec;
        t->u.e_byte_array_len.val_codec = c->u.byte_array_len.val_codec;
        if (cram_codec_decoder2encoder(fd, t->u.e_byte_array_len.len_codec) == -1 ||
            cram_codec_decoder2encoder(fd, t->u.e_byte_array_len.val_codec) == -1) {
            t->free(t);
            return -1;
        }

        // {len,val}_{encoding,dat} are undefined, but unused.
        // Leaving them unset here means we can test that assertion.
        *c = *t;
        free(t);
        break;
    }

    case E_BYTE_ARRAY_STOP:
        // shares struct with decode
        c->free   = cram_byte_array_stop_encode_free;
        c->store  = cram_byte_array_stop_encode_store;
        c->encode = cram_byte_array_stop_encode;
        break;

    default:
        return -1;
    }

    return 0;
}

int cram_codec_describe(cram_codec *c, kstring_t *ks) {
    if (c && c->describe)
        return c->describe(c, ks);
    else
        return ksprintf(ks, "?");
}
