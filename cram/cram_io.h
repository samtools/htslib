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

/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the low level CRAM I/O primitives.
 * This includes basic data types such as byte, int, ITF-8,
 * maps, bitwise I/O, etc.
 */

#ifndef CRAM_IO_H
#define CRAM_IO_H

#include <stdint.h>

#include "misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@{ ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
 * Also see the itf8_get and itf8_put macros.
 */

/*! INTERNAL: Converts two characters into an integer for use in switch{} */
#define CRAM_KEY(a,b) ((((unsigned char) a)<<8)|(((unsigned char) b)))

/*! Reads an integer in ITF-8 encoding from 'fd' and stores it in
 * *val.
 *
 * @return
 * Returns the number of bytes read on success;
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val);

extern const int itf8_bytes[16];
extern const int ltf8_bytes[256];

/*! Pushes a value in ITF8 format onto the end of a block.
 *
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * @return
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int32_t val);
int ltf8_put_blk(cram_block *blk, int64_t val);

/*! Pulls a literal 32-bit value from a block.
 *
 * @returns the number of bytes decoded;
 *         -1 on failure.
 */
int int32_get_blk(cram_block *b, int32_t *val);

/*! Pushes a literal 32-bit value onto the end of a block.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure.
 */
int int32_put_blk(cram_block *blk, int32_t val);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*! Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * @return
 * Returns block pointer on success;
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
                           int content_id);

/*! Reads a block from a cram file.
 *
 * @return
 * Returns cram_block pointer on success;
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd);

/*! Writes a CRAM block.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b);

/*! Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b);

/*! Uncompress a memory block using Zlib.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size);

/*! Uncompresses a CRAM block, if compressed.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_uncompress_block(cram_block *b);

/*! Compresses a block.
 *
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
                        int method, int level);
int cram_compress_block2(cram_fd *fd, cram_slice *s,
                         cram_block *b, cram_metrics *metrics,
                         int method, int level);

cram_metrics *cram_new_metrics(void);
char *cram_block_method2str(enum cram_block_method_int m);
char *cram_content_type2str(enum cram_content_type t);

/*
 * Find an external block by its content_id
 */

static inline cram_block *cram_get_block_by_id(cram_slice *slice, int id) {
  //fprintf(stderr, "%d\t%p\n", id, slice->block_by_id);
    uint32_t v = id;
    if (slice->block_by_id && v < 256) {
        return slice->block_by_id[v];
    } else {
        v = 256 + v % 251;
        if (slice->block_by_id &&
            slice->block_by_id[v] &&
            slice->block_by_id[v]->content_id == id)
            return slice->block_by_id[v];

        // Otherwise a linear search in case of collision
        int i;
        for (i = 0; i < slice->hdr->num_blocks; i++) {
            cram_block *b = slice->block[i];
            if (b && b->content_type == EXTERNAL && b->content_id == id)
                return b;
        }
    }
    return NULL;
}

/* --- Accessor macros for manipulating blocks on a byte by byte basis --- */

/* Block size and data pointer. */
#define BLOCK_SIZE(b) ((b)->byte)
#define BLOCK_DATA(b) ((b)->data)

/* Returns the address one past the end of the block */
#define BLOCK_END(b) (&(b)->data[(b)->byte])

/* Make block exactly 'l' bytes long */
static inline int block_resize_exact(cram_block *b, size_t len) {
    unsigned char *tmp = realloc(b->data, len);
    if (!tmp)
        return -1;
    b->alloc = len;
    b->data = tmp;
    return 0;
}

/* Request block to be at least 'l' bytes long */
static inline int block_resize(cram_block *b, size_t len) {
    if (b->alloc > len)
        return 0;

    size_t alloc = b->alloc+800;
    alloc = MAX(alloc + (alloc>>2), len);
    return block_resize_exact(b, alloc);
}


/* Ensure the block can hold at least another 'l' bytes */
static inline int block_grow(cram_block *b, size_t len) {
    return block_resize(b, BLOCK_SIZE(b) + len);
}

/* Append string 's' of length 'l'. */
static inline int block_append(cram_block *b, const void *s, size_t len) {
    if (block_grow(b, len) < 0)
        return -1;

    if (len) {
        memcpy(BLOCK_END(b), s, len);
        BLOCK_SIZE(b) += len;
    }

    return 0;
}

/* Append as single character 'c' */
static inline int block_append_char(cram_block *b, char c) {
    if (block_grow(b, 1) < 0)
        return -1;

    b->data[b->byte++] = c;
    return 0;
}

/* Append a single unsigned integer */
static inline unsigned char *append_uint32(unsigned char *cp, uint32_t i);
static inline int block_append_uint(cram_block *b, unsigned int i) {
    if (block_grow(b, 11) < 0)
        return -1;

    unsigned char *cp = &b->data[b->byte];
    b->byte += append_uint32(cp, i) - cp;
    return 0;
}

// Versions of above with built in goto block_err calls.
#define BLOCK_RESIZE_EXACT(b,l) if (block_resize_exact((b),(l))<0) goto block_err
#define BLOCK_RESIZE(b,l)       if (block_resize((b),(l))      <0) goto block_err
#define BLOCK_GROW(b,l)         if (block_grow((b),(l))        <0) goto block_err
#define BLOCK_APPEND(b,s,l)     if (block_append((b),(s),(l))  <0) goto block_err
#define BLOCK_APPEND_CHAR(b,c)  if (block_append_char((b),(c)) <0) goto block_err
#define BLOCK_APPEND_UINT(b,i)  if (block_append_uint((b),(i)) <0) goto block_err

static inline unsigned char *append_uint32(unsigned char *cp, uint32_t i) {
    uint32_t j;

    if (i == 0) {
        *cp++ = '0';
        return cp;
    }

    if (i < 100)        goto b1;
    if (i < 10000)      goto b3;
    if (i < 1000000)    goto b5;
    if (i < 100000000)  goto b7;

    if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
    if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7:if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
    if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5:if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
    if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3:if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
    if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1:if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
    if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

static inline unsigned char *append_sub32(unsigned char *cp, uint32_t i) {
    *cp++ = i / 100000000 + '0', i %= 100000000;
    *cp++ = i / 10000000  + '0', i %= 10000000;
    *cp++ = i / 1000000   + '0', i %= 1000000;
    *cp++ = i / 100000    + '0', i %= 100000;
    *cp++ = i / 10000     + '0', i %= 10000;
    *cp++ = i / 1000      + '0', i %= 1000;
    *cp++ = i / 100       + '0', i %= 100;
    *cp++ = i / 10        + '0', i %= 10;
    *cp++ = i             + '0';

    return cp;
}

static inline unsigned char *append_uint64(unsigned char *cp, uint64_t i) {
    uint64_t j;

    if (i <= 0xffffffff)
        return append_uint32(cp, i);

    if ((j = i/1000000000) > 1000000000) {
        cp = append_uint32(cp, j/1000000000);
        j %= 1000000000;
        cp = append_sub32(cp, j);
    } else {
        cp = append_uint32(cp, i / 1000000000);
    }
    cp = append_sub32(cp, i % 1000000000);

    return cp;
}

#define BLOCK_UPLEN(b) \
    (b)->comp_size = (b)->uncomp_size = BLOCK_SIZE((b))

/**@}*/
/**@{ ----------------------------------------------------------------------
 * Reference sequence handling
 */

/*! Loads a reference set from fn and stores in the cram_fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_load_reference(cram_fd *fd, char *fn);

/*! Generates a lookup table in refs based on the SQ headers in sam_hdr_t.
 *
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int refs2id(refs_t *r, sam_hdr_t *hdr);

void refs_free(refs_t *r);

/*! Returns a portion of a reference sequence from start to end inclusive.
 *
 * The returned pointer is owned by the cram_file fd and should not be freed
 * by the caller. It is valid only until the next cram_get_ref is called
 * with the same fd parameter (so is thread-safe if given multiple files).
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * @return
 * Returns reference on success;
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, hts_pos_t start, hts_pos_t end);
void cram_ref_incr(refs_t *r, int id);
void cram_ref_decr(refs_t *r, int id);
/**@}*/
/**@{ ----------------------------------------------------------------------
 * Containers
 */

/*! Creates a new container, specifying the maximum number of slices
 * and records permitted.
 *
 * @return
 * Returns cram_container ptr on success;
 *         NULL on failure
 */
cram_container *cram_new_container(int nrec, int nslice);
void cram_free_container(cram_container *c);

/*! Reads a container header.
 *
 * @return
 * Returns cram_container on success;
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd);

/*! Writes a container structure.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_container(cram_fd *fd, cram_container *h);

/*! Flushes a container to disk.
 *
 * Flushes a completely or partially full container to disk, writing
 * container structure, header and blocks. This also calls the encoder
 * functions.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_flush_container(cram_fd *fd, cram_container *c);
int cram_flush_container_mt(cram_fd *fd, cram_container *c);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * Compression headers; the first part of the container
 */

/*! Creates a new blank container compression header
 *
 * @return
 * Returns header ptr on success;
 *         NULL on failure
 */
cram_block_compression_hdr *cram_new_compression_header(void);

/*! Frees a cram_block_compression_hdr */
void cram_free_compression_header(cram_block_compression_hdr *hdr);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * Slices and slice headers
 */

/*! Frees a slice header */
void cram_free_slice_header(cram_block_slice_hdr *hdr);

/*! Frees a slice */
void cram_free_slice(cram_slice *s);

/*! Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * @return
 * Returns cram_slice ptr on success;
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs);

/*! Loads an entire slice.
 *
 * FIXME: In 1.0 the native unit of slices within CRAM is broken
 * as slices contain references to objects in other slices.
 * To work around this while keeping the slice oriented outer loop
 * we read all slices and stitch them together into a fake large
 * slice instead.
 *
 * @return
 * Returns cram_slice ptr on success;
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd);



/**@}*/
/**@{ ----------------------------------------------------------------------
 * CRAM file definition (header)
 */

/*! Reads a CRAM file definition structure.
 *
 * @return
 * Returns file_def ptr on success;
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd);

/*! Writes a cram_file_def structure to cram_fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def);

/*! Frees a cram_file_def structure. */
void cram_free_file_def(cram_file_def *def);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * SAM header I/O
 */

/*! Reads the SAM header from the first CRAM data block.
 *
 * Also performs minimal parsing to extract read-group
 * and sample information.
 *
 * @return
 * Returns SAM hdr ptr on success;
 *         NULL on failure
 */
sam_hdr_t *cram_read_SAM_hdr(cram_fd *fd);

/*! Writes a CRAM SAM header.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, sam_hdr_t *hdr);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*! Opens a CRAM file for read (mode "rb") or write ("wb").
 *
 * The filename may be "-" to indicate stdin or stdout.
 *
 * @return
 * Returns file handle on success;
 *         NULL on failure.
 */
cram_fd *cram_open(const char *filename, const char *mode);

/*! Opens an existing stream for reading or writing.
 *
 * @return
 * Returns file handle on success;
 *         NULL on failure.
 */
cram_fd *cram_dopen(struct hFILE *fp, const char *filename, const char *mode);

/*! Closes a CRAM file.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_close(cram_fd *fd);

/*
 * Seek within a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence);

/*
 * Flushes a CRAM file.
 * Useful for when writing to stdout without wishing to close the stream.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush(cram_fd *fd);

/*! Checks for end of file on a cram_fd stream.
 *
 * @return
 * Returns 0 if not at end of file
 *         1 if we hit an expected EOF (end of range or EOF block)
 *         2 for other EOF (end of stream without EOF block)
 */
int cram_eof(cram_fd *fd);

/*! Sets options on the cram_fd.
 *
 * See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum hts_fmt_option opt, ...);

/*! Sets options on the cram_fd.
 *
 * See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_voption(cram_fd *fd, enum hts_fmt_option opt, va_list args);

/*!
 * Attaches a header to a cram_fd.
 *
 * This should be used when creating a new cram_fd for writing where
 * we have an sam_hdr_t already constructed (eg from a file we've read
 * in).
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_header2(cram_fd *fd, const sam_hdr_t *hdr);

/*!
 * Returns the hFILE connected to a cram_fd.
 */
static inline struct hFILE *cram_hfile(cram_fd *fd) {
    return fd->fp;
}

#ifdef __cplusplus
}
#endif

#endif /* CRAM_IO_H */
