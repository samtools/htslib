/// @file htslib/bgzf.h
/// Low-level routines for direct BGZF operations.
/*
   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013, 2014, 2017, 2018-2019, 2022-2024 Genome Research Ltd

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

#ifndef HTSLIB_BGZF_H
#define HTSLIB_BGZF_H

#include <stdint.h>
#include <string.h>
#include <sys/types.h>

#include "hts_defs.h"

// Ensure ssize_t exists within this header. All #includes must precede this,
// and ssize_t must be undefined again at the end of this header.
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined ssize_t
#define HTSLIB_SSIZE_T
#define ssize_t intptr_t
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define BGZF_BLOCK_SIZE     0xff00 // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE 0x10000

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8
#define BGZF_ERR_MT     16 // stream cannot be multi-threaded
#define BGZF_ERR_CRC    32

struct hFILE;
struct hts_tpool;
struct kstring_t;
struct bgzf_mtaux_t;
typedef struct bgzidx_t bgzidx_t;
typedef struct bgzf_cache_t bgzf_cache_t;
struct z_stream_s;

struct BGZF {
    // Reserved bits should be written as 0; read as "don't care"
    unsigned errcode:16, reserved:1, is_write:1, no_eof_block:1, is_be:1;
    signed compress_level:9;
    unsigned last_block_eof:1, is_compressed:1, is_gzip:1;
    int cache_size;
    int block_length, block_clength, block_offset;
    int64_t block_address, uncompressed_address;
    void *uncompressed_block, *compressed_block;
    bgzf_cache_t *cache;
    struct hFILE *fp; // actual file handle
    struct bgzf_mtaux_t *mt; // only used for multi-threading
    bgzidx_t *idx;      // BGZF index
    int idx_build_otf;  // build index on the fly, set by bgzf_index_build_init()
    struct z_stream_s *gz_stream; // for gzip-compressed files
    int64_t seeked;     // virtual offset of last seek
};
#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif

    /******************
     * Basic routines *
     ******************/

    /**
     * Open an existing file descriptor for reading or writing.
     *
     * @param fd    file descriptor
     *              Note that the file must be opened in binary mode, or else
     *              there will be problems on platforms that make a difference
     *              between text and binary mode.
     * @param mode  mode matching /[rwag][u0-9]+/: 'r' for reading, 'w' for
     *              writing, 'a' for appending, 'g' for gzip rather than BGZF
     *              compression (with 'w' only), and digit specifies the zlib
     *              compression level.
     *              Note that there is a distinction between 'u' and '0': the
     *              first yields plain uncompressed output whereas the latter
     *              outputs uncompressed data wrapped in the zlib format.
     * @return      BGZF file handler; 0 on error
     */
    HTSLIB_EXPORT
    BGZF* bgzf_dopen(int fd, const char *mode);

    #define bgzf_fdopen(fd, mode) bgzf_dopen((fd), (mode)) // for backward compatibility

    /**
     * Open the specified file for reading or writing.
     */
    HTSLIB_EXPORT
    BGZF* bgzf_open(const char* path, const char *mode);

    /**
     * Open an existing hFILE stream for reading or writing.
     */
    HTSLIB_EXPORT
    BGZF* bgzf_hopen(struct hFILE *fp, const char *mode);

    /**
     * Close the BGZF and free all associated resources.
     *
     * @param fp    BGZF file handler
     * @return      0 on success and -1 on error
     */
    HTSLIB_EXPORT
    int bgzf_close(BGZF *fp);

    /**
     * Read up to _length_ bytes from the file storing into _data_.
     *
     * @param fp     BGZF file handler
     * @param data   data array to read into
     * @param length size of data to read
     * @return       number of bytes actually read; 0 on end-of-file and -1 on error
     */
    HTSLIB_EXPORT
    ssize_t bgzf_read(BGZF *fp, void *data, size_t length) HTS_RESULT_USED;

/**
 * bgzf_read optimised for small quantities, as a static inline
 * See bgzf_read() normal function for return values.
 */
static inline ssize_t bgzf_read_small(BGZF *fp, void *data, size_t length) {
    // A block length of 0 implies current block isn't loaded (see
    // bgzf_seek_common).  That gives negative available so careful on types
    if ((ssize_t)length < fp->block_length - fp->block_offset) {
        // Short cut the common and easy mode
        memcpy((uint8_t *)data,
               (uint8_t *)fp->uncompressed_block + fp->block_offset,
               length);
        fp->block_offset += length;
        fp->uncompressed_address += length;
        return length;
    } else {
        return bgzf_read(fp, data, length);
    }
}

    /**
     * Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
     * the complete _length_ bytes will be written (or queued for writing).
     *
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length size of data to write
     * @return       number of bytes written (i.e., _length_); negative on error
     */
    HTSLIB_EXPORT
    ssize_t bgzf_write(BGZF *fp, const void *data, size_t length) HTS_RESULT_USED;

/**
 * bgzf_write optimised for small quantities, as a static inline
 * See bgzf_write() normal function for return values.
 */
static inline
ssize_t bgzf_write_small(BGZF *fp, const void *data, size_t length) {
    if (fp->is_compressed
        && (size_t) (BGZF_BLOCK_SIZE - fp->block_offset) > length) {
        // Short cut the common and easy mode
        memcpy((uint8_t *)fp->uncompressed_block + fp->block_offset,
               data, length);
        fp->block_offset += length;
        return length;
    } else {
        return bgzf_write(fp, data, length);
    }
}

    /**
     * Write _length_ bytes from _data_ to the file, the index will be used to
     * decide the amount of uncompressed data to be written to each bgzip block.
     * If no I/O errors occur, the complete _length_ bytes will be written (or
     * queued for writing).
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length size of data to write
     * @return       number of bytes written (i.e., _length_); negative on error
     */
    HTSLIB_EXPORT
    ssize_t bgzf_block_write(BGZF *fp, const void *data, size_t length);

    /**
     * Returns the next byte in the file without consuming it.
     * @param fp     BGZF file handler
     * @return       -1 on EOF,
     *               -2 on error,
     *               otherwise the unsigned byte value.
     */
    HTSLIB_EXPORT
    int bgzf_peek(BGZF *fp);

    /**
     * Read up to _length_ bytes directly from the underlying stream without
     * decompressing.  Bypasses BGZF blocking, so must be used with care in
     * specialised circumstances only.
     *
     * @param fp     BGZF file handler
     * @param data   data array to read into
     * @param length number of raw bytes to read
     * @return       number of bytes actually read; 0 on end-of-file and -1 on error
     */
    HTSLIB_EXPORT
    ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length) HTS_RESULT_USED;

    /**
     * Write _length_ bytes directly to the underlying stream without
     * compressing.  Bypasses BGZF blocking, so must be used with care
     * in specialised circumstances only.
     *
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length number of raw bytes to write
     * @return       number of bytes actually written; -1 on error
     */
    HTSLIB_EXPORT
    ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length) HTS_RESULT_USED;

    /**
     * Write the data in the buffer to the file.
     *
     * @param fp     BGZF file handle
     * @return       0 on success and -1 on error
     */
    HTSLIB_EXPORT
    int bgzf_flush(BGZF *fp) HTS_RESULT_USED;

    /**
     * Return a virtual file pointer to the current location in the file.
     * No interpretation of the value should be made, other than a subsequent
     * call to bgzf_seek can be used to position the file at the same point.
     * Return value is non-negative on success.
     */
    #define bgzf_tell(fp) (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF))

    /**
     * Set the file to read from the location specified by _pos_.
     *
     * @param fp     BGZF file handler
     * @param pos    virtual file offset returned by bgzf_tell()
     * @param whence must be SEEK_SET
     * @return       0 on success and -1 on error
     *
     * @note It is not permitted to seek on files open for writing,
     * or files compressed with gzip (as opposed to bgzip).
     */
    HTSLIB_EXPORT
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence) HTS_RESULT_USED;

    /**
     * Check if the BGZF end-of-file (EOF) marker is present
     *
     * @param fp    BGZF file handler opened for reading
     * @return      1 if the EOF marker is present and correct;
     *              2 if it can't be checked, e.g., because fp isn't seekable;
     *              0 if the EOF marker is absent;
     *              -1 (with errno set) on error
     */
    HTSLIB_EXPORT
    int bgzf_check_EOF(BGZF *fp);

    /** Return the file's compression format
     *
     * @param fp  BGZF file handle
     * @return    A small integer matching the corresponding
     *            `enum htsCompression` value:
     *   - 0 / `no_compression` if the file is uncompressed
     *   - 1 / `gzip` if the file is plain GZIP-compressed
     *   - 2 / `bgzf` if the file is BGZF-compressed
     * @since 1.4
     */
    HTSLIB_EXPORT
    int bgzf_compression(BGZF *fp);

    /**
     * Check if a file is in the BGZF format
     *
     * @param fn    file name
     * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
     */
    HTSLIB_EXPORT
    int bgzf_is_bgzf(const char *fn) HTS_DEPRECATED("Use bgzf_compression() or hts_detect_format() instead");

    /*********************
     * Advanced routines *
     *********************/

    /**
     * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
     *
     * @param fp    BGZF file handler
     * @param size  size of cache in bytes; 0 to disable caching (default)
     */
    HTSLIB_EXPORT
    void bgzf_set_cache_size(BGZF *fp, int size);

    /**
     * Flush the file if the remaining buffer size is smaller than _size_
     * @return      0 if flushing succeeded or was not needed; negative on error
     */
    HTSLIB_EXPORT
    int bgzf_flush_try(BGZF *fp, ssize_t size) HTS_RESULT_USED;

    /**
     * Read one byte from a BGZF file. It is faster than bgzf_read()
     * @param fp     BGZF file handler
     * @return       byte read; -1 on end-of-file; <= -2 on error
     */
    HTSLIB_EXPORT
    int bgzf_getc(BGZF *fp);

    /**
     * Read one line from a BGZF file. It is faster than bgzf_getc()
     *
     * @param fp     BGZF file handler
     * @param delim  delimiter
     * @param str    string to write to; must be initialized
     * @return       length of the string (capped at INT_MAX);
     *               -1 on end-of-file; <= -2 on error
     */
    HTSLIB_EXPORT
    int bgzf_getline(BGZF *fp, int delim, struct kstring_t *str);

    /**
     * Read the next BGZF block.
     */
    HTSLIB_EXPORT
    int bgzf_read_block(BGZF *fp) HTS_RESULT_USED;

    /**
     * Enable multi-threading via a shared thread pool.  This means
     * both encoder and decoder can balance usage across a single pool
     * of worker jobs.
     *
     * @param fp          BGZF file handler
     * @param pool        The thread pool (see hts_create_threads)
     * @param qsize       The size of the job queue.  If 0 this is twice the
     *                    number of threads in the pool.
     */
    HTSLIB_EXPORT
    int bgzf_thread_pool(BGZF *fp, struct hts_tpool *pool, int qsize);

    /**
     * Enable multi-threading
     *
     * @param fp          BGZF file handler
     * @param n_threads   #threads used for reading / writing
     * @param n_sub_blks  Unused (was #blocks processed by each thread)
     */
    HTSLIB_EXPORT
    int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);

    /**
     * Compress a single BGZF block.
     *
     * @param dst    output buffer (must have size >= BGZF_MAX_BLOCK_SIZE)
     * @param dlen   size of output buffer; updated on return to the number
     *               of bytes actually written to dst
     * @param src    buffer to be compressed
     * @param slen   size of data to compress (must be <= BGZF_BLOCK_SIZE)
     * @param level  compression level
     * @return       0 on success and negative on error
     */
    HTSLIB_EXPORT
    int bgzf_compress(void *dst, size_t *dlen, const void *src, size_t slen, int level);

    /*******************
     * bgzidx routines *
     *******************/

    /**
     *  Position BGZF at the uncompressed offset
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *  @param uoffset      file offset in the uncompressed data
     *  @param where        must be SEEK_SET
     *
     *  Returns 0 on success and -1 on error.
     *
     *  @note It is not permitted to seek on files open for writing,
     *  or files compressed with gzip (as opposed to bgzip).
     */
    HTSLIB_EXPORT
    int bgzf_useek(BGZF *fp, off_t uoffset, int where) HTS_RESULT_USED;

    /**
     *  Position in uncompressed BGZF
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *
     *  Returns the current offset on success and -1 on error.
     */
    HTSLIB_EXPORT
    off_t bgzf_utell(BGZF *fp);

    /**
     * Tell BGZF to build index while compressing.
     *
     * @param fp          BGZF file handler; can be opened for reading or writing.
     *
     * Returns 0 on success and -1 on error.
     *
     * @note This function must be called before any data has been read or
     * written, and in particular before calling bgzf_mt() on the same
     * file handle (as threads may start reading data before the index
     * has been set up).
     */
    HTSLIB_EXPORT
    int bgzf_index_build_init(BGZF *fp);

    /// Load BGZF index
    /**
     * @param fp          BGZF file handler
     * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     * @return 0 on success and -1 on error.
     */
    HTSLIB_EXPORT
    int bgzf_index_load(BGZF *fp,
                        const char *bname, const char *suffix) HTS_RESULT_USED;

    /// Load BGZF index from an hFILE
    /**
     * @param fp   BGZF file handle
     * @param idx  hFILE to read from
     * @param name file name (for error reporting only; can be NULL)
     * @return 0 on success and -1 on error.
     *
     * Populates @p fp with index data read from the hFILE handle @p idx.
     * The file pointer to @idx should point to the start of the index
     * data when this function is called.
     *
     * The file name can optionally be passed in the @p name parameter.  This
     * is only used for printing error messages; if NULL the word "index" is
     * used instead.
     */
    HTSLIB_EXPORT
    int bgzf_index_load_hfile(BGZF *fp, struct hFILE *idx,
                              const char *name) HTS_RESULT_USED;

    /// Save BGZF index
    /**
     * @param fp          BGZF file handler
     * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     * @return 0 on success and -1 on error.
     */
    HTSLIB_EXPORT
    int bgzf_index_dump(BGZF *fp,
                        const char *bname, const char *suffix) HTS_RESULT_USED;

    /// Write a BGZF index to an hFILE
    /**
     * @param fp     BGZF file handle
     * @param idx    hFILE to write to
     * @param name   file name (for error reporting only, can be NULL)
     * @return 0 on success and -1 on error.
     *
     * Write index data from @p fp to the file @p idx.
     *
     * The file name can optionally be passed in the @p name parameter.  This
     * is only used for printing error messages; if NULL the word "index" is
     * used instead.
     */

    HTSLIB_EXPORT
    int bgzf_index_dump_hfile(BGZF *fp, struct hFILE *idx,
                              const char *name) HTS_RESULT_USED;

#ifdef __cplusplus
}
#endif

#ifdef HTSLIB_SSIZE_T
#undef HTSLIB_SSIZE_T
#undef ssize_t
#endif

#endif
