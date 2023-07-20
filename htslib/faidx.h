/// @file htslib/faidx.h
/// FASTA random access.
/*
   Copyright (C) 2008, 2009, 2013, 2014, 2016, 2017-2020, 2022-2023 Genome Research Ltd.

   Author: Heng Li <lh3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef HTSLIB_FAIDX_H
#define HTSLIB_FAIDX_H

#include <stdint.h>
#include "hts_defs.h"
#include "hts.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @file

  Index FASTA or FASTQ files and extract subsequence.

  The fai file index columns for FASTA are:
    - chromosome name
    - chromosome length: number of bases
    - offset: number of bytes to skip to get to the first base
        from the beginning of the file, including the length
        of the sequence description string (`>chr ..\n`)
    - line length: number of bases per line (excluding `\n`)
    - binary line length: number of bytes, including `\n`

   The index for FASTQ is similar to above:
    - chromosome name
    - chromosome length: number of bases
    - sequence offset: number of bytes to skip to get to the first base
        from the beginning of the file, including the length
        of the sequence description string (`@chr ..\n`)
    - line length: number of bases per line (excluding `\n`)
    - binary line length: number of bytes, including `\n`
    - quality offset: number of bytes to skip from the beginning of the file
        to get to the first quality value in the indexed entry.

    The FASTQ version of the index uses line length and binary line length
    for both the sequence and the quality values, so they must be line
    wrapped in the same way.
 */

struct faidx_t;
/// Opaque structure representing FASTA index
typedef struct faidx_t faidx_t;

/// Opaque structure; sole item needed from htslib/thread_pool.h
struct hts_tpool;

/// File format to be dealing with.
enum fai_format_options {
    FAI_NONE,
    FAI_FASTA,
    FAI_FASTQ
};

/// Build index for a FASTA or FASTQ or bgzip-compressed FASTA or FASTQ file.
/**  @param  fn  FASTA/FASTQ file name
     @param  fnfai Name of .fai file to build.
     @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
     @return     0 on success; or -1 on failure

If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
file will only be built if fn is bgzip-compressed.
*/
HTSLIB_EXPORT
int fai_build3(const char *fn, const char *fnfai, const char *fngzi) HTS_RESULT_USED;

/// Build index for a FASTA or FASTQ or bgzip-compressed FASTA or FASTQ file.
/** @param  fn  FASTA/FASTQ file name
    @return     0 on success; or -1 on failure

File "fn.fai" will be generated.  This function is equivalent to
fai_build3(fn, NULL, NULL);
*/
HTSLIB_EXPORT
int fai_build(const char *fn) HTS_RESULT_USED;

/// Destroy a faidx_t struct
HTSLIB_EXPORT
void fai_destroy(faidx_t *fai);

enum fai_load_options {
    FAI_CREATE = 0x01,
};

/// Load FASTA indexes.
/** @param  fn  File name of the FASTA file (can be compressed with bgzip).
    @param  fnfai File name of the FASTA index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @return Pointer to a faidx_t struct on success, NULL on failure.

If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
The bgzip index is only needed if fn is compressed.

If (flags & FAI_CREATE) is true, the index files will be built using
fai_build3() if they are not already present.

The struct returned by a successful call should be freed via fai_destroy()
when it is no longer needed.
*/
HTSLIB_EXPORT
faidx_t *fai_load3(const char *fn, const char *fnfai, const char *fngzi,
                   int flags);

/// Load index from "fn.fai".
/** @param  fn  File name of the FASTA file
    @return Pointer to a faidx_t struct on success, NULL on failure.

This function is equivalent to fai_load3(fn, NULL, NULL, FAI_CREATE|FAI_CACHE);
*/
HTSLIB_EXPORT
faidx_t *fai_load(const char *fn);

/// Load FASTA or FASTQ indexes.
/** @param  fn  File name of the FASTA/FASTQ file (can be compressed with bgzip).
    @param  fnfai File name of the FASTA/FASTQ index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @param  format FASTA or FASTQ file format
    @return Pointer to a faidx_t struct on success, NULL on failure.

If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
The bgzip index is only needed if fn is compressed.

If (flags & FAI_CREATE) is true, the index files will be built using
fai_build3() if they are not already present.

The struct returned by a successful call should be freed via fai_destroy()
when it is no longer needed.
*/
HTSLIB_EXPORT
faidx_t *fai_load3_format(const char *fn, const char *fnfai, const char *fngzi,
                   int flags, enum fai_format_options format);

/// Load index from "fn.fai".
/** @param  fn  File name of the FASTA/FASTQ file
    @param  format FASTA or FASTQ file format
    @return Pointer to a faidx_t struct on success, NULL on failure.

This function is equivalent to fai_load3_format(fn, NULL, NULL, FAI_CREATE|FAI_CACHE, format);
*/
HTSLIB_EXPORT
faidx_t *fai_load_format(const char *fn, enum fai_format_options format);

/// Fetch the sequence in a region
/** @param  fai  Pointer to the faidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @param  len  Length of the region; -2 if seq not present, -1 general error
    @return      Pointer to the sequence; `NULL` on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.

To work around ambiguous parsing issues, eg both "chr1" and "chr1:100-200"
are reference names, quote using curly braces.
Thus "{chr1}:100-200" and "{chr1:100-200}" disambiguate the above example.
*/
HTSLIB_EXPORT
char *fai_fetch(const faidx_t *fai, const char *reg, int *len);
HTSLIB_EXPORT
char *fai_fetch64(const faidx_t *fai, const char *reg, hts_pos_t *len);

/// Query the line-wrap length for a chromosome specified as part of a region
/** @param  fai  Pointer to the faidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @return      The line length (excluding newline),
                 negative on error.
*/
HTSLIB_EXPORT
hts_pos_t fai_line_length(const faidx_t *fai, const char *reg);

/// Fetch the quality string for a region for FASTQ files
/** @param  fai  Pointer to the faidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @param  len  Length of the region; -2 if seq not present, -1 general error
    @return      Pointer to the quality string; null on failure

The returned quality string is allocated by `malloc()` family and should be
destroyed by end users by calling `free()` on it.

Region names can be quoted with curly braces, as for fai_fetch().
*/
HTSLIB_EXPORT
char *fai_fetchqual(const faidx_t *fai, const char *reg, int *len);
HTSLIB_EXPORT
char *fai_fetchqual64(const faidx_t *fai, const char *reg, hts_pos_t *len);

/// Fetch the number of sequences
/** @param  fai  Pointer to the faidx_t struct
    @return      The number of sequences
*/
HTSLIB_EXPORT
int faidx_fetch_nseq(const faidx_t *fai) HTS_DEPRECATED("Please use faidx_nseq instead");

/// Fetch the sequence in a region
/** @param  fai  Pointer to the faidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
HTSLIB_EXPORT
char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);

/// Fetch the sequence in a region
/** @param  fai  Pointer to the faidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
HTSLIB_EXPORT
char *faidx_fetch_seq64(const faidx_t *fai, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/// Fetch the quality string in a region for FASTQ files
/** @param  fai  Pointer to the faidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
HTSLIB_EXPORT
char *faidx_fetch_qual(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);

/// Fetch the quality string in a region for FASTQ files
/** @param  fai  Pointer to the faidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
HTSLIB_EXPORT
char *faidx_fetch_qual64(const faidx_t *fai, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/// Query if sequence is present
/**   @param  fai  Pointer to the faidx_t struct
      @param  seq  Sequence name
      @return      1 if present or 0 if absent
*/
HTSLIB_EXPORT
int faidx_has_seq(const faidx_t *fai, const char *seq);

/// Return number of sequences in fai index
HTSLIB_EXPORT
int faidx_nseq(const faidx_t *fai);

/// Return name of i-th sequence
HTSLIB_EXPORT
const char *faidx_iseq(const faidx_t *fai, int i);

/// Return sequence length
/** @param  fai  Pointer to the faidx_t struct
    @param  seq  Name of the sequence
    @return Sequence length, or -1 if not present
*/
HTSLIB_EXPORT
hts_pos_t faidx_seq_len64(const faidx_t *fai, const char *seq);

/// Return sequence length
/** @param  fai  Pointer to the faidx_t struct
    @param  seq  Name of the sequence
    @return Sequence length, or -1 if not present

    @deprecated This funtion cannot handle very long sequences.
                Use faidx_seq_len64() instead.
*/
HTSLIB_EXPORT
int faidx_seq_len(const faidx_t *fai, const char *seq);

/// Parses a region string.
/** @param  fai   Pointer to the faidx_t struct
    @param  s     Region string
    @param  tid   Returns which i-th sequence is described in the region.
    @param  beg   Returns the start of the region (0 based)
    @param  end   Returns the one past last of the region (0 based)
    @param  flags Parsing method, see HTS_PARSE_* in hts.h.
    @return       Pointer to end of parsed s if successful, NULL if not.

    To work around ambiguous parsing issues, eg both "chr1" and "chr1:100-200"
    are reference names, quote using curly braces.
    Thus "{chr1}:100-200" and "{chr1:100-200}" disambiguate the above example.
*/
HTSLIB_EXPORT
const char *fai_parse_region(const faidx_t *fai, const char *s,
                             int *tid, hts_pos_t *beg, hts_pos_t *end,
                             int flags);

/// Adjust region to the actual sequence length
/** @param  fai   Pointer to the faidx_t struct
    @param  tid   Sequence index, as returned by fai_parse_region()
    @param  beg[in,out]   The start of the region (0 based)
    @param  end[in,out]   One past end of the region (0 based)
    @return 1, 2, or 3 if @p beg, @p end, or both are adjusted,
            0 if @p beg and @p end are unchanged
            -1 on error

    Looks up the length of @p tid, and then adjusts the values of @p beg
    and @p end if they fall outside the boundaries of the sequence.

    If @p beg > @p end, it will be set to @p end.

    The return value indicates which, if any, of the inputs have been
    adjusted.  -1 will be returned if @p tid is not a valid sequence index.
*/
HTSLIB_EXPORT
int fai_adjust_region(const faidx_t *fai, int tid,
                      hts_pos_t *beg, hts_pos_t *end);

/// Sets the cache size of the underlying BGZF compressed file
/** @param  fai         Pointer to the faidx_t struct
 *  @param  cache_size  Selected cache size in bytes
 */
HTSLIB_EXPORT
void fai_set_cache_size(faidx_t *fai, int cache_size);

/// Adds a thread pool to the underlying BGZF layer.
/** @param fai         FAI file handler
 *  @param pool        The thread pool (see hts_create_threads)
 *  @param qsize       The size of the job queue.  If 0 this is twice the
 *                     number of threads in the pool.
 */
HTSLIB_EXPORT
int fai_thread_pool(faidx_t *fai, struct hts_tpool *pool, int qsize);

/// Determines the path to the reference index file
/** @param  fa    String with the path to the reference file
 *  @return       String with the path to the reference index file, or NULL on failure

    If the reference path has the format reference.fa##idx##index.fa.fai,
    the index path is taken directly from it as index.fa.fai.
    If the reference file is local and the index file cannot be found, it
    will be created alongside the reference file.
    If the reference file is remote and the index file cannot be found,
    the method returns NULL.

    The returned string has to be freed by the user at the end of its scope.
 */
HTSLIB_EXPORT
char *fai_path(const char *fa);
#ifdef __cplusplus
}
#endif

#endif
