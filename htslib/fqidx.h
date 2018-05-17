/// @file htslib/fqidx.h
/// FASTQ random access.
/*
   Copyright (C) 2008, 2009, 2013, 2014, 2016, 2017, 2018 Genome Research Ltd.

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

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @file

  Index FASTQ files and extract subsequence.

  The fqi file index columns are:
    - chromosome name
    - chromosome length: number of bases
    - offset: number of bytes to skip to get to the first base
        from the beginning of the file, including the length
        of the sequence description string (`>chr ..\n`)
    - line length: number of bases per line (excluding `\n`)
    - binary line length: number of bytes, including `\n`
 */

struct __fqidx_t;
/// Opaque structure representing FASTQ index
typedef struct __fqidx_t fqidx_t;

/// Build index for a FASTQ or bgzip-compressed FASTQ file.
/**  @param  fn  FASTQ file name
     @param  fnfqi Name of .fqi file to build.
     @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
     @return     0 on success; or -1 on failure

If fnfqi is NULL, ".fqi" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
file will only be built if fn is bgzip-compressed.
*/
int fqi_build3(const char *fn, const char *fnfqi, const char *fngzi) HTS_RESULT_USED;

/// Build index for a FASTQ or bgzip-compressed FASTQ file.
/** @param  fn  FASTQ file name
    @return     0 on success; or -1 on failure

File "fn.fqi" will be generated.  This function is equivalent to
fqi_build3(fn, NULL, NULL);
*/
int fqi_build(const char *fn) HTS_RESULT_USED;

/// Destroy a fqidx_t struct
void fqi_destroy(fqidx_t *fqi);

enum fqi_load_options {
    FQI_CREATE = 0x01,
};

/// Load FASTQ indexes.
/** @param  fn  File name of the FASTQ file (can be compressed with bgzip).
    @param  fnfqi File name of the FASTQ index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @return Pointer to a fqidx_t struct on success, NULL on failure.

If fnfqi is NULL, ".fqi" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
The bgzip index is only needed if fn is compressed.

If (flags & FAI_CREATE) is true, the index files will be built using
fqi_build3() if they are not already present.
*/
fqidx_t *fqi_load3(const char *fn, const char *fnfqi, const char *fngzi,
                   int flags);

/// Load index from "fn.fqi".
/** @param  fn  File name of the FASTQ file
    @return Pointer to a fqidx_t struct on success, NULL on failure.

This function is equivalent to fqi_load3(fn, NULL, NULL, FAI_CREATE|FAI_CACHE);
*/
fqidx_t *fqi_load(const char *fn);

/// Fetch the sequence in a region
/** @param  fqi  Pointer to the fqidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @param  len  Length of the region; -2 if seq not present, -1 general error
    @return      Pointer to the sequence; `NULL` on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *fqi_fetch(const fqidx_t *fqi, const char* reg, int *len);
char *fqi_fetchqual(const fqidx_t *fqi, const char* reg, int *len);

/// Fetch the number of sequences
/** @param  fqi  Pointer to the fqidx_t struct
    @return      The number of sequences
*/
int fqidx_fetch_nseq(const fqidx_t *fqi) HTS_DEPRECATED("Please use fqidx_nseq instead");

/// Fetch the sequence in a region
/** @param  fqi  Pointer to the fqidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on failure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *fqidx_fetch_seq(const fqidx_t *fqi, const char *c_name, int p_beg_i, int p_end_i, int *len);

/// Query if sequence is present
/**   @param  fqi  Pointer to the fqidx_t struct
      @param  seq  Sequence name
      @return      1 if present or 0 if absent
*/
int fqidx_has_seq(const fqidx_t *fqi, const char *seq);

/// Return number of sequences in fqi index
int fqidx_nseq(const fqidx_t *fqi);

/// Return name of i-th sequence
const char *fqidx_iseq(const fqidx_t *fqi, int i);

/// Return sequence length, -1 if not present
int fqidx_seq_len(const fqidx_t *fqi, const char *seq);

#ifdef __cplusplus
}
#endif

#endif
