/*  hts_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2015-2016, 2018-2020, 2025 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_HTS_INTERNAL_H
#define HTSLIB_HTS_INTERNAL_H

#include <stddef.h>
#include <ctype.h>
#include <time.h>

#include "htslib/hts.h"
#include "textutils_internal.h"

#define HTS_MAX_EXT_LEN 9

#ifdef __cplusplus
extern "C" {
#endif

struct hFILE;

struct hts_json_token {
    char type;    ///< Token type
    char *str;    ///< Value as a C string (filled in for all token types)
    // TODO Add other fields to fill in for particular data types, e.g.
    // int inum;
    // float fnum;
};

struct cram_fd;

/*
 * Adjust CSI index parameters to support max_len_in bases
 *
 * @param max_len_in         Maximum position to be indexed
 * @param min_shift_[in,out] min_shift parameter
 * @param n_lvls_[in,out]    n_lvls parameter
 *
 * Adjusts *n_lvls_ (preferred) or *min_shift_ so that the resulting values
 * can be passed to hts_idx_init(, HTS_FMT_CSI, ...) in order to make an
 * index that can store positions up to max_len_in bases.
 */
void hts_adjust_csi_settings(int64_t max_len_in, int *min_shift_, int *n_lvls_);

/*
 * Check the existence of a local index file using part of the alignment file name.
 * The order is alignment.bam.csi, alignment.csi, alignment.bam.bai, alignment.bai
 * @param fn    - pointer to the file name
 * @param fnidx - pointer to the index file name placeholder
 * @return        1 for success, 0 for failure
 */
int hts_idx_check_local(const char *fn, int fmt, char **fnidx);

// Retrieve the name of the index file and also download it, if it is remote
char *hts_idx_getfn(const char *fn, const char *ext);

// Retrieve the name of the index file, but do not download it, if it is remote
char *hts_idx_locatefn(const char *fn, const char *ext);

// Used for on-the-fly indexing.  See the comments in hts.c.
void hts_idx_amend_last(hts_idx_t *idx, uint64_t offset);

int hts_idx_fmt(hts_idx_t *idx);

// Internal interface to save on-the-fly indexes.  The index file handle
// is kept open so hts_close() can close if after writing out the EOF
// block for its own file.
int hts_idx_save_but_not_close(hts_idx_t *idx, const char *fnidx, int fmt);

// Construct a unique filename based on fname and open it.
struct hFILE *hts_open_tmpfile(const char *fname, const char *mode, kstring_t *tmpname);

// Check that index is capable of storing items in range beg..end
int hts_idx_check_range(hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end);

// The CRAM implementation stores the loaded index within the cram_fd rather
// than separately as is done elsewhere in htslib.  So if p is a pointer to
// an hts_idx_t with p->fmt == HTS_FMT_CRAI, then it actually points to an
// hts_cram_idx_t and should be cast accordingly.
typedef struct hts_cram_idx_t {
    int fmt;
    struct cram_fd *cram;
} hts_cram_idx_t;

// Determine whether the string's contents appear to be UTF-16-encoded text.
// Returns 1 if they are, 2 if there is also a BOM, or 0 otherwise.
int hts_is_utf16_text(const kstring_t *str);

// Entry point to hFILE_multipart backend.
struct hFILE *hopen_htsget_redirect(struct hFILE *hfile, const char *mode);

struct hts_path_itr {
    kstring_t path, entry;
    void *dirv;  // DIR * privately
    const char *pathdir, *prefix, *suffix;
    size_t prefix_len, suffix_len, entry_dir_l;
};

void hts_path_itr_setup(struct hts_path_itr *itr, const char *path,
    const char *builtin_path, const char *prefix, size_t prefix_len,
    const char *suffix, size_t suffix_len);

const char *hts_path_itr_next(struct hts_path_itr *itr);

typedef void plugin_void_func(void);
plugin_void_func *load_plugin(void **pluginp, const char *filename, const char *symbol);
void *plugin_sym(void *plugin, const char *name, const char **errmsg);
plugin_void_func *plugin_func(void *plugin, const char *name, const char **errmsg);
void close_plugin(void *plugin);
const char *hts_plugin_path(void);

/*
 * Buffers up arguments to hts_idx_push for later use, once we've written all bar
 * this block.  This is necessary when multiple blocks are in flight (threading).
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf_idx_push(BGZF *fp, hts_idx_t *hidx, int tid, hts_pos_t beg, hts_pos_t end, uint64_t offset, int is_mapped);

static inline int find_file_extension(const char *fn, char ext_out[static HTS_MAX_EXT_LEN])
{
    const char *delim = fn ? strstr(fn, HTS_IDX_DELIM) : NULL, *ext;
    if (!fn) return -1;
    if (!delim) delim = fn + strlen(fn);
    for (ext = delim; ext > fn && *ext != '.' && *ext != '/'; --ext) {}
    if (*ext == '.' && ext > fn &&
        ((delim - ext == 3 && ext[1] == 'g' && ext[2] == 'z') || // permit .sam.gz as a valid file extension
        (delim - ext == 4 && ext[1] == 'b' && ext[2] == 'g' && ext[3] == 'z'))) // permit .vcf.bgz as a valid file extension
    {
        for (ext--; ext > fn && *ext != '.' && *ext != '/'; --ext) {}
    }
    if (*ext != '.' || delim - ext > HTS_MAX_EXT_LEN || delim - ext < 3)
        return -1;
    memcpy(ext_out, ext + 1, delim - ext - 1);
    ext_out[delim - ext - 1] = '\0';
    return 0;
}

static inline int hts_usleep(long long usec)
{
    struct timespec req = { usec / 1000000, (usec % 1000000) * 1000 };
    return nanosleep(&req, NULL);
}

/*!
  @abstract   Is SVLEN the reference length for a VCF ALT allele?
  @param alt  ALT allele
  @param size Length of @p alt; -1 if not known
  @return     1 if yes; 0 if not.

  This is used when reading VCF and in tabix to check if SVLEN should be taken
  into account when working out the reference length.  It should if the
  ALT allele is a symbolic one of type CNV, DEL, DUP or INV, plus
  sub-types like <CNV:TR> or <DEL:ME>.

  @p alt does not have to be NUL-terminated, but if not @p size should be
  greater than of equal to zero.  If @p is less than zero, @p alt must be
  NUL-terminated.
*/

static inline int svlen_on_ref_for_vcf_alt(const char *alt, int32_t size)
{
    size_t sz;
    if (*alt != '<') // Check if ALT is symbolic
        return 0;
    sz = size >= 0 ? (size_t) size : strlen(alt);
    if (sz < 5)      // Reject if not long enough
        return 0;
    if (alt[4] != '>' && alt[4] != ':')  // Reject if too long
        return 0;
    if (memcmp(alt, "<CNV", 4) != 0     // Copy-number variation
        && memcmp(alt, "<DEL", 4) != 0  // Deletion
        && memcmp(alt, "<DUP", 4) != 0  // Duplication
        && memcmp(alt, "<INV", 4) != 0) // Inversion
        return 0;
    return alt[sz - 1] == '>' ? 1 : 0; // Check symbolic allele ends correctly
}

#ifdef __cplusplus
}
#endif

#endif
