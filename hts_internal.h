/*  hts_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2015-2016, 2018-2019 Genome Research Ltd.

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

#include "htslib/hts.h"

#include "textutils_internal.h"

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
 * Check the existence of a local index file using part of the alignment file name.
 * The order is alignment.bam.csi, alignment.csi, alignment.bam.bai, alignment.bai
 * @param fn    - pointer to the file name
 * @param fnidx - pointer to the index file name placeholder
 * @return        1 for success, 0 for failure
 */
int hts_idx_check_local(const char *fn, int fmt, char **fnidx);

// Retrieve the name of the index file and also download it, if it is remote
char *hts_idx_getfn(const char *fn, const char *ext);

// Used for on-the-fly indexing.  See the comments in hts.c.
void hts_idx_amend_last(hts_idx_t *idx, uint64_t offset);

int hts_idx_fmt(hts_idx_t *idx);

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

void *load_plugin(void **pluginp, const char *filename, const char *symbol);
void *plugin_sym(void *plugin, const char *name, const char **errmsg);
void close_plugin(void *plugin);

/*
 * Buffers up arguments to hts_idx_push for later use, once we've written all bar
 * this block.  This is necessary when multiple blocks are in flight (threading).
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int bgzf_idx_push(BGZF *fp, hts_idx_t *hidx, int tid, hts_pos_t beg, hts_pos_t end, uint64_t offset, int is_mapped);

#ifdef __cplusplus
}
#endif

#endif
