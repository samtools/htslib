/*  hts_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2015-2016 Genome Research Ltd.

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

#ifdef __cplusplus
extern "C" {
#endif

// The <ctype.h> functions operate on ints such as are returned by fgetc(),
// i.e., characters represented as unsigned-char-valued ints, or EOF.
// To operate on plain chars (and to avoid warnings on some platforms),
// technically one must cast to unsigned char everywhere (see CERT STR37-C)
// or less painfully use these *_c() functions that operate on plain chars
// (but not EOF, which must be considered separately where it is applicable).
// TODO We may eventually wish to implement these functions directly without
// using their <ctype.h> equivalents, and thus make them immune to locales.
static inline int isalnum_c(char c) { return isalnum((unsigned char) c); }
static inline int isalpha_c(char c) { return isalpha((unsigned char) c); }
static inline int isdigit_c(char c) { return isdigit((unsigned char) c); }
static inline int isgraph_c(char c) { return isgraph((unsigned char) c); }
static inline int islower_c(char c) { return islower((unsigned char) c); }
static inline int isprint_c(char c) { return isprint((unsigned char) c); }
static inline int isspace_c(char c) { return isspace((unsigned char) c); }
static inline int isupper_c(char c) { return isupper((unsigned char) c); }
static inline char tolower_c(char c) { return tolower((unsigned char) c); }
static inline char toupper_c(char c) { return toupper((unsigned char) c); }


struct cram_fd;

char *hts_idx_getfn(const char *fn, const char *ext);

// The CRAM implementation stores the loaded index within the cram_fd rather
// than separately as is done elsewhere in htslib.  So if p is a pointer to
// an hts_idx_t with p->fmt == HTS_FMT_CRAI, then it actually points to an
// hts_cram_idx_t and should be cast accordingly.
typedef struct hts_cram_idx_t {
    int fmt;
    struct cram_fd *cram;
} hts_cram_idx_t;


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

#ifdef __cplusplus
}
#endif

#endif
