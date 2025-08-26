/*  bgzf_internal.h -- internal bgzf functions; not part of the public API.

    Copyright (C) 2025 Genome Research Ltd.

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

#include <assert.h>
#include "htslib/bgzf.h"

/*
 * BGZF private data interface
 * This exists so that we can pass BCF headers into interfaces that have
 * traditionally only taken a BGZF pointer without a corresponding bcf_hdr_t *,
 * notably the bcf_readrec() function used by BCF iterators.
 *
 * To preserve the BGZF API and ABI, this is tagged on to the existing
 * opaque bgzf_cache_t structure.  bgzf_cache_t is now defined here so we can
 * inline lookups.
 */

typedef void bgzf_private_data_cleanup_func(void *private_data);

struct kh_bgzf_cache_s;

struct bgzf_cache_t {
    struct kh_bgzf_cache_s *h;
    unsigned int last_pos;
    void *private_data;
    bgzf_private_data_cleanup_func *private_data_cleanup;
};

// Set private data.  cleanup will be called on bgzf_close() or
// bgzf_clear_private_data();

static inline void bgzf_set_private_data(BGZF *fp, void *private_data,
                                         bgzf_private_data_cleanup_func *fn) {
    assert(fp->cache != NULL);
    fp->cache->private_data = private_data;
    fp->cache->private_data_cleanup = fn;
}

static inline void bgzf_clear_private_data(BGZF *fp) {
    assert(fp->cache != NULL);
    if (fp->cache->private_data) {
        if (fp->cache->private_data_cleanup)
            fp->cache->private_data_cleanup(fp->cache->private_data);
        fp->cache->private_data = NULL;
    }
}

static inline void * bgzf_get_private_data(BGZF *fp) {
    assert(fp->cache != NULL);
    return fp->cache->private_data;
}
