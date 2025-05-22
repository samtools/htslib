/*  ref_files.c -- ref-cache reference file handling

    Copyright (C) 2025 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

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

#include <config.h>

#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "ref_files.h"
#include "misc.h"
#include "options.h"
#include "upstream.h"

#define HASH_SZ   0x10000
#define HASH_MASK (HASH_SZ - 1)

struct RefFile {
    char hexmd5[MD5_LEN];
    struct RefFile *prev_md5;
    struct RefFile *next_md5;
    off_t           size;
    off_t           available;
    unsigned int    ref_count;
    unsigned int    id;
    RefFileStatus   status;
    int             fd;
};

typedef struct RefFiles {
    RefFile *by_md5[HASH_SZ];
    unsigned int id;
} RefFiles;

RefFiles refs = { { NULL }, 0 };

static RefFile *get_ref_placeholder(const char *md5) {
    RefFile *r;
    unsigned int m5hash = (hexval(md5[0]) << 12 |
                           hexval(md5[1]) << 8  |
                           hexval(md5[2]) << 4  |
                           hexval(md5[3]) << 0) & HASH_MASK;
    /* See if it's already there */
    for (r = refs.by_md5[m5hash]; r != NULL; r = r->next_md5) {
        if (strncmp(md5, r->hexmd5, MD5_LEN) == 0) {
            r->ref_count++;
            return r;
        }
    }

    r = calloc(1, sizeof(RefFile));
    if (r == NULL)
        return NULL;

    r->status = REF_WAITING_UPSTREAM;
    r->fd = -1;
    r->ref_count = 1;
    r->id = ++refs.id;
    memcpy(r->hexmd5, md5, MD5_LEN);
    r->prev_md5 = NULL;
    if ((r->next_md5 = refs.by_md5[m5hash]) != NULL) {
        r->next_md5->prev_md5 = r;
    }
    refs.by_md5[m5hash] = r;

    return r;
}

RefFile *get_ref_file(const Options *opts, const char *md5, int upstream_fd) {
    RefFile *r = get_ref_placeholder(md5);
    char fname[MD5_LEN + 3];
    struct stat stat_buf;

    /* Already opened / created */
    if (r->ref_count > 1) return r;

    fname[0] = md5[0]; fname[1] = md5[1]; fname[2] = '/';
    fname[3] = md5[2]; fname[4] = md5[3]; fname[5] = '/';
    memcpy(fname + 6, md5 + 4, MD5_LEN - 4);
    fname[MD5_LEN + 2] = '\0';

    r->fd = openat(opts->cache_fd, fname, O_RDONLY);
    if (r->fd < 0) {
        if (errno == ENOENT) {
            if (upstream_fd >= 0) {
                /* Send request to download the file */
                if (upstream_send_cmd(upstream_fd, md5, r->id) != 0)
                    goto fail;
                r->status = REF_WAITING_UPSTREAM;
            } else {
                /* No upstream, set not found status */
                r->status = REF_NOT_FOUND;
            }
            return r;
        } else {
            goto fail;
        }
    }

    if (fstat(r->fd, &stat_buf) != 0) {
        fprintf(stderr, "Couldn't get length of %s/%s : %s\n",
                opts->cache_dir, fname, strerror(errno));
        goto fail;
    }

    r->size = r->available = stat_buf.st_size;
    r->status = REF_IS_COMPLETE;

    return r;

 fail:
    release_ref_file(r);
    return NULL;
}

RefFileStatus get_ref_status(const RefFile *ref) {
    return ref->status;
}

off_t get_ref_size(const RefFile *ref) {
    return ref->size;
}

off_t get_ref_available(const RefFile *ref) {
    return ref->available;
}

unsigned int get_ref_id(const RefFile *ref) {
    return ref->id;
}

int get_ref_complete(const RefFile *ref) {
    return ref->status == REF_IS_COMPLETE;
}

int get_ref_fd(const RefFile *ref) {
    return ref->fd;
}

void update_ref_download_started(RefFile *ref, int fd,
                                 int64_t size_if_complete) {
    ref->fd = fd;
    if (size_if_complete >= 0) {
        ref->status = REF_IS_COMPLETE;
        ref->size = ref->available = (off_t) size_if_complete;
    }
}

void update_ref_available(RefFile *ref, int64_t available) {
    assert(ref->available <= available);
    ref->available = available;
}

void update_ref_with_content_len(RefFile *ref, int64_t size) {
    ref->size = (off_t) size;
    if (ref->status < REF_DOWNLOAD_STARTED)
        ref->status = REF_DOWNLOAD_STARTED;
}

int set_ref_complete(RefFile *ref) {
    assert(ref != NULL);
    int no_content_length = ref->size == 0;
    ref->status = REF_IS_COMPLETE;
    ref->size = ref->available;
    return no_content_length;
}

int release_ref_file(RefFile *ref) {
    int res;

    if (--ref->ref_count > 0) return 0;

    /* Remove from md5 linked list */
    if (ref->prev_md5 == NULL) {
        unsigned int m5hash = (hexval(ref->hexmd5[0]) << 12 |
                               hexval(ref->hexmd5[1]) << 8  |
                               hexval(ref->hexmd5[2]) << 4  |
                               hexval(ref->hexmd5[3]) << 0) & HASH_MASK;
        refs.by_md5[m5hash] = ref->next_md5;
    } else {
        ref->prev_md5->next_md5 = ref->next_md5;
    }
    if (ref->next_md5 != NULL) ref->next_md5->prev_md5 = ref->prev_md5;

    /* Close the file (if open) */
    res = ref->fd >= 0 ? close(ref->fd) : 0;

    free(ref);

    return res;
}
