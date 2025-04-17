/*  ref_files.h -- ref-cache reference file handler

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

#ifndef REF_FILES_H
#define REF_FILES_H

#include <sys/types.h>
#include <stdint.h>

#include "../htslib/hts_defs.h"
#include "types.h"

#define MD5_LEN 32

typedef enum {
    REF_WAITING_UPSTREAM = 0,
    REF_NOT_FOUND        = 1,
    REF_DOWNLOAD_STARTED = 2,
    REF_IS_COMPLETE      = 3
} RefFileStatus;

RefFile * get_ref_file(const Options *opts, const char *md5, int upstream_fd)
HTS_ACCESS(read_only, 1)
HTS_ACCESS(read_only, 2);

RefFileStatus get_ref_status(const RefFile *ref)
HTS_ACCESS(read_only, 1);

off_t get_ref_size(const RefFile *ref)
HTS_ACCESS(read_only, 1);

off_t get_ref_available(const RefFile *ref)
HTS_ACCESS(read_only, 1);

unsigned int get_ref_id(const RefFile *ref)
HTS_ACCESS(read_only, 1);

int get_ref_complete(const RefFile *ref)
HTS_ACCESS(read_only, 1);

int get_ref_fd(const RefFile *ref)
HTS_ACCESS(read_only, 1);

void update_ref_download_started(RefFile *ref, int fd,
                                 int64_t size_if_complete);

void update_ref_available(RefFile *ref, int64_t available);

void update_ref_with_content_len(RefFile *ref, int64_t size);

int set_ref_complete(RefFile *ref);

int release_ref_file(RefFile *ref);

#endif /* REF_FILES_H */
