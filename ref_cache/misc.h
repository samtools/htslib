/*  misc.h -- ref-cache miscellaneous interfaces

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

#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include "../htslib/hts_defs.h"

extern const signed char hexvals[256];

#define hexval(C) (hexvals[(unsigned char) (C)])

static inline int setnonblock(int fd) {
    int val;

    if ((val = fcntl(fd, F_GETFL)) == -1) {
        perror("Couldn't get file descriptor flags");
        return -1;
    }

    if (fcntl(fd, F_SETFL, val | O_NONBLOCK)) {
        perror("Couldn't set socket to non-blocking mode");
        return -1;
    }
    return 0;
}

static inline
HTS_ACCESS(read_only, 2, 3)
ssize_t do_write_all(int fd, const void *buf, size_t count) {
    ssize_t res = 0;
    const unsigned char *ucbuf = buf;
    while (count > 0) {
        do {
            res = write(fd, ucbuf, count);
        } while (res < 0
                 && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
        if (res < 0) break;
        count -= (size_t) res;
        ucbuf += res;
    }
    return res >= 0 ? 0 : -1;
}

static inline
HTS_ACCESS(write_only, 2, 3)
ssize_t do_read_all(int fd, void *buf, size_t count) {
    ssize_t res = 0;
    ssize_t bytes = 0;
    unsigned char *ucbuf = buf;

    while ((size_t) bytes < count) {
        do {
            res = read(fd, ucbuf, count);
        } while (res < 0
                 && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
        if (res <= 0) break;
        bytes += res;
        ucbuf += res;
    }
    return res < 0 ? res : bytes;
}

static inline
HTS_ACCESS(read_only, 1, 2)
char * lim_strdup(const char *str, size_t len, size_t max_len) {
    char *out;

    if (len ==  0) return NULL;
    if (len < max_len) {
        out = malloc(len + 1);
        if (out == NULL) return NULL;
        memcpy(out, str, len);
        out[len] = '\0';
        return out;
    }
    out = malloc(max_len + 1);
    if (out == NULL) return NULL;
    memcpy(out, str, max_len - 3);
    memcpy(out + max_len - 3, "...", 4);
    return out;
}

#endif /* MISC_H_INCLUDED */
