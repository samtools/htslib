/*  sendfile_wrap.c -- wrapper around sendfile interfaces

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

#if defined(HAVE_LINUX_SENDFILE) && HAVE_LINUX_SENDFILE

#include <sys/sendfile.h>

#elif (defined(HAVE_FREEBSD_SENDFILE) && HAVE_FREEBSD_SENDFILE) \
  ||  (defined(HAVE_MACOS_SENDFILE) && HAVE_MACOS_SENDFILE)

// Need to remove these else the sendfile() definition may be hidden
#if defined(_XOPEN_SOURCE)
#  undef _XOPEN_SOURCE
#endif
#if defined(_POSIX_C_SOURCE)
#  undef _POSIX_C_SOURCE
#endif
#if defined(HAVE_MACOS_SENDFILE) && HAVE_MACOS_SENDFILE && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/uio.h>
#include <stddef.h>
#include <errno.h>
#endif

#include "sendfile_wrap.h"

#if defined(HAVE_LINUX_SENDFILE) && HAVE_LINUX_SENDFILE

ssize_t sendfile_wrap(int out_fd, int in_fd, off_t *offset, size_t count) {
    return sendfile(out_fd, in_fd, offset, count);
}

#elif defined(HAVE_FREEBSD_SENDFILE) && HAVE_FREEBSD_SENDFILE

ssize_t sendfile_wrap(int out_fd, int in_fd, off_t *offset, size_t count) {
    off_t sbytes = 0;
    int res;
    if (count == 0)
        return 0; // Else the whole file is sent
    res = sendfile(in_fd, out_fd, *offset, count, NULL, &sbytes, 0);
    *offset += sbytes;
    return res < 0 ? res : (ssize_t) sbytes;
}

#elif defined(HAVE_MACOS_SENDFILE) && HAVE_MACOS_SENDFILE

ssize_t sendfile_wrap(int out_fd, int in_fd, off_t *offset, size_t count) {
    off_t len = (off_t) count;
    int res;
    if (len == 0)
        return 0; // Else the whole file is sent
    res = sendfile(in_fd, out_fd, *offset, &len, NULL, 0);
    if (res == 0 || errno == EINTR || errno == EAGAIN)
        *offset += len;
    return res < 0 ? res : (ssize_t) len;
}

#else

// This should never be called
ssize_t sendfile_wrap(int out_fd, int in_fd, off_t *offset, size_t count) {
    // Stop complaints about unused variables
    if (out_fd >= 0 || in_fd >= 0 || offset || count)
        return -2;
    return -1;
}

#endif
