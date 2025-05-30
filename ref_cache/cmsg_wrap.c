/*  cmsg_wrap.c -- wrapper around cmsg interfaces

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

// We may disagree with some configure settings as certain definitions
// cause CMSG_LEN() to be hidden.  However, overriding them may also
// remove things we want for some compiler settings, so the code that
// needs it is isolated in this file.

#if defined(_XOPEN_SOURCE)
#  undef _XOPEN_SOURCE
#endif
#if defined(_POSIX_C_SOURCE)
#  undef _POSIX_C_SOURCE
#endif

#include <stddef.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>

#include "cmsg_wrap.h"

int make_scm_rights_cmsg(struct msghdr *msg, int fd,
                         char *buf, size_t buf_sz) {
    size_t needed = (size_t) CMSG_SPACE(sizeof(fd));
    struct cmsghdr *cmsg;
    unsigned char *fdptr;
    if (needed > buf_sz) {
        return -1;
    } else {
        msg->msg_control = buf;
    }
    msg->msg_controllen = (socklen_t) needed;
    memset(msg->msg_control, 0, (size_t) needed);
    cmsg = CMSG_FIRSTHDR(msg);
    cmsg->cmsg_level = SOL_SOCKET;
    cmsg->cmsg_type  = SCM_RIGHTS;
    cmsg->cmsg_len   = CMSG_LEN(sizeof(fd));
    fdptr = CMSG_DATA(cmsg);
    memcpy(fdptr, &fd, sizeof(fd));
    return 0;
}

int get_scm_rights_fd(struct msghdr *msg) {
    struct cmsghdr *cmsg;
    int fd = -1;

    for (cmsg = CMSG_FIRSTHDR(msg); cmsg != NULL; CMSG_NXTHDR(msg, cmsg)) {
        if (cmsg->cmsg_level == SOL_SOCKET && cmsg->cmsg_type == SCM_RIGHTS) {
            unsigned char *fdp = (unsigned char *) CMSG_DATA(cmsg);
            memcpy(&fd, fdp, sizeof(int));
            break;
        }
    }
    return fd;
}
