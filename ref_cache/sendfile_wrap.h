/*  sendfile_wrap.h -- wrapper around sendfile interfaces

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

#include <sys/types.h>

// Define HAVE_SENDFILE if configure found a sendfile interface
#if  (defined(HAVE_LINUX_SENDFILE) && HAVE_LINUX_SENDFILE) \
  || (defined(HAVE_FREEBSD_SENDFILE) && HAVE_FREEBSD_SENDFILE) \
  || (defined(HAVE_MACOS_SENDFILE) && HAVE_MACOS_SENDFILE)
#ifndef HAVE_SENDFILE
#  define HAVE_SENDFILE
#endif
#elif defined(HAVE_SENDFILE)
#undef HAVE_SENDFILE
#endif

ssize_t sendfile_wrap(int out_fd, int in_fd, off_t *offset, size_t count);
