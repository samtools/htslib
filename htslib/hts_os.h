/// @file hts_os.h
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

#ifndef HTSLIB_HTS_OS_H
#define HTSLIB_HTS_OS_H

extern void srand48(long seed);
extern double drand48(void);
extern long lrand48(void);

#ifdef _WIN32
/* Check if the fd is a cygwin/msys's pty. */
extern int is_cygpty(int fd);
#endif

#if defined(__MINGW32__)
#include <io.h>
#define mkdir(filename,mode) mkdir((filename))
#endif

#ifdef _WIN32
#include <stdlib.h>
#define srandom srand
#define random rand
#endif

#endif
