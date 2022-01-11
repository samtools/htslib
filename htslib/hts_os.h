/// @file hts_os.h
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017, 2019-2020 Genome Research Ltd.

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

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This is srand48_deterministic() on platforms that provide it, or srand48()
   otherwise (or our own POSIX srand48() on platforms that provide neither).
   Hence calling hts_srand48() will always set up the same POSIX-determined
   sequence of pseudo-random numbers on any platform, while calling srand48()
   may (e.g., on OpenBSD) set up a different non-deterministic sequence. */
HTSLIB_EXPORT
void hts_srand48(long seed);

HTSLIB_EXPORT
double hts_erand48(unsigned short xseed[3]);

HTSLIB_EXPORT
double hts_drand48(void);

HTSLIB_EXPORT
long hts_lrand48(void);

#if defined(_WIN32) && !defined(__CYGWIN__)
// Windows usually lacks *rand48(), but cygwin provides them.
#define srand48(S) hts_srand48((S))
#define erand48(X) hts_erand48((X))
#define drand48() hts_drand48()
#define lrand48() hts_lrand48()
#endif

#if 0  /* def _WIN32 - disabled for now, not currently used */
/* Check if the fd is a cygwin/msys's pty. */
extern int is_cygpty(int fd);
#endif

#ifdef __cplusplus
}
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

/* MSVC does not provide ssize_t in its <sys/types.h>. This ensures the type
   is available (unless suppressed by defining HTS_NO_SSIZE_T first). */
#if defined _MSC_VER && defined _INTPTR_T_DEFINED && !defined _SSIZE_T_DEFINED && !defined HTS_NO_SSIZE_T && !defined ssize_t
#define ssize_t intptr_t
#endif

#endif // HTSLIB_HTS_OS_H
