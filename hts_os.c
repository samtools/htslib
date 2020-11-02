/// @file hts_os.c
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017, 2019 Genome Research Ltd.

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
#include "htslib/hts_defs.h"

// Windows (maybe more) lack a drand48 implementation.
#ifndef HAVE_DRAND48
#include "os/rand.c"
#else
#include <stdlib.h>
HTSLIB_EXPORT
void hts_srand48(long seed)
{
#ifdef HAVE_SRAND48_DETERMINISTIC
    srand48_deterministic(seed);
#else
    srand48(seed);
#endif
}

HTSLIB_EXPORT
double hts_erand48(unsigned short xseed[3]) { return erand48(xseed); }

HTSLIB_EXPORT
double hts_drand48(void) { return drand48(); }

HTSLIB_EXPORT
long hts_lrand48(void) { return lrand48(); }
#endif

// // On Windows when using the MSYS or Cygwin terminals, isatty fails
// #ifdef _WIN32
// #define USE_FILEEXTD
// #include "os/iscygpty.c"
// #endif


#include <stdint.h>
#include <string.h>
#include "hts_internal.h"
#include "htslib/hts.h"
#include "htslib/hts_os.h"
#include "htslib/kstring.h"

unsigned int htslib_features(void) {
    unsigned int feat = 0;

#ifdef PACKAGE_URL
    feat |= HTS_FEATURE_CONFIGURE;
#endif

#ifdef ENABLE_PLUGINS
    feat |= HTS_FEATURE_PLUGINS;
#endif

#ifdef HAVE_LIBCURL
    feat |= HTS_FEATURE_LIBCURL;
#endif

#ifdef ENABLE_S3
    feat |= HTS_FEATURE_S3;
#endif

#ifdef ENABLE_GCS
    feat |= HTS_FEATURE_GCS;
#endif

#ifdef HAVE_LIBDEFLATE
    feat |= HTS_FEATURE_LIBDEFLATE;
#endif

#ifdef HAVE_LIBLZMA
    feat |= HTS_FEATURE_LZMA;
#endif

#ifdef HAVE_LIBBZ2
    feat |= HTS_FEATURE_BZIP2;
#endif

    return feat;
}

const char *htslib_test_feature(int id) {
    int feat = htslib_features();

    switch (id) {
    case HTS_FEATURE_CONFIGURE:
        return feat & HTS_FEATURE_CONFIGURE ? "yes" : NULL;
    case HTS_FEATURE_PLUGINS:
        return feat & HTS_FEATURE_PLUGINS ? "yes" : NULL;
    case HTS_FEATURE_LIBCURL:
        return feat & HTS_FEATURE_LIBCURL ? "yes" : NULL;
    case HTS_FEATURE_S3:
        return feat & HTS_FEATURE_S3 ? "yes" : NULL;
    case HTS_FEATURE_GCS:
        return feat & HTS_FEATURE_GCS ? "yes" : NULL;
    case HTS_FEATURE_LIBDEFLATE:
        return feat & HTS_FEATURE_LIBDEFLATE ? "yes" : NULL;
    case HTS_FEATURE_BZIP2:
        return feat & HTS_FEATURE_BZIP2 ? "yes" : NULL;
    case HTS_FEATURE_LZMA:
        return feat & HTS_FEATURE_LZMA ? "yes" : NULL;

    case HTS_FEATURE_CC:
        return HTS_CC;
    case HTS_FEATURE_CFLAGS:
        return HTS_CFLAGS;
    case HTS_FEATURE_LDFLAGS:
        return HTS_LDFLAGS;
    case HTS_FEATURE_CPPFLAGS:
        return HTS_CPPFLAGS;

    default:
        fprintf(stderr, "Unknown feature code: %d\n", id);
    }

    return NULL;
}

// Note this implementation also means we can just "strings" the library
// to find the configuration parameters.
const char *htslib_feature_string(void) {
    const char *fmt=

#ifdef PACKAGE_URL
    "build=configure "
#else
    "build=Makefile "
#endif

#ifdef ENABLE_PLUGINS
    "plugins=yes, plugin-path=%.1000s "
#else
    "plugins=no "
#endif

#ifdef HAVE_LIBCURL
    "libcurl=yes "
#else
    "libcurl=no "
#endif

#ifdef ENABLE_S3
    "S3=yes "
#else
    "S3=no "
#endif

#ifdef ENABLE_GCS
    "GCS=yes "
#else
    "GCS=no "
#endif

#ifdef HAVE_LIBDEFLATE
    "libdeflate=yes "
#else
    "libdeflate=no "
#endif

#ifdef HAVE_LIBLZMA
    "lzma=yes "
#else
    "lzma=no "
#endif

#ifdef HAVE_LIBBZ2
    "bzip2=yes ";
#else
    "bzip2=no ";
#endif

#ifdef ENABLE_PLUGINS
    static char config[1200];
    sprintf(config, fmt, htslib_plugin_path());
    return config;
#else
    return fmt;
#endif
}

// Plus hts_version here?
