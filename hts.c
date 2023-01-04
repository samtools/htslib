/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2023 Genome Research Ltd.
    Copyright (C) 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef HAVE_LIBLZMA
#ifdef HAVE_LZMA_H
#include <lzma.h>
#else
#include "os/lzma_stub.h"
#endif
#endif

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "version.h"
#include "config_vars.h"
#include "hts_internal.h"
#include "hfile_internal.h"
#include "sam_internal.h"
#include "htslib/hts_expr.h"
#include "htslib/hts_os.h" // drand48

#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/ksort.h"
#include "htslib/tbx.h"
#if defined(HAVE_EXTERNAL_LIBHTSCODECS)
#include <htscodecs/htscodecs.h>
#else
#include "htscodecs/htscodecs/htscodecs.h"
#endif

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif

KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

HTSLIB_EXPORT
int hts_verbose = HTS_LOG_WARNING;

const char *hts_version()
{
    return HTS_VERSION_TEXT;
}

unsigned int hts_features(void) {
    unsigned int feat = HTS_FEATURE_HTSCODECS; // Always present

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

const char *hts_test_feature(unsigned int id) {
    unsigned int feat = hts_features();

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

    case HTS_FEATURE_HTSCODECS:
        return htscodecs_version();

    case HTS_FEATURE_CC:
        return HTS_CC;
    case HTS_FEATURE_CFLAGS:
        return HTS_CFLAGS;
    case HTS_FEATURE_LDFLAGS:
        return HTS_LDFLAGS;
    case HTS_FEATURE_CPPFLAGS:
        return HTS_CPPFLAGS;

    default:
        fprintf(stderr, "Unknown feature code: %u\n", id);
    }

    return NULL;
}

// Note this implementation also means we can just "strings" the library
// to find the configuration parameters.
const char *hts_feature_string(void) {
    static char config[1200];
    const char *flags=

#ifdef PACKAGE_URL
    "build=configure "
#else
    "build=Makefile "
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
    "bzip2=yes "
#else
    "bzip2=no "
#endif

// "plugins=" must stay at the end as it is followed by "plugin-path="
#ifdef ENABLE_PLUGINS
    "plugins=yes";
#else
    "plugins=no";
#endif

#ifdef ENABLE_PLUGINS
    snprintf(config, sizeof(config),
             "%s plugin-path=%.1000s htscodecs=%.40s",
             flags, hts_plugin_path(), htscodecs_version());
#else
    snprintf(config, sizeof(config),
             "%s htscodecs=%.40s",
             flags, htscodecs_version());
#endif
    return config;
}


HTSLIB_EXPORT
const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

HTSLIB_EXPORT
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

HTSLIB_EXPORT
const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

/**********************
 *** Basic file I/O ***
 **********************/

static enum htsFormatCategory format_category(enum htsExactFormat fmt)
{
    switch (fmt) {
    case bam:
    case sam:
    case cram:
    case fastq_format:
    case fasta_format:
        return sequence_data;

    case vcf:
    case bcf:
        return variant_data;

    case bai:
    case crai:
    case csi:
    case fai_format:
    case fqi_format:
    case gzi:
    case tbi:
        return index_file;

    case bed:
    case d4_format:
        return region_list;

    case htsget:
    case hts_crypt4gh_format:
        return unknown_category;

    case unknown_format:
    case binary_format:
    case text_format:
    case empty_format:
    case format_maximum:
        break;
    }

    return unknown_category;
}

// Decompress several hundred bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static ssize_t
decompress_peek_gz(hFILE *fp, unsigned char *dest, size_t destsize)
{
    unsigned char buffer[2048];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return -1;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return -1;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}

#ifdef HAVE_LIBLZMA
// Similarly decompress a portion by peeking at the file, which must be
// positioned at the start of the file.
static ssize_t
decompress_peek_xz(hFILE *fp, unsigned char *dest, size_t destsize)
{
    unsigned char buffer[2048];
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);
    if (npeek < 0) return -1;

    lzma_stream ls = LZMA_STREAM_INIT;
    if (lzma_stream_decoder(&ls, lzma_easy_decoder_memusage(9), 0) != LZMA_OK)
        return -1;

    ls.next_in = buffer;
    ls.avail_in = npeek;
    ls.next_out = dest;
    ls.avail_out = destsize;

    int r = lzma_code(&ls, LZMA_RUN);
    if (! (r == LZMA_OK || r == LZMA_STREAM_END)) {
        lzma_end(&ls);
        return -1;
    }

    destsize = ls.total_out;
    lzma_end(&ls);

    return destsize;
}
#endif

// Parse "x.y" text, taking care because the string is not NUL-terminated
// and filling in major/minor only when the digits are followed by a delimiter,
// so we don't misread "1.10" as "1.1" due to reaching the end of the buffer.
static void
parse_version(htsFormat *fmt, const unsigned char *u, const unsigned char *ulim)
{
    const char *s    = (const char *) u;
    const char *slim = (const char *) ulim;
    short v;

    fmt->version.major = fmt->version.minor = -1;

    for (v = 0; s < slim && isdigit_c(*s); s++)
        v = 10 * v + *s - '0';

    if (s < slim) {
        fmt->version.major = v;
        if (*s == '.') {
            s++;
            for (v = 0; s < slim && isdigit_c(*s); s++)
                v = 10 * v + *s - '0';
            if (s < slim)
                fmt->version.minor = v;
        }
        else
            fmt->version.minor = 0;
    }
}

static int
cmp_nonblank(const char *key, const unsigned char *u, const unsigned char *ulim)
{
    const unsigned char *ukey = (const unsigned char *) key;

    while (*ukey)
        if (u >= ulim) return +1;
        else if (isspace_c(*u)) u++;
        else if (*u != *ukey) return (*ukey < *u)? -1 : +1;
        else u++, ukey++;

    return 0;
}

static int is_text_only(const unsigned char *u, const unsigned char *ulim)
{
    for (; u < ulim; u++)
        if (! (*u >= ' ' || *u == '\t' || *u == '\r' || *u == '\n'))
            return 0;

    return 1;
}

static int is_fastaq(const unsigned char *u, const unsigned char *ulim)
{
    const unsigned char *eol = memchr(u, '\n', ulim - u);

    // Check that the first line is entirely textual
    if (! is_text_only(u, eol? eol : ulim)) return 0;

    // If the first line is very long, consider the file to indeed be FASTA/Q
    if (eol == NULL) return 1;

    u = eol+1; // Now points to the first character of the second line

    // Scan over all base-encoding letters (including 'N' but not SEQ's '=')
    while (u < ulim && (seq_nt16_table[*u] != 15 || toupper(*u) == 'N')) {
        if (*u == '=') return 0;
        u++;
    }

    return (u == ulim || *u == '\r' || *u == '\n')? 1 : 0;
}

// Parse tab-delimited text, filling in a string of column types and returning
// the number of columns spotted (within [u,ulim), and up to column_len) or -1
// if non-printable characters were seen.  Column types:
//     i: integer, s: strand sign, C: CIGAR, O: SAM optional field, Z: anything
static int
parse_tabbed_text(char *columns, int column_len,
                  const unsigned char *u, const unsigned char *ulim,
                  int *complete)
{
    const char *str  = (const char *) u;
    const char *slim = (const char *) ulim;
    const char *s;
    int ncolumns = 0;

    enum { digit = 1, leading_sign = 2, cigar_operator = 4, other = 8 };
    unsigned seen = 0;
    *complete = 0;

    for (s = str; s < slim; s++)
        if (*s >= ' ') {
            if (isdigit_c(*s))
                seen |= digit;
            else if ((*s == '+' || *s == '-') && s == str)
                seen |= leading_sign;
            else if (strchr(BAM_CIGAR_STR, *s) && s > str && isdigit_c(s[-1]))
                seen |= cigar_operator;
            else
                seen |= other;
        }
        else if (*s == '\t' || *s == '\r' || *s == '\n') {
            size_t len = s - str;
            char type;

            if (seen == digit || seen == (leading_sign|digit)) type = 'i';
            else if (seen == (digit|cigar_operator)) type = 'C';
            else if (len == 1)
                switch (str[0]) {
                case '*': type = 'C'; break;
                case '+': case '-': case '.': type = 's'; break;
                default: type = 'Z'; break;
                }
            else if (len >= 5 && str[2] == ':' && str[4] == ':') type = 'O';
            else type = 'Z';

            columns[ncolumns++] = type;
            if (*s != '\t' || ncolumns >= column_len - 1) {
                *complete = 1; // finished the line or more columns than needed
                break;
            }

            str = s + 1;
            seen = 0;
        }
        else return -1;

    columns[ncolumns] = '\0';
    return ncolumns;
}

// Match COLUMNS as a prefix against PATTERN (so COLUMNS may run out first).
// Returns len(COLUMNS) (modulo '+'), or 0 if there is a mismatched entry.
static int colmatch(const char *columns, const char *pattern)
{
    int i;
    for (i = 0; columns[i] != '\0'; i++) {
        if (pattern[i] == '+') return i;
        if (! (columns[i] == pattern[i] || pattern[i] == 'Z')) return 0;
    }

    return i;
}

int hts_detect_format(hFILE *hfile, htsFormat *fmt)
{
    return hts_detect_format2(hfile, NULL, fmt);
}

int hts_detect_format2(hFILE *hfile, const char *fname, htsFormat *fmt)
{
    char extension[HTS_MAX_EXT_LEN], columns[24];
    unsigned char s[1024];
    int complete = 0;
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return -1;

    fmt->category = unknown_category;
    fmt->format = unknown_format;
    fmt->version.major = fmt->version.minor = -1;
    fmt->compression = no_compression;
    fmt->compression_level = -1;
    fmt->specific = NULL;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few records or lines.
        fmt->compression = gzip;
        if (len >= 18 && (s[3] & 4)) {
            if (memcmp(&s[12], "BC\2\0", 4) == 0)
                fmt->compression = bgzf;
            else if (memcmp(&s[12], "RAZF", 4) == 0)
                fmt->compression = razf_compression;
        }
        if (len >= 9 && s[2] == 8)
            fmt->compression_level = (s[8] == 2)? 9 : (s[8] == 4)? 1 : -1;

        len = decompress_peek_gz(hfile, s, sizeof s);
    }
    else if (len >= 10 && memcmp(s, "BZh", 3) == 0 &&
             (memcmp(&s[4], "\x31\x41\x59\x26\x53\x59", 6) == 0 ||
              memcmp(&s[4], "\x17\x72\x45\x38\x50\x90", 6) == 0)) {
        fmt->compression = bzip2_compression;
        fmt->compression_level = s[3] - '0';
        // Decompressing via libbz2 produces no output until it has a whole
        // block (of size 100Kb x level), which is too large for peeking.
        // So unfortunately we can recognise bzip2 but not the contents,
        // except that \x1772... magic indicates the stream is empty.
        if (s[4] == '\x31') return 0;
        else len = 0;
    }
    else if (len >= 6 && memcmp(s, "\xfd""7zXZ\0", 6) == 0) {
        fmt->compression = xz_compression;
#ifdef HAVE_LIBLZMA
        len = decompress_peek_xz(hfile, s, sizeof s);
#else
        // Without liblzma, we can't recognise the decompressed contents.
        return 0;
#endif
    }
    else if (len >= 4 && memcmp(s, "\x28\xb5\x2f\xfd", 4) == 0) {
        fmt->compression = zstd_compression;
        return 0;
    }
    else {
        len = hpeek(hfile, s, sizeof s);
    }
    if (len < 0) return -1;

    if (len == 0) {
        fmt->format = empty_format;
        return 0;
    }

    // We avoid using filename extensions wherever possible (as filenames are
    // not always available), but in a few cases they must be considered:
    //  - FASTA/Q indexes are simply tab-separated text; files that match these
    //    patterns but not the fai/fqi extension are usually generic BED files
    //  - GZI indexes have no magic numbers so can only be detected by filename
    if (fname && strcmp(fname, "-") != 0) {
        char *s;
        if (find_file_extension(fname, extension) < 0) extension[0] = '\0';
        for (s = extension; *s; s++) *s = tolower_c(*s);
    }
    else extension[0] = '\0';

    if (len >= 6 && memcmp(s,"CRAM",4) == 0 && s[4]>=1 && s[4]<=7 && s[5]<=7) {
        fmt->category = sequence_data;
        fmt->format = cram;
        fmt->version.major = s[4], fmt->version.minor = s[5];
        fmt->compression = custom;
        return 0;
    }
    else if (len >= 4 && s[3] <= '\4') {
        if (memcmp(s, "BAM\1", 4) == 0) {
            fmt->category = sequence_data;
            fmt->format = bam;
            // TODO Decompress enough to pick version from @HD-VN header
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BAI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = bai;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\4", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\2", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = s[3];
            fmt->version.minor = (len >= 5 && s[4] <= 2)? s[4] : 0;
            return 0;
        }
        else if (memcmp(s, "CSI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = csi;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "TBI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = tbi;
            return 0;
        }
        // GZI indexes have no magic numbers, so must be recognised solely by
        // filename extension.
        else if (strcmp(extension, "gzi") == 0) {
            fmt->category = index_file;
            fmt->format = gzi;
            return 0;
        }
    }
    else if (len >= 16 && memcmp(s, "##fileformat=VCF", 16) == 0) {
        fmt->category = variant_data;
        fmt->format = vcf;
        if (len >= 21 && s[16] == 'v')
            parse_version(fmt, &s[17], &s[len]);
        return 0;
    }
    else if (len >= 4 && s[0] == '@' &&
             (memcmp(s, "@HD\t", 4) == 0 || memcmp(s, "@SQ\t", 4) == 0 ||
              memcmp(s, "@RG\t", 4) == 0 || memcmp(s, "@PG\t", 4) == 0 ||
              memcmp(s, "@CO\t", 4) == 0)) {
        fmt->category = sequence_data;
        fmt->format = sam;
        // @HD-VN is not guaranteed to be the first tag, but then @HD is
        // not guaranteed to be present at all...
        if (len >= 9 && memcmp(s, "@HD\tVN:", 7) == 0)
            parse_version(fmt, &s[7], &s[len]);
        else
            fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }
    else if (len >= 8 && memcmp(s, "d4\xdd\xdd", 4) == 0) {
        fmt->category = region_list;
        fmt->format = d4_format;
        // How to decode the D4 Format Version bytes is not yet specified
        // so we don't try to set fmt->version.{major,minor}.
        return 0;
    }
    else if (cmp_nonblank("{\"htsget\":", s, &s[len]) == 0) {
        fmt->category = unknown_category;
        fmt->format = htsget;
        return 0;
    }
    else if (len > 8 && memcmp(s, "crypt4gh", 8) == 0) {
        fmt->category = unknown_category;
        fmt->format = hts_crypt4gh_format;
        return 0;
    }
    else if (len >= 1 && s[0] == '>' && is_fastaq(s, &s[len])) {
        fmt->category = sequence_data;
        fmt->format = fasta_format;
        return 0;
    }
    else if (len >= 1 && s[0] == '@' && is_fastaq(s, &s[len])) {
        fmt->category = sequence_data;
        fmt->format = fastq_format;
        return 0;
    }
    else if (parse_tabbed_text(columns, sizeof columns, s,
                               &s[len], &complete) > 0) {
        // A complete SAM line is at least 11 columns.  On unmapped long reads may
        // be missing two.  (On mapped long reads we must have an @ header so long
        // CIGAR is irrelevant.)
        if (colmatch(columns, "ZiZiiCZiiZZOOOOOOOOOOOOOOOOOOOO+")
            >= 9 + 2*complete) {
            fmt->category = sequence_data;
            fmt->format = sam;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (fmt->compression == gzip && colmatch(columns, "iiiiii") == 6) {
            fmt->category = index_file;
            fmt->format = crai;
            return 0;
        }
        else if (strstr(extension, "fqi") && colmatch(columns, "Ziiiii") == 6) {
            fmt->category = index_file;
            fmt->format = fqi_format;
            return 0;
        }
        else if (strstr(extension, "fai") && colmatch(columns, "Ziiii") == 5) {
            fmt->category = index_file;
            fmt->format = fai_format;
            return 0;
        }
        else if (colmatch(columns, "Zii+") >= 3) {
            fmt->category = region_list;
            fmt->format = bed;
            return 0;
        }
    }

    // Arbitrary text files can be read using hts_getline().
    if (is_text_only(s, &s[len])) fmt->format = text_format;

    // Nothing recognised: leave unset fmt-> fields as unknown.
    return 0;
}

char *hts_format_description(const htsFormat *format)
{
    kstring_t str = { 0, 0, NULL };

    switch (format->format) {
    case sam:   kputs("SAM", &str); break;
    case bam:   kputs("BAM", &str); break;
    case cram:  kputs("CRAM", &str); break;
    case fasta_format:  kputs("FASTA", &str); break;
    case fastq_format:  kputs("FASTQ", &str); break;
    case vcf:   kputs("VCF", &str); break;
    case bcf:
        if (format->version.major == 1) kputs("Legacy BCF", &str);
        else kputs("BCF", &str);
        break;
    case bai:   kputs("BAI", &str); break;
    case crai:  kputs("CRAI", &str); break;
    case csi:   kputs("CSI", &str); break;
    case fai_format:    kputs("FASTA-IDX", &str); break;
    case fqi_format:    kputs("FASTQ-IDX", &str); break;
    case gzi:   kputs("GZI", &str); break;
    case tbi:   kputs("Tabix", &str); break;
    case bed:   kputs("BED", &str); break;
    case d4_format:     kputs("D4", &str); break;
    case htsget: kputs("htsget", &str); break;
    case hts_crypt4gh_format: kputs("crypt4gh", &str); break;
    case empty_format:  kputs("empty", &str); break;
    default:    kputs("unknown", &str); break;
    }

    if (format->version.major >= 0) {
        kputs(" version ", &str);
        kputw(format->version.major, &str);
        if (format->version.minor >= 0) {
            kputc('.', &str);
            kputw(format->version.minor, &str);
        }
    }

    switch (format->compression) {
    case bzip2_compression:  kputs(" bzip2-compressed", &str); break;
    case razf_compression:   kputs(" legacy-RAZF-compressed", &str); break;
    case xz_compression:     kputs(" XZ-compressed", &str); break;
    case zstd_compression:   kputs(" Zstandard-compressed", &str); break;
    case custom: kputs(" compressed", &str); break;
    case gzip:   kputs(" gzip-compressed", &str); break;
    case bgzf:
        switch (format->format) {
        case bam:
        case bcf:
        case csi:
        case tbi:
            // These are by definition BGZF, so just use the generic term
            kputs(" compressed", &str);
            break;
        default:
            kputs(" BGZF-compressed", &str);
            break;
        }
        break;
    default: break;
    }

    switch (format->category) {
    case sequence_data: kputs(" sequence", &str); break;
    case variant_data:  kputs(" variant calling", &str); break;
    case index_file:    kputs(" index", &str); break;
    case region_list:   kputs(" genomic region", &str); break;
    default: break;
    }

    if (format->compression == no_compression)
        switch (format->format) {
        case text_format:
        case sam:
        case crai:
        case vcf:
        case bed:
        case fai_format:
        case fqi_format:
        case fasta_format:
        case fastq_format:
        case htsget:
            kputs(" text", &str);
            break;

        case empty_format:
            break;

        default:
            kputs(" data", &str);
            break;
        }
    else
        kputs(" data", &str);

    return ks_release(&str);
}

htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)
{
    char smode[101], *cp, *cp2, *mode_c;
    htsFile *fp = NULL;
    hFILE *hfile = NULL;
    char fmt_code = '\0';
    // see enum htsExactFormat in htslib/hts.h
    const char format_to_mode[] = "\0g\0\0b\0c\0\0b\0g\0\0\0\0\0Ff\0\0";

    strncpy(smode, mode, 99);
    smode[99]=0;
    if ((cp = strchr(smode, ',')))
        *cp = '\0';

    // Migrate format code (b or c) to the end of the smode buffer.
    for (cp2 = cp = smode; *cp; cp++) {
        if (*cp == 'b')
            fmt_code = 'b';
        else if (*cp == 'c')
            fmt_code = 'c';
        else
            *cp2++ = *cp;
    }
    mode_c = cp2;
    *cp2++ = fmt_code;
    *cp2++ = 0;

    // Set or reset the format code if opts->format is used
    if (fmt && fmt->format > unknown_format
        && fmt->format < sizeof(format_to_mode)) {
        *mode_c = format_to_mode[fmt->format];
    }

    // If we really asked for a compressed text format then mode_c above will
    // point to nul.  We set to 'z' to enable bgzf.
    if (strchr(mode, 'w') && fmt && fmt->compression == bgzf) {
        if (fmt->format == sam || fmt->format == vcf || fmt->format == text_format)
            *mode_c = 'z';
    }

    char *rmme = NULL, *fnidx = strstr(fn, HTS_IDX_DELIM);
    if ( fnidx ) {
        rmme = strdup(fn);
        if ( !rmme ) goto error;
        rmme[fnidx-fn] = 0;
        fn = rmme;
    }

    hfile = hopen(fn, smode);
    if (hfile == NULL) goto error;

    fp = hts_hopen(hfile, fn, smode);
    if (fp == NULL) goto error;

    // Compensate for the loss of exactness in htsExactFormat.
    // hts_hopen returns generics such as binary or text, but we
    // have been given something explicit here so use that instead.
    if (fp->is_write && fmt &&
        (fmt->format == bam || fmt->format == sam ||
         fmt->format == vcf || fmt->format == bcf ||
         fmt->format == bed || fmt->format == fasta_format ||
         fmt->format == fastq_format))
        fp->format.format = fmt->format;

    if (fmt && fmt->specific)
        if (hts_opt_apply(fp, fmt->specific) != 0)
            goto error;

    if ( rmme ) free(rmme);
    return fp;

error:
    hts_log_error("Failed to open file \"%s\"%s%s", fn,
                  errno ? " : " : "", errno ? strerror(errno) : "");
    if ( rmme ) free(rmme);

    if (hfile)
        hclose_abruptly(hfile);

    return NULL;
}

htsFile *hts_open(const char *fn, const char *mode) {
    return hts_open_format(fn, mode, NULL);
}

/*
 * Splits str into a prefix, delimiter ('\0' or delim), and suffix, writing
 * the prefix in lowercase into buf and returning a pointer to the suffix.
 * On return, buf is always NUL-terminated; thus assumes that the "keyword"
 * prefix should be one of several known values of maximum length buflen-2.
 * (If delim is not found, returns a pointer to the '\0'.)
 */
static const char *
scan_keyword(const char *str, char delim, char *buf, size_t buflen)
{
    size_t i = 0;
    while (*str && *str != delim) {
        if (i < buflen-1) buf[i++] = tolower_c(*str);
        str++;
    }

    buf[i] = '\0';
    return *str? str+1 : str;
}

/*
 * Parses arg and appends it to the option list.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int hts_opt_add(hts_opt **opts, const char *c_arg) {
    hts_opt *o, *t;
    char *val;

    /*
     * IMPORTANT!!!
     * If you add another string option here, don't forget to also add
     * it to the case statement in hts_opt_apply.
     */

    if (!c_arg)
        return -1;

    if (!(o =  malloc(sizeof(*o))))
        return -1;

    if (!(o->arg = strdup(c_arg))) {
        free(o);
        return -1;
    }

    if (!(val = strchr(o->arg, '=')))
        val = "1"; // assume boolean
    else
        *val++ = '\0';

    if (strcmp(o->arg, "decode_md") == 0 ||
        strcmp(o->arg, "DECODE_MD") == 0)
        o->opt = CRAM_OPT_DECODE_MD, o->val.i = atoi(val);

    else if (strcmp(o->arg, "verbosity") == 0 ||
             strcmp(o->arg, "VERBOSITY") == 0)
        o->opt = CRAM_OPT_VERBOSITY, o->val.i = atoi(val);

    else if (strcmp(o->arg, "seqs_per_slice") == 0 ||
             strcmp(o->arg, "SEQS_PER_SLICE") == 0)
        o->opt = CRAM_OPT_SEQS_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "bases_per_slice") == 0 ||
             strcmp(o->arg, "BASES_PER_SLICE") == 0)
        o->opt = CRAM_OPT_BASES_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "slices_per_container") == 0 ||
             strcmp(o->arg, "SLICES_PER_CONTAINER") == 0)
        o->opt = CRAM_OPT_SLICES_PER_CONTAINER, o->val.i = atoi(val);

    else if (strcmp(o->arg, "embed_ref") == 0 ||
             strcmp(o->arg, "EMBED_REF") == 0)
        o->opt = CRAM_OPT_EMBED_REF, o->val.i = atoi(val);

    else if (strcmp(o->arg, "no_ref") == 0 ||
             strcmp(o->arg, "NO_REF") == 0)
        o->opt = CRAM_OPT_NO_REF, o->val.i = atoi(val);

    else if (strcmp(o->arg, "pos_delta") == 0 ||
             strcmp(o->arg, "POS_DELTA") == 0)
        o->opt = CRAM_OPT_POS_DELTA, o->val.i = atoi(val);

    else if (strcmp(o->arg, "ignore_md5") == 0 ||
             strcmp(o->arg, "IGNORE_MD5") == 0)
        o->opt = CRAM_OPT_IGNORE_MD5, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_bzip2") == 0 ||
             strcmp(o->arg, "USE_BZIP2") == 0)
        o->opt = CRAM_OPT_USE_BZIP2, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_rans") == 0 ||
             strcmp(o->arg, "USE_RANS") == 0)
        o->opt = CRAM_OPT_USE_RANS, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_lzma") == 0 ||
             strcmp(o->arg, "USE_LZMA") == 0)
        o->opt = CRAM_OPT_USE_LZMA, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_tok") == 0 ||
             strcmp(o->arg, "USE_TOK") == 0)
        o->opt = CRAM_OPT_USE_TOK, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_fqz") == 0 ||
             strcmp(o->arg, "USE_FQZ") == 0)
        o->opt = CRAM_OPT_USE_FQZ, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_arith") == 0 ||
             strcmp(o->arg, "USE_ARITH") == 0)
        o->opt = CRAM_OPT_USE_ARITH, o->val.i = atoi(val);

    else if (strcmp(o->arg, "fast") == 0 ||
             strcmp(o->arg, "FAST") == 0)
        o->opt = HTS_OPT_PROFILE, o->val.i = HTS_PROFILE_FAST;

    else if (strcmp(o->arg, "normal") == 0 ||
             strcmp(o->arg, "NORMAL") == 0)
        o->opt = HTS_OPT_PROFILE, o->val.i = HTS_PROFILE_NORMAL;

    else if (strcmp(o->arg, "small") == 0 ||
             strcmp(o->arg, "SMALL") == 0)
        o->opt = HTS_OPT_PROFILE, o->val.i = HTS_PROFILE_SMALL;

    else if (strcmp(o->arg, "archive") == 0 ||
             strcmp(o->arg, "ARCHIVE") == 0)
        o->opt = HTS_OPT_PROFILE, o->val.i = HTS_PROFILE_ARCHIVE;

    else if (strcmp(o->arg, "reference") == 0 ||
             strcmp(o->arg, "REFERENCE") == 0)
        o->opt = CRAM_OPT_REFERENCE, o->val.s = val;

    else if (strcmp(o->arg, "version") == 0 ||
             strcmp(o->arg, "VERSION") == 0)
        o->opt = CRAM_OPT_VERSION, o->val.s =val;

    else if (strcmp(o->arg, "multi_seq_per_slice") == 0 ||
             strcmp(o->arg, "MULTI_SEQ_PER_SLICE") == 0)
        o->opt = CRAM_OPT_MULTI_SEQ_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "nthreads") == 0 ||
             strcmp(o->arg, "NTHREADS") == 0)
        o->opt = HTS_OPT_NTHREADS, o->val.i = atoi(val);

    else if (strcmp(o->arg, "cache_size") == 0 ||
             strcmp(o->arg, "CACHE_SIZE") == 0) {
        char *endp;
        o->opt = HTS_OPT_CACHE_SIZE;
        o->val.i = strtol(val, &endp, 0);
        // NB: Doesn't support floats, eg 1.5g
        // TODO: extend hts_parse_decimal? See also samtools sort.
        switch (*endp) {
        case 'g': case 'G': o->val.i *= 1024; // fall through
        case 'm': case 'M': o->val.i *= 1024; // fall through
        case 'k': case 'K': o->val.i *= 1024; break;
        case '\0': break;
        default:
            hts_log_error("Unrecognised cache size suffix '%c'", *endp);
            free(o->arg);
            free(o);
            return -1;
        }
    }

    else if (strcmp(o->arg, "required_fields") == 0 ||
             strcmp(o->arg, "REQUIRED_FIELDS") == 0)
        o->opt = CRAM_OPT_REQUIRED_FIELDS, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "lossy_names") == 0 ||
             strcmp(o->arg, "LOSSY_NAMES") == 0)
        o->opt = CRAM_OPT_LOSSY_NAMES, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "name_prefix") == 0 ||
             strcmp(o->arg, "NAME_PREFIX") == 0)
        o->opt = CRAM_OPT_PREFIX, o->val.s = val;

    else if (strcmp(o->arg, "store_md") == 0 ||
             strcmp(o->arg, "store_md") == 0)
        o->opt = CRAM_OPT_STORE_MD, o->val.i = atoi(val);

    else if (strcmp(o->arg, "store_nm") == 0 ||
             strcmp(o->arg, "store_nm") == 0)
        o->opt = CRAM_OPT_STORE_NM, o->val.i = atoi(val);

    else if (strcmp(o->arg, "block_size") == 0 ||
             strcmp(o->arg, "BLOCK_SIZE") == 0)
        o->opt = HTS_OPT_BLOCK_SIZE, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "level") == 0 ||
             strcmp(o->arg, "LEVEL") == 0)
        o->opt = HTS_OPT_COMPRESSION_LEVEL, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "filter") == 0 ||
             strcmp(o->arg, "FILTER") == 0)
        o->opt = HTS_OPT_FILTER, o->val.s = val;

    else if (strcmp(o->arg, "fastq_aux") == 0 ||
        strcmp(o->arg, "FASTQ_AUX") == 0)
        o->opt = FASTQ_OPT_AUX, o->val.s = val;

    else if (strcmp(o->arg, "fastq_barcode") == 0 ||
        strcmp(o->arg, "FASTQ_BARCODE") == 0)
        o->opt = FASTQ_OPT_BARCODE, o->val.s = val;

    else if (strcmp(o->arg, "fastq_rnum") == 0 ||
        strcmp(o->arg, "FASTQ_RNUM") == 0)
        o->opt = FASTQ_OPT_RNUM, o->val.i = 1;

    else if (strcmp(o->arg, "fastq_casava") == 0 ||
        strcmp(o->arg, "FASTQ_CASAVA") == 0)
        o->opt = FASTQ_OPT_CASAVA, o->val.i = 1;

    else if (strcmp(o->arg, "fastq_name2") == 0 ||
        strcmp(o->arg, "FASTQ_NAME2") == 0)
        o->opt = FASTQ_OPT_NAME2, o->val.i = 1;

    else {
        hts_log_error("Unknown option '%s'", o->arg);
        free(o->arg);
        free(o);
        return -1;
    }

    o->next = NULL;

    // Append; assumes small list.
    if (*opts) {
        t = *opts;
        while (t->next)
            t = t->next;
        t->next = o;
    } else {
        *opts = o;
    }

    return 0;
}

/*
 * Applies an hts_opt option list to a given htsFile.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int hts_opt_apply(htsFile *fp, hts_opt *opts) {
    hts_opt *last = NULL;

    for (; opts;  opts = (last=opts)->next) {
        switch (opts->opt) {
            case CRAM_OPT_REFERENCE:
                if (!(fp->fn_aux = strdup(opts->val.s)))
                    return -1;
                // fall through
            case CRAM_OPT_VERSION:
            case CRAM_OPT_PREFIX:
            case HTS_OPT_FILTER:
            case FASTQ_OPT_AUX:
            case FASTQ_OPT_BARCODE:
                if (hts_set_opt(fp,  opts->opt,  opts->val.s) != 0)
                    return -1;
                break;
            default:
                if (hts_set_opt(fp,  opts->opt,  opts->val.i) != 0)
                    return -1;
                break;
        }
    }

    return 0;
}

/*
 * Frees an hts_opt list.
 */
void hts_opt_free(hts_opt *opts) {
    hts_opt *last = NULL;
    while (opts) {
        opts = (last=opts)->next;
        free(last->arg);
        free(last);
    }
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_opt_list(htsFormat *fmt, const char *str) {
    while (str && *str) {
        const char *str_start;
        int len;
        char arg[8001];

        while (*str && *str == ',')
            str++;

        for (str_start = str; *str && *str != ','; str++);
        len = str - str_start;

        // Produce a nul terminated copy of the option
        strncpy(arg, str_start, len < 8000 ? len : 8000);
        arg[len < 8000 ? len : 8000] = '\0';

        if (hts_opt_add((hts_opt **)&fmt->specific, arg) != 0)
            return -1;

        if (*str)
            str++;
    }

    return 0;
}

/*
 * Accepts a string file format (sam, bam, cram, vcf, bam) optionally
 * followed by a comma separated list of key=value options and splits
 * these up into the fields of htsFormat struct.
 *
 * format is assumed to be already initialised, either to blank
 * "unknown" values or via previous hts_opt_add calls.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_format(htsFormat *format, const char *str) {
    char fmt[8];
    const char *cp = scan_keyword(str, ',', fmt, sizeof fmt);

    format->version.minor = 0; // unknown
    format->version.major = 0; // unknown

    if (strcmp(fmt, "sam") == 0) {
        format->category          = sequence_data;
        format->format            = sam;
        format->compression       = no_compression;
        format->compression_level = 0;
    } else if (strcmp(fmt, "sam.gz") == 0) {
        format->category          = sequence_data;
        format->format            = sam;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else if (strcmp(fmt, "bam") == 0) {
        format->category          = sequence_data;
        format->format            = bam;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else if (strcmp(fmt, "cram") == 0) {
        format->category          = sequence_data;
        format->format            = cram;
        format->compression       = custom;
        format->compression_level = -1;
    } else if (strcmp(fmt, "vcf") == 0) {
        format->category          = variant_data;
        format->format            = vcf;
        format->compression       = no_compression;
        format->compression_level = 0;
    } else if (strcmp(fmt, "bcf") == 0) {
        format->category          = variant_data;
        format->format            = bcf;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else if (strcmp(fmt, "fastq") == 0 || strcmp(fmt, "fq") == 0) {
        format->category          = sequence_data;
        format->format            = fastq_format;
        format->compression       = no_compression;
        format->compression_level = 0;
    } else if (strcmp(fmt, "fastq.gz") == 0 || strcmp(fmt, "fq.gz") == 0) {
        format->category          = sequence_data;
        format->format            = fastq_format;
        format->compression       = bgzf;
        format->compression_level = 0;
    } else if (strcmp(fmt, "fasta") == 0 || strcmp(fmt, "fa") == 0) {
        format->category          = sequence_data;
        format->format            = fasta_format;
        format->compression       = no_compression;
        format->compression_level = 0;
    } else if (strcmp(fmt, "fasta.gz") == 0 || strcmp(fmt, "fa.gz") == 0) {
        format->category          = sequence_data;
        format->format            = fasta_format;
        format->compression       = bgzf;
        format->compression_level = 0;
    } else {
        return -1;
    }

    return hts_parse_opt_list(format, cp);
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int hts_process_opts(htsFile *fp, const char *opts) {
    htsFormat fmt;

    fmt.specific = NULL;
    if (hts_parse_opt_list(&fmt, opts) != 0)
        return -1;

    if (hts_opt_apply(fp, fmt.specific) != 0) {
        hts_opt_free(fmt.specific);
        return -1;
    }

    hts_opt_free(fmt.specific);

    return 0;
}

static int hts_crypt4gh_redirect(const char *fn, const char *mode,
                                 hFILE **hfile_ptr, htsFile *fp) {
    hFILE *hfile1 = *hfile_ptr;
    hFILE *hfile2 = NULL;
    char fn_buf[512], *fn2 = fn_buf;
    const char *prefix = "crypt4gh:";
    size_t fn2_len = strlen(prefix) + strlen(fn) + 1;
    int ret = -1;

    if (fn2_len > sizeof(fn_buf)) {
        if (fn2_len >= INT_MAX) // Silence gcc format-truncation warning
            return -1;
        fn2 = malloc(fn2_len);
        if (!fn2) return -1;
    }

    // Reopen fn using the crypt4gh plug-in (if available)
    snprintf(fn2, fn2_len, "%s%s", prefix, fn);
    hfile2 = hopen(fn2, mode, "parent", hfile1, NULL);
    if (hfile2) {
        // Replace original hfile with the new one.  The original is now
        // enclosed within hfile2
        *hfile_ptr = hfile2;
        ret = 0;
    }

    if (fn2 != fn_buf)
        free(fn2);
    return ret;
}

htsFile *hts_hopen(hFILE *hfile, const char *fn, const char *mode)
{
    hFILE *hfile_orig = hfile;
    htsFile *fp = (htsFile*)calloc(1, sizeof(htsFile));
    char simple_mode[101], *cp, *opts;
    simple_mode[100] = '\0';

    if (fp == NULL) goto error;

    fp->fn = strdup(fn);
    fp->is_be = ed_is_big();

    // Split mode into simple_mode,opts strings
    if ((cp = strchr(mode, ','))) {
        strncpy(simple_mode, mode, cp-mode <= 100 ? cp-mode : 100);
        simple_mode[cp-mode] = '\0';
        opts = cp+1;
    } else {
        strncpy(simple_mode, mode, 100);
        opts = NULL;
    }

    if (strchr(simple_mode, 'r')) {
        const int max_loops = 5; // Should be plenty
        int loops = 0;
        if (hts_detect_format2(hfile, fn, &fp->format) < 0) goto error;

        // Deal with formats that re-direct an underlying file via a plug-in.
        // Loops as we may have crypt4gh served via htsget, or
        // crypt4gh-in-crypt4gh.
        while (fp->format.format == htsget ||
               fp->format.format == hts_crypt4gh_format) {
            // Ensure we don't get stuck in an endless redirect loop
            if (++loops > max_loops) {
                errno = ELOOP;
                goto error;
            }

            if (fp->format.format == htsget) {
                hFILE *hfile2 = hopen_htsget_redirect(hfile, simple_mode);
                if (hfile2 == NULL) goto error;

                hfile = hfile2;
            }
            else if (fp->format.format == hts_crypt4gh_format) {
                if (hts_crypt4gh_redirect(fn, simple_mode, &hfile, fp) < 0)
                    goto error;
            }

            // Re-detect format against the result of the redirection
            if (hts_detect_format2(hfile, fn, &fp->format) < 0) goto error;
        }
    }
    else if (strchr(simple_mode, 'w') || strchr(simple_mode, 'a')) {
        htsFormat *fmt = &fp->format;
        fp->is_write = 1;

        if (strchr(simple_mode, 'b')) fmt->format = binary_format;
        else if (strchr(simple_mode, 'c')) fmt->format = cram;
        else if (strchr(simple_mode, 'f')) fmt->format = fastq_format;
        else if (strchr(simple_mode, 'F')) fmt->format = fasta_format;
        else fmt->format = text_format;

        if (strchr(simple_mode, 'z')) fmt->compression = bgzf;
        else if (strchr(simple_mode, 'g')) fmt->compression = gzip;
        else if (strchr(simple_mode, 'u')) fmt->compression = no_compression;
        else {
            // No compression mode specified, set to the default for the format
            switch (fmt->format) {
            case binary_format: fmt->compression = bgzf; break;
            case cram: fmt->compression = custom; break;
            case fastq_format: fmt->compression = no_compression; break;
            case fasta_format: fmt->compression = no_compression; break;
            case text_format: fmt->compression = no_compression; break;
            default: abort();
            }
        }

        // Fill in category (if determinable; e.g. 'b' could be BAM or BCF)
        fmt->category = format_category(fmt->format);

        fmt->version.major = fmt->version.minor = -1;
        fmt->compression_level = -1;
        fmt->specific = NULL;
    }
    else { errno = EINVAL; goto error; }

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
        if (fp->fp.bgzf == NULL) goto error;
        fp->is_bin = fp->is_bgzf = 1;
        break;

    case cram:
        fp->fp.cram = cram_dopen(hfile, fn, simple_mode);
        if (fp->fp.cram == NULL) goto error;
        if (!fp->is_write)
            cram_set_option(fp->fp.cram, CRAM_OPT_DECODE_MD, -1); // auto
        fp->is_cram = 1;
        break;

    case empty_format:
    case text_format:
    case bed:
    case fasta_format:
    case fastq_format:
    case sam:
    case vcf:
        if (fp->format.compression != no_compression) {
            fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
            if (fp->fp.bgzf == NULL) goto error;
            fp->is_bgzf = 1;
        }
        else
            fp->fp.hfile = hfile;
        break;

    default:
        errno = EFTYPE;
        goto error;
    }

    if (opts)
        hts_process_opts(fp, opts);

    // If redirecting, close the original hFILE now (pedantically we would
    // instead close it in hts_close(), but this a simplifying optimisation)
    if (hfile != hfile_orig) hclose_abruptly(hfile_orig);

    return fp;

error:
    hts_log_error("Failed to open file %s", fn);

    // If redirecting, close the failed redirection hFILE that we have opened
    if (hfile != hfile_orig) hclose_abruptly(hfile);

    if (fp) {
        free(fp->fn);
        free(fp->fn_aux);
        free(fp);
    }
    return NULL;
}

int hts_close(htsFile *fp)
{
    int ret = 0, save;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        ret = bgzf_close(fp->fp.bgzf);
        break;

    case cram:
        if (!fp->is_write) {
            switch (cram_eof(fp->fp.cram)) {
            case 2:
                hts_log_warning("EOF marker is absent. The input is probably truncated");
                break;
            case 0:  /* not at EOF, but may not have wanted all seqs */
            default: /* case 1, expected EOF */
                break;
            }
        }
        ret = cram_close(fp->fp.cram);
        break;

    case empty_format:
    case text_format:
    case bed:
    case fasta_format:
    case fastq_format:
    case sam:
    case vcf:
        if (fp->format.format == sam)
            ret = sam_state_destroy(fp);
        else if (fp->format.format == fastq_format ||
                 fp->format.format == fasta_format)
            fastq_state_destroy(fp);

        if (fp->format.compression != no_compression)
            ret |= bgzf_close(fp->fp.bgzf);
        else
            ret |= hclose(fp->fp.hfile);
        break;

    default:
        ret = -1;
        break;
    }

    save = errno;
    sam_hdr_destroy(fp->bam_header);
    hts_idx_destroy(fp->idx);
    hts_filter_free(fp->filter);
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    errno = save;
    return ret;
}

int hts_flush(htsFile *fp)
{
    if (fp == NULL) return 0;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        return bgzf_flush(fp->fp.bgzf);

    case cram:
        return cram_flush(fp->fp.cram);

    case empty_format:
    case text_format:
    case bed:
    case fasta_format:
    case fastq_format:
    case sam:
    case vcf:
        if (fp->format.compression != no_compression)
            return bgzf_flush(fp->fp.bgzf);
        else
            return hflush(fp->fp.hfile);

    default:
        break;
    }

    return 0;
}

const htsFormat *hts_get_format(htsFile *fp)
{
    return fp? &fp->format : NULL;
}

const char *hts_format_file_extension(const htsFormat *format) {
    if (!format)
        return "?";

    switch (format->format) {
    case sam:  return "sam";
    case bam:  return "bam";
    case bai:  return "bai";
    case cram: return "cram";
    case crai: return "crai";
    case vcf:  return "vcf";
    case bcf:  return "bcf";
    case csi:  return "csi";
    case fai_format:   return "fai";
    case fqi_format:   return "fqi";
    case gzi:  return "gzi";
    case tbi:  return "tbi";
    case bed:  return "bed";
    case d4_format:    return "d4";
    case fasta_format: return "fa";
    case fastq_format: return "fq";
    default:   return "?";
    }
}

static hFILE *hts_hfile(htsFile *fp) {
    switch (fp->format.format) {
    case binary_format:// fall through
    case bcf:          // fall through
    case bam:          return bgzf_hfile(fp->fp.bgzf);
    case cram:         return cram_hfile(fp->fp.cram);
    case text_format:  return fp->fp.hfile;
    case vcf:          // fall through
    case fastq_format: // fall through
    case fasta_format: // fall through
    case sam:          return fp->format.compression != no_compression
                              ? bgzf_hfile(fp->fp.bgzf)
                              : fp->fp.hfile;
    default:           return NULL;
    }
}

int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    switch (opt) {
    case HTS_OPT_NTHREADS: {
        va_start(args, opt);
        int nthreads = va_arg(args, int);
        va_end(args);
        return hts_set_threads(fp, nthreads);
    }

    case HTS_OPT_BLOCK_SIZE: {
        hFILE *hf = hts_hfile(fp);

        if (hf) {
            va_start(args, opt);
            if (hfile_set_blksize(hf, va_arg(args, int)) != 0)
                hts_log_warning("Failed to change block size");
            va_end(args);
        }
        else {
            // To do - implement for vcf/bcf.
            hts_log_warning("Cannot change block size for this format");
        }

        return 0;
    }

    case HTS_OPT_THREAD_POOL: {
        va_start(args, opt);
        htsThreadPool *p = va_arg(args, htsThreadPool *);
        va_end(args);
        return hts_set_thread_pool(fp, p);
    }

    case HTS_OPT_CACHE_SIZE: {
        va_start(args, opt);
        int cache_size = va_arg(args, int);
        va_end(args);
        hts_set_cache_size(fp, cache_size);
        return 0;
    }

    case FASTQ_OPT_CASAVA:
    case FASTQ_OPT_RNUM:
    case FASTQ_OPT_NAME2:
        if (fp->format.format == fastq_format ||
            fp->format.format == fasta_format)
            return fastq_state_set(fp, opt);
        return 0;

    case FASTQ_OPT_AUX:
        if (fp->format.format == fastq_format ||
            fp->format.format == fasta_format) {
            va_start(args, opt);
            char *list = va_arg(args, char *);
            va_end(args);
            return fastq_state_set(fp, opt, list);
        }
        return 0;

    case FASTQ_OPT_BARCODE:
        if (fp->format.format == fastq_format ||
            fp->format.format == fasta_format) {
            va_start(args, opt);
            char *bc = va_arg(args, char *);
            va_end(args);
            return fastq_state_set(fp, opt, bc);
        }
        return 0;

    // Options below here flow through to cram_set_voption
    case HTS_OPT_COMPRESSION_LEVEL: {
        va_start(args, opt);
        int level = va_arg(args, int);
        va_end(args);
        if (fp->is_bgzf)
            fp->fp.bgzf->compress_level = level;
        else if (fp->format.format == cram)
            return cram_set_option(fp->fp.cram, opt, level);
        return 0;
    }

    case HTS_OPT_FILTER: {
        va_start(args, opt);
        char *expr = va_arg(args, char *);
        va_end(args);
        return hts_set_filter_expression(fp, expr);
    }

    case HTS_OPT_PROFILE: {
        va_start(args, opt);
        enum hts_profile_option prof = va_arg(args, int);
        va_end(args);
        if (fp->is_bgzf) {
            switch (prof) {
#ifdef HAVE_LIBDEFLATE
            case HTS_PROFILE_FAST:    fp->fp.bgzf->compress_level =  2; break;
            case HTS_PROFILE_NORMAL:  fp->fp.bgzf->compress_level = -1; break;
            case HTS_PROFILE_SMALL:   fp->fp.bgzf->compress_level = 10; break;
            case HTS_PROFILE_ARCHIVE: fp->fp.bgzf->compress_level = 12; break;
#else
            case HTS_PROFILE_FAST:    fp->fp.bgzf->compress_level =  1; break;
            case HTS_PROFILE_NORMAL:  fp->fp.bgzf->compress_level = -1; break;
            case HTS_PROFILE_SMALL:   fp->fp.bgzf->compress_level =  8; break;
            case HTS_PROFILE_ARCHIVE: fp->fp.bgzf->compress_level =  9; break;
#endif
            }
        } // else CRAM manages this in its own way
        break;
    }

    default:
        break;
    }

    if (fp->format.format != cram)
        return 0;

    va_start(args, opt);
    r = cram_set_voption(fp->fp.cram, opt, args);
    va_end(args);

    return r;
}

BGZF *hts_get_bgzfp(htsFile *fp);

int hts_set_threads(htsFile *fp, int n)
{
    if (fp->format.format == sam) {
        return sam_set_threads(fp, n);
    } else if (fp->format.compression == bgzf) {
        return bgzf_mt(hts_get_bgzfp(fp), n, 256/*unused*/);
    } else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_NTHREADS, n);
    }
    else return 0;
}

int hts_set_thread_pool(htsFile *fp, htsThreadPool *p) {
    if (fp->format.format == sam || fp->format.format == text_format) {
        return sam_set_thread_pool(fp, p);
    } else if (fp->format.compression == bgzf) {
        return bgzf_thread_pool(hts_get_bgzfp(fp), p->pool, p->qsize);
    } else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_THREAD_POOL, p);
    }
    else return 0;
}

void hts_set_cache_size(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf)
        bgzf_set_cache_size(hts_get_bgzfp(fp), n);
}

int hts_set_fai_filename(htsFile *fp, const char *fn_aux)
{
    free(fp->fn_aux);
    if (fn_aux) {
        fp->fn_aux = strdup(fn_aux);
        if (fp->fn_aux == NULL) return -1;
    }
    else fp->fn_aux = NULL;

    if (fp->format.format == cram)
        if (cram_set_option(fp->fp.cram, CRAM_OPT_REFERENCE, fp->fn_aux))
            return -1;

    return 0;
}

int hts_set_filter_expression(htsFile *fp, const char *expr)
{
    if (fp->filter)
        hts_filter_free(fp->filter);

    if (!expr)
        return 0;

    return (fp->filter = hts_filter_init(expr))
        ? 0 : -1;
}

hFILE *hts_open_tmpfile(const char *fname, const char *mode, kstring_t *tmpname)
{
    int pid = (int) getpid();
    unsigned ptr = (uintptr_t) tmpname;
    int n = 0;
    hFILE *fp = NULL;

    do {
        // Attempt to further uniquify the temporary filename
        unsigned t = ((unsigned) time(NULL)) ^ ((unsigned) clock()) ^ ptr;
        n++;

        ks_clear(tmpname);
        if (ksprintf(tmpname, "%s.tmp_%d_%d_%u", fname, pid, n, t) < 0) break;

        fp = hopen(tmpname->s, mode);
    } while (fp == NULL && errno == EEXIST && n < 100);

    return fp;
}

// For VCF/BCF backward sweeper. Not exposing these functions because their
// future is uncertain. Things will probably have to change with hFILE...
BGZF *hts_get_bgzfp(htsFile *fp)
{
    if (fp->is_bgzf)
        return fp->fp.bgzf;
    else
        return NULL;
}
int hts_useek(htsFile *fp, off_t uoffset, int where)
{
    if (fp->is_bgzf)
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
        return (hseek(fp->fp.hfile, uoffset, SEEK_SET) >= 0)? 0 : -1;
}
off_t hts_utell(htsFile *fp)
{
    if (fp->is_bgzf)
        return bgzf_utell(fp->fp.bgzf);
    else
        return htell(fp->fp.hfile);
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
    int ret;
    if (! (delimiter == KS_SEP_LINE || delimiter == '\n')) {
        hts_log_error("Unexpected delimiter %d", delimiter);
        abort();
    }

    switch (fp->format.compression) {
    case no_compression:
        str->l = 0;
        ret = kgetline2(str, (kgets_func2 *) hgetln, fp->fp.hfile);
        if (ret >= 0) ret = (str->l <= INT_MAX)? (int) str->l : INT_MAX;
        else if (herrno(fp->fp.hfile)) ret = -2, errno = herrno(fp->fp.hfile);
        else ret = -1;
        break;

    case gzip:
    case bgzf:
        ret = bgzf_getline(fp->fp.bgzf, '\n', str);
        break;

    default:
        abort();
    }

    ++fp->lineno;
    return ret;
}

char **hts_readlist(const char *string, int is_file, int *_n)
{
    unsigned int m = 0, n = 0;
    char **s = 0, **s_new;
    if ( is_file )
    {
        BGZF *fp = bgzf_open(string, "r");
        if ( !fp ) return NULL;

        kstring_t str;
        int ret;
        str.s = 0; str.l = str.m = 0;
        while ((ret = bgzf_getline(fp, '\n', &str)) >= 0)
        {
            if (str.l == 0) continue;
            if (hts_resize(char*, n + 1, &m, &s, 0) < 0)
                goto err;
            s[n] = strdup(str.s);
            if (!s[n])
                goto err;
            n++;
        }
        if (ret < -1) // Read error
            goto err;
        bgzf_close(fp);
        free(str.s);
    }
    else
    {
        const char *q = string, *p = string;
        while ( 1 )
        {
            if (*p == ',' || *p == 0)
            {
                if (hts_resize(char*, n + 1, &m, &s, 0) < 0)
                    goto err;
                s[n] = (char*)calloc(p - q + 1, 1);
                if (!s[n])
                    goto err;
                strncpy(s[n++], q, p - q);
                q = p + 1;
            }
            if ( !*p ) break;
            p++;
        }
    }
    // Try to shrink s to the minimum size needed
    s_new = (char**)realloc(s, n * sizeof(char*));
    if (!s_new)
        goto err;

    s = s_new;
    assert(n < INT_MAX); // hts_resize() should ensure this
    *_n = n;
    return s;

 err:
    for (m = 0; m < n; m++)
        free(s[m]);
    free(s);
    return NULL;
}

char **hts_readlines(const char *fn, int *_n)
{
    unsigned int m = 0, n = 0;
    char **s = 0, **s_new;
    BGZF *fp = bgzf_open(fn, "r");
    if ( fp ) { // read from file
        kstring_t str;
        int ret;
        str.s = 0; str.l = str.m = 0;
        while ((ret = bgzf_getline(fp, '\n', &str)) >= 0) {
            if (str.l == 0) continue;
            if (hts_resize(char *, n + 1, &m, &s, 0) < 0)
                goto err;
            s[n] = strdup(str.s);
            if (!s[n])
                goto err;
            n++;
        }
        if (ret < -1) // Read error
            goto err;
        bgzf_close(fp);
        free(str.s);
    } else if (*fn == ':') { // read from string
        const char *q, *p;
        for (q = p = fn + 1;; ++p)
            if (*p == ',' || *p == 0) {
                if (hts_resize(char *, n + 1, &m, &s, 0) < 0)
                    goto err;
                s[n] = (char*)calloc(p - q + 1, 1);
                if (!s[n])
                    goto err;
                strncpy(s[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
    } else return 0;
    // Try to shrink s to the minimum size needed
    s_new = (char**)realloc(s, n * sizeof(char*));
    if (!s_new)
        goto err;

    s = s_new;
    assert(n < INT_MAX); // hts_resize() should ensure this
    *_n = n;
    return s;

 err:
    for (m = 0; m < n; m++)
        free(s[m]);
    free(s);
    return NULL;
}

// DEPRECATED: To be removed in a future HTSlib release
int hts_file_type(const char *fname)
{
    int len = strlen(fname);
    if ( !strcasecmp(".vcf.gz",fname+len-7) ) return FT_VCF_GZ;
    if ( !strcasecmp(".vcf",fname+len-4) ) return FT_VCF;
    if ( !strcasecmp(".bcf",fname+len-4) ) return FT_BCF_GZ;
    if ( !strcmp("-",fname) ) return FT_STDIN;

    hFILE *f = hopen(fname, "r");
    if (f == NULL) return 0;

    htsFormat fmt;
    if (hts_detect_format2(f, fname, &fmt) < 0) { hclose_abruptly(f); return 0; }
    if (hclose(f) < 0) return 0;

    switch (fmt.format) {
    case vcf: return (fmt.compression == no_compression)? FT_VCF : FT_VCF_GZ;
    case bcf: return (fmt.compression == no_compression)? FT_BCF : FT_BCF_GZ;
    default:  return 0;
    }
}

int hts_check_EOF(htsFile *fp)
{
    if (fp->format.compression == bgzf)
        return bgzf_check_EOF(hts_get_bgzfp(fp));
    else if (fp->format.format == cram)
        return cram_check_EOF(fp->fp.cram);
    else
        return 3;
}


/****************
 *** Indexing ***
 ****************/

#define HTS_MIN_MARKER_DIST 0x10000

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

#define pair64_lt(a,b) ((a).u < (b).u)
#define pair64max_lt(a,b) ((a).u < (b).u || \
                           ((a).u == (b).u && (a).max < (b).max))

KSORT_INIT_STATIC(_off, hts_pair64_t, pair64_lt)
KSORT_INIT_STATIC(_off_max, hts_pair64_max_t, pair64max_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    hts_pos_t n, m;
    uint64_t *offset;
} lidx_t;

struct hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    int tbi_n, last_tbi_tid;
    struct {
        uint32_t last_bin, save_bin;
        hts_pos_t last_coor;
        int last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

static char * idx_format_name(int fmt) {
    switch (fmt) {
        case HTS_FMT_CSI: return "csi";
        case HTS_FMT_BAI: return "bai";
        case HTS_FMT_TBI: return "tbi";
        case HTS_FMT_CRAI: return "crai";
        default: return "unknown";
    }
}

#ifdef DEBUG_INDEX
static void idx_dump(const hts_idx_t *idx) {
    int i;
    int64_t j;

    if (!idx) fprintf(stderr, "Null index\n");

    fprintf(stderr, "format='%s', min_shift=%d, n_lvls=%d, n_bins=%d, l_meta=%u ",
            idx_format_name(idx->fmt), idx->min_shift, idx->n_lvls, idx->n_bins, idx->l_meta);
    fprintf(stderr, "n=%d, m=%d, n_no_coor=%"PRIu64"\n", idx->n, idx->m, idx->n_no_coor);
    for (i = 0; i < idx->n; i++) {
        bidx_t *bidx = idx->bidx[i];
        lidx_t *lidx = &idx->lidx[i];
        if (bidx) {
            fprintf(stderr, "======== BIN Index - tid=%d, n_buckets=%d, size=%d\n", i, bidx->n_buckets, bidx->size);
            int b;
            for (b = 0; b < META_BIN(idx); b++) {
                khint_t k;
                if ((k = kh_get(bin, bidx, b)) != kh_end(bidx)) {
                    bins_t *entries = &kh_value(bidx, k);
                    int l = hts_bin_level(b);
                    int64_t bin_width = 1LL << ((idx->n_lvls - l) * 3 + idx->min_shift);
                    fprintf(stderr, "\tbin=%d, level=%d, parent=%d, n_chunks=%d, loff=%"PRIu64", interval=[%"PRId64" - %"PRId64"]\n",
                        b, l, hts_bin_parent(b), entries->n, entries->loff, (b-hts_bin_first(l))*bin_width+1, (b+1-hts_bin_first(l))*bin_width);
                    for (j = 0; j < entries->n; j++)
                        fprintf(stderr, "\t\tchunk=%"PRId64", u=%"PRIu64", v=%"PRIu64"\n", j, entries->list[j].u, entries->list[j].v);
                }
            }
        }
        if (lidx) {
            fprintf(stderr, "======== LINEAR Index - tid=%d, n_values=%"PRId64"\n", i, lidx->n);
            for (j = 0; j < lidx->n; j++) {
                fprintf(stderr, "\t\tentry=%"PRId64", offset=%"PRIu64", interval=[%"PRId64" - %"PRId64"]\n",
                    j, lidx->offset[j], j*(1<<idx->min_shift)+1, (j+1)*(1<<idx->min_shift));
            }
        }
    }
}
#endif

static inline int insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
    khint_t k;
    bins_t *l;
    int absent;
    k = kh_put(bin, b, bin, &absent);
    if (absent < 0) return -1; // Out of memory
    l = &kh_value(b, k);
    if (absent) {
        l->m = 1; l->n = 0;
        l->list = (hts_pair64_t*)calloc(l->m, sizeof(hts_pair64_t));
        if (!l->list) {
            kh_del(bin, b, k);
            return -1;
        }
    } else if (l->n == l->m) {
        uint32_t new_m = l->m ? l->m << 1 : 1;
        hts_pair64_t *new_list = realloc(l->list, new_m * sizeof(hts_pair64_t));
        if (!new_list) return -1;
        l->list = new_list;
        l->m = new_m;
    }
    l->list[l->n].u = beg;
    l->list[l->n++].v = end;
    return 0;
}

static inline int insert_to_l(lidx_t *l, int64_t _beg, int64_t _end, uint64_t offset, int min_shift)
{
    int i;
    hts_pos_t beg, end;
    beg = _beg >> min_shift;
    end = (_end - 1) >> min_shift;
    if (l->m < end + 1) {
        size_t new_m = l->m * 2 > end + 1 ? l->m * 2 : end + 1;
        uint64_t *new_offset;

        new_offset = (uint64_t*)realloc(l->offset, new_m * sizeof(uint64_t));
        if (!new_offset) return -1;

        // fill unused memory with (uint64_t)-1
        memset(new_offset + l->m, 0xff, sizeof(uint64_t) * (new_m - l->m));
        l->m = new_m;
        l->offset = new_offset;
    }
    for (i = beg; i <= end; ++i) {
        if (l->offset[i] == (uint64_t)-1) l->offset[i] = offset;
    }
    if (l->n < end + 1) l->n = end + 1;
    return 0;
}

hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
{
    hts_idx_t *idx;
    idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
    if (idx == NULL) return NULL;
    idx->fmt = fmt;
    idx->min_shift = min_shift;
    idx->n_lvls = n_lvls;
    idx->n_bins = ((1<<(3 * n_lvls + 3)) - 1) / 7;
    idx->z.save_tid = idx->z.last_tid = -1;
    idx->z.save_bin = idx->z.last_bin = 0xffffffffu;
    idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
    idx->z.last_coor = 0xffffffffu;
    if (n) {
        idx->n = idx->m = n;
        idx->bidx = (bidx_t**)calloc(n, sizeof(bidx_t*));
        if (idx->bidx == NULL) { free(idx); return NULL; }
        idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
        if (idx->lidx == NULL) { free(idx->bidx); free(idx); return NULL; }
    }
    idx->tbi_n = -1;
    idx->last_tbi_tid = -1;
    return idx;
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
    bidx_t *bidx = idx->bidx[i];
    lidx_t *lidx = &idx->lidx[i];
    khint_t k;
    int l;
    // the last entry is always valid
    for (l=lidx->n-2; l >= 0; l--) {
        if (lidx->offset[l] == (uint64_t)-1)
            lidx->offset[l] = lidx->offset[l+1];
    }
    if (bidx == 0) return;
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) // set loff
        if (kh_exist(bidx, k))
        {
            if ( kh_key(bidx, k) < idx->n_bins )
            {
                int bot_bin = hts_bin_bot(kh_key(bidx, k), idx->n_lvls);
                // disable linear index if bot_bin out of bounds
                kh_val(bidx, k).loff = bot_bin < lidx->n ? lidx->offset[bot_bin] : 0;
            }
            else
                kh_val(bidx, k).loff = 0;
        }
    if (free_lidx) {
        free(lidx->offset);
        lidx->m = lidx->n = 0;
        lidx->offset = 0;
    }
}

static int compress_binning(hts_idx_t *idx, int i)
{
    bidx_t *bidx = idx->bidx[i];
    khint_t k;
    int l, m;
    if (bidx == 0) return 0;
    // merge a bin to its parent if the bin is too small
    for (l = idx->n_lvls; l > 0; --l) {
        unsigned start = hts_bin_first(l);
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
            bins_t *p, *q;
            if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins || kh_key(bidx, k) < start) continue;
            p = &kh_value(bidx, k);
            if (l < idx->n_lvls && p->n > 1) ks_introsort(_off, p->n, p->list);
            if ((p->list[p->n - 1].v>>16) - (p->list[0].u>>16) < HTS_MIN_MARKER_DIST) {
                khint_t kp;
                kp = kh_get(bin, bidx, hts_bin_parent(kh_key(bidx, k)));
                if (kp == kh_end(bidx)) continue;
                q = &kh_val(bidx, kp);
                if (q->n + p->n > q->m) {
                    uint32_t new_m = q->n + p->n;
                    hts_pair64_t *new_list;
                    kroundup32(new_m);
                    if (new_m > INT32_MAX) return -1; // Limited by index format
                    new_list = realloc(q->list, new_m * sizeof(*new_list));
                    if (!new_list) return -1;
                    q->m = new_m;
                    q->list = new_list;
                }
                memcpy(q->list + q->n, p->list, p->n * sizeof(hts_pair64_t));
                q->n += p->n;
                free(p->list);
                kh_del(bin, bidx, k);
            }
        }
    }
    k = kh_get(bin, bidx, 0);
    if (k != kh_end(bidx)) ks_introsort(_off, kh_val(bidx, k).n, kh_val(bidx, k).list);
    // merge adjacent chunks that start from the same BGZF block
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
        bins_t *p;
        if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins) continue;
        p = &kh_value(bidx, k);
        for (l = 1, m = 0; l < p->n; ++l) {
            if (p->list[m].v>>16 >= p->list[l].u>>16) {
                if (p->list[m].v < p->list[l].v) p->list[m].v = p->list[l].v;
            } else p->list[++m] = p->list[l];
        }
        p->n = m + 1;
    }
    return 0;
}

int hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
    int i, ret = 0;
    if (idx == NULL || idx->z.finished) return 0; // do not run this function on an empty index or multiple times
    if (idx->z.save_tid >= 0) {
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, final_offset);
        ret |= insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
    }
    for (i = 0; i < idx->n; ++i) {
        update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
        ret |= compress_binning(idx, i);
    }
    idx->z.finished = 1;

    return ret;
}

int hts_idx_check_range(hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end)
{
    int64_t maxpos = (int64_t) 1 << (idx->min_shift + idx->n_lvls * 3);
    if (tid < 0 || (beg <= maxpos && end <= maxpos))
        return 0;

    if (idx->fmt == HTS_FMT_CSI) {
        hts_log_error("Region %"PRIhts_pos"..%"PRIhts_pos" "
                      "cannot be stored in a csi index with these parameters. "
                      "Please use a larger min_shift or depth",
                      beg, end);
    } else {
        hts_log_error("Region %"PRIhts_pos"..%"PRIhts_pos
                      " cannot be stored in a %s index. Try using a csi index",
                      beg, end, idx_format_name(idx->fmt));
    }
    errno = ERANGE;
    return -1;
}

int hts_idx_push(hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end, uint64_t offset, int is_mapped)
{
    int bin;
    if (tid<0) beg = -1, end = 0;
    if (hts_idx_check_range(idx, tid, beg, end) < 0)
        return -1;
    if (tid >= idx->m) { // enlarge the index
        uint32_t new_m = idx->m * 2 > tid + 1 ? idx->m * 2 : tid + 1;
        bidx_t **new_bidx;
        lidx_t *new_lidx;

        new_bidx = (bidx_t**)realloc(idx->bidx, new_m * sizeof(bidx_t*));
        if (!new_bidx) return -1;
        idx->bidx = new_bidx;

        new_lidx = (lidx_t*) realloc(idx->lidx, new_m * sizeof(lidx_t));
        if (!new_lidx) return -1;
        idx->lidx = new_lidx;

        memset(&idx->bidx[idx->m], 0, (new_m - idx->m) * sizeof(bidx_t*));
        memset(&idx->lidx[idx->m], 0, (new_m - idx->m) * sizeof(lidx_t));
        idx->m = new_m;
    }
    if (idx->n < tid + 1) idx->n = tid + 1;
    if (idx->z.finished) return 0;
    if (idx->z.last_tid != tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
        if ( tid>=0 && idx->n_no_coor )
        {
            hts_log_error("NO_COOR reads not in a single block at the end %d %d", tid, idx->z.last_tid);
            return -1;
        }
        if (tid>=0 && idx->bidx[tid] != 0)
        {
            hts_log_error("Chromosome blocks not continuous");
            return -1;
        }
        idx->z.last_tid = tid;
        idx->z.last_bin = 0xffffffffu;
    } else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
        hts_log_error("Unsorted positions on sequence #%d: %"PRIhts_pos" followed by %"PRIhts_pos, tid+1, idx->z.last_coor+1, beg+1);
        return -1;
    }
    if (end < beg) {
        // Malformed ranges are errors. (Empty ranges (beg==end) are unusual but acceptable.)
        hts_log_error("Invalid record on sequence #%d: end %"PRId64" < begin %"PRId64, tid+1, end, beg+1);
        return -1;
    }
    if ( tid>=0 )
    {
        if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
        // shoehorn [-1,0) (VCF POS=0) into the leftmost bottom-level bin
        if (beg < 0)  beg = 0;
        if (end <= 0) end = 1;
        // idx->z.last_off points to the start of the current record
        if (insert_to_l(&idx->lidx[tid], beg, end,
                        idx->z.last_off, idx->min_shift) < 0) return -1;
    }
    else idx->n_no_coor++;
    bin = hts_reg2bin(beg, end, idx->min_shift, idx->n_lvls);
    if ((int)idx->z.last_bin != bin) { // then possibly write the binning index
        if (idx->z.save_bin != 0xffffffffu) { // save_bin==0xffffffffu only happens to the first record
            if (insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin,
                            idx->z.save_off, idx->z.last_off) < 0) return -1;
        }
        if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
            idx->z.off_end = idx->z.last_off;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.off_beg, idx->z.off_end) < 0) return -1;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.n_mapped, idx->z.n_unmapped) < 0) return -1;
            idx->z.n_mapped = idx->z.n_unmapped = 0;
            idx->z.off_beg = idx->z.off_end;
        }
        idx->z.save_off = idx->z.last_off;
        idx->z.save_bin = idx->z.last_bin = bin;
        idx->z.save_tid = tid;
    }
    if (is_mapped) ++idx->z.n_mapped;
    else ++idx->z.n_unmapped;
    idx->z.last_off = offset;
    idx->z.last_coor = beg;
    return 0;
}

// Needed for TBI only.  Ensure 'tid' with 'name' is in the index meta data.
// idx->meta needs to have been initialised first with an appropriate Tabix
// configuration via hts_idx_set_meta.
//
// NB number of references (first 4 bytes of tabix header) aren't in
// idx->meta, but held in idx->n instead.
int hts_idx_tbi_name(hts_idx_t *idx, int tid, const char *name) {
    // Horrid - we have to map incoming tid to a tbi alternative tid.
    // This is because TBI counts tids by "covered" refs while everything
    // else counts by Nth SQ/contig record in header.
    if (tid == idx->last_tbi_tid || tid < 0 || !name)
        return idx->tbi_n;

    uint32_t len = strlen(name)+1;
    uint8_t *tmp = (uint8_t *)realloc(idx->meta, idx->l_meta + len);
    if (!tmp)
        return -1;

    // Append name
    idx->meta = tmp;
    strcpy((char *)idx->meta + idx->l_meta, name);
    idx->l_meta += len;

    // Update seq length
    u32_to_le(le_to_u32(idx->meta+24)+len, idx->meta+24);

    idx->last_tbi_tid = tid;
    return ++idx->tbi_n;
}

// When doing samtools index we have a read_bam / hts_idx_push(bgzf_tell())
// loop.  idx->z.last_off is the previous bzgf_tell location, so we know
// the location the current bam record started at as well as where it ends.
//
// When building an index on the fly via a write_bam / hts_idx_push loop,
// this isn't quite identical as we may amend the virtual coord returned
// by bgzf_tell to the start of a new block if the next bam struct doesn't
// fit.  It's essentially the same thing, but for bit-identical indices
// we need to amend the idx->z.last_off when we know we're starting a new
// block.
void hts_idx_amend_last(hts_idx_t *idx, uint64_t offset)
{
    idx->z.last_off = offset;
}

void hts_idx_destroy(hts_idx_t *idx)
{
    khint_t k;
    int i;
    if (idx == 0) return;

    // For HTS_FMT_CRAI, idx actually points to a different type -- see sam.c
    if (idx->fmt == HTS_FMT_CRAI) {
        hts_cram_idx_t *cidx = (hts_cram_idx_t *) idx;
        cram_index_free(cidx->cram);
        free(cidx);
        return;
    }

    for (i = 0; i < idx->m; ++i) {
        bidx_t *bidx = idx->bidx[i];
        free(idx->lidx[i].offset);
        if (bidx == 0) continue;
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
            if (kh_exist(bidx, k))
                free(kh_value(bidx, k).list);
        kh_destroy(bin, bidx);
    }
    free(idx->bidx); free(idx->lidx); free(idx->meta);
    free(idx);
}

int hts_idx_fmt(hts_idx_t *idx) {
    return idx->fmt;
}

// The optimizer eliminates these ed_is_big() calls; still it would be good to
// TODO Determine endianness at configure- or compile-time

static inline ssize_t HTS_RESULT_USED idx_write_int32(BGZF *fp, int32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint32(BGZF *fp, uint32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint64(BGZF *fp, uint64_t x)
{
    if (ed_is_big()) x = ed_swap_8(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline void swap_bins(bins_t *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        ed_swap_8p(&p->list[i].u);
        ed_swap_8p(&p->list[i].v);
    }
}

static int idx_save_core(const hts_idx_t *idx, BGZF *fp, int fmt)
{
    int32_t i, j;

    #define check(ret) if ((ret) < 0) return -1

    // VCF TBI/CSI only writes IDs for non-empty bins (ie covered references)
    //
    // NOTE: CSI meta is undefined in spec, so this code has an assumption
    // that we're only using it for Tabix data.
    int nids = idx->n;
    if (idx->meta && idx->l_meta >= 4 && le_to_u32(idx->meta) == TBX_VCF) {
        for (i = nids = 0; i < idx->n; ++i) {
            if (idx->bidx[i])
                nids++;
        }
    }
    check(idx_write_int32(fp, nids));
    if (fmt == HTS_FMT_TBI && idx->l_meta)
        check(bgzf_write(fp, idx->meta, idx->l_meta));

    for (i = 0; i < idx->n; ++i) {
        khint_t k;
        bidx_t *bidx = idx->bidx[i];
        lidx_t *lidx = &idx->lidx[i];

        // write binning index
        if (nids == idx->n || bidx)
            check(idx_write_int32(fp, bidx? kh_size(bidx) : 0));
        if (bidx)
            for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
                if (kh_exist(bidx, k)) {
                    bins_t *p = &kh_value(bidx, k);
                    check(idx_write_uint32(fp, kh_key(bidx, k)));
                    if (fmt == HTS_FMT_CSI) check(idx_write_uint64(fp, p->loff));
                    //int j;for(j=0;j<p->n;++j)fprintf(stderr,"%d,%llx,%d,%llx:%llx\n",kh_key(bidx,k),kh_val(bidx, k).loff,j,p->list[j].u,p->list[j].v);
                    check(idx_write_int32(fp, p->n));
                    for (j = 0; j < p->n; ++j) {
                        //fprintf(stderr, "\t%ld\t%ld\n", p->list[j].u, p->list[j].v);
                        check(idx_write_uint64(fp, p->list[j].u));
                        check(idx_write_uint64(fp, p->list[j].v));
                    }
                }

        // write linear index
        if (fmt != HTS_FMT_CSI) {
            check(idx_write_int32(fp, lidx->n));
            for (j = 0; j < lidx->n; ++j)
                check(idx_write_uint64(fp, lidx->offset[j]));
        }
    }

    check(idx_write_uint64(fp, idx->n_no_coor));
#ifdef DEBUG_INDEX
    idx_dump(idx);
#endif

    return 0;
    #undef check
}

int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
{
    int ret, save;
    if (idx == NULL || fn == NULL) { errno = EINVAL; return -1; }
    char *fnidx = (char*)calloc(1, strlen(fn) + 5);
    if (fnidx == NULL) return -1;

    strcpy(fnidx, fn);
    switch (fmt) {
    case HTS_FMT_BAI: strcat(fnidx, ".bai"); break;
    case HTS_FMT_CSI: strcat(fnidx, ".csi"); break;
    case HTS_FMT_TBI: strcat(fnidx, ".tbi"); break;
    default: abort();
    }

    ret = hts_idx_save_as(idx, fn, fnidx, fmt);
    save = errno;
    free(fnidx);
    errno = save;
    return ret;
}

int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)
{
    BGZF *fp;

    #define check(ret) if ((ret) < 0) goto fail

    if (fnidx == NULL) return hts_idx_save(idx, fn, fmt);

    fp = bgzf_open(fnidx, (fmt == HTS_FMT_BAI)? "wu" : "w");
    if (fp == NULL) return -1;

    if (fmt == HTS_FMT_CSI) {
        check(bgzf_write(fp, "CSI\1", 4));
        check(idx_write_int32(fp, idx->min_shift));
        check(idx_write_int32(fp, idx->n_lvls));
        check(idx_write_uint32(fp, idx->l_meta));
        if (idx->l_meta) check(bgzf_write(fp, idx->meta, idx->l_meta));
    } else if (fmt == HTS_FMT_TBI) {
        check(bgzf_write(fp, "TBI\1", 4));
    } else if (fmt == HTS_FMT_BAI) {
        check(bgzf_write(fp, "BAI\1", 4));
    } else abort();

    check(idx_save_core(idx, fp, fmt));

    return bgzf_close(fp);
    #undef check

fail:
    bgzf_close(fp);
    return -1;
}

static int idx_read_core(hts_idx_t *idx, BGZF *fp, int fmt)
{
    int32_t i, n, is_be;
    is_be = ed_is_big();
    if (idx == NULL) return -4;
    for (i = 0; i < idx->n; ++i) {
        bidx_t *h;
        lidx_t *l = &idx->lidx[i];
        uint32_t key;
        int j, absent;
        bins_t *p;
        h = idx->bidx[i] = kh_init(bin);
        if (bgzf_read(fp, &n, 4) != 4) return -1;
        if (is_be) ed_swap_4p(&n);
        if (n < 0) return -3;
        for (j = 0; j < n; ++j) {
            khint_t k;
            if (bgzf_read(fp, &key, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&key);
            k = kh_put(bin, h, key, &absent);
            if (absent <  0) return -2; // No memory
            if (absent == 0) return -3; // Duplicate bin number
            p = &kh_val(h, k);
            if (fmt == HTS_FMT_CSI) {
                if (bgzf_read(fp, &p->loff, 8) != 8) return -1;
                if (is_be) ed_swap_8p(&p->loff);
            } else p->loff = 0;
            if (bgzf_read(fp, &p->n, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&p->n);
            if (p->n < 0) return -3;
            if ((size_t) p->n > SIZE_MAX / sizeof(hts_pair64_t)) return -2;
            p->m = p->n;
            p->list = (hts_pair64_t*)malloc(p->m * sizeof(hts_pair64_t));
            if (p->list == NULL) return -2;
            if (bgzf_read(fp, p->list, ((size_t) p->n)<<4) != ((size_t) p->n)<<4) return -1;
            if (is_be) swap_bins(p);
        }
        if (fmt != HTS_FMT_CSI) { // load linear index
            int j, k;
            uint32_t x;
            if (bgzf_read(fp, &x, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&x);
            l->n = x;
            if (l->n < 0) return -3;
            if ((size_t) l->n > SIZE_MAX / sizeof(uint64_t)) return -2;
            l->m = l->n;
            l->offset = (uint64_t*)malloc(l->n * sizeof(uint64_t));
            if (l->offset == NULL) return -2;
            if (bgzf_read(fp, l->offset, l->n << 3) != l->n << 3) return -1;
            if (is_be) for (j = 0; j < l->n; ++j) ed_swap_8p(&l->offset[j]);
            for (k = j = 0; j < l->n && l->offset[j] == 0; k = ++j); // stop at the first non-zero entry
            for (j = l->n-1; j > k; j--) // fill missing values; may happen given older samtools and tabix
                if (l->offset[j-1] == 0) l->offset[j-1] = l->offset[j];
            update_loff(idx, i, 0);
        }
    }
    if (bgzf_read(fp, &idx->n_no_coor, 8) != 8) idx->n_no_coor = 0;
    if (is_be) ed_swap_8p(&idx->n_no_coor);
#ifdef DEBUG_INDEX
    idx_dump(idx);
#endif

    return 0;
}

static hts_idx_t *idx_read(const char *fn)
{
    uint8_t magic[4];
    int i, is_be;
    hts_idx_t *idx = NULL;
    uint8_t *meta = NULL;
    BGZF *fp = bgzf_open(fn, "r");
    if (fp == NULL) return NULL;
    is_be = ed_is_big();
    if (bgzf_read(fp, magic, 4) != 4) goto fail;

    if (memcmp(magic, "CSI\1", 4) == 0) {
        uint32_t x[3], n;
        if (bgzf_read(fp, x, 12) != 12) goto fail;
        if (is_be) for (i = 0; i < 3; ++i) ed_swap_4p(&x[i]);
        if (x[2]) {
            if (SIZE_MAX - x[2] < 1) goto fail; // Prevent possible overflow
            if ((meta = (uint8_t*)malloc((size_t) x[2] + 1)) == NULL) goto fail;
            if (bgzf_read(fp, meta, x[2]) != x[2]) goto fail;
            // Prevent possible strlen past the end in tbx_index_load2
            meta[x[2]] = '\0';
        }
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        if (n > INT32_MAX) goto fail;
        if ((idx = hts_idx_init(n, HTS_FMT_CSI, 0, x[0], x[1])) == NULL) goto fail;
        idx->l_meta = x[2];
        idx->meta = meta;
        meta = NULL;
        if (idx_read_core(idx, fp, HTS_FMT_CSI) < 0) goto fail;
    }
    else if (memcmp(magic, "TBI\1", 4) == 0) {
        uint8_t x[8 * 4];
        uint32_t n;
        // Read file header
        if (bgzf_read(fp, x, sizeof(x)) != sizeof(x)) goto fail;
        n = le_to_u32(&x[0]); // location of n_ref
        if (n > INT32_MAX) goto fail;
        if ((idx = hts_idx_init(n, HTS_FMT_TBI, 0, 14, 5)) == NULL) goto fail;
        n = le_to_u32(&x[7*4]); // location of l_nm
        if (n > UINT32_MAX - 29) goto fail; // Prevent possible overflow
        idx->l_meta = 28 + n;
        if ((idx->meta = (uint8_t*)malloc(idx->l_meta + 1)) == NULL) goto fail;
        // copy format, col_seq, col_beg, col_end, meta, skip, l_nm
        // N.B. left in little-endian byte order.
        memcpy(idx->meta, &x[1*4], 28);
        // Read in sequence names.
        if (bgzf_read(fp, idx->meta + 28, n) != n) goto fail;
        // Prevent possible strlen past the end in tbx_index_load2
        idx->meta[idx->l_meta] = '\0';
        if (idx_read_core(idx, fp, HTS_FMT_TBI) < 0) goto fail;
    }
    else if (memcmp(magic, "BAI\1", 4) == 0) {
        uint32_t n;
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        if (n > INT32_MAX) goto fail;
        if ((idx = hts_idx_init(n, HTS_FMT_BAI, 0, 14, 5)) == NULL) goto fail;
        if (idx_read_core(idx, fp, HTS_FMT_BAI) < 0) goto fail;
    }
    else { errno = EINVAL; goto fail; }

    bgzf_close(fp);
    return idx;

fail:
    bgzf_close(fp);
    hts_idx_destroy(idx);
    free(meta);
    return NULL;
}

int hts_idx_set_meta(hts_idx_t *idx, uint32_t l_meta, uint8_t *meta,
                      int is_copy)
{
    uint8_t *new_meta = meta;
    if (is_copy) {
        size_t l = l_meta;
        if (l > SIZE_MAX - 1) {
            errno = ENOMEM;
            return -1;
        }
        new_meta = malloc(l + 1);
        if (!new_meta) return -1;
        memcpy(new_meta, meta, l);
        // Prevent possible strlen past the end in tbx_index_load2
        new_meta[l] = '\0';
    }
    if (idx->meta) free(idx->meta);
    idx->l_meta = l_meta;
    idx->meta = new_meta;
    return 0;
}

uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta)
{
    *l_meta = idx->l_meta;
    return idx->meta;
}

const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)
{
    if ( !idx || !idx->n )
    {
        *n = 0;
        return NULL;
    }

    int tid = 0, i;
    const char **names = (const char**) calloc(idx->n,sizeof(const char*));
    for (i=0; i<idx->n; i++)
    {
        bidx_t *bidx = idx->bidx[i];
        if ( !bidx ) continue;
        names[tid++] = getid(hdr,i);
    }
    *n = tid;
    return names;
}

int hts_idx_nseq(const hts_idx_t *idx) {
    if (!idx) return -1;
    return idx->n;
}

int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)
{
    if (!idx) return -1;
    if ( idx->fmt == HTS_FMT_CRAI ) {
        *mapped = 0; *unmapped = 0;
        return -1;
    }

    bidx_t *h = idx->bidx[tid];
    if (!h) return -1;
    khint_t k = kh_get(bin, h, META_BIN(idx));
    if (k != kh_end(h)) {
        *mapped = kh_val(h, k).list[1].u;
        *unmapped = kh_val(h, k).list[1].v;
        return 0;
    } else {
        *mapped = 0; *unmapped = 0;
        return -1;
    }
}

uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx)
{
    if (idx->fmt == HTS_FMT_CRAI) return 0;
    return idx->n_no_coor;
}

/****************
 *** Iterator ***
 ****************/

// Note: even with 32-bit hts_pos_t, end needs to be 64-bit here due to 1LL<<s.
static inline int reg2bins(int64_t beg, int64_t end, hts_itr_t *itr, int min_shift, int n_lvls)
{
    int l, t, s = min_shift + (n_lvls<<1) + n_lvls;
    if (beg >= end) return 0;
    if (end >= 1LL<<s) end = 1LL<<s;
    for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        hts_pos_t b, e;
        int n, i;
        b = t + (beg>>s); e = t + (end>>s); n = e - b + 1;
        if (itr->bins.n + n > itr->bins.m) {
            itr->bins.m = itr->bins.n + n;
            kroundup32(itr->bins.m);
            itr->bins.a = (int*)realloc(itr->bins.a, sizeof(int) * itr->bins.m);
        }
        for (i = b; i <= e; ++i) itr->bins.a[itr->bins.n++] = i;
    }
    return itr->bins.n;
}

static inline int reg2intervals(hts_itr_t *iter, const hts_idx_t *idx, int tid, int64_t beg, int64_t end, uint32_t interval, uint64_t min_off, uint64_t max_off, int min_shift, int n_lvls)
{
    int l, t, s;
    int i, j;
    hts_pos_t b, e;
    hts_pair64_max_t *off;
    bidx_t *bidx;
    khint_t k;
    int start_n_off = iter->n_off;

    if (!iter || !idx || (bidx = idx->bidx[tid]) == NULL || beg >= end)
        return -1;

    s = min_shift + (n_lvls<<1) + n_lvls;
    if (end >= 1LL<<s)
        end = 1LL<<s;

    for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        b = t + (beg>>s); e = t + (end>>s);

        for (i = b; i <= e; ++i) {
            if ((k = kh_get(bin, bidx, i)) != kh_end(bidx)) {
                bins_t *p = &kh_value(bidx, k);

                if (p->n) {
                    off = realloc(iter->off, (iter->n_off + p->n) * sizeof(*off));
                    if (!off)
                        return -2;

                    iter->off = off;
                    for (j = 0; j < p->n; ++j) {
                        if (p->list[j].v > min_off && p->list[j].u < max_off) {
                            iter->off[iter->n_off].u = min_off > p->list[j].u
                                ? min_off : p->list[j].u;
                            iter->off[iter->n_off].v = max_off < p->list[j].v
                                ? max_off : p->list[j].v;
                            // hts_pair64_max_t::max is now used to link
                            // file offsets to region list entries.
                            // The iterator can use this to decide if it
                            // can skip some file regions.
                            iter->off[iter->n_off].max = ((uint64_t) tid << 32) | interval;
                            iter->n_off++;
                        }
                    }
                }
            }
        }
    }

    if (iter->n_off - start_n_off > 1) {
        ks_introsort(_off_max, iter->n_off - start_n_off, iter->off + start_n_off);
        for (i = start_n_off, j = start_n_off + 1; j < iter->n_off; j++) {
            if (iter->off[i].v >= iter->off[j].u) {
                if (iter->off[i].v < iter->off[j].v)
                    iter->off[i].v = iter->off[j].v;
            } else {
                i++;
                if (i < j)
                    iter->off[i] = iter->off[j];
            }
        }
        iter->n_off = i + 1;
    }

    return iter->n_off;
}

static int compare_regions(const void *r1, const void *r2) {
    hts_reglist_t *reg1 = (hts_reglist_t *)r1;
    hts_reglist_t *reg2 = (hts_reglist_t *)r2;

    if (reg1->tid < 0 && reg2->tid >= 0)
        return 1;
    else if (reg1->tid >= 0 && reg2->tid < 0)
        return -1;
    else
        return reg1->tid - reg2->tid;
}

uint64_t hts_itr_off(const hts_idx_t* idx, int tid) {

    int i;
    bidx_t* bidx;
    uint64_t off0 = (uint64_t) -1;
    khint_t k;
    switch (tid) {
    case HTS_IDX_START:
        // Find the smallest offset, note that sequence ids may not be ordered sequentially
        for (i = 0; i < idx->n; i++) {
            bidx = idx->bidx[i];
            k = kh_get(bin, bidx, META_BIN(idx));
            if (k == kh_end(bidx))
                continue;

            if (off0 > kh_val(bidx, k).list[0].u)
                off0 = kh_val(bidx, k).list[0].u;
        }
        if (off0 == (uint64_t) -1 && idx->n_no_coor)
            off0 = 0;
        // only no-coor reads in this bam
        break;
    case HTS_IDX_NOCOOR:
        /* No-coor reads sort after all of the mapped reads.  The position
           is not stored in the index itself, so need to find the end
           offset for the last mapped read.  A loop is needed here in
           case references at the end of the file have no mapped reads,
           or sequence ids are not ordered sequentially.
           See issue samtools#568 and commits b2aab8, 60c22d and cc207d. */
        for (i = 0; i < idx->n; i++) {
            bidx = idx->bidx[i];
            k = kh_get(bin, bidx, META_BIN(idx));
            if (k != kh_end(bidx)) {
                if (off0 == (uint64_t) -1 || off0 < kh_val(bidx, k).list[0].v) {
                    off0 = kh_val(bidx, k).list[0].v;
                }
            }
        }
        if (off0 == (uint64_t) -1 && idx->n_no_coor)
            off0 = 0;
        // only no-coor reads in this bam
        break;
    case HTS_IDX_REST:
        off0 = 0;
        break;
    case HTS_IDX_NONE:
        off0 = 0;
        break;
    }

    return off0;
}

hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end, hts_readrec_func *readrec)
{
    int i, n_off, l, bin;
    hts_pair64_max_t *off;
    khint_t k;
    bidx_t *bidx;
    uint64_t min_off, max_off;
    hts_itr_t *iter;
    uint32_t unmapped = 0, rel_off;

    // It's possible to call this function with NULL idx iff
    // tid is one of the special values HTS_IDX_REST or HTS_IDX_NONE
    if (!idx && !(tid == HTS_IDX_REST || tid == HTS_IDX_NONE)) {
        errno = EINVAL;
        return NULL;
    }

    iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
    if (iter) {
        if (tid < 0) {
            uint64_t off = hts_itr_off(idx, tid);
            if (off != (uint64_t) -1) {
                iter->read_rest = 1;
                iter->curr_off = off;
                iter->readrec = readrec;
                if (tid == HTS_IDX_NONE)
                    iter->finished = 1;
            } else {
                free(iter);
                iter = NULL;
            }
        } else if (tid >= idx->n || (bidx = idx->bidx[tid]) == NULL) {
            iter->finished = 1;
        } else {
            if (beg < 0) beg = 0;
            if (end < beg) {
              free(iter);
              return NULL;
            }

            k = kh_get(bin, bidx, META_BIN(idx));
            if (k != kh_end(bidx))
                unmapped = kh_val(bidx, k).list[1].v;
            else
                unmapped = 1;

            iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;
            iter->readrec = readrec;

            if ( !kh_size(bidx) ) { iter->finished = 1; return iter; }

            rel_off = beg>>idx->min_shift;
            // compute min_off
            bin = hts_bin_first(idx->n_lvls) + rel_off;
            do {
                int first;
                k = kh_get(bin, bidx, bin);
                if (k != kh_end(bidx)) break;
                first = (hts_bin_parent(bin)<<3) + 1;
                if (bin > first) --bin;
                else bin = hts_bin_parent(bin);
            } while (bin);
            if (bin == 0) k = kh_get(bin, bidx, bin);
            min_off = k != kh_end(bidx)? kh_val(bidx, k).loff : 0;
            // min_off can be calculated more accurately if the
            // linear index is available
            if (idx->lidx[tid].offset
                && rel_off < idx->lidx[tid].n) {
                if (min_off < idx->lidx[tid].offset[rel_off])
                    min_off = idx->lidx[tid].offset[rel_off];
                if (unmapped) {
                    // unmapped reads are not covered by the linear index,
                    // so search backwards for a smaller offset
                    int tmp_off;
                    for (tmp_off = rel_off-1; tmp_off >= 0; tmp_off--) {
                        if (idx->lidx[tid].offset[tmp_off] < min_off) {
                            min_off = idx->lidx[tid].offset[tmp_off];
                            break;
                        }
                    }
                    // if the search went too far back or no satisfactory entry
                    // was found, revert to the bin index loff value
                    if (k != kh_end(bidx) && (min_off < kh_val(bidx, k).loff || tmp_off < 0))
                        min_off = kh_val(bidx, k).loff;
                }
            } else if (unmapped) { //CSI index
                if (k != kh_end(bidx))
                    min_off = kh_val(bidx, k).loff;
            }

            // compute max_off: a virtual offset from a bin to the right of end
            bin = hts_bin_first(idx->n_lvls) + ((end-1) >> idx->min_shift) + 1;
            if (bin >= idx->n_bins) bin = 0;
            while (1) {
                // search for an extant bin by moving right, but moving up to the
                // parent whenever we get to a first child (which also covers falling
                // off the RHS, which wraps around and immediately goes up to bin 0)
                while (bin % 8 == 1) bin = hts_bin_parent(bin);
                if (bin == 0) { max_off = (uint64_t)-1; break; }
                k = kh_get(bin, bidx, bin);
                if (k != kh_end(bidx) && kh_val(bidx, k).n > 0) { max_off = kh_val(bidx, k).list[0].u; break; }
                bin++;
            }

            // retrieve bins
            reg2bins(beg, end, iter, idx->min_shift, idx->n_lvls);

            for (i = n_off = 0; i < iter->bins.n; ++i)
                if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx))
                    n_off += kh_value(bidx, k).n;
            if (n_off == 0) {
                // No overlapping bins means the iterator has already finished.
                iter->finished = 1;
                return iter;
            }
            off = calloc(n_off, sizeof(*off));
            for (i = n_off = 0; i < iter->bins.n; ++i) {
                if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx)) {
                    int j;
                    bins_t *p = &kh_value(bidx, k);
                    for (j = 0; j < p->n; ++j)
                        if (p->list[j].v > min_off && p->list[j].u < max_off) {
                            off[n_off].u = min_off > p->list[j].u
                                ? min_off : p->list[j].u;
                            off[n_off].v = max_off < p->list[j].v
                                ? max_off : p->list[j].v;
                            // hts_pair64_max_t::max is now used to link
                            // file offsets to region list entries.
                            // The iterator can use this to decide if it
                            // can skip some file regions.
                            off[n_off].max = ((uint64_t) tid << 32) | j;
                            n_off++;
                        }
                }
            }

            if (n_off == 0) {
                free(off);
                iter->finished = 1;
                return iter;
            }
            ks_introsort(_off_max, n_off, off);
            // resolve completely contained adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i)
                if (off[l].v < off[i].v) off[++l] = off[i];
            n_off = l + 1;
            // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
            for (i = 1; i < n_off; ++i)
                if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
            // merge adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i) {
                if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
                else off[++l] = off[i];
            }
            n_off = l + 1;
            iter->n_off = n_off; iter->off = off;
        }
    }

    return iter;
}

int hts_itr_multi_bam(const hts_idx_t *idx, hts_itr_t *iter)
{
    int i, j, bin;
    khint_t k;
    bidx_t *bidx;
    uint64_t min_off, max_off, t_off = (uint64_t)-1;
    int tid;
    hts_pos_t beg, end;
    hts_reglist_t *curr_reg;
    uint32_t unmapped = 0, rel_off;

    if (!idx || !iter || !iter->multi)
        return -1;

    iter->i = -1;
    for (i=0; i<iter->n_reg; i++) {

        curr_reg = &iter->reg_list[i];
        tid = curr_reg->tid;

        if (tid < 0) {
            t_off = hts_itr_off(idx, tid);
            if (t_off != (uint64_t)-1) {
                switch (tid) {
                case HTS_IDX_NONE:
                    iter->finished = 1;
                    // fall through
                case HTS_IDX_START:
                case HTS_IDX_REST:
                    iter->curr_off = t_off;
                    iter->n_reg = 0;
                    iter->reg_list = NULL;
                    iter->read_rest = 1;
                    return 0;
                case HTS_IDX_NOCOOR:
                    iter->nocoor = 1;
                    iter->nocoor_off = t_off;
                }
            }
        } else {
            if (tid >= idx->n || (bidx = idx->bidx[tid]) == NULL || !kh_size(bidx))
                continue;

            k = kh_get(bin, bidx, META_BIN(idx));
            if (k != kh_end(bidx))
                unmapped = kh_val(bidx, k).list[1].v;
            else
                unmapped = 1;

            for(j=0; j<curr_reg->count; j++) {
                hts_pair32_t *curr_intv = &curr_reg->intervals[j];
                if (curr_intv->end < curr_intv->beg)
                    continue;

                beg = curr_intv->beg;
                end = curr_intv->end;
                rel_off = beg>>idx->min_shift;

                /* Compute 'min_off' by searching the lowest level bin containing 'beg'.
                       If the computed bin is not in the index, try the next bin to the
                       left, belonging to the same parent. If it is the first sibling bin,
                       try the parent bin. */
                bin = hts_bin_first(idx->n_lvls) + rel_off;
                do {
                    int first;
                    k = kh_get(bin, bidx, bin);
                    if (k != kh_end(bidx)) break;
                    first = (hts_bin_parent(bin)<<3) + 1;
                    if (bin > first) --bin;
                    else bin = hts_bin_parent(bin);
                } while (bin);
                if (bin == 0)
                    k = kh_get(bin, bidx, bin);
                min_off = k != kh_end(bidx)? kh_val(bidx, k).loff : 0;
                // min_off can be calculated more accurately if the
                // linear index is available
                if (idx->lidx[tid].offset
                    && rel_off < idx->lidx[tid].n) {
                    if (min_off < idx->lidx[tid].offset[rel_off])
                        min_off = idx->lidx[tid].offset[rel_off];
                    if (unmapped) {
                        int tmp_off;
                        for (tmp_off = rel_off-1; tmp_off >= 0; tmp_off--) {
                            if (idx->lidx[tid].offset[tmp_off] < min_off) {
                                min_off = idx->lidx[tid].offset[tmp_off];
                                break;
                            }
                        }

                        if (k != kh_end(bidx) && (min_off < kh_val(bidx, k).loff || tmp_off < 0))
                            min_off = kh_val(bidx, k).loff;
                    }
                } else if (unmapped) { //CSI index
                    if (k != kh_end(bidx))
                        min_off = kh_val(bidx, k).loff;
                }

                // compute max_off: a virtual offset from a bin to the right of end
                bin = hts_bin_first(idx->n_lvls) + ((end-1) >> idx->min_shift) + 1;
                if (bin >= idx->n_bins) bin = 0;
                while (1) {
                    // search for an extant bin by moving right, but moving up to the
                    // parent whenever we get to a first child (which also covers falling
                    // off the RHS, which wraps around and immediately goes up to bin 0)
                    while (bin % 8 == 1) bin = hts_bin_parent(bin);
                    if (bin == 0) { max_off = (uint64_t)-1; break; }
                    k = kh_get(bin, bidx, bin);
                    if (k != kh_end(bidx) && kh_val(bidx, k).n > 0) {
                        max_off = kh_val(bidx, k).list[0].u;
                        break;
                    }
                    bin++;
                }

                //convert coordinates to file offsets
                if (reg2intervals(iter, idx, tid, beg, end, j,
                                  min_off, max_off,
                                  idx->min_shift, idx->n_lvls) < 0) {
                    return -1;
                }
            }
        }
    }

    if (iter->n_off > 1)
        ks_introsort(_off_max, iter->n_off, iter->off);

    if(!iter->n_off && !iter->nocoor)
        iter->finished = 1;

    return 0;
}

int hts_itr_multi_cram(const hts_idx_t *idx, hts_itr_t *iter)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    int tid, i, n_off = 0;
    uint32_t j;
    hts_pos_t beg, end;
    hts_reglist_t *curr_reg;
    hts_pair32_t *curr_intv;
    hts_pair64_max_t *off = NULL, *tmp;
    cram_index *e = NULL;

    if (!cidx || !iter || !iter->multi)
        return -1;

    iter->is_cram = 1;
    iter->read_rest = 0;
    iter->off = NULL;
    iter->n_off = 0;
    iter->curr_off = 0;
    iter->i = -1;

    for (i=0; i<iter->n_reg; i++) {

        curr_reg = &iter->reg_list[i];
        tid = curr_reg->tid;

        if (tid >= 0) {
            tmp = realloc(off, (n_off + curr_reg->count) * sizeof(*off));
            if (!tmp)
                goto err;
            off = tmp;

            for (j=0; j < curr_reg->count; j++) {
                curr_intv = &curr_reg->intervals[j];
                if (curr_intv->end < curr_intv->beg)
                    continue;

                beg = curr_intv->beg;
                end = curr_intv->end;

/* First, fetch the container overlapping 'beg' and assign its file offset to u, then
 * find the container overlapping 'end' and assign the relative end of the slice to v.
 * The cram_ptell function will adjust with the container offset, which is not stored
 * in the index.
 */
                e = cram_index_query(cidx->cram, tid, beg+1, NULL);
                if (e) {
                    off[n_off].u = e->offset;
                    // hts_pair64_max_t::max is now used to link
                    // file offsets to region list entries.
                    // The iterator can use this to decide if it
                    // can skip some file regions.
                    off[n_off].max = ((uint64_t) tid << 32) | j;

                    if (end >= HTS_POS_MAX) {
                       e = cram_index_last(cidx->cram, tid, NULL);
                    } else {
                       e = cram_index_query_last(cidx->cram, tid, end+1);
                    }

                    if (e) {
                        off[n_off++].v = e->next
                            ? e->next
                            : e->offset + e->slice + e->len;
                    } else {
                        hts_log_warning("Could not set offset end for region %d:%"PRIhts_pos"-%"PRIhts_pos". Skipping", tid, beg, end);
                    }
                } else {
                    hts_log_warning("No index entry for region %d:%"PRIhts_pos"-%"PRIhts_pos"", tid, beg, end);
                }
            }
        } else {
            switch (tid) {
                case HTS_IDX_NOCOOR:
                    e = cram_index_query(cidx->cram, tid, 1, NULL);
                    if (e) {
                        iter->nocoor = 1;
                        iter->nocoor_off = e->offset;
                    } else {
                        hts_log_warning("No index entry for NOCOOR region");
                    }
                    break;
                case HTS_IDX_START:
                    e = cram_index_query(cidx->cram, tid, 1, NULL);
                    if (e) {
                        iter->read_rest = 1;
                        tmp = realloc(off, sizeof(*off));
                        if (!tmp)
                            goto err;
                        off = tmp;
                        off[0].u = e->offset;
                        off[0].v = 0;
                        n_off=1;
                    } else {
                        hts_log_warning("No index entries");
                    }
                    break;
                case HTS_IDX_REST:
                    break;
                case HTS_IDX_NONE:
                    iter->finished = 1;
                    break;
                default:
                    hts_log_error("Query with tid=%d not implemented for CRAM files", tid);
            }
        }
    }

    if (n_off) {
        ks_introsort(_off_max, n_off, off);
        iter->n_off = n_off; iter->off = off;
    }

    if(!n_off && !iter->nocoor)
        iter->finished = 1;

    return 0;

 err:
    free(off);
    return -1;
}

void hts_itr_destroy(hts_itr_t *iter)
{
    if (iter) {
        if (iter->multi) {
            hts_reglist_free(iter->reg_list, iter->n_reg);
        } else {
            free(iter->bins.a);
        }

        if (iter->off)
            free(iter->off);
        free(iter);
    }
}

static inline long long push_digit(long long i, char c)
{
    // ensure subtraction occurs first, avoiding overflow for >= MAX-48 or so
    int digit = c - '0';
    return 10 * i + digit;
}

long long hts_parse_decimal(const char *str, char **strend, int flags)
{
    long long n = 0;
    int digits = 0, decimals = 0, e = 0, lost = 0;
    char sign = '+', esign = '+';
    const char *s, *str_orig = str;

    while (isspace_c(*str)) str++;
    s = str;

    if (*s == '+' || *s == '-') sign = *s++;
    while (*s)
        if (isdigit_c(*s)) digits++, n = push_digit(n, *s++);
        else if (*s == ',' && (flags & HTS_PARSE_THOUSANDS_SEP)) s++;
        else break;

    if (*s == '.') {
        s++;
        while (isdigit_c(*s)) decimals++, digits++, n = push_digit(n, *s++);
    }

    switch (*s) {
    case 'e': case 'E':
        s++;
        if (*s == '+' || *s == '-') esign = *s++;
        while (isdigit_c(*s)) e = push_digit(e, *s++);
        if (esign == '-') e = -e;
        break;

    case 'k': case 'K': e += 3; s++; break;
    case 'm': case 'M': e += 6; s++; break;
    case 'g': case 'G': e += 9; s++; break;
    }

    e -= decimals;
    while (e > 0) n *= 10, e--;
    while (e < 0) lost += n % 10, n /= 10, e++;

    if (lost > 0) {
        hts_log_warning("Discarding fractional part of %.*s", (int)(s - str), str);
    }

    if (strend) {
        // Set to the original input str pointer if not valid number syntax
        *strend = (digits > 0)? (char *)s : (char *)str_orig;
    } else if (digits == 0) {
        hts_log_warning("Invalid numeric value %.8s[truncated]", str);
    } else if (*s) {
        if ((flags & HTS_PARSE_THOUSANDS_SEP) || (!(flags & HTS_PARSE_THOUSANDS_SEP) && *s != ','))
            hts_log_warning("Ignoring unknown characters after %.*s[%s]", (int)(s - str), str, s);
    }

    return (sign == '+')? n : -n;
}

static void *hts_memrchr(const void *s, int c, size_t n) {
    size_t i;
    unsigned char *u = (unsigned char *)s;
    for (i = n; i > 0; i--) {
        if (u[i-1] == c)
            return u+i-1;
    }

    return NULL;
}

/*
 * A variant of hts_parse_reg which is reference-id aware.  It uses
 * the iterator name2id callbacks to validate the region tokenisation works.
 *
 * This is necessary due to GRCh38 HLA additions which have reference names
 * like "HLA-DRB1*12:17".
 *
 * All parameters are mandatory.
 *
 * To work around ambiguous parsing issues, eg both "chr1" and "chr1:100-200"
 * are reference names, we may quote using curly braces.
 * Thus "{chr1}:100-200" and "{chr1:100-200}" disambiguate the above example.
 *
 * Flags are used to control how parsing works, and can be one of the below.
 *
 * HTS_PARSE_LIST:
 *     If present, the region is assmed to be a comma separated list and
 *     position parsing will not contain commas (this implicitly
 *     clears HTS_PARSE_THOUSANDS_SEP in the call to hts_parse_decimal).
 *     On success the return pointer will be the start of the next region, ie
 *     the character after the comma.  (If *ret != '\0' then the caller can
 *     assume another region is present in the list.)
 *
 *     If not set then positions may contain commas.  In this case the return
 *     value should point to the end of the string, or NULL on failure.
 *
 * HTS_PARSE_ONE_COORD:
 *     If present, X:100 is treated as the single base pair region X:100-100.
 *     In this case X:-100 is shorthand for X:1-100 and X:100- is X:100-<end>.
 *     (This is the standard bcftools region convention.)
 *
 *     When not set X:100 is considered to be X:100-<end> where <end> is
 *     the end of chromosome X (set to HTS_POS_MAX here).  X:100- and X:-100
 *     are invalid.
 *     (This is the standard samtools region convention.)
 *
 * Note the supplied string expects 1 based inclusive coordinates, but the
 * returned coordinates start from 0 and are half open, so pos0 is valid
 * for use in e.g. "for (pos0 = beg; pos0 < end; pos0++) {...}"
 *
 * On success a pointer to the byte after the end of the entire region
 *            specifier is returned (plus any trailing comma), and tid,
 *            beg & end will be set.
 * On failure NULL is returned.
 */
const char *hts_parse_region(const char *s, int *tid, hts_pos_t *beg,
                             hts_pos_t *end, hts_name2id_f getid, void *hdr,
                             int flags)
{
    if (!s || !tid || !beg || !end || !getid)
        return NULL;

    size_t s_len = strlen(s);
    kstring_t ks = { 0, 0, NULL };

    const char *colon = NULL, *comma = NULL;
    int quoted = 0;

    if (flags & HTS_PARSE_LIST)
        flags &= ~HTS_PARSE_THOUSANDS_SEP;
    else
        flags |= HTS_PARSE_THOUSANDS_SEP;

    const char *s_end = s + s_len;

    // Braced quoting of references is permitted to resolve ambiguities.
    if (*s == '{') {
        const char *close = memchr(s, '}', s_len);
        if (!close) {
            hts_log_error("Mismatching braces in \"%s\"", s);
            *tid = -1;
            return NULL;
        }
        s++;
        s_len--;
        if (close[1] == ':')
            colon = close+1;
        quoted = 1; // number of trailing characters to trim

        // Truncate to this item only, if appropriate.
        if (flags & HTS_PARSE_LIST) {
            comma = strchr(close, ',');
            if (comma) {
                s_len = comma-s;
                s_end = comma+1;
            }
        }
    } else {
        // Truncate to this item only, if appropriate.
        if (flags & HTS_PARSE_LIST) {
            comma = strchr(s, ',');
            if (comma) {
                s_len = comma-s;
                s_end = comma+1;
            }
        }

        colon = hts_memrchr(s, ':', s_len);
    }

    // No colon is simplest case; just check and return.
    if (colon == NULL) {
        *beg = 0; *end = HTS_POS_MAX;
        kputsn(s, s_len-quoted, &ks); // convert to nul terminated string
        if (!ks.s) {
            *tid = -2;
            return NULL;
        }

        *tid = getid(hdr, ks.s);
        free(ks.s);

        return *tid >= 0 ? s_end : NULL;
    }

    // Has a colon, but check whole name first.
    if (!quoted) {
        *beg = 0; *end = HTS_POS_MAX;
        kputsn(s, s_len, &ks); // convert to nul terminated string
        if (!ks.s) {
            *tid = -2;
            return NULL;
        }
        if ((*tid = getid(hdr, ks.s)) >= 0) {
            // Entire name matches, but also check this isn't
            // ambiguous.  eg we have ref chr1 and ref chr1:100-200
            // both present.
            ks.l = 0;
            kputsn(s, colon-s, &ks); // convert to nul terminated string
            if (!ks.s) {
                *tid = -2;
                return NULL;
            }
            if (getid(hdr, ks.s) >= 0) {
                free(ks.s);
                *tid = -1;
                hts_log_error("Range is ambiguous. "
                              "Use {%s} or {%.*s}%s instead",
                              s, (int)(colon-s), s, colon);
                return NULL;
            }
            free(ks.s);

            return s_end;
        }
        if (*tid < -1) // Failed to parse header
            return NULL;
    }

    // Quoted, or unquoted and whole string isn't a name.
    // Check the pre-colon part is valid.
    ks.l = 0;
    kputsn(s, colon-s-quoted, &ks); // convert to nul terminated string
    if (!ks.s) {
        *tid = -2;
        return NULL;
    }
    *tid = getid(hdr, ks.s);
    free(ks.s);
    if (*tid < 0)
        return NULL;

    // Finally parse the post-colon coordinates
    char *hyphen;
    *beg = hts_parse_decimal(colon+1, &hyphen, flags) - 1;
    if (*beg < 0) {
        if (*beg != -1 && *hyphen == '-' && colon[1] != '\0') {
            // User specified zero, but we're 1-based.
            hts_log_error("Coordinates must be > 0");
            return NULL;
        }
        if (isdigit_c(*hyphen) || *hyphen == '\0' || *hyphen == ',') {
            // interpret chr:-100 as chr:1-100
            *end = *beg==-1 ? HTS_POS_MAX : -(*beg+1);
            *beg = 0;
            return s_end;
        } else if (*beg < -1) {
            hts_log_error("Unexpected string \"%s\" after region", hyphen);
            return NULL;
        }
    }

    if (*hyphen == '\0' || ((flags & HTS_PARSE_LIST) && *hyphen == ',')) {
        *end = flags & HTS_PARSE_ONE_COORD ? *beg+1 : HTS_POS_MAX;
    } else if (*hyphen == '-') {
        *end = hts_parse_decimal(hyphen+1, &hyphen, flags);
        if (*hyphen != '\0' && *hyphen != ',') {
            hts_log_error("Unexpected string \"%s\" after region", hyphen);
            return NULL;
        }
    } else {
        hts_log_error("Unexpected string \"%s\" after region", hyphen);
        return NULL;
    }

    if (*end == 0)
        *end = HTS_POS_MAX; // interpret chr:100- as chr:100-<end>

    if (*beg >= *end) return NULL;

    return s_end;
}

// Next release we should mark this as deprecated?
// Use hts_parse_region above instead.
const char *hts_parse_reg64(const char *s, hts_pos_t *beg, hts_pos_t *end)
{
    char *hyphen;
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = HTS_POS_MAX;
        return s + strlen(s);
    }

    *beg = hts_parse_decimal(colon+1, &hyphen, HTS_PARSE_THOUSANDS_SEP) - 1;
    if (*beg < 0) *beg = 0;

    if (*hyphen == '\0') *end = HTS_POS_MAX;
    else if (*hyphen == '-') *end = hts_parse_decimal(hyphen+1, NULL, HTS_PARSE_THOUSANDS_SEP);
    else return NULL;

    if (*beg >= *end) return NULL;
    return colon;
}

const char *hts_parse_reg(const char *s, int *beg, int *end)
{
    hts_pos_t beg64 = 0, end64 = 0;
    const char *colon = hts_parse_reg64(s, &beg64, &end64);
    if (beg64 > INT_MAX) {
        hts_log_error("Position %"PRId64" too large", beg64);
        return NULL;
    }
    if (end64 > INT_MAX) {
        if (end64 == HTS_POS_MAX) {
            end64 = INT_MAX;
        } else {
            hts_log_error("Position %"PRId64" too large", end64);
            return NULL;
        }
    }
    *beg = beg64;
    *end = end64;
    return colon;
}

hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)
{
    int tid;
    hts_pos_t beg, end;

    if (strcmp(reg, ".") == 0)
        return itr_query(idx, HTS_IDX_START, 0, 0, readrec);
    else if (strcmp(reg, "*") == 0)
        return itr_query(idx, HTS_IDX_NOCOOR, 0, 0, readrec);

    if (!hts_parse_region(reg, &tid, &beg, &end, getid, hdr, HTS_PARSE_THOUSANDS_SEP))
        return NULL;

    return itr_query(idx, tid, beg, end, readrec);
}

hts_itr_t *hts_itr_regions(const hts_idx_t *idx, hts_reglist_t *reglist, int count, hts_name2id_f getid, void *hdr, hts_itr_multi_query_func *itr_specific, hts_readrec_func *readrec, hts_seek_func *seek, hts_tell_func *tell) {

    int i;

    if (!reglist)
        return NULL;

    hts_itr_t *itr = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
    if (itr) {
        itr->n_reg = count;
        itr->readrec = readrec;
        itr->seek = seek;
        itr->tell = tell;
        itr->reg_list = reglist;
        itr->finished = 0;
        itr->nocoor = 0;
        itr->multi = 1;

        for (i = 0; i < itr->n_reg; i++) {
            if (itr->reg_list[i].reg) {
                if (!strcmp(itr->reg_list[i].reg, ".")) {
                    itr->reg_list[i].tid = HTS_IDX_START;
                    continue;
                }

                if (!strcmp(itr->reg_list[i].reg, "*")) {
                    itr->reg_list[i].tid = HTS_IDX_NOCOOR;
                    continue;
                }

                itr->reg_list[i].tid = getid(hdr, reglist[i].reg);
                if (itr->reg_list[i].tid < 0) {
                    if (itr->reg_list[i].tid < -1) {
                        hts_log_error("Failed to parse header");
                        hts_itr_destroy(itr);
                        return NULL;
                    } else {
                        hts_log_warning("Region '%s' specifies an unknown reference name. Continue anyway", reglist[i].reg);
                    }
                }
            }
        }

        qsort(itr->reg_list, itr->n_reg, sizeof(hts_reglist_t), compare_regions);
        if (itr_specific(idx, itr) != 0) {
            hts_log_error("Failed to create the multi-region iterator!");
            hts_itr_destroy(itr);
            itr = NULL;
        }
    }

    return itr;
}

int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
{
    int ret, tid;
    hts_pos_t beg, end;
    if (iter == NULL || iter->finished) return -1;
    if (iter->read_rest) {
        if (iter->curr_off) { // seek to the start
            if (bgzf_seek(fp, iter->curr_off, SEEK_SET) < 0) {
                hts_log_error("Failed to seek to offset %"PRIu64"%s%s",
                              iter->curr_off,
                              errno ? ": " : "", strerror(errno));
                return -2;
            }
            iter->curr_off = 0; // only seek once
        }
        ret = iter->readrec(fp, data, r, &tid, &beg, &end);
        if (ret < 0) iter->finished = 1;
        iter->curr_tid = tid;
        iter->curr_beg = beg;
        iter->curr_end = end;
        return ret;
    }
    // A NULL iter->off should always be accompanied by iter->finished.
    assert(iter->off != NULL);
    for (;;) {
        if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
            if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
            if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
                if (bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET) < 0) {
                    hts_log_error("Failed to seek to offset %"PRIu64"%s%s",
                                  iter->off[iter->i+1].u,
                                  errno ? ": " : "", strerror(errno));
                    return -2;
                }
                iter->curr_off = bgzf_tell(fp);
            }
            ++iter->i;
        }
        if ((ret = iter->readrec(fp, data, r, &tid, &beg, &end)) >= 0) {
            iter->curr_off = bgzf_tell(fp);
            if (tid != iter->tid || beg >= iter->end) { // no need to proceed
                ret = -1; break;
            } else if (end > iter->beg && iter->end > beg) {
                iter->curr_tid = tid;
                iter->curr_beg = beg;
                iter->curr_end = end;
                return ret;
            }
        } else break; // end of file or error
    }
    iter->finished = 1;
    return ret;
}

int hts_itr_multi_next(htsFile *fd, hts_itr_t *iter, void *r)
{
    void *fp;
    int ret, tid, i, cr, ci;
    hts_pos_t beg, end;
    hts_reglist_t *found_reg;

    if (iter == NULL || iter->finished) return -1;

    if (iter->is_cram) {
        fp = fd->fp.cram;
    } else {
        fp = fd->fp.bgzf;
    }

    if (iter->read_rest) {
        if (iter->curr_off) { // seek to the start
            if (iter->seek(fp, iter->curr_off, SEEK_SET) < 0) {
                hts_log_error("Seek at offset %" PRIu64 " failed.", iter->curr_off);
                return -1;
            }
            iter->curr_off = 0; // only seek once
        }

        ret = iter->readrec(fp, fd, r, &tid, &beg, &end);
        if (ret < 0)
            iter->finished = 1;

        iter->curr_tid = tid;
        iter->curr_beg = beg;
        iter->curr_end = end;

        return ret;
    }
    // A NULL iter->off should always be accompanied by iter->finished.
    assert(iter->off != NULL || iter->nocoor != 0);

    int next_range = 0;
    for (;;) {
        // Note that due to the way bam indexing works, iter->off may contain
        // file chunks that are not actually needed as they contain data
        // beyond the end of the requested region.  These are filtered out
        // by comparing the tid and index into hts_reglist_t::intervals
        // (packed for reasons of convenience into iter->off[iter->i].max)
        // associated with the file region with iter->curr_tid and
        // iter->curr_intv.

        if (next_range
            || iter->curr_off == 0
            || iter->i >= iter->n_off
            || iter->curr_off >= iter->off[iter->i].v
            || (iter->off[iter->i].max >> 32 == iter->curr_tid
                && (iter->off[iter->i].max & 0xffffffff) < iter->curr_intv)) {

            // Jump to the next chunk.  It may be necessary to skip more
            // than one as the iter->off list can include overlapping entries.
            do {
                iter->i++;
            } while (iter->i < iter->n_off
                     && (iter->curr_off >= iter->off[iter->i].v
                         || (iter->off[iter->i].max >> 32 == iter->curr_tid
                             && (iter->off[iter->i].max & 0xffffffff) < iter->curr_intv)));

            if (iter->is_cram && iter->i < iter->n_off) {
                // Ensure iter->curr_reg is correct.
                //
                // We need this for CRAM as we shortcut some of the later
                // logic by getting an end-of-range and continuing to the
                // next offset.
                //
                // We cannot do this for BAM (and fortunately do not need to
                // either) because in BAM world a query to genomic positions
                // GX and GY leading to a seek offsets PX and PY may have
                // GX > GY and PX < PY.  (This is due to the R-tree and falling
                // between intervals, bumping up to a higher bin.)
                // CRAM strictly follows PX >= PY if GX >= GY, so this logic
                // works.
                int want_tid = iter->off[iter->i].max >> 32;
                if (!(iter->curr_reg < iter->n_reg &&
                      iter->reg_list[iter->curr_reg].tid == want_tid)) {
                    int j;
                    for (j = 0; j < iter->n_reg; j++)
                        if (iter->reg_list[j].tid == want_tid)
                            break;
                    if (j == iter->n_reg)
                        return -1;
                    iter->curr_reg = j;
                    iter->curr_tid = iter->reg_list[iter->curr_reg].tid;
                };
                iter->curr_intv = iter->off[iter->i].max & 0xffffffff;
            }

            if (iter->i >= iter->n_off) { // no more chunks, except NOCOORs
                if (iter->nocoor) {
                    next_range = 0;
                    if (iter->seek(fp, iter->nocoor_off, SEEK_SET) < 0) {
                        hts_log_error("Seek at offset %" PRIu64 " failed.", iter->nocoor_off);
                        return -1;
                    }
                    if (iter->is_cram) {
                        cram_range r = { HTS_IDX_NOCOOR };
                        cram_set_option(fp, CRAM_OPT_RANGE_NOSEEK, &r);
                    }

                    // The first slice covering the unmapped reads might
                    // contain a few mapped reads, so scroll
                    // forward until finding the first unmapped read.
                    do {
                        ret = iter->readrec(fp, fd, r, &tid, &beg, &end);
                    } while (tid >= 0 && ret >=0);

                    if (ret < 0)
                        iter->finished = 1;
                    else
                        iter->read_rest = 1;

                    iter->curr_off = 0; // don't seek any more
                    iter->curr_tid = tid;
                    iter->curr_beg = beg;
                    iter->curr_end = end;

                    return ret;
                } else {
                    ret = -1; break;
                }
            } else if (iter->i < iter->n_off) {
                // New chunk may overlap the last one, so ensure we
                // only seek forwards.
                if (iter->curr_off < iter->off[iter->i].u || next_range) {
                    iter->curr_off = iter->off[iter->i].u;

                    // CRAM has the capability of setting an end location.
                    // This means multi-threaded decodes can stop once they
                    // reach that point, rather than pointlessly decoding
                    // more slices than we'll be using.
                    //
                    // We have to be careful here.  Whenever we set the cram
                    // range we need a corresponding seek in order to ensure
                    // we can safely decode at that offset.  We use next_range
                    // var to ensure this is always true; this is set on
                    // end-of-range condition. It's never modified for BAM.
                    if (iter->is_cram) {
                        // Next offset.[uv] tuple, but it's already been
                        // included in our cram range, so don't seek and don't
                        // reset range so we can efficiently multi-thread.
                        if (next_range || iter->curr_off >= iter->end) {
                            if (iter->seek(fp, iter->curr_off, SEEK_SET) < 0) {
                                hts_log_error("Seek at offset %" PRIu64
                                        " failed.", iter->curr_off);
                                return -1;
                            }

                            // Find the genomic range matching this interval.
                            int j;
                            hts_reglist_t *rl = &iter->reg_list[iter->curr_reg];
                            cram_range r = {
                                    rl->tid,
                                    rl->intervals[iter->curr_intv].beg,
                                    rl->intervals[iter->curr_intv].end
                            };

                            // Expand it up to cover neighbouring intervals.
                            // Note we can only have a single chromosome in a
                            // range, so if we detect our blocks span chromosomes
                            // or we have a multi-ref mode slice, we just use
                            // HTS_IDX_START refid instead.  This doesn't actually
                            // seek (due to CRAM_OPT_RANGE_NOSEEK) and is simply
                            // and indicator of decoding with no end limit.
                            //
                            // That isn't as efficient as it could be, but it's
                            // no poorer than before and it works.
                            int tid = r.refid;
                            int64_t end = r.end;
                            int64_t v = iter->off[iter->i].v;
                            j = iter->i+1;
                            while (j < iter->n_off) {
                                if (iter->off[j].u > v)
                                    break;

                                uint64_t max = iter->off[j].max;
                                if ((max>>32) != tid)
                                    tid = HTS_IDX_START; // => no range limit

                                if (end < rl->intervals[max & 0xffffffff].end)
                                    end = rl->intervals[max & 0xffffffff].end;
                                if (v < iter->off[j].v)
                                    v = iter->off[j].v;
                                j++;
                            }
                            r.refid = tid;
                            r.end = end;

                            // Remember maximum 'v' here so we don't do
                            // unnecessary subsequent seeks for the next
                            // regions.  We can't change curr_off, but
                            // beg/end are used only by single region iterator so
                            // we cache it there to avoid changing the struct.
                            iter->end = v;

                            cram_set_option(fp, CRAM_OPT_RANGE_NOSEEK, &r);
                            next_range = 0;
                        }
                    } else { // Not CRAM
                        if (iter->seek(fp, iter->curr_off, SEEK_SET) < 0) {
                            hts_log_error("Seek at offset %" PRIu64 " failed.",
                                          iter->curr_off);
                            return -1;
                        }
                    }
                }
            }
        }

        ret = iter->readrec(fp, fd, r, &tid, &beg, &end);
        if (ret < 0) {
            if (iter->is_cram && cram_eof(fp)) {
                // Skip to end of range
                //
                // We should never be adjusting curr_off manually unless
                // we also can guarantee we'll be doing a seek after to
                // a new location.  Otherwise we'll be reading wrong offset
                // for the next container.
                //
                // We ensure this by adjusting our CRAM_OPT_RANGE
                // accordingly above, but to double check we also
                // set the skipped_block flag to enforce a seek also.
                iter->curr_off = iter->off[iter->i].v;
                next_range = 1;

                // Next region
                if (++iter->curr_intv >= iter->reg_list[iter->curr_reg].count){
                    if (++iter->curr_reg >= iter->n_reg)
                        break;
                    iter->curr_intv = 0;
                    iter->curr_tid = iter->reg_list[iter->curr_reg].tid;
                }
                continue;
            } else {
                break;
            }
        }

        iter->curr_off = iter->tell(fp);

        if (tid != iter->curr_tid) {
            hts_reglist_t key;
            key.tid = tid;

            found_reg = (hts_reglist_t *)bsearch(&key, iter->reg_list,
                                                 iter->n_reg,
                                                 sizeof(hts_reglist_t),
                                                 compare_regions);
            if (!found_reg)
                continue;

            iter->curr_reg = (found_reg - iter->reg_list);
            iter->curr_tid = tid;
            iter->curr_intv = 0;
        }

        cr = iter->curr_reg;
        ci = iter->curr_intv;

        for (i = ci; i < iter->reg_list[cr].count; i++) {
            if (end > iter->reg_list[cr].intervals[i].beg &&
                iter->reg_list[cr].intervals[i].end > beg) {
                iter->curr_beg = beg;
                iter->curr_end = end;
                iter->curr_intv = i;

                return ret;
            }

            // Check if the read starts beyond intervals[i].end
            // If so, the interval is finished so move on to the next.
            if (beg > iter->reg_list[cr].intervals[i].end)
                iter->curr_intv = i + 1;

            // No need to keep searching if the read ends before intervals[i].beg
            if (end < iter->reg_list[cr].intervals[i].beg)
                break;
        }
    }
    iter->finished = 1;

    return ret;
}

/**********************
 *** Retrieve index ***
 **********************/
// Local_fn and local_len will return a sub-region of 'fn'.
// Eg http://elsewhere/dir/foo.bam.bai?a=b may return
// foo.bam.bai via local_fn and local_len.
//
// Returns -1 if index couldn't be opened.
//         -2 on other errors
static int idx_test_and_fetch(const char *fn, const char **local_fn, int *local_len, int download)
{
    hFILE *remote_hfp = NULL;
    hFILE *local_fp = NULL;
    int save_errno;
    htsFormat fmt;
    kstring_t s = KS_INITIALIZE;
    kstring_t tmps = KS_INITIALIZE;

    if (hisremote(fn)) {
        const int buf_size = 1 * 1024 * 1024;
        int l;
        const char *p, *e;
        // Ignore ?# params: eg any file.fmt?param=val, except for S3 URLs
        e = fn + ((strncmp(fn, "s3://", 5) && strncmp(fn, "s3+http://", 10) && strncmp(fn, "s3+https://", 11)) ? strcspn(fn, "?#") : strcspn(fn, "?"));
        // Find the previous slash from there.
        p = e;
        while (p > fn && *p != '/') p--;
        if (*p == '/') p++;

        // Attempt to open local file first
        kputsn(p, e-p, &s);
        if (access(s.s, R_OK) == 0)
        {
            free(s.s);
            *local_fn = p;
            *local_len = e-p;
            return 0;
        }

        // Attempt to open remote file. Stay quiet on failure, it is OK to fail when trying first .csi then .bai or .tbi index.
        if ((remote_hfp = hopen(fn, "r")) == 0) {
            hts_log_info("Failed to open index file '%s'", fn);
            free(s.s);
            return -1;
        }
        if (hts_detect_format2(remote_hfp, fn, &fmt)) {
            hts_log_error("Failed to detect format of index file '%s'", fn);
            goto fail;
        }
        if (fmt.category != index_file || (fmt.format != bai &&  fmt.format != csi && fmt.format != tbi
                && fmt.format != crai && fmt.format != fai_format)) {
            hts_log_error("Format of index file '%s' is not supported", fn);
            goto fail;
        }

        if (download) {
            if ((local_fp = hts_open_tmpfile(s.s, "wx", &tmps)) == NULL) {
                hts_log_error("Failed to create file %s in the working directory", p);
                goto fail;
            }
            hts_log_info("Downloading file %s to local directory", fn);
            uint8_t *buf = (uint8_t*)calloc(buf_size, 1);
            if (!buf) {
                hts_log_error("%s", strerror(errno));
                goto fail;
            }
            while ((l = hread(remote_hfp, buf, buf_size)) > 0) {
                if (hwrite(local_fp, buf, l) != l) {
                    hts_log_error("Failed to write data to %s : %s",
                            fn, strerror(errno));
                    free(buf);
                    goto fail;
                }
            }
            free(buf);
            if (l < 0) {
                hts_log_error("Error reading \"%s\"", fn);
                goto fail;
            }
            if (hclose(local_fp) < 0) {
                hts_log_error("Error closing %s : %s", fn, strerror(errno));
                local_fp = NULL;
                goto fail;
            }
            local_fp = NULL;
            if (rename(tmps.s, s.s) < 0) {
                hts_log_error("Error renaming %s : %s", tmps.s, strerror(errno));
                goto fail;
            }
            ks_clear(&tmps);

            *local_fn = p;
            *local_len = e-p;
        } else {
            *local_fn = fn;
            *local_len = e-fn;
        }

        if (hclose(remote_hfp) != 0) {
            hts_log_error("Failed to close remote file %s", fn);
        }

        free(tmps.s);
        free(s.s);
        return 0;
    } else {
        hFILE *local_hfp;
        if ((local_hfp = hopen(fn, "r")) == 0) return -1;
        hclose_abruptly(local_hfp);
        *local_fn = fn;
        *local_len = strlen(fn);
        return 0;
    }

 fail:
    save_errno = errno;
    if (remote_hfp) hclose_abruptly(remote_hfp);
    if (local_fp) hclose_abruptly(local_fp);
    if (tmps.l > 0) unlink(tmps.s);
    free(tmps.s);
    free(s.s);
    errno = save_errno;
    return -2;
}

/*
 * Check the existence of a local index file using part of the alignment file name.
 * The order is alignment.bam.csi, alignment.csi, alignment.bam.bai, alignment.bai
 * @param fn    - pointer to the file name
 * @param fnidx - pointer to the index file name placeholder
 * @return        1 for success, 0 for failure
 */
int hts_idx_check_local(const char *fn, int fmt, char **fnidx) {
    int i, l_fn, l_ext;
    const char *fn_tmp = NULL;
    char *fnidx_tmp;
    char *csi_ext = ".csi";
    char *bai_ext = ".bai";
    char *tbi_ext = ".tbi";
    char *crai_ext = ".crai";
    char *fai_ext = ".fai";

    if (!fn)
        return 0;

    if (hisremote(fn)) {
        for (i = strlen(fn) - 1; i >= 0; --i)
            if (fn[i] == '/') {
                fn_tmp = (char *)&fn[i+1];
                break;
            }
    } else {
        // Borrowed from hopen_fd_fileuri()
        if (strncmp(fn, "file://localhost/", 17) == 0) fn_tmp = fn + 16;
        else if (strncmp(fn, "file:///", 8) == 0) fn_tmp = fn + 7;
        else fn_tmp = fn;
#if defined(_WIN32) || defined(__MSYS__)
        // For cases like C:/foo
        if (fn_tmp[0] == '/' && fn_tmp[1] && fn_tmp[2] == ':' && fn_tmp[3] == '/')
            fn_tmp++;
#endif
    }

    if (!fn_tmp) return 0;
    hts_log_info("Using alignment file '%s'", fn_tmp);
    l_fn = strlen(fn_tmp); l_ext = 5;
    fnidx_tmp = (char*)calloc(l_fn + l_ext + 1, 1);
    if (!fnidx_tmp) return 0;

    struct stat sbuf;

    // Try alignment.bam.csi first
    strcpy(fnidx_tmp, fn_tmp); strcpy(fnidx_tmp + l_fn, csi_ext);
    if(stat(fnidx_tmp, &sbuf) == 0) {
        *fnidx = fnidx_tmp;
        return 1;
    } else { // Then try alignment.csi
        for (i = l_fn - 1; i > 0; --i)
            if (fnidx_tmp[i] == '.') {
                strcpy(fnidx_tmp + i, csi_ext);
                if(stat(fnidx_tmp, &sbuf) == 0) {
                    *fnidx = fnidx_tmp;
                    return 1;
                }
                break;
            }
    }
    if (fmt == HTS_FMT_BAI) {
        // Next, try alignment.bam.bai
        strcpy(fnidx_tmp, fn_tmp); strcpy(fnidx_tmp + l_fn, bai_ext);
        if(stat(fnidx_tmp, &sbuf) == 0) {
            *fnidx = fnidx_tmp;
            return 1;
        } else { // And finally, try alignment.bai
            for (i = l_fn - 1; i > 0; --i)
                if (fnidx_tmp[i] == '.') {
                    strcpy(fnidx_tmp + i, bai_ext);
                    if(stat(fnidx_tmp, &sbuf) == 0) {
                        *fnidx = fnidx_tmp;
                        return 1;
                    }
                    break;
                }
        }
    } else if (fmt == HTS_FMT_TBI) { // Or .tbi
        strcpy(fnidx_tmp, fn_tmp); strcpy(fnidx_tmp + l_fn, tbi_ext);
        if(stat(fnidx_tmp, &sbuf) == 0) {
            *fnidx = fnidx_tmp;
            return 1;
        } else {
            for (i = l_fn - 1; i > 0; --i)
                if (fnidx_tmp[i] == '.') {
                    strcpy(fnidx_tmp + i, tbi_ext);
                    if(stat(fnidx_tmp, &sbuf) == 0) {
                        *fnidx = fnidx_tmp;
                        return 1;
                    }
                    break;
                }
        }
    } else if (fmt == HTS_FMT_CRAI) { // Or .crai
        strcpy(fnidx_tmp, fn_tmp); strcpy(fnidx_tmp + l_fn, crai_ext);
        if(stat(fnidx_tmp, &sbuf) == 0) {
            *fnidx = fnidx_tmp;
            return 1;
        } else {
            for (i = l_fn - 1; i > 0; --i)
                if (fnidx_tmp[i] == '.') {
                    strcpy(fnidx_tmp + i, crai_ext);
                    if(stat(fnidx_tmp, &sbuf) == 0) {
                        *fnidx = fnidx_tmp;
                        return 1;
                    }
                    break;
                }
        }
    } else if (fmt == HTS_FMT_FAI) { // Or .fai
        strcpy(fnidx_tmp, fn_tmp); strcpy(fnidx_tmp + l_fn, fai_ext);
        *fnidx = fnidx_tmp;
        if(stat(fnidx_tmp, &sbuf) == 0)
            return 1;
        else
            return 0;
    }

    free(fnidx_tmp);
    return 0;
}

static char *idx_filename(const char *fn, const char *ext, int download) {
    int ret, local_len;
    char *fnidx;
    const char *local_fn = NULL;
    kstring_t buffer = KS_INITIALIZE;

    // First try : append `ext` to `fn`
    if (!(fnidx = haddextension(&buffer, fn, 0, ext))) {
        free(buffer.s);
        return NULL;
    }
    if ((ret = idx_test_and_fetch(fnidx, &local_fn, &local_len, download)) == -1) {
        // Second try : replace suffix of `fn` with `ext`
        if (!(fnidx = haddextension(&buffer, fn, 1, ext))) {
            free(buffer.s);
            return NULL;
        }
        ret = idx_test_and_fetch(fnidx, &local_fn, &local_len, download);
    }

    if (ret < 0) {
        free(buffer.s);
        return NULL;
    }

    memmove(fnidx, local_fn, local_len);
    fnidx[local_len] = 0;
    return fnidx;
}

char *hts_idx_getfn(const char *fn, const char *ext)
{
    return idx_filename(fn, ext, HTS_IDX_SAVE_REMOTE);
}

char *hts_idx_locatefn(const char *fn, const char *ext)
{
    return idx_filename(fn, ext, 0);
}

static hts_idx_t *idx_find_and_load(const char *fn, int fmt, int flags)
{
    char *fnidx = strstr(fn, HTS_IDX_DELIM);
    hts_idx_t *idx;

    if ( fnidx ) {
        char *fn2 = strdup(fn);
        if (!fn2) {
            hts_log_error("%s", strerror(errno));
            return NULL;
        }
        fn2[fnidx - fn] = '\0';
        fnidx += strlen(HTS_IDX_DELIM);
        idx = hts_idx_load3(fn2, fnidx, fmt, flags);
        free(fn2);
        return idx;
    }

    if (hts_idx_check_local(fn, fmt, &fnidx) == 0 && hisremote(fn)) {
        if (flags & HTS_IDX_SAVE_REMOTE) {
            fnidx = idx_filename(fn, ".csi", HTS_IDX_SAVE_REMOTE);
            if (!fnidx) {
                switch (fmt) {
                case HTS_FMT_BAI: fnidx = idx_filename(fn, ".bai", HTS_IDX_SAVE_REMOTE); break;
                case HTS_FMT_TBI: fnidx = idx_filename(fn, ".tbi", HTS_IDX_SAVE_REMOTE); break;
                default: break;
                }
            }
        } else {
            fnidx = idx_filename(fn, ".csi", 0);
            if (!fnidx) {
                switch (fmt) {
                case HTS_FMT_BAI: fnidx = idx_filename(fn, ".bai", 0); break;
                case HTS_FMT_TBI: fnidx = idx_filename(fn, ".tbi", 0); break;
                default: break;
                }
            }
        }
    }
    if (!fnidx) {
        if (!(flags & HTS_IDX_SILENT_FAIL))
            hts_log_error("Could not retrieve index file for '%s'", fn);
        return 0;
    }

    if (flags & HTS_IDX_SAVE_REMOTE)
        idx = hts_idx_load3(fn, fnidx, fmt, flags);
    else
        idx = idx_read(fnidx);
    free(fnidx);
    return idx;
}

hts_idx_t *hts_idx_load(const char *fn, int fmt) {
    return idx_find_and_load(fn, fmt, 1);
}

hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)
{
    return hts_idx_load3(fn, fnidx, 0, 0);
}

hts_idx_t *hts_idx_load3(const char *fn, const char *fnidx, int fmt, int flags)
{
    const char *local_fn = NULL;
    char *local_fnidx = NULL;
    int local_len;
    if (!fnidx)
        return idx_find_and_load(fn, fmt, flags);

    // Check that the index file is up to date, the main file might have changed
    struct stat stat_idx,stat_main;
    int remote_fn = hisremote(fn), remote_fnidx = hisremote(fnidx);
    if ( !remote_fn && !remote_fnidx
         && !stat(fn, &stat_main) && !stat(fnidx, &stat_idx) )
    {
        if ( stat_idx.st_mtime < stat_main.st_mtime )
            hts_log_warning("The index file is older than the data file: %s", fnidx);
    }

    if (remote_fnidx && (flags & HTS_IDX_SAVE_REMOTE))
    {
        int ret = idx_test_and_fetch(fnidx, &local_fn, &local_len, 1);
        if (ret == 0) {
            local_fnidx = strdup(local_fn);
            if (local_fnidx) {
                local_fnidx[local_len] = '\0';
                fnidx = local_fnidx;
            }
        }
    }

    hts_idx_t *idx = idx_read(fnidx);
    if (!idx && !(flags & HTS_IDX_SILENT_FAIL))
        hts_log_error("Could not load local index file '%s'%s%s", fnidx,
                      errno ? " : " : "", errno ? strerror(errno) : "");


    free(local_fnidx);

    return idx;
}



/**********************
 ***     Memory     ***
 **********************/

/* For use with hts_expand macros *only* */
HTSLIB_EXPORT
size_t hts_realloc_or_die(size_t n, size_t m, size_t m_sz, size_t size,
                          int clear, void **ptr, const char *func) {
    /* If new_m and size are both below this limit, multiplying them
       together can't overflow */
    const size_t safe = (size_t) 1 << (sizeof(size_t) * 4);
    void *new_ptr;
    size_t bytes, new_m;

    new_m = n;
    kroundup_size_t(new_m);

    bytes = size * new_m;

    /* Check for overflow.  Both ensure that new_m will fit in m (we make the
       pessimistic assumption that m is signed), and that bytes has not
       wrapped around. */
    if (new_m > (((size_t) 1 << (m_sz * 8 - 1)) - 1)
        || ((size > safe || new_m > safe)
            && bytes / new_m != size)) {
        errno = ENOMEM;
        goto die;
    }

    new_ptr = realloc(*ptr, bytes);
    if (new_ptr == NULL) goto die;

    if (clear) {
        if (new_m > m) {
            memset((char *) new_ptr + m * size, 0, (new_m - m) * size);
        }
    }

    *ptr = new_ptr;

    return new_m;

 die:
    hts_log_error("%s", strerror(errno));
    exit(1);
}

/*
 * Companion to hts_resize() macro that does the actual allocation.
 *
 * Somewhat complicated as hts_resize() needs to write the new allocated
 * size back into *size_in_out, and the value pointed to may either be
 * int32_t, uint32_t or size_t depending on which array is being resized.
 * This is solved by making `size_in_out` a void pointer, getting the macro
 * to pass in the size of the item pointed to (in `size_sz`) and then using
 * an appropriate cast (based on the value of size_sz).  The function
 * ensures that the maximum size will be storable in a signed type of
 * the given size so storing to an int32_t should work correctly.
 *
 * Assumes that sizeof(uint32_t) and sizeof(int32_t) is 4,
 * sizeof(uint64_t) and sizeof(int64_t) is 8 and sizeof(size_t) is
 * either 4 or 8.  It also assumes casting from unsigned to signed will
 * work as long as the top bit isn't set.
 */

int hts_resize_array_(size_t item_size, size_t num, size_t size_sz,
                      void *size_in_out, void **ptr_in_out, int flags,
                      const char *func) {
    /* If new_size and item_size are both below this limit, multiplying them
       together can't overflow */
    const size_t safe = (size_t) 1 << (sizeof(size_t) * 4);
    void *new_ptr;
    size_t bytes, new_size;

    new_size = num;
    kroundup_size_t(new_size);
    bytes = item_size * new_size;

    /* Check for overflow.  Both ensure that alloc will fit in alloc_in_out (we
       make the pessimistic assumption that *alloc_in_out is signed), and that
       bytes has not wrapped around. */

    if ((new_size > (((size_t) 1 << (size_sz * 8 - 1)) - 1))
        || (((item_size > safe) || (new_size > safe))
            && bytes / new_size != item_size)) {
        hts_log(HTS_LOG_ERROR, func, "Memory allocation too large");
        errno = ENOMEM;
        return -1;
    }

    new_ptr = realloc(*ptr_in_out, bytes);
    if (new_ptr == NULL) {
        int save_errno = errno;
        hts_log(HTS_LOG_ERROR, func, "%s", strerror(errno));
        errno = save_errno;
        return -1;
    }

    if (flags & HTS_RESIZE_CLEAR) {
        size_t old_size;
        switch (size_sz) {
        case 4: old_size = *((uint32_t *) size_in_out); break;
        case 8: old_size = *((uint64_t *) size_in_out); break;
        default: abort();
        }
        if (new_size > old_size) {
            memset((char *) new_ptr + old_size * item_size, 0,
                   (new_size - old_size) * item_size);
        }
    }

    switch (size_sz) {
    case 4: *((uint32_t *) size_in_out) = new_size; break;
    case 8: *((uint64_t *) size_in_out) = new_size; break;
    default: abort();
    }

    *ptr_in_out = new_ptr;
    return 0;
}

void hts_lib_shutdown()
{
    hfile_shutdown(1);
}

void hts_free(void *ptr) {
    free(ptr);
}

void hts_set_log_level(enum htsLogLevel level)
{
    hts_verbose = level;
}

enum htsLogLevel hts_get_log_level()
{
    return hts_verbose;
}

static char get_severity_tag(enum htsLogLevel severity)
{
    switch (severity) {
    case HTS_LOG_ERROR:
        return 'E';
    case HTS_LOG_WARNING:
        return 'W';
    case HTS_LOG_INFO:
        return 'I';
    case HTS_LOG_DEBUG:
        return 'D';
    case HTS_LOG_TRACE:
        return 'T';
    default:
        break;
    }

    return '*';
}

void hts_log(enum htsLogLevel severity, const char *context, const char *format, ...)
{
    int save_errno = errno;
    if (severity <= hts_verbose) {
        va_list argptr;

        fprintf(stderr, "[%c::%s] ", get_severity_tag(severity), context);

        va_start(argptr, format);
        vfprintf(stderr, format, argptr);
        va_end(argptr);

        fprintf(stderr, "\n");
    }
    errno = save_errno;
}
