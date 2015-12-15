/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2015 Genome Research Ltd.
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

#include <config.h>

#include <zlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "cram/cram.h"
#include "htslib/hfile.h"
#include "version.h"
#include "hts_internal.h"

#include "htslib/kseq.h"
#define KS_BGZF 1
#if KS_BGZF
    // bgzf now supports gzip-compressed files, the gzFile branch can be removed
    KSTREAM_INIT2(, BGZF*, bgzf_read, 65536)
#else
    KSTREAM_INIT2(, gzFile, gzread, 16384)
#endif

#include "htslib/khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = 3;

const char *hts_version()
{
    return HTS_VERSION;
}

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

const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

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
        return sequence_data;

    case vcf:
    case bcf:
        return variant_data;

    case bai:
    case crai:
    case csi:
    case gzi:
    case tbi:
        return index_file;

    case bed:
        return region_list;

    case unknown_format:
    case binary_format:
    case text_format:
    case format_maximum:
        break;
    }

    return unknown_category;
}

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
    // Typically at most a couple of hundred bytes of input are required
    // to get a few bytes of output from inflate(), so hopefully this buffer
    // size suffices in general.
    unsigned char buffer[512];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return 0;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return 0;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}

// Parse "x.y" text, taking care because the string is not NUL-terminated
// and filling in major/minor only when the digits are followed by a delimiter,
// so we don't misread "1.10" as "1.1" due to reaching the end of the buffer.
static void
parse_version(htsFormat *fmt, const unsigned char *u, const unsigned char *ulim)
{
    const char *str  = (const char *) u;
    const char *slim = (const char *) ulim;
    const char *s;

    fmt->version.major = fmt->version.minor = -1;

    for (s = str; s < slim; s++) if (!isdigit(*s)) break;
    if (s < slim) {
        fmt->version.major = atoi(str);
        if (*s == '.') {
            str = &s[1];
            for (s = str; s < slim; s++) if (!isdigit(*s)) break;
            if (s < slim)
                fmt->version.minor = atoi(str);
        }
        else
            fmt->version.minor = 0;
    }
}

int hts_detect_format(hFILE *hfile, htsFormat *fmt)
{
    unsigned char s[21];
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return -1;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few bytes.
        fmt->compression = (len >= 18 && (s[3] & 4) &&
                            memcmp(&s[12], "BC\2\0", 4) == 0)? bgzf : gzip;
        len = decompress_peek(hfile, s, sizeof s);
    }
    else {
        fmt->compression = no_compression;
        len = hpeek(hfile, s, sizeof s);
    }
    if (len < 0) return -1;

    fmt->compression_level = -1;
    fmt->specific = NULL;

    if (len >= 6 && memcmp(s,"CRAM",4) == 0 && s[4]>=1 && s[4]<=3 && s[5]<=1) {
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
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
    }
    else if (len >= 16 && memcmp(s, "##fileformat=VCF", 16) == 0) {
        fmt->category = variant_data;
        fmt->format = vcf;
        if (len >= 21 && s[16] == 'v')
            parse_version(fmt, &s[17], &s[len]);
        else
            fmt->version.major = fmt->version.minor = -1;
        return 0;
    }
    else if (len >= 4 && s[0] == '@' &&
             (memcmp(s, "@HD\t", 4) == 0 || memcmp(s, "@SQ\t", 4) == 0 ||
              memcmp(s, "@RG\t", 4) == 0 || memcmp(s, "@PG\t", 4) == 0)) {
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
    else {
        // Various possibilities for tab-delimited text:
        // .crai   (gzipped tab-delimited six columns: seqid 5*number)
        // .bed    ([3..12] tab-delimited columns)
        // .bedpe  (>= 10 tab-delimited columns)
        // .sam    (tab-delimited >= 11 columns: seqid number seqid...)
        // FIXME For now, assume it's SAM
        fmt->category = sequence_data;
        fmt->format = sam;
        fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }

    fmt->category = unknown_category;
    fmt->format = unknown_format;
    fmt->version.major = fmt->version.minor = -1;
    fmt->compression = no_compression;
    return 0;
}

char *hts_format_description(const htsFormat *format)
{
    kstring_t str = { 0, 0, NULL };

    switch (format->format) {
    case sam:   kputs("SAM", &str); break;
    case bam:   kputs("BAM", &str); break;
    case cram:  kputs("CRAM", &str); break;
    case vcf:   kputs("VCF", &str); break;
    case bcf:
        if (format->version.major == 1) kputs("Legacy BCF", &str);
        else kputs("BCF", &str);
        break;
    case bai:   kputs("BAI", &str); break;
    case crai:  kputs("CRAI", &str); break;
    case csi:   kputs("CSI", &str); break;
    case tbi:   kputs("Tabix", &str); break;
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
        case sam:
        case crai:
        case vcf:
        case bed:
            kputs(" text", &str);
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
    char smode[102], *cp, *cp2, *mode_c;
    htsFile *fp = NULL;
    hFILE *hfile;
    char fmt_code = '\0';

    strncpy(smode, mode, 100);
    smode[100]=0;
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
    *cp2++ = 0;

    // Set or reset the format code if opts->format is used
    if (fmt && fmt->format != unknown_format)
        *mode_c = "\0g\0\0b\0c\0\0b\0g\0\0"[fmt->format];

    hfile = hopen(fn, smode);
    if (hfile == NULL) goto error;

    fp = hts_hopen(hfile, fn, smode);
    if (fp == NULL) goto error;

    if (fmt && fmt->specific)
        if (hts_opt_apply(fp, fmt->specific) != 0)
            goto error;

    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (hfile)
        hclose_abruptly(hfile);

    return NULL;
}

htsFile *hts_open(const char *fn, const char *mode) {
    return hts_open_format(fn, mode, NULL);
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

    else if (strcmp(o->arg, "slices_per_container") == 0 ||
             strcmp(o->arg, "SLICES_PER_CONTAINER") == 0)
        o->opt = CRAM_OPT_SLICES_PER_CONTAINER, o->val.i = atoi(val);

    else if (strcmp(o->arg, "embed_ref") == 0 ||
             strcmp(o->arg, "EMBED_REF") == 0)
        o->opt = CRAM_OPT_EMBED_REF, o->val.i = atoi(val);

    else if (strcmp(o->arg, "no_ref") == 0 ||
             strcmp(o->arg, "NO_REF") == 0)
        o->opt = CRAM_OPT_NO_REF, o->val.i = atoi(val);

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

    else if (strcmp(o->arg, "required_fields") == 0 ||
             strcmp(o->arg, "REQUIRED_FIELDS") == 0)
        o->opt = CRAM_OPT_REQUIRED_FIELDS, o->val.i = strtol(val, NULL, 0);

    else {
        fprintf(stderr, "Unknown option '%s'\n", o->arg);
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

    for (; opts;  opts = (last=opts)->next)
        if (hts_set_opt(fp,  opts->opt,  opts->val) != 0)
            return -1;

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
    const char *cp;

    if (!(cp = strchr(str, ',')))
        cp = str+strlen(str);

    format->version.minor = 0; // unknown
    format->version.major = 0; // unknown

    if (strncmp(str, "sam", cp-str) == 0) {
        format->category          = sequence_data;
        format->format            = sam;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strncmp(str, "bam", cp-str) == 0) {
        format->category          = sequence_data;
        format->format            = bam;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else if (strncmp(str, "cram", cp-str) == 0) {
        format->category          = sequence_data;
        format->format            = cram;
        format->compression       = custom;
        format->compression_level = -1;
    } else if (strncmp(str, "vcf", cp-str) == 0) {
        format->category          = variant_data;
        format->format            = vcf;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strncmp(str, "bcf", cp-str) == 0) {
        format->category          = variant_data;
        format->format            = bcf;
        format->compression       = bgzf;
        format->compression_level = -1;
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


htsFile *hts_hopen(struct hFILE *hfile, const char *fn, const char *mode)
{
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
        if (hts_detect_format(hfile, &fp->format) < 0) goto error;
    }
    else if (strchr(simple_mode, 'w') || strchr(simple_mode, 'a')) {
        htsFormat *fmt = &fp->format;
        fp->is_write = 1;

        if (strchr(simple_mode, 'b')) fmt->format = binary_format;
        else if (strchr(simple_mode, 'c')) fmt->format = cram;
        else fmt->format = text_format;

        if (strchr(simple_mode, 'z')) fmt->compression = bgzf;
        else if (strchr(simple_mode, 'g')) fmt->compression = gzip;
        else if (strchr(simple_mode, 'u')) fmt->compression = no_compression;
        else {
            // No compression mode specified, set to the default for the format
            switch (fmt->format) {
            case binary_format: fmt->compression = bgzf; break;
            case cram: fmt->compression = custom; break;
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
    else goto error;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
        if (fp->fp.bgzf == NULL) goto error;
        fp->is_bin = 1;
        break;

    case cram:
        fp->fp.cram = cram_dopen(hfile, fn, simple_mode);
        if (fp->fp.cram == NULL) goto error;
        if (!fp->is_write)
            cram_set_option(fp->fp.cram, CRAM_OPT_DECODE_MD, 1);
        fp->is_cram = 1;
        break;

    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = bgzf_hopen(hfile, simple_mode);
        #else
            // TODO Implement gzip hFILE adaptor
            hclose(hfile); // This won't work, especially for stdin
            gzFile gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
        #endif
            if (gzfp) fp->fp.voidp = ks_init(gzfp);
            else goto error;
        }
        else if (fp->format.compression != no_compression) {
            fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
            if (fp->fp.bgzf == NULL) goto error;
        }
        else
            fp->fp.hfile = hfile;
        break;

    default:
        goto error;
    }

    if (opts)
        hts_process_opts(fp, opts);

    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (fp) {
        free(fp->fn);
        free(fp->fn_aux);
        free(fp);
    }
    return NULL;
}

int hts_close(htsFile *fp)
{
    int ret, save;

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
                fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
                break;
            case 0:  /* not at EOF, but may not have wanted all seqs */
            default: /* case 1, expected EOF */
                break;
            }
        }
        ret = cram_close(fp->fp.cram);
        break;

    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = bgzf_close(gzfp);
        #else
            gzFile gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = gzclose(gzfp);
        #endif
            ks_destroy((kstream_t*)fp->fp.voidp);
        }
        else if (fp->format.compression != no_compression)
            ret = bgzf_close(fp->fp.bgzf);
        else
            ret = hclose(fp->fp.hfile);
        break;

    default:
        ret = -1;
        break;
    }

    save = errno;
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    errno = save;
    return ret;
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
    case gzi:  return "gzi";
    case tbi:  return "tbi";
    case bed:  return "bed";
    default:   return "?";
    }
}

int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    if (opt == HTS_OPT_NTHREADS) {
        va_start(args, opt);
        int nthreads = va_arg(args, int);
        va_end(args);
        return hts_set_threads(fp, nthreads);
    }

    if (fp->format.format != cram)
        return 0;

    va_start(args, opt);
    r = cram_set_voption(fp->fp.cram, opt, args);
    va_end(args);

    return r;
}

int hts_set_threads(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf) {
        return bgzf_mt(fp->fp.bgzf, n, 256);
    } else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_NTHREADS, n);
    }
    else return 0;
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

// For VCF/BCF backward sweeper. Not exposing these functions because their
// future is uncertain. Things will probably have to change with hFILE...
BGZF *hts_get_bgzfp(htsFile *fp)
{
    if ( fp->is_bin )
        return fp->fp.bgzf;
    else
        return ((kstream_t*)fp->fp.voidp)->f;
}
int hts_useek(htsFile *fp, long uoffset, int where)
{
    if ( fp->is_bin )
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
    {
        ks_rewind((kstream_t*)fp->fp.voidp);
        ((kstream_t*)fp->fp.voidp)->seek_pos = uoffset;
        return bgzf_useek(((kstream_t*)fp->fp.voidp)->f, uoffset, where);
    }
}
long hts_utell(htsFile *fp)
{
    if ( fp->is_bin )
        return bgzf_utell(fp->fp.bgzf);
    else
        return ((kstream_t*)fp->fp.voidp)->seek_pos;
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
    int ret, dret;
    ret = ks_getuntil((kstream_t*)fp->fp.voidp, delimiter, str, &dret);
    ++fp->lineno;
    return ret;
}

char **hts_readlist(const char *string, int is_file, int *_n)
{
    int m = 0, n = 0, dret;
    char **s = 0;
    if ( is_file )
    {
#if KS_BGZF
        BGZF *fp = bgzf_open(string, "r");
#else
        gzFile fp = gzopen(string, "r");
#endif
        if ( !fp ) return NULL;

        kstream_t *ks;
        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0)
        {
            if (str.l == 0) continue;
            n++;
            hts_expand(char*,n,m,s);
            s[n-1] = strdup(str.s);
        }
        ks_destroy(ks);
#if KS_BGZF
        bgzf_close(fp);
#else
        gzclose(fp);
#endif
        free(str.s);
    }
    else
    {
        const char *q = string, *p = string;
        while ( 1 )
        {
            if (*p == ',' || *p == 0)
            {
                n++;
                hts_expand(char*,n,m,s);
                s[n-1] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n-1], q, p - q);
                q = p + 1;
            }
            if ( !*p ) break;
            p++;
        }
    }
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
}

char **hts_readlines(const char *fn, int *_n)
{
    int m = 0, n = 0, dret;
    char **s = 0;
#if KS_BGZF
    BGZF *fp = bgzf_open(fn, "r");
#else
    gzFile fp = gzopen(fn, "r");
#endif
    if ( fp ) { // read from file
        kstream_t *ks;
        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
            if (str.l == 0) continue;
            if (m == n) {
                m = m? m<<1 : 16;
                s = (char**)realloc(s, m * sizeof(char*));
            }
            s[n++] = strdup(str.s);
        }
        ks_destroy(ks);
        #if KS_BGZF
            bgzf_close(fp);
        #else
            gzclose(fp);
        #endif
        s = (char**)realloc(s, n * sizeof(char*));
        free(str.s);
    } else if (*fn == ':') { // read from string
        const char *q, *p;
        for (q = p = fn + 1;; ++p)
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    s = (char**)realloc(s, m * sizeof(char*));
                }
                s[n] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
    } else return 0;
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
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
    if (hts_detect_format(f, &fmt) < 0) { hclose_abruptly(f); return 0; }
    if (hclose(f) < 0) return 0;

    switch (fmt.format) {
    case vcf: return (fmt.compression == no_compression)? FT_VCF : FT_VCF_GZ;
    case bcf: return (fmt.compression == no_compression)? FT_BCF : FT_BCF_GZ;
    default:  return 0;
    }
}

/****************
 *** Indexing ***
 ****************/

#define HTS_MIN_MARKER_DIST 0x10000

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

#define pair64_lt(a,b) ((a).u < (b).u)

#include "htslib/ksort.h"
KSORT_INIT(_off, hts_pair64_t, pair64_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta;
    struct {
        uint32_t last_bin, save_bin;
        int last_coor, last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

static inline void insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
    khint_t k;
    bins_t *l;
    int absent;
    k = kh_put(bin, b, bin, &absent);
    l = &kh_value(b, k);
    if (absent) {
        l->m = 1; l->n = 0;
        l->list = (hts_pair64_t*)calloc(l->m, sizeof(hts_pair64_t));
    }
    if (l->n == l->m) {
        l->m <<= 1;
        l->list = (hts_pair64_t*)realloc(l->list, l->m * sizeof(hts_pair64_t));
    }
    l->list[l->n].u = beg;
    l->list[l->n++].v = end;
}

static inline void insert_to_l(lidx_t *l, int64_t _beg, int64_t _end, uint64_t offset, int min_shift)
{
    int i, beg, end;
    beg = _beg >> min_shift;
    end = (_end - 1) >> min_shift;
    if (l->m < end + 1) {
        int old_m = l->m;
        l->m = end + 1;
        kroundup32(l->m);
        l->offset = (uint64_t*)realloc(l->offset, l->m * sizeof(uint64_t));
        memset(l->offset + old_m, 0xff, 8 * (l->m - old_m)); // fill l->offset with (uint64_t)-1
    }
    if (beg == end) { // to save a loop in this case
        if (l->offset[beg] == (uint64_t)-1) l->offset[beg] = offset;
    } else {
        for (i = beg; i <= end; ++i)
            if (l->offset[i] == (uint64_t)-1) l->offset[i] = offset;
    }
    if (l->n < end + 1) l->n = end + 1;
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
    idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
    idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
    idx->z.last_coor = 0xffffffffu;
    if (n) {
        idx->n = idx->m = n;
        idx->bidx = (bidx_t**)calloc(n, sizeof(bidx_t*));
        if (idx->bidx == NULL) { free(idx); return NULL; }
        idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
        if (idx->lidx == NULL) { free(idx->bidx); free(idx); return NULL; }
    }
    return idx;
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
    bidx_t *bidx = idx->bidx[i];
    lidx_t *lidx = &idx->lidx[i];
    khint_t k;
    int l;
    uint64_t offset0 = 0;
    if (bidx) {
        k = kh_get(bin, bidx, META_BIN(idx));
        if (k != kh_end(bidx))
            offset0 = kh_val(bidx, k).list[0].u;
        for (l = 0; l < lidx->n && lidx->offset[l] == (uint64_t)-1; ++l)
            lidx->offset[l] = offset0;
    } else l = 1;
    for (; l < lidx->n; ++l) // fill missing values
        if (lidx->offset[l] == (uint64_t)-1)
            lidx->offset[l] = lidx->offset[l-1];
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

static void compress_binning(hts_idx_t *idx, int i)
{
    bidx_t *bidx = idx->bidx[i];
    khint_t k;
    int l, m;
    if (bidx == 0) return;
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
                    q->m = q->n + p->n;
                    kroundup32(q->m);
                    q->list = (hts_pair64_t*)realloc(q->list, q->m * sizeof(hts_pair64_t));
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
}

void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
    int i;
    if (idx == NULL || idx->z.finished) return; // do not run this function on an empty index or multiple times
    if (idx->z.save_tid >= 0) {
        insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
    }
    for (i = 0; i < idx->n; ++i) {
        update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
        compress_binning(idx, i);
    }
    idx->z.finished = 1;
}

int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
{
    int bin;
    if (tid<0) beg = -1, end = 0;
    if (tid >= idx->m) { // enlarge the index
        int32_t oldm = idx->m;
        idx->m = idx->m? idx->m<<1 : 2;
        idx->bidx = (bidx_t**)realloc(idx->bidx, idx->m * sizeof(bidx_t*));
        idx->lidx = (lidx_t*) realloc(idx->lidx, idx->m * sizeof(lidx_t));
        memset(&idx->bidx[oldm], 0, (idx->m - oldm) * sizeof(bidx_t*));
        memset(&idx->lidx[oldm], 0, (idx->m - oldm) * sizeof(lidx_t));
    }
    if (idx->n < tid + 1) idx->n = tid + 1;
    if (idx->z.finished) return 0;
    if (idx->z.last_tid != tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
        if ( tid>=0 && idx->n_no_coor )
        {
            if (hts_verbose >= 1) fprintf(stderr,"[E::%s] NO_COOR reads not in a single block at the end %d %d\n", __func__, tid,idx->z.last_tid);
            return -1;
        }
        if (tid>=0 && idx->bidx[tid] != 0)
        {
            if (hts_verbose >= 1) fprintf(stderr, "[E::%s] chromosome blocks not continuous\n", __func__);
            return -1;
        }
        idx->z.last_tid = tid;
        idx->z.last_bin = 0xffffffffu;
    } else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
        if (hts_verbose >= 1) fprintf(stderr, "[E::%s] unsorted positions\n", __func__);
        return -1;
    }
    if ( tid>=0 )
    {
        if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
        if ( is_mapped)
            insert_to_l(&idx->lidx[tid], beg, end, idx->z.last_off, idx->min_shift); // last_off points to the start of the current record
    }
    else idx->n_no_coor++;
    bin = hts_reg2bin(beg, end, idx->min_shift, idx->n_lvls);
    if ((int)idx->z.last_bin != bin) { // then possibly write the binning index
        if (idx->z.save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
            insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, idx->z.last_off);
        if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
            idx->z.off_end = idx->z.last_off;
            insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, idx->z.off_end);
            insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
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

static inline long idx_write(int is_bgzf, void *fp, const void *buf, long l)
{
    if (is_bgzf) return bgzf_write((BGZF*)fp, buf, l);
    else return (long)fwrite(buf, 1, l, (FILE*)fp);
}

static inline void swap_bins(bins_t *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        ed_swap_8p(&p->list[i].u);
        ed_swap_8p(&p->list[i].v);
    }
}

static void hts_idx_save_core(const hts_idx_t *idx, void *fp, int fmt)
{
    int32_t i, size, is_be;
    int is_bgzf = (fmt != HTS_FMT_BAI);
    is_be = ed_is_big();
    if (is_be) {
        uint32_t x = idx->n;
        idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
    } else idx_write(is_bgzf, fp, &idx->n, 4);
    if (fmt == HTS_FMT_TBI && idx->l_meta) idx_write(is_bgzf, fp, idx->meta, idx->l_meta);
    for (i = 0; i < idx->n; ++i) {
        khint_t k;
        bidx_t *bidx = idx->bidx[i];
        lidx_t *lidx = &idx->lidx[i];
        // write binning index
        size = bidx? kh_size(bidx) : 0;
        if (is_be) { // big endian
            uint32_t x = size;
            idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
        } else idx_write(is_bgzf, fp, &size, 4);
        if (bidx == 0) goto write_lidx;
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
            bins_t *p;
            if (!kh_exist(bidx, k)) continue;
            p = &kh_value(bidx, k);
            if (is_be) { // big endian
                uint32_t x;
                x = kh_key(bidx, k); idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
                if (fmt == HTS_FMT_CSI) {
                    uint64_t y = kh_val(bidx, k).loff;
                    idx_write(is_bgzf, fp, ed_swap_4p(&y), 8);
                }
                x = p->n; idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
                swap_bins(p);
                idx_write(is_bgzf, fp, p->list, 16 * p->n);
                swap_bins(p);
            } else {
                idx_write(is_bgzf, fp, &kh_key(bidx, k), 4);
                if (fmt == HTS_FMT_CSI) idx_write(is_bgzf, fp, &kh_val(bidx, k).loff, 8);
                //int j;for(j=0;j<p->n;++j)fprintf(stderr,"%d,%llx,%d,%llx:%llx\n",kh_key(bidx,k),kh_val(bidx, k).loff,j,p->list[j].u,p->list[j].v);
                idx_write(is_bgzf, fp, &p->n, 4);
                idx_write(is_bgzf, fp, p->list, p->n << 4);
            }
        }
write_lidx:
        if (fmt != HTS_FMT_CSI) {
            if (is_be) {
                int32_t x = lidx->n;
                idx_write(is_bgzf, fp, ed_swap_4p(&x), 4);
                for (x = 0; x < lidx->n; ++x) ed_swap_8p(&lidx->offset[x]);
                idx_write(is_bgzf, fp, lidx->offset, lidx->n << 3);
                for (x = 0; x < lidx->n; ++x) ed_swap_8p(&lidx->offset[x]);
            } else {
                idx_write(is_bgzf, fp, &lidx->n, 4);
                idx_write(is_bgzf, fp, lidx->offset, lidx->n << 3);
            }
        }
    }
    if (is_be) { // write the number of reads without coordinates
        uint64_t x = idx->n_no_coor;
        idx_write(is_bgzf, fp, &x, 8);
    } else idx_write(is_bgzf, fp, &idx->n_no_coor, 8);
}

int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
{
    int ret, save;
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
    if (fnidx == NULL) return hts_idx_save(idx, fn, fmt);

    if (fmt == HTS_FMT_CSI) {
        BGZF *fp;
        uint32_t x[3];
        int is_be, i;
        is_be = ed_is_big();
        fp = bgzf_open(fnidx, "w");
        if (fp == NULL) return -1;
        bgzf_write(fp, "CSI\1", 4);
        x[0] = idx->min_shift; x[1] = idx->n_lvls; x[2] = idx->l_meta;
        if (is_be) {
            for (i = 0; i < 3; ++i)
                bgzf_write(fp, ed_swap_4p(&x[i]), 4);
        } else bgzf_write(fp, &x, 12);
        if (idx->l_meta) bgzf_write(fp, idx->meta, idx->l_meta);
        hts_idx_save_core(idx, fp, HTS_FMT_CSI);
        bgzf_close(fp);
    } else if (fmt == HTS_FMT_TBI) {
        BGZF *fp = bgzf_open(fnidx, "w");
        if (fp == NULL) return -1;
        bgzf_write(fp, "TBI\1", 4);
        hts_idx_save_core(idx, fp, HTS_FMT_TBI);
        bgzf_close(fp);
    } else if (fmt == HTS_FMT_BAI) {
        FILE *fp = fopen(fnidx, "w");
        if (fp == NULL) return -1;
        fwrite("BAI\1", 1, 4, fp);
        hts_idx_save_core(idx, fp, HTS_FMT_BAI);
        fclose(fp);
    } else abort();

    return 0;
}

static int hts_idx_load_core(hts_idx_t *idx, BGZF *fp, int fmt)
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
        for (j = 0; j < n; ++j) {
            khint_t k;
            if (bgzf_read(fp, &key, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&key);
            k = kh_put(bin, h, key, &absent);
            if (absent <= 0) return -3; // Duplicate bin number
            p = &kh_val(h, k);
            if (fmt == HTS_FMT_CSI) {
                if (bgzf_read(fp, &p->loff, 8) != 8) return -1;
                if (is_be) ed_swap_8p(&p->loff);
            } else p->loff = 0;
            if (bgzf_read(fp, &p->n, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&p->n);
            p->m = p->n;
            p->list = (hts_pair64_t*)malloc(p->m * sizeof(hts_pair64_t));
            if (p->list == NULL) return -2;
            if (bgzf_read(fp, p->list, p->n<<4) != p->n<<4) return -1;
            if (is_be) swap_bins(p);
        }
        if (fmt != HTS_FMT_CSI) { // load linear index
            int j;
            if (bgzf_read(fp, &l->n, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&l->n);
            l->m = l->n;
            l->offset = (uint64_t*)malloc(l->n * sizeof(uint64_t));
            if (l->offset == NULL) return -2;
            if (bgzf_read(fp, l->offset, l->n << 3) != l->n << 3) return -1;
            if (is_be) for (j = 0; j < l->n; ++j) ed_swap_8p(&l->offset[j]);
            for (j = 1; j < l->n; ++j) // fill missing values; may happen given older samtools and tabix
                if (l->offset[j] == 0) l->offset[j] = l->offset[j-1];
            update_loff(idx, i, 1);
        }
    }
    if (bgzf_read(fp, &idx->n_no_coor, 8) != 8) idx->n_no_coor = 0;
    if (is_be) ed_swap_8p(&idx->n_no_coor);
    return 0;
}

static hts_idx_t *hts_idx_load_local(const char *fn)
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
            if ((meta = (uint8_t*)malloc(x[2])) == NULL) goto fail;
            if (bgzf_read(fp, meta, x[2]) != x[2]) goto fail;
        }
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        if ((idx = hts_idx_init(n, HTS_FMT_CSI, 0, x[0], x[1])) == NULL) goto fail;
        idx->l_meta = x[2];
        idx->meta = meta;
        meta = NULL;
        if (hts_idx_load_core(idx, fp, HTS_FMT_CSI) < 0) goto fail;
    }
    else if (memcmp(magic, "TBI\1", 4) == 0) {
        uint32_t x[8];
        if (bgzf_read(fp, x, 32) != 32) goto fail;
        if (is_be) for (i = 0; i < 8; ++i) ed_swap_4p(&x[i]);
        if ((idx = hts_idx_init(x[0], HTS_FMT_TBI, 0, 14, 5)) == NULL) goto fail;
        idx->l_meta = 28 + x[7];
        if ((idx->meta = (uint8_t*)malloc(idx->l_meta)) == NULL) goto fail;
        memcpy(idx->meta, &x[1], 28);
        if (bgzf_read(fp, idx->meta + 28, x[7]) != x[7]) goto fail;
        if (hts_idx_load_core(idx, fp, HTS_FMT_TBI) < 0) goto fail;
    }
    else if (memcmp(magic, "BAI\1", 4) == 0) {
        uint32_t n;
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        idx = hts_idx_init(n, HTS_FMT_BAI, 0, 14, 5);
        if (hts_idx_load_core(idx, fp, HTS_FMT_BAI) < 0) goto fail;
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

void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)
{
    if (idx->meta) free(idx->meta);
    idx->l_meta = l_meta;
    if (is_copy) {
        idx->meta = (uint8_t*)malloc(l_meta);
        memcpy(idx->meta, meta, l_meta);
    } else idx->meta = meta;
}

uint8_t *hts_idx_get_meta(hts_idx_t *idx, int *l_meta)
{
    *l_meta = idx->l_meta;
    return idx->meta;
}

const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)
{
    if ( !idx->n )
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

int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)
{
    if ( idx->fmt == HTS_FMT_CRAI ) {
        *mapped = 0; *unmapped = 0;
        return -1;
    }

    bidx_t *h = idx->bidx[tid];
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
    return idx->n_no_coor;
}

/****************
 *** Iterator ***
 ****************/

static inline int reg2bins(int64_t beg, int64_t end, hts_itr_t *itr, int min_shift, int n_lvls)
{
    int l, t, s = min_shift + (n_lvls<<1) + n_lvls;
    if (beg >= end) return 0;
    if (end >= 1LL<<s) end = 1LL<<s;
    for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        int b, e, n, i;
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

hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
{
    int i, n_off, l, bin;
    hts_pair64_t *off;
    khint_t k;
    bidx_t *bidx;
    uint64_t min_off;
    hts_itr_t *iter = 0;
    if (tid < 0) {
        int finished0 = 0;
        uint64_t off0 = (uint64_t)-1;
        khint_t k;
        switch (tid) {
        case HTS_IDX_START:
            // Find the smallest offset, note that sequence ids may not be ordered sequentially
            for (i=0; i<idx->n; i++)
            {
                bidx = idx->bidx[i];
                k = kh_get(bin, bidx, META_BIN(idx));
                if (k == kh_end(bidx)) continue;
                if ( off0 > kh_val(bidx, k).list[0].u ) off0 = kh_val(bidx, k).list[0].u;
            }
            if ( off0==(uint64_t)-1 && idx->n_no_coor ) off0 = 0; // only no-coor reads in this bam
            break;

        case HTS_IDX_NOCOOR:
            if ( idx->n>0 )
            {
                bidx = idx->bidx[idx->n - 1];
                k = kh_get(bin, bidx, META_BIN(idx));
                if (k != kh_end(bidx)) off0 = kh_val(bidx, k).list[0].v;
            }
            if ( off0==(uint64_t)-1 && idx->n_no_coor ) off0 = 0; // only no-coor reads in this bam
            break;

        case HTS_IDX_REST:
            off0 = 0;
            break;

        case HTS_IDX_NONE:
            finished0 = 1;
            off0 = 0;
            break;

        default:
            return 0;
        }
        if (off0 != (uint64_t)-1) {
            iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
            iter->read_rest = 1;
            iter->finished = finished0;
            iter->curr_off = off0;
            iter->readrec = readrec;
            return iter;
        } else return 0;
    }

    if (beg < 0) beg = 0;
    if (end < beg) return 0;
    if (tid >= idx->n || (bidx = idx->bidx[tid]) == NULL) return 0;

    iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
    iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;
    iter->readrec = readrec;

    // compute min_off
    bin = hts_bin_first(idx->n_lvls) + (beg>>idx->min_shift);
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
    // retrieve bins
    reg2bins(beg, end, iter, idx->min_shift, idx->n_lvls);
    for (i = n_off = 0; i < iter->bins.n; ++i)
        if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx))
            n_off += kh_value(bidx, k).n;
    if (n_off == 0) return iter;
    off = (hts_pair64_t*)calloc(n_off, sizeof(hts_pair64_t));
    for (i = n_off = 0; i < iter->bins.n; ++i) {
        if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx)) {
            int j;
            bins_t *p = &kh_value(bidx, k);
            for (j = 0; j < p->n; ++j)
                if (p->list[j].v > min_off) off[n_off++] = p->list[j];
        }
    }
    if (n_off == 0) {
        free(off); return iter;
    }
    ks_introsort(_off, n_off, off);
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
    return iter;
}

void hts_itr_destroy(hts_itr_t *iter)
{
    if (iter) { free(iter->off); free(iter->bins.a); free(iter); }
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
    int decimals = 0, e = 0, lost = 0;
    char sign = '+', esign = '+';
    const char *s;

    while (isspace(*str)) str++;
    s = str;

    if (*s == '+' || *s == '-') sign = *s++;
    while (*s)
        if (isdigit(*s)) n = push_digit(n, *s++);
        else if (*s == ',' && (flags & HTS_PARSE_THOUSANDS_SEP)) s++;
        else break;

    if (*s == '.') {
        s++;
        while (isdigit(*s)) decimals++, n = push_digit(n, *s++);
    }

    if (*s == 'E' || *s == 'e') {
        s++;
        if (*s == '+' || *s == '-') esign = *s++;
        while (isdigit(*s)) e = push_digit(e, *s++);
        if (esign == '-') e = -e;
    }

    e -= decimals;
    while (e > 0) n *= 10, e--;
    while (e < 0) lost += n % 10, n /= 10, e++;

    if (lost > 0 && hts_verbose >= 3)
        fprintf(stderr, "[W::%s] discarding fractional part of %.*s\n",
                __func__, (int)(s - str), str);

    if (strend) *strend = (char *) s;
    else if (*s && hts_verbose >= 2)
        fprintf(stderr, "[W::%s] ignoring unknown characters after %.*s[%s]\n",
                __func__, (int)(s - str), str, s);

    return (sign == '+')? n : -n;
}

const char *hts_parse_reg(const char *s, int *beg, int *end)
{
    char *hyphen;
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = INT_MAX;
        return s + strlen(s);
    }

    *beg = hts_parse_decimal(colon+1, &hyphen, HTS_PARSE_THOUSANDS_SEP) - 1;
    if (*beg < 0) *beg = 0;

    if (*hyphen == '\0') *end = INT_MAX;
    else if (*hyphen == '-') *end = hts_parse_decimal(hyphen+1, NULL, HTS_PARSE_THOUSANDS_SEP);
    else return NULL;

    if (*beg >= *end) return NULL;
    return colon;
}

hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)
{
    int tid, beg, end;
    const char *q;

    if (strcmp(reg, ".") == 0)
        return itr_query(idx, HTS_IDX_START, 0, 0, readrec);
    else if (strcmp(reg, "*") == 0)
        return itr_query(idx, HTS_IDX_NOCOOR, 0, 0, readrec);

    q = hts_parse_reg(reg, &beg, &end);
    if (q) {
        char *tmp = (char*)alloca(q - reg + 1);
        strncpy(tmp, reg, q - reg);
        tmp[q - reg] = 0;
        tid = getid(hdr, tmp);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        tid = getid(hdr, reg);
        beg = 0; end = INT_MAX;
    }

    if (tid < 0) return NULL;
    return itr_query(idx, tid, beg, end, readrec);
}

int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
{
    int ret, tid, beg, end;
    if (iter == NULL || iter->finished) return -1;
    if (iter->read_rest) {
        if (iter->curr_off) { // seek to the start
            bgzf_seek(fp, iter->curr_off, SEEK_SET);
            iter->curr_off = 0; // only seek once
        }
        ret = iter->readrec(fp, data, r, &tid, &beg, &end);
        if (ret < 0) iter->finished = 1;
        iter->curr_tid = tid;
        iter->curr_beg = beg;
        iter->curr_end = end;
        return ret;
    }
    if (iter->off == 0) return -1;
    for (;;) {
        if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
            if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
            if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
                bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
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

/**********************
 *** Retrieve index ***
 **********************/

static char *test_and_fetch(const char *fn)
{
    FILE *fp;
    if (hisremote(fn)) {
        const int buf_size = 1 * 1024 * 1024;
        hFILE *fp_remote;
        uint8_t *buf;
        int l;
        const char *p;
        for (p = fn + strlen(fn) - 1; p >= fn; --p)
            if (*p == '/') break;
        ++p; // p now points to the local file name
        // Attempt to open local file first
        if ((fp = fopen((char*)p, "rb")) != 0)
        {
            fclose(fp);
            return (char*)p;
        }
        // Attempt to open remote file. Stay quiet on failure, it is OK to fail when trying first .csi then .tbi index.
        if ((fp_remote = hopen(fn, "r")) == 0) return 0;
        if ((fp = fopen(p, "w")) == 0) {
            if (hts_verbose >= 1) fprintf(stderr, "[E::%s] fail to create file '%s' in the working directory\n", __func__, p);
            hclose_abruptly(fp_remote);
            return 0;
        }
        if (hts_verbose >= 3) fprintf(stderr, "[M::%s] downloading file '%s' to local directory\n", __func__, fn);
        buf = (uint8_t*)calloc(buf_size, 1);
        while ((l = hread(fp_remote, buf, buf_size)) > 0) fwrite(buf, 1, l, fp);
        free(buf);
        fclose(fp);
        if (hclose(fp_remote) != 0) fprintf(stderr, "[E::%s] fail to close remote file '%s'\n", __func__, fn);
        return (char*)p;
    } else {
        if ((fp = fopen(fn, "rb")) == 0) return 0;
        fclose(fp);
        return (char*)fn;
    }
}

char *hts_idx_getfn(const char *fn, const char *ext)
{
    int i, l_fn, l_ext;
    char *fnidx, *ret;
    l_fn = strlen(fn); l_ext = strlen(ext);
    fnidx = (char*)calloc(l_fn + l_ext + 1, 1);
    strcpy(fnidx, fn); strcpy(fnidx + l_fn, ext);
    if ((ret = test_and_fetch(fnidx)) == 0) {
        for (i = l_fn - 1; i > 0; --i)
            if (fnidx[i] == '.') break;
        strcpy(fnidx + i, ext);
        ret = test_and_fetch(fnidx);
    }
    if (ret == 0) {
        free(fnidx);
        return 0;
    }
    l_fn = strlen(ret);
    memmove(fnidx, ret, l_fn + 1);
    return fnidx;
}

hts_idx_t *hts_idx_load(const char *fn, int fmt)
{
    char *fnidx;
    hts_idx_t *idx;
    fnidx = hts_idx_getfn(fn, ".csi");
    if (! fnidx) fnidx = hts_idx_getfn(fn, fmt == HTS_FMT_BAI? ".bai" : ".tbi");
    if (fnidx == 0) return 0;

    idx = hts_idx_load2(fn, fnidx);
    free(fnidx);
    return idx;
}

hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)
{
    // Check that the index file is up to date, the main file might have changed
    struct stat stat_idx,stat_main;
    if ( !stat(fn, &stat_main) && !stat(fnidx, &stat_idx) )
    {
        if ( stat_idx.st_mtime < stat_main.st_mtime )
            fprintf(stderr, "Warning: The index file is older than the data file: %s\n", fnidx);
    }

    return hts_idx_load_local(fnidx);
}
