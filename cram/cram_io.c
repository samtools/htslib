/*
Copyright (c) 2012-2023 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * CRAM I/O primitives.
 *
 * - ITF8 encoding and decoding.
 * - Block based I/O
 * - Zlib inflating and deflating (memory)
 * - CRAM basic data structure reading and writing
 * - File opening / closing
 * - Reference sequence handling
 */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif
#ifdef HAVE_LIBLZMA
#ifdef HAVE_LZMA_H
#include <lzma.h>
#else
#include "../os/lzma_stub.h"
#endif
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>

#ifdef HAVE_LIBDEFLATE
#include <libdeflate.h>
#define crc32(a,b,c) libdeflate_crc32((a),(b),(c))
#endif

#include "cram.h"
#include "os.h"
#include "../htslib/hts.h"
#include "open_trace_file.h"

#if defined(HAVE_EXTERNAL_LIBHTSCODECS)
#include <htscodecs/rANS_static.h>
#include <htscodecs/rANS_static4x16.h>
#include <htscodecs/arith_dynamic.h>
#include <htscodecs/tokenise_name3.h>
#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/varint.h> // CRAM v4.0 variable-size integers
#else
#include "../htscodecs/htscodecs/rANS_static.h"
#include "../htscodecs/htscodecs/rANS_static4x16.h"
#include "../htscodecs/htscodecs/arith_dynamic.h"
#include "../htscodecs/htscodecs/tokenise_name3.h"
#include "../htscodecs/htscodecs/fqzcomp_qual.h"
#include "../htscodecs/htscodecs/varint.h"
#endif

//#define REF_DEBUG

#ifdef REF_DEBUG
#include <sys/syscall.h>
#define gettid() (int)syscall(SYS_gettid)

#define RP(...) fprintf (stderr, __VA_ARGS__)
#else
#define RP(...)
#endif

#include "../htslib/hfile.h"
#include "../htslib/bgzf.h"
#include "../htslib/faidx.h"
#include "../hts_internal.h"

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX
#endif

#define TRIAL_SPAN 70
#define NTRIALS 3

#define CRAM_DEFAULT_LEVEL 5

/* ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
 * Also see the itf8_get and itf8_put macros in cram_io.h
 */

/*
 * LEGACY: consider using itf8_decode_crc.
 *
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val_p) {
    static int nbytes[16] = {
        0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
        1,1,1,1,                                        // 1000xxxx - 1011xxxx
        2,2,                                            // 1100xxxx - 1101xxxx
        3,                                              // 1110xxxx
        4,                                              // 1111xxxx
    };

    static int nbits[16] = {
        0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
        0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
        0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
        0x0f,                                           // 1110xxxx
        0x0f,                                           // 1111xxxx
    };

    int32_t val = hgetc(fd->fp);
    if (val == -1)
        return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
        *val_p = val;
        return 1;

    case 1:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 2;

    case 2:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 3;

    case 3:
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
        return 4;

    case 4: // really 3.5 more, why make it different?
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<4) | (((unsigned char)hgetc(fd->fp)) & 0x0f);
        *val_p = val;
    }

    return 5;
}

int itf8_decode_crc(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    static int nbytes[16] = {
        0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
        1,1,1,1,                                        // 1000xxxx - 1011xxxx
        2,2,                                            // 1100xxxx - 1101xxxx
        3,                                              // 1110xxxx
        4,                                              // 1111xxxx
    };

    static int nbits[16] = {
        0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
        0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
        0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
        0x0f,                                           // 1110xxxx
        0x0f,                                           // 1111xxxx
    };
    unsigned char c[5];

    int32_t val = hgetc(fd->fp);
    if (val == -1)
        return -1;

    c[0]=val;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    if (i > 0) {
        if (hread(fd->fp, &c[1], i) < i)
            return -1;
    }

    switch(i) {
    case 0:
        *val_p = val;
        *crc = crc32(*crc, c, 1);
        return 1;

    case 1:
        val = (val<<8) | c[1];
        *val_p = val;
        *crc = crc32(*crc, c, 2);
        return 2;

    case 2:
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        *val_p = val;
        *crc = crc32(*crc, c, 3);
        return 3;

    case 3:
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        *val_p = val;
        *crc = crc32(*crc, c, 4);
        return 4;

    case 4: // really 3.5 more, why make it different?
        {
            uint32_t uv = val;
            uv = (uv<<8) |  c[1];
            uv = (uv<<8) |  c[2];
            uv = (uv<<8) |  c[3];
            uv = (uv<<4) | (c[4] & 0x0f);
            // Avoid implementation-defined behaviour on negative values
            *val_p = uv < 0x80000000UL ? (int32_t) uv : -((int32_t) (0xffffffffUL - uv)) - 1;
            *crc = crc32(*crc, c, 5);
        }
    }

    return 5;
}

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
static inline int itf8_put(char *cp, int32_t val) {
    unsigned char *up = (unsigned char *)cp;
    if        (!(val & ~0x00000007f)) { // 1 byte
        *up = val;
        return 1;
    } else if (!(val & ~0x00003fff)) { // 2 byte
        *up++ = (val >> 8 ) | 0x80;
        *up   = val & 0xff;
        return 2;
    } else if (!(val & ~0x01fffff)) { // 3 byte
        *up++ = (val >> 16) | 0xc0;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 3;
    } else if (!(val & ~0x0fffffff)) { // 4 byte
        *up++ = (val >> 24) | 0xe0;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 4;
    } else {                           // 5 byte
        *up++ = 0xf0 | ((val>>28) & 0xff);
        *up++ = (val >> 20) & 0xff;
        *up++ = (val >> 12) & 0xff;
        *up++ = (val >> 4 ) & 0xff;
        *up = val & 0x0f;
        return 5;
    }
}


/* 64-bit itf8 variant */
static inline int ltf8_put(char *cp, int64_t val) {
    unsigned char *up = (unsigned char *)cp;
    if        (!(val & ~((1LL<<7)-1))) {
        *up = val;
        return 1;
    } else if (!(val & ~((1LL<<(6+8))-1))) {
        *up++ = (val >> 8 ) | 0x80;
        *up   = val & 0xff;
        return 2;
    } else if (!(val & ~((1LL<<(5+2*8))-1))) {
        *up++ = (val >> 16) | 0xc0;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 3;
    } else if (!(val & ~((1LL<<(4+3*8))-1))) {
        *up++ = (val >> 24) | 0xe0;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 4;
    } else if (!(val & ~((1LL<<(3+4*8))-1))) {
        *up++ = (val >> 32) | 0xf0;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 5;
    } else if (!(val & ~((1LL<<(2+5*8))-1))) {
        *up++ = (val >> 40) | 0xf8;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 6;
    } else if (!(val & ~((1LL<<(1+6*8))-1))) {
        *up++ = (val >> 48) | 0xfc;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 7;
    } else if (!(val & ~((1LL<<(7*8))-1))) {
        *up++ = (val >> 56) | 0xfe;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 8;
    } else {
        *up++ = 0xff;
        *up++ = (val >> 56) & 0xff;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 9;
    }
}

/*
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return hwrite(fd->fp, buf, len) == len ? 0 : -1;
}

const int itf8_bytes[16] = {
    1, 1, 1, 1,  1, 1, 1, 1,
    2, 2, 2, 2,  3, 3, 4, 5
};

const int ltf8_bytes[256] = {
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,

    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,
    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    5, 5, 5, 5,  5, 5, 5, 5,  6, 6, 6, 6,  7, 7, 8, 9
};

/*
 * LEGACY: consider using ltf8_decode_crc.
 */
int ltf8_decode(cram_fd *fd, int64_t *val_p) {
    int c = hgetc(fd->fp);
    int64_t val = (unsigned char)c;
    if (c == -1)
        return -1;

    if (val < 0x80) {
        *val_p =   val;
        return 1;

    } else if (val < 0xc0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & (((1LL<<(6+8)))-1);
        return 2;

    } else if (val < 0xe0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(5+2*8))-1);
        return 3;

    } else if (val < 0xf0) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(4+3*8))-1);
        return 4;

    } else if (val < 0xf8) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(3+4*8))-1);
        return 5;

    } else if (val < 0xfc) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(2+5*8))-1);
        return 6;

    } else if (val < 0xfe) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(1+6*8))-1);
        return 7;

    } else if (val < 0xff) {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val & ((1LL<<(7*8))-1);
        return 8;

    } else {
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        val = (val<<8) | (unsigned char)hgetc(fd->fp);
        *val_p = val;
    }

    return 9;
}

int ltf8_decode_crc(cram_fd *fd, int64_t *val_p, uint32_t *crc) {
    unsigned char c[9];
    int64_t val = hgetc(fd->fp);
    if (val < 0)
        return -1;

    c[0] = val;

    if (val < 0x80) {
        *val_p =   val;
        *crc = crc32(*crc, c, 1);
        return 1;

    } else if (val < 0xc0) {
        int v = hgetc(fd->fp);
        if (v < 0)
            return -1;
        val = (val<<8) | (c[1]=v);
        *val_p = val & (((1LL<<(6+8)))-1);
        *crc = crc32(*crc, c, 2);
        return 2;

    } else if (val < 0xe0) {
        if (hread(fd->fp, &c[1], 2) < 2)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        *val_p = val & ((1LL<<(5+2*8))-1);
        *crc = crc32(*crc, c, 3);
        return 3;

    } else if (val < 0xf0) {
        if (hread(fd->fp, &c[1], 3) < 3)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        *val_p = val & ((1LL<<(4+3*8))-1);
        *crc = crc32(*crc, c, 4);
        return 4;

    } else if (val < 0xf8) {
        if (hread(fd->fp, &c[1], 4) < 4)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        *val_p = val & ((1LL<<(3+4*8))-1);
        *crc = crc32(*crc, c, 5);
        return 5;

    } else if (val < 0xfc) {
        if (hread(fd->fp, &c[1], 5) < 5)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        val = (val<<8) | c[5];
        *val_p = val & ((1LL<<(2+5*8))-1);
        *crc = crc32(*crc, c, 6);
        return 6;

    } else if (val < 0xfe) {
        if (hread(fd->fp, &c[1], 6) < 6)
            return -1;
        val = (val<<8) | c[1];
        val = (val<<8) | c[2];
        val = (val<<8) | c[3];
        val = (val<<8) | c[4];
        val = (val<<8) | c[5];
        val = (val<<8) | c[6];
        *val_p = val & ((1LL<<(1+6*8))-1);
        *crc = crc32(*crc, c, 7);
        return 7;

    } else if (val < 0xff) {
        uint64_t uval = val;
        if (hread(fd->fp, &c[1], 7) < 7)
            return -1;
        uval = (uval<<8) | c[1];
        uval = (uval<<8) | c[2];
        uval = (uval<<8) | c[3];
        uval = (uval<<8) | c[4];
        uval = (uval<<8) | c[5];
        uval = (uval<<8) | c[6];
        uval = (uval<<8) | c[7];
        *val_p = uval & ((1ULL<<(7*8))-1);
        *crc = crc32(*crc, c, 8);
        return 8;

    } else {
        uint64_t uval;
        if (hread(fd->fp, &c[1], 8) < 8)
            return -1;
        uval =             c[1];
        uval = (uval<<8) | c[2];
        uval = (uval<<8) | c[3];
        uval = (uval<<8) | c[4];
        uval = (uval<<8) | c[5];
        uval = (uval<<8) | c[6];
        uval = (uval<<8) | c[7];
        uval = (uval<<8) | c[8];
        *crc = crc32(*crc, c, 9);
        // Avoid implementation-defined behaviour on negative values
        *val_p = c[1] < 0x80 ? (int64_t) uval : -((int64_t) (0xffffffffffffffffULL - uval)) - 1;
    }

    return 9;
}

/*
 * Pushes a value in ITF8 format onto the end of a block.
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int32_t val) {
    char buf[5];
    int sz;

    sz = itf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

int ltf8_put_blk(cram_block *blk, int64_t val) {
    char buf[9];
    int sz;

    sz = ltf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int64_t safe_itf8_get(char **cp, const char *endp, int *err) {
    const unsigned char *up = (unsigned char *)*cp;

    if (endp && endp - *cp < 5 &&
        (*cp >= endp || endp - *cp < itf8_bytes[up[0]>>4])) {
        if (err) *err = 1;
        return 0;
    }

    if (up[0] < 0x80) {
        (*cp)++;
        return up[0];
    } else if (up[0] < 0xc0) {
        (*cp)+=2;
        return ((up[0] <<8) |  up[1])                           & 0x3fff;
    } else if (up[0] < 0xe0) {
        (*cp)+=3;
        return ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
    } else if (up[0] < 0xf0) {
        (*cp)+=4;
        uint32_t uv = (((uint32_t)up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
        return (int32_t)uv;
    } else {
        (*cp)+=5;
        uint32_t uv = (((uint32_t)up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
        return (int32_t)uv;
    }
}

static int64_t safe_ltf8_get(char **cp, const char *endp, int *err) {
    unsigned char *up = (unsigned char *)*cp;

    if (endp && endp - *cp < 9 &&
        (*cp >= endp || endp - *cp < ltf8_bytes[up[0]])) {
        if (err) *err = 1;
        return 0;
    }

    if (up[0] < 0x80) {
        (*cp)++;
        return up[0];
    } else if (up[0] < 0xc0) {
        (*cp)+=2;
        return (((uint64_t)up[0]<< 8) |
                 (uint64_t)up[1]) & (((1LL<<(6+8)))-1);
    } else if (up[0] < 0xe0) {
        (*cp)+=3;
        return (((uint64_t)up[0]<<16) |
                ((uint64_t)up[1]<< 8) |
                 (uint64_t)up[2]) & ((1LL<<(5+2*8))-1);
    } else if (up[0] < 0xf0) {
        (*cp)+=4;
        return (((uint64_t)up[0]<<24) |
                ((uint64_t)up[1]<<16) |
                ((uint64_t)up[2]<< 8) |
                 (uint64_t)up[3]) & ((1LL<<(4+3*8))-1);
    } else if (up[0] < 0xf8) {
        (*cp)+=5;
        return (((uint64_t)up[0]<<32) |
                ((uint64_t)up[1]<<24) |
                ((uint64_t)up[2]<<16) |
                ((uint64_t)up[3]<< 8) |
                 (uint64_t)up[4]) & ((1LL<<(3+4*8))-1);
    } else if (up[0] < 0xfc) {
        (*cp)+=6;
        return (((uint64_t)up[0]<<40) |
                ((uint64_t)up[1]<<32) |
                ((uint64_t)up[2]<<24) |
                ((uint64_t)up[3]<<16) |
                ((uint64_t)up[4]<< 8) |
                 (uint64_t)up[5]) & ((1LL<<(2+5*8))-1);
    } else if (up[0] < 0xfe) {
        (*cp)+=7;
        return (((uint64_t)up[0]<<48) |
                ((uint64_t)up[1]<<40) |
                ((uint64_t)up[2]<<32) |
                ((uint64_t)up[3]<<24) |
                ((uint64_t)up[4]<<16) |
                ((uint64_t)up[5]<< 8) |
                 (uint64_t)up[6]) & ((1LL<<(1+6*8))-1);
    } else if (up[0] < 0xff) {
        (*cp)+=8;
        return (((uint64_t)up[1]<<48) |
                ((uint64_t)up[2]<<40) |
                ((uint64_t)up[3]<<32) |
                ((uint64_t)up[4]<<24) |
                ((uint64_t)up[5]<<16) |
                ((uint64_t)up[6]<< 8) |
                 (uint64_t)up[7]) & ((1LL<<(7*8))-1);
    } else {
        (*cp)+=9;
        return (((uint64_t)up[1]<<56) |
                ((uint64_t)up[2]<<48) |
                ((uint64_t)up[3]<<40) |
                ((uint64_t)up[4]<<32) |
                ((uint64_t)up[5]<<24) |
                ((uint64_t)up[6]<<16) |
                ((uint64_t)up[7]<< 8) |
                 (uint64_t)up[8]);
    }
}

// Wrapper for now
static int safe_itf8_put(char *cp, char *cp_end, int32_t val) {
    return itf8_put(cp, val);
}

static int safe_ltf8_put(char *cp, char *cp_end, int64_t val) {
    return ltf8_put(cp, val);
}

static int itf8_size(int64_t v) {
    return ((!((v)&~0x7f))?1:(!((v)&~0x3fff))?2:(!((v)&~0x1fffff))?3:(!((v)&~0xfffffff))?4:5);
}

//-----------------------------------------------------------------------------

// CRAM v4.0 onwards uses a different variable sized integer encoding
// that is size agnostic.

// Local interface to varint.h inline version, so we can use in func ptr.
// Note a lot of these use the unsigned interface but take signed int64_t.
// This is because the old CRAM ITF8 inteface had signed -1 as unsigned
// 0xffffffff.
static int uint7_size(int64_t v) {
    return var_size_u64(v);
}

static int64_t uint7_get_32(char **cp, const char *endp, int *err) {
    uint32_t val;
    int nb = var_get_u32((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t sint7_get_32(char **cp, const char *endp, int *err) {
    int32_t val;
    int nb = var_get_s32((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t uint7_get_64(char **cp, const char *endp, int *err) {
    uint64_t val;
    int nb = var_get_u64((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int64_t sint7_get_64(char **cp, const char *endp, int *err) {
    int64_t val;
    int nb = var_get_s64((uint8_t *)(*cp), (const uint8_t *)endp, &val);
    (*cp) += nb;
    if (!nb && err) *err = 1;
    return val;
}

static int uint7_put_32(char *cp, char *endp, int32_t val) {
    return var_put_u32((uint8_t *)cp, (uint8_t *)endp, val);
}

static int sint7_put_32(char *cp, char *endp, int32_t val) {
    return var_put_s32((uint8_t *)cp, (uint8_t *)endp, val);
}

static int uint7_put_64(char *cp, char *endp, int64_t val) {
    return var_put_u64((uint8_t *)cp, (uint8_t *)endp, val);
}

static int sint7_put_64(char *cp, char *endp, int64_t val) {
    return var_put_s64((uint8_t *)cp, (uint8_t *)endp, val);
}

// Put direct to to cram_block
static int uint7_put_blk_32(cram_block *blk, int32_t v) {
    uint8_t buf[10];
    int sz = var_put_u32(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int sint7_put_blk_32(cram_block *blk, int32_t v) {
    uint8_t buf[10];
    int sz = var_put_s32(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int uint7_put_blk_64(cram_block *blk, int64_t v) {
    uint8_t buf[10];
    int sz = var_put_u64(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

static int sint7_put_blk_64(cram_block *blk, int64_t v) {
    uint8_t buf[10];
    int sz = var_put_s64(buf, buf+10, v);
    BLOCK_APPEND(blk, buf, sz);
    return sz;

 block_err:
    return -1;
}

// Decode 32-bits with CRC update from cram_fd
static int uint7_decode_crc32(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    uint8_t b[5], i = 0;
    int c;
    uint32_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u32(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 5 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = v;
    return i;
}

// Decode 32-bits with CRC update from cram_fd
static int sint7_decode_crc32(cram_fd *fd, int32_t *val_p, uint32_t *crc) {
    uint8_t b[5], i = 0;
    int c;
    uint32_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u32(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 5 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = (v>>1) ^ -(v&1);
    return i;
}


// Decode 64-bits with CRC update from cram_fd
static int uint7_decode_crc64(cram_fd *fd, int64_t *val_p, uint32_t *crc) {
    uint8_t b[10], i = 0;
    int c;
    uint64_t v = 0;

#ifdef VARINT2
    b[0] = hgetc(fd->fp);
    if (b[0] < 177) {
    } else if (b[0] < 241) {
        b[1] = hgetc(fd->fp);
    } else if (b[0] < 249) {
        b[1] = hgetc(fd->fp);
        b[2] = hgetc(fd->fp);
    } else {
        int n = b[0]+2, z = 1;
        while (n-- >= 249)
            b[z++] = hgetc(fd->fp);
    }
    i = var_get_u64(b, NULL, &v);
#else
//    // Little endian
//    int s = 0;
//    do {
//        b[i++] = c = hgetc(fd->fp);
//        if (c < 0)
//            return -1;
//        v |= (c & 0x7f) << s;
//      s += 7;
//    } while (i < 10 && (c & 0x80));

    // Big endian, see also htscodecs/varint.h
    do {
        b[i++] = c = hgetc(fd->fp);
        if (c < 0)
            return -1;
        v = (v<<7) | (c & 0x7f);
    } while (i < 5 && (c & 0x80));
#endif
    *crc = crc32(*crc, b, i);

    *val_p = v;
    return i;
}

//-----------------------------------------------------------------------------

/*
 * Decodes a 32-bit little endian value from fd and stores in val.
 *
 * Returns the number of bytes read on success
 *         -1 on failure
 */
static int int32_decode(cram_fd *fd, int32_t *val) {
    int32_t i;
    if (4 != hread(fd->fp, &i, 4))
        return -1;

    *val = le_int4(i);
    return 4;
}

/*
 * Encodes a 32-bit little endian value 'val' and writes to fd.
 *
 * Returns the number of bytes written on success
 *         -1 on failure
 */
static int int32_encode(cram_fd *fd, int32_t val) {
    uint32_t v = le_int4(val);
    if (4 != hwrite(fd->fp, &v, 4))
        return -1;

    return 4;
}

/* As int32_decoded/encode, but from/to blocks instead of cram_fd */
int int32_get_blk(cram_block *b, int32_t *val) {
    if (b->uncomp_size - BLOCK_SIZE(b) < 4)
        return -1;

    uint32_t v =
         ((uint32_t) b->data[b->byte  ])        |
        (((uint32_t) b->data[b->byte+1]) <<  8) |
        (((uint32_t) b->data[b->byte+2]) << 16) |
        (((uint32_t) b->data[b->byte+3]) << 24);
    // Avoid implementation-defined behaviour on negative values
    *val = v < 0x80000000U ? (int32_t) v : -((int32_t) (0xffffffffU - v)) - 1;
    BLOCK_SIZE(b) += 4;
    return 4;
}

/* As int32_decoded/encode, but from/to blocks instead of cram_fd */
int int32_put_blk(cram_block *b, int32_t val) {
    unsigned char cp[4];
    uint32_t v = val;
    cp[0] = ( v      & 0xff);
    cp[1] = ((v>>8)  & 0xff);
    cp[2] = ((v>>16) & 0xff);
    cp[3] = ((v>>24) & 0xff);

    BLOCK_APPEND(b, cp, 4);
    return 0;

 block_err:
    return -1;
}

#ifdef HAVE_LIBDEFLATE
/* ----------------------------------------------------------------------
 * libdeflate compression code, with interface to match
 * zlib_mem_{in,de}flate for simplicity elsewhere.
 */

// Named the same as the version that uses zlib as we always use libdeflate for
// decompression when available.
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_decompressor failed");
        return NULL;
    }

    uint8_t *data = NULL, *new_data;
    if (!*size)
        *size = csize*2;
    for(;;) {
        new_data = realloc(data, *size);
        if (!new_data) {
            hts_log_error("Memory allocation failure");
            goto fail;
        }
        data = new_data;

        int ret = libdeflate_gzip_decompress(z, cdata, csize, data, *size, size);

        // Auto grow output buffer size if needed and try again.
        // Fortunately for all bar one call of this we know the size already.
        if (ret == LIBDEFLATE_INSUFFICIENT_SPACE) {
            (*size) *= 1.5;
            continue;
        }

        if (ret != LIBDEFLATE_SUCCESS) {
            hts_log_error("Inflate operation failed: %d", ret);
            goto fail;
        } else {
            break;
        }
    }

    libdeflate_free_decompressor(z);
    return (char *)data;

 fail:
    libdeflate_free_decompressor(z);
    free(data);
    return NULL;
}

// Named differently as we use both zlib/libdeflate for compression.
static char *libdeflate_deflate(char *data, size_t size, size_t *cdata_size,
                                int level, int strat) {
    level = level > 0 ? level : 6; // libdeflate doesn't honour -1 as default
    level *= 1.23;     // NB levels go up to 12 here; 5 onwards is +1
    level += level>=8; // 5,6,7->6,7,8  8->10  9->12
    if (level > 12) level = 12;

    if (strat == Z_RLE) // not supported by libdeflate
        level = 1;

    struct libdeflate_compressor *z = libdeflate_alloc_compressor(level);
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_compressor failed");
        return NULL;
    }

    unsigned char *cdata = NULL; /* Compressed output */
    size_t cdata_alloc;
    cdata = malloc(cdata_alloc = size*1.05+100);
    if (!cdata) {
        hts_log_error("Memory allocation failure");
        libdeflate_free_compressor(z);
        return NULL;
    }

    *cdata_size = libdeflate_gzip_compress(z, data, size, cdata, cdata_alloc);
    libdeflate_free_compressor(z);

    if (*cdata_size == 0) {
        hts_log_error("Call to libdeflate_gzip_compress failed");
        free(cdata);
        return NULL;
    }

    return (char *)cdata;
}

#else

/* ----------------------------------------------------------------------
 * zlib compression code - from Gap5's tg_iface_g.c
 * They're static here as they're only used within the cram_compress_block
 * and cram_uncompress_block functions, which are the external interface.
 */
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, and scale after that */
    data = malloc(data_alloc = csize*1.2+100);
    if (!data)
        return NULL;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    //err = inflateInit(&s);
    err = inflateInit2(&s, 15 + 32);
    if (err != Z_OK) {
        hts_log_error("Call to zlib inflateInit failed: %s", s.msg);
        free(data);
        return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
        unsigned char *data_tmp;
        int alloc_inc;

        s.next_out = &data[s.total_out];
        err = inflate(&s, Z_NO_FLUSH);
        if (err == Z_STREAM_END)
            break;

        if (err != Z_OK) {
            hts_log_error("Call to zlib inflate failed: %s", s.msg);
            free(data);
            inflateEnd(&s);
            return NULL;
        }

        /* More to come, so realloc based on growth so far */
        alloc_inc = (double)s.avail_in/s.total_in * s.total_out + 100;
        data = realloc((data_tmp = data), data_alloc += alloc_inc);
        if (!data) {
            free(data_tmp);
            inflateEnd(&s);
            return NULL;
        }
        s.avail_out += alloc_inc;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return (char *)data;
}
#endif

#if !defined(HAVE_LIBDEFLATE) || LIBDEFLATE_VERSION_MAJOR < 1 || (LIBDEFLATE_VERSION_MAJOR ==  1 && LIBDEFLATE_VERSION_MINOR <= 8)
static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size,
                              int level, int strat) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

    cdata = malloc(cdata_alloc = size*1.05+100);
    if (!cdata)
        return NULL;
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.avail_in = size;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, level, Z_DEFLATED, 15|16, 9, strat);
    if (err != Z_OK) {
        hts_log_error("Call to zlib deflateInit2 failed: %s", s.msg);
        return NULL;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
        s.next_out = &cdata[cdata_pos];
        s.avail_out = cdata_alloc - cdata_pos;
        if (cdata_alloc - cdata_pos <= 0) {
            hts_log_error("Deflate produced larger output than expected");
            return NULL;
        }
        err = deflate(&s, Z_NO_FLUSH);
        cdata_pos = cdata_alloc - s.avail_out;
        if (err != Z_OK) {
            hts_log_error("Call to zlib deflate failed: %s", s.msg);
            break;
        }
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
        hts_log_error("Call to zlib deflate failed: %s", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
        hts_log_error("Call to zlib deflate failed: %s", s.msg);
    }
    return (char *)cdata;
}
#endif

#ifdef HAVE_LIBLZMA
/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using liblzma (xz)
 *
 * On a test set this shrunk the main db from 136157104 bytes to 114796168, but
 * caused tg_index to grow from 2m43.707s to 15m3.961s. Exporting as bfastq
 * went from 18.3s to 36.3s. So decompression suffers too, but not as bad
 * as compression times.
 *
 * For now we disable this functionality. If it's to be reenabled make sure you
 * improve the mem_inflate implementation as it's just a test hack at the
 * moment.
 */

static char *lzma_mem_deflate(char *data, size_t size, size_t *cdata_size,
                              int level) {
    char *out;
    size_t out_size = lzma_stream_buffer_bound(size);
    *cdata_size = 0;

    out = malloc(out_size);

    /* Single call compression */
    if (LZMA_OK != lzma_easy_buffer_encode(level, LZMA_CHECK_CRC32, NULL,
                                           (uint8_t *)data, size,
                                           (uint8_t *)out, cdata_size,
                                           out_size))
        return NULL;

    return out;
}

static char *lzma_mem_inflate(char *cdata, size_t csize, size_t *size) {
    lzma_stream strm = LZMA_STREAM_INIT;
    size_t out_size = 0, out_pos = 0;
    char *out = NULL, *new_out;
    int r;

    /* Initiate the decoder */
    if (LZMA_OK != lzma_stream_decoder(&strm, lzma_easy_decoder_memusage(9), 0))
        return NULL;

    /* Decode loop */
    strm.avail_in = csize;
    strm.next_in = (uint8_t *)cdata;

    for (;strm.avail_in;) {
        if (strm.avail_in > out_size - out_pos) {
            out_size += strm.avail_in * 4 + 32768;
            new_out = realloc(out, out_size);
            if (!new_out)
                goto fail;
            out = new_out;
        }
        strm.avail_out = out_size - out_pos;
        strm.next_out = (uint8_t *)&out[out_pos];

        r = lzma_code(&strm, LZMA_RUN);
        if (LZMA_OK != r && LZMA_STREAM_END != r) {
            hts_log_error("LZMA decode failure (error %d)", r);
            goto fail;
        }

        out_pos = strm.total_out;

        if (r == LZMA_STREAM_END)
            break;
    }

    /* finish up any unflushed data; necessary? */
    r = lzma_code(&strm, LZMA_FINISH);
    if (r != LZMA_OK && r != LZMA_STREAM_END) {
        hts_log_error("Call to lzma_code failed with error %d", r);
        goto fail;
    }

    new_out = realloc(out, strm.total_out > 0 ? strm.total_out : 1);
    if (new_out)
        out = new_out;
    *size = strm.total_out;

    lzma_end(&strm);

    return out;

 fail:
    lzma_end(&strm);
    free(out);
    return NULL;
}
#endif

/* ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*
 * Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * Returns block pointer on success
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
                           int content_id) {
    cram_block *b = malloc(sizeof(*b));
    if (!b)
        return NULL;
    b->method = b->orig_method = RAW;
    b->content_type = content_type;
    b->content_id = content_id;
    b->comp_size = 0;
    b->uncomp_size = 0;
    b->data = NULL;
    b->alloc = 0;
    b->byte = 0;
    b->bit = 7; // MSB
    b->crc32 = 0;
    b->idx = 0;
    b->m = NULL;

    return b;
}

/*
 * Reads a block from a cram file.
 * Returns cram_block pointer on success.
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd) {
    cram_block *b = malloc(sizeof(*b));
    unsigned char c;
    uint32_t crc = 0;
    if (!b)
        return NULL;

    //fprintf(stderr, "Block at %d\n", (int)ftell(fd->fp));

    if (-1 == (b->method      = hgetc(fd->fp))) { free(b); return NULL; }
    c = b->method; crc = crc32(crc, &c, 1);
    if (-1 == (b->content_type= hgetc(fd->fp))) { free(b); return NULL; }
    c = b->content_type; crc = crc32(crc, &c, 1);
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->content_id, &crc))  { free(b); return NULL; }
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->comp_size, &crc))   { free(b); return NULL; }
    if (-1 == fd->vv.varint_decode32_crc(fd, &b->uncomp_size, &crc)) { free(b); return NULL; }

    //fprintf(stderr, "  method %d, ctype %d, cid %d, csize %d, ucsize %d\n",
    //      b->method, b->content_type, b->content_id, b->comp_size, b->uncomp_size);

    if (b->method == RAW) {
        if (b->uncomp_size < 0 || b->comp_size != b->uncomp_size) {
            free(b);
            return NULL;
        }
        b->alloc = b->uncomp_size;
        if (!(b->data = malloc(b->uncomp_size))){ free(b); return NULL; }
        if (b->uncomp_size != hread(fd->fp, b->data, b->uncomp_size)) {
            free(b->data);
            free(b);
            return NULL;
        }
    } else {
        if (b->comp_size < 0 || b->uncomp_size < 0) {
            free(b);
            return NULL;
        }
        b->alloc = b->comp_size;
        if (!(b->data = malloc(b->comp_size)))  { free(b); return NULL; }
        if (b->comp_size != hread(fd->fp, b->data, b->comp_size)) {
            free(b->data);
            free(b);
            return NULL;
        }
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        if (-1 == int32_decode(fd, (int32_t *)&b->crc32)) {
            free(b->data);
            free(b);
            return NULL;
        }

        b->crc32_checked = fd->ignore_md5;
        b->crc_part = crc;
    } else {
        b->crc32_checked = 1; // CRC not present
    }

    b->orig_method = b->method;
    b->idx = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}


/*
 * Computes the size of a cram block, including the block
 * header itself.
 */
uint32_t cram_block_size(cram_block *b) {
    unsigned char dat[100], *cp = dat;;
    uint32_t sz;

    *cp++ = b->method;
    *cp++ = b->content_type;
    cp += itf8_put((char*)cp, b->content_id);
    cp += itf8_put((char*)cp, b->comp_size);
    cp += itf8_put((char*)cp, b->uncomp_size);

    sz = cp-dat + 4;
    sz += b->method == RAW ? b->uncomp_size : b->comp_size;

    return sz;
}

/*
 * Writes a CRAM block.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b) {
    char vardata[100];
    int vardata_o = 0;

    assert(b->method != RAW || (b->comp_size == b->uncomp_size));

    if (hputc(b->method,       fd->fp)  == EOF) return -1;
    if (hputc(b->content_type, fd->fp)  == EOF) return -1;
    vardata_o += fd->vv.varint_put32(vardata          , vardata+100, b->content_id);
    vardata_o += fd->vv.varint_put32(vardata+vardata_o, vardata+100, b->comp_size);
    vardata_o += fd->vv.varint_put32(vardata+vardata_o, vardata+100, b->uncomp_size);
    if (vardata_o != hwrite(fd->fp, vardata, vardata_o))
        return -1;

    if (b->data) {
        if (b->method == RAW) {
            if (b->uncomp_size != hwrite(fd->fp, b->data, b->uncomp_size))
                return -1;
        } else {
            if (b->comp_size != hwrite(fd->fp, b->data, b->comp_size))
                return -1;
        }
    } else {
        // Absent blocks should be size 0
        assert(b->method == RAW && b->uncomp_size == 0);
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        char dat[100], *cp = (char *)dat;
        uint32_t crc;

        *cp++ = b->method;
        *cp++ = b->content_type;
        cp += fd->vv.varint_put32(cp, dat+100, b->content_id);
        cp += fd->vv.varint_put32(cp, dat+100, b->comp_size);
        cp += fd->vv.varint_put32(cp, dat+100, b->uncomp_size);
        crc = crc32(0L, (uc *)dat, cp-dat);

        if (b->method == RAW) {
            b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->uncomp_size);
        } else {
            b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->comp_size);
        }

        if (-1 == int32_encode(fd, b->crc32))
            return -1;
    }

    return 0;
}

/*
 * Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b) {
    if (!b)
        return;
    if (b->data)
        free(b->data);
    free(b);
}

/*
 * Uncompresses a CRAM block, if compressed.
 */
int cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    if (b->crc32_checked == 0) {
        uint32_t crc = crc32(b->crc_part, b->data ? b->data : (uc *)"", b->alloc);
        b->crc32_checked = 1;
        if (crc != b->crc32) {
            hts_log_error("Block CRC32 failure");
            return -1;
        }
    }

    if (b->uncomp_size == 0) {
        // blank block
        b->method = RAW;
        return 0;
    }
    assert(b->uncomp_size >= 0); // cram_read_block should ensure this

    switch (b->method) {
    case RAW:
        return 0;

    case GZIP:
        uncomp_size = b->uncomp_size;
        uncomp = zlib_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);

        if (!uncomp)
            return -1;
        if (uncomp_size != b->uncomp_size) {
            free(uncomp);
            return -1;
        }
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = uncomp_size;
        b->method = RAW;
        break;

#ifdef HAVE_LIBBZ2
    case BZIP2: {
        unsigned int usize = b->uncomp_size;
        if (!(uncomp = malloc(usize)))
            return -1;
        if (BZ_OK != BZ2_bzBuffToBuffDecompress(uncomp, &usize,
                                                (char *)b->data, b->comp_size,
                                                0, 0)) {
            free(uncomp);
            return -1;
        }
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = usize;
        b->method = RAW;
        b->uncomp_size = usize; // Just in case it differs
        break;
    }
#else
    case BZIP2:
        hts_log_error("Bzip2 compression is not compiled into this version. Please rebuild and try again");
        return -1;
#endif

#ifdef HAVE_LIBLZMA
    case LZMA:
        uncomp = lzma_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);
        if (!uncomp)
            return -1;
        if (uncomp_size != b->uncomp_size) {
            free(uncomp);
            return -1;
        }
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = uncomp_size;
        b->method = RAW;
        break;
#else
    case LZMA:
        hts_log_error("Lzma compression is not compiled into this version. Please rebuild and try again");
        return -1;
        break;
#endif

    case RANS: {
        unsigned int usize = b->uncomp_size, usize2;
        uncomp = (char *)rans_uncompress(b->data, b->comp_size, &usize2);
        if (!uncomp)
            return -1;
        if (usize != usize2) {
            free(uncomp);
            return -1;
        }
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = usize2;
        b->method = RAW;
        b->uncomp_size = usize2; // Just in case it differs
        //fprintf(stderr, "Expanded %d to %d\n", b->comp_size, b->uncomp_size);
        break;
    }

    case FQZ: {
        uncomp_size = b->uncomp_size;
        uncomp = fqz_decompress((char *)b->data, b->comp_size, &uncomp_size, NULL, 0);
        if (!uncomp)
            return -1;
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = uncomp_size;
        b->method = RAW;
        b->uncomp_size = uncomp_size;
        break;
    }

    case RANS_PR0: {
        unsigned int usize = b->uncomp_size, usize2;
        uncomp = (char *)rans_uncompress_4x16(b->data, b->comp_size, &usize2);
        if (!uncomp)
            return -1;
        if (usize != usize2) {
            free(uncomp);
            return -1;
        }
        b->orig_method = RANS_PR0 + (b->data[0]&1)
            + 2*((b->data[0]&0x40)>0) + 4*((b->data[0]&0x80)>0);
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = usize2;
        b->method = RAW;
        b->uncomp_size = usize2; // Just incase it differs
        //fprintf(stderr, "Expanded %d to %d\n", b->comp_size, b->uncomp_size);
        break;
    }

    case ARITH_PR0: {
        unsigned int usize = b->uncomp_size, usize2;
        uncomp = (char *)arith_uncompress_to(b->data, b->comp_size, NULL, &usize2);
        if (!uncomp)
            return -1;
        if (usize != usize2) {
            free(uncomp);
            return -1;
        }
        b->orig_method = ARITH_PR0 + (b->data[0]&1)
            + 2*((b->data[0]&0x40)>0) + 4*((b->data[0]&0x80)>0);
        free(b->data);
        b->data = (unsigned char *)uncomp;
        b->alloc = usize2;
        b->method = RAW;
        b->uncomp_size = usize2; // Just incase it differs
        //fprintf(stderr, "Expanded %d to %d\n", b->comp_size, b->uncomp_size);
        break;
    }

    case TOK3: {
        uint32_t out_len;
        uint8_t *cp = tok3_decode_names(b->data, b->comp_size, &out_len);
        if (!cp)
            return -1;
        b->orig_method = TOK3;
        b->method = RAW;
        free(b->data);
        b->data = cp;
        b->alloc = out_len;
        b->uncomp_size = out_len;
        break;
    }

    default:
        return -1;
    }

    return 0;
}

static char *cram_compress_by_method(cram_slice *s, char *in, size_t in_size,
                                     int content_id, size_t *out_size,
                                     enum cram_block_method_int method,
                                     int level, int strat) {
    switch (method) {
    case GZIP:
    case GZIP_RLE:
    case GZIP_1:
        // Read names bizarrely benefit from zlib over libdeflate for
        // mid-range compression levels.  Focusing purely of ratio or
        // speed, libdeflate still wins.  It also seems to win for
        // other data series too.
        //
        // Eg RN at level 5;  libdeflate=55.9MB  zlib=51.6MB
#ifdef HAVE_LIBDEFLATE
#  if (LIBDEFLATE_VERSION_MAJOR < 1 || (LIBDEFLATE_VERSION_MAJOR == 1 && LIBDEFLATE_VERSION_MINOR <= 8))
        if (content_id == DS_RN && level >= 4 && level <= 7)
            return zlib_mem_deflate(in, in_size, out_size, level, strat);
        else
#  endif
            return libdeflate_deflate(in, in_size, out_size, level, strat);
#else
        return zlib_mem_deflate(in, in_size, out_size, level, strat);
#endif

    case BZIP2: {
#ifdef HAVE_LIBBZ2
        unsigned int comp_size = in_size*1.01 + 600;
        char *comp = malloc(comp_size);
        if (!comp)
            return NULL;

        if (BZ_OK != BZ2_bzBuffToBuffCompress(comp, &comp_size,
                                              in, in_size,
                                              level, 0, 30)) {
            free(comp);
            return NULL;
        }
        *out_size = comp_size;
        return comp;
#else
        return NULL;
#endif
    }

    case FQZ:
    case FQZ_b:
    case FQZ_c:
    case FQZ_d: {
        // Extract the necessary portion of the slice into an fqz_slice struct.
        // These previously were the same thing, but this permits us to detach
        // the codec from the rest of this CRAM implementation.
        fqz_slice *f = malloc(2*s->hdr->num_records * sizeof(uint32_t) + sizeof(fqz_slice));
        if (!f)
            return NULL;
        f->num_records = s->hdr->num_records;
        f->len = (uint32_t *)(((char *)f) + sizeof(fqz_slice));
        f->flags = f->len + s->hdr->num_records;
        int i;
        for (i = 0; i < s->hdr->num_records; i++) {
            f->flags[i] = s->crecs[i].flags;
            f->len[i] = (i+1 < s->hdr->num_records
                         ? s->crecs[i+1].qual - s->crecs[i].qual
                         : s->block[DS_QS]->uncomp_size - s->crecs[i].qual);
        }
        char *comp = fqz_compress(strat & 0xff /* cram vers */, f,
                                  in, in_size, out_size, strat >> 8, NULL);
        free(f);
        return comp;
    }

    case LZMA:
#ifdef HAVE_LIBLZMA
        return lzma_mem_deflate(in, in_size, out_size, level);
#else
        return NULL;
#endif

    case RANS0:
    case RANS1: {
        unsigned int out_size_i;
        unsigned char *cp;
        cp = rans_compress((unsigned char *)in, in_size, &out_size_i,
                           method == RANS0 ? 0 : 1);
        *out_size = out_size_i;
        return (char *)cp;
    }

    case RANS_PR0:
    case RANS_PR1:
    case RANS_PR64:
    case RANS_PR9:
    case RANS_PR128:
    case RANS_PR129:
    case RANS_PR192:
    case RANS_PR193: {
        unsigned int out_size_i;
        unsigned char *cp;

        // see enum cram_block. We map RANS_* methods to order bit-fields
        static int methmap[] = { 1, 64,9, 128,129, 192,193 };

        cp = rans_compress_4x16((unsigned char *)in, in_size, &out_size_i,
                                method == RANS_PR0 ? 0 : methmap[method - RANS_PR1]);
        *out_size = out_size_i;
        return (char *)cp;
    }

    case ARITH_PR0:
    case ARITH_PR1:
    case ARITH_PR64:
    case ARITH_PR9:
    case ARITH_PR128:
    case ARITH_PR129:
    case ARITH_PR192:
    case ARITH_PR193: {
        unsigned int out_size_i;
        unsigned char *cp;

        // see enum cram_block. We map ARITH_* methods to order bit-fields
        static int methmap[] = { 1, 64,9, 128,129, 192,193 };

        cp = arith_compress_to((unsigned char *)in, in_size, NULL, &out_size_i,
                               method == ARITH_PR0 ? 0 : methmap[method - ARITH_PR1]);
        *out_size = out_size_i;
        return (char *)cp;
    }

    case TOK3:
    case TOKA: {
        int out_len;
        int lev = level;
        if (method == TOK3 && lev > 3)
            lev = 3;
        uint8_t *cp = tok3_encode_names(in, in_size, lev, strat, &out_len, NULL);
        *out_size = out_len;
        return (char *)cp;
    }

    case RAW:
        break;

    default:
        return NULL;
    }

    return NULL;
}


/*
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 *
 * Method and level -1 implies defaults, as specified in cram_fd.
 */
int cram_compress_block2(cram_fd *fd, cram_slice *s,
                         cram_block *b, cram_metrics *metrics,
                         int method, int level) {

    if (!b)
        return 0;

    char *comp = NULL;
    size_t comp_size = 0;
    int strat;

    // Internally we have parameterised methods that externally map
    // to the same CRAM method value.
    // See enum_cram_block_method_int in cram_structs.h.
    int methmap[] = {
        // Externally defined values
        RAW, GZIP, BZIP2, LZMA, RANS, RANSPR, ARITH, FQZ, TOK3,

        // Reserved for possible expansion
        0, 0,

        // Internally parameterised versions matching back to above
        // external values
        GZIP, GZIP,
        FQZ, FQZ, FQZ,
        RANS,
        RANSPR, RANSPR, RANSPR, RANSPR, RANSPR, RANSPR, RANSPR,
        TOK3,
        ARITH,  ARITH,  ARITH,  ARITH,  ARITH,  ARITH,  ARITH,
    };

    if (b->method != RAW) {
        // Maybe already compressed if s->block[0] was compressed and
        // we have e.g. s->block[DS_BA] set to s->block[0] due to only
        // one base type present and hence using E_HUFFMAN on block 0.
        // A second explicit attempt to compress the same block then
        // occurs.
        return 0;
    }

    if (method == -1) {
        method = 1<<GZIP;
        if (fd->use_bz2)
            method |= 1<<BZIP2;
        if (fd->use_lzma)
            method |= 1<<LZMA;
    }

    if (level == -1)
        level = fd->level;

    //fprintf(stderr, "IN: block %d, sz %d\n", b->content_id, b->uncomp_size);

    if (method == RAW || level == 0 || b->uncomp_size == 0) {
        b->method = RAW;
        b->comp_size = b->uncomp_size;
        //fprintf(stderr, "Skip block id %d\n", b->content_id);
        return 0;
    }

#ifndef ABS
#    define ABS(a) ((a)>=0?(a):-(a))
#endif

    if (metrics) {
        pthread_mutex_lock(&fd->metrics_lock);
        // Sudden changes in size trigger a retrial.  These are mainly
        // triggered when switching to sorted / unsorted, where the number
        // of elements in a slice radically changes.
        //
        // We also get large fluctuations based on genome coordinate for
        // e.g. SA:Z and SC series, but we consider the typical scale of
        // delta between blocks and use this to look for abnormality.
        if (metrics->input_avg_sz &&
            (b->uncomp_size + 1000 > 4*(metrics->input_avg_sz+1000) ||
             b->uncomp_size + 1000 < (metrics->input_avg_sz+1000)/4) &&
            ABS(b->uncomp_size-metrics->input_avg_sz)
                > 10*metrics->input_avg_delta) {
            metrics->next_trial = 0;
        }

        if (metrics->trial > 0 || --metrics->next_trial <= 0) {
            int m, unpackable = metrics->unpackable;
            size_t sz_best = b->uncomp_size;
            size_t sz[CRAM_MAX_METHOD] = {0};
            int method_best = 0; // RAW
            char *c_best = NULL, *c = NULL;

            metrics->input_avg_delta =
                0.9 * (metrics->input_avg_delta +
                       ABS(b->uncomp_size - metrics->input_avg_sz));

            metrics->input_avg_sz += b->uncomp_size*.2;
            metrics->input_avg_sz *= 0.8;

            if (metrics->revised_method)
                method = metrics->revised_method;
            else
                metrics->revised_method = method;

            if (metrics->next_trial <= 0) {
                metrics->next_trial = TRIAL_SPAN;
                metrics->trial = NTRIALS;
                for (m = 0; m < CRAM_MAX_METHOD; m++)
                    metrics->sz[m] /= 2;
                metrics->unpackable = 0;
            }

            // Compress this block using the best method
            if (unpackable && CRAM_MAJOR_VERS(fd->version) > 3) {
                // No point trying bit-pack if 17+ symbols.
                if (method & (1<<RANS_PR128))
                    method = (method|(1<<RANS_PR0))&~(1<<RANS_PR128);
                if (method & (1<<RANS_PR129))
                    method = (method|(1<<RANS_PR1))&~(1<<RANS_PR129);
                if (method & (1<<RANS_PR192))
                    method = (method|(1<<RANS_PR64))&~(1<<RANS_PR192);
                if (method & (1<<RANS_PR193))
                    method = (method|(1<<RANS_PR64)|(1<<RANS_PR1))&~(1<<RANS_PR193);

                if (method & (1<<ARITH_PR128))
                    method = (method|(1<<ARITH_PR0))&~(1<<ARITH_PR128);
                if (method & (1<<ARITH_PR129))
                    method = (method|(1<<ARITH_PR1))&~(1<<ARITH_PR129);
                if (method & (1<<ARITH_PR192))
                    method = (method|(1<<ARITH_PR64))&~(1<<ARITH_PR192);
                if (method & (1u<<ARITH_PR193))
                    method = (method|(1<<ARITH_PR64)|(1<<ARITH_PR1))&~(1u<<ARITH_PR193);
            }

            // Libdeflate doesn't have a Z_RLE strategy.
            // We treat it as level 1, but iff we haven't also
            // explicitly listed that in the method list.
#ifdef HAVE_LIBDEFLATE
            if ((method & (1<<GZIP_RLE)) && (method & (1<<GZIP_1)))
                method &= ~(1<<GZIP_RLE);
#endif

            pthread_mutex_unlock(&fd->metrics_lock);

            for (m = 0; m < CRAM_MAX_METHOD; m++) {
                if (method & (1u<<m)) {
                    int lvl = level;
                    switch (m) {
                    case GZIP:     strat = Z_FILTERED; break;
                    case GZIP_1:   strat = Z_DEFAULT_STRATEGY; lvl = 1; break;
                    case GZIP_RLE: strat = Z_RLE; break;
                    case FQZ:      strat = CRAM_MAJOR_VERS(fd->version); break;
                    case FQZ_b:    strat = CRAM_MAJOR_VERS(fd->version)+256; break;
                    case FQZ_c:    strat = CRAM_MAJOR_VERS(fd->version)+2*256; break;
                    case FQZ_d:    strat = CRAM_MAJOR_VERS(fd->version)+3*256; break;
                    case TOK3:     strat = 0; break;
                    case TOKA:     strat = 1; break;
                    default:       strat = 0;
                    }

                    c = cram_compress_by_method(s, (char *)b->data, b->uncomp_size,
                                                b->content_id, &sz[m], m, lvl, strat);

                    if (c && sz_best > sz[m]) {
                        sz_best = sz[m];
                        method_best = m;
                        if (c_best)
                            free(c_best);
                        c_best = c;
                    } else if (c) {
                        free(c);
                    } else {
                        sz[m] = b->uncomp_size*2+1000; // arbitrarily worse than raw
                    }
                } else {
                    sz[m] = b->uncomp_size*2+1000; // arbitrarily worse than raw
                }
            }

            if (c_best) {
                free(b->data);
                b->data = (unsigned char *)c_best;
                b->method = method_best; // adjusted to methmap[method_best] later
                b->comp_size = sz_best;
            }

            // Accumulate stats for all methods tried
            pthread_mutex_lock(&fd->metrics_lock);
            for (m = 0; m < CRAM_MAX_METHOD; m++)
                // don't be overly sure on small blocks.
                // +2000 means eg bzip2 vs gzip (1.07 to 1.04) or gz vs rans1
                // needs to be at least 60 bytes smaller to overcome the
                // fixed size addition.
                metrics->sz[m] += sz[m]+2000;

            // When enough trials performed, find the best on average
            if (--metrics->trial == 0) {
                int best_method = RAW;
                int best_sz = INT_MAX;

                // Relative costs of methods. See enum_cram_block_method_int
                // and methmap
                double meth_cost[32] = {
                    // Externally defined methods
                    1,    // 0  raw
                    1.04, // 1  gzip (Z_FILTERED)
                    1.07, // 2  bzip2
                    1.08, // 3  lzma
                    1.00, // 4  rans    (O0)
                    1.00, // 5  ranspr  (O0)
                    1.04, // 6  arithpr (O0)
                    1.05, // 7  fqz
                    1.05, // 8  tok3 (rans)
                    1.00, 1.00, // 9,10 reserved

                    // Paramterised versions of above
                    1.01, // gzip rle
                    1.01, // gzip -1

                    1.05, 1.05, 1.05, // FQZ_b,c,d

                    1.01, // rans O1

                    1.01, // rans_pr1
                    1.00, // rans_pr64; if smaller, usually fast
                    1.03, // rans_pr65/9
                    1.00, // rans_pr128
                    1.01, // rans_pr129
                    1.00, // rans_pr192
                    1.01, // rans_pr193

                    1.07, // tok3 arith

                    1.04, // arith_pr1
                    1.04, // arith_pr64
                    1.04, // arith_pr9
                    1.03, // arith_pr128
                    1.04, // arith_pr129
                    1.04, // arith_pr192
                    1.04, // arith_pr193
                };

                // Scale methods by cost based on compression level
                if (fd->level <= 1) {
                    for (m = 0; m < CRAM_MAX_METHOD; m++)
                        metrics->sz[m] *= 1+(meth_cost[m]-1)*4;
                } else if (fd->level <= 3) {
                    for (m = 0; m < CRAM_MAX_METHOD; m++)
                        metrics->sz[m] *= 1+(meth_cost[m]-1);
                } else if (fd->level <= 6) {
                    for (m = 0; m < CRAM_MAX_METHOD; m++)
                        metrics->sz[m] *= 1+(meth_cost[m]-1)/2;
                } else if (fd->level <= 7) {
                    for (m = 0; m < CRAM_MAX_METHOD; m++)
                        metrics->sz[m] *= 1+(meth_cost[m]-1)/3;
                } // else cost is ignored

                // Ensure these are never used; BSC and ZSTD
                metrics->sz[9] = metrics->sz[10] = INT_MAX;

                for (m = 0; m < CRAM_MAX_METHOD; m++) {
                    if ((!metrics->sz[m]) || (!(method & (1u<<m))))
                        continue;

                    if (best_sz > metrics->sz[m])
                        best_sz = metrics->sz[m], best_method = m;
                }

                if (best_method != metrics->method) {
                    //metrics->trial = (NTRIALS+1)/2; // be sure
                    //metrics->next_trial /= 1.5;
                    metrics->consistency = 0;
                } else {
                    metrics->next_trial *= MIN(2, 1+metrics->consistency/4.0);
                    metrics->consistency++;
                }

                metrics->method = best_method;
                switch (best_method) {
                case GZIP:     strat = Z_FILTERED; break;
                case GZIP_1:   strat = Z_DEFAULT_STRATEGY; break;
                case GZIP_RLE: strat = Z_RLE; break;
                case FQZ:      strat = CRAM_MAJOR_VERS(fd->version); break;
                case FQZ_b:    strat = CRAM_MAJOR_VERS(fd->version)+256; break;
                case FQZ_c:    strat = CRAM_MAJOR_VERS(fd->version)+2*256; break;
                case FQZ_d:    strat = CRAM_MAJOR_VERS(fd->version)+3*256; break;
                case TOK3:     strat = 0; break;
                case TOKA:     strat = 1; break;
                default:       strat = 0;
                }
                metrics->strat  = strat;

                // If we see at least MAXFAIL trials in a row for a specific
                // compression method with more than MAXDELTA aggregate
                // size then we drop this from the list of methods used
                // for this block type.
#define MAXDELTA 0.20
#define MAXFAILS 4
                for (m = 0; m < CRAM_MAX_METHOD; m++) {
                    if (best_method == m) {
                        metrics->cnt[m] = 0;
                        metrics->extra[m] = 0;
                    } else if (best_sz < metrics->sz[m]) {
                        double r = (double)metrics->sz[m] / best_sz - 1;
                        int mul = 1+(fd->level>=7);
                        if (++metrics->cnt[m] >= MAXFAILS*mul &&
                            (metrics->extra[m] += r) >= MAXDELTA*mul)
                            method &= ~(1u<<m);

                        // Special case for fqzcomp as it rarely changes
                        if (m == FQZ || m == FQZ_b || m == FQZ_c || m == FQZ_d) {
                            if (metrics->sz[m] > best_sz)
                                method &= ~(1u<<m);
                        }
                    }
                }

                //if (fd->verbose > 1 && method != metrics->revised_method)
                //    fprintf(stderr, "%d: revising method from %x to %x\n",
                //          b->content_id, metrics->revised_method, method);
                metrics->revised_method = method;
            }
            pthread_mutex_unlock(&fd->metrics_lock);
        } else {
            metrics->input_avg_delta =
                0.9 * (metrics->input_avg_delta +
                       ABS(b->uncomp_size - metrics->input_avg_sz));

            metrics->input_avg_sz += b->uncomp_size*.2;
            metrics->input_avg_sz *= 0.8;

            strat = metrics->strat;
            method = metrics->method;

            pthread_mutex_unlock(&fd->metrics_lock);
            comp = cram_compress_by_method(s, (char *)b->data, b->uncomp_size,
                                           b->content_id, &comp_size, method,
                                           method == GZIP_1 ? 1 : level,
                                           strat);
            if (!comp)
                return -1;

            if (comp_size < b->uncomp_size) {
                free(b->data);
                b->data = (unsigned char *)comp;
                b->comp_size = comp_size;
                b->method = method;
            } else {
                free(comp);
            }
        }

    } else {
        // no cached metrics, so just do zlib?
        comp = cram_compress_by_method(s, (char *)b->data, b->uncomp_size,
                                       b->content_id, &comp_size, GZIP, level, Z_FILTERED);
        if (!comp) {
            hts_log_error("Compression failed!");
            return -1;
        }

        if (comp_size < b->uncomp_size) {
            free(b->data);
            b->data = (unsigned char *)comp;
            b->comp_size = comp_size;
            b->method = GZIP;
        } else {
            free(comp);
        }
        strat = Z_FILTERED;
    }

    hts_log_info("Compressed block ID %d from %d to %d by method %s",
                 b->content_id, b->uncomp_size, b->comp_size,
                 cram_block_method2str(b->method));

    b->method = methmap[b->method];

    return 0;
}
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
                        int method, int level) {
    return cram_compress_block2(fd, NULL, b, metrics, method, level);
}

cram_metrics *cram_new_metrics(void) {
    cram_metrics *m = calloc(1, sizeof(*m));
    if (!m)
        return NULL;
    m->trial = NTRIALS-1;
    m->next_trial = TRIAL_SPAN/2; // learn quicker at start
    m->method = RAW;
    m->strat = 0;
    m->revised_method = 0;
    m->unpackable = 0;

    return m;
}

char *cram_block_method2str(enum cram_block_method_int m) {
    switch(m) {
    case RAW:         return "RAW";
    case GZIP:        return "GZIP";
    case BZIP2:       return "BZIP2";
    case LZMA:        return "LZMA";
    case RANS0:       return "RANS0";
    case RANS1:       return "RANS1";
    case GZIP_RLE:    return "GZIP_RLE";
    case GZIP_1:      return "GZIP_1";
    case FQZ:         return "FQZ";
    case FQZ_b:       return "FQZ_b";
    case FQZ_c:       return "FQZ_c";
    case FQZ_d:       return "FQZ_d";
    case RANS_PR0:    return "RANS_PR0";
    case RANS_PR1:    return "RANS_PR1";
    case RANS_PR64:   return "RANS_PR64";
    case RANS_PR9:    return "RANS_PR9";
    case RANS_PR128:  return "RANS_PR128";
    case RANS_PR129:  return "RANS_PR129";
    case RANS_PR192:  return "RANS_PR192";
    case RANS_PR193:  return "RANS_PR193";
    case TOK3:        return "TOK3_R";
    case TOKA:        return "TOK3_A";
    case ARITH_PR0:   return "ARITH_PR0";
    case ARITH_PR1:   return "ARITH_PR1";
    case ARITH_PR64:  return "ARITH_PR64";
    case ARITH_PR9:   return "ARITH_PR9";
    case ARITH_PR128: return "ARITH_PR128";
    case ARITH_PR129: return "ARITH_PR129";
    case ARITH_PR192: return "ARITH_PR192";
    case ARITH_PR193: return "ARITH_PR193";
    case BM_ERROR: break;
    }
    return "?";
}

char *cram_content_type2str(enum cram_content_type t) {
    switch (t) {
    case FILE_HEADER:         return "FILE_HEADER";
    case COMPRESSION_HEADER:  return "COMPRESSION_HEADER";
    case MAPPED_SLICE:        return "MAPPED_SLICE";
    case UNMAPPED_SLICE:      return "UNMAPPED_SLICE";
    case EXTERNAL:            return "EXTERNAL";
    case CORE:                return "CORE";
    case CT_ERROR:            break;
    }
    return "?";
}

/* ----------------------------------------------------------------------
 * Reference sequence handling
 *
 * These revolve around the refs_t structure, which may potentially be
 * shared between multiple cram_fd.
 *
 * We start with refs_create() to allocate an empty refs_t and then
 * populate it with @SQ line data using refs_from_header(). This is done on
 * cram_open().  Also at start up we can call cram_load_reference() which
 * is used with "scramble -r foo.fa". This replaces the fd->refs with the
 * new one specified. In either case refs2id() is then called which
 * maps ref_entry names to @SQ ids (refs_t->ref_id[]).
 *
 * Later, possibly within a thread, we will want to know the actual ref
 * seq itself, obtained by calling cram_get_ref().  This may use the
 * UR: or M5: fields or the filename specified in the original
 * cram_load_reference() call.
 *
 * Given the potential for multi-threaded reference usage, we have
 * reference counting (sorry for the confusing double use of "ref") to
 * track the number of callers interested in any specific reference.
 */

/*
 * Frees/unmaps a reference sequence and associated file handles.
 */
static void ref_entry_free_seq(ref_entry *e) {
    if (e->mf)
        mfclose(e->mf);
    if (e->seq && !e->mf)
        free(e->seq);

    e->seq = NULL;
    e->mf = NULL;
}

void refs_free(refs_t *r) {
    RP("refs_free()\n");

    if (--r->count > 0)
        return;

    if (!r)
        return;

    if (r->pool)
        string_pool_destroy(r->pool);

    if (r->h_meta) {
        khint_t k;

        for (k = kh_begin(r->h_meta); k != kh_end(r->h_meta); k++) {
            ref_entry *e;

            if (!kh_exist(r->h_meta, k))
                continue;
            if (!(e = kh_val(r->h_meta, k)))
                continue;
            ref_entry_free_seq(e);
            free(e);
        }

        kh_destroy(refs, r->h_meta);
    }

    if (r->ref_id)
        free(r->ref_id);

    if (r->fp)
        bgzf_close(r->fp);

    pthread_mutex_destroy(&r->lock);

    free(r);
}

static refs_t *refs_create(void) {
    refs_t *r = calloc(1, sizeof(*r));

    RP("refs_create()\n");

    if (!r)
        return NULL;

    if (!(r->pool = string_pool_create(8192)))
        goto err;

    r->ref_id = NULL; // see refs2id() to populate.
    r->count = 1;
    r->last = NULL;
    r->last_id = -1;

    if (!(r->h_meta = kh_init(refs)))
        goto err;

    pthread_mutex_init(&r->lock, NULL);

    return r;

 err:
    refs_free(r);
    return NULL;
}

/*
 * Opens a reference fasta file as a BGZF stream, allowing for
 * compressed files.  It automatically builds a .fai file if
 * required and if compressed a .gzi bgzf index too.
 *
 * Returns a BGZF handle on success;
 *         NULL on failure.
 */
static BGZF *bgzf_open_ref(char *fn, char *mode, int is_md5) {
    BGZF *fp;

    if (!is_md5 && !hisremote(fn)) {
        char fai_file[PATH_MAX];

        snprintf(fai_file, PATH_MAX, "%s.fai", fn);
        if (access(fai_file, R_OK) != 0)
            if (fai_build(fn) != 0)
                return NULL;
    }

    if (!(fp = bgzf_open(fn, mode))) {
        perror(fn);
        return NULL;
    }

    if (fp->is_compressed == 1 && bgzf_index_load(fp, fn, ".gzi") < 0) {
        hts_log_error("Unable to load .gzi index '%s.gzi'", fn);
        bgzf_close(fp);
        return NULL;
    }

    return fp;
}

/*
 * Loads a FAI file for a reference.fasta.
 * "is_err" indicates whether failure to load is worthy of emitting an
 * error message. In some cases (eg with embedded references) we
 * speculatively load, just in case, and silently ignore errors.
 *
 * Returns the refs_t struct on success (maybe newly allocated);
 *         NULL on failure
 */
static refs_t *refs_load_fai(refs_t *r_orig, const char *fn, int is_err) {
    hFILE *fp = NULL;
    char fai_fn[PATH_MAX];
    char line[8192];
    refs_t *r = r_orig;
    size_t fn_l = strlen(fn);
    int id = 0, id_alloc = 0;

    RP("refs_load_fai %s\n", fn);

    if (!r)
        if (!(r = refs_create()))
            goto err;

    if (r->fp)
        if (bgzf_close(r->fp) != 0)
            goto err;
    r->fp = NULL;

    /* Look for a FASTA##idx##FAI format */
    char *fn_delim = strstr(fn, HTS_IDX_DELIM);
    if (fn_delim) {
        if (!(r->fn = string_ndup(r->pool, fn, fn_delim - fn)))
            goto err;
        fn_delim += strlen(HTS_IDX_DELIM);
        snprintf(fai_fn, PATH_MAX, "%s", fn_delim);
    } else {
        /* An index file was provided, instead of the actual reference file */
        if (fn_l > 4 && strcmp(&fn[fn_l-4], ".fai") == 0) {
            if (!r->fn) {
                if (!(r->fn = string_ndup(r->pool, fn, fn_l-4)))
                    goto err;
            }
            snprintf(fai_fn, PATH_MAX, "%s", fn);
        } else {
        /* Only the reference file provided. Get the index file name from it */
            if (!(r->fn = string_dup(r->pool, fn)))
                goto err;
            sprintf(fai_fn, "%.*s.fai", PATH_MAX-5, fn);
        }
    }

    if (!(r->fp = bgzf_open_ref(r->fn, "r", 0))) {
        hts_log_error("Failed to open reference file '%s'", r->fn);
        goto err;
    }

    if (!(fp = hopen(fai_fn, "r"))) {
        hts_log_error("Failed to open index file '%s'", fai_fn);
        if (is_err)
            perror(fai_fn);
        goto err;
    }
    while (hgets(line, 8192, fp) != NULL) {
        ref_entry *e = malloc(sizeof(*e));
        char *cp;
        int n;
        khint_t k;

        if (!e)
            return NULL;

        // id
        for (cp = line; *cp && !isspace_c(*cp); cp++)
            ;
        *cp++ = 0;
        e->name = string_dup(r->pool, line);

        // length
        while (*cp && isspace_c(*cp))
            cp++;
        e->length = strtoll(cp, &cp, 10);

        // offset
        while (*cp && isspace_c(*cp))
            cp++;
        e->offset = strtoll(cp, &cp, 10);

        // bases per line
        while (*cp && isspace_c(*cp))
            cp++;
        e->bases_per_line = strtol(cp, &cp, 10);

        // line length
        while (*cp && isspace_c(*cp))
            cp++;
        e->line_length = strtol(cp, &cp, 10);

        // filename
        e->fn = r->fn;

        e->count = 0;
        e->seq = NULL;
        e->mf = NULL;
        e->is_md5 = 0;
        e->validated_md5 = 0;

        k = kh_put(refs, r->h_meta, e->name, &n);
        if (-1 == n)  {
            free(e);
            return NULL;
        }

        if (n) {
            kh_val(r->h_meta, k) = e;
        } else {
            ref_entry *re = kh_val(r->h_meta, k);
            if (re && (re->count != 0 || re->length != 0)) {
                /* Keep old */
                free(e);
            } else {
                /* Replace old */
                if (re)
                    free(re);
                kh_val(r->h_meta, k) = e;
            }
        }

        if (id >= id_alloc) {
            ref_entry **new_refs;
            int x;

            id_alloc = id_alloc ?id_alloc*2 : 16;
            new_refs = realloc(r->ref_id, id_alloc * sizeof(*r->ref_id));
            if (!new_refs)
                goto err;
            r->ref_id = new_refs;

            for (x = id; x < id_alloc; x++)
                r->ref_id[x] = NULL;
        }
        r->ref_id[id] = e;
        r->nref = ++id;
    }

    if(hclose(fp) < 0)
        goto err;
    return r;

 err:
    if (fp)
        hclose_abruptly(fp);

    if (!r_orig)
        refs_free(r);

    return NULL;
}

/*
 * Verifies that the CRAM @SQ lines and .fai files match.
 */
static void sanitise_SQ_lines(cram_fd *fd) {
    int i;

    if (!fd->header || !fd->header->hrecs)
        return;

    if (!fd->refs || !fd->refs->h_meta)
        return;

    for (i = 0; i < fd->header->hrecs->nref; i++) {
        const char *name = fd->header->hrecs->ref[i].name;
        khint_t k = kh_get(refs, fd->refs->h_meta, name);
        ref_entry *r;

        // We may have @SQ lines which have no known .fai, but do not
        // in themselves pose a problem because they are unused in the file.
        if (k == kh_end(fd->refs->h_meta))
            continue;

        if (!(r = (ref_entry *)kh_val(fd->refs->h_meta, k)))
            continue;

        if (r->length && r->length != fd->header->hrecs->ref[i].len) {
            assert(strcmp(r->name, fd->header->hrecs->ref[i].name) == 0);

            // Should we also check MD5sums here to ensure the correct
            // reference was given?
            hts_log_warning("Header @SQ length mismatch for ref %s, %"PRIhts_pos" vs %d",
                            r->name, fd->header->hrecs->ref[i].len, (int)r->length);

            // Fixing the parsed @SQ header will make MD:Z: strings work
            // and also stop it producing N for the sequence.
            fd->header->hrecs->ref[i].len = r->length;
        }
    }
}

/*
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int refs2id(refs_t *r, sam_hdr_t *hdr) {
    int i;
    sam_hrecs_t *h = hdr->hrecs;

    if (r->ref_id)
        free(r->ref_id);
    if (r->last)
        r->last = NULL;

    r->ref_id = calloc(h->nref, sizeof(*r->ref_id));
    if (!r->ref_id)
        return -1;

    r->nref = h->nref;
    for (i = 0; i < h->nref; i++) {
        khint_t k = kh_get(refs, r->h_meta, h->ref[i].name);
        if (k != kh_end(r->h_meta)) {
            r->ref_id[i] = kh_val(r->h_meta, k);
        } else {
            hts_log_warning("Unable to find ref name '%s'", h->ref[i].name);
        }
    }

    return 0;
}

/*
 * Generates refs_t entries based on @SQ lines in the header.
 * Returns 0 on success
 *         -1 on failure
 */
static int refs_from_header(cram_fd *fd) {
    if (!fd)
        return -1;

    refs_t *r = fd->refs;
    if (!r)
        return -1;

    sam_hdr_t *h = fd->header;
    if (!h)
        return 0;

    if (!h->hrecs) {
        if (-1 == sam_hdr_fill_hrecs(h))
            return -1;
    }

    if (h->hrecs->nref == 0)
        return 0;

    //fprintf(stderr, "refs_from_header for %p mode %c\n", fd, fd->mode);

    /* Existing refs are fine, as long as they're compatible with the hdr. */
    ref_entry **new_ref_id = realloc(r->ref_id, (r->nref + h->hrecs->nref) * sizeof(*r->ref_id));
    if (!new_ref_id)
        return -1;
    r->ref_id = new_ref_id;

    int i, j;
    /* Copy info from h->ref[i] over to r */
    for (i = 0, j = r->nref; i < h->hrecs->nref; i++) {
        sam_hrec_type_t *ty;
        sam_hrec_tag_t *tag;
        khint_t k;
        int n;

        k = kh_get(refs, r->h_meta, h->hrecs->ref[i].name);
        if (k != kh_end(r->h_meta))
            // Ref already known about
            continue;

        if (!(r->ref_id[j] = calloc(1, sizeof(ref_entry))))
            return -1;

        if (!h->hrecs->ref[i].name)
            return -1;

        r->ref_id[j]->name = string_dup(r->pool, h->hrecs->ref[i].name);
        if (!r->ref_id[j]->name) return -1;
        r->ref_id[j]->length = 0; // marker for not yet loaded

        /* Initialise likely filename if known */
        if ((ty = sam_hrecs_find_type_id(h->hrecs, "SQ", "SN", h->hrecs->ref[i].name))) {
            if ((tag = sam_hrecs_find_key(ty, "M5", NULL))) {
                r->ref_id[j]->fn = string_dup(r->pool, tag->str+3);
                //fprintf(stderr, "Tagging @SQ %s / %s\n", r->ref_id[h]->name, r->ref_id[h]->fn);
            }
        }

        k = kh_put(refs, r->h_meta, r->ref_id[j]->name, &n);
        if (n <= 0) // already exists or error
            return -1;
        kh_val(r->h_meta, k) = r->ref_id[j];

        j++;
    }
    r->nref = j;

    return 0;
}

/*
 * Attaches a header to a cram_fd.
 *
 * This should be used when creating a new cram_fd for writing where
 * we have a header already constructed (eg from a file we've read
 * in).
 */
int cram_set_header2(cram_fd *fd, const sam_hdr_t *hdr) {
    if (!fd || !hdr )
        return -1;

    if (fd->header != hdr) {
        if (fd->header)
            sam_hdr_destroy(fd->header);
        fd->header = sam_hdr_dup(hdr);
        if (!fd->header)
            return -1;
    }
    return refs_from_header(fd);
}

int cram_set_header(cram_fd *fd, sam_hdr_t *hdr) {
    return cram_set_header2(fd, hdr);
}

/*
 * Returns whether the path refers to a directory.
 */
static int is_directory(char *fn) {
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

/*
 * Converts a directory and a filename into an expanded path, replacing %s
 * in directory with the filename and %[0-9]+s with portions of the filename
 * Any remaining parts of filename are added to the end with /%s.
 */
static int expand_cache_path(char *path, char *dir, const char *fn) {
    char *cp, *start = path;
    size_t len;
    size_t sz = PATH_MAX;

    while ((cp = strchr(dir, '%'))) {
        if (cp-dir >= sz) return -1;
        strncpy(path, dir, cp-dir);
        path += cp-dir;
        sz -= cp-dir;

        if (*++cp == 's') {
            len = strlen(fn);
            if (len >= sz) return -1;
            strcpy(path, fn);
            path += len;
            sz -= len;
            fn += len;
            cp++;
        } else if (*cp >= '0' && *cp <= '9') {
            char *endp;
            long l;

            l = strtol(cp, &endp, 10);
            l = MIN(l, strlen(fn));
            if (*endp == 's') {
                if (l >= sz) return -1;
                strncpy(path, fn, l);
                path += l;
                fn += l;
                sz -= l;
                *path = 0;
                cp = endp+1;
            } else {
                if (sz < 3) return -1;
                *path++ = '%';
                *path++ = *cp++;
            }
        } else {
            if (sz < 3) return -1;
            *path++ = '%';
            *path++ = *cp++;
        }
        dir = cp;
    }

    len = strlen(dir);
    if (len >= sz) return -1;
    strcpy(path, dir);
    path += len;
    sz -= len;

    len = strlen(fn) + ((*fn && path > start && path[-1] != '/') ? 1 : 0);
    if (len >= sz) return -1;
    if (*fn && path > start && path[-1] != '/')
        *path++ = '/';
    strcpy(path, fn);
    return 0;
}

/*
 * Make the directory containing path and any prefix directories.
 */
static void mkdir_prefix(char *path, int mode) {
    char *cp = strrchr(path, '/');
    if (!cp)
        return;

    *cp = 0;
    if (is_directory(path)) {
        *cp = '/';
        return;
    }

    if (mkdir(path, mode) == 0) {
        chmod(path, mode);
        *cp = '/';
        return;
    }

    mkdir_prefix(path, mode);
    mkdir(path, mode);
    chmod(path, mode);
    *cp = '/';
}

/*
 * Return the cache directory to use, based on the first of these
 * environment variables to be set to a non-empty value.
 */
static const char *get_cache_basedir(const char **extra) {
    char *base;

    *extra = "";

    base = getenv("XDG_CACHE_HOME");
    if (base && *base) return base;

    base = getenv("HOME");
    if (base && *base) { *extra = "/.cache"; return base; }

    base = getenv("TMPDIR");
    if (base && *base) return base;

    base = getenv("TEMP");
    if (base && *base) return base;

    return "/tmp";
}

/*
 * Queries the M5 string from the header and attempts to populate the
 * reference from this using the REF_PATH environment.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_populate_ref(cram_fd *fd, int id, ref_entry *r) {
    char *ref_path = getenv("REF_PATH");
    sam_hrec_type_t *ty;
    sam_hrec_tag_t *tag;
    char path[PATH_MAX];
    kstring_t path_tmp = KS_INITIALIZE;
    char cache[PATH_MAX], cache_root[PATH_MAX];
    char *local_cache = getenv("REF_CACHE");
    mFILE *mf;
    int local_path = 0;

    hts_log_info("Running cram_populate_ref on fd %p, id %d", (void *)fd, id);

    cache_root[0] = '\0';

    if (!ref_path || *ref_path == '\0') {
        /*
         * If we have no ref path, we use the EBI server.
         * However to avoid spamming it we require a local ref cache too.
         */
        ref_path = "https://www.ebi.ac.uk/ena/cram/md5/%s";
        if (!local_cache || *local_cache == '\0') {
            const char *extra;
            const char *base = get_cache_basedir(&extra);
            snprintf(cache_root, PATH_MAX, "%s%s/hts-ref", base, extra);
            snprintf(cache,PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
            local_cache = cache;
            hts_log_info("Populating local cache: %s", local_cache);
        }
    }

    if (!r->name)
        return -1;

    if (!(ty = sam_hrecs_find_type_id(fd->header->hrecs, "SQ", "SN", r->name)))
        return -1;

    if (!(tag = sam_hrecs_find_key(ty, "M5", NULL)))
        goto no_M5;

    hts_log_info("Querying ref %s", tag->str+3);

    /* Use cache if available */
    if (local_cache && *local_cache) {
        if (expand_cache_path(path, local_cache, tag->str+3) == 0)
            local_path = 1;
    }

#ifndef HAVE_MMAP
    char *path2;
    /* Search local files in REF_PATH; we can open them and return as above */
    if (!local_path && (path2 = find_path(tag->str+3, ref_path))) {
        int len = snprintf(path, PATH_MAX, "%s", path2);
        free(path2);
        if (len > 0 && len < PATH_MAX) // in case it's too long
            local_path = 1;
    }
#endif

    /* Found via REF_CACHE or local REF_PATH file */
    if (local_path) {
        struct stat sb;
        BGZF *fp;

        if (0 == stat(path, &sb)
            && S_ISREG(sb.st_mode)
            && (fp = bgzf_open(path, "r"))) {
            r->length = sb.st_size;
            r->offset = r->line_length = r->bases_per_line = 0;

            r->fn = string_dup(fd->refs->pool, path);

            if (fd->refs->fp)
                if (bgzf_close(fd->refs->fp) != 0)
                    return -1;
            fd->refs->fp = fp;
            fd->refs->fn = r->fn;
            r->is_md5 = 1;
            r->validated_md5 = 1;

            // Fall back to cram_get_ref() where it'll do the actual
            // reading of the file.
            return 0;
        }
    }


    /* Otherwise search full REF_PATH; slower as loads entire file */
    if ((mf = open_path_mfile(tag->str+3, ref_path, NULL))) {
        size_t sz;
        r->seq = mfsteal(mf, &sz);
        if (r->seq) {
            r->mf = NULL;
        } else {
            // keep mf around as we couldn't detach
            r->seq = mf->data;
            r->mf = mf;
        }
        r->length = sz;
        r->is_md5 = 1;
        r->validated_md5 = 1;
    } else {
        refs_t *refs;
        const char *fn;

    no_M5:
        /* Failed to find in search path or M5 cache, see if @SQ UR: tag? */
        if (!(tag = sam_hrecs_find_key(ty, "UR", NULL)))
            return -1;

        fn = (strncmp(tag->str+3, "file:", 5) == 0)
            ? tag->str+8
            : tag->str+3;

        if (fd->refs->fp) {
            if (bgzf_close(fd->refs->fp) != 0)
                return -1;
            fd->refs->fp = NULL;
        }
        if (!(refs = refs_load_fai(fd->refs, fn, 0)))
            return -1;
        sanitise_SQ_lines(fd);

        fd->refs = refs;
        if (fd->refs->fp) {
            if (bgzf_close(fd->refs->fp) != 0)
                return -1;
            fd->refs->fp = NULL;
        }

        if (!fd->refs->fn)
            return -1;

        if (-1 == refs2id(fd->refs, fd->header))
            return -1;
        if (!fd->refs->ref_id || !fd->refs->ref_id[id])
            return -1;

        // Local copy already, so fall back to cram_get_ref().
        return 0;
    }

    /* Populate the local disk cache if required */
    if (local_cache && *local_cache) {
        hFILE *fp;

        if (*cache_root && !is_directory(cache_root)) {
            hts_log_warning("Creating reference cache directory %s\n"
                            "This may become large; see the samtools(1) manual page REF_CACHE discussion",
                            cache_root);
        }

        if (expand_cache_path(path, local_cache, tag->str+3) < 0) {
            return 0; // Not fatal - we have the data already so keep going.
        }
        hts_log_info("Writing cache file '%s'", path);
        mkdir_prefix(path, 01777);

        fp = hts_open_tmpfile(path, "wx", &path_tmp);
        if (!fp) {
            perror(path_tmp.s);
            free(path_tmp.s);

            // Not fatal - we have the data already so keep going.
            return 0;
        }

        // Check md5sum
        hts_md5_context *md5;
        char unsigned md5_buf1[16];
        char md5_buf2[33];

        if (!(md5 = hts_md5_init())) {
            hclose_abruptly(fp);
            unlink(path_tmp.s);
            free(path_tmp.s);
            return -1;
        }
        hts_md5_update(md5, r->seq, r->length);
        hts_md5_final(md5_buf1, md5);
        hts_md5_destroy(md5);
        hts_md5_hex(md5_buf2, md5_buf1);

        if (strncmp(tag->str+3, md5_buf2, 32) != 0) {
            hts_log_error("Mismatching md5sum for downloaded reference");
            hclose_abruptly(fp);
            unlink(path_tmp.s);
            free(path_tmp.s);
            return -1;
        }

        ssize_t length_written = hwrite(fp, r->seq, r->length);
        if (hclose(fp) < 0 || length_written != r->length ||
            chmod(path_tmp.s, 0444) < 0 ||
            rename(path_tmp.s, path) < 0) {
            hts_log_error("Creating reference at %s failed: %s",
                          path, strerror(errno));
            unlink(path_tmp.s);
        }
    }

    free(path_tmp.s);
    return 0;
}

static void cram_ref_incr_locked(refs_t *r, int id) {
    RP("%d INC REF %d, %d %p\n", gettid(), id,
       (int)(id>=0 && r->ref_id[id]?r->ref_id[id]->count+1:-999),
       id>=0 && r->ref_id[id]?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id] || !r->ref_id[id]->seq)
        return;

    if (r->last_id == id)
        r->last_id = -1;

    ++r->ref_id[id]->count;
}

void cram_ref_incr(refs_t *r, int id) {
    pthread_mutex_lock(&r->lock);
    cram_ref_incr_locked(r, id);
    pthread_mutex_unlock(&r->lock);
}

static void cram_ref_decr_locked(refs_t *r, int id) {
    RP("%d DEC REF %d, %d %p\n", gettid(), id,
       (int)(id>=0 && r->ref_id[id]?r->ref_id[id]->count-1:-999),
       id>=0 && r->ref_id[id]?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id] || !r->ref_id[id]->seq) {
        return;
    }

    if (--r->ref_id[id]->count <= 0) {
        assert(r->ref_id[id]->count == 0);
        if (r->last_id >= 0) {
            if (r->ref_id[r->last_id]->count <= 0 &&
                r->ref_id[r->last_id]->seq) {
                RP("%d FREE REF %d (%p)\n", gettid(),
                   r->last_id, r->ref_id[r->last_id]->seq);
                ref_entry_free_seq(r->ref_id[r->last_id]);
                if (r->ref_id[r->last_id]->is_md5) r->ref_id[r->last_id]->length = 0;
            }
        }
        r->last_id = id;
    }
}

void cram_ref_decr(refs_t *r, int id) {
    pthread_mutex_lock(&r->lock);
    cram_ref_decr_locked(r, id);
    pthread_mutex_unlock(&r->lock);
}

/*
 * Used by cram_ref_load and cram_get_ref. The file handle will have
 * already been opened, so we can catch it. The ref_entry *e informs us
 * of whether this is a multi-line fasta file or a raw MD5 style file.
 * Either way we create a single contiguous sequence.
 *
 * Returns all or part of a reference sequence on success (malloced);
 *         NULL on failure.
 */
static char *load_ref_portion(BGZF *fp, ref_entry *e, int start, int end) {
    off_t offset, len;
    char *seq;

    if (end < start)
        end = start;

    /*
     * Compute locations in file. This is trivial for the MD5 files, but
     * is still necessary for the fasta variants.
     *
     * Note the offset here, as with faidx, has the assumption that white-
     * space (the diff between line_length and bases_per_line) only occurs
     * at the end of a line of text.
     */
    offset = e->line_length
        ? e->offset + (start-1)/e->bases_per_line * e->line_length +
          (start-1) % e->bases_per_line
        : start-1;

    len = (e->line_length
           ? e->offset + (end-1)/e->bases_per_line * e->line_length +
             (end-1) % e->bases_per_line
           : end-1) - offset + 1;

    if (bgzf_useek(fp, offset, SEEK_SET) < 0) {
        perror("bgzf_useek() on reference file");
        return NULL;
    }

    if (len == 0 || !(seq = malloc(len))) {
        return NULL;
    }

    if (len != bgzf_read(fp, seq, len)) {
        perror("bgzf_read() on reference file");
        free(seq);
        return NULL;
    }

    /* Strip white-space if required. */
    if (len != end-start+1) {
        hts_pos_t i, j;
        char *cp = seq;
        char *cp_to;

        // Copy up to the first white-space, and then repeatedly just copy
        // bases_per_line verbatim, and use the slow method to end again.
        //
        // This may seem excessive, but this code can be a significant
        // portion of total CRAM decode CPU time for shallow data sets.
        for (i = j = 0; i < len; i++) {
            if (!isspace_c(cp[i]))
                cp[j++] = cp[i] & ~0x20;
            else
                break;
        }
        while (i < len && isspace_c(cp[i]))
            i++;
        while (i < len - e->line_length) {
            hts_pos_t j_end = j + e->bases_per_line;
            while (j < j_end)
                cp[j++] = cp[i++] & ~0x20; // toupper equiv
            i += e->line_length - e->bases_per_line;
        }
        for (; i < len; i++) {
            if (!isspace_c(cp[i]))
                cp[j++] = cp[i] & ~0x20;
        }

        cp_to = cp+j;

        if (cp_to - seq != end-start+1) {
            hts_log_error("Malformed reference file");
            free(seq);
            return NULL;
        }
    } else {
        int i;
        for (i = 0; i < len; i++) {
            seq[i] = toupper_c(seq[i]);
        }
    }

    return seq;
}

/*
 * Load the entire reference 'id'.
 * This also increments the reference count by 1.
 *
 * Returns ref_entry on success;
 *         NULL on failure
 */
ref_entry *cram_ref_load(refs_t *r, int id, int is_md5) {
    ref_entry *e = r->ref_id[id];
    int start = 1, end = e->length;
    char *seq;

    if (e->seq) {
        return e;
    }

    assert(e->count == 0);

    if (r->last) {
#ifdef REF_DEBUG
        int idx = 0;
        for (idx = 0; idx < r->nref; idx++)
            if (r->last == r->ref_id[idx])
                break;
        RP("%d cram_ref_load DECR %d\n", gettid(), idx);
#endif
        assert(r->last->count > 0);
        if (--r->last->count <= 0) {
            RP("%d FREE REF %d (%p)\n", gettid(), id, r->ref_id[id]->seq);
            if (r->last->seq)
                ref_entry_free_seq(r->last);
        }
    }

    if (!r->fn)
        return NULL;

    /* Open file if it's not already the current open reference */
    if (strcmp(r->fn, e->fn) || r->fp == NULL) {
        if (r->fp)
            if (bgzf_close(r->fp) != 0)
                return NULL;
        r->fn = e->fn;
        if (!(r->fp = bgzf_open_ref(r->fn, "r", is_md5)))
            return NULL;
    }

    RP("%d Loading ref %d (%d..%d)\n", gettid(), id, start, end);

    if (!(seq = load_ref_portion(r->fp, e, start, end))) {
        return NULL;
    }

    RP("%d Loaded ref %d (%d..%d) = %p\n", gettid(), id, start, end, seq);

    RP("%d INC REF %d, %"PRId64"\n", gettid(), id, (e->count+1));
    e->seq = seq;
    e->mf = NULL;
    e->count++;

    /*
     * Also keep track of last used ref so incr/decr loops on the same
     * sequence don't cause load/free loops.
     */
    RP("%d cram_ref_load INCR %d => %"PRId64"\n", gettid(), id, e->count+1);
    r->last = e;
    e->count++;

    return e;
}

/*
 * Returns a portion of a reference sequence from start to end inclusive.
 * The returned pointer is owned by either the cram_file fd or by the
 * internal refs_t structure and should not be freed  by the caller.
 *
 * The difference is whether or not this refs_t is in use by just the one
 * cram_fd or by multiples, or whether we have multiple threads accessing
 * references. In either case fd->shared will be true and we start using
 * reference counting to track the number of users of a specific reference
 * sequence.
 *
 * Otherwise the ref seq returned is allocated as part of cram_fd itself
 * and will be freed up on the next call to cram_get_ref or cram_close.
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * To cease using a reference, call cram_ref_decr().
 *
 * Returns reference on success,
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, int start, int end) {
    ref_entry *r;
    char *seq;
    int ostart = start;

    if (id == -1 || start < 1)
        return NULL;

    /* FIXME: axiomatic query of r->seq being true?
     * Or shortcut for unsorted data where we load once and never free?
     */

    //fd->shared_ref = 1; // hard code for now to simplify things

    pthread_mutex_lock(&fd->ref_lock);

    RP("%d cram_get_ref on fd %p, id %d, range %d..%d\n", gettid(), fd, id, start, end);

    /*
     * Unsorted data implies we want to fetch an entire reference at a time.
     * We just deal with this at the moment by claiming we're sharing
     * references instead, which has the same requirement.
     */
    if (fd->unsorted)
        fd->shared_ref = 1;


    /* Sanity checking: does this ID exist? */
    if (id >= fd->refs->nref) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (!fd->refs || !fd->refs->ref_id[id]) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (!(r = fd->refs->ref_id[id])) {
        hts_log_error("No reference found for id %d", id);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }


    /*
     * It has an entry, but may not have been populated yet.
     * Any manually loaded .fai files have their lengths known.
     * A ref entry computed from @SQ lines (M5 or UR field) will have
     * r->length == 0 unless it's been loaded once and verified that we have
     * an on-disk filename for it.
     *
     * 19 Sep 2013: Moved the lock here as the cram_populate_ref code calls
     * open_path_mfile and libcurl, which isn't multi-thread safe unless I
     * rewrite my code to have one curl handle per thread.
     */
    pthread_mutex_lock(&fd->refs->lock);
    if (r->length == 0) {
        if (fd->ref_fn)
            hts_log_warning("Reference file given, but ref '%s' not present",
                            r->name);
        if (cram_populate_ref(fd, id, r) == -1) {
            hts_log_warning("Failed to populate reference for id %d", id);
            pthread_mutex_unlock(&fd->refs->lock);
            pthread_mutex_unlock(&fd->ref_lock);
            return NULL;
        }
        r = fd->refs->ref_id[id];
        if (fd->unsorted)
            cram_ref_incr_locked(fd->refs, id);
    }


    /*
     * We now know that we the filename containing the reference, so check
     * for limits. If it's over half the reference we'll load all of it in
     * memory as this will speed up subsequent calls.
     */
    if (end < 1)
        end = r->length;
    if (end >= r->length)
        end  = r->length;

    if (end - start >= 0.5*r->length || fd->shared_ref) {
        start = 1;
        end = r->length;
    }

    /*
     * Maybe we have it cached already? If so use it.
     *
     * Alternatively if we don't have the sequence but we're sharing
     * references and/or are asking for the entire length of it, then
     * load the full reference into the refs structure and return
     * a pointer to that one instead.
     */
    if (fd->shared_ref || r->seq || (start == 1 && end == r->length)) {
        char *cp;

        if (id >= 0) {
            if (r->seq) {
                cram_ref_incr_locked(fd->refs, id);
            } else {
                ref_entry *e;
                if (!(e = cram_ref_load(fd->refs, id, r->is_md5))) {
                    pthread_mutex_unlock(&fd->refs->lock);
                    pthread_mutex_unlock(&fd->ref_lock);
                    return NULL;
                }

                /* unsorted data implies cache ref indefinitely, to avoid
                 * continually loading and unloading.
                 */
                if (fd->unsorted)
                    cram_ref_incr_locked(fd->refs, id);
            }

            fd->ref = NULL; /* We never access it directly */
            fd->ref_start = 1;
            fd->ref_end   = r->length;
            fd->ref_id    = id;

            cp = fd->refs->ref_id[id]->seq + ostart-1;
        } else {
            fd->ref = NULL;
            cp = NULL;
        }

        RP("%d cram_get_ref returning for id %d, count %d\n", gettid(), id, (int)r->count);

        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return cp;
    }

    /*
     * Otherwise we're not sharing, we don't have a copy of it already and
     * we're only asking for a small portion of it.
     *
     * In this case load up just that segment ourselves, freeing any old
     * small segments in the process.
     */

    /* Unmapped ref ID */
    if (id < 0 || !fd->refs->fn) {
        if (fd->ref_free) {
            free(fd->ref_free);
            fd->ref_free = NULL;
        }
        fd->ref = NULL;
        fd->ref_id = id;
        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    /* Open file if it's not already the current open reference */
    if (strcmp(fd->refs->fn, r->fn) || fd->refs->fp == NULL) {
        if (fd->refs->fp)
            if (bgzf_close(fd->refs->fp) != 0)
                return NULL;
        fd->refs->fn = r->fn;
        if (!(fd->refs->fp = bgzf_open_ref(fd->refs->fn, "r", r->is_md5))) {
            pthread_mutex_unlock(&fd->refs->lock);
            pthread_mutex_unlock(&fd->ref_lock);
            return NULL;
        }
    }

    if (!(fd->ref = load_ref_portion(fd->refs->fp, r, start, end))) {
        pthread_mutex_unlock(&fd->refs->lock);
        pthread_mutex_unlock(&fd->ref_lock);
        return NULL;
    }

    if (fd->ref_free)
        free(fd->ref_free);

    fd->ref_id    = id;
    fd->ref_start = start;
    fd->ref_end   = end;
    fd->ref_free = fd->ref;
    seq = fd->ref;

    pthread_mutex_unlock(&fd->refs->lock);
    pthread_mutex_unlock(&fd->ref_lock);

    return seq ? seq + ostart - start : NULL;
}

/*
 * If fd has been opened for reading, it may be permitted to specify 'fn'
 * as NULL and let the code auto-detect the reference by parsing the
 * SAM header @SQ lines.
 */
int cram_load_reference(cram_fd *fd, char *fn) {
    int ret = 0;

    if (fn) {
        fd->refs = refs_load_fai(fd->refs, fn,
                                 !(fd->embed_ref>0 && fd->mode == 'r'));
        fn = fd->refs ? fd->refs->fn : NULL;
        if (!fn)
            ret = -1;
        sanitise_SQ_lines(fd);
    }
    fd->ref_fn = fn;

    if ((!fd->refs || (fd->refs->nref == 0 && !fn)) && fd->header) {
        if (fd->refs)
            refs_free(fd->refs);
        if (!(fd->refs = refs_create()))
            return -1;
        if (-1 == refs_from_header(fd))
            return -1;
    }

    if (fd->header)
        if (-1 == refs2id(fd->refs, fd->header))
            return -1;

    return ret;
}

/* ----------------------------------------------------------------------
 * Containers
 */

/*
 * Creates a new container, specifying the maximum number of slices
 * and records permitted.
 *
 * Returns cram_container ptr on success
 *         NULL on failure
 */
cram_container *cram_new_container(int nrec, int nslice) {
    cram_container *c = calloc(1, sizeof(*c));
    enum cram_DS_ID id;

    if (!c)
        return NULL;

    c->curr_ref = -2;

    c->max_c_rec = nrec * nslice;
    c->curr_c_rec = 0;

    c->max_rec = nrec;
    c->record_counter = 0;
    c->num_bases = 0;
    c->s_num_bases = 0;

    c->max_slice = nslice;
    c->curr_slice = 0;

    c->pos_sorted = 1;
    c->max_apos   = 0;
    c->multi_seq  = 0;
    c->qs_seq_orient = 1;
    c->no_ref = 0;
    c->embed_ref = -1; // automatic selection

    c->bams = NULL;

    if (!(c->slices = calloc(nslice != 0 ? nslice : 1, sizeof(cram_slice *))))
        goto err;
    c->slice = NULL;

    if (!(c->comp_hdr = cram_new_compression_header()))
        goto err;
    c->comp_hdr_block = NULL;

    for (id = DS_RN; id < DS_TN; id++)
        if (!(c->stats[id] = cram_stats_create())) goto err;

    //c->aux_B_stats = cram_stats_create();

    if (!(c->tags_used = kh_init(m_tagmap)))
        goto err;
    c->refs_used = 0;
    c->ref_free = 0;

    return c;

 err:
    if (c) {
        if (c->slices)
            free(c->slices);
        free(c);
    }
    return NULL;
}

void cram_free_container(cram_container *c) {
    enum cram_DS_ID id;
    int i;

    if (!c)
        return;

    if (c->refs_used)
        free(c->refs_used);

    if (c->landmark)
        free(c->landmark);

    if (c->comp_hdr)
        cram_free_compression_header(c->comp_hdr);

    if (c->comp_hdr_block)
        cram_free_block(c->comp_hdr_block);

    // Free the slices; filled out by encoder only
    if (c->slices) {
        for (i = 0; i < c->max_slice; i++) {
            if (c->slices[i])
                cram_free_slice(c->slices[i]);
            if (c->slices[i] == c->slice)
                c->slice = NULL;
        }
        free(c->slices);
    }

    // Free the current slice; set by both encoder & decoder
    if (c->slice) {
        cram_free_slice(c->slice);
        c->slice = NULL;
    }

    for (id = DS_RN; id < DS_TN; id++)
        if (c->stats[id]) cram_stats_free(c->stats[id]);

    //if (c->aux_B_stats) cram_stats_free(c->aux_B_stats);

    if (c->tags_used) {
        khint_t k;

        for (k = kh_begin(c->tags_used); k != kh_end(c->tags_used); k++) {
            if (!kh_exist(c->tags_used, k))
                continue;

            cram_tag_map *tm = (cram_tag_map *)kh_val(c->tags_used, k);
            if (tm) {
                cram_codec *c = tm->codec;

                if (c) c->free(c);
                free(tm);
            }
        }

        kh_destroy(m_tagmap, c->tags_used);
    }

    if (c->ref_free)
        free(c->ref);

    free(c);
}

/*
 * Reads a container header.
 *
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd) {
    cram_container c2, *c;
    int i, s;
    size_t rd = 0;
    uint32_t crc = 0;

    fd->err = 0;
    fd->eof = 0;

    memset(&c2, 0, sizeof(c2));
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        if ((s = fd->vv.varint_decode32_crc(fd, &c2.length, &crc)) == -1) {
            fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
    } else if (CRAM_MAJOR_VERS(fd->version) < 4) {
        uint32_t len;
        if ((s = int32_decode(fd, &c2.length)) == -1) {
            if (CRAM_MAJOR_VERS(fd->version) == 2 &&
                CRAM_MINOR_VERS(fd->version) == 0)
                fd->eof = 1; // EOF blocks arrived in v2.1
            else
                fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
        len = le_int4(c2.length);
        crc = crc32(0L, (unsigned char *)&len, 4);
    } else {
        if ((s = fd->vv.varint_decode32_crc(fd, &c2.length, &crc))   == -1) {
            fd->eof = fd->empty_container ? 1 : 2;
            return NULL;
        } else {
            rd+=s;
        }
    }
    if ((s = fd->vv.varint_decode32s_crc(fd, &c2.ref_seq_id, &crc))   == -1) return NULL; else rd+=s;
    if (CRAM_MAJOR_VERS(fd->version) >= 4) {
        int64_t i64;
        if ((s = fd->vv.varint_decode64_crc(fd, &i64, &crc))== -1) return NULL; else rd+=s;
        c2.ref_seq_start = i64;
        if ((s = fd->vv.varint_decode64_crc(fd, &i64, &crc)) == -1) return NULL; else rd+=s;
        c2.ref_seq_span = i64;
    } else {
        int32_t i32;
        if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc))== -1) return NULL; else rd+=s;
        c2.ref_seq_start = i32;
        if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc)) == -1) return NULL; else rd+=s;
        c2.ref_seq_span = i32;
    }
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_records, &crc))  == -1) return NULL; else rd+=s;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        c2.record_counter = 0;
        c2.num_bases = 0;
    } else {
        if (CRAM_MAJOR_VERS(fd->version) >= 3) {
            if ((s = fd->vv.varint_decode64_crc(fd, &c2.record_counter, &crc)) == -1)
                return NULL;
            else
                rd += s;
        } else {
            int32_t i32;
            if ((s = fd->vv.varint_decode32_crc(fd, &i32, &crc)) == -1)
                return NULL;
            else
                rd += s;
            c2.record_counter = i32;
        }

        if ((s = fd->vv.varint_decode64_crc(fd, &c2.num_bases, &crc))== -1)
            return NULL;
        else
            rd += s;
    }
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_blocks, &crc))   == -1)
        return NULL;
    else
        rd+=s;
    if ((s = fd->vv.varint_decode32_crc(fd, &c2.num_landmarks, &crc))== -1)
        return NULL;
    else
        rd+=s;

    if (c2.num_landmarks < 0 || c2.num_landmarks >= SIZE_MAX / sizeof(int32_t))
        return NULL;

    if (!(c = calloc(1, sizeof(*c))))
        return NULL;

    *c = c2;

    if (c->num_landmarks && !(c->landmark = malloc(c->num_landmarks * sizeof(int32_t)))) {
        fd->err = errno;
        cram_free_container(c);
        return NULL;
    }
    for (i = 0; i < c->num_landmarks; i++) {
        if ((s = fd->vv.varint_decode32_crc(fd, &c->landmark[i], &crc)) == -1) {
            cram_free_container(c);
            return NULL;
        } else {
            rd += s;
        }
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        if (-1 == int32_decode(fd, (int32_t *)&c->crc32)) {
            cram_free_container(c);
            return NULL;
        } else {
            rd+=4;
        }

        if (crc != c->crc32) {
            hts_log_error("Container header CRC32 failure");
            cram_free_container(c);
            return NULL;
        }
    }

    c->offset = rd;
    c->slices = NULL;
    c->slice = NULL;
    c->curr_slice = 0;
    c->max_slice = c->num_landmarks;
    c->slice_rec = 0;
    c->curr_rec = 0;
    c->max_rec = 0;

    if (c->ref_seq_id == -2) {
        c->multi_seq = 1;
        fd->multi_seq = 1;
    }

    fd->empty_container =
        (c->num_records == 0 &&
         c->ref_seq_id == -1 &&
         c->ref_seq_start == 0x454f46 /* EOF */) ? 1 : 0;

    return c;
}


/* MAXIMUM storage size needed for the container. */
int cram_container_size(cram_container *c) {
    return 55 + 5*c->num_landmarks;
}


/*
 * Stores the container structure in dat and returns *size as the
 * number of bytes written to dat[].  The input size of dat is also
 * held in *size and should be initialised to cram_container_size(c).
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_store_container(cram_fd *fd, cram_container *c, char *dat, int *size)
{
    char *cp = (char *)dat;
    int i;

    // Check the input buffer is large enough according to our stated
    // requirements. (NOTE: it may actually take less.)
    if (cram_container_size(c) > *size)
        return -1;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        cp += itf8_put(cp, c->length);
    } else {
        *(int32_t *)cp = le_int4(c->length);
        cp += 4;
    }
    if (c->multi_seq) {
        cp += fd->vv.varint_put32(cp, NULL, -2);
        cp += fd->vv.varint_put32(cp, NULL, 0);
        cp += fd->vv.varint_put32(cp, NULL, 0);
    } else {
        cp += fd->vv.varint_put32s(cp, NULL, c->ref_seq_id);
        if (CRAM_MAJOR_VERS(fd->version) >= 4) {
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_span);
        } else {
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_span);
        }
    }
    cp += fd->vv.varint_put32(cp, NULL, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
        cp += fd->vv.varint_put64(cp, NULL, c->record_counter);
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        cp += fd->vv.varint_put32(cp, NULL, c->record_counter);
    }
    cp += fd->vv.varint_put64(cp, NULL, c->num_bases);
    cp += fd->vv.varint_put32(cp, NULL, c->num_blocks);
    cp += fd->vv.varint_put32(cp, NULL, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
        cp += fd->vv.varint_put32(cp, NULL, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        c->crc32 = crc32(0L, (uc *)dat, cp-dat);
        cp[0] =  c->crc32        & 0xff;
        cp[1] = (c->crc32 >>  8) & 0xff;
        cp[2] = (c->crc32 >> 16) & 0xff;
        cp[3] = (c->crc32 >> 24) & 0xff;
        cp += 4;
    }

    *size = cp-dat; // actual used size

    return 0;
}


/*
 * Writes a container structure.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_container(cram_fd *fd, cram_container *c) {
    char buf_a[1024], *buf = buf_a, *cp;
    int i;

    if (61 + c->num_landmarks * 10 >= 1024) {
        buf = malloc(61 + c->num_landmarks * 10);
        if (!buf)
            return -1;
    }
    cp = buf;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        cp += itf8_put(cp, c->length);
    } else if (CRAM_MAJOR_VERS(fd->version) <= 3) {
        *(int32_t *)cp = le_int4(c->length);
        cp += 4;
    } else {
        cp += fd->vv.varint_put32(cp, NULL, c->length);
    }
    if (c->multi_seq) {
        cp += fd->vv.varint_put32(cp, NULL, (uint32_t)-2);
        cp += fd->vv.varint_put32(cp, NULL, 0);
        cp += fd->vv.varint_put32(cp, NULL, 0);
    } else {
        cp += fd->vv.varint_put32s(cp, NULL, c->ref_seq_id);
        if (CRAM_MAJOR_VERS(fd->version) >= 4) {
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put64(cp, NULL, c->ref_seq_span);
        } else {
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_start);
            cp += fd->vv.varint_put32(cp, NULL, c->ref_seq_span);
        }
    }
    cp += fd->vv.varint_put32(cp, NULL, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) >= 3)
        cp += fd->vv.varint_put64(cp, NULL, c->record_counter);
    else
        cp += fd->vv.varint_put32(cp, NULL, c->record_counter);
    cp += fd->vv.varint_put64(cp, NULL, c->num_bases);
    cp += fd->vv.varint_put32(cp, NULL, c->num_blocks);
    cp += fd->vv.varint_put32(cp, NULL, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
        cp += fd->vv.varint_put32(cp, NULL, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
        c->crc32 = crc32(0L, (uc *)buf, cp-buf);
        cp[0] =  c->crc32        & 0xff;
        cp[1] = (c->crc32 >>  8) & 0xff;
        cp[2] = (c->crc32 >> 16) & 0xff;
        cp[3] = (c->crc32 >> 24) & 0xff;
        cp += 4;
    }

    if (cp-buf != hwrite(fd->fp, buf, cp-buf)) {
        if (buf != buf_a)
            free(buf);
        return -1;
    }

    if (buf != buf_a)
        free(buf);

    return 0;
}

// common component shared by cram_flush_container{,_mt}
static int cram_flush_container2(cram_fd *fd, cram_container *c) {
    int i, j;

    if (c->curr_slice > 0 && !c->slices)
        return -1;

    //fprintf(stderr, "Writing container %d, sum %u\n", c->record_counter, sum);

    off_t c_offset = htell(fd->fp); // File offset of container

    /* Write the container struct itself */
    if (0 != cram_write_container(fd, c))
        return -1;

    off_t hdr_size = htell(fd->fp) - c_offset;

    /* And the compression header */
    if (0 != cram_write_block(fd, c->comp_hdr_block))
        return -1;

    /* Followed by the slice blocks */
    off_t file_offset = htell(fd->fp);
    for (i = 0; i < c->curr_slice; i++) {
        cram_slice *s = c->slices[i];
        off_t spos = file_offset - c_offset - hdr_size;

        if (0 != cram_write_block(fd, s->hdr_block))
            return -1;

        for (j = 0; j < s->hdr->num_blocks; j++) {
            if (0 != cram_write_block(fd, s->block[j]))
                return -1;
        }

        file_offset = htell(fd->fp);
        off_t sz = file_offset - c_offset - hdr_size - spos;

        if (fd->idxfp) {
            if (cram_index_slice(fd, c, s, fd->idxfp, c_offset, spos, sz) < 0)
                return -1;
        }
    }

    return 0;
}

/*
 * Flushes a completely or partially full container to disk, writing
 * container structure, header and blocks. This also calls the encoder
 * functions.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush_container(cram_fd *fd, cram_container *c) {
    /* Encode the container blocks and generate compression header */
    if (0 != cram_encode_container(fd, c))
        return -1;

    return cram_flush_container2(fd, c);
}

typedef struct {
    cram_fd *fd;
    cram_container *c;
} cram_job;

void *cram_flush_thread(void *arg) {
    cram_job *j = (cram_job *)arg;

    /* Encode the container blocks and generate compression header */
    if (0 != cram_encode_container(j->fd, j->c)) {
        hts_log_error("Call to cram_encode_container failed");
        return NULL;
    }

    return arg;
}

static int cram_flush_result(cram_fd *fd) {
    int i, ret = 0;
    hts_tpool_result *r;
    cram_container *lc = NULL;

    // NB: we can have one result per slice, not per container,
    // so we need to free the container only after all slices
    // within it have been freed.  (Automatic via reference counting.)
    while ((r = hts_tpool_next_result(fd->rqueue))) {
        cram_job *j = (cram_job *)hts_tpool_result_data(r);
        cram_container *c;

        if (!j) {
            hts_tpool_delete_result(r, 0);
            return -1;
        }

        fd = j->fd;
        c = j->c;

        if (fd->mode == 'w')
            if (0 != cram_flush_container2(fd, c))
                return -1;

        // Free the slices; filled out by encoder only
        if (c->slices) {
            for (i = 0; i < c->max_slice; i++) {
                if (c->slices[i])
                    cram_free_slice(c->slices[i]);
                if (c->slices[i] == c->slice)
                    c->slice = NULL;
                c->slices[i] = NULL;
            }
        }

        // Free the current slice; set by both encoder & decoder
        if (c->slice) {
            cram_free_slice(c->slice);
            c->slice = NULL;
        }
        c->curr_slice = 0;

        // Our jobs will be in order, so we free the last
        // container when our job has switched to a new one.
        if (c != lc) {
            if (lc) {
                if (fd->ctr == lc)
                    fd->ctr = NULL;
                if (fd->ctr_mt == lc)
                    fd->ctr_mt = NULL;
                cram_free_container(lc);
            }
            lc = c;
        }

        hts_tpool_delete_result(r, 1);
    }
    if (lc) {
        if (fd->ctr == lc)
            fd->ctr = NULL;
        if (fd->ctr_mt == lc)
            fd->ctr_mt = NULL;
        cram_free_container(lc);
    }

    return ret;
}

// Note: called while metrics_lock is held.
// Will be left in this state too, but may temporarily unlock.
void reset_metrics(cram_fd *fd) {
    int i;

    if (fd->pool) {
        // If multi-threaded we have multiple blocks being
        // compressed already and several on the to-do list
        // (fd->rqueue->pending).  It's tricky to reset the
        // metrics exactly the correct point, so instead we
        // just flush the pool, reset, and then continue again.

        // Don't bother starting a new trial before then though.
        for (i = 0; i < DS_END; i++) {
            cram_metrics *m = fd->m[i];
            if (!m)
                continue;
            m->next_trial = 999;
        }

        pthread_mutex_unlock(&fd->metrics_lock);
        hts_tpool_process_flush(fd->rqueue);
        pthread_mutex_lock(&fd->metrics_lock);
    }

    for (i = 0; i < DS_END; i++) {
        cram_metrics *m = fd->m[i];
        if (!m)
            continue;

        m->trial = NTRIALS;
        m->next_trial = TRIAL_SPAN;
        m->revised_method = 0;
        m->unpackable = 0;

        memset(m->sz, 0, sizeof(m->sz));
    }
}

int cram_flush_container_mt(cram_fd *fd, cram_container *c) {
    cram_job *j;

    // At the junction of mapped to unmapped data the compression
    // methods may need to change due to very different statistical
    // properties; particularly BA if minhash sorted.
    //
    // However with threading we'll have several in-flight blocks
    // arriving out of order.
    //
    // So we do one trial reset of NThreads to last for NThreads
    // duration to get us over this transition period, followed
    // by another retrial of the usual ntrials & trial span.
    pthread_mutex_lock(&fd->metrics_lock);
    if (c->n_mapped < 0.3*c->curr_rec &&
        fd->last_mapped > 0.7*c->max_rec) {
        reset_metrics(fd);
    }
    fd->last_mapped = c->n_mapped * (c->max_rec+1)/(c->curr_rec+1) ;
    pthread_mutex_unlock(&fd->metrics_lock);

    if (!fd->pool)
        return cram_flush_container(fd, c);

    if (!(j = malloc(sizeof(*j))))
        return -1;
    j->fd = fd;
    j->c = c;

    // Flush the job.  Note our encoder queue may be full, so we
    // either have to keep trying in non-blocking mode (what we do) or
    // use a dedicated separate thread for draining the queue.
    for (;;) {
        errno = 0;
        hts_tpool_dispatch2(fd->pool, fd->rqueue, cram_flush_thread, j, 1);
        int pending = (errno == EAGAIN);
        if (cram_flush_result(fd) != 0)
            return -1;
        if (!pending)
            break;

        usleep(1000);
    }

    return 0;
}

/* ----------------------------------------------------------------------
 * Compression headers; the first part of the container
 */

/*
 * Creates a new blank container compression header
 *
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_new_compression_header(void) {
    cram_block_compression_hdr *hdr = calloc(1, sizeof(*hdr));
    if (!hdr)
        return NULL;

    if (!(hdr->TD_blk = cram_new_block(CORE, 0))) {
        free(hdr);
        return NULL;
    }

    if (!(hdr->TD_hash = kh_init(m_s2i))) {
        cram_free_block(hdr->TD_blk);
        free(hdr);
        return NULL;
    }

    if (!(hdr->TD_keys = string_pool_create(8192))) {
        kh_destroy(m_s2i, hdr->TD_hash);
        cram_free_block(hdr->TD_blk);
        free(hdr);
        return NULL;
    }

    return hdr;
}

void cram_free_compression_header(cram_block_compression_hdr *hdr) {
    int i;

    if (hdr->landmark)
        free(hdr->landmark);

    if (hdr->preservation_map)
        kh_destroy(map, hdr->preservation_map);

    for (i = 0; i < CRAM_MAP_HASH; i++) {
        cram_map *m, *m2;
        for (m = hdr->rec_encoding_map[i]; m; m = m2) {
            m2 = m->next;
            if (m->codec)
                m->codec->free(m->codec);
            free(m);
        }
    }

    for (i = 0; i < CRAM_MAP_HASH; i++) {
        cram_map *m, *m2;
        for (m = hdr->tag_encoding_map[i]; m; m = m2) {
            m2 = m->next;
            if (m->codec)
                m->codec->free(m->codec);
            free(m);
        }
    }

    for (i = 0; i < DS_END; i++) {
        if (hdr->codecs[i])
            hdr->codecs[i]->free(hdr->codecs[i]);
    }

    if (hdr->TL)
        free(hdr->TL);
    if (hdr->TD_blk)
        cram_free_block(hdr->TD_blk);
    if (hdr->TD_hash)
        kh_destroy(m_s2i, hdr->TD_hash);
    if (hdr->TD_keys)
        string_pool_destroy(hdr->TD_keys);

    free(hdr);
}


/* ----------------------------------------------------------------------
 * Slices and slice headers
 */

void cram_free_slice_header(cram_block_slice_hdr *hdr) {
    if (!hdr)
        return;

    if (hdr->block_content_ids)
        free(hdr->block_content_ids);

    free(hdr);

    return;
}

void cram_free_slice(cram_slice *s) {
    if (!s)
        return;

    if (s->hdr_block)
        cram_free_block(s->hdr_block);

    if (s->block) {
        int i;

        if (s->hdr) {
            for (i = 0; i < s->hdr->num_blocks; i++) {
                if (i > 0 && s->block[i] == s->block[0])
                    continue;
                cram_free_block(s->block[i]);
            }
        }
        free(s->block);
    }

    if (s->block_by_id)
        free(s->block_by_id);

    if (s->hdr)
        cram_free_slice_header(s->hdr);

    if (s->seqs_blk)
        cram_free_block(s->seqs_blk);

    if (s->qual_blk)
        cram_free_block(s->qual_blk);

    if (s->name_blk)
        cram_free_block(s->name_blk);

    if (s->aux_blk)
        cram_free_block(s->aux_blk);

    if (s->base_blk)
        cram_free_block(s->base_blk);

    if (s->soft_blk)
        cram_free_block(s->soft_blk);

    if (s->cigar)
        free(s->cigar);

    if (s->crecs)
        free(s->crecs);

    if (s->features)
        free(s->features);

    if (s->TN)
        free(s->TN);

    if (s->pair_keys)
        string_pool_destroy(s->pair_keys);

    if (s->pair[0])
        kh_destroy(m_s2i, s->pair[0]);
    if (s->pair[1])
        kh_destroy(m_s2i, s->pair[1]);

    if (s->aux_block)
        free(s->aux_block);

    free(s);
}

/*
 * Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs) {
    cram_slice *s = calloc(1, sizeof(*s));
    if (!s)
        return NULL;

    if (!(s->hdr = (cram_block_slice_hdr *)calloc(1, sizeof(*s->hdr))))
        goto err;
    s->hdr->content_type = type;

    s->hdr_block = NULL;
    s->block = NULL;
    s->block_by_id = NULL;
    s->last_apos = 0;
    if (!(s->crecs = malloc(nrecs * sizeof(cram_record))))  goto err;
    s->cigar_alloc = 1024;
    if (!(s->cigar = malloc(s->cigar_alloc * sizeof(*s->cigar)))) goto err;
    s->ncigar = 0;

    if (!(s->seqs_blk = cram_new_block(EXTERNAL, 0)))       goto err;
    if (!(s->qual_blk = cram_new_block(EXTERNAL, DS_QS)))   goto err;
    if (!(s->name_blk = cram_new_block(EXTERNAL, DS_RN)))   goto err;
    if (!(s->aux_blk  = cram_new_block(EXTERNAL, DS_aux)))  goto err;
    if (!(s->base_blk = cram_new_block(EXTERNAL, DS_IN)))   goto err;
    if (!(s->soft_blk = cram_new_block(EXTERNAL, DS_SC)))   goto err;

    s->features = NULL;
    s->nfeatures = s->afeatures = 0;

#ifndef TN_external
    s->TN = NULL;
    s->nTN = s->aTN = 0;
#endif

    // Volatile keys as we do realloc in dstring
    if (!(s->pair_keys = string_pool_create(8192))) goto err;
    if (!(s->pair[0] = kh_init(m_s2i)))             goto err;
    if (!(s->pair[1] = kh_init(m_s2i)))             goto err;

#ifdef BA_external
    s->BA_len = 0;
#endif

    return s;

 err:
    if (s)
        cram_free_slice(s);

    return NULL;
}

/*
 * Loads an entire slice.
 * FIXME: In 1.0 the native unit of slices within CRAM is broken
 * as slices contain references to objects in other slices.
 * To work around this while keeping the slice oriented outer loop
 * we read all slices and stitch them together into a fake large
 * slice instead.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd) {
    cram_block *b = cram_read_block(fd);
    cram_slice *s = calloc(1, sizeof(*s));
    int i, n, max_id, min_id;

    if (!b || !s)
        goto err;

    s->hdr_block = b;
    switch (b->content_type) {
    case MAPPED_SLICE:
    case UNMAPPED_SLICE:
        if (!(s->hdr = cram_decode_slice_header(fd, b)))
            goto err;
        break;

    default:
        hts_log_error("Unexpected block of type %s",
                      cram_content_type2str(b->content_type));
        goto err;
    }

    if (s->hdr->num_blocks < 1) {
        hts_log_error("Slice does not include any data blocks");
        goto err;
    }

    s->block = calloc(n = s->hdr->num_blocks, sizeof(*s->block));
    if (!s->block)
        goto err;

    for (max_id = i = 0, min_id = INT_MAX; i < n; i++) {
        if (!(s->block[i] = cram_read_block(fd)))
            goto err;

        if (s->block[i]->content_type == EXTERNAL) {
            if (max_id < s->block[i]->content_id)
                max_id = s->block[i]->content_id;
            if (min_id > s->block[i]->content_id)
                min_id = s->block[i]->content_id;
        }
    }

    if (!(s->block_by_id = calloc(512, sizeof(s->block[0]))))
        goto err;

    for (i = 0; i < n; i++) {
        if (s->block[i]->content_type != EXTERNAL)
            continue;
        uint32_t v = s->block[i]->content_id;
        if (v >= 256)
            v = 256 + v % 251;
        s->block_by_id[v] = s->block[i];
    }

    /* Initialise encoding/decoding tables */
    s->cigar_alloc = 1024;
    if (!(s->cigar = malloc(s->cigar_alloc * sizeof(*s->cigar)))) goto err;
    s->ncigar = 0;

    if (!(s->seqs_blk = cram_new_block(EXTERNAL, 0)))      goto err;
    if (!(s->qual_blk = cram_new_block(EXTERNAL, DS_QS)))  goto err;
    if (!(s->name_blk = cram_new_block(EXTERNAL, DS_RN)))  goto err;
    if (!(s->aux_blk  = cram_new_block(EXTERNAL, DS_aux))) goto err;
    if (!(s->base_blk = cram_new_block(EXTERNAL, DS_IN)))  goto err;
    if (!(s->soft_blk = cram_new_block(EXTERNAL, DS_SC)))  goto err;

    s->crecs = NULL;

    s->last_apos = s->hdr->ref_seq_start;
    s->decode_md = fd->decode_md;

    return s;

 err:
    if (b)
        cram_free_block(b);
    if (s) {
        s->hdr_block = NULL;
        cram_free_slice(s);
    }
    return NULL;
}


/* ----------------------------------------------------------------------
 * CRAM file definition (header)
 */

/*
 * Reads a CRAM file definition structure.
 * Returns file_def ptr on success
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd) {
    cram_file_def *def = malloc(sizeof(*def));
    if (!def)
        return NULL;

    if (26 != hread(fd->fp, &def->magic[0], 26)) {
        free(def);
        return NULL;
    }

    if (memcmp(def->magic, "CRAM", 4) != 0) {
        free(def);
        return NULL;
    }

    if (def->major_version > 4) {
        hts_log_error("CRAM version number mismatch. Expected 1.x, 2.x, 3.x or 4.x, got %d.%d",
                      def->major_version, def->minor_version);
        free(def);
        return NULL;
    }

    fd->first_container += 26;
    fd->curr_position = fd->first_container;
    fd->last_slice = 0;

    return def;
}

/*
 * Writes a cram_file_def structure to cram_fd.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def) {
    return (hwrite(fd->fp, &def->magic[0], 26) == 26) ? 0 : -1;
}

void cram_free_file_def(cram_file_def *def) {
    if (def) free(def);
}

/* ----------------------------------------------------------------------
 * SAM header I/O
 */


/*
 * Reads the SAM header from the first CRAM data block.
 * Also performs minimal parsing to extract read-group
 * and sample information.

 * Returns SAM hdr ptr on success
 *         NULL on failure
 */
sam_hdr_t *cram_read_SAM_hdr(cram_fd *fd) {
    int32_t header_len;
    char *header;
    sam_hdr_t *hdr;

    /* 1.1 onwards stores the header in the first block of a container */
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        /* Length */
        if (-1 == int32_decode(fd, &header_len))
            return NULL;

        /* Alloc and read */
        if (header_len < 0 || NULL == (header = malloc((size_t) header_len+1)))
            return NULL;

        if (header_len != hread(fd->fp, header, header_len)) {
            free(header);
            return NULL;
        }
        header[header_len] = '\0';

        fd->first_container += 4 + header_len;
    } else {
        cram_container *c = cram_read_container(fd);
        cram_block *b;
        int i;
        int64_t len;

        if (!c)
            return NULL;

        fd->first_container += c->length + c->offset;
        fd->curr_position = fd->first_container;

        if (c->num_blocks < 1) {
            cram_free_container(c);
            return NULL;
        }

        if (!(b = cram_read_block(fd))) {
            cram_free_container(c);
            return NULL;
        }
        if (cram_uncompress_block(b) != 0) {
            cram_free_container(c);
            cram_free_block(b);
            return NULL;
        }

        len = b->comp_size + 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
            fd->vv.varint_size(b->content_id) +
            fd->vv.varint_size(b->uncomp_size) +
            fd->vv.varint_size(b->comp_size);

        /* Extract header from 1st block */
        if (-1 == int32_get_blk(b, &header_len) ||
            header_len < 0 || /* Spec. says signed...  why? */
            b->uncomp_size - 4 < header_len) {
            cram_free_container(c);
            cram_free_block(b);
            return NULL;
        }
        if (NULL == (header = malloc((size_t) header_len+1))) {
            cram_free_container(c);
            cram_free_block(b);
            return NULL;
        }
        memcpy(header, BLOCK_END(b), header_len);
        header[header_len] = '\0';
        cram_free_block(b);

        /* Consume any remaining blocks */
        for (i = 1; i < c->num_blocks; i++) {
            if (!(b = cram_read_block(fd))) {
                cram_free_container(c);
                free(header);
                return NULL;
            }
            len += b->comp_size + 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
                fd->vv.varint_size(b->content_id) +
                fd->vv.varint_size(b->uncomp_size) +
                fd->vv.varint_size(b->comp_size);
            cram_free_block(b);
        }

        if (c->length > 0 && len > 0 && c->length > len) {
            // Consume padding
            char *pads = malloc(c->length - len);
            if (!pads) {
                cram_free_container(c);
                free(header);
                return NULL;
            }

            if (c->length - len != hread(fd->fp, pads, c->length - len)) {
                cram_free_container(c);
                free(header);
                free(pads);
                return NULL;
            }
            free(pads);
        }

        cram_free_container(c);
    }

    /* Parse */
    hdr = sam_hdr_init();
    if (!hdr) {
        free(header);
        return NULL;
    }

    if (-1 == sam_hdr_add_lines(hdr, header, header_len)) {
        free(header);
        sam_hdr_destroy(hdr);
        return NULL;
    }

    hdr->l_text = header_len;
    hdr->text = header;

    return hdr;

}

/*
 * Converts 'in' to a full pathname to store in out.
 * Out must be at least PATH_MAX bytes long.
 */
static void full_path(char *out, char *in) {
    size_t in_l = strlen(in);
    if (hisremote(in)) {
        if (in_l > PATH_MAX) {
            hts_log_error("Reference path is longer than %d", PATH_MAX);
            return;
        }
        strncpy(out, in, PATH_MAX-1);
        out[PATH_MAX-1] = 0;
        return;
    }
    if (*in == '/' ||
        // Windows paths
        (in_l > 3 && toupper_c(*in) >= 'A'  && toupper_c(*in) <= 'Z' &&
         in[1] == ':' && (in[2] == '/' || in[2] == '\\'))) {
        strncpy(out, in, PATH_MAX-1);
        out[PATH_MAX-1] = 0;
    } else {
        int len;

        // unable to get dir or out+in is too long
        if (!getcwd(out, PATH_MAX) ||
            (len = strlen(out))+1+strlen(in) >= PATH_MAX) {
            strncpy(out, in, PATH_MAX-1);
            out[PATH_MAX-1] = 0;
            return;
        }

        sprintf(out+len, "/%.*s", PATH_MAX - 2 - len, in);

        // FIXME: cope with `pwd`/../../../foo.fa ?
    }
}

/*
 * Writes a CRAM SAM header.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, sam_hdr_t *hdr) {
    size_t header_len;
    int blank_block = (CRAM_MAJOR_VERS(fd->version) >= 3);

    /* Write CRAM MAGIC if not yet written. */
    if (fd->file_def->major_version == 0) {
        fd->file_def->major_version = CRAM_MAJOR_VERS(fd->version);
        fd->file_def->minor_version = CRAM_MINOR_VERS(fd->version);
        if (0 != cram_write_file_def(fd, fd->file_def))
            return -1;
    }

    /* 1.0 requires an UNKNOWN read-group */
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        if (!sam_hrecs_find_rg(hdr->hrecs, "UNKNOWN"))
            if (sam_hdr_add_line(hdr, "RG",
                            "ID", "UNKNOWN", "SM", "UNKNOWN", NULL))
                return -1;
    }

    if (-1 == refs_from_header(fd))
        return -1;
    if (-1 == refs2id(fd->refs, fd->header))
        return -1;

    /* Fix M5 strings */
    if (fd->refs && !fd->no_ref && fd->embed_ref <= 1) {
        int i;
        for (i = 0; i < hdr->hrecs->nref; i++) {
            sam_hrec_type_t *ty;
            char *ref;

            if (!(ty = sam_hrecs_find_type_id(hdr->hrecs, "SQ", "SN", hdr->hrecs->ref[i].name)))
                return -1;

            if (!sam_hrecs_find_key(ty, "M5", NULL)) {
                char unsigned buf[16];
                char buf2[33];
                int rlen;
                hts_md5_context *md5;

                if (!fd->refs ||
                    !fd->refs->ref_id ||
                    !fd->refs->ref_id[i]) {
                    return -1;
                }
                rlen = fd->refs->ref_id[i]->length;
                ref = cram_get_ref(fd, i, 1, rlen);
                if (NULL == ref) {
                    if (fd->embed_ref == -1) {
                        // auto embed-ref
                        hts_log_warning("No M5 tags present and could not "
                                        "find reference");
                        hts_log_warning("Enabling embed_ref=2 option");
                        hts_log_warning("NOTE: the CRAM file will be bigger "
                                        "than using an external reference");
                        pthread_mutex_lock(&fd->ref_lock);
                        fd->embed_ref = 2;
                        pthread_mutex_unlock(&fd->ref_lock);
                        break;
                    }
                    return -1;
                }
                rlen = fd->refs->ref_id[i]->length; /* In case it just loaded */
                if (!(md5 = hts_md5_init()))
                    return -1;
                hts_md5_update(md5, ref, rlen);
                hts_md5_final(buf, md5);
                hts_md5_destroy(md5);
                cram_ref_decr(fd->refs, i);

                hts_md5_hex(buf2, buf);
                fd->refs->ref_id[i]->validated_md5 = 1;
                if (sam_hdr_update_line(hdr, "SQ", "SN", hdr->hrecs->ref[i].name, "M5", buf2, NULL))
                    return -1;
            }

            if (fd->ref_fn) {
                char ref_fn[PATH_MAX];
                full_path(ref_fn, fd->ref_fn);
                if (sam_hdr_update_line(hdr, "SQ", "SN", hdr->hrecs->ref[i].name, "UR", ref_fn, NULL))
                    return -1;
            }
        }
    }

    /* Length */
    header_len = sam_hdr_length(hdr);
    if (header_len > INT32_MAX) {
        hts_log_error("Header is too long for CRAM format");
        return -1;
    }
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        if (-1 == int32_encode(fd, header_len))
            return -1;

        /* Text data */
        if (header_len != hwrite(fd->fp, sam_hdr_str(hdr), header_len))
            return -1;
    } else {
        /* Create block(s) inside a container */
        cram_block *b = cram_new_block(FILE_HEADER, 0);
        cram_container *c = cram_new_container(0, 0);
        int padded_length;
        char *pads;
        int is_cram_3 = (CRAM_MAJOR_VERS(fd->version) >= 3);

        if (!b || !c) {
            if (b) cram_free_block(b);
            if (c) cram_free_container(c);
            return -1;
        }

        if (int32_put_blk(b, header_len) < 0)
            return -1;
        if (header_len)
            BLOCK_APPEND(b, sam_hdr_str(hdr), header_len);
        BLOCK_UPLEN(b);

        // Compress header block if V3.0 and above
        if (CRAM_MAJOR_VERS(fd->version) >= 3)
            if (cram_compress_block(fd, b, NULL, -1, -1) < 0)
                return -1;

        if (blank_block) {
            c->length = b->comp_size + 2 + 4*is_cram_3 +
                fd->vv.varint_size(b->content_id)   +
                fd->vv.varint_size(b->uncomp_size)  +
                fd->vv.varint_size(b->comp_size);

            c->num_blocks = 2;
            c->num_landmarks = 2;
            if (!(c->landmark = malloc(2*sizeof(*c->landmark)))) {
                cram_free_block(b);
                cram_free_container(c);
                return -1;
            }
            c->landmark[0] = 0;
            c->landmark[1] = c->length;

            // Plus extra storage for uncompressed secondary blank block
            padded_length = MIN(c->length*.5, 10000);
            c->length += padded_length + 2 + 4*is_cram_3 +
                fd->vv.varint_size(b->content_id) +
                fd->vv.varint_size(padded_length)*2;
        } else {
            // Pad the block instead.
            c->num_blocks = 1;
            c->num_landmarks = 1;
            if (!(c->landmark = malloc(sizeof(*c->landmark))))
                return -1;
            c->landmark[0] = 0;

            padded_length = MAX(c->length*1.5, 10000) - c->length;

            c->length = b->comp_size + padded_length +
                2 + 4*is_cram_3 +
                fd->vv.varint_size(b->content_id)   +
                fd->vv.varint_size(b->uncomp_size)  +
                fd->vv.varint_size(b->comp_size);

            if (NULL == (pads = calloc(1, padded_length))) {
                cram_free_block(b);
                cram_free_container(c);
                return -1;
            }
            BLOCK_APPEND(b, pads, padded_length);
            BLOCK_UPLEN(b);
            free(pads);
        }

        if (-1 == cram_write_container(fd, c)) {
            cram_free_block(b);
            cram_free_container(c);
            return -1;
        }

        if (-1 == cram_write_block(fd, b)) {
            cram_free_block(b);
            cram_free_container(c);
            return -1;
        }

        if (blank_block) {
            BLOCK_RESIZE(b, padded_length);
            memset(BLOCK_DATA(b), 0, padded_length);
            BLOCK_SIZE(b) = padded_length;
            BLOCK_UPLEN(b);
            b->method = RAW;
            if (-1 == cram_write_block(fd, b)) {
                cram_free_block(b);
                cram_free_container(c);
                return -1;
            }
        }

        cram_free_block(b);
        cram_free_container(c);
    }

    if (0 != hflush(fd->fp))
        return -1;

    RP("=== Finishing saving header ===\n");

    return 0;

 block_err:
    return -1;
}

/* ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*
 * Sets CRAM variable sized integer decode function tables.
 * CRAM 1, 2, and 3.x all used ITF8 for uint32 and UTF8 for uint64.
 * CRAM 4.x uses the same encoding mechanism for 32-bit and 64-bit
 * (or anything inbetween), but also now supports signed values.
 *
 * Version is the CRAM major version number.
 * vv is the vector table (probably &cram_fd->vv)
 */
static void cram_init_varint(varint_vec *vv, int version) {
    if (version >= 4) {
        vv->varint_get32 = uint7_get_32; // FIXME: varint.h API should be size agnostic
        vv->varint_get32s = sint7_get_32;
        vv->varint_get64 = uint7_get_64;
        vv->varint_get64s = sint7_get_64;
        vv->varint_put32 = uint7_put_32;
        vv->varint_put32s = sint7_put_32;
        vv->varint_put64 = uint7_put_64;
        vv->varint_put64s = sint7_put_64;
        vv->varint_put32_blk = uint7_put_blk_32;
        vv->varint_put32s_blk = sint7_put_blk_32;
        vv->varint_put64_blk = uint7_put_blk_64;
        vv->varint_put64s_blk = sint7_put_blk_64;
        vv->varint_size = uint7_size;
        vv->varint_decode32_crc = uint7_decode_crc32;
        vv->varint_decode32s_crc = sint7_decode_crc32;
        vv->varint_decode64_crc = uint7_decode_crc64;
    } else {
        vv->varint_get32 = safe_itf8_get;
        vv->varint_get32s = safe_itf8_get;
        vv->varint_get64 = safe_ltf8_get;
        vv->varint_get64s = safe_ltf8_get;
        vv->varint_put32 = safe_itf8_put;
        vv->varint_put32s = safe_itf8_put;
        vv->varint_put64 = safe_ltf8_put;
        vv->varint_put64s = safe_ltf8_put;
        vv->varint_put32_blk = itf8_put_blk;
        vv->varint_put32s_blk = itf8_put_blk;
        vv->varint_put64_blk = ltf8_put_blk;
        vv->varint_put64s_blk = ltf8_put_blk;
        vv->varint_size = itf8_size;
        vv->varint_decode32_crc = itf8_decode_crc;
        vv->varint_decode32s_crc = itf8_decode_crc;
        vv->varint_decode64_crc = ltf8_decode_crc;
    }
}

/*
 * Initialises the lookup tables. These could be global statics, but they're
 * clumsy to setup in a multi-threaded environment unless we generate
 * verbatim code and include that.
 */
static void cram_init_tables(cram_fd *fd) {
    int i;

    memset(fd->L1, 4, 256);
    fd->L1['A'] = 0; fd->L1['a'] = 0;
    fd->L1['C'] = 1; fd->L1['c'] = 1;
    fd->L1['G'] = 2; fd->L1['g'] = 2;
    fd->L1['T'] = 3; fd->L1['t'] = 3;

    memset(fd->L2, 5, 256);
    fd->L2['A'] = 0; fd->L2['a'] = 0;
    fd->L2['C'] = 1; fd->L2['c'] = 1;
    fd->L2['G'] = 2; fd->L2['g'] = 2;
    fd->L2['T'] = 3; fd->L2['t'] = 3;
    fd->L2['N'] = 4; fd->L2['n'] = 4;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
        for (i = 0; i < 0x200; i++) {
            int f = 0;

            if (i & CRAM_FPAIRED)      f |= BAM_FPAIRED;
            if (i & CRAM_FPROPER_PAIR) f |= BAM_FPROPER_PAIR;
            if (i & CRAM_FUNMAP)       f |= BAM_FUNMAP;
            if (i & CRAM_FREVERSE)     f |= BAM_FREVERSE;
            if (i & CRAM_FREAD1)       f |= BAM_FREAD1;
            if (i & CRAM_FREAD2)       f |= BAM_FREAD2;
            if (i & CRAM_FSECONDARY)   f |= BAM_FSECONDARY;
            if (i & CRAM_FQCFAIL)      f |= BAM_FQCFAIL;
            if (i & CRAM_FDUP)         f |= BAM_FDUP;

            fd->bam_flag_swap[i]  = f;
        }

        for (i = 0; i < 0x1000; i++) {
            int g = 0;

            if (i & BAM_FPAIRED)           g |= CRAM_FPAIRED;
            if (i & BAM_FPROPER_PAIR)  g |= CRAM_FPROPER_PAIR;
            if (i & BAM_FUNMAP)        g |= CRAM_FUNMAP;
            if (i & BAM_FREVERSE)      g |= CRAM_FREVERSE;
            if (i & BAM_FREAD1)        g |= CRAM_FREAD1;
            if (i & BAM_FREAD2)        g |= CRAM_FREAD2;
            if (i & BAM_FSECONDARY)    g |= CRAM_FSECONDARY;
            if (i & BAM_FQCFAIL)       g |= CRAM_FQCFAIL;
            if (i & BAM_FDUP)          g |= CRAM_FDUP;

            fd->cram_flag_swap[i] = g;
        }
    } else {
        /* NOP */
        for (i = 0; i < 0x1000; i++)
            fd->bam_flag_swap[i] = i;
        for (i = 0; i < 0x1000; i++)
            fd->cram_flag_swap[i] = i;
    }

    memset(fd->cram_sub_matrix, 4, 32*32);
    for (i = 0; i < 32; i++) {
        fd->cram_sub_matrix[i]['A'&0x1f]=0;
        fd->cram_sub_matrix[i]['C'&0x1f]=1;
        fd->cram_sub_matrix[i]['G'&0x1f]=2;
        fd->cram_sub_matrix[i]['T'&0x1f]=3;
        fd->cram_sub_matrix[i]['N'&0x1f]=4;
    }
    for (i = 0; i < 20; i+=4) {
        int j;
        for (j = 0; j < 20; j++) {
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
            fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
        }
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+0]&0x1f]=0;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+1]&0x1f]=1;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+2]&0x1f]=2;
        fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+3]&0x1f]=3;
    }

    cram_init_varint(&fd->vv, CRAM_MAJOR_VERS(fd->version));
}

// Default version numbers for CRAM
static int major_version = 3;
static int minor_version = 0;

/*
 * Opens a CRAM file for read (mode "rb") or write ("wb").
 * The filename may be "-" to indicate stdin or stdout.
 *
 * Returns file handle on success
 *         NULL on failure.
 */
cram_fd *cram_open(const char *filename, const char *mode) {
    hFILE *fp;
    cram_fd *fd;
    char fmode[3]= { mode[0], '\0', '\0' };

    if (strlen(mode) > 1 && (mode[1] == 'b' || mode[1] == 'c')) {
        fmode[1] = 'b';
    }

    fp = hopen(filename, fmode);
    if (!fp)
        return NULL;

    fd = cram_dopen(fp, filename, mode);
    if (!fd)
        hclose_abruptly(fp);

    return fd;
}

/* Opens an existing stream for reading or writing.
 *
 * Returns file handle on success;
 *         NULL on failure.
 */
cram_fd *cram_dopen(hFILE *fp, const char *filename, const char *mode) {
    int i;
    char *cp;
    cram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
        return NULL;

    fd->level = CRAM_DEFAULT_LEVEL;
    for (i = 0; mode[i]; i++) {
        if (mode[i] >= '0' && mode[i] <= '9') {
            fd->level = mode[i] - '0';
            break;
        }
    }

    fd->fp = fp;
    fd->mode = *mode;
    fd->first_container = 0;
    fd->curr_position = 0;

    if (fd->mode == 'r') {
        /* Reader */

        if (!(fd->file_def = cram_read_file_def(fd)))
            goto err;

        fd->version = fd->file_def->major_version * 256 +
            fd->file_def->minor_version;

        cram_init_tables(fd);

        if (!(fd->header = cram_read_SAM_hdr(fd))) {
            cram_free_file_def(fd->file_def);
            goto err;
        }

    } else {
        /* Writer */
        cram_file_def *def = calloc(1, sizeof(*def));
        if (!def)
            return NULL;

        fd->file_def = def;

        def->magic[0] = 'C';
        def->magic[1] = 'R';
        def->magic[2] = 'A';
        def->magic[3] = 'M';
        def->major_version = 0; // Indicator to write file def later.
        def->minor_version = 0;
        memset(def->file_id, 0, 20);
        strncpy(def->file_id, filename, 20);

        fd->version = major_version * 256 + minor_version;
        cram_init_tables(fd);

        /* SAM header written later along with this file_def */
    }

    fd->prefix = strdup((cp = strrchr(filename, '/')) ? cp+1 : filename);
    if (!fd->prefix)
        goto err;
    fd->first_base = fd->last_base = -1;
    fd->record_counter = 0;

    fd->ctr = NULL;
    fd->ctr_mt = NULL;
    fd->refs  = refs_create();
    if (!fd->refs)
        goto err;
    fd->ref_id = -2;
    fd->ref = NULL;

    fd->decode_md = 0;
    fd->seqs_per_slice = SEQS_PER_SLICE;
    fd->bases_per_slice = BASES_PER_SLICE;
    fd->slices_per_container = SLICE_PER_CNT;
    fd->embed_ref = -1; // automatic selection
    fd->no_ref = 0;
    fd->ap_delta = 0;
    fd->ignore_md5 = 0;
    fd->lossy_read_names = 0;
    fd->use_bz2 = 0;
    fd->use_rans = (CRAM_MAJOR_VERS(fd->version) >= 3);
    fd->use_tok = (CRAM_MAJOR_VERS(fd->version) >= 3) && (CRAM_MINOR_VERS(fd->version) >= 1);
    fd->use_lzma = 0;
    fd->multi_seq = -1;
    fd->multi_seq_user = -1;
    fd->unsorted   = 0;
    fd->shared_ref = 0;
    fd->store_md = 0;
    fd->store_nm = 0;
    fd->last_RI_count = 0;

    fd->index       = NULL;
    fd->own_pool    = 0;
    fd->pool        = NULL;
    fd->rqueue      = NULL;
    fd->job_pending = NULL;
    fd->ooc         = 0;
    fd->required_fields = INT_MAX;

    for (i = 0; i < DS_END; i++) {
        fd->m[i] = cram_new_metrics();
        if (!fd->m[i])
            goto err;
    }

    if (!(fd->tags_used = kh_init(m_metrics)))
        goto err;

    fd->range.refid = -2; // no ref.
    fd->eof = 1;          // See samtools issue #150
    fd->ref_fn = NULL;

    fd->bl = NULL;

    /* Initialise dummy refs from the @SQ headers */
    if (-1 == refs_from_header(fd))
        goto err;

    return fd;

 err:
    if (fd)
        free(fd);

    return NULL;
}

/*
 * Seek within a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence) {
    char buf[65536];

    fd->ooc = 0;

    cram_drain_rqueue(fd);

    if (hseek(fd->fp, offset, whence) >= 0) {
        return 0;
    }

    if (!(whence == SEEK_CUR && offset >= 0))
        return -1;

    /* Couldn't fseek, but we're in SEEK_CUR mode so read instead */
    while (offset > 0) {
        int len = MIN(65536, offset);
        if (len != hread(fd->fp, buf, len))
            return -1;
        offset -= len;
    }

    return 0;
}

/*
 * Flushes a CRAM file.
 * Useful for when writing to stdout without wishing to close the stream.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush(cram_fd *fd) {
    if (!fd)
        return -1;

    if (fd->mode == 'w' && fd->ctr) {
        if(fd->ctr->slice)
            cram_update_curr_slice(fd->ctr, fd->version);

        if (-1 == cram_flush_container_mt(fd, fd->ctr))
            return -1;
    }

    return 0;
}

/*
 * Writes an EOF block to a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_eof_block(cram_fd *fd) {
    // EOF block is a container with special values to aid detection
    if (CRAM_MAJOR_VERS(fd->version) >= 2) {
        // Empty container with
        //   ref_seq_id -1
        //   start pos 0x454f46 ("EOF")
        //   span 0
        //   nrec 0
        //   counter 0
        //   nbases 0
        //   1 block (landmark 0)
        //   (CRC32)
        cram_container c;
        memset(&c, 0, sizeof(c));
        c.ref_seq_id = -1;
        c.ref_seq_start = 0x454f46; // "EOF"
        c.ref_seq_span = 0;
        c.record_counter = 0;
        c.num_bases = 0;
        c.num_blocks = 1;
        int32_t land[1] = {0};
        c.landmark = land;

        // An empty compression header block with
        //   method raw (0)
        //   type comp header (1)
        //   content id 0
        //   block contents size 6
        //   raw size 6
        //     empty preservation map (01 00)
        //     empty data series map (01 00)
        //     empty tag map (01 00)
        //   block CRC
        cram_block_compression_hdr ch;
        memset(&ch, 0, sizeof(ch));
        c.comp_hdr_block = cram_encode_compression_header(fd, &c, &ch, 0);

        c.length = c.comp_hdr_block->byte            // Landmark[0]
            + 5                                      // block struct
            + 4*(CRAM_MAJOR_VERS(fd->version) >= 3); // CRC
        if (cram_write_container(fd, &c) < 0 ||
            cram_write_block(fd, c.comp_hdr_block) < 0) {
            cram_close(fd);
            cram_free_block(c.comp_hdr_block);
            return -1;
        }
        if (ch.preservation_map)
            kh_destroy(map, ch.preservation_map);
        cram_free_block(c.comp_hdr_block);

        // V2.1 bytes
        // 0b 00 00 00 ff ff ff ff 0f // Cont HDR: size, ref seq id
        // e0 45 4f 46 00 00 00       // Cont HDR: pos, span, nrec, counter
        // 00 01 00                   // Cont HDR: nbase, nblk, landmark
        // 00 01 00 06 06             // Comp.HDR blk
        // 01 00 01 00 01 00          // Comp.HDR blk

        // V3.0 bytes:
        // 0f 00 00 00 ff ff ff ff 0f // Cont HDR: size, ref seq id
        // e0 45 4f 46 00 00 00       // Cont HDR: pos, span, nrec, counter
        // 00 01 00                   // Cont HDR: nbase, nblk, landmark
        // 05 bd d9 4f                // CRC32
        // 00 01 00 06 06             // Comp.HDR blk
        // 01 00 01 00 01 00          // Comp.HDR blk
        // ee 63 01 4b                // CRC32

        // V4.0 bytes:
        // 0f 00 00 00 8f ff ff ff    // Cont HDR: size, ref seq id
        // 82 95 9e 46 00 00 00       // Cont HDR: pos, span, nrec, counter
        // 00 01 00                   // Cont HDR: nbase, nblk, landmark
        // ac d6 05 bc                // CRC32
        // 00 01 00 06 06             // Comp.HDR blk
        // 01 00 01 00 01 00          // Comp.HDR blk
        // ee 63 01 4b                // CRC32
    }

    return 0;
}
/*
 * Closes a CRAM file.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_close(cram_fd *fd) {
    spare_bams *bl, *next;
    int i;

    if (!fd)
        return -1;

    if (fd->mode == 'w' && fd->ctr) {
        if(fd->ctr->slice)
            cram_update_curr_slice(fd->ctr, fd->version);

        if (-1 == cram_flush_container_mt(fd, fd->ctr))
            return -1;
    }

    if (fd->mode != 'w')
        cram_drain_rqueue(fd);

    if (fd->pool && fd->eof >= 0 && fd->rqueue) {
        hts_tpool_process_flush(fd->rqueue);

        if (0 != cram_flush_result(fd))
            return -1;

        if (fd->mode == 'w')
            fd->ctr = NULL; // prevent double freeing

        pthread_mutex_destroy(&fd->metrics_lock);
        pthread_mutex_destroy(&fd->ref_lock);
        pthread_mutex_destroy(&fd->bam_list_lock);

        //fprintf(stderr, "CRAM: destroy queue %p\n", fd->rqueue);

        hts_tpool_process_destroy(fd->rqueue);
    }

    if (fd->mode == 'w') {
        /* Write EOF block */
        if (0 != cram_write_eof_block(fd))
            return -1;
    }

    for (bl = fd->bl; bl; bl = next) {
        int i, max_rec = fd->seqs_per_slice * fd->slices_per_container;

        next = bl->next;
        for (i = 0; i < max_rec; i++) {
            if (bl->bams[i])
                bam_free(bl->bams[i]);
        }
        free(bl->bams);
        free(bl);
    }

    if (hclose(fd->fp) != 0)
        return -1;

    if (fd->file_def)
        cram_free_file_def(fd->file_def);

    if (fd->header)
        sam_hdr_destroy(fd->header);

    free(fd->prefix);

    if (fd->ctr)
        cram_free_container(fd->ctr);

    if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
        cram_free_container(fd->ctr_mt);

    if (fd->refs)
        refs_free(fd->refs);
    if (fd->ref_free)
        free(fd->ref_free);

    for (i = 0; i < DS_END; i++)
        if (fd->m[i])
            free(fd->m[i]);

    if (fd->tags_used) {
        khint_t k;

        for (k = kh_begin(fd->tags_used); k != kh_end(fd->tags_used); k++) {
            if (kh_exist(fd->tags_used, k))
                free(kh_val(fd->tags_used, k));
        }

        kh_destroy(m_metrics, fd->tags_used);
    }

    if (fd->index)
        cram_index_free(fd);

    if (fd->own_pool && fd->pool)
        hts_tpool_destroy(fd->pool);

    if (fd->idxfp)
        if (bgzf_close(fd->idxfp) < 0)
            return -1;

    free(fd);
    return 0;
}

/*
 * Returns 1 if we hit an EOF while reading.
 */
int cram_eof(cram_fd *fd) {
    return fd->eof;
}


/*
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    va_start(args, opt);
    r = cram_set_voption(fd, opt, args);
    va_end(args);

    return r;
}

/*
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_voption(cram_fd *fd, enum hts_fmt_option opt, va_list args) {
    refs_t *refs;

    if (!fd) {
        errno = EBADF;
        return -1;
    }

    switch (opt) {
    case CRAM_OPT_DECODE_MD:
        fd->decode_md = va_arg(args, int);
        break;

    case CRAM_OPT_PREFIX:
        if (fd->prefix)
            free(fd->prefix);
        if (!(fd->prefix = strdup(va_arg(args, char *))))
            return -1;
        break;

    case CRAM_OPT_VERBOSITY:
        break;

    case CRAM_OPT_SEQS_PER_SLICE:
        fd->seqs_per_slice = va_arg(args, int);
        if (fd->bases_per_slice == BASES_PER_SLICE)
            fd->bases_per_slice = fd->seqs_per_slice * 500;
        break;

    case CRAM_OPT_BASES_PER_SLICE:
        fd->bases_per_slice = va_arg(args, int);
        break;

    case CRAM_OPT_SLICES_PER_CONTAINER:
        fd->slices_per_container = va_arg(args, int);
        break;

    case CRAM_OPT_EMBED_REF:
        fd->embed_ref = va_arg(args, int);
        break;

    case CRAM_OPT_NO_REF:
        fd->no_ref = va_arg(args, int);
        break;

    case CRAM_OPT_POS_DELTA:
        fd->ap_delta = va_arg(args, int);
        break;

    case CRAM_OPT_IGNORE_MD5:
        fd->ignore_md5 = va_arg(args, int);
        break;

    case CRAM_OPT_LOSSY_NAMES:
        fd->lossy_read_names = va_arg(args, int);
        // Currently lossy read names required paired (attached) reads.
        // TLEN 0 or being 1 out causes read pairs to be detached, breaking
        // the lossy read name compression, so we have extra options to
        // slacken the exact TLEN round-trip checks.
        fd->tlen_approx = fd->lossy_read_names;
        fd->tlen_zero = fd->lossy_read_names;
        break;

    case CRAM_OPT_USE_BZIP2:
        fd->use_bz2 = va_arg(args, int);
        break;

    case CRAM_OPT_USE_RANS:
        fd->use_rans = va_arg(args, int);
        break;

    case CRAM_OPT_USE_TOK:
        fd->use_tok = va_arg(args, int);
        break;

    case CRAM_OPT_USE_FQZ:
        fd->use_fqz = va_arg(args, int);
        break;

    case CRAM_OPT_USE_ARITH:
        fd->use_arith = va_arg(args, int);
        break;

    case CRAM_OPT_USE_LZMA:
        fd->use_lzma = va_arg(args, int);
        break;

    case CRAM_OPT_SHARED_REF:
        fd->shared_ref = 1;
        refs = va_arg(args, refs_t *);
        if (refs != fd->refs) {
            if (fd->refs)
                refs_free(fd->refs);
            fd->refs = refs;
            fd->refs->count++;
        }
        break;

    case CRAM_OPT_RANGE: {
        int r = cram_seek_to_refpos(fd, va_arg(args, cram_range *));
        pthread_mutex_lock(&fd->range_lock);
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        pthread_mutex_unlock(&fd->range_lock);
        return r;
    }

    case CRAM_OPT_RANGE_NOSEEK: {
        // As per CRAM_OPT_RANGE, but no seeking
        pthread_mutex_lock(&fd->range_lock);
        cram_range *r = va_arg(args, cram_range *);
        fd->range = *r;
        if (r->refid == HTS_IDX_NOCOOR) {
            fd->range.refid = -1;
            fd->range.start = 0;
        } else if (r->refid == HTS_IDX_START || r->refid == HTS_IDX_REST) {
            fd->range.refid = -2; // special case in cram_next_slice
        }
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        fd->ooc = 0;
        fd->eof = 0;
        pthread_mutex_unlock(&fd->range_lock);
        return 0;
    }

    case CRAM_OPT_REFERENCE:
        return cram_load_reference(fd, va_arg(args, char *));

    case CRAM_OPT_VERSION: {
        int major, minor;
        char *s = va_arg(args, char *);
        if (2 != sscanf(s, "%d.%d", &major, &minor)) {
            hts_log_error("Malformed version string %s", s);
            return -1;
        }
        if (!((major == 1 &&  minor == 0) ||
              (major == 2 && (minor == 0 || minor == 1)) ||
              (major == 3 && (minor == 0 || minor == 1)) ||
              (major == 4 &&  minor == 0))) {
            hts_log_error("Unknown version string; use 1.0, 2.0, 2.1, 3.0, 3.1 or 4.0");
            errno = EINVAL;
            return -1;
        }

        if (major > 3 || (major == 3 && minor > 0)) {
            hts_log_warning(
                "CRAM version %s is still a draft and subject to change.\n"
                "This is a technology demonstration that should not be "
                "used for archival data.", s);
        }

        fd->version = major*256 + minor;

        fd->use_rans = (CRAM_MAJOR_VERS(fd->version) >= 3) ? 1 : 0;

        fd->use_tok = ((CRAM_MAJOR_VERS(fd->version) == 3 &&
                        CRAM_MINOR_VERS(fd->version) >= 1) ||
                        CRAM_MAJOR_VERS(fd->version) >= 4) ? 1 : 0;
        cram_init_tables(fd);

        break;
    }

    case CRAM_OPT_MULTI_SEQ_PER_SLICE:
        fd->multi_seq_user = fd->multi_seq = va_arg(args, int);
        break;

    case CRAM_OPT_NTHREADS: {
        int nthreads =  va_arg(args, int);
        if (nthreads >= 1) {
            if (!(fd->pool = hts_tpool_init(nthreads)))
                return -1;

            fd->rqueue = hts_tpool_process_init(fd->pool, nthreads*2, 0);
            pthread_mutex_init(&fd->metrics_lock, NULL);
            pthread_mutex_init(&fd->ref_lock, NULL);
            pthread_mutex_init(&fd->range_lock, NULL);
            pthread_mutex_init(&fd->bam_list_lock, NULL);
            fd->shared_ref = 1;
            fd->own_pool = 1;
        }
        break;
    }

    case CRAM_OPT_THREAD_POOL: {
        htsThreadPool *p = va_arg(args, htsThreadPool *);
        fd->pool = p ? p->pool : NULL;
        if (fd->pool) {
            fd->rqueue = hts_tpool_process_init(fd->pool,
                                                p->qsize ? p->qsize : hts_tpool_size(fd->pool)*2,
                                                0);
            pthread_mutex_init(&fd->metrics_lock, NULL);
            pthread_mutex_init(&fd->ref_lock, NULL);
            pthread_mutex_init(&fd->range_lock, NULL);
            pthread_mutex_init(&fd->bam_list_lock, NULL);
        }
        fd->shared_ref = 1; // Needed to avoid clobbering ref between threads
        fd->own_pool = 0;

        //fd->qsize = 1;
        //fd->decoded = calloc(fd->qsize, sizeof(cram_container *));
        //hts_tpool_dispatch(fd->pool, cram_decoder_thread, fd);
        break;
    }

    case CRAM_OPT_REQUIRED_FIELDS:
        fd->required_fields = va_arg(args, int);
        if (fd->range.refid != -2)
            fd->required_fields |= SAM_POS;
        break;

    case CRAM_OPT_STORE_MD:
        fd->store_md = va_arg(args, int);
        break;

    case CRAM_OPT_STORE_NM:
        fd->store_nm = va_arg(args, int);
        break;

    case HTS_OPT_COMPRESSION_LEVEL:
        fd->level = va_arg(args, int);
        break;

    case HTS_OPT_PROFILE: {
        enum hts_profile_option prof = va_arg(args, int);
        switch (prof) {
        case HTS_PROFILE_FAST:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 1;
            fd->use_tok = 0;
            fd->seqs_per_slice = 10000;
            break;

        case HTS_PROFILE_NORMAL:
            break;

        case HTS_PROFILE_SMALL:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 6;
            fd->use_bz2 = 1;
            fd->use_fqz = 1;
            fd->seqs_per_slice = 25000;
            break;

        case HTS_PROFILE_ARCHIVE:
            if (fd->level == CRAM_DEFAULT_LEVEL) fd->level = 7;
            fd->use_bz2 = 1;
            fd->use_fqz = 1;
            fd->use_arith = 1;
            if (fd->level > 7)
                fd->use_lzma = 1;
            fd->seqs_per_slice = 100000;
            break;
        }

        if (fd->bases_per_slice == BASES_PER_SLICE)
            fd->bases_per_slice = fd->seqs_per_slice * 500;
        break;
    }

    default:
        hts_log_error("Unknown CRAM option code %d", opt);
        errno = EINVAL;
        return -1;
    }

    return 0;
}

int cram_check_EOF(cram_fd *fd)
{
    // Byte 9 in these templates is & with 0x0f to resolve differences
    // between ITF-8 interpretations between early Java and C
    // implementations of CRAM
    static const unsigned char TEMPLATE_2_1[30] = {
        0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
        0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
        0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00
    };
    static const unsigned char TEMPLATE_3[38] = {
        0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0,
        0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x05,
        0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00,
        0x01, 0x00, 0x01, 0x00, 0xee, 0x63, 0x01, 0x4b
    };

    unsigned char buf[38]; // max(sizeof TEMPLATE_*)

    uint8_t major = CRAM_MAJOR_VERS(fd->version);
    uint8_t minor = CRAM_MINOR_VERS(fd->version);

    const unsigned char *template;
    ssize_t template_len;
    if ((major < 2) ||
        (major == 2 && minor == 0)) {
        return 3; // No EOF support in cram versions less than 2.1
    } else if (major == 2 && minor == 1) {
        template = TEMPLATE_2_1;
        template_len = sizeof TEMPLATE_2_1;
    } else {
        template = TEMPLATE_3;
        template_len = sizeof TEMPLATE_3;
    }

    off_t offset = htell(fd->fp);
    if (hseek(fd->fp, -template_len, SEEK_END) < 0) {
        if (errno == ESPIPE) {
            hclearerr(fd->fp);
            return 2;
        }
        else {
            return -1;
        }
    }
    if (hread(fd->fp, buf, template_len) != template_len) return -1;
    if (hseek(fd->fp, offset, SEEK_SET) < 0) return -1;
    buf[8] &= 0x0f;
    return (memcmp(template, buf, template_len) == 0)? 1 : 0;
}
