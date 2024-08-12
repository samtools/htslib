/*  sam_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2019-2020, 2023-2024 Genome Research Ltd.

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

#ifndef HTSLIB_SAM_INTERNAL_H
#define HTSLIB_SAM_INTERNAL_H

#include <errno.h>
#include <stdint.h>

#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

// Used internally in the SAM format multi-threading.
int sam_state_destroy(samFile *fp);
int sam_set_thread_pool(htsFile *fp, htsThreadPool *p);
int sam_set_threads(htsFile *fp, int nthreads);

// Fastq state
int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...);
void fastq_state_destroy(samFile *fp);

// bam1_t data (re)allocation
int sam_realloc_bam_data(bam1_t *b, size_t desired);

static inline int realloc_bam_data(bam1_t *b, size_t desired)
{
    if (desired <= b->m_data) return 0;
    return sam_realloc_bam_data(b, desired);
}

static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    size_t new_len = (size_t) b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < bytes) { // Too big or overflow
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return sam_realloc_bam_data(b, new_len);
}

/*
 * Convert a nibble encoded BAM sequence to a string of bases.
 *
 * We do this 2 bp at a time for speed. Equiv to:
 *
 * for (i = 0; i < len; i++)
 *    seq[i] = seq_nt16_str[bam_seqi(nib, i)];
 */
static inline void nibble2base_default(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len/2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i*2], &code2base[(size_t)nib[i]*2], 2);

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

#if defined HAVE_ATTRIBUTE_CONSTRUCTOR && \
    ((defined __x86_64__ && defined HAVE_ATTRIBUTE_TARGET && defined HAVE_BUILTIN_CPU_SUPPORT_SSSE3) || \
     (defined __ARM_NEON))
#define BUILDING_SIMD_NIBBLE2BASE
#endif

static inline void nibble2base(uint8_t *nib, char *seq, int len) {
#ifdef BUILDING_SIMD_NIBBLE2BASE
    extern void (*htslib_nibble2base)(uint8_t *nib, char *seq, int len);
    htslib_nibble2base(nib, seq, len);
#else
    nibble2base_default(nib, seq, len);
#endif
}

#ifdef __cplusplus
}
#endif

#endif
