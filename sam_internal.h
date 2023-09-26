/*  sam_internal.h -- internal functions; not part of the public API.

    Copyright (C) 2019-2020 Genome Research Ltd.

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
#ifdef __SSSE3__
#include "tmmintrin.h"
#endif
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
static inline void nibble2base(uint8_t *nib, char *seq, int len) {
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

    seq[0] = 0;
    const char *seq_end_ptr = seq + len;
    char *seq_cursor = seq;
    const uint8_t *nibble_cursor = nib;
    const char *seq_end_ptr_twoatatime = seq + (len & (~1ULL));
    #ifdef __SSSE3__
    const char *seq_vec_end_ptr = seq_end_ptr - (2 * sizeof(__m128i));
    __m128i first_upper_shuffle = _mm_setr_epi8(
        0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7, 0xff);
    __m128i first_lower_shuffle = _mm_setr_epi8(
        0xff, 0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7);
    __m128i second_upper_shuffle = _mm_setr_epi8(
        8, 0xff, 9, 0xff, 10, 0xff, 11, 0xff, 12, 0xff, 13, 0xff, 14, 0xff, 15, 0xff);
    __m128i second_lower_shuffle = _mm_setr_epi8(
        0xff, 8, 0xff, 9, 0xff, 10, 0xff, 11, 0xff, 12, 0xff, 13, 0xff, 14, 0xff, 15);
    __m128i nuc_lookup_vec = _mm_lddqu_si128((__m128i *)seq_nt16_str);
    /* Work on 16 encoded characters at the time resulting in 32 decoded characters
       Examples are given for 8 encoded characters A until H to keep it readable.
        Encoded stored as |AB|CD|EF|GH|
        Shuffle into |AB|00|CD|00|EF|00|GH|00| and
                     |00|AB|00|CD|00|EF|00|GH|
        Shift upper to the right resulting into
                     |0A|B0|0C|D0|0E|F0|0G|H0| and
                     |00|AB|00|CD|00|EF|00|GH|
        Merge with or resulting into (X stands for garbage)
                     |0A|XB|0C|XD|0E|XF|0G|XH|
        Bitwise and with 0b1111 leads to:
                     |0A|0B|0C|0D|0E|0F|0G|0H|
        We can use the resulting 4-bit integers as indexes for the shuffle of
        the nucleotide lookup. */
    while (seq_cursor < seq_vec_end_ptr) {
        __m128i encoded = _mm_lddqu_si128((__m128i *)nibble_cursor);

        __m128i first_upper = _mm_shuffle_epi8(encoded, first_upper_shuffle);
        __m128i first_lower = _mm_shuffle_epi8(encoded, first_lower_shuffle);
        __m128i shifted_first_upper = _mm_srli_epi64(first_upper, 4);
        __m128i first_merged = _mm_or_si128(shifted_first_upper, first_lower);
        __m128i first_indexes = _mm_and_si128(first_merged, _mm_set1_epi8(0b1111));
        __m128i first_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, first_indexes);
        _mm_storeu_si128((__m128i *)seq_cursor, first_nucleotides);

        __m128i second_upper = _mm_shuffle_epi8(encoded, second_upper_shuffle);
        __m128i second_lower = _mm_shuffle_epi8(encoded, second_lower_shuffle);
        __m128i shifted_second_upper = _mm_srli_epi64(second_upper, 4);
        __m128i second_merged = _mm_or_si128(shifted_second_upper, second_lower);
        __m128i second_indexes = _mm_and_si128(second_merged, _mm_set1_epi8(0b1111));
        __m128i second_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, second_indexes);
        _mm_storeu_si128((__m128i *)(seq_cursor + 16), second_nucleotides);

        nibble_cursor += sizeof(__m128i);
        seq_cursor += 2 * sizeof(__m128i);
    }
    #endif
    while (seq_cursor < seq_end_ptr_twoatatime) {
        // Note size_t cast helps gcc optimiser.
        memcpy(seq_cursor, code2base + ((size_t)*nibble_cursor * 2), 2);
        seq_cursor += 2;
        nibble_cursor += 1;
    }
    if (seq_cursor != seq_end_ptr) {
        /* There is a single encoded nuc left */
        uint8_t nibble_c = *nibble_cursor;
        uint8_t upper_nuc_index = nibble_c >> 4;
        seq_cursor[0] = seq_nt16_str[upper_nuc_index];
    }
}

#ifdef __cplusplus
}
#endif

#endif
