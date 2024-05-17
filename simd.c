/*  simd.c -- SIMD optimised versions of various internal functions.

    Copyright (C) 2024 Genome Research Ltd.

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

#include "htslib/sam.h"
#include "sam_internal.h"

#include "immintrin.h"

#ifdef BUILDING_SIMD_NIBBLE2BASE

/*
 * Convert a nibble encoded BAM sequence to a string of bases.
 *
 * Using SSSE3 instructions, 16 codepoints that hold 2 bases each can be
 * unpacked into 32 indexes from 0-15. Using the pshufb instruction these can
 * be converted to the IUPAC characters.
 * It falls back on the nibble2base_default function for the remainder.
 */

__attribute__((target("ssse3")))
static void nibble2base_ssse3(uint8_t *nib, char *seq, int len) {
    seq[0] = 0;
    const char *seq_end_ptr = seq + len;
    char *seq_cursor = seq;
    uint8_t *nibble_cursor = nib;
    const char *seq_vec_end_ptr = seq_end_ptr - (2 * sizeof(__m128i));
    __m128i first_upper_shuffle = _mm_setr_epi8(
        0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7, -1);
    __m128i first_lower_shuffle = _mm_setr_epi8(
        -1, 0, -1, 1, -1, 2, -1, 3, -1, 4, -1, 5, -1, 6, -1, 7);
    __m128i second_upper_shuffle = _mm_setr_epi8(
        8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15, -1);
    __m128i second_lower_shuffle = _mm_setr_epi8(
        -1, 8, -1, 9, -1, 10, -1, 11, -1, 12, -1, 13, -1, 14, -1, 15);
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
        __m128i first_indexes = _mm_and_si128(first_merged, _mm_set1_epi8(15));
        __m128i first_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, first_indexes);
        _mm_storeu_si128((__m128i *)seq_cursor, first_nucleotides);

        __m128i second_upper = _mm_shuffle_epi8(encoded, second_upper_shuffle);
        __m128i second_lower = _mm_shuffle_epi8(encoded, second_lower_shuffle);
        __m128i shifted_second_upper = _mm_srli_epi64(second_upper, 4);
        __m128i second_merged = _mm_or_si128(shifted_second_upper, second_lower);
        __m128i second_indexes = _mm_and_si128(second_merged, _mm_set1_epi8(15));
        __m128i second_nucleotides = _mm_shuffle_epi8(nuc_lookup_vec, second_indexes);
        _mm_storeu_si128((__m128i *)(seq_cursor + 16), second_nucleotides);

        nibble_cursor += sizeof(__m128i);
        seq_cursor += 2 * sizeof(__m128i);
    }
    nibble2base_default(nibble_cursor, seq_cursor, seq_end_ptr - seq_cursor);
}

void (*htslib_nibble2base)(uint8_t *nib, char *seq, int len) = nibble2base_default;

__attribute__((constructor))
static void nibble2base_resolve(void) {
    if (__builtin_cpu_supports("ssse3")) {
        htslib_nibble2base = nibble2base_ssse3;
    }
}

#endif // BUILDING_SIMD_NIBBLE2BASE

// Potentially useful diagnostic, and prevents "empty translation unit" errors
const char htslib_simd[] =
    "SIMD functions present:"
#ifdef BUILDING_SIMD_NIBBLE2BASE
    " nibble2base"
#endif
    ".";