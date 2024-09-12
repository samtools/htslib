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

// These must be defined before the first system include to ensure that legacy
// BSD types needed by <sys/sysctl.h> remain defined when _XOPEN_SOURCE is set.
#if defined __APPLE__
#define _DARWIN_C_SOURCE
#elif defined __NetBSD__
#define _NETBSD_SOURCE
#endif

#include "htslib/sam.h"
#include "sam_internal.h"

#if defined __x86_64__
#include <immintrin.h>
#elif defined __ARM_NEON
#include <arm_neon.h>
#endif

#if defined __arm__ || defined __aarch64__

#if defined __linux__ || defined __FreeBSD__
#include <sys/auxv.h>
#elif defined __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#elif defined __NetBSD__
#include <stddef.h>
#include <sys/param.h>
#include <sys/sysctl.h>
#ifdef __aarch64__
#include <aarch64/armreg.h>
#else
#include <arm/armreg.h>
#endif
#elif defined _WIN32
#include <processthreadsapi.h>
#endif

static inline int cpu_supports_neon(void) {
#if defined __linux__ && defined __arm__ && defined HWCAP_NEON
    return (getauxval(AT_HWCAP) & HWCAP_NEON) != 0;
#elif defined __linux__ && defined __arm__ && defined HWCAP_ARM_NEON
    return (getauxval(AT_HWCAP) & HWCAP_ARM_NEON) != 0;
#elif defined __linux__ && defined __aarch64__ && defined HWCAP_ASIMD
    return (getauxval(AT_HWCAP) & HWCAP_ASIMD) != 0;
#elif defined __APPLE__ && defined __aarch64__
    int32_t ctl;
    size_t ctlsize = sizeof ctl;
    if (sysctlbyname("hw.optional.AdvSIMD", &ctl, &ctlsize, NULL, 0) != 0) return 0;
    if (ctlsize != sizeof ctl) return 0;
    return ctl;
#elif defined __FreeBSD__ && defined __arm__ && defined HWCAP_NEON
    unsigned long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_NEON) != 0;
#elif defined __FreeBSD__ && defined __aarch64__ && defined HWCAP_ASIMD
    unsigned long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_ASIMD) != 0;
#elif defined __NetBSD__ && defined __arm__ && defined ARM_MVFR0_ASIMD_MASK
    uint32_t buf[16];
    size_t buflen = sizeof buf;
    if (sysctlbyname("machdep.id_mvfr", buf, &buflen, NULL, 0) != 0) return 0;
    if (buflen < sizeof(uint32_t)) return 0;
    return (buf[0] & ARM_MVFR0_ASIMD_MASK) == 0x00000002;
#elif defined __NetBSD__ && defined __aarch64__ && defined ID_AA64PFR0_EL1_ADVSIMD
    struct aarch64_sysctl_cpu_id buf;
    size_t buflen = sizeof buf;
    if (sysctlbyname("machdep.cpu0.cpu_id", &buf, &buflen, NULL, 0) != 0) return 0;
    if (buflen < offsetof(struct aarch64_sysctl_cpu_id, ac_aa64pfr0) + sizeof(uint64_t)) return 0;
    return (buf.ac_aa64pfr0 & ID_AA64PFR0_EL1_ADVSIMD & 0x00e00000) == 0;
#elif defined _WIN32
    return IsProcessorFeaturePresent(PF_ARM_V8_INSTRUCTIONS_AVAILABLE) != 0;
#else
    return 0;
#endif
}

#endif

#ifdef BUILDING_SIMD_NIBBLE2BASE

void (*htslib_nibble2base)(uint8_t *nib, char *seq, int len) = nibble2base_default;

#if defined __x86_64__

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
    const char *seq_end_ptr = seq + len;
    char *seq_cursor = seq;
    uint8_t *nibble_cursor = nib;
    const char *seq_vec_end_ptr = seq_end_ptr - (2 * sizeof(__m128i) - 1);
    __m128i nuc_lookup_vec = _mm_lddqu_si128((__m128i *)seq_nt16_str);
    /* Nucleotides are encoded 4-bits per nucleotide and stored in 8-bit bytes
       as follows: |AB|CD|EF|GH|. The 4-bit codes (going from 0-15) can be used
       together with the pshufb instruction as a lookup table. The most efficient
       way is to use bitwise AND and shift to create two vectors. One with all
       the upper codes (|A|C|E|G|) and one with the lower codes (|B|D|F|H|).
       The lookup can then be performed and the resulting vectors can be
       interleaved again using the unpack instructions. */
    while (seq_cursor < seq_vec_end_ptr) {
        __m128i encoded = _mm_lddqu_si128((__m128i *)nibble_cursor);
        __m128i encoded_upper = _mm_srli_epi64(encoded, 4);
        encoded_upper = _mm_and_si128(encoded_upper, _mm_set1_epi8(15));
        __m128i encoded_lower = _mm_and_si128(encoded, _mm_set1_epi8(15));
        __m128i nucs_upper = _mm_shuffle_epi8(nuc_lookup_vec, encoded_upper);
        __m128i nucs_lower = _mm_shuffle_epi8(nuc_lookup_vec, encoded_lower);
        __m128i first_nucleotides = _mm_unpacklo_epi8(nucs_upper, nucs_lower);
        __m128i second_nucleotides = _mm_unpackhi_epi8(nucs_upper, nucs_lower);
        _mm_storeu_si128((__m128i *)seq_cursor, first_nucleotides);
        _mm_storeu_si128((__m128i *)(seq_cursor + sizeof(__m128i)),
                         second_nucleotides);
        nibble_cursor += sizeof(__m128i);
        seq_cursor += 2 * sizeof(__m128i);
    }
    nibble2base_default(nibble_cursor, seq_cursor, seq_end_ptr - seq_cursor);
}

__attribute__((constructor))
static void nibble2base_resolve(void) {
    if (__builtin_cpu_supports("ssse3")) {
        htslib_nibble2base = nibble2base_ssse3;
    }
}

#elif defined __ARM_NEON

static void nibble2base_neon(uint8_t *nib, char *seq0, int len) {
    uint8x16_t low_nibbles_mask = vdupq_n_u8(0x0f);
    uint8x16_t nuc_lookup_vec = vld1q_u8((const uint8_t *) seq_nt16_str);
#ifndef __aarch64__
    uint8x8x2_t nuc_lookup_vec2 = {{ vget_low_u8(nuc_lookup_vec), vget_high_u8(nuc_lookup_vec) }};
#endif

    uint8_t *seq = (uint8_t *) seq0;
    int blocks;

    for (blocks = len / 32; blocks > 0; --blocks) {
        uint8x16_t encoded = vld1q_u8(nib);
        nib += 16;

        /* Translate the high and low nibbles to nucleotide letters separately,
           then interleave them back together via vzipq for writing. */

        uint8x16_t high_nibbles = vshrq_n_u8(encoded, 4);
        uint8x16_t low_nibbles  = vandq_u8(encoded, low_nibbles_mask);

#ifdef __aarch64__
        uint8x16_t high_nucleotides = vqtbl1q_u8(nuc_lookup_vec, high_nibbles);
        uint8x16_t low_nucleotides  = vqtbl1q_u8(nuc_lookup_vec, low_nibbles);
#else
        uint8x8_t high_low  = vtbl2_u8(nuc_lookup_vec2, vget_low_u8(high_nibbles));
        uint8x8_t high_high = vtbl2_u8(nuc_lookup_vec2, vget_high_u8(high_nibbles));
        uint8x16_t high_nucleotides = vcombine_u8(high_low, high_high);

        uint8x8_t low_low  = vtbl2_u8(nuc_lookup_vec2, vget_low_u8(low_nibbles));
        uint8x8_t low_high = vtbl2_u8(nuc_lookup_vec2, vget_high_u8(low_nibbles));
        uint8x16_t low_nucleotides = vcombine_u8(low_low, low_high);
#endif

#ifdef __aarch64__
        vst1q_u8_x2(seq, vzipq_u8(high_nucleotides, low_nucleotides));
#else
        // Avoid vst1q_u8_x2 as GCC erroneously omits it on 32-bit ARM
        uint8x16x2_t nucleotides = {{ high_nucleotides, low_nucleotides }};
        vst2q_u8(seq, nucleotides);
#endif
        seq += 32;
    }

    if (len % 32 != 0)
        nibble2base_default(nib, (char *) seq, len % 32);
}

static __attribute__((constructor)) void nibble2base_resolve(void) {
    if (cpu_supports_neon()) htslib_nibble2base = nibble2base_neon;
}

#endif

#endif // BUILDING_SIMD_NIBBLE2BASE

// Potentially useful diagnostic, and prevents "empty translation unit" errors
const char htslib_simd[] =
    "SIMD functions present:"
#ifdef BUILDING_SIMD_NIBBLE2BASE
    " nibble2base"
#endif
    ".";
