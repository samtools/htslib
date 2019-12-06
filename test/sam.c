/*  test/sam.c -- SAM/BAM/CRAM API test cases.

    Copyright (C) 2014-2019 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>

// Suppress message for faidx_fetch_nseq(), which we're intentionally testing
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/hts_log.h"

KHASH_SET_INIT_STR(keep)
typedef khash_t(keep) *keephash_t;

#ifndef HTS_VERSION
#error HTS_VERSION not defined
#endif
#if HTS_VERSION < 100900
#error HTS_VERSION comparison incorrect
#endif

int status;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) fail(const char *fmt, ...)
{
    va_list args;

    fprintf(stderr, "Failed: ");
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");

    status = EXIT_FAILURE;
}

uint8_t *check_bam_aux_get(const bam1_t *aln, const char *tag, char type)
{
    uint8_t *p = bam_aux_get(aln, tag);
    if (p) {
        if (*p == type) return p;
        else fail("%s field of type '%c', expected '%c'", tag, *p, type);
    }
    else fail("can't find %s field", tag);

    return NULL;
}

static void check_int_B_array(bam1_t *aln, char *tag,
                             uint32_t nvals, int64_t *vals) {
    uint8_t *p;
    if ((p = check_bam_aux_get(aln, tag, 'B')) != NULL) {
        uint32_t i;

        if (bam_auxB_len(p) != nvals)
            fail("Wrong length reported for %s field, got %u, expected %u",
                 tag, bam_auxB_len(p), nvals);

        for (i = 0; i < nvals; i++) {
            if (bam_auxB2i(p, i) != vals[i]) {
                fail("Wrong value from bam_auxB2i for %s field index %u, "
                     "got %"PRId64" expected %"PRId64,
                     tag, i, bam_auxB2i(p, i), vals[i]);
            }
            if (bam_auxB2f(p, i) != (double) vals[i]) {
                fail("Wrong value from bam_auxB2f for %s field index %u, "
                     "got %f expected %f",
                     tag, i, bam_auxB2f(p, i), (double) vals[i]);
            }
        }
    }
}

#define PI 3.141592653589793
#define E  2.718281828459045
#define HELLO "Hello, world!"
#define NEW_HELLO "Yo, dude"
#define BEEF "DEADBEEF"

#define str(x) #x
#define xstr(x) str(x)

#define NELE(x) (sizeof(x)/sizeof(x[0]))

static int test_update_int(bam1_t *aln,
                           const char target_id[2], int64_t target_val,
                           char expected_type,
                           const char next_id[2], int64_t next_val,
                           char next_type) {
    uint8_t *p;

    // Try updating target
    if (bam_aux_update_int(aln, target_id, target_val) < 0) {
        fail("update %.2s tag", target_id);
        return -1;
    }

    // Check it's there and has the right type and value
    p = bam_aux_get(aln, target_id);
    if (!p) {
        fail("find  %.2s tag", target_id);
        return -1;
    }
    if (*p != expected_type || bam_aux2i(p) != target_val) {
        fail("%.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, *p, bam_aux2i(p), expected_type, target_val);
        return -1;
    }

    // If given, check that the next tag hasn't been clobbered by the
    // update above.
    if (!*next_id) return 0;
    p = bam_aux_get(aln, next_id);
    if (!p) {
        fail("find  %.2s tag after updating %.2s", next_id, target_id);
        return -1;
    }
    if (*p != next_type || bam_aux2i(p) != next_val) {
        fail("after updating %.2s to %"PRId64":"
             " %.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, target_val,
             next_id, *p, bam_aux2i(p), next_type, next_val);
        return -1;
    }
    return 0;
}

#define CHECK_ARRAY_VALS(T, GET_VAL, FMT1, FMT2) do {                   \
        T * vals = (T *) data;                                          \
    uint32_t i;                                                         \
    for (i = 0; i < nitems; i++) {                                      \
        if (GET_VAL(p, i) != vals[i]) {                                 \
            fail("Wrong value from %s for %.2s field index %u, "        \
                 "got %" FMT1 " expected %" FMT2,                       \
                 xstr(GET_VAL), target_id, i, GET_VAL(p, i), vals[i]);  \
            return -1;                                                  \
        }                                                               \
    }                                                                   \
} while (0)

static int test_update_array(bam1_t *aln, const char target_id[2],
                             uint8_t type, uint32_t nitems, void *data,
                             const char next_id[2], int64_t next_val,
                             char next_type)
{
    uint8_t *p;

    // Try updating target
    if (bam_aux_update_array(aln, target_id, type, nitems, data) < 0) {
        fail("update %2.s tag", target_id);
        return -1;
    }

    // Check values
    p = bam_aux_get(aln, target_id);
    if (!p) {
        fail("find  %.2s tag", target_id);
        return -1;
    }
    switch (type) {
        case 'c':
            CHECK_ARRAY_VALS(int8_t, bam_auxB2i, PRId64, PRId8); break;
        case 'C':
            CHECK_ARRAY_VALS(uint8_t, bam_auxB2i, PRId64, PRIu8); break;
        case 's':
            CHECK_ARRAY_VALS(int16_t, bam_auxB2i, PRId64, PRId16); break;
        case 'S':
            CHECK_ARRAY_VALS(uint16_t, bam_auxB2i, PRId64, PRIu16); break;
        case 'i':
            CHECK_ARRAY_VALS(int32_t, bam_auxB2i, PRId64, PRId32); break;
        case 'I':
            CHECK_ARRAY_VALS(uint32_t, bam_auxB2i, PRId64, PRIu32); break;
        case 'f':
            CHECK_ARRAY_VALS(float, bam_auxB2f, "e", "e"); break;
    }

    // If given, check that the next tag hasn't been clobbered by the
    // update above.
    if (!*next_id) return 0;
    p = bam_aux_get(aln, next_id);
    if (!p) {
        fail("find  %.2s tag after updating %.2s", next_id, target_id);
        return -1;
    }
    if (*p != next_type || bam_aux2i(p) != next_val) {
        fail("after updating %.2s:"
             " %.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, next_id, *p, bam_aux2i(p), next_type, next_val);
        return -1;
    }

    return 0;
}

// This function uses bam_hdr_t etc as a check ensuring the legacy typedef
// and functions continue to compile successfully.
static int aux_fields1(void)
{
    static const char sam[] = "data:,"
"@SQ\tSN:one\tLN:1000\n"
"@SQ\tSN:two\tLN:500\n"
"r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:" xstr(PI) "\tXd:d:" xstr(E) "\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,+2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:d:2.46801\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295\n"
"r2\t0x8D\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq\n"
;

    // Canonical form of the alignment records above, as output by sam_format1()
    static const char r1[] = "r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXi:i:37\tXf:f:3.14159\tXd:d:2.71828\tXZ:Z:" NEW_HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:f:9.8765\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295\tN0:i:-1234\tN1:i:1234\tN2:i:-2\tN3:i:3\tF1:f:4.5678\tN4:B:S,65535,32768,1,0\tN5:i:4242";
    static const char r2[] = "r2\t141\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    uint8_t *p;
    kstring_t ks = { 0, 0, NULL };
    int64_t b0vals[5] = { -2147483648LL,-1,0,1,2147483647LL }; // i
    int64_t b1vals[4] = { 0,1,2147483648LL,4294967295LL };     // I
    int64_t b2vals[5] = { -32768,-1,0,1,32767 };           // s
    int64_t b3vals[4] = { 0,1,32768,65535 };               // S
    int64_t b4vals[5] = { -128,-1,0,1,127 };               // c
    int64_t b5vals[4] = { 0,1,127,255 };                   // C
    // NB: Floats not doubles below!
    // See https://randomascii.wordpress.com/2012/06/26/doubles-are-not-floats-so-dont-compare-them/
    float bfvals[2] = { -3.14159f, 2.71828f };

    int8_t n4v1[] = { -128, -64, -32, -16, -8, -4, -2, -1,
                      0, 1, 2, 4, 8, 16, 32, 64, 127 };
    uint32_t n4v2[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1234, 5678, 1U << 31, 0 };
    int16_t n4v3[] = { -32768, -1, 0, 1, 32767 };
    float n4v4[] = { 0, 1, 2, 10, 20, 30, 1.5, -2.5 };
    uint8_t n4v5[] = { 0, 255 };
    int32_t n4v6[] = { -2147483647 - 1, 10, -1, 0, 1, 2147483647 };
    uint16_t n4v7[] = { 65535, 32768, 1, 0 };

    int32_t ival = -1234;
    uint32_t uval = 1234;
    float f1 = 4.5678;
    float f2 = 9.8765;

    size_t nvals, i;

    if (sam_read1(in, header, aln) >= 0) {
        if ((p = check_bam_aux_get(aln, "XA", 'A')) && bam_aux2A(p) != 'k')
            fail("XA field is '%c', expected 'k'", bam_aux2A(p));

        bam_aux_del(aln,p);
        if (bam_aux_get(aln,"XA"))
            fail("XA field was not deleted");

        if ((p = check_bam_aux_get(aln, "Xi", 'C')) && bam_aux2i(p) != 37)
            fail("Xi field is %"PRId64", expected 37", bam_aux2i(p));

        if ((p = check_bam_aux_get(aln, "Xf", 'f')) && fabs(bam_aux2f(p) - PI) > 1E-6)
            fail("Xf field is %.12f, expected pi", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "Xd", 'd')) && fabs(bam_aux2f(p) - E) > 1E-6)
            fail("Xf field is %.12f, expected e", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "XZ", 'Z')) && strcmp(bam_aux2Z(p), HELLO) != 0)
            fail("XZ field is \"%s\", expected \"%s\"", bam_aux2Z(p), HELLO);

        bam_aux_update_str(aln,"XZ",strlen(NEW_HELLO)+1,NEW_HELLO);
        if ((p = check_bam_aux_get(aln, "XZ", 'Z')) && strcmp(bam_aux2Z(p), NEW_HELLO) != 0)
            fail("XZ field is \"%s\", expected \"%s\"", bam_aux2Z(p), NEW_HELLO);


        if ((p = check_bam_aux_get(aln, "XH", 'H')) && strcmp(bam_aux2Z(p), BEEF) != 0)
            fail("XH field is \"%s\", expected \"%s\"", bam_aux2Z(p), BEEF);

        if ((p = check_bam_aux_get(aln, "XB", 'B'))
            && ! (memcmp(p, "Bc", 2) == 0
                  && memcmp(p + 2, "\x03\x00\x00\x00\xfe\x00\x02", 7) == 0))
            fail("XB field is %c,..., expected c,-2,0,+2", p[1]);

        check_int_B_array(aln, "B0", NELE(b0vals), b0vals);
        check_int_B_array(aln, "B1", NELE(b1vals), b1vals);
        check_int_B_array(aln, "B2", NELE(b2vals), b2vals);
        check_int_B_array(aln, "B3", NELE(b3vals), b3vals);
        check_int_B_array(aln, "B4", NELE(b4vals), b4vals);
        check_int_B_array(aln, "B5", NELE(b5vals), b5vals);

        nvals = NELE(bfvals);
        if ((p = check_bam_aux_get(aln, "Bf", 'B')) != NULL) {
            if (bam_auxB_len(p) != nvals)
                fail("Wrong length reported for Bf field, got %d, expected %zd",
                     bam_auxB_len(p), nvals);

            for (i = 0; i < nvals; i++) {
                if (bam_auxB2f(p, i) != bfvals[i]) {
                    fail("Wrong value from bam_auxB2f for Bf field index %zd, "
                         "got %f expected %f",
                         i, bam_auxB2f(p, i), bfvals[i]);
                }
            }
        }

        if ((p = check_bam_aux_get(aln, "ZZ", 'I')) && bam_aux2i(p) != 1000000)
            fail("ZZ field is %"PRId64", expected 1000000", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y1")) && bam_aux2i(p) != -2147483647-1)
            fail("Y1 field is %"PRId64", expected -2^31", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y2")) && bam_aux2i(p) != -2147483647)
            fail("Y2 field is %"PRId64", expected -2^31+1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y3")) && bam_aux2i(p) != -1)
            fail("Y3 field is %"PRId64", expected -1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y4")) && bam_aux2i(p) != 0)
            fail("Y4 field is %"PRId64", expected 0", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y5")) && bam_aux2i(p) != 1)
            fail("Y5 field is %"PRId64", expected 1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y6")) && bam_aux2i(p) != 2147483647)
            fail("Y6 field is %"PRId64", expected 2^31-1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y7")) && bam_aux2i(p) != 2147483648LL)
            fail("Y7 field is %"PRId64", expected 2^31", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y8")) && bam_aux2i(p) != 4294967295LL)
            fail("Y8 field is %"PRId64", expected 2^32-1", bam_aux2i(p));

        // Try appending some new tags
        if (bam_aux_append(aln, "N0", 'i', sizeof(ival), (uint8_t *) &ival) != 0)
            fail("Failed to append N0:i tag");

        if ((p = bam_aux_get(aln, "N0")) && bam_aux2i(p) != ival)
            fail("N0 field is %"PRId64", expected %d", bam_aux2i(p), ival);

        if (bam_aux_append(aln, "N1", 'I', sizeof(uval), (uint8_t *) &uval) != 0)
            fail("failed to append N1:I tag");

        if ((p = bam_aux_get(aln, "N1")) && bam_aux2i(p) != uval)
            fail("N1 field is %"PRId64", expected %u", bam_aux2i(p), uval);

        // Append tags with bam_aux_update_int()
        if (bam_aux_update_int(aln, "N2", -2) < 0)
            fail("failed to append N2:c tag");

        if (bam_aux_update_int(aln, "N3", 3) < 0)
            fail("failed to append N3:C tag");

        p = bam_aux_get(aln, "N2");
        if (!p)
            fail("failed to retrieve N2 tag");
        else if (*p != 'c' || bam_aux2i(p) != -2)
            fail("N2 field is %c:%"PRId64", expected c:-2", *p, bam_aux2i(p));

        p = bam_aux_get(aln, "N3");
        if (!p)
            fail("failed to retrieve N3 tag");
        else if (*p != 'C' || bam_aux2i(p) != 3)
            fail("N3 field is %c:%"PRId64", expected C:3", *p, bam_aux2i(p));

        // Try changing values with bam_aux_update_int()
        i = test_update_int(aln, "N2", 2, 'C', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", 1234, 'S', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", -1, 's', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", 4294967295U, 'I', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", -2, 'i', "N3", 3, 'C');

        // Append a value with bam_aux_update_float()
        if (bam_aux_update_float(aln, "F1", f1) < 0)
            fail("append F1:f tag");

        p = bam_aux_get(aln, "F1");
        if (!p)
            fail("retrieve F1 tag");
        else if (*p != 'f' || bam_aux2f(p) != f1)
            fail("F1 field is %c:%e, expected f:%e", *p, bam_aux2f(p), f1);

        // Change a double tag to a float
        if (bam_aux_update_float(aln, "F2", f2) < 0)
            fail("update F2 tag");

        p = bam_aux_get(aln, "F2");
        if (!p)
            fail("retrieve F2 tag");
        else if (*p != 'f' || bam_aux2f(p) != f2)
            fail("F2 field is %c:%e, expected f:%e", *p, bam_aux2f(p), f2);

        // Check the next one is intact too
        p = bam_aux_get(aln, "Y1");
        if (!p)
            fail("retrieve Y1 tag");
        else if (*p != 'i' && bam_aux2i(p) != -2147483647-1)
            fail("Y1 field is %"PRId64", expected -2^31", bam_aux2i(p));

        // bam_aux_update_array tests
        // append a new array
        i = test_update_array(aln, "N4", 'c', NELE(n4v1), n4v1, "\0\0", 0, 0);

        // Add a sentinal to check resizes work
        if (i == 0) i = test_update_int(aln, "N5", 4242, 'S', "\0\0", 0, 0);

        // alter the array tag a few times
        if (i == 0)
            i = test_update_array(aln, "N4", 'I', NELE(n4v2), n4v2,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 's', NELE(n4v3), n4v3,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'f', NELE(n4v4), n4v4,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'c', NELE(n4v5), n4v5,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'i', NELE(n4v6), n4v6,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'S', NELE(n4v7), n4v7,
                                  "N5", 4242, 'S');

        if (sam_format1(header, aln, &ks) < 0)
            fail("can't format record");

        if (strcmp(ks.s, r1) != 0)
            fail("record formatted incorrectly: \"%s\"", ks.s);
    }
    else fail("can't read record");

    if (sam_read1(in, header, aln) >= 0) {
        if (sam_format1(header, aln, &ks) < 0)
            fail("can't format record r2");

        if (aln->core.flag != 0x8D)
            fail("r2 flag value is 0x%X, expected 0x8D", aln->core.flag);

        if (strcmp(ks.s, r2) != 0)
            fail("record r2 formatted incorrectly: \"%s\"", ks.s);
    }
    else fail("can't read record r2");

    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
    free(ks.s);

    return 1;
}

static void set_qname(void)
{
    static const char sam[] = "data:,"
"@SQ\tSN:one\tLN:1000\n"
"@SQ\tSN:two\tLN:500\n"
"r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:" xstr(PI) "\tXd:d:" xstr(E) "\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,+2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:d:2.46801\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295\n"
"r22\t0x8D\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq\n"
"r12345678\t0x8D\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq\n"
;

    // Canonical form of the alignment records above, as output by sam_format1()
    static const char r1[] = "r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:3.14159\tXd:d:2.71828\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:d:2.46801\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295";
    static const char r2[] = "r234\t141\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq";
    static const char r3[] = "xyz\t141\t*\t0\t0\t*\t*\t0\t0\tATGC\tqqqq";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    kstring_t ks = { 0, 0, NULL };

    if (sam_read1(in, header, aln) >= 0) {
        bam_set_qname(aln, "r1");
        if (sam_format1(header, aln, &ks) < 0) fail("can't format record");
        if (strcmp(ks.s, r1) != 0) fail("record formatted incorrectly:\nGot: \"%s\"\nExp: \"%s\"\n", ks.s, r1);
    }
    else fail("can't read record");

    if (sam_read1(in, header, aln) >= 0) {
        bam_set_qname(aln, "r234");
        if (sam_format1(header, aln, &ks) < 0) fail("can't format record");
        if (strcmp(ks.s, r2) != 0) fail("record formatted incorrectly:\nGot: \"%s\"\nExp: \"%s\"\n", ks.s, r2);
    }
    else fail("can't read record");

    if (sam_read1(in, header, aln) >= 0) {
        bam_set_qname(aln, "xyz");
        if (sam_format1(header, aln, &ks) < 0) fail("can't format record");
        if (strcmp(ks.s, r3) != 0) fail("record formatted incorrectly:\nGot: \"%s\"\nExp: \"%s\"\n", ks.s, r3);
    }
    else fail("can't read record");

    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
    free(ks.s);
}

static void iterators1(void)
{
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0));
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_NONE, 0, 0));
}

// This function uses bam_hdr_t etc as a check ensuring the legacy typedef
// and functions continue to compile successfully.
static void copy_check_alignment(const char *infname, const char *informat,
    const char *outfname, const char *outmode, const char *outref)
{
    samFile *in = sam_open(infname, "r");
    samFile *out = sam_open(outfname, outmode);
    bam1_t *aln = bam_init1();
    bam_hdr_t *header = NULL;
    int res;

    if (!in) {
        fail("couldn't open %s", infname);
        goto err;
    }
    if (!out) {
        fail("couldn't open %s with mode %s", outfname, outmode);
        goto err;
    }
    if (!aln) {
        fail("bam_init1() failed");
        goto err;
    }

    if (outref) {
        if (hts_set_opt(out, CRAM_OPT_REFERENCE, outref) < 0) {
            fail("setting reference %s for %s", outref, outfname);
            goto err;
        }
    }

    header = sam_hdr_read(in);
    if (!header) {
        fail("reading header from %s", infname);
        goto err;
    }
    if (sam_hdr_write(out, header) < 0) fail("writing headers to %s", outfname);

    while ((res = sam_read1(in, header, aln)) >= 0) {
        int mod4 = ((intptr_t) bam_get_cigar(aln)) % 4;
        if (mod4 != 0)
            fail("%s CIGAR not 4-byte aligned; offset is 4k+%d for \"%s\"",
                 informat, mod4, bam_get_qname(aln));

        if (sam_write1(out, header, aln) < 0) fail("writing to %s", outfname);
    }
    if (res < -1) {
        fail("failed to read alignment from %s", infname);
    }

 err:
    bam_destroy1(aln);
    aln = NULL;
    bam_hdr_destroy(header);
    header = NULL;
    if (in) sam_close(in);
    if (out) sam_close(out);
}

static int check_target_names(sam_hdr_t *header, int expected_n_targets,
                              const char **expected_targets,
                              const int   *expected_lengths) {
    int i;

    // Check consistency of target_names array
    if (!header->target_name) {
        fail("target_name is NULL");
        return -1;
    }
    if (!header->target_len) {
        fail("target_len is NULL");
        return -1;
    }
    if (header->n_targets != expected_n_targets) {
        fail("header->n_targets (%d) != expected_n_targets (%d)",
             header->n_targets, expected_n_targets);
        return -1;
    }
    for (i = 0; i < expected_n_targets; i++) {
        if (!header->target_name[i]
            || strcmp(header->target_name[i], expected_targets[i]) != 0) {
            fail("header->target_name[%d] (%s) != \"%s\"",
                 i, header->target_name[i] ? header->target_name[i] : "NULL",
                 expected_targets[i]);
            return -1;
        }
        if (header->target_len[i] != expected_lengths[i]) {
            fail("header->target_len[%d] (%d) != %d",
                 i, header->target_len[i], expected_lengths[i]);
            return -1;
        }
    }
    return 0;
}

static void use_header_api(void) {
    static const char header_text[] = "data:,"
            "@HD\tVN:1.4\tGO:group\tSS:coordinate:queryname\n"
            "@SQ\tSN:ref0\tLN:100\n"
            "@CO\tThis line below will be updated\n"
            "@SQ\tSN:ref1\tLN:5001\tM5:983dalu9ue2\n"
            "@SQ\tSN:ref1.5\tLN:5001\n"
            "@CO\tThis line is good\n"
            "@SQ\tSN:ref2\tLN:5002\n";

    static const char rg_line[] =
            { '@', 'R', 'G', '\t', 'I', 'D', ':', 'r', 'u', 'n', '1' };

    static const char expected[] =
            "@HD\tVN:1.5\tSO:coordinate\n"
            "@CO\tThis line below will be updated\n"
            "@SQ\tSN:ref1\tLN:5001\tM5:kja8u34a2q3\n"
            "@CO\tThis line is good\n"
            "@SQ\tSN:ref2\tLN:5002\n"
            "@SQ\tSN:ref3\tLN:5003\n"
            "@PG\tID:samtools\tPN:samtools\tVN:1.9\n"
            "@RG\tID:run1\n"
            "@RG\tID:run4\n";

    static const char *expected_targets[] = { "ref1", "ref2", "ref3" };
    static const int   expected_lengths[] = { 5001, 5002, 5003 };
    const int expected_n_targets = sizeof(expected_targets) / sizeof(char *);

    const char outfname[] = "test/sam_header.tmp.sam_";
    const char outmode[] = "w";
    FILE *inf = NULL;
    char buffer[sizeof(expected) + 1024];

    samFile *in = sam_open(header_text, "r");
    samFile *out = sam_open(outfname, outmode);
    sam_hdr_t *header = NULL;
    kstring_t ks = { 0, 0, NULL };
    size_t bytes;
    int r;
    const char *name;

    if (!in) {
        fail("couldn't open file");
        goto err;
    }
    if (!out) {
        fail("couldn't open %s with mode %s", outfname, outmode);
        goto err;
    }

    header = sam_hdr_read(in);
    if (!header) {
        fail("reading header from file");
        goto err;
    }
    r = sam_hdr_remove_tag_id(header, "HD", NULL, NULL, "GO");
    if (r != 1) { fail("sam_hdr_remove_tag_id"); goto err; }

    r = sam_hdr_update_hd(header, "VN", "1.5");
    if (r != 0) { fail("sam_hdr_update_hd"); goto err; }

    r = sam_hdr_add_line(header, "SQ", "SN", "ref3", "LN", "5003", NULL);
    if (r < 0) { fail("sam_hdr_add_line"); goto err; }

    r = sam_hdr_update_line(header, "SQ", "SN", "ref1",
                             "M5", "kja8u34a2q3", NULL);
    if (r != 0) { fail("sam_hdr_update_line SQ"); goto err; }

    r = sam_hdr_add_pg(header, "samtools", "VN", "1.9", NULL);
    if (r != 0) { fail("sam_hdr_add_pg"); goto err; }

    // Test addition with no newline or trailing NUL
    r = sam_hdr_add_lines(header, rg_line, sizeof(rg_line));
    if (r != 0) { fail("sam_hdr_add_lines rg_line"); goto err; }

    // Test header line removal
    r = sam_hdr_add_line(header, "RG", "ID", "run2", NULL);
    if (r < 0) { fail("sam_hdr_add_line"); goto err; }

    r = sam_hdr_add_line(header, "RG", "ID", "run3", NULL);
    if (r < 0) { fail("sam_hdr_add_line"); goto err; }

    r = sam_hdr_add_line(header, "RG", "ID", "run4", NULL);
    if (r < 0) { fail("sam_hdr_add_line"); goto err; }

    r = sam_hdr_line_index(header, "RG", "run4");
    if (r != 3) { fail("sam_hdr_line_index - run4~3"); goto err; }

    r = sam_hdr_line_index(header, "RG", "run5");
    if (r != -1) { fail("sam_hdr_line_index - run5~-1"); goto err; }

    name = sam_hdr_line_name(header, "RG", 2);
    if (!name || strcmp(name, "run3")) { fail("sam_hdr_line_name - 2~run3"); goto err; }

    name = sam_hdr_line_name(header, "RG", 10);
    if (name) { fail("sam_hdr_line_name - 10~NULL"); goto err; }

    r = sam_hdr_remove_line_id(header, "RG", "ID", "run2");
    if (r < 0) { fail("sam_hdr_remove_line_id"); goto err; }

    r = sam_hdr_find_tag_id(header, "RG", "ID", "run3", "ID", &ks);
    if (r < 0 || !ks.s || strcmp(ks.s, "run3") != 0) {
        fail("sam_hdr_find_tag_id() expected \"run3\" got \"%s\"",
             r == 0 && ks.s ? ks.s : "NULL");
        goto err;
    }

    r = sam_hdr_remove_line_pos(header, "RG", 1); // Removes run3
    if (r < 0) { fail("sam_hdr_remove_line_pos"); goto err; }

    r = sam_hdr_remove_line_id(header, "SQ", "SN", "ref0");
    if (r < 0) { fail("sam_hdr_remove_line_id"); goto err; }

    r = sam_hdr_remove_line_pos(header, "SQ", 1); // Removes ref1.5
    if (r < 0) { fail("sam_hdr_remove_line_pos"); goto err; }

    r = sam_hdr_find_tag_id(header, "SQ", "SN", "ref1", "M5", &ks);
    if (r < 0 || !ks.s || strcmp(ks.s, "kja8u34a2q3") != 0) {
        fail("sam_hdr_find_tag_id() expected \"kja8u34a2q3\" got \"%s\"",
             r == 0 && ks.s ? ks.s : "NULL");
        goto err;
    }

    r = sam_hdr_line_index(header, "RG", "run4");
    if (r != 1) { fail("sam_hdr_line_index - run4~1"); goto err; }

    name = sam_hdr_line_name(header, "RG", 2);
    if (name) { fail("sam_hdr_line_name - 2~NULL"); goto err; }

    r = sam_hdr_remove_tag_hd(header, "SS");
    if (r < 0) {
        fail("sam_hdr_remove_tag_hd");
    }

    r = sam_hdr_find_hd(header, &ks);
    if (r < 0 || !ks.s || strcmp(ks.s, "@HD\tVN:1.5") != 0) {
        fail("sam_hdr_find_hd() expected \"@HD\tVN:1.5\" got \"%s\"",
             r == 0 && ks.s ? ks.s : "NULL");
    }

    r = sam_hdr_find_tag_hd(header, "VN", &ks);
    if (r < 0 || !ks.s || strcmp(ks.s, "1.5") != 0) {
        fail("sam_hdr_find_tag_hd() expected \"1.5\" got \"%s\"",
             r == 0 && ks.s ? ks.s : "NULL");
    }

    r = sam_hdr_update_hd(header, "SO", "coordinate");
    if (r < 0) {
        fail("sam_hdr_update_hd");
    }

    if (check_target_names(header, expected_n_targets, expected_targets,
                           expected_lengths) < 0) {
        goto err;
    }

    if ((r = sam_hdr_count_lines(header, "HD")) != 1) {
        fail("incorrect HD line count - expected 1, got %d", r);
        goto err;
    }
    if ((r = sam_hdr_count_lines(header, "SQ")) != 3) {
        fail("incorrect SQ line count - expected 3, got %d", r);
        goto err;
    }
    if ((r = sam_hdr_count_lines(header, "PG")) != 1) {
        fail("incorrect PG line count - expected 1, got %d", r);
        goto err;
    }
    if ((r = sam_hdr_count_lines(header, "RG")) != 2) {
        fail("incorrect RG line count - expected 2, got %d", r);
        goto err;
    }
    if ((r = sam_hdr_count_lines(header, "CO")) != 2) {
        fail("incorrect CO line count - expected 2, got %d", r);
        goto err;
    }

    if (sam_hdr_write(out, header) < 0) {
        fail("writing headers to \"%s\"", outfname);
        goto err;
    }
    r = sam_close(out);
    out = NULL;
    if (r < 0) {
        fail("close \"%s\"", outfname);
        goto err;
    }

    inf = fopen(outfname, "r");
    if (!inf) {
        fail("Opening written header \"%s\"", outfname);
        goto err;
    }
    bytes = fread(buffer, 1, sizeof(buffer), inf);
    if (bytes != sizeof(expected) - 1 || memcmp(buffer, expected, bytes) != 0) {
        fail("edited header does not match expected version");
        fprintf(stderr,
                "---------- Expected:\n%.*s\n"
                "++++++++++ Got:\n%.*s\n"
                "====================\n",
                (int) sizeof(expected), expected,
                (int) bytes, buffer);
        goto err;
    }

    free(ks_release(&ks));

 err:
    sam_hdr_destroy(header);
    header = NULL;
    if (in) sam_close(in);
    if (out) sam_close(out);
    if (inf) fclose(inf);
    free(ks_release(&ks));
}

static void test_header_pg_lines(void) {
    static const char header_text[] = "data:,"
        "@HD\tVN:1.5\n"
        "@PG\tID:prog1\tPN:prog1\n"
        "@PG\tID:prog2\tPN:prog2\tPP:prog1\n";

    static const char expected[] =
        "@HD\tVN:1.5\n"
        "@PG\tID:prog1\tPN:prog1\n"
        "@PG\tID:prog2\tPN:prog2\tPP:prog1\n"
        "@PG\tID:prog3\tPN:prog3\tPP:prog2\n"
        "@PG\tID:prog4\tPN:prog4\tPP:prog1\n"
        "@PG\tID:prog5\tPN:prog5\tPP:prog2\n"
        "@PG\tID:prog6\tPN:prog6\tPP:prog3\n"
        "@PG\tID:prog6.1\tPN:prog6\tPP:prog4\n"
        "@PG\tID:prog6.2\tPN:prog6\tPP:prog5\n"
        "@PG\tPN:prog7\tID:my_id\tPP:prog6\n";

    samFile *in = sam_open(header_text, "r");
    sam_hdr_t *header = NULL;
    const char *text = NULL;
    enum htsLogLevel old_log_level;
    int r;

    if (!in) {
        fail("couldn't open file");
        goto err;
    }

    header = sam_hdr_read(in);
    if (!header) {
        fail("reading header from file");
        goto err;
    }

    r = sam_hdr_add_pg(header, "prog3", NULL);
    if (r != 0) { fail("sam_hdr_add_pg prog3"); goto err; }


    r = sam_hdr_add_pg(header, "prog4", "PP", "prog1", NULL);
    if (r != 0) { fail("sam_hdr_add_pg prog4"); goto err; }

    r = sam_hdr_add_line(header, "PG", "ID",
                         "prog5", "PN", "prog5", "PP", "prog2", NULL);
    if (r != 0) { fail("sam_hdr_add_line @PG ID:prog5"); goto err; }

    r = sam_hdr_add_pg(header, "prog6", NULL);
    if (r != 0) { fail("sam_hdr_add_pg prog6"); goto err; }

    r = sam_hdr_add_pg(header, "prog7", "ID", "my_id", "PP", "prog6", NULL);
    if (r != 0) { fail("sam_hdr_add_pg prog7"); goto err; }

    text = sam_hdr_str(header);
    if (!text) { fail("sam_hdr_str"); goto err; }

    // These should fail
    old_log_level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_OFF);

    r = sam_hdr_add_pg(header, "prog8", "ID", "my_id", NULL);
    if (r == 0) { fail("sam_hdr_add_pg prog8 (unexpected success)"); goto err; }

    r = sam_hdr_add_pg(header, "prog9", "PP", "non-existent", NULL);
    if (r == 0) { fail("sam_hdr_add_pg prog9 (unexpected success)"); goto err; }

    hts_set_log_level(old_log_level);
    // End failing tests

    if (strcmp(text, expected) != 0) {
        fail("edited header does not match expected version");
        fprintf(stderr,
                "---------- Expected:\n%s\n"
                "++++++++++ Got:\n%s\n"
                "====================\n",
                expected, text);
        goto err;
    }

 err:
    sam_hdr_destroy(header);
    header = NULL;
    if (in) sam_close(in);
    return;
}

static void test_header_updates(void) {
    static const char header_text[] =
        "@HD\tVN:1.4\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@SQ\tSN:chr2\tLN:200\n"
        "@SQ\tSN:chr3\tLN:300\n"
        "@RG\tID:run1\n"
        "@RG\tID:run2\n"
        "@RG\tID:run3\n"
        "@PG\tID:prog1\tPN:prog1\n";

    static const char expected[] =
        "@HD\tVN:1.4\n"
        "@SQ\tSN:1\tLN:100\n"
        "@SQ\tSN:chr2\tLN:2000\n"
        "@SQ\tSN:chr3\tLN:300\n"
        "@RG\tID:run1\tDS:hello\n"
        "@RG\tID:aliquot2\n"
        "@RG\tID:run3\n"
        "@PG\tID:prog1\tPN:prog1\n";

    static const char *expected_targets[] = { "1", "chr2", "chr3" };
    static const int   expected_lengths[] = { 100, 2000, 300 };
    const int expected_n_targets = sizeof(expected_targets) / sizeof(char *);

    sam_hdr_t *header = sam_hdr_parse(sizeof(header_text) - 1, header_text);
    const char *hdr_str;
    int r, i, old_log_level;

    if (!header) {
        fail("creating sam header");
        goto err;
    }

    if (sam_hdr_name2tid(header, "chr1") != 0) { // Should now be unknown
        fail("sam_hdr_name2tid(\"chr1\") != 0");
        goto err;
    }

    r = sam_hdr_update_line(header, "SQ", "SN", "chr2", "LN", "2000", NULL);
    if (r != 0) { fail("sam_hdr_update_line SQ SN chr2 LN 2000"); goto err; }
    r = sam_hdr_update_line(header, "SQ", "SN", "chr1", "SN", "1", NULL);
    if (r != 0) { fail("sam_hdr_update_line SQ SN chr1 SN 1"); goto err; }
    r = sam_hdr_update_line(header, "RG", "ID", "run1", "DS", "hello", NULL);
    if (r != 0) { fail("sam_hdr_update_line RG ID run1 DS hello"); goto err; }
    r = sam_hdr_update_line(header, "RG", "ID", "run2", "ID", "aliquot2", NULL);
    if (r != 0) { fail("sam_hdr_update_line RG ID run2 ID aliquot2"); goto err; }

    // These should fail
    old_log_level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_OFF);

    r = sam_hdr_update_line(header, "PG", "ID", "prog1", "ID", "prog2", NULL);
    if (r == 0) { fail("sam_hdr_update_line PG ID prog1 ID prog2"); goto err; }

    r = sam_hdr_update_line(header, "SQ", "SN", "chr3", "SN", "chr2", NULL);
    if (r == 0) { fail("sam_hdr_update_line SQ SN chr3 SN chr2"); goto err; }

    r = sam_hdr_update_line(header, "RG", "ID", "run3", "ID", "run1", NULL);
    if (r == 0) { fail("sam_hdr_update_line RG ID run3 ID run1"); goto err; }

    hts_set_log_level(old_log_level);
    // End failing tests

    if (check_target_names(header, expected_n_targets, expected_targets,
                           expected_lengths) < 0) {
        goto err;
    }

    for (i = 0; i < expected_n_targets; i++) {
        if (sam_hdr_name2tid(header, expected_targets[i]) != i) {
            fail("sam_hdr_name2tid unexpected result");
            goto err;
        }
    }
    if (sam_hdr_name2tid(header, "chr1") != -1) { // Should now be unknown
        fail("sam_hdr_name2tid(\"chr1\") != -1");
        goto err;
    }

    hdr_str = sam_hdr_str(header);
    if (!hdr_str || strcmp(hdr_str, expected) != 0) {
        fail("edited header does not match expected version");
        fprintf(stderr,
                "---------- Expected:\n%s\n"
                "++++++++++ Got:\n%s\n"
                "====================\n",
                expected, hdr_str ? hdr_str : "<NULL>");
        goto err;
    }

 err:
    sam_hdr_destroy(header);
}

static void test_header_remove_lines(void) {
    static const char header_text[] =
        "@HD\tVN:1.4\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@SQ\tSN:chr2\tLN:200\n"
        "@SQ\tSN:chr3\tLN:300\n"
        "@RG\tID:run1\n"
        "@RG\tID:run2\n"
        "@RG\tID:run3\n"
        "@PG\tID:prog1\tPN:prog1\n";

    static const char expected[] =
        "@HD\tVN:1.4\n"
        "@SQ\tSN:chr1\tLN:100\n"
        "@SQ\tSN:chr3\tLN:300\n"
        "@PG\tID:prog1\tPN:prog1\n";

    sam_hdr_t *header = sam_hdr_parse(sizeof(header_text) - 1, header_text);
    keephash_t rh = kh_init(keep);
    khint_t k;
    const char *hdr_str;
    int r = 0;

    if (!header) {
        fail("creating sam header");
        goto err;
    }
    if (!rh) {
        fail("creating keep hash table");
        goto err;
    }

    kh_put(keep, rh, strdup("chr3"), &r);
    if (r < 0) { fail("adding chr3 to hash table"); goto err; }
    kh_put(keep, rh, strdup("chr1"), &r);
    if (r < 0) { fail("adding chr1 to hash table"); goto err; }

    r = sam_hdr_remove_lines(header, "SQ", "SN", rh);
    if (r != 0) { fail("sam_hdr_remove_lines SQ SN rh"); goto err; }

    r = sam_hdr_remove_lines(header, "RG", "ID", NULL);
    if (r != 0) { fail("sam_hdr_remove_lines RG ID NULL"); goto err; }

    hdr_str = sam_hdr_str(header);
    if (!hdr_str || strcmp(hdr_str, expected) != 0) {
        fail("edited header does not match expected version");
        fprintf(stderr,
                "---------- Expected:\n%s\n"
                "++++++++++ Got:\n%s\n"
                "====================\n",
                expected, hdr_str ? hdr_str : "<NULL>");
        goto err;
    }

 err:
    if (rh) {
        for (k = 0; k < kh_end(rh); ++k)
            if (kh_exist(rh, k)) free((char*)kh_key(rh, k));
        kh_destroy(keep, rh);
    }
    if (header) sam_hdr_destroy(header);
}

static void check_ref_lookup(sam_hdr_t *header, const char *msg, ...) {
    const char *name;
    va_list args;
    va_start(args, msg);
    while ((name = va_arg(args, const char *)) != NULL) {
        int exp = va_arg(args, int);
        int tid = sam_hdr_name2tid(header, name);
        if (tid != exp)
            fail("%s: altname \"%s\" => %d (expected %d)", msg, name, tid, exp);
    }
    va_end(args);
}

static void test_header_ref_altnames(void) {
    static const char initial_header[] =
        "@SQ\tSN:1\tLN:100\tAN:chr1\n"
        "@SQ\tSN:chr2\tAN:2\tLN:200\n"
        "@SQ\tSN:3\tLN:300\n"
        "@SQ\tSN:chrMT\tLN:16569\tAN:MT,chrM,M\n";

    sam_hdr_t *header = sam_hdr_init();
    if (header == NULL) { fail("sam_hdr_init"); return; }

    if (sam_hdr_add_lines(header, initial_header, 0) < 0)
        fail("sam_hdr_add_lines() for altnames");

    check_ref_lookup(header, "initial",
        "1", 0, "chr1", 0, "2", 1, "chr2", 1, "3", 2,
        "chrMT", 3, "chrM", 3, "M", 3, "fred", -1, "barney", -1,
        NULL);

    if (sam_hdr_add_line(header, "SQ", "AN", "fred", "LN", "500", "SN", "barney", NULL) < 0)
        fail("sam_hdr_add_line() for altnames");

    check_ref_lookup(header, "barney added",
        "1", 0, "chr1", 0, "2", 1, "chr2", 1, "3", 2,
        "chrMT", 3, "chrM", 3, "M", 3, "fred", 4, "barney", 4,
        NULL);

    if (sam_hdr_remove_line_id(header, "SQ", "SN", "chr2") < 0)
        fail("sam_hdr_remove_line_id() for altnames");

    check_ref_lookup(header, "chr2 removed",
        "1", 0, "chr1", 0, "2", -1, "chr2", -1, "3", 1,
        "chrMT", 2, "chrM", 2, "M", 2, "fred", 3, "barney", 3,
        NULL);

    if (sam_hdr_remove_tag_id(header, "SQ", "SN", "1", "AN") < 0)
        fail("sam_hdr_remove_tag_id() for altnames");

    check_ref_lookup(header, "1's AN removed",
        "1", 0, "chr1", -1, "CM000663", -1, "2", -1, "chr2", -1, "3", 1,
        "chrMT", 2, "chrM", 2, "M", 2, "fred", 3, "barney", 3,
        NULL);

    sam_hdr_destroy(header);

    static const char initial_header_duplicates[] =
        "@SQ\tSN:1\tLN:100\tAN:foo,2\n"
        "@SQ\tSN:2\tLN:200\tAN:bar\n"
        "@SQ\tSN:3\tLN:300\tAN:baz,3\n";

    header = sam_hdr_init();
    if (header == NULL) { fail("sam_hdr_init"); return; }

    int old_log_level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_ERROR); // Silence "Duplicate entry AN:2" warning

    if (sam_hdr_add_lines(header, initial_header_duplicates, 0) < 0)
        fail("sam_hdr_add_lines() for altnames with duplicates");

    hts_set_log_level(old_log_level);

    // Check "2" is SN:2 and not AN:2
    check_ref_lookup(header, "initial_header_duplicates",
                     "1", 0, "foo", 0,
                     "2", 1, "bar", 1,
                     "3", 2, "baz", 2, NULL);

    if (sam_hdr_remove_tag_id(header, "SQ", "SN", "1", "AN") < 0)
        fail("sam_hdr_remove_tag_id() for duplicate altnames SN:1");

    // Check "2" still works and "foo" does not
    check_ref_lookup(header, "initial_header_duplicates",
                     "1", 0, "foo", -1,
                     "2", 1, "bar", 1,
                     "3", 2, "baz", 2, NULL);

    if (sam_hdr_remove_tag_id(header, "SQ", "SN", "3", "AN") < 0)
        fail("sam_hdr_remove_tag_id() for duplicate altnames SN:3");

    // Check "3" still works and "baz" does not
    check_ref_lookup(header, "initial_header_duplicates",
                     "1", 0, "foo", -1,
                     "2", 1, "bar", 1,
                     "3", 2, "baz", -1, NULL);

    sam_hdr_destroy(header);
}

#define ABC50   "abcdefghijklmnopqrstuvwxyabcdefghijklmnopqrstuvwxy"
#define ABC250  ABC50 ABC50 ABC50 ABC50 ABC50

static void samrecord_layout(void)
{
    static const char qnames[] = "data:,"
"@SQ\tSN:CHROMOSOME_II\tLN:5000\n"
    "a\t0\tCHROMOSOME_II\t100\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
   "bc\t0\tCHROMOSOME_II\t200\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
  "def\t0\tCHROMOSOME_II\t300\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
 "ghij\t0\tCHROMOSOME_II\t400\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
"klmno\t0\tCHROMOSOME_II\t500\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
     ABC250 "\t0\tCHROMOSOME_II\t600\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
    ABC250 "1\t0\tCHROMOSOME_II\t650\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
   ABC250 "12\t0\tCHROMOSOME_II\t700\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
  ABC250 "123\t0\tCHROMOSOME_II\t750\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
 ABC250 "1234\t0\tCHROMOSOME_II\t800\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
;

    size_t bam1_t_size, bam1_t_size2;

    assert(sizeof(hts_pos_t) == 8 || sizeof(hts_pos_t) == 4);
    int core_size = sizeof(hts_pos_t) == 8 ? 48 : 36;
    bam1_t_size = (core_size + sizeof(int) + sizeof(char *) + sizeof(uint64_t)
                   + 2 * sizeof(uint32_t));
    bam1_t_size2 = bam1_t_size + 4;  // Account for padding on some platforms

    if (sizeof (bam1_core_t) != core_size)
        fail("sizeof bam1_core_t is %zu, expected %d",
             sizeof (bam1_core_t), core_size);

    if (sizeof (bam1_t) != bam1_t_size && sizeof (bam1_t) != bam1_t_size2)
        fail("sizeof bam1_t is %zu, expected either %zu or %zu",
             sizeof(bam1_t), bam1_t_size, bam1_t_size2);

    copy_check_alignment(qnames, "SAM",
                         "test/sam_alignment.tmp.bam", "wb", NULL);
    copy_check_alignment("test/sam_alignment.tmp.bam", "BAM",
                         "test/sam_alignment.tmp.cram", "wc", "test/ce.fa");
    copy_check_alignment("test/sam_alignment.tmp.cram", "CRAM",
                         "test/sam_alignment.tmp.sam_", "w", NULL);
}

static int check_ref_lengths(const sam_hdr_t *header,
                             const hts_pos_t *expected_lengths,
                             int num_refs, const char *hdr_name)
{
    int i;
    for (i = 0; i < num_refs; i++) {
        hts_pos_t ln = sam_hdr_tid2len(header, i);
        if (ln != expected_lengths[i]) {
            fail("Wrong length for %s ref %d : "
                 "expected %"PRIhts_pos" got %"PRIhts_pos"\n",
                 hdr_name, i, expected_lengths[i], ln);
            return -1;
        }
    }
    return 0;
}

static void check_big_ref(int parse_header)
{
    static const char sam_text[] = "data:,"
        "@HD\tVN:1.4\n"
        "@SQ\tSN:large#1\tLN:5000000000\n"
        "@SQ\tSN:small#1\tLN:100\n"
        "@SQ\tSN:large#2\tLN:9223372034707292158\n"
        "@SQ\tSN:small#2\tLN:1\n"
        "r1\t0\tlarge#1\t4999999000\t50\t8M\t*\t0\t0\tACGTACGT\tabcdefgh\n"
        "r2\t0\tsmall#1\t1\t50\t8M\t*\t0\t0\tACGTACGT\tabcdefgh\n"
        "r3\t0\tlarge#2\t9223372034707292000\t50\t8M\t*\t0\t0\tACGTACGT\tabcdefgh\n"
        "p1\t99\tlarge#2\t1\t50\t8M\t=\t9223372034707292150\t9223372034707292158\tACGTACGT\tabcdefgh\n"
        "p1\t147\tlarge#2\t9223372034707292150\t50\t8M\t=\t1\t-9223372034707292158\tACGTACGT\tabcdefgh\n"
        "r4\t0\tsmall#2\t2\t50\t8M\t*\t0\t0\tACGTACGT\tabcdefgh\n";
    const hts_pos_t expected_lengths[] = {
        5000000000LL, 100LL, 9223372034707292158LL, 1LL
    };
    const int expected_tids[] = {
        0, 1, 2, 2, 2, 3
    };
    const int expected_mtid[] = {
        -1, -1, -1, 2, 2, -1
    };
    const hts_pos_t expected_positions[] = {
        4999999000LL - 1, 1LL - 1, 9223372034707292000LL - 1, 1LL - 1,
        9223372034707292150LL - 1, 2LL - 1
    };
    const hts_pos_t expected_mpos[] = {
        -1, -1, -1, 9223372034707292150LL - 1, 1LL - 1, -1
    };
    samFile *in = NULL, *out = NULL;
    sam_hdr_t *header = NULL, *dup_header = NULL;
    bam1_t *aln = bam_init1();
    const int num_refs = sizeof(expected_lengths) / sizeof(expected_lengths[0]);
    const int num_align = sizeof(expected_tids) / sizeof(expected_tids[0]);
    const char *outfname = "test/sam_big_ref.tmp.sam_";
    int i, r;
    char buffer[sizeof(sam_text) + 1024];
    FILE *inf = NULL;
    size_t bytes;

    if (!aln) {
        fail("Out of memory");
        goto cleanup;
    }

    in = sam_open(sam_text, "r");
    if (!in) {
        fail("Opening SAM file");
        goto cleanup;
    }
    out = sam_open(outfname, "w");
    if (!out) {
        fail("Opening output SAM file \"%s\"", outfname);
        goto cleanup;
    }
    header = sam_hdr_read(in);
    if (!header) {
        fail("Reading SAM header");
        goto cleanup;
    }
    if (parse_header) {
        // This will force the header to be parsed
        if (sam_hdr_count_lines(header, "SQ") != num_refs) {
            fail("Wrong number of SQ lines in header");
            goto cleanup;
        }
    }
    if (check_ref_lengths(header, expected_lengths, num_refs, "header") < 0)
        goto cleanup;

    dup_header = sam_hdr_dup(header);
    if (!dup_header) {
        fail("Failed to duplicate header");
    }

    if (check_ref_lengths(dup_header, expected_lengths,
                          num_refs, "duplicate header") < 0)
        goto cleanup;

    if (sam_hdr_count_lines(dup_header, "SQ") != num_refs) {
        fail("Wrong number of SQ lines in duplicate header");
        goto cleanup;
    }

    if (check_ref_lengths(dup_header, expected_lengths,
                          num_refs, "parsed duplicate header") < 0)
        goto cleanup;

    if (sam_hdr_write(out, header) < 0) {
        fail("Failed to write SAM header");
        goto cleanup;
    }
    i = 0;
    while ((r = sam_read1(in, header, aln)) >= 0) {
        if (i >= num_align) {
            fail("Too many alignment records.\n");
            goto cleanup;
        }
        if (aln->core.tid != expected_tids[i]) {
            fail("Wrong tid for record %d : expected %d got %d\n",
                 i, expected_tids[i], aln->core.tid);
            goto cleanup;
        }
        if (aln->core.mtid != expected_mtid[i]) {
            fail("Wrong mate tid for record %d : expected %d got %d\n",
                 i, expected_mtid[i], aln->core.mtid);
            goto cleanup;
        }
        if (aln->core.pos != expected_positions[i]) {
            fail("Wrong position for record %d : "
                 "expected %"PRIhts_pos" got %"PRIhts_pos"\n",
                 i, expected_positions[i], aln->core.pos);
        }
        if (aln->core.mpos != expected_mpos[i]) {
            fail("Wrong mate position for record %d : "
                 "expected %"PRIhts_pos" got %"PRIhts_pos"\n",
                 i, expected_mpos[i], aln->core.mpos);
        }
        if (sam_write1(out, header, aln) < 0) {
            fail("Failed to write alignment record %d\n", i);
            goto cleanup;
        }
        i++;
    }
    if (r < -1) {
        fail("Error reading SAM alignment\n");
        goto cleanup;
    }
    if (i < num_align) {
        fail("Not enough alignment records\n");
        goto cleanup;
    }
    r = sam_close(in); in = NULL;
    if (r < 0) {
        fail("sam_close(in)");
        goto cleanup;
    }
    r = sam_close(out); out = NULL;
    if (r < 0) {
        fail("sam_close(out)");
        goto cleanup;
    }

    inf = fopen(outfname, "r");
    if (!inf) {
        fail("Opening \"%s\"", outfname);
        goto cleanup;
    }
    bytes = fread(buffer, 1, sizeof(buffer), inf);
    if (bytes != sizeof(sam_text) - 7
        || memcmp(buffer, sam_text + 6, bytes - 7) != 0) {
        fail("Output file does not match original version");
        fprintf(stderr,
                "---------- Expected:\n%.*s\n"
                "++++++++++ Got:\n%.*s\n"
                "====================\n",
                (int) sizeof(sam_text) - 7, sam_text + 6,
                (int) bytes, buffer);
        goto cleanup;
    }

 cleanup:
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_hdr_destroy(dup_header);
    if (in) sam_close(in);
    if (out) sam_close(out);
    if (inf) fclose(inf);
    unlink(outfname);
    return;
}

static void faidx1(const char *filename)
{
    int n, n_exp = 0, n_fq_exp = 0;
    char tmpfilename[FILENAME_MAX], line[500];
    FILE *fin, *fout;
    faidx_t *fai;

    fin = fopen(filename, "rb");
    if (fin == NULL) fail("can't open %s", filename);
    sprintf(tmpfilename, "%s.tmp", filename);
    fout = fopen(tmpfilename, "wb");
    if (fout == NULL) fail("can't create temporary %s", tmpfilename);
    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '>') n_exp++;
        if (line[0] == '+' && line[1] == '\n') n_fq_exp++;
        fputs(line, fout);
    }
    fclose(fin);
    fclose(fout);

    if (n_exp == 0 && n_fq_exp != 0) {
        // probably a fastq file
        n_exp = n_fq_exp;
    }

    if (fai_build(tmpfilename) < 0) fail("can't index %s", tmpfilename);
    fai = fai_load(tmpfilename);
    if (fai == NULL) { fail("can't load faidx file %s", tmpfilename); return; }

    n = faidx_fetch_nseq(fai);
    if (n != n_exp)
        fail("%s: faidx_fetch_nseq returned %d, expected %d", filename, n, n_exp);

    n = faidx_nseq(fai);
    if (n != n_exp)
        fail("%s: faidx_nseq returned %d, expected %d", filename, n, n_exp);

    fai_destroy(fai);
}

static void test_empty_sam_file(const char *filename)
{
    samFile *in = sam_open(filename, "r");
    if (in) {
        enum htsExactFormat format = hts_get_format(in)->format;
        bam1_t *aln = bam_init1();
        sam_hdr_t *header = sam_hdr_read(in);
        int ret = sam_read1(in, header, aln);

        if (format != empty_format)
            fail("detected %s as %d (expected empty_format)", filename, format);
        if (header)
            fail("sam_hdr_read() from %s should fail", filename);
        if (ret >= -1)
            fail("sam_read1() from %s returned %d but should fail", filename, ret);

        bam_destroy1(aln);
        sam_hdr_destroy(header);
        sam_close(in);
    }
    else fail("can't open %s to read as SAM", filename);
}

static void test_text_file(const char *filename, int nexp)
{
    htsFile *in = hts_open(filename, "r");
    if (in) {
        kstring_t str = KS_INITIALIZE;
        int ret, n = 0;
        while ((ret = hts_getline(in, '\n', &str)) >= 0) n++;
        if (ret != -1) fail("hts_getline got an error from %s", filename);
        if (n != nexp) fail("hts_getline read %d lines from %s (expected %d)", n, filename, nexp);

        hts_close(in);
        free(str.s);
    }
    else fail("can't open %s to read as text", filename);
}

static void check_enum1(void)
{
    // bgzf_compression() returns int, but enjoys this correspondence
    if (no_compression != 0) fail("no_compression is %d", no_compression);
    if (gzip != 1) fail("gzip is %d", gzip);
    if (bgzf != 2) fail("bgzf is %d", bgzf);
}

static void check_cigar_tab(void)
{
    int i, n_neg = 0;

    for (i = 0; i < 256; ++i)
        if (bam_cigar_table[i] < 0) n_neg++;

    if (n_neg + strlen(BAM_CIGAR_STR) != 256)
        fail("bam_cigar_table has %d unset entries", n_neg);

    for (i = 0; BAM_CIGAR_STR[i]; ++i)
        if (bam_cigar_table[(unsigned char) BAM_CIGAR_STR[i]] != i)
            fail("bam_cigar_table['%c'] is not %d", BAM_CIGAR_STR[i], i);
}

#define MAX_RECS 1000
#define SEQ_LEN 100
#define REC_LENGTH 150 // Undersized so some won't fit.

static int generator(const char *name)
{
    FILE *f = fopen(name, "w");
    char *ref = NULL;
    char qual[101];
    size_t i;
    uint32_t lfsr = 0xbadcafe;
    int res = -1;

    if (!f) {
        fail("Couldn't open \"%s\"", name);
        return -1;
    }

    ref = malloc(MAX_RECS + SEQ_LEN + 1);
    if (!ref) goto cleanup;
    for (i = 0; i < MAX_RECS + SEQ_LEN; i++) {
        // Linear-feedback shift register to make random reference
        lfsr ^= lfsr << 13;
        lfsr ^= lfsr >> 17;
        lfsr ^= lfsr << 5;
        ref[i] = "ACGT"[lfsr & 3];
    }
    ref[MAX_RECS + SEQ_LEN] = '\0';
    for (i = 0; i < SEQ_LEN; i++) {
        qual[i] = 'A' + (i & 0xf);
    }

    if (fputs("@HD\tVN:1.4\n", f) < 0) goto cleanup;
    if (fprintf(f, "@SQ\tSN:ref1\tLN:%u\n", MAX_RECS + SEQ_LEN) < 0)
        goto cleanup;
    for (i = 0; i < MAX_RECS; i++) {
        if (fprintf(f, "read%zu\t0\tref1\t%zu\t64\t100M\t*\t0\t0\t%.*s\t%.*s\n",
                    i + 1, i + 1, SEQ_LEN, ref + i, SEQ_LEN, qual) < 0)
            goto cleanup;
    }

    if (fclose(f) == 0)
        res = 0;
    f = NULL;

 cleanup:
    if (f) fclose(f);
    free(ref);
    return res;
}

static int read_data_block(const char *in_name, samFile *fp_in,
                           const char *out_name, samFile *fp_out,
                           sam_hdr_t *header, bam1_t *recs, size_t max_recs,
                           uint8_t *buffer, size_t bufsz, size_t *nrecs_out) {
    size_t buff_used = 0, nrecs;
    uint32_t new_m_data;
    int ret = -1, res = -1;

    for (nrecs = 0; nrecs < max_recs; nrecs++) {
        bam_set_mempolicy(&recs[nrecs],
                          BAM_USER_OWNS_STRUCT|BAM_USER_OWNS_DATA);

        recs[nrecs].data = &buffer[buff_used];
        recs[nrecs].m_data = bufsz - buff_used;

        res = sam_read1(fp_in, header, &recs[nrecs]);
        if (res < 0) break; // EOF or error

        if (fp_out) {
            if (sam_write1(fp_out, header, &recs[nrecs]) < 0) {
                nrecs++; // To return correct count
                fail("sam_write1() to \"%s\"", out_name);
                goto out;
            }
        }

        if ((bam_get_mempolicy(&recs[nrecs]) & BAM_USER_OWNS_DATA) == 0) {
            continue;  // Data not put in buffer
        }

        new_m_data = ((uint32_t) recs[nrecs].l_data + 7) & (~7U);
        if (new_m_data < recs[nrecs].m_data) recs[nrecs].m_data = new_m_data;

        buff_used += recs[nrecs].m_data;
    }
    if (res < -1) {
        fail("sam_read1() from \"%s\" failed", in_name);
    } else {
        ret = 0;
    }

 out:
    *nrecs_out = nrecs;
    return ret;
}

static void test_mempolicy(void)
{
    size_t bufsz = MAX_RECS * REC_LENGTH, nrecs = 0, i;
    bam1_t *recs = calloc(MAX_RECS, sizeof(bam1_t));
    uint8_t *buffer = malloc(bufsz);
    const char *fname = "test/sam_alignment.tmp.sam";
    const char *bam_name = "test/sam_alignment.tmp.bam";
    const char *cram_name = "test/sam_alignment.tmp.cram";
    const char *tag_text =
        "lengthy text ... lengthy text ... lengthy text ... lengthy text ... "
        "lengthy text ... lengthy text ... lengthy text ... lengthy text ... "
        "lengthy text ... lengthy text ... lengthy text ... lengthy text ... "
        "lengthy text ... lengthy text ... lengthy text ... lengthy text ... "
        "lengthy text ... lengthy text ... lengthy text ... lengthy text ... ";
    int res = 0;
    samFile *fp = NULL, *bam_fp = NULL, *cram_fp = NULL;
    htsFormat cram_fmt;
    sam_hdr_t *header = NULL;

    if (!recs || !buffer) {
        fail("Allocating buffer");
        goto cleanup;
    }

    memset(&cram_fmt, 0, sizeof(cram_fmt));

    // Make test file
    if (generator(fname) < 0)
        goto cleanup;

    // Open and read header
    fp = sam_open(fname, "r");
    if (!fp) {
        fail("sam_open(\"%s\")", fname);
        goto cleanup;
    }

    bam_fp = sam_open(bam_name, "wb");
    if (!fp) {
        fail("sam_open(\"%s\")", bam_name);
        goto cleanup;
    }

    header = sam_hdr_read(fp);
    if (!header) {
        fail("read header from \"%s\"", fname);
        goto cleanup;
    }

    if (sam_hdr_write(bam_fp, header) < 0) {
        fail("sam_hdr_write() to \"%s\"", bam_name);
        goto cleanup;
    }

    if (read_data_block(fname, fp, bam_name, bam_fp, header, recs,
                        MAX_RECS, buffer, bufsz, &nrecs) < 0)
        goto cleanup;

    res = sam_close(bam_fp);
    bam_fp = NULL;
    if (res < 0) {
        fail("sam_close(\"%s\")", bam_name);
        goto cleanup;
    }

    // Add a big tag to some records so they no longer fit in the allocated
    // buffer space.
    for (i = 0; i < MAX_RECS; i += 11) {
        if (bam_aux_update_str(&recs[i], "ZZ",
                               sizeof(tag_text) - 1, tag_text) < 0) {
            fail("bam_aux_update_str()");
            goto cleanup;
        }
    }

    // Delete all the records.  bam_destroy1() should free the data
    // for the ones that were expanded.
    for (i = 0; i < nrecs; i++) {
        bam_destroy1(&recs[i]);
    }

    res = sam_close(fp);
    fp = NULL;
    if (res < 0) {
        fail("sam_close(\"%s\")", fname);
        goto cleanup;
    }

    // Same test but reading BAM, writing CRAM
    nrecs = 0;
    sam_hdr_destroy(header);
    header = NULL;

    bam_fp = sam_open(bam_name, "r");
    if (!bam_fp) {
        fail("sam_open(\"%s\", \"r\")", bam_name);
        goto cleanup;
    }

    if (hts_parse_format(&cram_fmt, "cram,no_ref") < 0) {
        fail("hts_parse_format");
        goto cleanup;
    }
    cram_fp = hts_open_format(cram_name, "wc", &cram_fmt);
    if (!cram_fp) {
        fail("hts_open_format(\"%s\", \"wc\")", cram_name);
        goto cleanup;
    }

    header = sam_hdr_read(bam_fp);
    if (!header) {
        fail("read header from \"%s\"", bam_name);
        goto cleanup;
    }

    if (sam_hdr_write(cram_fp, header) < 0) {
        fail("sam_hdr_write() to \"%s\"", cram_name);
        goto cleanup;
    }

    if (read_data_block(bam_name, bam_fp, cram_name, cram_fp, header, recs,
                        MAX_RECS, buffer, bufsz, &nrecs) < 0)
        goto cleanup;

    res = sam_close(cram_fp);
    cram_fp = NULL;
    if (res < 0) {
        fail("sam_close(\"%s\")", cram_name);
        goto cleanup;
    }

    for (i = 0; i < MAX_RECS; i += 11) {
        if (bam_aux_update_str(&recs[i], "ZZ",
                               sizeof(tag_text) - 1, tag_text) < 0) {
            fail("bam_aux_update_str()");
            goto cleanup;
        }
    }

    for (i = 0; i < nrecs; i++) {
        bam_destroy1(&recs[i]);
    }

    // Now try reading the cram file
    nrecs = 0;
    sam_hdr_destroy(header);
    header = NULL;

    cram_fp = sam_open(cram_name, "r");
    if (!cram_fp) {
        fail("sam_open(\"%s\", \"r\")", cram_name);
        goto cleanup;
    }

    header = sam_hdr_read(cram_fp);
    if (!header) {
        fail("read header from \"%s\"", cram_name);
        goto cleanup;
    }

    if (read_data_block(cram_name, cram_fp, NULL, NULL, header, recs,
                        MAX_RECS, buffer, bufsz, &nrecs) < 0)
        goto cleanup;

    for (i = 0; i < MAX_RECS; i += 11) {
        if (bam_aux_update_str(&recs[i], "ZZ",
                               sizeof(tag_text) - 1, tag_text) < 0) {
            fail("bam_aux_update_str()");
            goto cleanup;
        }
    }

 cleanup:
    sam_hdr_destroy(header);
    if (fp) sam_close(fp);
    if (bam_fp) sam_close(bam_fp);
    if (cram_fp) sam_close(cram_fp);

    for (i = 0; i < nrecs; i++) {
        bam_destroy1(&recs[i]);
    }
    free(buffer);
    free(recs);
    if (cram_fmt.specific) {
        hts_opt_free(cram_fmt.specific);
    }
}

int main(int argc, char **argv)
{
    int i;

    status = EXIT_SUCCESS;

    aux_fields1();
    iterators1();
    samrecord_layout();
    use_header_api();
    test_header_pg_lines();
    test_header_updates();
    test_header_remove_lines();
    test_header_ref_altnames();
    test_empty_sam_file("test/emptyfile");
    test_text_file("test/emptyfile", 0);
    test_text_file("test/xx#pair.sam", 7);
    test_text_file("test/xx.fa", 7);
    test_text_file("test/fastqs.fq", 500);
    check_enum1();
    check_cigar_tab();
    check_big_ref(0);
    check_big_ref(1);
    test_mempolicy();
    set_qname();
    for (i = 1; i < argc; i++) faidx1(argv[i]);

    return status;
}
