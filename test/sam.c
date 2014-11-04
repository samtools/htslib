/*  test/sam.c -- SAM/BAM/CRAM API test cases.

    Copyright (C) 2014 Genome Research Ltd.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "htslib/sam.h"
#include "htslib/kstring.h"

int status;

static void fail(const char *fmt, ...)
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
        else fail("%s field of type '%c', expected '%c'\n", tag, *p, type);
    }
    else fail("can't find %s field\n", tag);

    return NULL;
}

#define PI 3.141592653589793
#define E  2.718281828459045
#define HELLO "Hello, world!"
#define BEEF "DEADBEEF"

#define str(x) #x
#define xstr(x) str(x)

static int aux_fields1(void)
{
    static const char sam[] = "data:"
"@SQ\tSN:one\tLN:1000\n"
"@SQ\tSN:two\tLN:500\n"
"r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:" xstr(PI) "\tXd:d:" xstr(E) "\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,+2\tZZ:i:1000000\n";

    // Canonical form of the alignment record above, as output by sam_format1()
    static const char r1[] = "r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:3.14159\tXd:d:2.71828\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,2\tZZ:i:1000000";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    uint8_t *p;
    uint32_t n;
    kstring_t ks = { 0, 0, NULL };

    if (sam_read1(in, header, aln) >= 0) {
        if ((p = check_bam_aux_get(aln, "XA", 'A')) && bam_aux2A(p) != 'k')
            fail("XA field is '%c', expected 'k'", bam_aux2A(p));

        if ((p = check_bam_aux_get(aln, "Xi", 'C')) && bam_aux2i(p) != 37)
            fail("Xi field is %d, expected 37", bam_aux2i(p));

        if ((p = check_bam_aux_get(aln, "Xf", 'f')) && fabs(bam_aux2f(p) - PI) > 1E-6)
            fail("Xf field is %.12f, expected pi", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "Xd", 'd')) && fabs(bam_aux2f(p) - E) > 1E-6)
            fail("Xf field is %.12f, expected e", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "XZ", 'Z')) && strcmp(bam_aux2Z(p), HELLO) != 0)
            fail("XZ field is \"%s\", expected \"%s\"", bam_aux2Z(p), HELLO);

        if ((p = check_bam_aux_get(aln, "XH", 'H')) && strcmp(bam_aux2Z(p), BEEF) != 0)
            fail("XH field is \"%s\", expected \"%s\"", bam_aux2Z(p), BEEF);

        // TODO Invent and use bam_aux2B()
        if ((p = check_bam_aux_get(aln, "XB", 'B')) && ! (memcmp(p, "Bc", 2) == 0 && (memcpy(&n, p+2, 4), n) == 3 && memcmp(p+6, "\xfe\x00\x02", 3) == 0))
            fail("XB field is %c,..., expected c,-2,0,+2", p[1]);

        if ((p = check_bam_aux_get(aln, "ZZ", 'I')) && bam_aux2i(p) != 1000000)
            fail("ZZ field is %d, expected 1000000", bam_aux2i(p));

        if (sam_format1(header, aln, &ks) < 0)
            fail("can't format record");

        if (strcmp(ks.s, r1) != 0)
            fail("record formatted incorrectly: \"%s\"", ks.s);

        free(ks.s);
    }
    else fail("can't read record");

    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);

    return 1;
}

static void iterators1(void)
{
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0));
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_NONE, 0, 0));
}

int main(void)
{
    status = EXIT_SUCCESS;

    aux_fields1();
    iterators1();

    return status;
}
