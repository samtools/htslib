/*  test/test_nibbles.c -- Test SIMD optimised function implementations.

    Copyright (C) 2024 Centre for Population Genomics.

    Author: John Marshall <jmarshall@hey.com>

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_CLOCK_GETTIME_CPUTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "../htslib/sam.h"
#include "../sam_internal.h"

long long gettime(void) {
#ifdef HAVE_CLOCK_GETTIME_CPUTIME
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec * 1000000000LL + ts.tv_nsec;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000000LL + tv.tv_usec;
#endif
}

char *fmttime(long long elapsed) {
    static char buf[64];

#ifdef HAVE_CLOCK_GETTIME_CPUTIME
    long long sec = elapsed / 1000000000;
    long long nsec = elapsed % 1000000000;
    sprintf(buf, "%lld.%09lld processor seconds", sec, nsec);
#else
    long long sec = elapsed / 1000000;
    long long usec = elapsed % 1000000;
    sprintf(buf, "%lld.%06lld wall-time seconds", sec, usec);
#endif

    return buf;
}

void nibble2base_single(uint8_t *nib, char *seq, int len) {
    int i;
    for (i = 0; i < len; i++)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

unsigned char nibble[5000];
char buf[10000];

int validate_nibble2base(void) {
    char defbuf[500];
    int i, start, len;
    unsigned long long total = 0, failed = 0;

    for (i = 0; i < sizeof nibble; i++)
        nibble[i] = i % 256;

    for (start = 0; start < 80; start++)
        for (len = 0; len < 400; len++) {
            memset(defbuf, '\0', sizeof defbuf);
            nibble2base_single(&nibble[start], defbuf, len);

            memset(buf, '\0', sizeof defbuf);
            nibble2base(&nibble[start], buf, len);

            total++;
            if (strcmp(defbuf, buf) != 0) {
                printf("%s expected\n%s FAIL\n\n", defbuf, buf);
                failed++;
            }
        }

    if (failed > 0) {
        fprintf(stderr, "Failures: %llu (out of %llu tests)\n", failed, total);
        return 1;
    }

    return 0;
}

int time_nibble2base(int length, unsigned long count) {
    unsigned long i, total = 0;

    for (i = 0; i < length; i++)
        nibble[i] = i % 256;

    printf("Timing %lu nibble2base iterations with read length %d...\n", count, length);
    long long start = gettime();

    for (i = 0; i < count; i++) {
        nibble2base(nibble, buf, length);
        total += buf[i % length];
    }

    long long stop = gettime();
    printf("%s (summing to %lu)\n", fmttime(stop - start), total);
    return 0;
}

int main(int argc, char **argv) {
    int readlen = 5000;
    unsigned long count = 1000000;
    int status = 0;
    int c;

    if (argc == 1)
        printf(
"Usage: test_nibbles [-c NUM] [-r NUM] [-n|-v]...\n"
"Options:\n"
"  -c NUM  Specify number of iterations [%lu]\n"
"  -n      Run nibble2base speed tests\n"
"  -r NUM  Specify read length [%d]\n"
"  -v      Run all validation tests\n"
"", count, readlen);

    while ((c = getopt(argc, argv, "c:nr:v")) >= 0)
        switch (c) {
        case 'c':
            count = strtoul(optarg, NULL, 0);
            break;

        case 'n':
            status += time_nibble2base(readlen, count);
            break;

        case 'r':
            readlen = atoi(optarg);
            break;

        case 'v':
            status += validate_nibble2base();
            break;
        }

    return status;
}
