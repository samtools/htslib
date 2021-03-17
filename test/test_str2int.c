/* test/test_str2int.c -- Test integer string conversion (and safe printing)

   Copyright (C) 2019-2020 Genome Research Ltd.

   Author: Rob Davies <rmd@sanger.ac.uk>

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
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../textutils_internal.h"

// Test hts_str2int() and hts_str2uint() on various values around the
// maximum (or minimum for negative numbers) allowed for the given
// number of bits.  Ensures that the failed flag is set when the output
// isn't going to fit, that the correct value is returned and that
// 'end' points to the character following the number.
static int check_str2int(int verbose) {
    char buffer[64], *end;
    int64_t val;
    uint64_t num, uval;
    int failed = 0, efail, i, offset;
    const char sentinel = '#';

    // Positive value (unsigned)
    for (i = 1; i < 64; i++) {
        num = (1ULL << i) - 1;
        for (offset = i < 5 ? -(1LL << (i - 1)) : -16; offset <= 30; offset++) {
            efail = (offset > 0);
            snprintf(buffer, sizeof(buffer), "%" PRIu64 "%c",
                     num + offset, sentinel);

            uval = hts_str2uint(buffer, &end, i, &failed);
            if (failed != efail || uval != (!efail ? num + offset : num)
                || *end != sentinel) {
                fprintf(stderr, "hts_str2uint failed: %d bit "
                        "%s %"PRIu64" '%c' %d (%d)\n",
                        i, buffer, uval, *end, failed, efail);
                return -1;
            } else if (verbose) {
                fprintf(stderr, "hts_str2uint OK: %d bit "
                        "%s %"PRIu64" '%c' %d (%d)\n",
                        i, buffer, uval, *end, failed, efail);
            }
            failed = 0;
        }

        // Positive value (signed)
        for (offset = i < 5 ? -(1LL << (i - 1)) : -16; offset <= 30; offset++) {
            efail = (offset > 0);
            snprintf(buffer, sizeof(buffer), "%" PRIu64 "%c",
                     num + offset, sentinel);

            val = hts_str2int(buffer, &end, i + 1, &failed);
            if (failed != efail || val != (!efail ? num + offset : num)
                || *end != sentinel) {
                fprintf(stderr,
                        "hts_str2int  failed: %d bit "
                        "%s %"PRId64" '%c' %d (%d)\n",
                        i + 1, buffer, val, *end, failed, efail);
                return -1;
            } else if (verbose) {
                fprintf(stderr, "hts_str2int  OK: %d bit "
                        "%s %"PRId64" '%c' %d (%d)\n",
                        i + 1, buffer, val, *end, failed, efail);
            }
            failed = 0;
        }

        // Negative value (signed)
        for (offset = i < 5 ? -(1LL << (i - 1)) : -16; offset <= 30; offset++) {
            efail = (offset > 0);
            snprintf(buffer, sizeof(buffer), "-%" PRIu64 "%c",
                     num + offset + 1, sentinel);

            val = hts_str2int(buffer, &end, i + 1, &failed);
            // Cast of val to unsigned in this comparison avoids undefined
            // behaviour when checking INT64_MIN.
            if (failed != efail
                || -((uint64_t) val) != (!efail ? num + offset + 1 : num + 1)
                || *end != sentinel) {
                fprintf(stderr,
                        "hts_str2int  failed: %d bit "
                        "%s %"PRId64" '%c' %d (%d)\n",
                        i + 1, buffer, val, *end, failed, efail);
                return -1;
            } else if (verbose) {
                fprintf(stderr, "hts_str2int  OK: %d bit "
                        "%s %"PRId64" '%c' %d (%d)\n",
                        i + 1, buffer, val, *end, failed, efail);
            }
            failed = 0;
        }
    }

    // Special case for UINT64_MAX
    for (offset = 0; offset <= 999; offset++) {
        efail = offset > 615;
        snprintf(buffer, sizeof(buffer), "18446744073709551%03d%c",
                 offset, sentinel);
        uval = hts_str2uint(buffer, &end, 64, &failed);
        if (failed != efail
            || uval != (efail ? UINT64_MAX : 18446744073709551000ULL + offset)
            || *end != sentinel) {
            fprintf(stderr, "hts_str2uint failed: 64 bit %s "
                    "%"PRIu64" '%c' %d (%d)\n",
                    buffer, uval, *end, failed, efail);
            return -1;
        } else if (verbose) {
            fprintf(stderr, "hts_str2uint OK: 64 bit "
                    "%s %"PRIu64" '%c' %d (%d)\n",
                    buffer, uval, *end, failed, efail);
        }
    }
    return 0;
}

static int
check_strprint2(int verbose, const char *str, size_t len, size_t destlen,
                char quote, const char *expect) {
    char buf[100];
    hts_strprint(buf, destlen, quote, str, len);
    if (strcmp(buf, expect) != 0) {
        fprintf(stderr, "hts_strprint failed: length %zu: got \"%.*s\", "
                "expected \"%s\"\n", destlen, (int) destlen, buf, expect);
        return -1;
    }
    else if (verbose) {
        fprintf(stderr, "hts_strprint OK: length %zu: got \"%s\"\n",
                destlen, expect);
    }
    return 0;
}

static int
check_strprint1(int v, const char *str, size_t destlen, const char *expect) {
    return check_strprint2(v, str, SIZE_MAX, destlen, '\0', expect);
}

static int
check_strprintq(int v, const char *str, size_t destlen, char quote,
                const char *expect) {
    return check_strprint2(v, str, SIZE_MAX, destlen, quote, expect);
}

static int check_strprint(int v) {
    int res = 0;

    res |= check_strprint1(v, "chr10", 9, "chr10");
    res |= check_strprint1(v, "chr10", 6, "chr10");
    res |= check_strprint1(v, "chr10", 5, "c...");
    res |= check_strprint1(v, "chr10", 4, "...");
    res |= check_strprint1(v, "tab\twxyz",10, "tab\\twxyz");
    res |= check_strprint1(v, "tab\twxyz", 9, "tab\\t...");
    res |= check_strprint1(v, "tab\twxyz", 8, "tab\\...");
    res |= check_strprint1(v, "tab\twxyz", 7, "tab...");
    res |= check_strprint1(v, "tab\twxyz", 6, "ta...");
    res |= check_strprint1(v, "\xab", 5, "\\xAB");
    res |= check_strprint1(v, "\xab", 4, "...");
    res |= check_strprint1(v, "hello\xff", 40, "hello\\xFF");
    res |= check_strprint1(v, "hello\xff", 10, "hello\\xFF");
    res |= check_strprint1(v, "hello\xff", 9, "hello...");
    res |= check_strprint1(v, "hello\t", 40, "hello\\t");
    res |= check_strprint1(v, "hello\t", 8, "hello\\t");
    res |= check_strprint1(v, "hello\t", 7, "hel...");
    res |= check_strprint1(v, "\t", 40, "\\t");
    res |= check_strprint1(v, "", 40, "");

    res |= check_strprintq(v, "chr10", 9, '\'', "'chr10'");
    res |= check_strprintq(v, "chr10", 8, '\'', "'chr10'");
    res |= check_strprintq(v, "chr10", 7, '\'', "'c'...");
    res |= check_strprintq(v, "chr10", 6, '\'', "''...");
    res |= check_strprintq(v, "quo'wxyz",12, '\'', "'quo\\'wxyz'");
    res |= check_strprintq(v, "quo'wxyz",11, '\'', "'quo\\''...");
    res |= check_strprintq(v, "quo'wxyz",10, '\'', "'quo\\'...");

    res |= check_strprint2(v, "foo\0bar", SIZE_MAX, 10, '\0', "foo");
    res |= check_strprint2(v, "foo\0bar", 7,10, '\0', "foo\\0bar");
    res |= check_strprint2(v, "foo\0bar", 7, 9, '\0', "foo\\0bar");
    res |= check_strprint2(v, "foo\0bar", 7, 8, '\0', "foo\\...");

    return res;
}

int main(int argc, char **argv) {
    int verbose = 0, opt, res;

    while ((opt = getopt(argc, argv, "v")) != -1) {
        switch (opt) {
        case 'v':
            verbose = 1;
            break;
        default:
            fprintf(stderr, "Usage: %s [-v]\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

    res = check_str2int(verbose);
    res |= check_strprint(verbose);
    return res ? EXIT_FAILURE : EXIT_SUCCESS;
}
