/* test_parse_int.c -- Test integer string conversion

   Copyright (C) 2019 Genome Research Ltd.

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
#include <getopt.h>
#include "textutils_internal.h"

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
    const char sentinal = '#';

    // Positive value (unsigned)
    for (i = 1; i < 64; i++) {
        num = (1ULL << i) - 1;
        for (offset = i < 5 ? -(1LL << (i - 1)) : -16; offset <= 30; offset++) {
            efail = (offset > 0);
            snprintf(buffer, sizeof(buffer), "%" PRIu64 "%c",
                     num + offset, sentinal);

            uval = hts_str2uint(buffer, &end, i, &failed);
            if (failed != efail || uval != (!efail ? num + offset : num)
                || *end != sentinal) {
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
                     num + offset, sentinal);

            val = hts_str2int(buffer, &end, i + 1, &failed);
            if (failed != efail || val != (!efail ? num + offset : num)
                || *end != sentinal) {
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
                     num + offset + 1, sentinal);

            val = hts_str2int(buffer, &end, i + 1, &failed);
            // Cast of val to unsigned in this comparison avoids undefined
            // behaviour when checking INT64_MIN.
            if (failed != efail
                || -((uint64_t) val) != (!efail ? num + offset + 1 : num + 1)
                || *end != sentinal) {
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
                 offset, sentinal);
        uval = hts_str2uint(buffer, &end, 64, &failed);
        if (failed != efail
            || uval != (efail ? UINT64_MAX : 18446744073709551000ULL + offset)
            || *end != sentinal) {
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
    return res ? EXIT_FAILURE : EXIT_SUCCESS;
}
