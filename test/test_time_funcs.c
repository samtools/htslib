/*  test_time_compat.c -- Test time functions

    Copyright (C) 2022 Genome Research Ltd.

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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

#include "../hts_time_funcs.h"

int test_normalised(time_t start, time_t end, time_t incr) {
    time_t i, j;
    struct tm *utc;

    for (i = start; i < end; i += incr) {
        utc = gmtime(&i);
        j = hts_time_gm(utc);
        if (i != j) {
            fprintf(stderr,
                    "hts_time_gm() failed, got %"PRId64" expected %"PRId64"\n",
                    (int64_t) j, (int64_t) i);
            return 1;
        }
    }
    return 0;
}

int test_specific(int year, int mon, int mday, int hour, int min, int sec,
                  time_t expected) {
    struct tm utc = { sec, min, hour, mday, mon - 1, year - 1900, 0, 0, 0 };
    time_t res = hts_time_gm(&utc);
    if (res != expected) {
        fprintf(stderr,
                "hts_time_gm() failed for %4d/%02d/%02d %02d:%02d:%02d :"
                " got %"PRId64" expected %"PRId64"\n",
                year, mon, mday, hour, min, sec,
                (int64_t) res, (int64_t) expected);
        return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    int res = 0;

    if (test_normalised(0, INT_MAX - 1000, 1000) != 0)
        return EXIT_FAILURE;
    if (sizeof(time_t) >= 8) {
        if (test_normalised(INT_MAX - 1000,
                            (time_t)((int64_t) INT_MAX * 2), 1000) != 0)
            return EXIT_FAILURE;
    }

    // 2022-06-14 12:32:10
    res |= test_specific(2022, 6, 14, 12, 32, 10, 1655209930);
    // 2022-06-14 12:32:10
    res |= test_specific(1993, 9, 10514, 12, 32, 10, 1655209930);
    // 2022-02-28 12:00:00
    res |= test_specific(2020, 2, 28, 12, 0, 0, 1582891200);
    // 2022-02-29 12:00:00
    res |= test_specific(2020, 2, 29, 12, 0, 0, 1582977600);
    // 2022-03-01 12:00:00
    res |= test_specific(2020, 2, 30, 12, 0, 0, 1583064000);
    // 2022-02-29 12:00:00
    res |= test_specific(2020, 3, 0, 12, 0, 0, 1582977600);
    // 2020-02-01 12:00:00
    res |= test_specific(2019, 14, 1, 12, 0, 0, 1580558400);
    // 2020-03-01 12:00:00
    res |= test_specific(2019, 15, 1, 12, 0, 0, 1583064000);
    // 2021-03-01 12:00:00
    res |= test_specific(2019, 27, 1, 12, 0, 0, 1614600000);
    // 2024-02-01 12:00:00
    res |= test_specific(2019, 62, 1, 12, 0, 0, 1706788800);
    // 2024-03-01 12:00:00
    res |= test_specific(2019, 63, 1, 12, 0, 0, 1709294400);
    // 2020-12-31 23:59:59
    res |= test_specific(2021, 0, 31, 23, 59, 59, 1609459199);
    // 2020-03-01 12:00:00
    res |= test_specific(2021, -9, 1, 12, 0, 0, 1583064000);
    // 2020-02-01 12:00:00
    res |= test_specific(2021, -10, 1, 12, 0, 0, 1580558400);
    // 2019-02-01 12:00:00
    res |= test_specific(2021, -22, 1, 12, 0, 0, 1549022400);
    // 1970-01-01 00:00:00
    res |= test_specific(1970, 1, 1, 0, 0, 0, 0);
    // 2038-01-19 03:14:07
    res |= test_specific(1970, 1, 1, 0, 0, INT_MAX, INT_MAX);
    // 2038-01-19 03:14:07
    res |= test_specific(2038, 1, 19, 3, 14, 7, INT_MAX);
    if (sizeof(time_t) < 8) {
        // 2038-01-19 03:14:08
        res |= test_specific(2038, 1, 19, 3, 14, 8, (time_t) -1);
    } else {
        // 2038-01-19 03:14:08
        res |= test_specific(2038, 1, 19, 3, 14, 8,
                             (time_t)((int64_t) INT_MAX + 1));
    }

    return res == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
