/*  test_kstring.c -- kstring unit tests

    Copyright (C) 2018 Genome Research Ltd.

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
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>

#include <htslib/kstring.h>

static inline void clamp(int64_t *val, int64_t min, int64_t max) {
    if (*val < min) *val = min;
    if (*val > max) *val = max;
}

static int test_kputuw_from_to(kstring_t *str, unsigned int s, unsigned int e) {
    unsigned int i = s;

    for (;;) {
        str->l = 0;
        memset(str->s, 0xff, str->m);
        if (kputuw(i, str) < 0 || !str->s) {
            perror("kputuw");
            return -1;
        }
        if (str->l >= str->m || str->s[str->l] != '\0') {
            fprintf(stderr, "No NUL termination on string from kputuw\n");
            return -1;
        }
        if (i != strtoul(str->s, NULL, 10)) {
            fprintf(stderr,
                    "kputuw wrote the wrong value, expected %u, got %s\n",
                    i, str->s);
            return -1;
        }
        if (i >= e) break;
        i++;
    }
    return 0;
}

static int test_kputuw(int64_t start, int64_t end) {
    kstring_t str = { 0, 0, NULL };
    int64_t val;

    str.s = malloc(2);
    if (!str.s) {
        perror("malloc");
        return -1;
    }
    str.m = 2;

    for (val = 0; val < UINT_MAX; val = val == 0 ? 1 : val * 10) {
        unsigned int s = val == 0 ? 0 : val - 5;
        unsigned int e = val + 5;

        if (test_kputuw_from_to(&str, s, e) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    if (test_kputuw_from_to(&str, UINT_MAX - 5, UINT_MAX) < 0) {
        free(ks_release(&str));
        return -1;
    }

    str.m = 1; // Force a resize
    clamp(&start, 0, UINT_MAX);
    clamp(&end,   0, UINT_MAX);

    if (test_kputuw_from_to(&str, start, end) < 0) {
        free(ks_release(&str));
        return -1;
    }

    free(ks_release(&str));

    return 0;
}

static int test_kputw_from_to(kstring_t *str, int s, int e) {
    int i = s;

    for (;;) {
        str->l = 0;
        memset(str->s, 0xff, str->m);
        if (kputw(i, str) < 0 || !str->s) {
            perror("kputw");
            return -1;
        }
        if (str->l >= str->m || str->s[str->l] != '\0') {
            fprintf(stderr, "No NUL termination on string from kputw\n");
            return -1;
        }
        if (i != strtol(str->s, NULL, 10)) {
            fprintf(stderr,
                    "kputw wrote the wrong value, expected %u, got %s\n",
                    i, str->s);
            return -1;
        }
        if (i >= e) break;
        i++;
    }
    return 0;
}

static int test_kputw(int64_t start, int64_t end) {
    kstring_t str = { 0, 0, NULL };
    int64_t val;

    str.s = malloc(2);
    if (!str.s) {
        perror("malloc");
        return -1;
    }
    str.m = 2;

    for (val = 1; val < INT_MAX; val *= 10) {
        if (test_kputw_from_to(&str, val > 5 ? val - 5 : 0, val + 5) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    for (val = -1; val > INT_MIN; val *= 10) {
        if (test_kputw_from_to(&str, val - 5, val < -5 ? val + 5 : 0) < 0) {
            free(ks_release(&str));
            return -1;
        }
    }

    if (test_kputw_from_to(&str, INT_MAX - 5, INT_MAX) < 0) {
        free(ks_release(&str));
        return -1;
    }

    if (test_kputw_from_to(&str, INT_MIN, INT_MIN + 5) < 0) {
        free(ks_release(&str));
        return -1;
    }

    str.m = 1; // Force a resize
    clamp(&start, INT_MIN, INT_MAX);
    clamp(&end,   INT_MIN, INT_MAX);

    if (test_kputw_from_to(&str, start, end) < 0) {
        free(ks_release(&str));
        return -1;
    }

    free(ks_release(&str));

    return 0;
}

int main(int argc, char **argv) {
    int opt, res = EXIT_SUCCESS;
    int64_t start = 0;
    int64_t end   = 0;
    char *test    = NULL;

    while ((opt = getopt(argc, argv, "e:s:t:")) != -1) {
        switch (opt) {
        case 's':
            start = strtoll(optarg, NULL, 0);
            break;
        case 'e':
            end = strtoll(optarg, NULL, 0);
            break;
        case 't':
            test = optarg;
            break;
        default:
            fprintf(stderr, "Usage : %s [-s <num>] [-e <num>] [-t <test>]\n",
                    argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (!test || strcmp(test, "kputuw") == 0)
        if (test_kputuw(start, end) != 0) res = EXIT_FAILURE;

    if (!test || strcmp(test, "kputw") == 0)
        if (test_kputw(start, end) != 0) res = EXIT_FAILURE;

    return res;
}
