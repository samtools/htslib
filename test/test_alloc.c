/*  test_alloc.c -- hts_alloc unit tests

    Copyright (C) 2026 Genome Research Ltd.

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
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "../htslib/hts_alloc.h"

typedef struct {
    size_t      a;
    size_t      b;
    size_t      expected;
} test_case_t;

static int try_add_sat2(size_t a, size_t b, size_t expected, int verbose) {
    size_t r_ab = hts_add_sat2(a, b);
    size_t r_ba = hts_add_sat2(b, a);
    int fails = 0;
    if (r_ab != expected) {
        fails++;
        fprintf(stderr,
                "Failed: hts_add_sat2(%zu, %zu) got %zu expected %zu\n",
                a, b, r_ab, expected);
    } else if (verbose) {
        fprintf(stderr, "hts_add_sat2(%zu, %zu) = %zu\n", a, b, r_ab);
    }
    if (r_ba != expected) {
        fails++;
        fprintf(stderr,
                "Failed: hts_add_sat2(%zu, %zu) got %zu expected %zu\n",
                b, a, r_ba, expected);
    } else if (verbose) {
        fprintf(stderr, "hts_add_sat2(%zu, %zu) = %zu\n", b, a, r_ba);
    }
    return fails;
}

static int test_add_sat2(int verbose) {
    const test_case_t tests[] = {
        { 0, 0, 0 },
        { 0, 1, 1 },
        { 1, 1, 2 },
        { 0, SIZE_MAX - 1, SIZE_MAX - 1 },
        { 1, SIZE_MAX - 1, SIZE_MAX },
        { 2, SIZE_MAX - 1, SIZE_MAX },
        { 0, SIZE_MAX, SIZE_MAX },
        { 0, SIZE_MAX, SIZE_MAX },
        { SIZE_MAX, SIZE_MAX, SIZE_MAX },
    };
    size_t n;
    int fails = 0;

    for (n = 0; n < sizeof(tests)/sizeof(tests[0]); n++) {
        fails += try_add_sat2(tests[n].a, tests[n].b, tests[n].expected, verbose);
    }
    return fails;
}

static int try_prod_sat2(size_t a, size_t b, size_t expected, int verbose) {
    size_t r_ab = hts_prod_sat2(a, b);
    size_t r_ba = hts_prod_sat2(b, a);
    int fails = 0;
    if (r_ab != expected) {
        fails++;
        fprintf(stderr,
                "Failed: hts_prod_sat2(%zu, %zu) got %zu expected %zu\n",
                a, b, r_ab, expected);
    } else if (verbose) {
        fprintf(stderr, "hts_prod_sat2(%zu, %zu) = %zu\n", a, b, r_ab);
    }
    if (r_ba != expected) {
        fails++;
        fprintf(stderr,
                "Failed: hts_prod_sat2(%zu, %zu) got %zu expected %zu\n",
                b, a, r_ba, expected);
    } else if (verbose) {
        fprintf(stderr, "hts_prod_sat2(%zu, %zu) = %zu\n", b, a, r_ba);
    }
    return fails;
}

static int test_prod_sat2(int verbose) {
    const test_case_t tests[] = {
        { 0, 0, 0 },
        { 0, 1, 0 },
        { 1, 1, 1 },
        { 2, 3, 6 },
        { 0, SIZE_MAX, 0 },
        { 1, SIZE_MAX - 1, SIZE_MAX - 1 },
        { 2, SIZE_MAX - 1, SIZE_MAX },
        { 0, SIZE_MAX, 0 },
        { 1, SIZE_MAX, SIZE_MAX },
        { 2, SIZE_MAX, SIZE_MAX },
        { SIZE_MAX, SIZE_MAX, SIZE_MAX },
    };
    size_t n, size_bits = sizeof(size_t) * 8;
    int fails = 0;

    for (n = 0; n < sizeof(tests)/sizeof(tests[0]); n++) {
        fails += try_prod_sat2(tests[n].a, tests[n].b, tests[n].expected,
                               verbose);
    }

    for (n = 0; n < size_bits; n++) {
        size_t a = ((size_t) 1) << n;
        size_t b = ((size_t) 1) << (size_bits - 1 - n);
        // Should not overflow
        fails += try_prod_sat2(a, b, ((size_t) 1) << (size_bits - 1), verbose);
        if (n > 0) {
            // Should overflow
            fails += try_prod_sat2(a, b * 2, SIZE_MAX, verbose);
        }
    }

    return fails;
}

static void ptr2hex(char *out, size_t sz, void *ptr) {
    if (ptr) {
        snprintf(out, sz, "0x%p", ptr);
    } else {
        snprintf(out, sz, "NULL");
    }
}

static int alloc_result(const char *prefix, size_t a, size_t b,
                        const char *suffix, void *ptr, size_t should_work,
                        int verbose) {
    char pbuf[32];
    int failed = 0;
    if ((ptr != NULL) ^ should_work) {
        failed = 1;
        ptr2hex(pbuf, sizeof(pbuf), ptr);
        fprintf(stderr, "Failed: %s%zu, %zu%s got %s expected %s\n",
                prefix, a, b, suffix, pbuf, should_work ? "pointer" : "NULL");
    } else if (verbose) {
        ptr2hex(pbuf, sizeof(pbuf), ptr);
        fprintf(stderr, "%s%zu, %zu%s = %s\n", prefix, a, b, suffix, pbuf);
    }
    return failed;
}

static int test_malloc(unsigned long max_mb, int verbose) {
    const size_t size_limit = PTRDIFF_MAX < SIZE_MAX ? PTRDIFF_MAX : SIZE_MAX;
    const size_t max_bytes = hts_prod_sat2(1024*1024, max_mb);
    test_case_t tests[] = {
        { 1, 1, 1 },
        { max_bytes, 1, 1 },
        { 1, max_bytes, 1 },
        { max_bytes / 2, 2, 1 },
        { 2, max_bytes / 2, 1 },
        { size_limit, 1, 0 },
        { 1, size_limit, 0 },
        { 2, SIZE_MAX / 2 + 100,  0 },
        { SIZE_MAX / 2 + 100, 2, 0 },
    };
    size_t n;
    int fails = 0, i;
    for (n = 0; n < sizeof(tests)/sizeof(tests[0]); n++) {
        void *ptr = hts_malloc_p(tests[n].a, tests[n].b);
        fails += alloc_result("hts_malloc_p(", tests[n].a, tests[n].b, ")",
                              ptr, tests[n].expected, verbose);
        free(ptr);

        if (tests[n].b > 100) {
            char suffix[32];
            ptr = hts_malloc_ps(tests[n].a, tests[n].b - 100, 100);
            fails += alloc_result("hts_malloc_ps(",
                                  tests[n].a, tests[n].b - 100, ", 100)",
                                  ptr, tests[n].expected, verbose);
            free(ptr);

            ptr = hts_calloc_ps(tests[n].a, tests[n].b - 100, 100);
            fails += alloc_result("hts_calloc_ps(",
                                  tests[n].a, tests[n].b - 100, ", 100)",
                                  ptr, tests[n].expected, verbose);
            free(ptr);

            snprintf(suffix, sizeof(suffix), ", %zu)", tests[n].a * 100);
            ptr = hts_malloc_pse(tests[n].a, tests[n].b - 100, 0,
                                 tests[n].a * 100);
            fails += alloc_result("hts_malloc_pse(",
                                  tests[n].a, tests[n].b - 100, suffix,
                                  ptr, tests[n].expected, verbose);
            free(ptr);

            ptr = hts_calloc_pse(tests[n].a, tests[n].b - 100, 0,
                                 tests[n].a * 100);
            fails += alloc_result("hts_calloc_pse(",
                                  tests[n].a, tests[n].b - 100, suffix,
                                  ptr, tests[n].expected, verbose);
            free(ptr);
        }

        if (tests[n].expected) {
            ptr = hts_malloc_ps(tests[n].a, tests[n].b, size_limit);
            fails += alloc_result("hts_malloc_ps(",
                                  tests[n].a, tests[n].b, ", size_limit)",
                                  ptr, 0, verbose);
            free(ptr);

            ptr = hts_malloc_pse(tests[n].a, tests[n].b, 0, size_limit);
            fails += alloc_result("hts_malloc_pse(",
                                  tests[n].a, tests[n].b, ", 0, size_limit)",
                                  ptr, 0, verbose);
            free(ptr);
        }

        for (i = 0; i < 2; i++) {
            // i == 0 start with a NULL pointer
            // i == 1 start with an allocted pointer
            void *tmp = i ? malloc(1) : NULL;
            const char *prefix = i ? "hts_realloc_p(ptr, " : "hts_realloc_p(NULL, ";
            ptr = hts_realloc_p(tmp, tests[n].a, tests[n].b);
            fails += alloc_result(prefix, tests[n].a, tests[n].b, ")",
                                  ptr, tests[n].expected, verbose);
            free(ptr ? ptr : tmp);

            if (tests[n].b > 100) {
                char suffix[32];
                tmp = i ? malloc(1) : NULL;
                prefix = i ? "hts_realloc_ps(ptr, " : "hts_realloc_ps(NULL, ";
                ptr = hts_realloc_ps(tmp, tests[n].a, tests[n].b - 100, 100);
                fails += alloc_result(prefix,
                                      tests[n].a, tests[n].b - 100, ", 100)",
                                      ptr, tests[n].expected, verbose);
                free(ptr ? ptr : tmp);

                tmp = i ? malloc(1) : NULL;
                prefix = i ? "hts_realloc_pse(ptr, " : "hts_realloc_pse(NULL, ";
                snprintf(suffix, sizeof(suffix), ", %zu)", tests[n].a * 100);
                ptr = hts_realloc_pse(tmp, tests[n].a, tests[n].b - 100, 0,
                                      tests[n].a * 100);
                fails += alloc_result(prefix,
                                      tests[n].a, tests[n].b - 100, suffix,
                                      ptr, tests[n].expected, verbose);
                free(ptr ? ptr : tmp);
            }

            if (tests[n].expected) {
                tmp = i ? malloc(1) : NULL;
                prefix = i ? "hts_realloc_ps(ptr, " : "hts_realloc_ps(NULL, ";
                ptr = hts_realloc_ps(tmp, tests[n].a, tests[n].b, size_limit);
                fails += alloc_result("hts_realloc_ps(NULL, ",
                                      tests[n].a, tests[n].b, ", size_limit)",
                                      ptr, 0, verbose);
                free(ptr ? ptr : tmp);

                tmp = i ? malloc(1) : NULL;
                prefix = i ? "hts_realloc_pse(ptr, " : "hts_realloc_pse(NULL, ";
                ptr = hts_realloc_pse(tmp, tests[n].a, tests[n].b, 0, size_limit);
                fails += alloc_result("hts_realloc_pse(NULL, ",
                                      tests[n].a, tests[n].b, ", 0, size_limit)",
                                      ptr, 0, verbose);
                free(ptr ? ptr : tmp);
            }
        }
    }
    return fails;
}

void show_help(FILE *fp) {
    fprintf(fp, "Usage: test_alloc [-m <megabytes>] [-v]\n");
}

int main(int argc, char **argv) {
    int fails = 0;
    int verbose = 0;
    unsigned long max_mb = 10, max_max_mb = sizeof(size_t) > 4 ? 50000 : 2100;
    int opt;

    while ((opt = getopt(argc, argv, "m:v")) > 0) {
        switch (opt) {
        case 'm':
            max_mb = strtoul(optarg, NULL, 0);
            if (max_mb < 1 || max_mb > max_max_mb) {
                fprintf(stderr, "test_alloc: -m should be between 1 and %lu\n",
                        max_max_mb);
                show_help(stderr);
                return EXIT_FAILURE;
            }
            break;
        case 'v':
            verbose = 1;
            break;
        case 'h':
            show_help(stdout);
            return EXIT_SUCCESS;
        default:
            fprintf(stderr, "test_alloc: Unknown option '%c'\n", opt);
            show_help(stderr);
            return EXIT_FAILURE;
        }
    }

    fails += test_add_sat2(verbose);
    fails += test_prod_sat2(verbose);
    fails += test_malloc(max_mb, verbose);
    if (fails > 0) {
        fprintf(stderr, "test_alloc: %d tests failed\n", fails);
    }
    return fails == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
