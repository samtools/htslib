/*  test_ksort.c -- ksort unit tests

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
#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_CLOCK_GETTIME_CPUTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "../htslib/hts_alloc.h"
#include "../htslib/ksort.h"
#include "../htslib/khash.h"

typedef enum {
    MERGESORT = 0,
    INTROSORT,
    COMBSORT,
    HEAPSORT,
    NUM_ALGOS
} algorithm;

typedef enum {
    UNSTABLE,
    STABLE
} stability;

typedef struct alg_info {
    const char *name;
    stability is_stable;
} alg_info;

const static alg_info sort_algorithms[NUM_ALGOS] = {
    { "mergesort", STABLE   },
    { "introsort", UNSTABLE },
    { "combsort",  UNSTABLE },
    { "heapsort",  UNSTABLE }
};

typedef struct item {
    size_t val;
    size_t idx;
} item;

static unsigned long long ncomp = 0;

#define lt_item(a, b) (++ncomp, (a).val < (b).val)

KSORT_INIT_STATIC(item, item, lt_item)

static long long get_time(void) {
#ifdef HAVE_CLOCK_GETTIME_CPUTIME
    struct timespec ts;
    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts) < 0) {
        perror("clock_gettime");
        return -1;
    }
    return ts.tv_sec * 1000000000LL + ts.tv_nsec;
#else
    struct timeval tv;
    if (gettimeofday(&tv, NULL) < 0) {
        perror("gettimeofday");
        return -1;
    }
    return tv.tv_sec * 1000000LL + tv.tv_usec;
#endif
}

static item *copy(size_t n, item *a) {
    item *dupe = hts_malloc_p(sizeof(*a), n);
    if (!dupe) {
        fprintf(stderr, "Failed to allocate %zu bytes\n", sizeof(*a) * n);
        return NULL;
    }
    memcpy(dupe, a, n * sizeof(*a));
    return dupe;
}

static int check_order_unstable(size_t n, item *a, const char *test,
                                const char *input_order) {
    size_t i;
    for (i = 1; i < n; i++) {
        if (a[i-1].val > a[i].val) {
            fprintf(stderr, "%s %zu %s: items %zu,%zu out of order\n",
                    test, n, input_order, i-1, i);
            return 1;
        }
    }
    return 0;
}

static int check_order_stable(size_t n, item *a, const char *test,
                              const char *input_order) {
    size_t i;
    for (i = 1; i < n; i++) {
        if ((a[i-1].val > a[i].val)
            || (a[i-1].val == a[i].val && a[i-1].idx >= a[i].idx)) {
            fprintf(stderr, "%s %zu %s: items %zu,%zu out of order\n",
                    test, n, input_order, i-1, i);
            return 1;
        }
    }
    return 0;
}

static void print_timing(algorithm algo, const char *input_order, size_t n,
                         unsigned long long ncomp, long long elapsed) {

#ifdef HAVE_CLOCK_GETTIME_CPUTIME
    long long sec = elapsed / 1000000000;
    long long nsec = elapsed % 1000000000;
    printf("%s %zu %s: %llu comparisons, %lld.%09lld processor seconds\n",
           sort_algorithms[algo].name, n, input_order, ncomp, sec, nsec);
#else
    long long sec = elapsed / 1000000;
    long long usec = elapsed % 1000000;
    printf("%s %zu %s: %llu comparisons, %lld.%06lld wall-time seconds\n",
           sort_algorithms[algo].name, n, input_order, ncomp, sec, usec);
#endif
}

static int test_sort(algorithm algo, size_t n, item *a,
                     const char *input_order, int benchmark) {
    item *c = copy(n, a);
    int retval;
    long long t0, t1;
    if (!c)
        return 1;

    ncomp = 0;
    t0 = get_time();
    switch (algo) {
      case MERGESORT: ks_mergesort(item, n, c, NULL); break;
      case INTROSORT: ks_introsort(item, n, c); break;
      case COMBSORT:  ks_combsort(item, n, c); break;
      case HEAPSORT:
          ks_heapmake(item, n, c);
          ks_heapsort(item, n, c);
          break;
      default:
          fprintf(stderr, "test_sort: Unknown algorithm number: %d\n",
                  (int) algo);
          return 1;
    }
    t1 = get_time();
    if (sort_algorithms[algo].is_stable == STABLE) {
        retval = check_order_stable(n, c, sort_algorithms[algo].name,
                                    input_order);
    } else {
        retval = check_order_unstable(n, c, sort_algorithms[algo].name,
                                      input_order);
    }
    free(c);
    if (benchmark)
        print_timing(algo, input_order, n, ncomp, t1 - t0);
    return retval;
}

int run_sort_tests(const size_t n, const unsigned int algo_mask,
                   const int benchmark) {
    item *a = hts_malloc_p(sizeof(*a), n);
    size_t i;
    int ret = 0, alg;

    // Increasing vals
    for (i = 0; i < n; i++) a[i].val = a[i].idx = i;
    for (alg = 0; alg < NUM_ALGOS; alg++)
        if ((algo_mask & (1 << alg)) != 0)
            ret |= test_sort(alg, n, a, "increasing", benchmark);

    // Decreasing vals
    for (i = 0; i < n; i++) a[i].val = n - i;
    for (alg = 0; alg < NUM_ALGOS; alg++)
        if ((algo_mask & (1 << alg)) != 0)
            ret |= test_sort(alg, n, a, "decreasing", benchmark);

    // Randomised small vals
    for (i = 0; i < n; i++)
        a[i].val = (__ac_Wang_hash((khint_t) (i & 0xffffffffU))) & 0xffff;
    for (alg = 0; alg < NUM_ALGOS; alg++)
        if ((algo_mask & (1 << alg)) != 0)
            ret |= test_sort(alg, n, a, "random", benchmark);

    free(a);

    return ret;
}

static void show_usage(FILE *out, char *prog) {
    fprintf(out, "Usage : %s [-b]\n", prog);
    fprintf(out, " Options:\n");
    fprintf(out, "  -a    Algorithm to test (may be used more than once)\n");
    fprintf(out, "  -b    Enable benchmarking\n");
    fprintf(out, "  -h    Show this help\n");
    fprintf(out, "  -n    Number of items to sort\n");
}

int main(int argc, char **argv) {
    int opt, ret = 0;
    size_t n = 0;

    int benchmark = 0;
    unsigned int algo_mask = ~0;

    while ((opt = getopt(argc, argv, "a:bhn:")) != -1) {
        switch (opt) {
        case 'a': {
            int alg;
            unsigned int found = 0;
            if (algo_mask & (1U << NUM_ALGOS))
                algo_mask = 0;
            for (alg = 0; alg < NUM_ALGOS; alg++) {
                if (strcmp(sort_algorithms[alg].name, optarg) == 0)
                    found = (1U << alg);
            }
            if (!found) {
                fprintf(stderr, "Unknown sort algorithm: %s\n", optarg);
                return EXIT_FAILURE;
            }
            algo_mask |= found;
            break;
        }
        case 'b':
            benchmark = 1;
            break;
        case 'h':
            show_usage(stdout, argv[0]);
            return EXIT_SUCCESS;
        case 'n':
            n = strtoul(optarg, NULL, 0);
            break;
        default:
            show_usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (n > 0)
        return run_sort_tests(n, algo_mask, benchmark);

    for (n = 0; n <= 16; n++) {
        ret |= run_sort_tests(n, algo_mask, benchmark);
    }

    // Easily divisible size
    ret |= run_sort_tests(131072, algo_mask, benchmark);
    // Not easily divisible size
    ret |= run_sort_tests(524287, algo_mask, benchmark);
    return ret;
}
