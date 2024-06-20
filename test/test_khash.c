/*  test_khash.c -- khash unit tests

    Copyright (C) 2024 Genome Research Ltd.
    Copyright (C) 2024 Centre for Population Genomics.

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
#ifdef HAVE_CLOCK_GETTIME_CPUTIME
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <htslib/khash.h>
#include <htslib/kroundup.h>

#define MAX_ENTRIES 99999999

KHASH_MAP_INIT_STR(str2int, int)

static void write_stats_str2int(khash_t(str2int) *h) {
    khint_t empty = 0, deleted = 0, hist_size = 0, *hist = NULL;

    if (kh_stats(str2int, h, &empty, &deleted, &hist_size, &hist) == 0) {
        khint_t i;
        printf("n_buckets = %u\n",
                kh_n_buckets(h));
        printf("empty     = %u\n", empty);
        printf("deleted   = %u\n", deleted);
        for (i = 0; i < hist_size; i++) {
            printf("dist[ %8u ] = %u\n", i, hist[i]);
        }
        free(hist);
    }
}

char * make_keys(size_t num, size_t kl) {
    size_t i;
    char *keys;

    if (num > MAX_ENTRIES) return NULL;
    keys = malloc(kl * num);
    if (!keys) {
        perror(NULL);
        return NULL;
    }
    for (i = 0; i < num; i++) {
        if (snprintf(keys + kl * i, kl, "test%zu", i) >= kl) {
            free(keys);
            return NULL;
        }
    }

    return keys;
}

static int add_str2int_entry(khash_t(str2int) *h, char *key, khint_t val) {
    int ret = 0;
    khint_t k = kh_put(str2int, h, key, &ret);

    if (ret != 1 && ret != 2) {
        fprintf(stderr, "Unexpected return from kh_put(%s) : %d\n", key, ret);
        return -1;
    }
    kh_val(h, k) = val;
    return 0;
}

static int check_str2int_entry(khash_t(str2int) *h, char *key, khint_t val,
                               uint8_t is_deleted) {
    khint_t k = kh_get(str2int, h, key);
    if (is_deleted) {
        if (k < kh_end(h)) {
            fprintf(stderr, "Found deleted entry %s in hash table\n", key);
            return -1;
        } else {
            return 0;
        }
    }

    if (k >= kh_end(h)) {
        fprintf(stderr, "Couldn't find %s in hash table\n", key);
        return -1;
    }
    if (strcmp(kh_key(h, k), key) != 0) {
        fprintf(stderr, "Wrong key in hash table, expected %s got %s\n",
                key, kh_key(h, k));
        return -1;
    }
    if (kh_val(h, k) != val) {
        fprintf(stderr, "Wrong value in hash table, expected %u got %u\n",
                val, kh_val(h, k));
        return -1;
    }
    return 0;
}

static int del_str2int_entry(khash_t(str2int) *h, char *key) {
    khint_t k = kh_get(str2int, h, key);
    if (k >= kh_end(h)) {
        fprintf(stderr, "Couldn't find %s to delete from hash table\n", key);
        return -1;
    }
    kh_del(str2int, h, k);
    return 0;
}

static int test_str2int(size_t max, size_t to_del, int show_stats) {
    const size_t kl = 16;
    size_t mask = max;
    char *keys = make_keys(max, kl);
    uint8_t *flags = NULL;
    khash_t(str2int) *h;
    khint_t i;
    uint32_t r = 0x533d;

    if (!keys) return -1;

    h = kh_init(str2int);
    if (!h) goto memfail;

    // Add some entries
    for (i = 0; i < max; i++) {
        if (add_str2int_entry(h, keys + i * kl, i) != 0)
            goto fail;
    }

    // Check they exist
    for (i = 0; i < max; i++) {
        if (check_str2int_entry(h, keys + i * kl, i, 0) != 0)
            goto fail;
    }

    if (show_stats) {
        printf("Initial fill:\n");
        write_stats_str2int(h);
    }

    // Delete a random selection
    flags = calloc(max, sizeof(*flags));
    if (!flags) {
        perror("");
        goto fail;
    }

    kroundup_size_t(mask);
    --mask;

    // Note that this method may become slow for a high %age removed
    // as it searches for the last available entries.  Despite this, it
    // seems to be acceptable for the number of entries allowed.
    for (i = 0; i < to_del; i++) {
        khint_t victim;
        // LFSR, see http://users.ece.cmu.edu/~koopman/lfsr/index.html
        do {
            r = (r >> 1) ^ ((r & 1) * 0x80000057U);
            victim = (r & mask) - 1;
        } while (victim >= max || flags[victim]);
        if (del_str2int_entry(h, keys + victim * kl) != 0)
            goto fail;
        flags[victim] = 1;
    }

    // Check correct entries are present
    for (i = 0; i < max; i++) {
        if (check_str2int_entry(h, keys + i * kl, i, flags[i]) != 0)
            goto fail;
    }

    if (show_stats) {
        printf("\nAfter deletion:\n");
        write_stats_str2int(h);
    }

    // Re-insert deleted entries
    for (i = 0; i < max; i++) {
        if (flags[i] && add_str2int_entry(h, keys + i * kl, i) != 0)
            goto fail;
    }

    // Ensure they're all back
    for (i = 0; i < max; i++) {
        if (check_str2int_entry(h, keys + i * kl, i, 0) != 0)
            goto fail;
    }

    if (show_stats) {
        printf("\nAfter re-insert:\n");
        write_stats_str2int(h);
    }

    kh_destroy(str2int, h);
    free(keys);
    free(flags);

    return 0;

 memfail:
    perror(NULL);
 fail:
    kh_destroy(str2int, h);
    free(keys);
    free(flags);
    return -1;
}

static size_t read_keys(const char *keys_file, char **keys_out,
                        char ***key_locations_out) {
    FILE *in = fopen(keys_file, "r");
    char *keys = NULL, *key, *end;
    size_t keys_size = 1000000;
    size_t keys_used = 0;
    size_t avail, got, nkeys = 0;
    char **key_locations = NULL;
    struct stat fileinfo = { 0 };

    if (!in)
        return 0;

    // Slurp entire file
    if (fstat(fileno(in), &fileinfo) < 0) {
        if (fileinfo.st_size > keys_size)
            keys_size = (size_t) fileinfo.st_size;
    }

    keys = malloc(keys_size + 1);
    if (!keys)
        goto fail;

    do {
        avail = keys_size - keys_used;
        if (avail == 0) {
            size_t new_size = keys_size + 1000000;
            char *new_keys = realloc(keys, new_size + 1);
            if (!new_keys)
                goto fail;
            keys = new_keys;
            keys_size = new_size;
            avail = keys_size - keys_used;
        }
        got = fread(keys + keys_used, 1, avail, in);
        keys_used += got;
    } while (got == avail);
    keys[keys_used] = '\0';

    if (ferror(in))
        goto fail;
    if (fclose(in) < 0)
        goto fail;
    in = NULL;

    // Split by line
    end = keys + keys_used;
    for (key = keys; key != NULL; key = memchr(key, '\n', end - key)) {
        while (*key == '\n') key++;
        if (key < end) nkeys++;
    }

    key_locations = malloc(nkeys * sizeof(*key_locations));
    if (!key_locations)
        goto fail;

    nkeys = 0;
    for (key = keys; key != NULL; key = memchr(key, '\n', end - key)) {
        while (*key == '\n') *key++ = '\0';
        if (key < end) {
            key_locations[nkeys++] = key;
        }
    }
    *keys_out = keys;
    *key_locations_out = key_locations;
    return nkeys;

 fail:
    if (in)
        fclose(in);
    free(keys);
    *keys_out = NULL;
    *key_locations_out = NULL;
    return 0;
}

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

static char * fmt_time(long long elapsed) {
    static char buf[64];
#ifdef HAVE_CLOCK_GETTIME_CPUTIME
    long long sec = elapsed / 1000000000;
    long long nsec = elapsed % 1000000000;
    snprintf(buf, sizeof(buf), "%lld.%09lld processor seconds", sec, nsec);
#else
    long long sec = elapsed / 1000000;
    long long usec = elapsed % 1000000;
    snprintf(buf, sizeof(buf), "%lld.%06lld wall-time seconds", sec, usec);
#endif
    return buf;
}

static int benchmark(const char *keys_file) {
    const size_t kl = 16;
    size_t max = 50000000;
    size_t i;
    char *keys = NULL;
    char **key_locations = NULL;
    khash_t(str2int) *h;
    long long start, end;

    if (keys_file) {
        max = read_keys(keys_file, &keys, &key_locations);
    } else {
        keys = make_keys(max, kl);
    }

    if (!keys) return -1;

    h = kh_init(str2int);
    if (!h) goto fail;

    if ((start = get_time()) < 0)
        goto fail;

    if (keys_file) {
        for (i = 0; i < max; i++) {
            int ret;
            khint_t k = kh_put(str2int, h, key_locations[i], &ret);
            if (ret < 0) {
                fprintf(stderr, "Unexpected return from kh_put(%s) : %d\n",
                        key_locations[i], ret);
                goto fail;
            }
            kh_val(h, k) = i;
        }
    } else {
        for (i = 0; i < max; i++) {
            int ret;
            khint_t k = kh_put(str2int, h, keys + i * kl, &ret);
            if (ret <= 0) {
                fprintf(stderr, "Unexpected return from kh_put(%s) : %d\n",
                        keys + i * kl, ret);
                goto fail;
            }
            kh_val(h, k) = i;
        }
    }

    if ((end = get_time()) < 0)
        goto fail;

    printf("Insert %zu %s\n", max, fmt_time(end - start));

    if ((start = get_time()) < 0)
        goto fail;

    if (keys_file) {
        for (i = 0; i < max; i++) {
            khint_t k = kh_get(str2int, h, key_locations[i]);
            if (k >= kh_end(h)) {
                fprintf(stderr, "Couldn't find %s in hash table\n",
                        key_locations[i]);
                goto fail;
            }
        }
    } else {
        for (i = 0; i < max; i++) {
            khint_t k = kh_get(str2int, h, keys + i * kl);
            if (k >= kh_end(h)) {
                fprintf(stderr, "Couldn't find %s in hash table\n",
                        keys + i * kl);
                goto fail;
            }
        }
    }

    if ((end = get_time()) < 0)
        goto fail;

    printf("Lookup %zu %s\n", max, fmt_time(end - start));

    write_stats_str2int(h);

    kh_destroy(str2int, h);
    free(keys);
    free(key_locations);

    return 0;
 fail:
    kh_destroy(str2int, h);
    free(keys);
    return -1;
}

static void show_usage(FILE *out, char *prog) {
    fprintf(out, "Usage : %s [-t <test>] [-i <file>]\n", prog);
    fprintf(out, " Options:\n");
    fprintf(out, "  -t <TEST>   Test to run (str2int, benchmark)\n");
    fprintf(out, "  -i <FILE>   Optional input file for benchmark\n");
    fprintf(out, "  -n <INT>    Number of items to add\n");
    fprintf(out, "  -f <FRAC>   Fraction to delete and re-insert\n");
    fprintf(out, "  -d          Dump hash table stats\n");
    fprintf(out, "  -h          Show this help\n");
}

int main(int argc, char **argv) {
    int opt, res = EXIT_SUCCESS;
    char *test = NULL;
    char *input_file = NULL;
    size_t max = 1000;
    double del_frac = 0.25;
    int show_stats = 0;

    while ((opt = getopt(argc, argv, "df:hi:n:t:")) != -1) {
        switch (opt) {
        case 'd':
            show_stats = 1;
            break;
        case 'f':
            del_frac = strtod(optarg, NULL);
            if (del_frac < 0 || del_frac > 1.0) {
                fprintf(stderr, "Error: -d must be between 0.0 and 1.0\n");
                return EXIT_FAILURE;
            }
            break;
        case 'h':
            show_usage(stdout, argv[0]);
            return EXIT_SUCCESS;
        case 'i':
            input_file = optarg;
            break;
        case 'n':
            max = strtoul(optarg, NULL, 0);
            if (max == 0 || max > 99999999) {
                fprintf(stderr, "Error: -n must be between 1 and %u\n",
                        MAX_ENTRIES);
                return EXIT_FAILURE;
            }
            break;
        case 't':
            test = optarg;
            break;
        default:
            show_usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (!test || strcmp(test, "str2int") == 0) {
        if (test_str2int(max, (size_t) (max * del_frac), show_stats) != 0)
            res = EXIT_FAILURE;
    }

    if (test && strcmp(test, "benchmark") == 0) {
        if (benchmark(input_file) != 0)
            res = EXIT_FAILURE;
    }

    return res;
}
