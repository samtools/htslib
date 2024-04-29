/*  test/plugins-dlhts.c -- Test plugins with dynamically loaded libhts.

    Copyright (C) 2020 University of Glasgow.

    Author: John Marshall <John.W.Marshall@glasgow.ac.uk>

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

#if defined _WIN32 || defined __CYGWIN__ || defined __MSYS__
#define SKIP "running on Windows"
#elif !defined ENABLE_PLUGINS
#define SKIP "plugins being disabled"
#endif

#ifndef SKIP

#include <dlfcn.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>

#ifndef EPROTONOSUPPORT
#define EPROTONOSUPPORT ENOSYS
#endif

void *sym(void *htslib, const char *name)
{
    void *ptr = dlsym(htslib, name);
    if (ptr == NULL) {
        fprintf(stderr, "Can't find symbol \"%s\": %s\n", name, dlerror());
        exit(EXIT_FAILURE);
    }
    return ptr;
}

typedef void void_func(void);
void_func *func(void *htslib, const char *name) {
    void_func *fptr;
    *(void **) &fptr = sym(htslib, name);
    return fptr;
}

int errors = 0;
int verbose = 0;

struct hFILE;
typedef struct hFILE *hopen_func(const char *fname, const char *mode, ...);
typedef void hclose_abruptly_func(struct hFILE *fp);

hopen_func *hopen_p;
hclose_abruptly_func *hclose_abruptly_p;

void test_hopen(const char *fname, int expected)
{
    struct hFILE *fp = hopen_p(fname, "r");
    if (fp) {
        hclose_abruptly_p(fp);
        fprintf(stderr, "Opening \"%s\" actually succeeded\n", fname);
        errors++;
        return;
    }

    int supported = (errno != EPROTONOSUPPORT);
    if (supported != expected) {
        fprintf(stderr, "Opening \"%s\" failed badly: %s\n", fname, strerror(errno));
        errors++;
    }
    else if (verbose)
        printf("Opening \"%s\" produces %s\n", fname, strerror(errno));
}

void verbose_log(const char *message)
{
    fflush(stderr);
    if (verbose) puts(message);
    fflush(stdout);
}

int main(int argc, char **argv)
{
    int dlflags = RTLD_NOW;
    int skip = 0;
    int c;

    while ((c = getopt(argc, argv, "glv")) >= 0)
        switch (c) {
        case 'g': dlflags |= RTLD_GLOBAL; break;
        case 'l': dlflags |= RTLD_LOCAL;  break;
        case 'v': verbose++; break;
        }

    if (optind >= argc) {
        fprintf(stderr, "Usage: plugins-dlhts [-glv] LIBHTSFILE\n");
        return EXIT_FAILURE;
    }

    void *htslib = dlopen(argv[optind], dlflags);
    if (htslib == NULL) {
        fprintf(stderr, "Can't dlopen \"%s\": %s\n", argv[optind], dlerror());
        return EXIT_FAILURE;
    }

    if (verbose) {
        int *hts_verbosep = sym(htslib, "hts_verbose");
        *hts_verbosep += verbose;

        typedef const char *cstr_func(void);
        printf("Loaded HTSlib %s\n", ((cstr_func *) func(htslib, "hts_version"))());
    }

    hopen_p = (hopen_func *) func(htslib, "hopen");
    hclose_abruptly_p = (hclose_abruptly_func *) func(htslib, "hclose_abruptly");

    test_hopen("bad-scheme:unsupported", 0);

#ifdef __APPLE__
    /* Skip -l tests as we don't link plugins back to libhts on macOS, as this
       would conflict with a statically linked libhts.a on this platform. */
    skip = (dlflags & RTLD_LOCAL) != 0;
#endif

    if (! skip) {
#ifdef HAVE_LIBCURL
        test_hopen("https://localhost:99999/invalid_port", 1);
#endif
#ifdef ENABLE_GCS
        test_hopen("gs:invalid", 1);
#endif
#ifdef ENABLE_S3
        test_hopen("s3:invalid", 1);
#endif
    }
    else
        verbose_log("Skipping most tests");

    verbose_log("Calling hts_lib_shutdown()");
    (func(htslib, "hts_lib_shutdown"))();

    verbose_log("Calling dlclose(htslib)");
    if (dlclose(htslib) < 0) {
        fprintf(stderr, "Can't dlclose \"%s\": %s\n", argv[optind], dlerror());
        errors++;
    }

    verbose_log("Returning from main()");

    if (errors > 0) {
        printf("FAILED: %d errors\n", errors);
        return EXIT_FAILURE;
    }

    if (verbose) printf("All tests passed\n");
    return EXIT_SUCCESS;
}

#else

int main(void)
{
    printf("Tests skipped due to " SKIP "\n");
    return EXIT_SUCCESS;
}

#endif
