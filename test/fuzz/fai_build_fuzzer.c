/*  test/fuzz/fai_build_fuzzer.c -- Fuzz driver for the FASTA/FASTQ index parser.

    Copyright (C) 2026 Ada Logics Ltd.

    Author: Arthur Chan <arthur.chan@adalogics.com>

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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../../htslib/faidx.h"

// Exercises the FASTA/FASTQ body parser (fai_build_core) and the .fai/.gzi
// index parser (fai_read), neither of which hts_open_fuzzer reaches.

#define TEMP_DIR "/tmp/"
#define TEMP_PREFIX "htsfai."

// Whether paths of the form "/proc/self/fd/<fd>" can be used to refer to
// files that have been unlinked but are still open (a Linux feature).
// Determined once, on the first call.
static int use_proc_fd = -1;

// Fetches (at most) the first 64 sequences, to exercise the random-access
// seek paths.
static void fetch_a_few(faidx_t *fai) {
    int n = faidx_nseq(fai);
    if (n > 64)
        n = 64;
    for (int i = 0; i < n; i++) {
        const char *name = faidx_iseq(fai, i);
        if (!name)
            continue;
        int len = 0;
        char *seq = fai_fetch(fai, name, &len);
        free(seq);
    }
}

// Creates a temporary file under TEMP_DIR (the only location the OSS-Fuzz
// sandbox allows us to write to) and returns an open descriptor for it.
//
// Where possible the file is unlinked immediately, keeping the descriptor
// open so it can still be referred to via its /proc/self/fd/<fd> name.
// The kernel reclaims unlinked files when the process exits, even if it is
// killed, so no clean-up code is needed and no stray files are ever left
// behind - unlike named temporary files, which survive if the fuzzer is
// killed before it can tidy up (e.g. on AFL++ timeouts).  As the files are
// anonymous, they cannot be confused with, or clobber, files made by other
// processes when fuzzers are run in parallel.
//
// On systems without /proc, the path is kept and removed on a best-effort
// basis by remove_temp_file().
static int create_temp_file(char *path, size_t path_size, const char *suffix) {
    char tmpl[64];
    int fd;

    snprintf(tmpl, sizeof tmpl, TEMP_DIR TEMP_PREFIX "%ld.%s.XXXXXX",
             (long) getpid(), suffix);
    fd = mkstemp(tmpl);
    if (fd < 0)
        return -1;

    if (use_proc_fd) {
        unlink(tmpl);
        snprintf(path, path_size, "/proc/self/fd/%d", fd);
    } else {
        snprintf(path, path_size, "%s", tmpl);
    }

    return fd;
}

static void remove_temp_file(const char *path) {
    if (!use_proc_fd)
        unlink(path);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    char fa_path[96], fai_path[96], gzi_path[96];
    int fd_fa, fd_fai, fd_gzi;
    faidx_t *fai;

    if (size == 0 || size > (1 << 20))
        return 0;

    if (use_proc_fd < 0) {
        struct stat st;
        use_proc_fd = (stat("/proc/self/fd", &st) == 0);
    }

    fd_fa = create_temp_file(fa_path, sizeof fa_path, "fa");
    if (fd_fa < 0)
        return 0;
    if (write(fd_fa, data, size) != (ssize_t) size) {
        close(fd_fa);
        remove_temp_file(fa_path);
        return 0;
    }

    // The index files are created empty, so build the index explicitly with
    // fai_build3() (fai_load3() with FAI_CREATE would find the empty files
    // and skip the build), then load it back; fetching a few sequences
    // exercises the random-access seek paths.
    fd_fai = create_temp_file(fai_path, sizeof fai_path, "fai");
    fd_gzi = create_temp_file(gzi_path, sizeof gzi_path, "gzi");
    if (fd_fai >= 0 && fd_gzi >= 0
        && fai_build3(fa_path, fai_path, gzi_path) == 0) {
        // Load as both FASTA and FASTQ, to exercise both index parsers
        // (fai_read).  Failed loads simply return NULL.
        fai = fai_load3(fa_path, fai_path, gzi_path, 0);
        if (fai) {
            fetch_a_few(fai);
            fai_destroy(fai);
        }

        fai = fai_load3_format(fa_path, fai_path, gzi_path, 0, FAI_FASTQ);
        if (fai) {
            fetch_a_few(fai);
            fai_destroy(fai);
        }
    }

    close(fd_fa);
    if (fd_fai >= 0)
        close(fd_fai);
    if (fd_gzi >= 0)
        close(fd_gzi);
    remove_temp_file(fa_path);
    remove_temp_file(fai_path);
    remove_temp_file(gzi_path);
    return 0;
}

