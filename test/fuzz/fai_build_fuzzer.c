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

#include <dirent.h>
#include <errno.h>
#include <signal.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "../../htslib/faidx.h"

// Exercises the FASTA/FASTQ body parser (fai_build_core) and the .fai/.gzi
// index parser (fai_read), neither of which hts_open_fuzzer reaches.

#define TEMP_DIR "/tmp/"
#define TEMP_PREFIX "htsfai."

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

// Removes leftover temporary files from runs that were killed before they
// could clean up after themselves.
static void unlink_stale_temp_files(void) {
    DIR *dir = opendir(TEMP_DIR);
    struct dirent *entry;

    if (dir == NULL)
        return;

    while ((entry = readdir(dir)) != NULL) {
        const char *name = entry->d_name;
        char *end;
        long pid;
        char path[512];

        if (strncmp(name, TEMP_PREFIX, sizeof TEMP_PREFIX - 1) != 0)
            continue;
        pid = strtol(name + sizeof TEMP_PREFIX - 1, &end, 10);
        if (end == name + sizeof TEMP_PREFIX - 1 || pid <= 0
            || (*end != '\0' && *end != '.'))
            continue;
        if (kill((pid_t) pid, 0) == 0 || errno != ESRCH)
            continue;

        if (snprintf(path, sizeof path, TEMP_DIR "%s", name)
            < (int) sizeof path)
            unlink(path);
    }

    closedir(dir);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    char fa[64], idx[80];
    int fd;

    if (size == 0 || size > (1 << 20))
        return 0;

    // Clean up any files left behind by killed runs, then (re)create ours.
    unlink_stale_temp_files();

    snprintf(fa, sizeof fa, TEMP_DIR TEMP_PREFIX "%ld.XXXXXX",
             (long) getpid());
    fd = mkstemp(fa);
    if (fd < 0)
        return 0;
    if (write(fd, data, size) != (ssize_t) size) {
        close(fd);
        unlink(fa);
        return 0;
    }
    close(fd);

    // fai_load builds the index by parsing the FASTA/FASTQ body, then reads
    // it back; fetching a few sequences exercises the random-access seek
    // paths.
    faidx_t *fai = fai_load(fa);
    if (fai) {
        fetch_a_few(fai);
        fai_destroy(fai);
    }

    snprintf(idx, sizeof idx, "%s.fai", fa);
    unlink(idx);
    snprintf(idx, sizeof idx, "%s.gzi", fa);
    unlink(idx);
    unlink(fa);
    return 0;
}
