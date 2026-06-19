/*  test/fuzz/fai_build_fuzzer.c -- Fuzz driver for the FASTA/FASTQ index parser.

    Copyright (C) 2026 Genome Research Ltd.

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
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "../../htslib/faidx.h"

// Exercises the FASTA/FASTQ body parser (fai_build_core) and the .fai/.gzi
// index parser (fai_read), neither of which hts_open_fuzzer reaches.
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size == 0 || size > (1 << 20))
        return 0;

    char fa[] = "/tmp/htsfai.XXXXXX";
    int fd = mkstemp(fa);
    if (fd < 0)
        return 0;
    if (write(fd, data, size) != (ssize_t) size) {
        close(fd);
        unlink(fa);
        return 0;
    }
    close(fd);

    // fai_load builds the index by parsing the FASTA/FASTQ body, then reads it
    // back; fetching a few sequences exercises the random-access seek paths.
    faidx_t *fai = fai_load(fa);
    if (fai) {
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
        fai_destroy(fai);
    }

    char idx[64];
    snprintf(idx, sizeof idx, "%s.fai", fa);
    unlink(idx);
    snprintf(idx, sizeof idx, "%s.gzi", fa);
    unlink(idx);
    unlink(fa);
    return 0;
}
