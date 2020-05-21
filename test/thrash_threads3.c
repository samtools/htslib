/* The MIT/Expat License

Copyright (C) 2017 Genome Research Ltd.

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
DEALINGS IN THE SOFTWARE.
 */
// Simple open,read,close thrash.

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "../htslib/bgzf.h"

int main(int argc, char *argv[]) {
    char buf[1000000];
    int i;

    if (argc <= 1) {
        fprintf(stderr, "Usage: thrash_threads3 input.bam\n");
        exit(1);
    }

    for (i = 0; i < 10000; i++) {
        printf("i=%d\n", i);
        BGZF *fpin  = bgzf_open(argv[1], "r");
        if (bgzf_read(fpin, buf, i*10) < 0) abort();
        bgzf_mt(fpin, 8, 256);
        if (bgzf_read(fpin, buf, i*10) < 0) abort();
        if (bgzf_close(fpin) < 0) abort();
    }
    return 0;
}
