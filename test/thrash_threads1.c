// Test extreme rapid turnover of readers, to check for
// race conditions between reader thread launching and file close.

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "htslib/bgzf.h"

int main(int argc, char *argv[]) {
    if (argc <= 1) {
	fprintf(stderr, "Usage: thrash_threads1 input.bam\n");
	exit(1);
    }

    int i;
    for (i = 0; i < 10000; i++) {
	printf("i=%d\n", i);
	BGZF *fpin  = bgzf_open(argv[1], "r");
	bgzf_mt(fpin, 2, 256);
	if (bgzf_close(fpin) < 0) abort();
    }
    return 0;
}
