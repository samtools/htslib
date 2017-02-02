// Simple open,read,close thrash.

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "htslib/bgzf.h"

int main(int argc, char *argv[]) {
    char buf[1000];
    int i;

    if (argc <= 1) {
	fprintf(stderr, "Usage: thrash_threads3 input.bam\n");
	exit(1);
    }

    for (i = 0; i < 10000; i++) {
	printf("i=%d\n", i);
	BGZF *fpin  = bgzf_open(argv[1], "r");
	bgzf_mt(fpin, 8, 256);
	if (bgzf_read(fpin, buf, 1000) < 0) abort();
	if (bgzf_close(fpin) < 0) abort();
    }
    return 0;
}
