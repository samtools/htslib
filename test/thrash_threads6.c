// Spam seeks
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

int main(int argc, char *argv[]) {
    if (argc <= 1) {
	fprintf(stderr, "Usage: thrash_threads4 input.bam\n");
	exit(1);
    }

    // Find a valid seek location ~64M into the file
    int i;
    BGZF *fpin  = bgzf_open(argv[1], "r");
    char buf[100000];
    for (i = 0; i < 100; i++)
	if (bgzf_read(fpin, buf, 65536) < 0)
	    abort();
    int64_t pos = bgzf_tell(fpin);
    while ((i = bgzf_read(fpin, buf, 65536)) > 0)
	continue;
    if (i < 0) abort();
    int64_t end = bgzf_tell(fpin);
    bgzf_close(fpin);

#define N 1000

    // Spam random seeks & reads
    for (i = 0; i < 1000; i++) {
	printf("i=%d\t", i);
	fpin  = bgzf_open(argv[1], "r");
	bgzf_mt(fpin, 8, 256);
	int j, eof = 0;
	for (j = 0; j < 80; j++) {
	    int n = rand() % 6;
	    putchar('0'+n); fflush(stdout);
	    switch (n) {
	    case 0: // start
		if (bgzf_seek(fpin, 0LL, SEEK_SET) < 0) puts("!");//abort();
		eof = 0;
		break;
	    case 1: // mid
		if (bgzf_seek(fpin, pos, SEEK_SET) < 0) puts("!");//abort();
		eof = 0;
		break;
	    case 2: // eof
		if (bgzf_seek(fpin, end, SEEK_SET) < 0) puts("!");//abort();
		eof = 1;
		break;
	    case 3: case 4: {
		int l = rand()%(n==3?100000:100);
		if (bgzf_read(fpin, buf, l) != l*(1-eof)) abort();
		break;
	    }
	    case 5:
		usleep(N);
		break;
	    }
	}
	printf("\n");
	if (bgzf_close(fpin))
	    abort();
    }

    return 0;
}
