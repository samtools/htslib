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
    char buf[65536];
    for (i = 0; i < 1000; i++)
	if (bgzf_read(fpin, buf, 65536) < 0)
	    abort();
    int64_t pos = bgzf_tell(fpin);
    bgzf_close(fpin);

#define N 1000

    // Spam seeks
    for (i = 0; i < 1000; i++) {
	printf("i=%d\n", i);
	fpin  = bgzf_open(argv[1], "r");
	bgzf_mt(fpin, 8, 256);
	if (bgzf_seek(fpin, pos, SEEK_SET) < 0) abort();
	usleep(N);
	//if (bgzf_read(fpin, buf, 65536) < 0) abort();
	//write(1, buf, 65536);
	if (bgzf_seek(fpin, 0LL, SEEK_SET) < 0) abort();
	usleep(N);
	//if (bgzf_read(fpin, buf, 65536) < 0) abort();
	//write(1, buf, 65536);
	if (bgzf_close(fpin))
	    abort();
    }

    return 0;
}
