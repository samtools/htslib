// Test extreme rapid turnover of writers, to check for
// race conditions between reader thread launching and file close.

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

int main(int argc, char *argv[]) {
    int i;
    for (i = 0; i < 1000; i++) {
	printf("i=%d\n", i);
	BGZF *fp  = bgzf_open("/dev/null", "w");
	bgzf_mt(fp, 8, 256);
	if (bgzf_close(fp))
	    abort();
    }

    return 0;
}
