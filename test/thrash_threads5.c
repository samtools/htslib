// A basic 'zcat filename [N-threads]'

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

#define N 1000
int main(int argc, char *argv[]) {
    char buf[N];
    ssize_t l, t = 0;
    BGZF *fpin  = bgzf_open(argv[1], "r");
    hts_tpool *p = NULL;
    if (argc > 2) {
        p = hts_tpool_init(atoi(argv[2]));
        bgzf_thread_pool(fpin,  p, 0);
    }
    int n = rand()%(N-1)+1;
    while ((l = bgzf_read(fpin, buf, n)) > 0) {
        if (l != write(1, buf, l)) abort();
	t += l;
        if (l != n) {
            fprintf(stderr, "expected %d bytes, got %d\n", n, (int)l);
            break;
        }
        n = rand()%(N-1)+1;
    }
    fprintf(stderr, "close=%d\n", (int)bgzf_close(fpin));
    if (p) hts_tpool_destroy(p);

    fprintf(stderr, "wrote %d bytes\n", (int)t);

    return 0;
}
