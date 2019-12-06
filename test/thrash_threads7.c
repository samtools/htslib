/* The MIT/Expat License

Copyright (C) 2017-2018 Genome Research Ltd.

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
/*
 * Test for thread lock-ups caused by a race condition on the queue list
 * where the process tpool_worker is working on could get detached just
 * after it finished running a job.  This would result on the pointer
 * to the next process to be searched for work being set to NULL, which
 * stopped all the workers from finding anything to do.
 */


#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <errno.h>
#include "htslib/thread_pool.h"


void *job(void *v) {
    unsigned int *usecs = (unsigned int *) v;
    usleep(*usecs);
    return NULL;
}

int main(int argc, char *argv[]) {
    int run_for_secs = 120;
    int num_threads = 8;
    int num_jobs = 8, count = 0, n_proc = 8, i;
    struct timeval end, now;
    hts_tpool *p = NULL;
    hts_tpool_process *q[n_proc];

    p = hts_tpool_init(num_threads);
    if (!p) {
        perror("hts_tpool_init");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n_proc; i++) {
        q[i] = hts_tpool_process_init(p, 10, 1);
        if (!q[i]) {
            perror("hts_tpool_process_init");
            exit(EXIT_FAILURE);
        }
    }

    if (gettimeofday(&end, NULL) != 0) {
        perror("gettimeofday");
        exit(EXIT_FAILURE);
    }

    end.tv_sec += run_for_secs;

    do {
        unsigned int *t;
        int qnum = rand() % n_proc;
        t = malloc(num_jobs * sizeof(*t));
        if (!t) {
            perror("malloc");
            exit(EXIT_FAILURE);
        }
        if ((count++ & 15) == 0) {
            fprintf(stderr, "\r%d ", count);
            alarm(10);
        }
        for (i = 0; i < num_jobs; i++) {
            t[i] = 1000;
            if (hts_tpool_dispatch(p, q[qnum], job, &t[i]) < 0) {
                perror("hts_tpool_dispatch");
                exit(EXIT_FAILURE);
            }
        }
        hts_tpool_process_flush(q[qnum]);
        hts_tpool_process_destroy(q[qnum]);
        free(t);
        q[qnum] = hts_tpool_process_init(p, 10, 1);
        if (!q[qnum]) {
            perror("hts_tpool_process_init");
            exit(EXIT_FAILURE);
        }

        if (gettimeofday(&now, NULL) != 0) {
            perror("gettimeofday");
            exit(EXIT_FAILURE);
        }
    } while (now.tv_sec < end.tv_sec
             || (now.tv_sec == end.tv_sec && now.tv_usec < end.tv_usec));
    for (i = 0; i < n_proc; i++) {
        hts_tpool_process_flush(q[i]);
        hts_tpool_process_destroy(q[i]);
    }
    hts_tpool_destroy(p);
    fprintf(stderr, "\n");

    return EXIT_SUCCESS;
}
