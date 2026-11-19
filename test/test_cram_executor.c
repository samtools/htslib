/*  test/test_cram_executor.c -- test cram_set_container_executor().

    Copyright (C) 2026 Nils Homer.

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

/*
 * Writes the same SAM file as CRAM twice -- once with htslib's own encoding,
 * once with a container executor whose worker threads encode containers
 * concurrently while records are still being added -- and checks that both
 * decode to identical records.
 *
 * Usage: test_cram_executor in.sam ref.fa out_prefix
 */

#include <config.h>

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../htslib/cram.h"
#include "../htslib/hts.h"
#include "../htslib/hts_log.h"
#include "../htslib/kstring.h"
#include "../htslib/sam.h"

#define NWORKERS 4
#define SEQS_PER_SLICE 20 // small, so the test file spans many containers

enum state { QUEUED, ENCODING, ENCODED, FAILED };

typedef struct {
    cram_container *c;
    enum state st;
} job_t;

typedef struct {
    cram_fd *fd;
    pthread_mutex_t lock;
    pthread_cond_t cond;
    job_t *jobs;
    size_t njobs, maxjobs, next_write;
    int shutdown;
} executor_t;

static int dispatch(void *arg, cram_fd *fd, cram_container *c) {
    executor_t *ex = arg;
    pthread_mutex_lock(&ex->lock);
    if (ex->njobs == ex->maxjobs) {
        size_t n = ex->maxjobs ? ex->maxjobs * 2 : 16;
        job_t *j = realloc(ex->jobs, n * sizeof(*j));
        if (!j) {
            pthread_mutex_unlock(&ex->lock);
            return -1;
        }
        ex->jobs = j;
        ex->maxjobs = n;
    }
    ex->jobs[ex->njobs].c = c;
    ex->jobs[ex->njobs].st = QUEUED;
    ex->njobs++;
    pthread_cond_broadcast(&ex->cond);
    pthread_mutex_unlock(&ex->lock);
    return 0;
}

static void *worker(void *arg) {
    executor_t *ex = arg;
    pthread_mutex_lock(&ex->lock);
    for (;;) {
        size_t i;
        // Take the newest queued container, so encodes finish out of order.
        for (i = ex->njobs; i > 0 && ex->jobs[i-1].st != QUEUED; i--)
            ;
        if (i == 0) {
            if (ex->shutdown)
                break;
            pthread_cond_wait(&ex->cond, &ex->lock);
            continue;
        }
        job_t *j = &ex->jobs[i-1];
        j->st = ENCODING;
        cram_container *c = j->c;
        pthread_mutex_unlock(&ex->lock);
        int ret = cram_encode_container(ex->fd, c);
        pthread_mutex_lock(&ex->lock);
        ex->jobs[i-1].st = ret == 0 ? ENCODED : FAILED;
        pthread_cond_broadcast(&ex->cond);
    }
    pthread_mutex_unlock(&ex->lock);
    return NULL;
}

// Write every encoded container at the head of the queue, in order.  With
// wait set, keep going until all handed-over containers are written.
static int write_ready(executor_t *ex, int wait) {
    pthread_mutex_lock(&ex->lock);
    for (;;) {
        if (ex->next_write == ex->njobs)
            break;
        job_t *j = &ex->jobs[ex->next_write];
        if (j->st == FAILED) {
            pthread_mutex_unlock(&ex->lock);
            return -1;
        }
        if (j->st != ENCODED) {
            if (!wait)
                break;
            pthread_cond_wait(&ex->cond, &ex->lock);
            continue;
        }
        cram_container *c = j->c;
        ex->next_write++;
        pthread_mutex_unlock(&ex->lock);
        if (cram_flush_encoded_container(ex->fd, c) != 0)
            return -1;
        pthread_mutex_lock(&ex->lock);
    }
    pthread_mutex_unlock(&ex->lock);
    return 0;
}

static int write_cram(const char *in_fn, const char *ref, const char *out_fn,
                      int use_executor, size_t *ncontainers) {
    int ret = -1, i, started = 0;
    pthread_t threads[NWORKERS];
    executor_t ex;
    memset(&ex, 0, sizeof(ex));
    pthread_mutex_init(&ex.lock, NULL);
    pthread_cond_init(&ex.cond, NULL);

    samFile *in = sam_open(in_fn, "r");
    samFile *out = sam_open(out_fn, "wc");
    sam_hdr_t *h = in ? sam_hdr_read(in) : NULL;
    bam1_t *b = bam_init1();
    if (!in || !out || !h || !b)
        goto out;
    if (hts_set_fai_filename(out, ref) < 0
        || hts_set_opt(out, CRAM_OPT_SEQS_PER_SLICE, SEQS_PER_SLICE) < 0)
        goto out;

    if (use_executor) {
        ex.fd = out->fp.cram;
        if (cram_set_container_executor(ex.fd, dispatch, &ex) < 0)
            goto out;
        // A thread pool is refused while an executor is set.
        enum htsLogLevel level = hts_get_log_level();
        hts_set_log_level(HTS_LOG_OFF);
        int refused = hts_set_opt(out, CRAM_OPT_NTHREADS, 2) != 0;
        hts_set_log_level(level);
        if (!refused) {
            fprintf(stderr, "CRAM_OPT_NTHREADS accepted with an executor\n");
            goto out;
        }
        for (i = 0; i < NWORKERS; i++) {
            if (pthread_create(&threads[i], NULL, worker, &ex) != 0)
                goto out;
            started++;
        }
    }

    if (sam_hdr_write(out, h) < 0)
        goto out;
    int r;
    while ((r = sam_read1(in, h, b)) >= 0) {
        if (sam_write1(out, h, b) < 0)
            goto out;
        if (use_executor && write_ready(&ex, 0) < 0)
            goto out;
    }
    if (r < -1)
        goto out;
    if (use_executor) {
        if (cram_flush(ex.fd) < 0 || write_ready(&ex, 1) < 0)
            goto out;
        *ncontainers = ex.njobs;
    }
    ret = 0;

 out:
    if (started) {
        pthread_mutex_lock(&ex.lock);
        ex.shutdown = 1;
        pthread_cond_broadcast(&ex.cond);
        pthread_mutex_unlock(&ex.lock);
        for (i = 0; i < started; i++)
            pthread_join(threads[i], NULL);
    }
    // Containers never written (error paths) are still ours.
    for (size_t k = ex.next_write; k < ex.njobs; k++)
        cram_free_container(ex.jobs[k].c);
    free(ex.jobs);
    if (b) bam_destroy1(b);
    if (h) sam_hdr_destroy(h);
    if (in) sam_close(in);
    if (out && sam_close(out) < 0)
        ret = -1;
    return ret;
}

// Compare two CRAM files record by record (as SAM text).
static int compare(const char *a_fn, const char *b_fn, const char *ref) {
    int ret = -1, ra, rb;
    long n = 0;
    kstring_t sa = KS_INITIALIZE, sb = KS_INITIALIZE;
    samFile *a = sam_open(a_fn, "r"), *b = sam_open(b_fn, "r");
    sam_hdr_t *ha = NULL, *hb = NULL;
    bam1_t *ba = bam_init1(), *bb = bam_init1();
    if (!a || !b || !ba || !bb
        || hts_set_fai_filename(a, ref) < 0 || hts_set_fai_filename(b, ref) < 0
        || !(ha = sam_hdr_read(a)) || !(hb = sam_hdr_read(b)))
        goto out;
    for (;;) {
        ra = sam_read1(a, ha, ba);
        rb = sam_read1(b, hb, bb);
        if (ra < -1 || rb < -1)
            goto out;
        if (ra != rb) {
            fprintf(stderr, "Record counts differ after %ld records\n", n);
            goto out;
        }
        if (ra == -1)
            break;
        if (sam_format1(ha, ba, &sa) < 0 || sam_format1(hb, bb, &sb) < 0)
            goto out;
        if (strcmp(sa.s, sb.s) != 0) {
            fprintf(stderr, "Record %ld differs:\n%s\n%s\n", n, sa.s, sb.s);
            goto out;
        }
        n++;
    }
    if (n == 0) {
        fprintf(stderr, "No records read\n");
        goto out;
    }
    ret = 0;
 out:
    ks_free(&sa);
    ks_free(&sb);
    if (ba) bam_destroy1(ba);
    if (bb) bam_destroy1(bb);
    if (ha) sam_hdr_destroy(ha);
    if (hb) sam_hdr_destroy(hb);
    if (a) sam_close(a);
    if (b) sam_close(b);
    return ret;
}

static int collect(void *arg, cram_fd *fd, cram_container *c) {
    executor_t *ex = arg;
    (void)fd;
    if (ex->njobs == ex->maxjobs) {
        size_t n = ex->maxjobs ? ex->maxjobs * 2 : 16;
        job_t *j = realloc(ex->jobs, n * sizeof(*j));
        if (!j) return -1;
        ex->jobs = j;
        ex->maxjobs = n;
    }
    ex->jobs[ex->njobs++].c = c;
    return 0;
}

// Hand containers to an executor that writes only the first: the close must
// fail and leave no EOF block, so the file cannot pass for complete.
static int abandoned_write(const char *in_fn, const char *ref, const char *out_fn) {
    int ret = -1, closed = 0;
    executor_t ex;
    memset(&ex, 0, sizeof(ex));
    samFile *in = sam_open(in_fn, "r");
    samFile *out = sam_open(out_fn, "wc");
    sam_hdr_t *h = in ? sam_hdr_read(in) : NULL;
    bam1_t *b = bam_init1();
    if (!in || !out || !h || !b
        || hts_set_fai_filename(out, ref) < 0
        || hts_set_opt(out, CRAM_OPT_SEQS_PER_SLICE, SEQS_PER_SLICE) < 0
        || cram_set_container_executor(out->fp.cram, collect, &ex) < 0
        || sam_hdr_write(out, h) < 0)
        goto out;
    while (sam_read1(in, h, b) >= 0)
        if (sam_write1(out, h, b) < 0)
            goto out;
    // Hand over the final container too, so only the unwritten count (not
    // an unflushed container) can stop the close from writing an EOF block.
    if (cram_flush(out->fp.cram) < 0)
        goto out;
    if (ex.njobs < 2) {
        fprintf(stderr, "Expected several containers, got %zu\n", ex.njobs);
        goto out;
    }
    if (cram_encode_container(out->fp.cram, ex.jobs[0].c) < 0
        || cram_flush_encoded_container(out->fp.cram, ex.jobs[0].c) < 0)
        goto out;
    for (size_t k = 1; k < ex.njobs; k++)
        cram_free_container(ex.jobs[k].c);
    ex.njobs = 0;
    enum htsLogLevel level = hts_get_log_level();
    hts_set_log_level(HTS_LOG_OFF);
    int close_ret = sam_close(out);
    hts_set_log_level(level);
    out = NULL;
    closed = 1;
    if (close_ret >= 0) {
        fprintf(stderr, "Closing with unwritten containers succeeded\n");
        goto out;
    }
    samFile *check = sam_open(out_fn, "r");
    if (!check) goto out;
    hts_set_log_level(HTS_LOG_OFF);
    int has_eof = hts_check_EOF(check);
    hts_set_log_level(level);
    sam_close(check);
    if (has_eof != 0) {
        fprintf(stderr, "An abandoned CRAM has an EOF block (hts_check_EOF %d)\n", has_eof);
        goto out;
    }
    ret = 0;
 out:
    for (size_t k = 0; k < ex.njobs; k++)
        cram_free_container(ex.jobs[k].c);
    free(ex.jobs);
    if (b) bam_destroy1(b);
    if (h) sam_hdr_destroy(h);
    if (in) sam_close(in);
    if (out && !closed) sam_close(out);
    return ret;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s in.sam ref.fa out_prefix\n", argv[0]);
        return EXIT_FAILURE;
    }
    kstring_t plain = KS_INITIALIZE, exec = KS_INITIALIZE, gone = KS_INITIALIZE;
    size_t ncontainers = 0;
    int ret = EXIT_FAILURE;
    if (ksprintf(&plain, "%s.plain.cram", argv[3]) < 0
        || ksprintf(&exec, "%s.executor.cram", argv[3]) < 0
        || ksprintf(&gone, "%s.abandoned.cram", argv[3]) < 0)
        goto out;
    if (write_cram(argv[1], argv[2], plain.s, 0, NULL) < 0) {
        fprintf(stderr, "Failed to write %s\n", plain.s);
        goto out;
    }
    if (write_cram(argv[1], argv[2], exec.s, 1, &ncontainers) < 0) {
        fprintf(stderr, "Failed to write %s with an executor\n", exec.s);
        goto out;
    }
    if (ncontainers < 2) {
        fprintf(stderr, "Expected several containers, got %zu\n", ncontainers);
        goto out;
    }
    if (compare(plain.s, exec.s, argv[2]) < 0)
        goto out;
    if (abandoned_write(argv[1], argv[2], gone.s) < 0) {
        fprintf(stderr, "Abandoned write check failed\n");
        goto out;
    }
    ret = EXIT_SUCCESS;
 out:
    ks_free(&plain);
    ks_free(&exec);
    ks_free(&gone);
    return ret;
}
