/*
Copyright (c) 2013 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <config.h>

#include <stdlib.h>

#include <signal.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <stdarg.h>

#include "htslib/thread_pool.h"

//#define DEBUG
//#define DEBUG_TIME

#ifdef DEBUG
static int worker_id(t_pool *p) {
    int i;
    pthread_t s = pthread_self();
    for (i = 0; i < p->tsize; i++) {
	if (pthread_equal(s, p->t[i].tid))
	    return i;
    }
    return -1;
}

int DBG_OUT(FILE *fp, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return vfprintf(fp, fmt, args);
}
#else
#define DBG_OUT(...) do{}while(0)
#endif

/* ----------------------------------------------------------------------------
 * A queue to hold results from the thread pool.
 *
 * Each thread pool may have jobs of multiple types being queued up and
 * interleaved, so we allow several results queue per pool.
 *
 * The jobs themselves are expected to push their results onto their
 * appropriate results queue.
 */

/*
 * Adds a result to the end of the result queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int t_pool_add_result(t_pool_job *j, void *data) {
    t_results_queue *q = j->q;
    t_pool_result *r;

    DBG_OUT(stderr, "%d: Adding result to queue %p, serial %d, %d of %d\n",
	    worker_id(j->p), q, j->serial, q->queue_len+1, q->qsize);

    /* No results queue is fine if we don't want any results back */
    if (!q)
	return 0;

    if (!(r = malloc(sizeof(*r))))
	return -1;

    r->next = NULL;
    r->data = data;
    r->serial = j->serial;

    pthread_mutex_lock(&q->p->pool_m);

    if (q->result_tail) {
	q->result_tail->next = r;
	q->result_tail = r;
    } else {
	q->result_head = q->result_tail = r;
    }
    q->queue_len++;
    q->pending--;

    DBG_OUT(stderr, "%d: Broadcasting result_avail (id %d)\n",
	    worker_id(j->p), r->serial);
    pthread_cond_signal(&q->result_avail_c);
    DBG_OUT(stderr, "%d: Broadcast complete\n", worker_id(j->p));

    pthread_mutex_unlock(&q->p->pool_m);

    return 0;
}

/*
 * Returns true or false depending on whether the results queue is full.
 */
int t_pool_results_queue_full(t_results_queue *q) {
    int full;

    pthread_mutex_lock(&q->p->pool_m);
    full = (q->qsize && q->queue_len >= q->qsize);
    pthread_mutex_unlock(&q->p->pool_m);
    return full;
}

static void wake_next_worker(t_pool *p, int locked);

/* Core of t_pool_next_result() */
static t_pool_result *t_pool_next_result_locked(t_results_queue *q) {
    t_pool_result *r, *last;

    for (last = NULL, r = q->result_head; r; last = r, r = r->next) {
	if (r->serial == q->next_serial)
	    break;
    }

    if (r) {
	if (q->result_head == r)
	    q->result_head = r->next;
	else
	    last->next = r->next;

	if (q->result_tail == r)
	    q->result_tail = last;

	if (!q->result_head)
	    q->result_tail = NULL;

	q->next_serial++;
	q->queue_len--;

	if (q->qsize && q->queue_len < q->qsize) {
	    pthread_cond_signal(&q->p->full_c);
	    wake_next_worker(q->p, 1);
	}
    }

    return r;
}

/*
 * Pulls a result off the head of the result queue. Caller should
 * free it (and any internals as appropriate) after use. This doesn't
 * wait for a result to be present.
 *
 * Results will be returned in strict order.
 * 
 * Returns t_pool_result pointer if a result is ready.
 *         NULL if not.
 */
t_pool_result *t_pool_next_result(t_results_queue *q) {
    t_pool_result *r;

    DBG_OUT(stderr, "Requesting next result on queue %p\n", q);

    pthread_mutex_lock(&q->p->pool_m);
    r = t_pool_next_result_locked(q);
    pthread_mutex_unlock(&q->p->pool_m);

    DBG_OUT(stderr, "(q=%p) Found %p\n", q, r);

    return r;
}

t_pool_result *t_pool_next_result_wait(t_results_queue *q) {
    t_pool_result *r;

    pthread_mutex_lock(&q->p->pool_m);
    while (!(r = t_pool_next_result_locked(q))) {
	/* Possible race here now avoided via _locked() call, but incase... */
	struct timeval now;
	struct timespec timeout;

	if (q->shutdown) {
	    pthread_mutex_unlock(&q->p->pool_m);
	    return NULL;
	}

	gettimeofday(&now, NULL);
	timeout.tv_sec = now.tv_sec + 10;
	timeout.tv_nsec = now.tv_usec * 1000;

	pthread_cond_timedwait(&q->result_avail_c, &q->p->pool_m, &timeout);

	if (q->shutdown) {
	    pthread_mutex_unlock(&q->p->pool_m);
	    return NULL;
	}
    }
    pthread_mutex_unlock(&q->p->pool_m);

    return r;
}

/*
 * Returns true if there are no items on the finished results queue and
 * also none still pending.
 */
int t_pool_results_queue_empty(t_results_queue *q) {
    int empty;

    pthread_mutex_lock(&q->p->pool_m);
    empty = q->queue_len == 0 && q->pending == 0;
    pthread_mutex_unlock(&q->p->pool_m);

    return empty;
}

/*
 * Returns the number of completed jobs on the results queue.
 */
int t_pool_results_queue_len(t_results_queue *q) {
    int len;

    pthread_mutex_lock(&q->p->pool_m);
    len = q->queue_len;
    pthread_mutex_unlock(&q->p->pool_m);

    return len;
}

int t_pool_results_queue_sz(t_results_queue *q) {
    int len;

    pthread_mutex_lock(&q->p->pool_m);
    len = q->queue_len + q->pending;
    pthread_mutex_unlock(&q->p->pool_m);

    return len;
}

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
void t_pool_delete_result(t_pool_result *r, int free_data) {
    if (!r)
	return;

    if (free_data && r->data)
	free(r->data);

    free(r);
}

/*
 * Initialises a results queue.
 *
 * Results queue pointer on success;
 *         NULL on failure
 */
t_results_queue *t_results_queue_init(t_pool *p, int qsize) {
    t_results_queue *q = malloc(sizeof(*q));

    pthread_cond_init(&q->result_avail_c, NULL);

    q->p           = p;
    q->result_head = NULL;
    q->result_tail = NULL;
    q->next_serial = 0;
    q->curr_serial = 0;
    q->queue_len   = 0;
    q->pending     = 0;
    q->qsize       = qsize;
    q->shutdown    = 0;

    return q;
}

/* Sends a shutdown signal to the results queue, preventing
 * new data from appearing and permitting draining of results.
 */
void t_results_queue_shutdown(t_results_queue *q) {
    pthread_mutex_lock(&q->p->pool_m);
    pthread_cond_signal(&q->result_avail_c);
    q->shutdown = 1;
    pthread_mutex_unlock(&q->p->pool_m);
}

/* Deallocates memory for a results queue.
 * Must be called before the thread pool is destroyed.
 */
void t_results_queue_destroy(t_results_queue *q) {
    DBG_OUT(stderr, "Destroying results queue %p\n", q);

    if (!q)
	return;

    pthread_cond_destroy(&q->result_avail_c);

    memset(q, 0xbb, sizeof(*q));
    free(q);

    DBG_OUT(stderr, "Destroyed results queue %p\n", q);
}

/* ----------------------------------------------------------------------------
 * The thread pool.
 */

#define TDIFF(t2,t1) ((t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec)

/*
 * A worker thread.
 *
 * Each thread waits for the pool to be non-empty.
 * As soon as this applies, one of them succeeds in getting the lock
 * and then executes the job.
 */
static void *t_pool_worker(void *arg) {
    t_pool_worker_t *w = (t_pool_worker_t *)arg;
    t_pool *p = w->p;
    t_pool_job *j;
#ifdef DEBUG_TIME
    struct timeval t1, t2, t3;
#endif

    for (;;) {
	// Pop an item off the pool queue
#ifdef DEBUG_TIME
	gettimeofday(&t1, NULL);
#endif

	pthread_mutex_lock(&p->pool_m);

#ifdef DEBUG_TIME
	gettimeofday(&t2, NULL);
	p->wait_time += TDIFF(t2,t1);
	w->wait_time += TDIFF(t2,t1);
#endif

	// If there is something on the job list and a higher priority
	// thread waiting, let it handle this instead.
//	while (p->head && p->t_stack_top != -1 && p->t_stack_top < w->idx) {
//	    pthread_mutex_unlock(&p->pool_m);
//	    pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
//	    pthread_mutex_lock(&p->pool_m);
//	}

	while ((!p->head && !p->shutdown)
	       || (p->head && t_pool_results_queue_full(p->head->q))) {

	    // Push beyond the queue somewhat as it's been shown that once we
	    // start running it's benefifical to keep running, so a different queue
	    // size for stopping vs starting works well.  Needed? Sufficient?
	    if (p->head && p->head->q->qsize &&
		p->head->serial < p->head->q->next_serial + p->head->q->qsize*2)
		break;

	    p->nwaiting++;

	    if (p->njobs == 0)
		pthread_cond_signal(&p->empty_c);
#ifdef DEBUG_TIME
	    gettimeofday(&t2, NULL);
#endif

	    // Push this thread to the top of the waiting stack
	    if (p->t_stack_top == -1 || p->t_stack_top > w->idx)
		p->t_stack_top = w->idx;

	    p->t_stack[w->idx] = 1;
	    pthread_cond_wait(&w->pending_c, &p->pool_m);
	    p->t_stack[w->idx] = 0;

	    /* Find new t_stack_top */
	    {
		int i;
		p->t_stack_top = -1;
		for (i = 0; i < p->tsize; i++) {
		    if (p->t_stack[i]) {
			p->t_stack_top = i;
			break;
		    }
		}
	    }

#ifdef DEBUG_TIME
	    gettimeofday(&t3, NULL);
	    p->wait_time += TDIFF(t3,t2);
	    w->wait_time += TDIFF(t3,t2);
#endif
	    p->nwaiting--;
	}

	if (p->shutdown) {
#ifdef DEBUG_TIME
	    p->total_time += TDIFF(t3,t1);
#endif
#ifdef DEBUG
	    fprintf(stderr, "%d: Shutting down\n", worker_id(p));
#endif
	    pthread_mutex_unlock(&p->pool_m);
	    pthread_exit(NULL);
	}

	j = p->head;
	if (!(p->head = j->next))
	    p->tail = NULL;

	if (p->njobs-- >= p->qsize)
	    pthread_cond_signal(&p->full_c);

	if (p->njobs == 0)
	    pthread_cond_signal(&p->empty_c);

	pthread_mutex_unlock(&p->pool_m);
	    
	// We have job 'j' - now execute it.
	t_pool_add_result(j, j->func(j->arg));	
#ifdef DEBUG_TIME
	pthread_mutex_lock(&p->pool_m);
	gettimeofday(&t3, NULL);
	p->total_time += TDIFF(t3,t1);
	pthread_mutex_unlock(&p->pool_m);
#endif
	memset(j, 0xbb, sizeof(*j));
	free(j);
    }

    return NULL;
}

static void wake_next_worker(t_pool *p, int locked) {
    if (!locked)
	pthread_mutex_lock(&p->pool_m);

    /*
    ERR!

    && p->head->q->qsize - p->head->q->queue_len >= p->njobs
       64 - 128 >= 8

    Why is queue_len so large? That's -64 >= 8

Basically we've got a lot of data in the output queue, too much to
process (64 per encoder / decoder) so aren't going to wake up anything
more.

However nothing is reading from the decode output queue; the main
thread is stuck trying to dispatch an encode job.
    */


    if (p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting
	&& p->head->q->qsize - p->head->q->queue_len >= p->njobs)
	pthread_cond_signal(&p->t[p->t_stack_top].pending_c);

    if (!locked)
	pthread_mutex_unlock(&p->pool_m);
}

/*
 * Creates a worker pool of length qsize with tsize worker threads.
 *
 * Returns pool pointer on success;
 *         NULL on failure
 */
t_pool *t_pool_init(int qsize, int tsize) {
    int i;
    t_pool *p = malloc(sizeof(*p));
    p->qsize = qsize;
    p->tsize = tsize;
    p->njobs = 0;
    p->nwaiting = 0;
    p->shutdown = 0;
    p->head = p->tail = NULL;
    p->t_stack = NULL;
#ifdef DEBUG_TIME
    p->total_time = p->wait_time = 0;
#endif

    p->t = malloc(tsize * sizeof(p->t[0]));

    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE_NP);
    
    pthread_mutex_init(&p->pool_m, &attr);
    pthread_cond_init(&p->empty_c, NULL);
    pthread_cond_init(&p->full_c, NULL);

    pthread_mutexattr_destroy(&attr);

    pthread_mutex_lock(&p->pool_m);

    if (!(p->t_stack = malloc(tsize * sizeof(*p->t_stack))))
	return NULL;
    p->t_stack_top = -1;

    for (i = 0; i < tsize; i++) {
	t_pool_worker_t *w = &p->t[i];
	p->t_stack[i] = 0;
	w->p = p;
	w->idx = i;
	w->wait_time = 0;
	pthread_cond_init(&w->pending_c, NULL);
	if (0 != pthread_create(&w->tid, NULL, t_pool_worker, w))
	    return NULL;
    }

    pthread_mutex_unlock(&p->pool_m);

    return p;
}

/*
 * Adds an item to the work pool.
 *
 * FIXME: Maybe return 1,0,-1 and distinguish between job dispathed vs
 * result returned. Ie rather than blocking on full queue we're permitted
 * to return early on "result available" event too.
 * Caller would then have a while loop around t_pool_dispatch.
 * Or, return -1 and set errno to EAGAIN to indicate job not yet submitted.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int t_pool_dispatch(t_pool *p, t_results_queue *q,
		    void *(*func)(void *arg), void *arg) {
    return t_pool_dispatch2(p, q, func, arg, 0);
}

/*
 * As above but optional non-block flag.
 *
 * nonblock  0 => block if input queue is full
 * nonblock +1 => don't block if input queue is full, but do not add task
 * nonblock -1 => add task regardless of whether queue is full (over-size)
 */
int t_pool_dispatch2(t_pool *p, t_results_queue *q,
		     void *(*func)(void *arg), void *arg, int nonblock) {
    t_pool_job *j;

    DBG_OUT(stderr, "Dispatching job for queue %p, serial %d\n", q, q->curr_serial);

    pthread_mutex_lock(&p->pool_m);

    if (p->njobs >= p->qsize && nonblock == 1) {
	pthread_mutex_unlock(&p->pool_m);
	errno = EAGAIN;
	return -1;
    }

    if (!(j = malloc(sizeof(*j))))
	return -1;
    j->func = func;
    j->arg = arg;
    j->next = NULL;
    j->p = p;
    j->q = q;
    if (q) {
	j->serial = q->curr_serial++;
	q->pending++;
    } else {
	j->serial = 0;
    }

    if (nonblock == 0) {
	while (p->njobs >= p->qsize)
	    pthread_cond_wait(&p->full_c, &p->pool_m);
    }

    p->njobs++;
    
//    if (q->curr_serial % 100 == 0)
//	fprintf(stderr, "p->njobs = %d    p->qsize = %d\n", p->njobs, p->qsize);

    if (p->tail) {
	p->tail->next = j;
	p->tail = j;
    } else {
	p->head = p->tail = j;
    }

    DBG_OUT(stderr, "Dispatched (serial %d)\n", j->serial);

    // Let a worker know we have data.
    // Keep incoming queue at 1 per running thread, so there is always
    // something waiting when they end their current task.  If we go above
    // this signal to start more threads (if available). This has the effect
    // of concentrating jobs to fewer cores when we are I/O bound, which in
    // turn benefits systems with auto CPU frequency scaling.
    wake_next_worker(p, 1);

    pthread_mutex_unlock(&p->pool_m);

    return 0;
}

/*
 * Flushes the pool, but doesn't exit. This simply drains the queue and
 * ensures all worker threads have finished their current task.
 *
 * NOTE: we may need to have one input queue per type of job, meaning we
 * keep track of the total number of jobs per queue and flush is a wait on the
 * output number matching the input number (with none left in the input).
 * (Currently we cannot flush, say, BAM encode jobs without also flushing
 * BAM decode jobs.)
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int t_pool_flush(t_pool *p) {
    int i;

    DBG_OUT(stderr, "Flushing pool %p\n", p);

    // Drains the queue
    pthread_mutex_lock(&p->pool_m);

    // Wake up everything for the final sprint!
    for (i = 0; i < p->tsize; i++)
	if (p->t_stack[i])
	    pthread_cond_signal(&p->t[i].pending_c);

    while (p->njobs || p->nwaiting != p->tsize)
	pthread_cond_wait(&p->empty_c, &p->pool_m);

    pthread_mutex_unlock(&p->pool_m);

    DBG_OUT(stderr, "Flushed complete for pool %p, njobs=%d, nwaiting=%d\n",
	    p, p->njobs, p->nwaiting);

    return 0;
}


/*
 * Destroys a thread pool. If 'kill' is true the threads are terminated now,
 * otherwise they are joined into the main thread so they will finish their
 * current work load.
 *
 * Use t_pool_destroy(p,0) after a t_pool_flush(p) on a normal shutdown or
 * t_pool_destroy(p,1) to quickly exit after a fatal error.
 */
void t_pool_destroy(t_pool *p, int kill) {
    int i;
    
    DBG_OUT(stderr, "Destroying pool %p, kill=%d\n", p, kill);

    /* Send shutdown message to worker threads */
    if (!kill) {
	pthread_mutex_lock(&p->pool_m);
	p->shutdown = 1;

	DBG_OUT(stderr, "Sending shutdown request\n");

	for (i = 0; i < p->tsize; i++)
	    pthread_cond_signal(&p->t[i].pending_c);

	pthread_mutex_unlock(&p->pool_m);

	DBG_OUT(stderr, "Shutdown complete\n");

	for (i = 0; i < p->tsize; i++)
	    pthread_join(p->t[i].tid, NULL);
    } else {
	for (i = 0; i < p->tsize; i++)
	    pthread_kill(p->t[i].tid, SIGINT);
    }

    pthread_mutex_destroy(&p->pool_m);
    pthread_cond_destroy(&p->empty_c);
    pthread_cond_destroy(&p->full_c);
    for (i = 0; i < p->tsize; i++)
	pthread_cond_destroy(&p->t[i].pending_c);

#ifdef DEBUG_TIME
    fprintf(stderr, "Total time=%f\n", p->total_time / 1000000.0);
    fprintf(stderr, "Wait  time=%f\n", p->wait_time  / 1000000.0);
    fprintf(stderr, "%d%% utilisation\n",
	    (int)(100 - ((100.0 * p->wait_time) / p->total_time + 0.5)));
    for (i = 0; i < p->tsize; i++)
	fprintf(stderr, "%d: Wait time=%f\n", i,
		p->t[i].wait_time / 1000000.0);
#endif

    if (p->t_stack)
	free(p->t_stack);

    free(p->t);
    free(p);

    DBG_OUT(stderr, "Destroyed pool %p\n", p);
}


/*-----------------------------------------------------------------------------
 * Test app.
 */

#ifdef TEST_MAIN

#include <stdio.h>
#include <math.h>

void *doit(void *arg) {
    int i, k, x = 0;
    int job = *(int *)arg;
    int *res;

    printf("Worker: execute job %d\n", job);

    usleep(random() % 1000000); // to coerce job completion out of order
    if (0) {
	for (k = 0; k < 100; k++) {
	    for (i = 0; i < 100000; i++) {
		x++;
		x += x * sin(i);
		x += x * cos(x);
	    }
	}
	x *= 100;
	x += job;
    } else {
	x = job*job;
    }

    printf("Worker: job %d terminating, x=%d\n", job, x);

    free(arg);

    res = malloc(sizeof(*res));
    *res = x;

    return res;
}

#define NTHREADS 8

int main(int argc, char **argv) {
    t_pool *p = t_pool_init(NTHREADS*2, NTHREADS);
    t_results_queue *q = t_results_queue_init();
    int i;
    t_pool_result *r;

    // Dispatch jobs
    for (i = 0; i < 20; i++) {
	int *ip = malloc(sizeof(*ip));
	*ip = i;
	printf("Submitting %d\n", i);
	t_pool_dispatch(p, q, doit, ip);
	
	// Check for results
	if ((r = t_pool_next_result(q))) {
	    printf("RESULT: %d\n", *(int *)r->data);
	    t_pool_delete_result(r, 1);
	}
    }

    t_pool_flush(p);

    while ((r = t_pool_next_result(q))) {
	printf("RESULT: %d\n", *(int *)r->data);
	t_pool_delete_result(r, 1);
    }

    t_pool_destroy(p, 0);
    t_results_queue_destroy(q);

    return 0;
}
#endif
