/// @file htslib/thread_pool.h
/// Thread pool for multi-threading applications.
/*
    Copyright (c) 2013-2017, 2019 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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
 * This file implements a thread pool for multi-threading applications.  It
 * consists of two distinct interfaces: thread pools and thread process
 * queues (a queue of both jobs to-do and of the results of completed jobs).
 * Do not confuse "process" here with a unix PID; rather it is analogous to a
 * program reading a stream of data blocks, processing them in some manner,
 * and outputting a stream of new data blocks.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling except to pick tasks for the
 * processes that have room to store the result.
 *
 * Upon completion, the return value from the function pointer is
 * added to back to the process result queue if required.  We may have
 * multiple "processes" in use for the one pool.
 *
 * To see example usage, please look at the #ifdef TEST_MAIN code in
 * thread_pool.c.
 */

#ifndef HTSLIB_THREAD_POOL_H
#define HTSLIB_THREAD_POOL_H

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 * Opaque data types.
 *
 * Actual definitions are in thread_pool_internal.h, but these should only
 * be used by thread_pool.c itself.
 */

/*
 * An hts_tpool_process implements a queue of input jobs to process and a
 * queue of resulting output  post-processing.  Internally it consists of two
 * buffered queues, analogous to the pipes in a unix pipeline:
 *    ...input | process |  output...
 *
 * Both input and output queues have size limits to prevent either queue from
 * growing too large and serial numbers to ensure sequential consumption of
 * the output.
 *
 * The thread pool may have many heterogeneous tasks, each using its own
 * process mixed into the same thread pool.
 */
typedef struct hts_tpool_process hts_tpool_process;

/*
 * The single pool structure itself.
 *
 * This knows nothing about the nature of the jobs or where their output is
 * going, but it maintains a list of process-queues associated with this pool
 * from which the jobs are taken.
 */
typedef struct hts_tpool hts_tpool;

/*
 * An output, after job has executed.
 */
typedef struct hts_tpool_result hts_tpool_result;


/*-----------------------------------------------------------------------------
 * Thread pool external functions
 */


/*
 * Creates a worker pool with n worker threads.
 *
 * Returns pool pointer on success;
 *         NULL on failure
 *
 * The hts_tpool struct returned by a successful call should be freed
 * via hts_tpool_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
hts_tpool *hts_tpool_init(int n);


/*
 * Returns the number of requested threads for a pool.
 */
HTSLIB_EXPORT
int hts_tpool_size(hts_tpool *p);


/// Add an item to the work pool.
/**
 * @param p     Thread pool
 * @param q     Process queue
 * @param func  Function run by the thread pool
 * @param arg   Data for use by func()
 * @return 0 on success
 *        -1 on failure
 */
// FIXME: should this drop the hts_tpool*p argument? It's just q->p
HTSLIB_EXPORT
int hts_tpool_dispatch(hts_tpool *p, hts_tpool_process *q,
                       void *(*func)(void *arg), void *arg);

/// Add an item to the work pool, with nonblocking option.
/**
 * @param p         Thread pool
 * @param q         Process queue
 * @param func      Function run by the thread pool
 * @param arg       Data for use by func()
 * @param nonblock  Non-blocking flag (see description)
 * @return 0 on success
 *        -1 on failure
 *
 * The @p nonblock parameter can take one of the following values:
 *      0 => block if input queue is full
 *     +1 => don't block if input queue is full, but do not add task
 *     -1 => add task regardless of whether queue is full (over-size)
 *
 * If @p nonblock is +1 and the queue is full, -1 will be returned and
 * `errno` is set to `EAGAIN`.
 */
HTSLIB_EXPORT
int hts_tpool_dispatch2(hts_tpool *p, hts_tpool_process *q,
                        void *(*func)(void *arg), void *arg, int nonblock);

/// Add an item to the work pool, with nonblocking and cleanup callbacks.
/**
 * @param p               Thread pool
 * @param q               Process queue
 * @param exec_func       Function run by the thread pool
 * @param arg             Data for use by func()
 * @param job_cleanup     Callback to clean up when discarding jobs
 * @param result_cleanup  Callback to clean up when discarding result data
 * @param nonblock        Non-blocking flag (see description)
 * @return 0 on success
 *        -1 on failure
 *
 * The @p nonblock parameter can take one of the following values:
 *      0 => block if input queue is full
 *     +1 => don't block if input queue is full, but do not add task
 *     -1 => add task regardless of whether queue is full (over-size)
 *
 * If @p nonblock is +1 and the queue is full, -1 will be returned and
 * `errno` is set to `EAGAIN`.
 *
 * The job_cleanup() and result_cleanup() callbacks are used when discarding
 * data from a queue, for example when calling hts_tpool_process_reset()
 * or hts_tpool_process_destroy().
 *
 * If not NULL, job_cleanup() will be called for each pending job with the
 * value of @p arg that was set for that job.  This can be used to free
 * any data associated with @p arg, and also @p arg itself.
 *
 * Similarly, result_cleanup() can be used to free any results left by
 * jobs that had started before hts_tpool_process_reset() was called.
 * The argument passed to result_cleanup() is the pointer that would
 * have been returned by calling hts_tpool_result_data() on the result
 * when pulled from the queue.
 *
 * job_cleanup() and result_cleanup() are only called when discarding jobs.
 * For jobs that are processed normally, it is the resposibility of
 * exec_func() and / or consumers of any results to do any cleaning up
 * necessary.
 */
HTSLIB_EXPORT
int hts_tpool_dispatch3(hts_tpool *p, hts_tpool_process *q,
                        void *(*exec_func)(void *arg), void *arg,
                        void (*job_cleanup)(void *arg),
                        void (*result_cleanup)(void *data),
                        int nonblock);

/*
 * Wakes up a single thread stuck in dispatch and make it return with
 * errno EAGAIN.
 */
HTSLIB_EXPORT
void hts_tpool_wake_dispatch(hts_tpool_process *q);

/*
 * Flushes the process-queue, but doesn't exit. This simply drains the queue
 * and ensures all worker threads have finished their current tasks
 * associated with this process.
 *
 * NOT: This does not mean the worker threads are not executing jobs in
 * another process-queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
HTSLIB_EXPORT
int hts_tpool_process_flush(hts_tpool_process *q);

/*
 * Resets a process to the initial state.
 *
 * This removes any queued up input jobs, disables any notification of
 * new results/output, flushes what is left and then discards any
 * queued output.  Anything consumer stuck in a wait on results to
 * appear should stay stuck and will only wake up when new data is
 * pushed through the queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
HTSLIB_EXPORT
int hts_tpool_process_reset(hts_tpool_process *q, int free_results);

/* Returns the process queue size */
HTSLIB_EXPORT
int hts_tpool_process_qsize(hts_tpool_process *q);


/*
 * Destroys a thread pool.  The threads are joined into the main
 * thread so they will finish their current work load.
 */
HTSLIB_EXPORT
void hts_tpool_destroy(hts_tpool *p);

/*
 * Destroys a thread pool without waiting on jobs to complete.
 * Use hts_tpool_kill(p) to quickly exit after a fatal error.
 */
HTSLIB_EXPORT
void hts_tpool_kill(hts_tpool *p);

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This doesn't wait for a
 * result to be present.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL if not.
 */
HTSLIB_EXPORT
hts_tpool_result *hts_tpool_next_result(hts_tpool_process *q);

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This will wait for
 * a result to be present if none are currently available.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL on error or during shutdown.
 */
HTSLIB_EXPORT
hts_tpool_result *hts_tpool_next_result_wait(hts_tpool_process *q);

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
HTSLIB_EXPORT
void hts_tpool_delete_result(hts_tpool_result *r, int free_data);

/*
 * Returns the data portion of a hts_tpool_result, corresponding
 * to the actual "result" itself.
 */
HTSLIB_EXPORT
void *hts_tpool_result_data(hts_tpool_result *r);

/*
 * Initialises a thread process-queue.
 *
 * In_only, if true, indicates that the process generates does not need to
 * hold any output.  Otherwise an output queue is used to store the results
 * of processing each input job.
 *
 * Results hts_tpool_process pointer on success;
 *         NULL on failure
 *
 * The hts_tpool_process struct returned by a successful call should be freed
 * via hts_tpool_process_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
hts_tpool_process *hts_tpool_process_init(hts_tpool *p, int qsize, int in_only);


/* Deallocates memory for a thread process-queue.
 * Must be called before the thread pool is destroyed.
 */
HTSLIB_EXPORT
void hts_tpool_process_destroy(hts_tpool_process *q);

/*
 * Returns true if there are no items in the process results queue and
 * also none still pending.
 */
HTSLIB_EXPORT
int hts_tpool_process_empty(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue.
 */
HTSLIB_EXPORT
int hts_tpool_process_len(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue plus the
 * number running and queued up to run.
 */
HTSLIB_EXPORT
int hts_tpool_process_sz(hts_tpool_process *q);

/*
 * Shutdown a process.
 *
 * This sets the shutdown flag and wakes any threads waiting on process
 * condition variables.
 */
HTSLIB_EXPORT
void hts_tpool_process_shutdown(hts_tpool_process *q);

/*
 * Attach and detach a thread process-queue with / from the thread pool
 * scheduler.
 *
 * We need to do attach after making a thread process, but may also wish
 * to temporarily detach if we wish to stop running jobs on a specific
 * process while permitting other process to continue.
 */
HTSLIB_EXPORT
void hts_tpool_process_attach(hts_tpool *p, hts_tpool_process *q);

HTSLIB_EXPORT
void hts_tpool_process_detach(hts_tpool *p, hts_tpool_process *q);

/*
 * Increment and decrement the reference count in a process-queue.
 * If the queue is being driven from two external (non thread-pool)
 * threads, eg "main" and a "reader", this permits each end to
 * decrement its use of the process-queue independently.
 */
HTSLIB_EXPORT
void hts_tpool_process_ref_incr(hts_tpool_process *q);

HTSLIB_EXPORT
void hts_tpool_process_ref_decr(hts_tpool_process *q);

#ifdef __cplusplus
}
#endif

#endif
