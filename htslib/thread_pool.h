/*
Copyright (c) 2013-2016 Genome Research Ltd.
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

#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <pthread.h>
#include <inttypes.h>

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
 * The thread pool may have many hetergeneous tasks, each using its own
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
 */
hts_tpool *hts_tpool_init(int n);


/*
 * Returns the number of requested threads for a pool.
 */
int hts_tpool_size(hts_tpool *p);


/*
 * Adds an item to the work pool.
 *
 * FIXME: permit q to be NULL, indicating a global/default pool held by
 * the thread pool itself?  This pool would be for jobs that have no
 * output, so fire and forget only with..
 *
 * Returns 0 on success
 *        -1 on failure
 */

// FIXME: should this drop the hts_tpool*p argument? It's just q->p
int hts_tpool_dispatch(hts_tpool *p, hts_tpool_process *q,
                       void *(*func)(void *arg), void *arg);
int hts_tpool_dispatch2(hts_tpool *p, hts_tpool_process *q,
                        void *(*func)(void *arg), void *arg, int nonblock);

/*
 * Wakes up a single thread stuck in dispatch and make it return with
 * errno EAGAIN.
 */
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
int hts_tpool_process_flush(hts_tpool_process *q);

/*
 * Resets a process to the intial state.
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
int hts_tpool_process_reset(hts_tpool_process *q, int free_results);

/* Returns the process queue size */
int hts_tpool_process_qsize(hts_tpool_process *q);


/*
 * Destroys a thread pool.  The threads are joined into the main
 * thread so they will finish their current work load. 
 */
void hts_tpool_destroy(hts_tpool *p);

/*
 * Destroys a thread pool without waiting on jobs to complete.
 * Use hts_tpool_kill(p) to quickly exit after a fatal error.
 */
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
hts_tpool_result *hts_tpool_next_result_wait(hts_tpool_process *q);

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
void hts_tpool_delete_result(hts_tpool_result *r, int free_data);

/*
 * Returns the data portion of a hts_tpool_result, corresponding
 * to the actual "result" itself.
 */
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
 */
hts_tpool_process *hts_tpool_process_init(hts_tpool *p, int qsize, int in_only);


/* Deallocates memory for a thread process-queue.
 * Must be called before the thread pool is destroyed.
 */
void hts_tpool_process_destroy(hts_tpool_process *q);

/*
 * Flushes the thread pool, but doesn't exit. This simply drains the
 * process-queue and ensures all worker threads have finished their current
 * task if associated with this process.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int hts_tpool_process_flush(hts_tpool_process *q);

/*
 * Returns true if there are no items in the process results queue and
 * also none still pending.
 */
int hts_tpool_process_empty(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue.
 */
int hts_tpool_process_len(hts_tpool_process *q);

/*
 * Returns the number of completed jobs in the process results queue plus the
 * number running and queued up to run.
 */
int hts_tpool_process_sz(hts_tpool_process *q);

/*
 * Shutdown a process.
 *
 * This sets the shutdown flag and wakes any threads waiting on process
 * condition variables.
 */
void hts_tpool_process_shutdown(hts_tpool_process *q);

/*
 * Attach and detach a thread process-queue with / from the thread pool
 * scheduler.
 *
 * We need to do attach after making a thread process, but may also wish
 * to temporarily detach if we wish to stop running jobs on a specific
 * process while permitting other process to continue.
 */
void hts_tpool_process_attach(hts_tpool *p, hts_tpool_process *q);
void hts_tpool_process_detach(hts_tpool *p, hts_tpool_process *q);

#ifdef __cplusplus
}
#endif

#endif /* _THREAD_POOL_H_ */
