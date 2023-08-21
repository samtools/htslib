/*  qtask_ordered.c --  showcases the htslib api usage

    Copyright (C) 2024 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE

*/

/* The purpose of this code is to demonstrate the library apis and need proper error handling and optimisation */

#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

typedef struct data {
    int count;                  //used up size
    int maxsize;                //max size per data chunk
    bam1_t **bamarray;          //bam1_t array for optimal queueing
    struct data *next;          //pointer to next one - to reuse earlier allocations
} data;

typedef struct datacache
{
    pthread_mutex_t lock;       //synchronizes the access to cache
    data *list;                 //data storage
} datacache;

typedef struct orderedwrite {
    samFile *outfile;           //output file handle
    sam_hdr_t *samhdr;          //header used to write data
    hts_tpool_process *queue;   //queue from which results to be retrieved
    datacache *cache;           //to re-use allocated storage
    int result;                 //result code returned by writer thread
} orderedwrite;

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: qtask_ordered infile threadcount outdir [chunksize]\n\
Calculates GC ratio - sum(G,C) / sum(A,T,C,G) - and adds to each alignment\n\
as xr:f aux tag. Output is saved in outdir.\n\
chunksize [4096] sets the number of alignments clubbed together to process.\n");
    return;
}

/// getbamstorage - allocates storage for alignments to queue
/** @param chunk number of bam data to allocate
 * @param bamcache cached storage
returns already allocated data storage if one is available, otherwise allocates new
*/
data* getbamstorage(int chunk, datacache *bamcache)
{
    int i = 0;
    data *bamdata = NULL;

    if (!bamcache) {
        return NULL;
    }
    //get from cache if there is an already allocated storage
    if (pthread_mutex_lock(&bamcache->lock)) {
        return NULL;
    }
    if (bamcache->list) {                   //available
        bamdata = bamcache->list;
        bamcache->list = bamdata->next;     //remove and set next one as available
        bamdata->next = NULL;               //remove link
        bamdata->count = 0;
        goto end;
    }
    //allocate and use
    if (!(bamdata = malloc(sizeof(data)))) {
        goto end;
    }
    bamdata->bamarray = malloc(chunk * sizeof(bam1_t*));
    if (!bamdata->bamarray) {
        free(bamdata);
        bamdata = NULL;
        goto end;
    }
    for (i = 0; i < chunk; ++i) {
        bamdata->bamarray[i] = bam_init1();
    }
    bamdata->maxsize = chunk;
    bamdata->count = 0;
    bamdata->next = NULL;

end:
    pthread_mutex_unlock(&bamcache->lock);
    return bamdata;
}

/// cleanup_bamstorage - frees a bamdata struct plus contents
/** @param arg Pointer to data to free
    @p arg has type void * so it can be used as a callback passed
    to hts_tpool_dispatch3().
 */
void cleanup_bamstorage(void *arg)
{
    data *bamdata = (data *) arg;
    if (!bamdata)
        return;
    if (bamdata->bamarray) {
        int i;
        for (i = 0; i < bamdata->maxsize; i++) {
            bam_destroy1(bamdata->bamarray[i]);
        }
        free(bamdata->bamarray);
    }
    free(bamdata);
}

/// thread_ordered_proc - does the processing of task in queue and queues the output back
/** @param args pointer to set of data to be processed
returns the processed data
the processing could be in any order based on the number of threads in use but read of output
from queue will be in order
a null data indicates the end of input and a null is returned to be added back to result queue
*/
void *thread_ordered_proc(void *args)
{
    int i = 0, pos = 0;
    data *bamdata = (data*)args;
    float gcratio = 0;
    uint8_t *data = NULL;

    if (bamdata == NULL)
        return NULL; // Indicates no more input

    for ( i = 0; i < bamdata->count; ++i) {
        //add count
        uint64_t count[16] = {0};
        data = bam_get_seq(bamdata->bamarray[i]);
        for (pos = 0; pos < bamdata->bamarray[i]->core.l_qseq; ++pos) {
            count[bam_seqi(data,pos)]++;
        }
        /*it is faster to count all and use offset to get required counts rather than select
        require ones inside the loop*/
        gcratio = (count[2] /*C*/ + count[4] /*G*/) / (float) (count[1] /*A*/ + count[8] /*T*/ + count[2] + count[4]);

        if (bam_aux_append(bamdata->bamarray[i], "xr", 'f', sizeof(gcratio), (const uint8_t*)&gcratio) < 0) {
            fprintf(stderr, "Failed to add aux tag xr, errno: %d\n", errno);
            break;
        }
    }
    return bamdata;
}

/// threadfn_orderedwrite - thread that read the output from queue and writes
/** @param args pointer to data specific for the thread
returns NULL
*/
void *threadfn_orderedwrite(void *args)
{
    orderedwrite *tdata = (orderedwrite*)args;
    hts_tpool_result *r = NULL;
    data *bamdata = NULL;
    int i = 0;

    tdata->result = 0;

    //get result and write; wait if no result is in queue - until shutdown of queue
    while (tdata->result == 0 &&
           (r = hts_tpool_next_result_wait(tdata->queue)) != NULL) {
        bamdata = (data*) hts_tpool_result_data(r);

        if (bamdata == NULL) {
            // Indicator for no more input. Time to stop.
            hts_tpool_delete_result(r, 0);
            break;
        }

        for (i = 0; i < bamdata->count; ++i) {
            if (sam_write1(tdata->outfile, tdata->samhdr, bamdata->bamarray[i]) < 0) {
                fprintf(stderr, "Failed to write output data\n");
                tdata->result = -1;
                break;
            }
        }
        hts_tpool_delete_result(r, 0);              //release the result memory

        pthread_mutex_lock(&tdata->cache->lock);
        bamdata->next = tdata->cache->list;         //make current list as next
        tdata->cache->list = bamdata;               //set as current to reuse
        pthread_mutex_unlock(&tdata->cache->lock);
    }

    // Shut down the process queue.  If we stopped early due to a write failure,
    // this will signal to the other end that something has gone wrong.
    hts_tpool_process_shutdown(tdata->queue);

    return NULL;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *outdir = NULL;
    char *file = NULL;
    int c = 0, ret = EXIT_FAILURE, cnt = 0, started_thread = 0, chunk = 0;
    size_t size = 0;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    pthread_t thread;
    orderedwrite twritedata = {0};
    hts_tpool *pool = NULL;
    hts_tpool_process *queue = NULL;
    htsThreadPool tpool = {NULL, 0};
    data *bamdata = NULL;
    datacache bamcache = {PTHREAD_MUTEX_INITIALIZER, NULL};

    //qtask infile threadcount outdir [chunksize]
    if (argc != 4 && argc != 5) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    cnt = atoi(argv[2]);
    outdir = argv[3];
    if (argc == 5) {    //chunk size present
        chunk = atoi(argv[4]);
    }
    if (cnt < 1) {      //set proper thread count
        cnt = 1;
    }
    if (chunk < 1) {    //set valid  chunk size
        chunk = 4096;
    }

    //allocate space for output
    size = (strlen(outdir) + sizeof("/out.bam") + 1);   //space for output file name and null termination
    if (!(file = malloc(size))) {
        fprintf(stderr, "Failed to set output path\n");
        goto end;
    }
    snprintf(file, size, "%s/out.bam", outdir);         //output file name
    if (!(pool = hts_tpool_init(cnt))) {                //thread pool
        fprintf(stderr, "Failed to create thread pool\n");
        goto end;
    }
    tpool.pool = pool;      //to share the pool for file read and write as well
    //queue to use with thread pool, for task and results
    if (!(queue = hts_tpool_process_init(pool, cnt * 2, 0))) {
        fprintf(stderr, "Failed to create queue\n");
        goto end;
    }
    //open input file - r reading
    if (!(infile = sam_open(inname, "r"))) {
        fprintf(stderr, "Could not open %s\n", inname);
        goto end;
    }
    //open output files - w write as SAM, wb  write as BAM
    if (!(outfile = sam_open(file, "wb"))) {
        fprintf(stderr, "Could not open output file\n");
        goto end;
    }
    //share the thread pool with i/o files
    if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
          hts_set_opt(outfile, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        fprintf(stderr, "Failed to set threads to i/o files\n");
        goto end;
    }
    //read header, required to resolve the target names to proper ids
    if (!(in_samhdr = sam_hdr_read(infile))) {
        fprintf(stderr, "Failed to read header from file!\n");
        goto end;
    }
    //write header
    if ((sam_hdr_write(outfile, in_samhdr) == -1)) {
        fprintf(stderr, "Failed to write header\n");
        goto end;
    }

    /* tasks are queued, worker threads get them and process in parallel;
    the results are queued and they are to be removed in parallel as well */

    // start output writer thread for ordered processing
    twritedata.outfile = outfile;
    twritedata.samhdr  = in_samhdr;
    twritedata.result  = 0;
    twritedata.queue   = queue;
    twritedata.cache   = &bamcache;
    if (pthread_create(&thread, NULL, threadfn_orderedwrite, &twritedata)) {
        fprintf(stderr, "Failed to create writer thread\n");
        goto end;
    }
    started_thread = 1;

    c = 0;
    while (c >= 0) {
        if (!(bamdata = getbamstorage(chunk, &bamcache))) {
            fprintf(stderr, "Failed to allocate memory\n");
            break;
        }
        //read alignments, upto max size for this lot
        for (cnt = 0; cnt < bamdata->maxsize; ++cnt) {
            c = sam_read1(infile, in_samhdr, bamdata->bamarray[cnt]);
            if (c < 0) {
                break;      // EOF or failure
            }
        }
        if (c >= -1 ) {
            //max size data or reached EOF
            bamdata->count = cnt;
            // Queue the data for processing.  hts_tpool_dispatch3() is
            // used here as it allows in-flight data to be cleaned up
            // properly when stopping early due to errors.
            if (hts_tpool_dispatch3(pool, queue, thread_ordered_proc, bamdata,
                                    cleanup_bamstorage, cleanup_bamstorage,
                                    0) == -1) {
                fprintf(stderr, "Failed to schedule processing\n");
                goto end;
            }
            bamdata = NULL;
        } else {
            fprintf(stderr, "Error in reading data\n");
            break;
        }
    }

    ret = EXIT_SUCCESS;

 end:
    // Tidy up after having dispatched all of the data.

    // Note that the order here is important.  In particular, we need
    // to join the thread that was started earlier before freeing anything
    // to avoid any use-after-free errors.

    // It's also possible to get here early due to various error conditions,
    // so we need to carefully check which parts of the program state have
    // been created before trying to clean them up.

    if (queue) {
        if (-1 == c) {
            // EOF read, send a marker to tell the threadfn_orderedwrite()
            // function to shut down.
            if (hts_tpool_dispatch(pool, queue, thread_ordered_proc,
                                   NULL) == -1) {
                fprintf(stderr, "Failed to schedule sentinel job\n");
                ret = EXIT_FAILURE;
            }
        } else {
            // Error or we never wrote anything.  Shut down the queue to
            // ensure threadfn_orderedwrite() wakes up and terminates.
            hts_tpool_process_shutdown(queue);
        }
    }

    // Wait for threadfn_orderedwrite to finish.
    if (started_thread) {
        pthread_join(thread, NULL);

        // Once the writer thread has finished, check the result it sent back
        if (twritedata.result != 0) {
            ret = EXIT_FAILURE;
        }
    }

    if (queue) {
        // Once threadfn_orderedwrite has stopped, the queue can be
        // cleaned up.
        hts_tpool_process_destroy(queue);
    }

    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        if (sam_close(infile) != 0) {
            ret = EXIT_FAILURE;
        }
    }
    if (outfile) {
        if (sam_close(outfile) != 0) {
            ret = EXIT_FAILURE;
        }
    }

    pthread_mutex_lock(&bamcache.lock);
    if (bamcache.list) {
        struct data *tmp = NULL;
        while (bamcache.list) {
            tmp = bamcache.list;
            bamcache.list = bamcache.list->next;
            cleanup_bamstorage(tmp);
        }
    }
    pthread_mutex_unlock(&bamcache.lock);

    if (file) {
        free(file);
    }
    if (pool) {
        hts_tpool_destroy(pool);
    }
    return ret;
}
