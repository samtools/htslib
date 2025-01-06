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

struct datacache;

typedef struct basecount {
    uint64_t counts[16];        //count of all bases
} basecount;

typedef struct data {
    int count;                  //used up size
    int maxsize;                //max size per data chunk
    bam1_t **bamarray;          //bam1_t array for optimal queueing

    struct datacache *cache;
    basecount *bases;           //count of all possible bases
    struct data *next;          //pointer to next one - to reuse earlier allocations
} data;

typedef struct datacache
{
    pthread_mutex_t lock;       //synchronizes the access to cache
    data *list;                 //data storage
} datacache;

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: qtask_unordered infile threadcount [chunksize]\n\
Shows the base counts and calculates GC ratio - sum(G,C) / sum(A,T,C,G)\n\
chunksize [4096] sets the number of alignments clubbed together to process.\n");
    return;
}

/// getbamstorage - allocates storage for alignments to queue
/** @param chunk number of bam data to allocate
 * @param bases storage of result
 * @param bamcache cached storage
returns already allocated data storage if one is available, otherwise allocates new
*/
data* getbamstorage(int chunk, basecount *bases, datacache *bamcache)
{
    int i = 0;
    data *bamdata = NULL;

    if (!bamcache || !bases) {
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

        bamdata->bases = bases;
        bamdata->cache = bamcache;
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

    bamdata->bases = bases;
    bamdata->cache = bamcache;

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

/// thread_unordered_proc - does the processing of task in queue and updates result
/** @param args pointer to set of data to be processed
returns NULL
the processing could be in any order based on the number of threads in use
*/
void *thread_unordered_proc(void *args)
{
    int i = 0;
    data *bamdata = (data*)args;
    uint64_t pos = 0;
    uint8_t *data = NULL;
    uint64_t counts[16] = {0};
    for ( i = 0; i < bamdata->count; ++i) {
        data = bam_get_seq(bamdata->bamarray[i]);
        for (pos = 0; pos < bamdata->bamarray[i]->core.l_qseq; ++pos) {
            /* it is faster to count all bases and select required ones later
            compared to select and count here */
            counts[bam_seqi(data, pos)]++;
        }
    }
    //update result and add the memory block for reuse
    pthread_mutex_lock(&bamdata->cache->lock);
    for (i = 0; i < 16; i++) {
        bamdata->bases->counts[i] += counts[i];
    }

    bamdata->next = bamdata->cache->list;
    bamdata->cache->list = bamdata;
    pthread_mutex_unlock(&bamdata->cache->lock);

    return NULL;
}

/// main - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL;
    int c = 0, ret = EXIT_FAILURE, cnt = 0, chunk = 0;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    hts_tpool *pool = NULL;
    hts_tpool_process *queue = NULL;
    htsThreadPool tpool = {NULL, 0};
    data *bamdata = NULL;
    basecount gccount = {{0}};
    datacache bamcache = {PTHREAD_MUTEX_INITIALIZER, NULL};

    //qtask infile threadcount [chunksize]
    if (argc != 3 && argc != 4) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    cnt = atoi(argv[2]);
    if (argc == 4) {
        chunk = atoi(argv[3]);
    }
    if (cnt < 1) {
        cnt = 1;
    }
    if (chunk < 1) {
        chunk = 4096;
    }

    if (!(pool = hts_tpool_init(cnt))) {
        fprintf(stderr, "Failed to create thread pool\n");
        goto end;
    }
    tpool.pool = pool;      //to share the pool for file read and write as well
    //queue to use with thread pool, for tasks
    if (!(queue = hts_tpool_process_init(pool, cnt * 2, 1))) {
        fprintf(stderr, "Failed to create queue\n");
        goto end;
    }
    //open input file - r reading
    if (!(infile = sam_open(inname, "r"))) {
        fprintf(stderr, "Could not open %s\n", inname);
        goto end;
    }
    //share the thread pool with i/o files
    if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        fprintf(stderr, "Failed to set threads to i/o files\n");
        goto end;
    }
    //read header, required to resolve the target names to proper ids
    if (!(in_samhdr = sam_hdr_read(infile))) {
        fprintf(stderr, "Failed to read header from file!\n");
        goto end;
    }

    /*tasks are queued, worker threads get them and process in parallel;
    all bases are counted instead of counting atcg alone as it is faster*/

    c = 0;
    while (c >= 0) {
        //use cached storage to avoid allocate/deallocate overheads
        if (!(bamdata = getbamstorage(chunk, &gccount, &bamcache))) {
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
            if (hts_tpool_dispatch3(pool, queue, thread_unordered_proc, bamdata,
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

     if (-1 == c) {
        // EOF read, ensure all are processed, waits for all to finish
        if (hts_tpool_process_flush(queue) == -1) {
            fprintf(stderr, "Failed to flush queue\n");
        } else { //all done
            //refer seq_nt16_str to find position of required bases
            fprintf(stdout, "GCratio: %f\nBase counts:\n",
                (gccount.counts[2] /*C*/ + gccount.counts[4] /*G*/) / (float)
                    (gccount.counts[1] /*A*/ + gccount.counts[8] /*T*/ +
                        gccount.counts[2] + gccount.counts[4]));

            for (cnt = 0; cnt < 16; ++cnt) {
                fprintf(stdout, "%c: %"PRIu64"\n", seq_nt16_str[cnt], gccount.counts[cnt]);
            }

            ret = EXIT_SUCCESS;
        }
    }
 end:
    if (queue) {
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

    if (pool) {
        hts_tpool_destroy(pool);
    }
    return ret;
}
