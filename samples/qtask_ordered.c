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

typedef struct orderedwrite {
    samFile *outfile;           //output file handle
    sam_hdr_t *samhdr;          //header used to write data
    hts_tpool_process *queue;   //queue from which results to be retrieved

    pthread_cond_t *done;       //indicates exit condition
    pthread_mutex_t *lock;      //to synchronise queue access
} orderedwrite;

typedef struct data {
    int count;                  //used up size
    int size;                   //max size
    bam1_t **bamarray;          //bam1_t array for optimal queueing
} data;

#define WAIT 1                  //1 sec
/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: qtask_ordered infile threadcount outdir [chunksize]\n\
Calculates GC ratio - sum(G,C) / sum(A,T,C,G,N) - and adds to each alignment\n\
as xr:f aux tag. Output is saved in outdir.\n\
chunksize [100] sets the number of alignments clubbed together to process.\n");
    return;
}

/// addcount - calculates and adds the aux tag
/** @param bamdata pointer to bam data
returns 0 on success and -1 for failure
*/
int addcount(bam1_t *bamdata)
{
    int pos = 0, gc = 0, ret = 0;
    float gcratio = 0;
    uint8_t *data = bam_get_seq(bamdata);
    for (pos = 0; pos < bamdata->core.l_qseq; ++pos) {
        switch(seq_nt16_str[bam_seqi(data, pos)]) {
            case 'G':   //fall through
            case 'C':
                gc++;
            break;
        }
    }
    gcratio = gc / (float) bamdata->core.l_qseq;

    if (bam_aux_append(bamdata, "xr", 'f', sizeof(gcratio), (const uint8_t*)&gcratio) < 0) {
        printf("Failed to add aux tag xr, errno: %d\n", errno);
        ret = -1;
    }
    return ret;
}


/// getbamstorage - allocates storage for alignments to queue
/** @param chunk no of alignments to be queued together
returns allocated data.
*/
data* getbamstorage(int chunk)
{
    int i = 0;
    data *bamdata = malloc(sizeof(data));
    if (!bamdata) {
        return NULL;
    }
    bamdata->bamarray = malloc(chunk * sizeof(bam1_t*));
    if (!bamdata->bamarray) {
        return NULL;
    }
    for (i = 0; i < chunk; ++i) {
        bamdata->bamarray[i] = bam_init1();
    }
    bamdata->count = 0;
    bamdata->size = chunk;

    return bamdata;
}

/// thread_ordered_proc - does the processing of task in queue and queues the output back
/** @param args pointer to set of data to be processed
returns the processed data
the processing could be in any order based on the number of threads in use but read of output
from queue will be in order
*/
void *thread_ordered_proc(void *args)
{
    int i = 0;
    data *bamdata = (data*)args;
    for ( i = 0; i < bamdata->count; ++i) {
        //add count
        if (addcount(bamdata->bamarray[i]) < 0) {
            printf("Failed to calculate gc data\n");
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

    struct timeval now;
    struct timespec timeout;
    long usec = 0;

    pthread_mutex_lock(tdata->lock);                //lock to check the exit condition
    do {
        //get result and write; wait if no result is in queue - until shutdown of queue
        while ((r = hts_tpool_next_result(tdata->queue))) {
            bamdata = (data*) hts_tpool_result_data(r);
            for (i = 0; i < bamdata->count; ++i) {
                if (sam_write1(tdata->outfile, tdata->samhdr, bamdata->bamarray[i]) < 0) {
                    printf("Failed to write output data\n");
                    break;
                }
            }
            hts_tpool_delete_result(r, 0);          //release the result memory
            if (bamdata) {
                for (i = 0; i < bamdata->size; ++i) {
                    bam_destroy1(bamdata->bamarray[i]); //clear the bamdata;
                }
                free(bamdata->bamarray);
                free(bamdata);
            }
        }
        //no more data in queues, check and wait till exit is triggered
        gettimeofday(&now, NULL);
        usec = now.tv_usec + 100000;    //+100msec
        if (usec >= 1000000) {          //overflow
            usec %= 1000000;
            now.tv_sec++;
        }
        timeout.tv_sec = now.tv_sec;
        timeout.tv_nsec = usec * 1000;

    } while (pthread_cond_timedwait(tdata->done, tdata->lock, &timeout) == ETIMEDOUT);

    pthread_mutex_unlock(tdata->lock);

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
    int c = 0, ret = EXIT_FAILURE, cnt = 0, clearthread = 0, chunk = 0;
    size_t size = 0;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    pthread_mutex_t stopcondsynch = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t stopcond = PTHREAD_COND_INITIALIZER;
    pthread_t thread;
    orderedwrite twritedata = {0};
    hts_tpool *pool = NULL;
    hts_tpool_process *queue = NULL;
    htsThreadPool tpool = {NULL, 0};
    data *bamdata = NULL;

    //qtask infile threadcount outdir [chunksize]
    if (argc != 4 && argc != 5) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    cnt = atoi(argv[2]);
    outdir = argv[3];
    if (argc == 5) {
        chunk = atoi(argv[4]);
    }

    if (cnt < 1) {
        cnt = 1;
    }
    if (chunk < 1) {
        chunk = 10000;
    }

    //allocate space for output
    size = (strlen(outdir) + sizeof("/out.sam") + 1); //space for output file name and null termination
    file = malloc(size);
    if (!file) {
        printf("Failed to set output path\n");
        goto end;
    }
    snprintf(file, size, "%s/out.sam", outdir);         //output file name
    if (!(pool = hts_tpool_init(cnt))) {                //thread pool
        printf("Failed to create thread pool\n");
        goto end;
    }
    tpool.pool = pool;      //to share the pool for file read and write as well
    //queue to use with thread pool, for task and results
    if (!(queue = hts_tpool_process_init(pool, cnt * 2, 0))) {
        printf("Failed to create queue\n");
        goto end;
    }
    //open input file - r reading
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    //open output files - w write as SAM, wb  write as BAM
    if (!(outfile = sam_open(file, "w"))) {
        printf("Could not open output file\n");
        goto end;
    }
    //share the thread pool with i/o files
    if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
          hts_set_opt(outfile, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        printf("Failed to set threads to i/o files\n");
        goto end;
    }
    //read header, required to resolve the target names to proper ids
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    /*tasks are queued, worker threads get them and processes in parallel;
    the results are queued and they are to be removed in parallel as well*/
    //start output writer thread for ordered processing
    twritedata.outfile = outfile;
    twritedata.queue = queue;
    twritedata.done = &stopcond;
    twritedata.lock = &stopcondsynch;
    twritedata.samhdr = in_samhdr;
    if (pthread_create(&thread, NULL, threadfn_orderedwrite, &twritedata)) {
        printf("Failed to create writer thread\n");
        goto end;
    }
    clearthread = 1;
    //write header
    if ((sam_hdr_write(outfile, in_samhdr) == -1)) {
        printf("Failed to write header\n");
        goto end;
    }

    c = 0;
    while (c >= 0) {
        bamdata = getbamstorage(chunk);
        //read alignments, upto max size for this lot
        for (cnt = 0; cnt < bamdata->size; ++cnt) {
            c = sam_read1(infile, in_samhdr, bamdata->bamarray[cnt]);
            if (c >= 0) {
                continue;   //read next
            }
            else {
                break;      //failure
            }
        }
        if (c >= -1 ) {
            //max size data or reached EOF
            bamdata->count = ( c >= 0 )? bamdata->size : cnt;
            //queue the lot for processing
            if (hts_tpool_dispatch(pool, queue, thread_ordered_proc,
                bamdata) == -1) {
                printf("Failed to schedule processing\n");
                goto end;
            }
            bamdata = NULL;
        }
        else {
            printf("Error in reading data\n");
            break;
        }
    }
    if (-1 == c) {
        //EOF read, trigger processing for anything pending, NOTE: will be blocked until queue is cleared
        if (hts_tpool_process_flush(queue) == -1) {
            printf("Failed to flush queues\n");
            goto end;
        }
        //all tasks done, check for processing completion
        while (1) {
            if (!hts_tpool_process_empty(queue)) {
                usleep(WAIT * 1000000);     //results yet to be empty, check again
                continue;
            }
            break;
        }
        ret = EXIT_SUCCESS;
    }

    //trigger exit for ordered write thread
    pthread_mutex_lock(twritedata.lock);
    pthread_cond_signal(twritedata.done);
    pthread_mutex_unlock(twritedata.lock);
end:
    //cleanup
    if (clearthread) {
        pthread_join(thread, NULL);
    }
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (bamdata) {
        for (cnt = 0; cnt < bamdata->size; ++cnt) {
            bam_destroy1(bamdata->bamarray[cnt]);
        }
        free(bamdata);
    }
    if (file) {
        free(file);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (queue) {
        hts_tpool_process_destroy(queue);
    }
    if (pool) {
        hts_tpool_destroy(pool);
    }
    return ret;
}
