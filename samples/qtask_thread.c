/*  split_thread3.c --  showcases the htslib api usage

    Copyright (C) 2023 Genome Research Ltd.

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

/* The pupose of this code is to demonstrate the library apis and need proper error handling and optimization */

#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

//thread specific data
typedef struct td_unordered {
    samFile *outfile;           //output file handle
    sam_hdr_t *samhdr;          //header used to write data
    bam1_t *bamdata;            //data to be written

    pthread_mutex_t *lock;      //to synchronize write to file
} td_unordered;

typedef struct td_orderedwrite {
    samFile *outfile;           //output file handle
    sam_hdr_t *samhdr;          //header used to write data
    hts_tpool_process *queue;   //queue from which results to be retrieved

    pthread_cond_t *done;       //indicates exit condition
    pthread_mutex_t *lock;      //to synchronize queue access
} td_orderedwrite;

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: qtask infile threadcount outdir\n\
Splits the input file alignments to read1 and read2 and saves as 1.sam and 2.bam in given directory\n\
GC ratio - sum(G,C) / sum(A,T,C,G,N) - is calculated and added on each alignment as xr:f aux tag.\n");
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


/// thread_unordered_proc - does the processing of task in queue1 (Read1) and writes output
/** @param args pointer to data specific for the thread
returns none
the output could be unordered based on the number of threads in use
*/
void *thread_unordered_proc(void *args)
{
    td_unordered *tdata = (td_unordered*)args;

    //add count
    if (addcount(tdata->bamdata) < 0) {
        printf("Failed to calculate gc data\n");
    }

    //lock to have exclusive write access to outfile
    pthread_mutex_lock(tdata->lock);
    if (sam_write1(tdata->outfile, tdata->samhdr, tdata->bamdata) < 0) {
        printf("Failed to write output data\n");
    }
    pthread_mutex_unlock(tdata->lock);

    return NULL;
}

/// thread_ordered_proc - does the processing of task in queue2 (Read2) and queues the output back
/** @param argss pointer to data specific for the thread
returns the processed bamdata
the processing could be unordered based on the number of threads in use but read of output from queue
will be in order
*/
void *thread_ordered_proc(void *args)
{
    bam1_t *bamdata = (bam1_t*)args;
    //add count
    if (addcount(bamdata) < 0) {
        printf("Failed to calculate gc data\n");
    }
    return bamdata;
}

/// thread_ordered_proc - thread that read the output from queue2 (Read2) and writes
/** @param argss pointer to data specific for the thread
returns NULL
*/
void *threadfn_orderedwrite(void *args)
{
    td_orderedwrite *tdata = (td_orderedwrite*)args;
    hts_tpool_result *r = NULL;
    bam1_t *bamdata = NULL;

    struct timeval now;
    struct timespec timeout;

    //get time to set the wait time on exit condition check
    gettimeofday(&now, NULL);
    timeout.tv_sec = now.tv_sec + 0;
    timeout.tv_nsec = 0;

    pthread_mutex_lock(tdata->lock);                //lock to check the exit condition
    while (pthread_cond_timedwait(tdata->done, tdata->lock, &timeout) == ETIMEDOUT) {
        //exit not set, get result and write; wait if no result is in queue - until shutdown of queue
        while ((r = hts_tpool_next_result_wait(tdata->queue))) {
            bamdata = (bam1_t*) hts_tpool_result_data(r);
            if (bamdata && sam_write1(tdata->outfile, tdata->samhdr, bamdata) < 0) {
                printf("Failed to write output data\n");
            }
            hts_tpool_delete_result(r, 0);          //release the result memory
            if (bamdata) {
                bam_destroy1(bamdata);              //clear the bamdata;
                bamdata = NULL;
            }
        }
        //no more data in queues, check and wait till exit is triggered
        gettimeofday(&now, NULL);
        timeout.tv_sec = now.tv_sec + 1;            //atmost 1 sec wait
        timeout.tv_nsec = 0;
    }
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
    char *file1 = NULL, *file2 = NULL;
    int c = 0, ret = EXIT_FAILURE, size = 0, cnt = 0;
    samFile *infile = NULL, *outfile1 = NULL, *outfile2 = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    pthread_mutex_t filesynch = PTHREAD_MUTEX_INITIALIZER, stopcondsynch = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t stopcond = PTHREAD_COND_INITIALIZER;
    pthread_t thread = 0;
    td_orderedwrite twritedata = {0};
    hts_tpool *pool = NULL;
    hts_tpool_process *queue1 = NULL, *queue2 = NULL;

    //split_t3 infile threadcount outdir
    if (argc != 4) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    cnt = atoi(argv[2]);
    outdir = argv[3];

    if (cnt < 1) {
        cnt = 1;
    }

    //allocate space for output
    size = sizeof(char) * (strlen(outdir) + sizeof("/1.sam") + 1); //space for output file name and null termination
    file1 = malloc(size); file2 = malloc(size);
    if (!file1 || !file2) {
        printf("Failed to set output path\n");
        goto end;
    }

    //output file names
    snprintf(file1, size, "%s/1.sam", outdir);  //for SAM output
    snprintf(file2, size, "%s/2.bam", outdir);  //for BAM output
    //bam data storage
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    //thread pool
    if (!(pool = hts_tpool_init(cnt))) {
        printf("Failed to create thread pool\n");
        goto end;
    }
    //queue to use with thread pool, no queuing of results for Q1 and results queued for Q2
    if (!(queue1 = hts_tpool_process_init(pool, cnt * 2, 1)) || !(queue2 = hts_tpool_process_init(pool, cnt * 2, 0))) {
        printf("Failed to create queue\n");
        goto end;
    }
    //open input file - r reading
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }

    //open output files - w write as SAM, wb  write as BAM
    outfile1 = sam_open(file1, "w"); outfile2 = sam_open(file2, "wb");
    if (!outfile1 || !outfile2) {
        printf("Could not open output file\n");
        goto end;
    }
    //read header, required to resolve the target names to proper ids
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    //start output writer thread for ordered processing
    twritedata.outfile = outfile2; twritedata.queue = queue2; twritedata.done = &stopcond; twritedata.lock = &stopcondsynch;
    twritedata.samhdr = in_samhdr;
    if (pthread_create(&thread, NULL, threadfn_orderedwrite, &twritedata)) {
        printf("Failed to create writer thread\n");
        goto end;
    }
    //write header
    if ((sam_hdr_write(outfile1, in_samhdr) == -1) || (sam_hdr_write(outfile2, in_samhdr) == -1)) {
        printf("Failed to write header\n");
        goto end;
    }
    //check flags and schedule
    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (bamdata->core.flag & BAM_FREAD1) {
            //read1, scheduling for unordered output, using Q1; processing occurs in order of thread execution
            td_unordered *tdata = calloc(1, sizeof(td_unordered));
            tdata->bamdata = bamdata; tdata->outfile = outfile1; tdata->lock = &filesynch; tdata->samhdr = in_samhdr;
            if (hts_tpool_dispatch3(pool, queue1, thread_unordered_proc, tdata, (void (*)(void *))&free, (void (*)(void*))&bam_destroy1, 0) == -1) {
                printf("Failed to schedule unordered processing\n");
                goto end;
            }
        }
        else if (bamdata->core.flag & BAM_FREAD2) {
            //read2, scheduling for ordered output, using Q2; proces and result queueing occurs in order of thread execution and result retrieval
            //showcases the threaded execution and ordered result handling on read2
            if (hts_tpool_dispatch3(pool, queue2, thread_ordered_proc, bamdata, NULL, (void (*)(void*))&bam_destroy1, 0) == -1) {
                printf("Failed to schedule ordered processing\n");
                goto end;
            }
        }
        bamdata = bam_init1();                  //for next read
    }
    if (-1 == c) {
        //EOF read
        //clear the queues of any tasks; NOTE: will be blocked until queue is cleared
        if (hts_tpool_process_flush(queue1) == -1 || hts_tpool_process_flush(queue2) == -1) {
            printf("Failed to flush queues\n");
            goto end;
        }
        ret = EXIT_SUCCESS;
    }
    else {
        printf("Error in reading data\n");
    }

    //trigger exit for ordered write thread
    if (thread) {
        //shutown queues to exit the result wait
        hts_tpool_process_shutdown(queue1);
        hts_tpool_process_shutdown(queue2);

        pthread_mutex_lock(twritedata.lock);
        pthread_cond_signal(twritedata.done);
        pthread_mutex_unlock(twritedata.lock);
        pthread_join(thread, NULL);
        thread = NULL;
    }
end:
    //cleanup
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (file1) {
        free(file1);
    }
    if (file2) {
        free(file2);
    }
    if (outfile1) {
        sam_close(outfile1);
    }
    if (outfile2) {
        sam_close(outfile2);
    }
    if (queue1) {
        hts_tpool_process_destroy(queue1);
    }
    if (queue2) {
        hts_tpool_process_destroy(queue2);
    }
    if (pool) {
        hts_tpool_destroy(pool);
    }
    if (thread) {
        pthread_join(thread, NULL);
    }
    return ret;
}
