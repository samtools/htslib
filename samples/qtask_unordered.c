/*  qtask_unordered.c --  showcases the htslib api usage

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
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

typedef struct beddata {
    char *name;                 //chromosome name
    hts_pos_t start;            //start position
    hts_pos_t end;              //end position
} beddata;

typedef struct splitdata {
    beddata *region;            //region information
    const char *infile;         //input file
    const char *outdir;         //output path
} splitdata;

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: qtask_unordered infile threadcount outdir bedfile\n\
Splits the input to files specific for regions given in bedfile. Output is\n\
saved in outdir. Expects the index file to be present along with infile.\n");
    return;
}

/// readbedfile - read the bedfile and return region array
/** @param file bed file name
 *  @param data - output, pointer to array of beddata data, to hold bed regions
returns number of regions in data
*/
int readbedfile(const char *file, struct beddata **data)
{
    int ret = -1, cnt = 0 , max = 0;
    kstring_t line = KS_INITIALIZE;
    htsFile *bedfile = NULL;
    char *sptr = NULL, *token = NULL;
    const char *sep ="\t";
    beddata *bedtoks = NULL;

    if (!data || *data) {
        printf("Invalid argument\n");
        goto fail;
    }
    if (!(bedfile = hts_open(file, "r"))) {
        printf("Failed to open bedfile\n");
        goto fail;
    }
    //get lines one by one and get region details
    while (!(ret = kgetline(&line, (kgets_func*)hgets, bedfile->fp.hfile))) {
        if (!line.l) {
            continue;               //skip empty line
        }
        if  (line.s[0] == '#' || !strncmp("track ", line.s, sizeof("track ") - 1)) {
            ks_clear(&line);
            continue;               //ignore track lines and comments
        }
        token = strtok_r(line.s, sep, &sptr);
        if (token) {                //allocate memory for regions
            if ((cnt+1) > max) {
                max += 100;         //for another 100 regions
                bedtoks = realloc(bedtoks, sizeof(beddata) * max);
            }
        }
        else {
            break;
        }
        bedtoks[cnt].name = strdup(token);  //chromosome name
        token = strtok_r(NULL, sep, &sptr);
        if (!token) {
            break;
        }                                   //start position
        bedtoks[cnt].start = token ? atoll(token) : 0;
        token = strtok_r(NULL, sep, &sptr);
        if (!token) {
            break;
        }                                   //end position
        bedtoks[cnt++].end = token ? atoll(token) : 0;
        ks_clear(&line);
    }
    if (ret != EOF) {
        goto fail;
    }
    ret = cnt;
    if (bedfile) {
        hts_close(bedfile);
    }
    ks_free(&line);
    if (!cnt) {
        goto fail;
    }
    *data = bedtoks;
    return ret;
fail:
    if (bedfile) {
        hts_close(bedfile);
    }
    ks_free(&line);
    if (bedtoks) {
        for (max = cnt, cnt = 0; cnt < max; ++cnt) {
            free(bedtoks[cnt].name);
        }
        free(bedtoks);
    }
    bedtoks = NULL;
    return 0;
}

/// splittoregions - saves the relevant data to separate file
/** @param args pointer to set of data to be processed
returns NULL
the processing could be in any order based on the number of threads in use
*/
void * splittoregions(void *args)
{
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    hts_itr_t *iter = NULL;
    hts_idx_t *idx = NULL;
    splitdata *data = (splitdata*)args;
    char *file = NULL, *region = NULL;
    int size = 0, ret = 0;
    if (!(infile = sam_open(data->infile, "r"))) {
        printf("Failed to open input file\n");
        goto end;
    }
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header data\n");
        goto end;
    }
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    size = strlen(data->region->name) + 50;         //region specification
    if (!(region = malloc(size))) {
        printf("Failed to allocate memory\n");
        goto end;
    }
    snprintf(region, size, "%s:%"PRIhts_pos"-%"PRIhts_pos, data->region->name, data->region->start, data->region->end);
    size += strlen(data->outdir);                   //output file with path
    if (!(file = malloc(size))) {
        printf("Failed to allocate memory\n");
        goto end;
    }
    snprintf(file, size, "%s/%s_%"PRIhts_pos"_%"PRIhts_pos".sam", data->outdir, data->region->name, data->region->start, data->region->end);
    if (!(idx = sam_index_load(infile, data->infile))) {
        printf("Failed to load index\n");
        goto end;
    }
    if (!(iter = sam_itr_querys(idx, in_samhdr, region))) {
        printf("Failed to create iterator\n");
        goto end;
    }
    if (!(outfile = sam_open(file, "w"))) {
        printf("Failed to open output file\n");
        goto end;
    }
    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        printf("Failed to write header\n");
        goto end;
    }
    while ((ret = sam_itr_next(infile, iter, bamdata)) >= 0) {  //read and write relevant data
        if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
            printf("Failed to write data\n");
            goto end;
        }
    }
    if (ret != -1) {
        printf("Failed to get all data\n");
    }

end:
    free(data);
    if (infile) {
        sam_close(infile);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (iter) {
        sam_itr_destroy(iter);
    }
    if (idx) {
        hts_idx_destroy(idx);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (file) {
        free(file);
    }
    if (region) {
        free(region);
    }
    return NULL;
}

/// main - splits the data to region specific files
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *outdir = NULL, *bedfile = NULL;
    int c = 0, ret = EXIT_FAILURE, cnt = 0, regcnt = 0;
    hts_tpool *pool = NULL;
    hts_tpool_process *queue = NULL;
    beddata *regions = NULL;

    //qtask infile threadcount outdir [chunksize]
    if (argc != 5) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    cnt = atoi(argv[2]);
    outdir = argv[3];
    bedfile = argv[4];
    //get regions from bedfile
    if ((regcnt = readbedfile(bedfile, &regions)) <= 0) {
        printf("Failed to get bed data\n");
        goto end;
    }
    if (cnt < 1) {
        cnt = 1;
    }
    if (!(pool = hts_tpool_init(cnt))) {                //thread pool
        printf("Failed to create thread pool\n");
        goto end;
    }
    //queue to use with thread pool, for tasks
    if (!(queue = hts_tpool_process_init(pool, cnt * 2, 1))) {
        printf("Failed to create queue\n");
        goto end;
    }
    for (c = 0; c < regcnt; ++c) {
        struct splitdata *task = malloc(sizeof(splitdata));
        task->infile = inname;
        task->outdir = outdir;
        task->region = regions + c;
        //schedule jobs to run in parallel
        if (hts_tpool_dispatch(pool, queue, splittoregions, task) < 0) {
            printf("Failed to schedule processing\n");
            goto end;
        }
    }
    //trigger processing for anything pending, NOTE: will be blocked until queue is cleared
    if (hts_tpool_process_flush(queue) == -1) {
        printf("Failed to flush queues\n");
        goto end;
    }
    ret = EXIT_SUCCESS;

    //shutdown queues to exit the result wait
    hts_tpool_process_shutdown(queue);

end:
    //cleanup
    for (c = 0; c < regcnt; ++c) {
        free(regions[c].name);
    }
    free(regions);

    if (queue) {
        hts_tpool_process_destroy(queue);
    }
    if (pool) {
        hts_tpool_destroy(pool);
    }
    return ret;
}
