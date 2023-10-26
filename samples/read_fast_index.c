/*  read_fast_index.c --  showcases the htslib api usage

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

/* The purpose of this code is to demonstrate the library apis and need proper error handling and optimisation */

#include <getopt.h>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: read_fast_i <infile> A/Q 0/1 regiondef\n\
Reads the fasta/fastq file using index and shows the content.\n\
For fasta files use A and Q for fastq files.\n\
Region can be 1 or more of <reference name>[:start-end] entries separated by comma.\n\
For single region, give regcount as 0 and non 0 for multi-regions.\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *region = NULL, *data = NULL, *remaining = NULL;
    int ret = EXIT_FAILURE, tid = -1, usemulti = 0;
    faidx_t *idx = NULL;
    enum fai_format_options fmt = FAI_FASTA;
    hts_pos_t len = 0, beg = 0, end = 0;

    //read_fast_i infile A/Q regcount region
    if (argc != 5) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    if (argv[2][0] == 'Q') {
        fmt = FAI_FASTQ;
    }
    usemulti = atoi(argv[3]);
    region = argv[4];

    //load index
    if (!(idx = fai_load3_format(inname, NULL, NULL, FAI_CREATE, fmt))) {
        printf("Failed to load index\n");
        goto end;
    }

    if (!usemulti) {
        //get data from given region
        if (!(data = fai_fetch64(idx, region, &len))) {
            if (-1 == len) {
                printf("Failed to get data\n");                 //failure
                goto end;
            }
            else {
                printf("Data not found for given region\n");    //no data
            }
        }
        else {
            printf("Data: %"PRId64" %s\n", len, data);
            free((void*)data);
            //get quality for fastq type
            if (fmt == FAI_FASTQ) {
                if (!(data = fai_fetchqual64(idx, region, &len))) {
                    if (len == -1) {
                        printf("Failed to get data\n");
                        goto end;
                    }
                    else {
                        printf("Data not found for given region\n");
                    }
                }
                else {
                    printf("Qual: %"PRId64" %s\n", len, data);
                    free((void*)data);
                }
            }
        }
    }
    else {
        //parse, get each region and get data for each
        while ((remaining = fai_parse_region(idx, region, &tid, &beg, &end, HTS_PARSE_LIST))) {     //here expects regions as csv
            //parsed the region, correct end points based on actual data
            if (fai_adjust_region(idx, tid, &beg, &end) == -1) {
                printf("Error in adjusting region for tid %d\n", tid);
                goto end;
            }
            //get data for given region
            if (!(data = faidx_fetch_seq64(idx, faidx_iseq(idx, tid), beg, end, &len))) {
                if (len == -1) {
                    printf("Failed to get data\n");                 //failure
                    goto end;
                }
                else {
                    printf("No data found for given region\n");     //no data
                }
            }
            else {
                printf("Data: %"PRIhts_pos" %s\n", len, data);
                free((void*)data);
                data = NULL;

                //get quality data for fastq
                if (fmt == FAI_FASTQ) {
                    if (!(data = faidx_fetch_qual64(idx, faidx_iseq(idx, tid), beg, end, &len))) {
                        if (len == -1) {
                            printf("Failed to get qual data\n");
                            goto end;
                        }
                        else {
                            printf("No data found for given region\n");
                        }
                    }
                    else {
                        printf("Qual: %"PRIhts_pos" %s\n", len, data);
                        free((void*)data);
                        data = NULL;
                    }
                }
            }
            region = remaining;                                     //parse remaining region defs
        }
    }

    ret = EXIT_SUCCESS;
end:
    //clean up
    if (idx) {
        fai_destroy(idx);
    }
    return ret;
}
