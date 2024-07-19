/*  write_fast.c --  showcases the htslib api usage

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
#include <time.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: write_fast <file> <sequence> [<qualities]\n\
Appends a fasta/fastq file.\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *outname = NULL;             //output file name
    int ret = EXIT_FAILURE;
    samFile *outfile = NULL;                //sam file
    bam1_t *bamdata = NULL;                 //to hold the read data
    char mode[4] = "a";
    const char *data = NULL, *qual = NULL;  //ref data and quality
    char name[256] = {0};

    if (argc > 4 || argc < 3) {
        print_usage(stdout);
        goto end;
    }
    outname = argv[1];
    data = argv[2];
    if (argc == 4) {    //fastq data
        qual = argv[3];
        if (strlen(data) != strlen(qual)) {     //check for proper length of data and quality values
            printf("Incorrect reference and quality data\n");
            goto end;
        }
    }

    //initialize
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    if (sam_open_mode(mode + 1, outname, NULL) < 0) {
        printf("Invalid file name\n");
        goto end;
    }
    //open output file
    if (!(outfile = sam_open(outname, mode))) {         //expects the name to have correct extension!
        printf("Could not open %s\n", outname);
        goto end;
    }
    /* if the file name extension is not appropriate to the content, inconsistent data will be present in output.
    if required, htsFormat and sam_open_format can be explicitly used to ensure appropriateness of content.
    htsFormat fmt = {sequence_data, fastq_format / fasta_format};
    sam_open_format(outname, mode, fmt);
    */

    snprintf(name, sizeof(name), "Test_%ld", (long) time(NULL));
    //data
    if (bam_set1(bamdata, strlen(name), name, BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, strlen(data), data, qual, 0) < 0) {
        printf("Failed to set data\n");
        goto end;
    }
    //as we write only FASTA/FASTQ, we can get away without providing headers
    if (sam_write1(outfile, NULL, bamdata) < 0) {
        printf("Failed to write data\n");
        goto end;
    }
    ret = EXIT_SUCCESS;
end:
    //clean up
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    return ret;
}
