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

/* The pupose of this code is to demonstrate the library apis and need proper error handling and optimization */

#include <getopt.h>
#include <unistd.h>
#include <htslib/sam.h>

/// print_usage - show flags_demo usage
/** @param fp pointer to the file / terminal to which demo_usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: write_fast <file>\n\
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
    sam_hdr_t *out_samhdr = NULL;           //header of file
    bam1_t *bamdata = NULL;                 //to hold the read data
    char mode[4] = "a";

    if (argc != 2) {
        print_usage(stdout);
        goto end;
    }
    outname = argv[1];

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
    if (!(outfile = sam_open(outname, mode))) {
        printf("Could not open %s\n", outname);
        goto end;
    }
    //dummy data
    if (bam_set1(bamdata, sizeof("test"), "test", BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, 10, "AACTGACTGA", "1234567890", 0) < 0) {
        printf("Failed to set data\n");
        goto end;
    }
    if (sam_write1(outfile, out_samhdr, bamdata) < 0) {
        printf("Failed to write data\n");
        goto end;
    }

    ret = EXIT_SUCCESS;
end:
    //clean up
    if (out_samhdr) {
        sam_hdr_destroy(out_samhdr);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    return ret;
}
