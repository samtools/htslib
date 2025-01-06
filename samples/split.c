/*  split.c --  showcases the htslib api usage

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

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: split infile outdir\n\
Splits the input file alignments to read1 and read2 and saves as 1.sam and 2.bam in given directory\n\
Shows the basic writing of output\n");
    return;
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
    int c = 0, ret = EXIT_FAILURE, size = 0;
    samFile *infile = NULL, *outfile1 = NULL, *outfile2 = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;

    if (argc != 3) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    outdir = argv[2];

    //allocate space for output
    size = sizeof(char) * (strlen(outdir) + sizeof("/1.sam") + 1); //space for output file name and null termination
    file1 = malloc(size);
    file2 = malloc(size);
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
    //open input file - r reading
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    //open output files - w write as SAM, wb  write as BAM
    outfile1 = sam_open(file1, "w");            //as SAM
    outfile2 = sam_open(file2, "wb");           //as BAM
    if (!outfile1 || !outfile2) {
        printf("Could not open output file\n");
        goto end;
    }

    //read header, required to resolve the target names to proper ids
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    //write header
    if ((sam_hdr_write(outfile1, in_samhdr) == -1) || (sam_hdr_write(outfile2, in_samhdr) == -1)) {
        printf("Failed to write header\n");
        goto end;
    }

    //check flags and write
    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (bamdata->core.flag & BAM_FREAD1) {
            if (sam_write1(outfile1, in_samhdr, bamdata) < 0) {
                printf("Failed to write output data\n");
                goto end;
            }
        }
        else if (bamdata->core.flag & BAM_FREAD2) {
            if (sam_write1(outfile2, in_samhdr, bamdata) < 0) {
                printf("Failed to write output data\n");
                goto end;
            }
        }
    }
    if (-1 == c) {
        //EOF
        ret = EXIT_SUCCESS;
    }
    else {
        printf("Error in reading data\n");
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
    return ret;
}
