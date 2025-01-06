/*  split2.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: split2 infile outdir\n\
Splits the input file alignments to read1 and read2 and saves as 1.sam and 2.bam in given directory\n\
Shows file type selection through name and format api\n");
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
    char *file1 = NULL, *file2 = NULL, mode1[5] = "w", mode2[5] = "w";
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
    size = sizeof(char) * (strlen(outdir) + sizeof("/1.sam.gz") + 1); //space for output file name and null termination
    file1 = malloc(size);
    file2 = malloc(size);
    if (!file1 || !file2) {
        printf("Failed to set output path\n");
        goto end;
    }

    //output file names
    snprintf(file1, size, "%s/1.sam.gz", outdir);   //name of Read1 file
    snprintf(file2, size, "%s/2.sam", outdir);      //name of Read2 file
    //bam data storage
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    //set file open mode based on file name for 1st and as explicit for 2nd
    if ((sam_open_mode(mode1+1, file1, NULL) == -1) || (sam_open_mode(mode2+1, file2, "sam.gz") == -1)) {
        printf("Failed to set open mode\n");
        goto end;
    }
    //open input file
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    //open output files
    outfile1 = sam_open(file1, mode1);                          //as compressed SAM through sam_open
    outfile2 = sam_open_format(file2, mode2, NULL);             //as compressed SAM through sam_open_format
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
