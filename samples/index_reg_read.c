/*  index_reg_read.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: read_reg infile idxfile region\n\
Reads alignments matching to a specific region\n\
\\. from start of file\n\
\\* only unmapped reads\n\
REFNAME all reads referring REFNAME\n\
REFNAME:S all reads referring REFNAME and overlapping from S onwards\n\
REFNAME:S-E all reads referring REFNAME overlapping from S to E\n\
REFNAME:-E all reads referring REFNAME overlapping upto E\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *region = NULL;
    char *idxfile = NULL;
    int c = 0, ret = EXIT_FAILURE;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;

    //readreg infile indexfile region
    if (argc != 4) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    idxfile = argv[2];
    region = argv[3];

    //initialize bam data storage
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }

    //open files
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open input file\n");
        goto end;
    }
    //using stdout as output file for ease of dumping data
    if (!(outfile = sam_open("-", "w"))) {
        printf("Could not open out file\n");
        goto end;
    }
    //load index file
    if (!(idx = sam_index_load2(infile, inname, idxfile))) {
        printf("Failed to load the index\n");
        goto end;
    }
    //can use sam_index_load if the index file is present in same location and follows standard naming conventions (i.e. <filename.ext>.<index type ext crai/bai/csi/tbi>)

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    //create iterator
    if (!(iter = sam_itr_querys(idx, in_samhdr, region))) {
        printf("Failed to get iterator\n");
        goto end;
    }
    //read using iterator
    while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
        //write to output
        if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
            printf("Failed to write output\n");
            goto end;
        }
    }
    if (c != -1) {
        printf("Error during read\n");
        goto end;
    }
    ret = EXIT_SUCCESS;

end:
    //cleanup
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (iter) {
        sam_itr_destroy(iter);
    }
    if (idx) {
        hts_idx_destroy(idx);
    }
    return ret;
}
