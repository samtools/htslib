/*  index_write.c --  showcases the htslib api usage

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
#include <libgen.h>
#include <unistd.h>
#include <htslib/sam.h>

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: idx_on_write infile shiftsize outdir\n\
Creates compressed sam file and index file for it in given directory\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *outdir = NULL;
    char *inname = NULL, *fileidx = NULL, *outname = NULL, outmode[4] = "w";
    int c = 0, ret = EXIT_FAILURE, size = 0;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;

    //idx_on_write infile sizeshift outputdirectory
    if (argc != 4) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    size = atoi(argv[2]);
    outdir = argv[3];

    //allocate space for output name - outdir/filename.ext.idxextNUL
    c = strlen(basename(inname)) + strlen(outdir) + 10;
    fileidx = malloc(sizeof(char) * c);
    outname = malloc(sizeof(char) * c);
    if (!fileidx || !outname) {
        printf("Couldnt allocate memory\n");
        goto end;
    }
    //initialize bam storage
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }

    //open files
    if ((infile = sam_open(inname, "r"))) {
        //get file type and create output names
        if (infile->format.format == cram) {
            //set as crai
            snprintf(fileidx, c, "%s/%s.crai", outdir, basename(inname));
            snprintf(outname, c, "%s/%s", outdir, basename(inname));
        }
        else {
            //set as either bai or csi based on interval
            if (infile->format.format == sam && infile->format.compression == no_compression) {
                //create as gzip compressed
                snprintf(outname, c, "%s/%s.gz", outdir, basename(inname));
                snprintf(fileidx, c, "%s/%s.gz.%s", outdir, basename(inname), !size ? "bai" : "csi");
            }
            else {
                //with same name as input
                snprintf(outname, c, "%s/%s", outdir, basename(inname));
                snprintf(fileidx, c, "%s/%s.%s", outdir, basename(inname), !size ? "bai" : "csi");
            }
        }
    }
    c = 0;
    sam_open_mode(outmode + 1, outname, NULL);          //set extra write options based on name
    outfile = sam_open(outname, outmode);
    if (!outfile || !infile) {
        printf("Could not open files\n");
        goto end;
    }

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    //write header
    if (sam_hdr_write(outfile, in_samhdr)) {
        printf("Failed to write header\n");
        goto end;
    }

    // initialize indexing, before start of write
    if (sam_idx_init(outfile, in_samhdr, size, fileidx)) {
        printf("idx initialization failed\n");
        goto end;
    }
    //read and write alignments
    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
            printf("Failed to write data\n");
            goto end;
        }
    }
    if (c != -1) {
        printf("Error in reading data\n");
        goto end;
    }
    //else EOF, save index
    if (sam_idx_save(outfile)) {
        printf("Could not save index\n");
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
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (fileidx) {
        free(fileidx);
    }
    if (outname) {
        free(outname);
    }
    if (outfile) {
        sam_close(outfile);
    }
    return ret;
}
