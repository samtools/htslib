/*  add_header.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: add_header infile\n\
Adds new header lines of SQ, RG, PG and CO types\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, sq[] = "@SQ\tSN:TR1\tLN:100\n@SQ\tSN:TR2\tLN:50";
    int c = 0, ret = EXIT_FAILURE;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    kstring_t data = KS_INITIALIZE;

    //update_header infile header idval tag value
    if (argc != 2) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];

    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    if (!(outfile = sam_open("-", "w"))) {                  //use stdout as the output file for ease of display of update
        printf("Could not open stdout\n");
        goto end;
    }

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    //dump command line arguments for PG line
    for (c = 0;  c < argc; ++c) {
        kputs(argv[c], &data);
        kputc(' ', &data);
    }

    //add SQ line with SN as TR1 and TR2
    if (sam_hdr_add_lines(in_samhdr, &sq[0], 0)) {          //length as 0 for NULL terminated data
        printf("Failed to add SQ lines\n");
        goto end;
    }

    //add RG line with ID as RG1
    if (sam_hdr_add_line(in_samhdr, "RG", "ID", "RG1", "LB", "Test", "SM", "S1", NULL)) {
        printf("Failed to add RG line\n");
        goto end;
    }

    //add pg line
    if (sam_hdr_add_pg(in_samhdr, "add_header", "VN", "Test", "CL", data.s, NULL)) {            //NULL is to indicate end of args
        printf("Failed to add PG line\n");
        goto end;
    }

    if (sam_hdr_add_line(in_samhdr, "CO", "Test data", NULL)) {                                 //NULL is to indicate end of args
        printf("Failed to add PG line\n");
        goto end;
    }

    //write output
    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        printf("Failed to write output\n");
        goto end;
    }
    ret = EXIT_SUCCESS;
    //bam data write to follow....
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
    ks_free(&data);
    return ret;
}
