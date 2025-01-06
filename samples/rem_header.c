/*  rem_header.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: rem_header infile header [id]\n\
Removes header line of given type and id\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *header = NULL, *idval = NULL;
    char *id = NULL;
    int ret = EXIT_FAILURE;
    samFile *infile = NULL, *outfile = NULL;
    sam_hdr_t *in_samhdr = NULL;

    //update_header infile header idval tag value
    if (argc <3 || argc > 4) {
        //3 & 4 are ok, 3-> all of given header type, 4->given id of given header type to be removed
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    header = argv[2];
    if (argc == 4) {
        idval = argv[3];
    }

    //unique identifier for each of the header types
    if (header[0] == 'H' && header[1] == 'D') {
        id = NULL;
    }
    else if (header[0] == 'S' && header[1] == 'Q') {
        id = "SN";
    }
    else if (header[0] == 'R' && header[1] == 'G') {
        id = "ID";
    }
    else if (header[0] == 'P' && header[1] == 'G') {
        id = "ID";
    }
    else if (header[0] == 'C' && header[1] == 'O') {
        //CO field can be removed using the position of it using sam_hdr_remove_line_pos
        id = "";
    }
    else {
        printf("Invalid header type\n");
        goto end;
    }

    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    if (!(outfile = sam_open("-", "w"))) {      //use stdout as the output file for ease of display of update
        printf("Could not open stdout\n");
        goto end;
    }

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }
    if (idval) {
        //remove specific line
        if (sam_hdr_remove_line_id(in_samhdr, header, id, idval)) {
            printf("Failed to remove header line\n");
            goto end;
        }
    }
    else {
        //remove multiple lines of a header type
        if (sam_hdr_remove_lines(in_samhdr, header, id, NULL)) {
            printf("Failed to remove header line\n");
            goto end;
        }
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
    return ret;
}
