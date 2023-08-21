/*  mod_aux.c --  showcases the htslib api usage

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
#include <strings.h>
#include <unistd.h>
#include <htslib/sam.h>

/// print_usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: mod_aux infile QNAME tag type val\n\
Add/update the given aux tag to all alignments\n\
type A-char C-int F-float Z-string\n");
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *tag = NULL, *qname = NULL, *val = NULL;
    char type = '\0';
    int ret = EXIT_FAILURE, ret_r = 0, length = 0;
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL, *outfile = NULL;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;

    //mod_aux infile QNAME tag type val
    if (argc != 6) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    qname = argv[2];
    tag = argv[3];
    type = argv[4][0];
    val = argv[5];

    if (!(bamdata = bam_init1())) {
        printf("Failed to allocate data memory!\n");
        goto end;
    }

    //open input file
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    //open output file
    if (!(outfile = sam_open("-", "w"))) {
        printf("Could not open std output\n");
        goto end;
    }

    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    if (sam_hdr_write(outfile, in_samhdr) == -1) {
        printf("Failed to write header\n");
        goto end;
    }

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (strcasecmp(bam_get_qname(bamdata), qname)) {
            if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                printf("Failed to write output\n");
                goto end;
            }
            continue;   //not matching
        }

        errno = 0;
        //matched to qname, update aux
        if (!(data = bam_aux_get(bamdata, tag))) {
            int i = 0; float f = 0;
            //tag not present append
            switch (type) {
                case 'f':
                case 'd':
                    length = sizeof(float);
                    f = atof(val);
                    val = (const char*) &f;
                    type = 'f';
                break;
                case 'C':
                case 'S':
                case 'I':
                    length = sizeof(int);
                    i = atoi(val);
                    val = (const char*) &i;
                break;
                case 'Z':
                    length = strlen(val) + 1;   //1 for NUL termination
                break;
                case 'A':
                    length = 1;
                break;
                default:
                    printf("Invalid type mentioned\n");
                    goto end;
                break;
            }
            if (bam_aux_append(bamdata, tag, type, length, (const uint8_t*)val)) {
                printf("Failed to append aux data, errno: %d\n", errno);
                goto end;
            }
        }
        else {
            char auxtype = bam_aux_type(data);
            //update the tag with newer value
            switch (type) {
                case 'f':
                case 'd':
                    if (auxtype != 'f' && auxtype != 'd') {
                        printf("Invalid aux type passed\n");
                        goto end;
                    }
                    if (bam_aux_update_float(bamdata, tag, atof(val))) {
                        printf("Failed to update float data, errno: %d\n", errno);
                        goto end;
                    }
                break;
                case 'C':
                case 'S':
                case 'I':
                    if (auxtype != 'c' && auxtype != 'C' && auxtype != 's' && auxtype != 'S' && auxtype != 'i' && auxtype != 'I') {
                        printf("Invalid aux type passed\n");
                        goto end;
                    }
                    if (bam_aux_update_int(bamdata, tag, atoll(val))) {
                        printf("Failed to update int data, errno: %d\n", errno);
                        goto end;
                    }
                break;
                case 'Z':
                    if (auxtype != 'Z') {
                        printf("Invalid aux type passed\n");
                        goto end;
                    }
                    length = strlen(val) + 1;   //1 for NUL termination
                    if (bam_aux_update_str(bamdata, tag, length, val)) {
                        //with length as -1, length will be detected based on null terminated val data
                        printf("Failed to update string data, errno: %d\n", errno);
                        goto end;
                    }
                break;
                case 'A':
                    if (auxtype != 'A') {
                        printf("Invalid aux type passed\n");
                        goto end;
                    }
                    //update the char data directly on buffer
                    *(data+1) = val[0];
                break;
                default:
                    printf("Invalid data type\n");
                    goto end;
                break;
            }
        }
        if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
            printf("Failed to write output\n");
            goto end;
        }
    }
    if (ret_r < -1) {
        //read error
        printf("Failed to read data\n");
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
    return ret;
}
