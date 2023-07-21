/*  read_header.c --  showcases the htslib api usage

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

/// print_usage - print the demo_usage
/** @param fp pointer to the file / terminal to which demo_usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: read_header infile header [id val] [tag]\n\
This shows given tag from given header or the whole line\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *header = NULL, *tag = NULL, *idval = NULL;
    char *id = NULL;
    int c = 0, ret = EXIT_FAILURE, linecnt = 0;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    kstring_t data = KS_INITIALIZE;

    //read_header infile header tag
    if (argc < 3 || argc > 6) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    header = argv[2];
    if (argc == 4) {            //header and tag
        tag = argv[3];
        //find unique identifier field name for requested header type
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
            id = "";
        }
        else {
            printf("Invalid header type\n");
            goto end;
        }
    }
    else if (argc == 5) {       //header id val
        id = argv[3];
        idval = argv[4];
    }
    else if (argc == 6) {       //header id val tag
        id = argv[3];
        idval = argv[4];
        tag = argv[5];
    }

    //open input files
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    if (id && idval) {
        if (tag) {
            ret = sam_hdr_find_tag_id(in_samhdr, header, id, idval, tag, &data);
        }
        else {
            ret = sam_hdr_find_line_id(in_samhdr, header, id, idval, &data);
        }

        if (ret == 0) {
            printf("%s\n", data.s);
        }
        else  if (ret == -1) {
            printf("No matching tag found\n");
            goto end;
        }
        else {
            printf("Failed to find header line\n");
            goto end;
        }
    }
    else {
        //get count of given header type
        linecnt = sam_hdr_count_lines(in_samhdr, header);
        if (linecnt == 0) {
            printf("No matching line found\n");
            goto end;
        }
        for (c = 0; c < linecnt; ++c ) {
            if (tag) {
                //non CO, get the tag requested
                ret = sam_hdr_find_tag_pos(in_samhdr, header, c, tag, &data);
            }
            else {
                //CO header, there are no tags but the whole line
                ret = sam_hdr_find_line_pos(in_samhdr, header, c, &data);
            }

            if (ret == 0) {
                printf("%s\n", data.s);
                continue;
            }
            else if (ret == -1) {
                printf("Tag not present\n");
                continue;
            }
            else {
                printf("Failed to get tag\n");
                goto end;
            }
        }
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
    ks_free(&data);
    return ret;
}
