/*  read_refname.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: read_refname infile minsize\n\
This shows name of references which has length above the given size\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *id = NULL;
    int c = 0, ret = EXIT_FAILURE, linecnt = 0, pos = 0;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    kstring_t data = KS_INITIALIZE;
    int64_t minsize = 0, size = 0;

    if (argc != 3 && argc != 2) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    if (argc == 3) {
        minsize = atoll(argv[2]);
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

    linecnt = sam_hdr_count_lines(in_samhdr, "SQ"); //get reference count
    if (linecnt <= 0) {
        if (!linecnt) {
            printf("No reference line present\n");
        }
        else {
            printf("Failed to get reference line count\n");
        }
        goto end;
    }
    //iterate and check each reference's length
    for (pos = 1, c = 0; c < linecnt; ++c) {
        if ((ret = sam_hdr_find_tag_pos(in_samhdr, "SQ", c, "LN", &data) == -2)) {
            printf("Failed to get length\n");
            goto end;
        }
        else if (ret == -1) {
            //length not present, ignore
            continue;
        }
        //else have length
        size = atoll(data.s);
        if (size < minsize) {
            //not required
            continue;
        }
        if (!(id = sam_hdr_line_name(in_samhdr, "SQ", c))) {    //sam_hdr_find_tag_pos(in_samhdr, "SQ", c, "SN", &data) can also do the same!
            printf("Failed to get id for reference data\n");
            goto end;
        }
        printf("%d,%s,%s\n", pos, id, data.s);
        pos++;
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
