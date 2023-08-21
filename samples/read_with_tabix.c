/*  read_with_tabix.c --  showcases the htslib api usage

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
#include <htslib/tbx.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: read_with_tbx <infile> <shift> <region>\n\
Reads sam.gz file using tabix index.\n");
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
    int ret = EXIT_FAILURE, built = 0, shift = 0, c = 0;
    tbx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    samFile *infile = NULL;
    kstring_t data = KS_INITIALIZE;

    //read_with_tbx infile shift region
    if (argc != 4) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];
    shift = atoi(argv[2]);
    region = argv[3];

    //open file for read
    if (!(infile = sam_open(inname, "r"))) {
        printf("Failed to open file\n");
        goto end;
    }
load:
    //load index
    if (!(idx = tbx_index_load3(inname, NULL, HTS_IDX_SILENT_FAIL))) {
        if (!built) {                   //try to build index and then load
            built = 1;                  //reset flag to avoid repeated attempt
            if (tbx_index_build3(inname, NULL, shift, 1, &tbx_conf_sam) == -1) {
                printf("Failed to create index\n");
            }
            else {
                goto load;              //built, try again to load
            }
        /*tbx_conf_sam for compressed sam, tbx_conf_vcf for compressed vcf,
        tbx_conf_bed for bed ...*/
        }
        printf("Failed to load index\n");
        goto end;
    }

    //read using index and region
    if (!(iter = tbx_itr_querys(idx, region))) {
        printf("Failed to get iterator\n");
        goto end;
    }
    while ((c = tbx_itr_next(infile, idx, iter, &data)) >= 0) {
        printf("%s\n", data.s);
    }
    if (c != -1) {
        printf("Failed to read all data\n");
        goto end;
    }
    //tabix doesnt support multi region description but application can explicitly extract region
    //from such queries and use on tbx_itr_querys
    /*while (rem = hts_parse_region(..region.)) {
        if(rem[-1] == ',') rem[-1]='\0';
        iter = tbx_itr_querys(..region);
        ....
        region = rem;
    }*/

    ret = EXIT_SUCCESS;
end:
    //clean up
    ks_free(&data);

    if (idx) {
        tbx_destroy(idx);
    }
    if (infile) {
        sam_close(infile);
    }
    return ret;
}
