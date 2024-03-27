/*  mod_aux_ba.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: mod_aux_ba infile\n\
Updates the count of bases as an aux array on all alignments\n\
BA:B:I,count of ACTGN\n");
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL;
    int i = 0, ret = EXIT_FAILURE, ret_r = 0;
    uint32_t cnt[5] = {0};  //A C G T N
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL, *outfile = NULL;
    bam1_t *bamdata = NULL;

    //mod_aux infile
    if (argc != 2) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];

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
        errno = 0;
        memset(cnt, 0, sizeof(cnt));
        for (i = 0; i < bamdata->core.l_qseq; ++i) {
            switch (seq_nt16_str[bam_seqi(bam_get_seq(bamdata),i)]) {
                case 'A':
                ++cnt[0];
                break;
                case 'C':
                ++cnt[1];
                break;
                case 'G':
                ++cnt[2];
                break;
                case 'T':
                ++cnt[3];
                break;
                default:    //N
                ++cnt[4];
                break;
            }
        }

        if (bam_aux_update_array(bamdata, "BA", 'I', sizeof(cnt)/sizeof(cnt[0]), cnt)) {
            printf("Failed to update base array, errno %d", errno);
            goto end;
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
