/*  modstate.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: modstate infile option\n\
Shows the base modifications on the alignment\n\
Option can be 1 or 2 to select the api to use\n");
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL;
    int ret = EXIT_FAILURE;
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL;

    int ret_r = 0, i = 0 , r = 0, j = 0, pos = 0, opt = 0, k = 0, cnt = 0, *bm = NULL;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;
    hts_base_mod_state *ms = NULL;


    //modstate infile 1/2
    if (argc != 3) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    opt = atoi(argv[2]) - 1;    //option 1 or 2?

    if (!(bamdata = bam_init1())) {
        printf("Failed to allocate data memory!\n");
        goto end;
    }

    if (!(ms = hts_base_mod_state_alloc())) {
        printf("Failed to allocate state memory\n");
        goto end;
    }

    //open input file
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
        i = 0;
        data = bam_get_seq(bamdata);
        if (bam_parse_basemod(bamdata, ms)) {
            printf("Failed to parse the base mods\n");
            goto end;
        }
        //dump the modifications
        printf("Modifications:");
        bm = bam_mods_recorded(ms, &cnt);
        for (k = 0; k < cnt; ++k) {
            printf("%c", bm[k]);
        }
        printf("\n");
        hts_base_mod mod[5] = {0};  //for ATCGN
        if (opt) {
            //option 1
            for (; i < bamdata->core.l_qseq; ++i) {
                if ((r = bam_mods_at_next_pos(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]))) <= -1) {
                    printf("Failed to get modifications\n");
                    goto end;
                }
                else if (r > (sizeof(mod) / sizeof(mod[0]))) {
                    printf("More modifications than this app can handle, update the app\n");
                    goto end;
                }
                else if (!r) {
                    //no modification at this pos
                    printf("%c", seq_nt16_str[bam_seqi(data, i)]);
                }
                //modifications
                for (j = 0; j < r; ++j) {
                    printf("%c%c%c", mod[j].canonical_base, mod[j].strand ? '-' : '+', mod[j].modified_base);
                }
            }
        }
        else {
            //option 2
            while ((r = bam_next_basemod(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]), &pos)) >= 0) {
                for (; i < bamdata->core.l_qseq && i < pos; ++i) {
                    printf("%c", seq_nt16_str[bam_seqi(data, i)]);
                }
                //modifications
                for (j = 0; j < r; ++j) {
                    printf("%c%c%c", mod[j].canonical_base, mod[j].strand ? '-' : '+', mod[j].modified_base);
                }
                if (i == pos)
                    i++;        //skip the modification already displayed
                if (!r) {
                    for (; i < bamdata->core.l_qseq; ++i) {
                        printf("%c", seq_nt16_str[bam_seqi(data, i)]);
                    }
                    break;
                }
            }
            if (r <= -1) {
                printf("Failed to get modifications\n");
                goto end;
            }
        }
        printf("\n");
    }

    if (ret_r == -1) {
        //check last alignment's base modification
        int strand = 0, impl = 0;
        char canonical = 0, modification[] = "mhfcgebaon";      //possible modifications
        printf("\n\nLast alignment has \n");
        for (k = 0; k < sizeof(modification) - 1; ++k) {        //avoiding NUL termination
            if (bam_mods_query_type(ms, modification[k], &strand, &impl, &canonical)) {
                printf ("No modification of %c type\n", modification[k]);
            }
            else {
                printf("%s strand has %c modified with %c, can %sassume unlisted as unmodified\n", strand?"-/bottom/reverse":"+/top/forward", canonical, modification[k], impl?"" : "not " );
            }
        }
        // no error!
        ret = EXIT_SUCCESS;
    }
    else {
        printf("Failed to read data\n");
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

    if (ms) {
        hts_base_mod_state_free(ms);
    }
    return ret;
}
