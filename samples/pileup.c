/*  pileup.c --  showcases the htslib api usage

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
#include <ctype.h>
#include <htslib/sam.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: pileup infile\n\
Shows the pileup api usage.\n");
    return;
}

typedef struct plpconf {
    char *inname;
    samFile *infile;
    sam_hdr_t *in_samhdr;
} plpconf;

/// @brief plpconstructor
/// @param data client data?
/// @param b bam being loaded
///     @param cd client data
/// @return
int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    /*plpconf *conf= (plpconf*)data;
    can access the data passed to pileup init from data
    can do any alignment specific allocation / data storage here in param cd
    it can hold either a float, 64 bit int or a pointer
    when using cd, initialize and use as it will be reused after destructor*/
    return 0;
}

int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    /*plpconf *conf= (plpconf*)data;
    can access the data passed to pileup init from data
    deallocate any alignment specific allocation made in constructor and stored in param cd*/
    return 0;
}

/// @brief bam_plp_auto_f reads alignment data for pileup operation
/// @param data client callback data holding alignment file handle
/// @param b bamdata read
/// @return same as sam_read1
int readdata(void *data, bam1_t *b)
{
    plpconf *conf = (plpconf*)data;
    if (!conf || !conf->infile) {
        return -2;  //cant read data
    }

    //read alignment and send
    return sam_read1(conf->infile, conf->infile->bam_header, b);
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    int ret = EXIT_FAILURE;
    bam1_t *bamdata = NULL;
    plpconf conf = {0};
    bam_plp_t plpiter = NULL;
    int tid = -1, n = -1, j = 0, k = 0;
    int refpos = -1;
    const bam_pileup1_t *plp = NULL;

    //infile
    if (argc != 2) {
        print_usage(stderr);
        goto end;
    }
    conf.inname = argv[1];

    //initialize
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    //open input files
    if (!(conf.infile = sam_open(conf.inname, "r"))) {
        printf("Could not open %s\n", conf.inname);
        goto end;
    }
    //read header
    if (!(conf.in_samhdr = sam_hdr_read(conf.infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    if (!(plpiter = bam_plp_init(readdata, &conf))) {
        printf("Failed to initialize pileup data\n");
        goto end;
    }

    //set constructor destructor callbacks
    bam_plp_constructor(plpiter, plpconstructor);
    bam_plp_destructor(plpiter, plpdestructor);

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &n))) {
        printf("%d\t%d\t", tid+1, refpos+1);

        for (j = 0; j < n; ++j) {
            //doesnt detect succeeding insertion and deletion together here, only insertion is identified
            //deletion is detected in plp->is_del as and when pos reaches the position
            //if detection ahead is required, use bam_plp_insertion here which gives deletion length along with insertion
            if (plp[j].is_del || plp[j].is_refskip) {
                printf("*");
                continue;
            }
            //start and end are displayed in UPPER and rest on LOWER
            printf("%c", plp[j].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                            (plp[j].is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) : tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)])));
            if (plp[j].indel > 0) {
                //insertions, anyway not start or end
                printf("+%d", plp[j].indel);
                for (k = 0; k < plp[j].indel; ++k) {
                    printf("%c", tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos + k + 1)]));
                }
            }
            else if (plp[j].indel < 0) {
                printf("%d", plp[j].indel);
                for (k = 0; k < -plp[j].indel; ++k) {
                    printf("?");
                }
            }
            printf(" ");
        }
        printf("\n");
        fflush(stdout);
    }

    ret = EXIT_SUCCESS;
end:
    //clean up
    if (conf.in_samhdr) {
        sam_hdr_destroy(conf.in_samhdr);
    }
    if (conf.infile) {
        sam_close(conf.infile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (plpiter) {
        bam_plp_destroy(plpiter);
    }
    return ret;
}
