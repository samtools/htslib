/*  pileup_mod.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: pileup_mod infile\n\
Shows the pileup api usage with base modification.\n");
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
    //plpconf *conf= (plpconf*)data; can use this to access anything required from the data in pileup init

    //when using cd, initialize and use as it will be reused after destructor
    cd->p = hts_base_mod_state_alloc();
    if (!cd->p) {
        printf("Failed to allocate base modification state\n");
        return 1;
    }

    //parse the bam data and gather modification data from MM tags
    return (-1 == bam_parse_basemod(b, (hts_base_mod_state*)cd->p)) ? 1 : 0;
}

int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    if (cd->p) {
        hts_base_mod_state_free((hts_base_mod_state *)cd->p);
        cd->p = NULL;
    }
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
    int tid = -1, depth = -1, j = 0, k = 0, inslen = 0, dellen = 0, modlen = 0;
    #define NMODS 5
    hts_base_mod mods[NMODS] = {0};     //ACGT N
    int refpos = -1;
    const bam_pileup1_t *plp = NULL;
    kstring_t insdata = KS_INITIALIZE;

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

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &depth))) {
        memset(&mods, 0, sizeof(mods));
        printf("%d\t%d\t", tid+1, refpos+1);

        for (j = 0; j < depth; ++j) {
            dellen = 0;

            if (plp[j].is_del || plp[j].is_refskip) {
                printf("*");
                continue;
            }
            /*invoke bam_mods_at_qpos before bam_plp_insertion_mod that the base modification
            is retrieved before change in pileup pos thr' plp_insertion_mod call*/
            if ((modlen = bam_mods_at_qpos(plp[j].b, plp[j].qpos, plp[j].cd.p, mods, NMODS)) == -1) {
                printf("Failed to get modifications\n");
                goto end;
            }

            //use plp_insertion/_mod to get insertion and del at the same position
            if ((inslen = bam_plp_insertion_mod(&plp[j], (hts_base_mod_state*)plp[j].cd.p, &insdata, &dellen)) == -1) {
                printf("Failed to get insertion status\n");
                goto end;
            }

            //start and end are displayed in UPPER and rest on LOWER, only 1st modification considered
            //base and modification
            printf("%c%c%c", plp[j].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                (plp[j].is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                    tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)])),
                modlen > 0 ? mods[0].strand ? '-' : '+' : '\0',
                modlen > 0 ? mods[0].modified_base : '\0');
            //insertion and deletions
            if (plp[j].indel > 0) {
                //insertion
                /*insertion data from plp_insertion_mod, note this shows the quality value as well
                which is different from base and modification above;the lower case display is not attempted either*/
                printf("+%d%s", plp[j].indel, insdata.s);
                //handle deletion if any
                if (dellen) {
                    printf("-%d", dellen);
                    for (k = 0; k < dellen; ++k) {
                        printf("?");
                    }
                }
            }
            else if (plp[j].indel < 0) {
                //deletion
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
    ks_free(&insdata);
    return ret;
}
