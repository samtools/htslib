/*  mpileup.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: mpileup infile ...\n\
Shows the mpileup api usage.\n");
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
    return 0;
}

int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
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
    plpconf** conf = NULL;
    bam_mplp_t mplpiter = NULL;
    int tid = -1, input = 0, k = 0, dpt = 0, *depth = NULL;
    hts_pos_t refpos = -1;
    const bam_pileup1_t **plp = NULL;

    //infile ...
    if (argc < 2) {
        print_usage(stderr);
        goto end;
    }
    if ((conf = calloc(argc - 1, sizeof(plpconf*)))) {
        for (input = 0; input < argc - 1; ++input) {
            conf[input] = calloc(1, sizeof(plpconf));
        }
    }
    depth = calloc(argc - 1, sizeof(int));
    plp = calloc(argc - 1, sizeof(bam_pileup1_t*));
    if (!conf || !depth || !plp) {
        printf("Failed to allocate memory\n");
        goto end;
    }
    for (input = 0; input < argc - 1; ++input) {
        conf[input]->inname = argv[input+1];
    }

    //initialize
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    //open input files
    for(input = 0; input < argc - 1; ++input) {
        if (!(conf[input]->infile = sam_open(conf[input]->inname, "r"))) {
            printf("Could not open %s\n", conf[input]->inname);
            goto end;
        }
        //read header
        if (!(conf[input]->in_samhdr = sam_hdr_read(conf[input]->infile))) {
            printf("Failed to read header from file!\n");
            goto end;
        }
    }

    if (!(mplpiter = bam_mplp_init(argc - 1, readdata, (void**) conf))) {
        printf("Failed to initialize mpileup data\n");
        goto end;
    }

    //set constructor destructor callbacks
    bam_mplp_constructor(mplpiter, plpconstructor);
    bam_mplp_destructor(mplpiter, plpdestructor);

    while (bam_mplp64_auto(mplpiter, &tid, &refpos, depth, plp) > 0) {
        printf("%d\t%"PRIhts_pos"\t", tid+1, refpos+1);

        for (input = 0; input < argc - 1; ++input) {
            for (dpt = 0; dpt  < depth[input]; ++dpt) {
                if (plp[input][dpt].is_del || plp[input][dpt].is_refskip) {
                    printf("*");
                    continue;
                }
                //start and end are displayed in UPPER and rest on LOWER
                printf("%c", plp[input][dpt].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b), plp[input][dpt].qpos)]) :
                                (plp[input]->is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b), plp[input][dpt].qpos)]) : tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b), plp[input][dpt].qpos)])));
                if (plp[input][dpt].indel > 0) {
                    //insertions, anyway not start or end
                    printf("+%d", plp[input][dpt].indel);
                    for (k = 0; k < plp[input][dpt].indel; ++k) {
                        printf("%c", tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b), plp[input][dpt].qpos + k + 1)]));
                    }
                }
                else if (plp[input][dpt].indel < 0) {
                    printf("%d", plp[input][dpt].indel);
                    for (k = 0; k < -plp[input][dpt].indel; ++k) {
                        printf("?");
                    }
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
    if (conf) {
        for (input = 0; input < argc - 1; ++input) {
            if (conf[input] && conf[input]->in_samhdr) {
                sam_hdr_destroy(conf[input]->in_samhdr);
            }
            if (conf[input] && conf[input]->infile) {
                sam_close(conf[input]->infile);
            }
            if (conf[input]) {
                free(conf[input]);
            }
        }
        free(conf);
    }

    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (mplpiter) {
        bam_mplp_destroy(mplpiter);
    }
    if (depth) {
        free(depth);
    }
    if (plp) {
        free(plp);
    }
    return ret;
}
