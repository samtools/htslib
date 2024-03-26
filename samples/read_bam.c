/*  read_bam.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: read_bam infile\n\
Shows the alignment data from file\n");
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *tidname = NULL, *flags = NULL;
    int ret = EXIT_FAILURE;
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL;

    int ret_r = 0, i = 0;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;
    uint32_t *cigar = NULL;


    //read_bam infile
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
    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
        //QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG:TYPE:VALUE]…
        printf("NAME: %s\n", bam_get_qname(bamdata));                                   //get the query name using the macro
        flags = bam_flag2str(bamdata->core.flag);                                       //flags as string
        printf("FLG: %d - %s\n", bamdata->core.flag, flags);                            //flag is available in core structure
        free((void*)flags);
        tidname = sam_hdr_tid2name(in_samhdr, bamdata->core.tid);
        printf("RNAME/TID: %d - %s\n", bamdata->core.tid, tidname? tidname: "" );       //retrieves the target name using the value in bam and by referring the header
        printf("POS: %"PRIhts_pos"\n", bamdata->core.pos + 1);                          //internally position is 0 based and on text output / SAM it is 1 based
        printf("MQUAL: %d\n", bamdata->core.qual);                                      //map quality value

        cigar = bam_get_cigar(bamdata);                                                 //retrieves the cigar data
        printf("CGR: ");
        for (i = 0; i < bamdata->core.n_cigar; ++i) {                                   //no. of cigar data entries
            printf("%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));       //the macros gives the count of operation and the symbol of operation for given cigar entry
        }
        printf("\nTLEN/ISIZE: %"PRIhts_pos"\n", bamdata->core.isize);

        data = bam_get_seq(bamdata);                                                    //get the sequence data
        if (bamdata->core.l_qseq != bam_cigar2qlen(bamdata->core.n_cigar, cigar)) {     //checks the length with CIGAR and query
            printf("\nLength doesnt matches to cigar data\n");
            goto end;
        }

        printf("SEQ: ");
        for (i = 0; i < bamdata->core.l_qseq ; ++i) {                                   //sequence length
            printf("%c", seq_nt16_str[bam_seqi(data, i)]);                              //retrieves the base from (internal compressed) sequence data
        }
        printf("\nQUAL: ");
        for (int i = 0; i < bamdata->core.l_qseq ; ++i) {
            printf("%c", bam_get_qual(bamdata)[i]+33);                                  //retrives the quality value
        }
        printf("\n\n");
    }

    if (ret_r == -1) {
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
    return ret;
}
