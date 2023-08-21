/*  mod_bam.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: mod_bam infile QNAME fieldpos newval\n\
Modifies the alignment data field\n\
fieldpos - 1 QNAME 2 FLAG 3 RNAME 4 POS 5 MAPQ 6 CIGAR 7 RNEXT 8 PNEXT 9 TLEN 10 SEQ 11 QUAL\n");
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *qname = NULL;
    char *val = NULL;
    int c = 0, ret = EXIT_FAILURE, field = 0;
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL, *outfile = NULL;
    int ret_r = 0, i = 0;
    bam1_t *bamdata = NULL;

    //mod_bam infile QNAME fieldpos newval
    if (argc != 5) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    qname = argv[2];
    //1 QNAME 2 FLAG 3 RNAME 4 POS 5 MAPQ 6 CIGAR 7 RNEXT 8 PNEXT 9 TLEN 10 SEQ 11 QUAL
    field = atoi(argv[3]);
    val = argv[4];

    if (!(bamdata = bam_init1())) {
        printf("Failed to allocate data memory!\n");
        goto end;
    }

    //open input file
    if (!(infile = sam_open(inname, "r")) || !(outfile = sam_open("-", "w"))) {
        printf("Could not open input/output\n");
        goto end;
    }
    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    if (sam_hdr_write(outfile, in_samhdr) == -1) {
        printf("Failed to write header\n");
        goto end;
    }

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
        //QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG:TYPE:VALUE]…
        ret = 0;
        if (!strcasecmp(qname, bam_get_qname(bamdata))) {
            //the required one
            switch(field) {
            case 1:// QNAME
                ret = bam_set_qname(bamdata, val);
            break;
            case 2:// FLAG
                bamdata->core.flag = atol(val) & 0xFFFF;
            break;
            case 3:// RNAME
            case 7:// RNEXT
                if ((ret = sam_hdr_name2tid(in_samhdr, val)) < 0) {
                    printf("Invalid reference name\n");
                    ret = -1;
                    break;
                }
                if (field == 3) {
                    //reference
                    bamdata->core.tid = ret;
                }
                else {
                    //mate reference
                    bamdata->core.mtid = ret;
                }
            break;
            case 4:// POS
                bamdata->core.pos = atoll(val);
            break;
            case 5:// MAPQ
                bamdata->core.qual = atoi(val) & 0x0FF;
            break;
            case 6:// CIGAR
            {
                uint32_t *cigar = NULL;
                size_t size = 0;
                ssize_t ncigar = 0;
                bam1_t *newbam = bam_init1();
                if (!newbam) {
                    printf("Failed to create new bam data\n");
                    ret = -1;
                    break;
                }
                //get cigar array and set all data in new bam record
                if ((ncigar = sam_parse_cigar(val, NULL, &cigar, &size)) < 0) {
                    printf("Failed to parse cigar\n");
                    ret = -1;
                    break;
                }
                if (bam_set1(newbam, bamdata->core.l_qname, bam_get_qname(bamdata), bamdata->core.flag, bamdata->core.tid, bamdata->core.pos, bamdata->core.qual,
                    ncigar, cigar, bamdata->core.mtid, bamdata->core.mpos, bamdata->core.isize, bamdata->core.l_qseq, (const char*)bam_get_seq(bamdata), (const char*)bam_get_qual(bamdata), bam_get_l_aux(bamdata)) < 0) {
                    printf("Failed to set bamdata\n");
                    ret = -1;
                    break;
                }
                //correct sequence data as input is expected in ascii format and not as compressed inside bam!
                memcpy(bam_get_seq(newbam), bam_get_seq(bamdata), (bamdata->core.l_qseq + 1) / 2);
                //copy the aux data
                memcpy(bam_get_aux(newbam), bam_get_aux(bamdata), bam_get_l_aux(bamdata));

                bam_destroy1(bamdata);
                bamdata = newbam;
            }
            break;
            case 8:// PNEXT
                bamdata->core.mpos = atoll(val);
            break;
            case 9:// TLEN
                bamdata->core.isize = atoll(val);
            break;
            case 10:// SEQ
                i = strlen(val);
                if (bamdata->core.l_qseq != i) {
                    printf("SEQ length different\n");
                    ret = -1;
                    //as it is different, have to update quality data and cigar data as well and more info is required for it, which is not handled in this sample
                    //accessing raw memory and moving is one option; creating and using new bam1_t object is another option.
                    break;
                }
                for( c = 0; c < i; ++c) {
                    bam_set_seqi(bam_get_seq(bamdata), c, seq_nt16_table[(unsigned char)val[c]]);
                }
            break;
            case 11:// QUAL
                i = strlen(val);
                if (i != bamdata->core.l_qseq) {
                    printf("Qual length different than sequence\n");
                    ret = -1;
                    break;
                }
                for (c = 0; c < i; ++c) {
                    val[c] -= 33;               //phred score from ascii value
                }
                memcpy(bam_get_qual(bamdata), val, i);
            break;
            default:
                printf("Invalid input\n");
                goto end;
            break;
            }
            if (ret < 0) {
                printf("Failed to set new data\n");
                ret = EXIT_FAILURE;
                goto end;
            }
        }
        if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
            printf("Failed to write bam data\n");
            ret = EXIT_FAILURE;
            goto end;
        }
    }

    if (ret_r == -1 || ret != EXIT_FAILURE) {
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
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    return ret;
}
