/*  read_aux.c --  showcases the htslib api usage

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
    fprintf(fp, "Usage: read_aux infile tag\n\
Read the given aux tag from alignments either as SAM string or as raw data\n");
}

/// printauxdata - prints aux data
/** @param fp - file to which it to be printed - stdout or null
 *  @param type - aux type
 *  @param idx - index in array, -1 when not an array type
 *  @param data - data
 *  recurses when the data is array type
returns 1 on failure 0 on success
*/
int printauxdata(FILE *fp, char type, int32_t idx, const uint8_t *data)
{
    uint32_t auxBcnt = 0;
    int i = 0;
    char auxBType = 'Z';

    //the tag is already queried and ensured to exist and the type is retrieved from the tag data, also iterated within index for arrays, so no error is expected here.
    //when these apis are used explicitly, these error conditions needs to be handled based on return value and errno
    switch(type) {
    case 'A':
        fprintf(fp, "%c", bam_aux2A(data));                                                 //byte data
        break;
    case 'c':
        fprintf(fp, "%d", (int8_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));    //signed 1 byte data; bam_auxB2i - from array or bam_aux2i - non array data
        break;
    case 'C':
        fprintf(fp, "%u", (uint8_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));   //unsigned 1 byte data
        break;
    case 's':
        fprintf(fp, "%d", (int16_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));   //signed 2 byte data
        break;
    case 'S':
        fprintf(fp, "%u", (uint16_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));  //unsigned 2 byte data
        break;
    case 'i':
        fprintf(fp, "%d", (int32_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));   //signed 4 byte data
        break;
    case 'I':
        fprintf(fp, "%u", (uint32_t)(idx > -1 ? bam_auxB2i(data, idx) : bam_aux2i(data)));  //unsigned 4 byte data
        break;
    case 'f':
    case 'd':
        fprintf(fp, "%g", (float)(idx > -1 ? bam_auxB2f(data, idx) : bam_aux2f(data)));     //floating point data, 4 bytes
        break;
    case 'H':
    case 'Z':
        fprintf(fp, "%s", bam_aux2Z(data));                                                 //array of char or hex data
        break;
    case 'B':                                                                               //array of char/int/float
        auxBcnt = bam_auxB_len(data);                                                       //length of array
        auxBType = bam_aux_type(data + 1);                                                  //type of element in array
        fprintf(fp, "%c", auxBType);
        for (i = 0; i < auxBcnt; ++i) {                                                     //iterate the array
            fprintf(fp, ",");
            //calling recursively  with index to reuse a few lines
            if (printauxdata(fp, auxBType, i, data) == EXIT_FAILURE) {
                return EXIT_FAILURE;
            }
        }
        break;
    default:
        printf("Invalid aux tag?\n");
        return EXIT_FAILURE;
        break;
    }
    return EXIT_SUCCESS;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL, *tag = NULL;
    int c = 0, ret = EXIT_FAILURE, ret_r = 0, i = 0;
    sam_hdr_t *in_samhdr = NULL;
    samFile *infile = NULL;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;
    kstring_t sdata = KS_INITIALIZE;

    //read_aux infile tag
    if (argc != 3) {
        print_usage(stderr);
        goto end;
    }
    inname = argv[1];
    tag = argv[2];

    if (!(bamdata = bam_init1())) {
        printf("Failed to allocate data memory!\n");
        goto end;
    }

    //open input file
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }

    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf("Failed to read header from file!\n");
        goto end;
    }

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        errno = 0; i++;
        ks_clear(&sdata);
        if (i % 2) {    //use options alternatively to demonstrate both
            //option 1 - get data as string with tag and type
            if ((c = bam_aux_get_str(bamdata, tag, &sdata)) == 1) {
                printf("%s\n",sdata.s);
            }
            else if (c == 0 && errno == ENOENT) {
                //tag not present
                printf("Tag not present\n");
            }
            else {
                //error
                printf("Failed to get tag\n");
                goto end;
            }
        }
        else {
            //option 2 - get raw data
            if (!(data = bam_aux_get(bamdata, tag))) {
                //tag data not returned, errno gives the reason
                if (errno == ENOENT) {
                    printf("Tag not present\n");
                }
                else {
                    printf("Invalid aux data\n");
                }
            }
            else {
                //got the tag, read and print
                if (printauxdata(stdout, bam_aux_type(data), -1, data) == EXIT_FAILURE) {
                    printf("Failed to read aux data\n");
                    goto end;
                }
                printf("\n");
            }
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
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    ks_free(&sdata);
    return ret;
}
