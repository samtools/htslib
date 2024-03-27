/*  index_fasta.c --  showcases the htslib api usage

    Copyright (C) 2024 Genome Research Ltd.

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
#include <time.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: index_fasta <file>\n\
Indexes a fasta/fastq file and saves along with source.\n");
    return;
}

/// main - indexes fasta/fastq file
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *filename = NULL;             //file name
    int ret = EXIT_FAILURE;

    if (argc != 2) {
        print_usage(stdout);
        goto end;
    }
    filename = argv[1];

    // index the file
    if (fai_build3(filename, NULL, NULL) == -1) {
        printf("Indexing failed with %d\n", errno);
        goto end;
    }
    //this creates an .fai file. If the file is bgzipped, a .gzi file will be created along with .fai
    ret = EXIT_SUCCESS;
end:
    //clean up
    return ret;
}
