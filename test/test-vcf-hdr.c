/*  test/test-vcf-api.c -- VCF test harness.

    Copyright (C) 2015 EMBL - European Bioinformatics Institute

    Author: Cristina Yenyxe Gonzalez Garcia <cyenyxe@ebi.ac.uk>

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
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>


int main(int argc, char **argv)
{
    char *in_fname = argc>1 ? argv[1] : "test/test-vcf-hdr-in.vcf";
    
    // Parse header with trailing spaces
    htsFile *file = hts_open(in_fname, "r");
    bcf_hdr_t *header = bcf_hdr_read(file);
        
    if (header->nhrec != 15)
    {
        fprintf(stderr, "The header contains %d lines, 15 were expected\n", header->nhrec);
        hts_close(file);
        exit(1);
    }
    
    // Write results
    char *out_fname = "test/test-vcf-hdr.out.new";
    htsFile *out = hts_open(out_fname, "wu");
    vcf_hdr_write(out, header);
    
    int ret;
    if (ret = hts_close(out))
    {
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", out_fname, ret);
        exit(ret);
    }
    
    if (ret = hts_close(file))
    {
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", in_fname, ret);
        exit(ret);
    }
    
    return 0;
}

