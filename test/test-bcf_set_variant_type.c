/*  test/test-bcf_set_variant_type.c -- bcf_set_variant_type test harness.

    Copyright (C) 2022 Genome Research Ltd.

    Author: Martin Pollard <mp15@sanger.ac.uk>

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

#include <config.h>

#include <string.h>

#include "../htslib/hts.h"
#include "../vcf.c"

void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    if (strrchr(format, '\n') == NULL) fputc('\n', stderr);
    exit(-1);
}

static void test_bcf_set_variant_type(void)
{
    // Test SNVs
    bcf_variant_t var1;
    bcf_set_variant_type("A", "T", &var1);
    if ( var1.type != VCF_SNP)
    {
        error("A -> T was not detected as a SNP");
    }

    // Test INDEL
    bcf_variant_t var2a;
    bcf_set_variant_type("A", "AA", &var2a);
    if ( var2a.type != (VCF_INDEL|VCF_INS) )
    {
        error("A -> AA was not detected as an INDEL");
    }
    bcf_variant_t var2b;
    bcf_set_variant_type("AA", "A", &var2b);
    if ( var2b.type != (VCF_INDEL|VCF_DEL) )
    {
        error("AA -> A was not detected as a INDEL");
    }

    // Test breakends
    bcf_variant_t var3a;
    bcf_set_variant_type("N", "N]16:33625444]", &var3a);
    if ( var3a.type != VCF_BND)
    {
        error("N]16:33625444] was not detected as a breakend");
    }

    bcf_variant_t var3b;
    bcf_set_variant_type("N", "N[16:33625444[", &var3b);
    if (var3b.type != VCF_BND)
    {
        error("N[16:33625444[ was not detected as a breakend");
    }

    bcf_variant_t var3c;
    bcf_set_variant_type("N", "]16:33625444]N", &var3c);
    if ( var3c.type != VCF_BND)
    {
        error("]16:33625444]N was not detected as a breakend");
    }

    bcf_variant_t var3d;
    bcf_set_variant_type("N", "[16:33625444[N", &var3d);
    if ( var3d.type != VCF_BND)
    {
        error("[16:33625444[N was not detected as a breakend");
    }
    // Test special reference alleles
    bcf_variant_t var4a;
    bcf_set_variant_type("A", "<NON_REF>", &var4a);
    if ( var4a.type != VCF_REF)
    {
        error("<NON_REF> was not detected as a special reference allele");
    }
    bcf_variant_t var4b;
    bcf_set_variant_type("A", "<*>", &var4b);
    if ( var4b.type != VCF_REF)
    {
        error("<*> was not detected as a special reference allele");
    }
    // Test MNP
    bcf_variant_t var5;
    bcf_set_variant_type("AA", "TT", &var5);
    if ( var5.type != VCF_MNP)
    {
        error("AA->TT was not detected as a MNP");
    }
    // Test Overlapping allele
    bcf_variant_t var6;
    bcf_set_variant_type("A", "*", &var6);
    if ( var6.type != VCF_OVERLAP)
    {
        error("A->* was not detected as an overlap");
    }
    // Test .
    bcf_variant_t var7;
    bcf_set_variant_type("A", ".", &var7);
    if ( var7.type != VCF_REF)
    {
        error("A->. was not detected as a special reference allele");
    }
}

int main(int argc, char **argv)
{
    test_bcf_set_variant_type();
    return 0;
}

