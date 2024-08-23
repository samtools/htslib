/*
    Copyright (C) 2017, 2020, 2023 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Test bcf synced reader allele pairing
*/

#include <config.h>

#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include <stdarg.h>
#include <inttypes.h>
#include <strings.h>
#include <errno.h>

#include "../htslib/hts_defs.h"
#include "../htslib/synced_bcf_reader.h"
#include "../htslib/hts.h"
#include "../htslib/vcf.h"

void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

void HTS_NORETURN usage(int exit_code)
{
    fprintf(stderr, "Usage: test-bcf-sr [OPTIONS] vcf-list.txt\n");
    fprintf(stderr, "       test-bcf-sr [OPTIONS] -args file1.bcf [...]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "       --args                   pass filenames directly in argument list\n");
    fprintf(stderr, "       --no-index               allow streaming\n");
    fprintf(stderr, "   -o, --output <file>          output file (stdout if not set)\n");
    fprintf(stderr, "   -O, --output-fmt <fmt>       fmt: vcf,bcf,summary\n");
    fprintf(stderr, "   -p, --pair <logic[+ref]>     logic: snps,indels,both,snps+ref,indels+ref,both+ref,exact,some,all\n");
    fprintf(stderr, "   -r, --regions <reg_list>     comma-separated list of regions\n");
    fprintf(stderr, "   -t, --targets <reg_list>     comma-separated list of targets\n");
    fprintf(stderr, "\n");
    exit(exit_code);
}

void write_summary_format(bcf_srs_t *sr, FILE *out)
{
    int n, i, j;
    while ((n = bcf_sr_next_line(sr)) > 0) {
        for (i=0; i<sr->nreaders; i++)
        {
            if ( !bcf_sr_has_line(sr,i) ) continue;
            bcf1_t *rec = bcf_sr_get_line(sr, i);
            if (!rec) error("bcf_sr_get_line() unexpectedly returned NULL\n");
            fprintf(out, "%s:%"PRIhts_pos,
                    bcf_seqname_safe(bcf_sr_get_header(sr,i),rec),rec->pos+1);
            break;
        }

        for (i=0; i<sr->nreaders; i++)
        {
            fprintf(out, "\t");

            if ( !bcf_sr_has_line(sr,i) )
            {
                fprintf(out, "%s","-");
                continue;
            }

            bcf1_t *rec = bcf_sr_get_line(sr, i);
            if (!rec) error("bcf_sr_get_line() unexpectedly returned NULL\n");
            fprintf(out, "%s", rec->n_allele > 1 ? rec->d.allele[1] : ".");
            for (j=2; j<rec->n_allele; j++)
            {
                fprintf(out, ",%s", rec->d.allele[j]);
            }
        }
        fprintf(out, "\n");
    }
}

void write_vcf_bcf_format(bcf_srs_t *sr, bcf_hdr_t *hdr, vcfFile *vcf_out,
                          const char *fmt_type)
{
    int i, n;
    if (bcf_hdr_write(vcf_out, hdr) != 0)
        error("Couldn't write %s header\n", fmt_type);

    while ((n = bcf_sr_next_line(sr)) > 0) {
        for (i=0; i<sr->nreaders; i++)
        {
            if ( !bcf_sr_has_line(sr,i) ) continue;
            bcf1_t *rec = bcf_sr_get_line(sr, i);
            if (!rec) error("bcf_sr_get_line() unexpectedly returned NULL\n");
            if (vcf_write(vcf_out, hdr, rec) < 0)
                error("vcf_write() failed\n");
        }
    }
}

int main(int argc, char *argv[])
{
    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"output-fmt",required_argument,NULL,'O'},
        {"pair",required_argument,NULL,'p'},
        {"regions",required_argument,NULL,'r'},
        {"targets",required_argument,NULL,'t'},
        {"no-index",no_argument,NULL,1000},
        {"args",no_argument,NULL,1001},
        {NULL,0,NULL,0}
    };

    int c, pair = 0, use_index = 1, use_fofn = 1;
    enum htsExactFormat out_fmt = text_format; // for original pos + alleles
    const char *out_fn = NULL, *regions = NULL, *targets = NULL;
    while ((c = getopt_long(argc, argv, "o:O:p:r:t:h", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'o':
                out_fn = optarg;
                break;
            case 'O':
                if (!strcasecmp(optarg, "vcf"))          out_fmt = vcf;
                else if (!strcasecmp(optarg, "bcf"))     out_fmt = bcf;
                else if (!strcasecmp(optarg, "summary")) out_fmt = text_format;
                else error("Unknown output format \"%s\"\n", optarg);
                break;
            case 'p':
                if ( !strcmp(optarg,"snps") )            pair |= BCF_SR_PAIR_SNPS;
                else if ( !strcmp(optarg,"snp+ref") )    pair |= BCF_SR_PAIR_SNPS|BCF_SR_PAIR_SNP_REF;
                else if ( !strcmp(optarg,"snps+ref") )   pair |= BCF_SR_PAIR_SNPS|BCF_SR_PAIR_SNP_REF;
                else if ( !strcmp(optarg,"indels") )     pair |= BCF_SR_PAIR_INDELS;
                else if ( !strcmp(optarg,"indel+ref") )  pair |= BCF_SR_PAIR_INDELS|BCF_SR_PAIR_INDEL_REF;
                else if ( !strcmp(optarg,"indels+ref") ) pair |= BCF_SR_PAIR_INDELS|BCF_SR_PAIR_INDEL_REF;
                else if ( !strcmp(optarg,"both") )       pair |= BCF_SR_PAIR_BOTH;
                else if ( !strcmp(optarg,"both+ref") )   pair |= BCF_SR_PAIR_BOTH_REF;
                else if ( !strcmp(optarg,"any") )        pair |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"all") )        pair |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"some") )       pair |= BCF_SR_PAIR_SOME;
                else if ( !strcmp(optarg,"exact") )      pair  = BCF_SR_PAIR_EXACT;
                else error("The --pair logic \"%s\" not recognised.\n", optarg);
                break;
            case 'r':
                regions = optarg;
                break;
            case 't':
                targets = optarg;
                break;
            case 1000:
                use_index = 0;
                break;
            case 1001:
                use_fofn = 0;
                break;
            case 'h':
                usage(EXIT_SUCCESS);
            default: usage(EXIT_FAILURE);
        }
    }
    if ( !pair ) pair = BCF_SR_PAIR_EXACT;
    if ( optind == argc ) usage(EXIT_FAILURE);

    int i, nvcf;
    char **vcfs = NULL;
    if (use_fofn) {
        vcfs = hts_readlist(argv[optind], 1, &nvcf);
        if ( !vcfs ) error("Could not parse %s\n", argv[optind]);
    } else {
        vcfs = &argv[optind];
        nvcf = argc - optind;
    }

    bcf_srs_t *sr = bcf_sr_init();
    if (!sr) error("bcf_sr_init() failed\n");
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, pair);
    if (use_index) {
        bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    } else {
        bcf_sr_set_opt(sr, BCF_SR_ALLOW_NO_IDX);
    }

    if (regions)
    {
        if (bcf_sr_set_regions(sr, regions, 0) != 0)
            error("Failed to set regions\n");
    }

    if (targets)
    {
        if (bcf_sr_set_targets(sr, targets, 0, 0) != 0)
            error("Failed to set targets\n");
    }

    for (i=0; i<nvcf; i++)
        if ( !bcf_sr_add_reader(sr,vcfs[i]) ) error("Failed to open %s: %s\n", vcfs[i],bcf_sr_strerror(sr->errnum));

    if (!sr->readers || sr->nreaders < 1)
        error("No readers set, even though one was added\n");

    if (out_fmt == text_format) {
        FILE *out = stdout;
        if (out_fn)
        {
            out = fopen(out_fn, "w");
            if (!out) error("Couldn't open \"%s\" for writing: %s\n",
                            out_fn, strerror(errno));
        }
        write_summary_format(sr, out);
        if (out_fn)
        {
            if (fclose(out) != 0)
                error("Error on closing %s : %s\n",
                      out_fn, strerror(errno));
        }
    } else {
        const char *fmt_type = out_fmt == vcf ? "VCF" : "BCF";

        bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
        if (!hdr) error("%s output, but don't have a header\n", fmt_type);

        if (!out_fn) { out_fn = "-"; }
        vcfFile *vcf_out = vcf_open(out_fn, out_fmt == vcf ? "w" : "wb");
        if (!vcf_out) error("Couldn't open \"%s\" for writing: %s\n",
                            out_fn, strerror(errno));
        write_vcf_bcf_format(sr, hdr, vcf_out, fmt_type);
        if (vcf_close(vcf_out) != 0)
            error("Error on closing \"%s\"\n", out_fn);
    }

    if (sr->errnum) error("Synced reader error: %s\n",
                          bcf_sr_strerror(sr->errnum));

    bcf_sr_destroy(sr);
    if (use_fofn)
    {
        for (i=0; i<nvcf; i++)
            free(vcfs[i]);
        free(vcfs);
    }

    return 0;
}

