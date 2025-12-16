/*  test/test-vcf-api.c -- VCF test harness.

    Copyright (C) 2013, 2014, 2017-2021, 2023, 2025 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "../htslib/hts.h"
#include "../htslib/vcf.h"
#include "../htslib/vcfutils.h"
#include "../htslib/kbitset.h"
#include "../htslib/kstring.h"
#include "../htslib/kseq.h"

void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    if (strrchr(format, '\n') == NULL) fputc('\n', stderr);
    exit(-1);
}

#define STRINGIFY(x) #x
#define check0(x) ((x) == 0 ? (void) 0 : error("Failed: %s", STRINGIFY(x)))

static int check_alleles(bcf1_t *rec, const char **alleles, int num) {
    int i;
    if (rec->n_allele !=  num) {
        fprintf(stderr, "Wrong number of alleles - expected %d, got %d\n",
                num, rec->n_allele);
        return -1;
    }
    if (bcf_unpack(rec, BCF_UN_STR) != 0)
        return -1;
    for (i = 0; i < num; i++) {
        if (0 != strcmp(alleles[i], rec->d.allele[i])) {
            fprintf(stderr,
                    "Mismatch for allele %d : expected '%s' got '%s'\n",
                    i, alleles[i], rec->d.allele[i]);
            return -1;
        }
    }
    return 0;
}

static void test_update_alleles(bcf_hdr_t *hdr, bcf1_t *rec)
{
    // Exercise bcf_update_alleles() a bit
    const char *alleles1[2] = { "G", "A" };
    const char *alleles2[3] = { "C", "TGCA", "CATG" };
#define rep10(x) x x x x x x x x x x
    const char *alleles3[3] = { rep10("ATTCTAGATC"), "TGCA",
                                rep10("CTATTATCTCTAATGACATG") };
#undef rep10
    const char *alleles4[3] = { alleles3[2], NULL, alleles3[0] };
    // Add some alleles
    check0(bcf_update_alleles(hdr, rec, alleles1, 2));
    check0(check_alleles(rec, alleles1, 2));
    // Erase them
    check0(bcf_update_alleles(hdr, rec, NULL, 0));
    check0(check_alleles(rec, NULL, 0));
    // Expand to three
    check0(bcf_update_alleles(hdr, rec, alleles2, 3));
    check0(check_alleles(rec, alleles2, 3));
    // Now try some bigger ones (should force a realloc)
    check0(bcf_update_alleles(hdr, rec, alleles3, 3));
    check0(check_alleles(rec, alleles3, 3));
    // Ensure it works even if one of the alleles points into the
    // existing structure
    alleles4[1] = rec->d.allele[1];
    check0(bcf_update_alleles(hdr, rec, alleles4, 3));
    alleles4[1] = alleles3[1]; // Will have been clobbered by the update
    check0(check_alleles(rec, alleles4, 3));
    // Ensure it works when the alleles point into the existing data,
    // rec->d.allele is used to define the input array and the
    // order of the entries is changed.  The result of this should
    // be the same as alleles2.
    char *tmp = rec->d.allele[0] + strlen(rec->d.allele[0]) - 4;
    rec->d.allele[0] = rec->d.allele[2] + strlen(rec->d.allele[2]) - 1;
    rec->d.allele[2] = tmp;
    check0(bcf_update_alleles(hdr, rec, (const char **) rec->d.allele, 3));
    check0(check_alleles(rec, alleles2, 3));
}

void write_bcf(char *fname)
{
    // Init
    htsFile *fp    = hts_open(fname,"wb");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) error("bcf_hdr_init : %s", strerror(errno));
    bcf1_t *rec    = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));

    // Check no-op on fresh bcf1_t
    check0(bcf_update_alleles(hdr, rec, NULL, 0));

    // Create VCF header
    kstring_t str = {0,0,0};
    check0(bcf_hdr_append(hdr, "##fileDate=20090805"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=UF,Number=1,Type=Integer,Description=\"Unused FORMAT\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=UI,Number=1,Type=Integer,Description=\"Unused INFO\">"));
    check0(bcf_hdr_append(hdr, "##FILTER=<ID=Flt,Description=\"Unused FILTER\">"));
    check0(bcf_hdr_append(hdr, "##unused=<XX=AA,Description=\"Unused generic\">"));
    check0(bcf_hdr_append(hdr, "##unused=<ID=BB,Description=\"Unused generic with ID\">"));
    check0(bcf_hdr_append(hdr, "##unused=unformatted text 1"));
    check0(bcf_hdr_append(hdr, "##unused=unformatted text 2"));
    check0(bcf_hdr_append(hdr, "##contig=<ID=Unused,length=1>"));
    check0(bcf_hdr_append(hdr, "##source=myImputationProgramV3.1"));
    check0(bcf_hdr_append(hdr, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta"));
    check0(bcf_hdr_append(hdr, "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"));
    check0(bcf_hdr_append(hdr, "##phasing=partial"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=NEG,Number=.,Type=Integer,Description=\"Test -ve Numbers\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"));
    check0(bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">"));
    check0(bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than half of samples have data\">"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String 1\">"));

    // Try a few header modifications
    bcf_hdr_remove(hdr, BCF_HL_CTG, "Unused");
    check0(bcf_hdr_append(hdr, "##contig=<ID=Unused,length=62435964>"));
    bcf_hdr_remove(hdr, BCF_HL_FMT, "TS");
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String\">"));
    bcf_hdr_remove(hdr, BCF_HL_INFO, "NEG");
    check0(bcf_hdr_append(hdr, "##INFO=<ID=NEG,Number=.,Type=Integer,Description=\"Test Negative Numbers\">"));
    bcf_hdr_remove(hdr, BCF_HL_FLT, "s50");
    check0(bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">"));

    check0(bcf_hdr_add_sample(hdr, "NA00001"));
    check0(bcf_hdr_add_sample(hdr, "NA00002"));
    check0(bcf_hdr_add_sample(hdr, "NA00003"));
    check0(bcf_hdr_add_sample(hdr, NULL));      // to update internal structures
    if ( bcf_hdr_write(fp, hdr)!=0 ) error("Failed to write to %s\n", fname);


    // Add a record
    // 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;NEG=-127;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    // .. CHROM
    rec->rid = bcf_hdr_name2id(hdr, "20");
    // .. POS
    rec->pos = 14369;
    // .. ID
    check0(bcf_update_id(hdr, rec, "rs6054257"));
    // .. REF and ALT
    test_update_alleles(hdr, rec);
    const char *alleles[2] = { "G", "A" };
    check0(bcf_update_alleles_str(hdr, rec, "G,A"));
    check0(check_alleles(rec, alleles, 2));
    // .. QUAL
    rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    check0(bcf_update_filter(hdr, rec, &tmpi, 1));
    // .. INFO
    tmpi = 3;
    check0(bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1));
    tmpi = 500;
    check0(bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1));
    tmpi = 100000;
    check0(bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1));
    tmpi = 14;
    check0(bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1));
    tmpi = -127;
    check0(bcf_update_info_int32(hdr, rec, "NEG", &tmpi, 1));
    float tmpf = 0.5;
    check0(bcf_update_info_float(hdr, rec, "AF", &tmpf, 1));
    check0(bcf_update_info_flag(hdr, rec, "DB", NULL, 1));
    check0(bcf_update_info_flag(hdr, rec, "H2", NULL, 1));
    // .. FORMAT
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_phased(0);
    tmpia[1] = bcf_gt_phased(0);
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_gt_phased(0);
    tmpia[4] = bcf_gt_unphased(1);
    tmpia[5] = bcf_gt_unphased(1);
    check0(bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2));
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    check0(bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr)));
    tmpia[0] = 0;
    tmpia[1] = 0;
    tmpia[2] = 1;
    check0(bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr)));
    tmpia[0] = 1;
    tmpia[1] = 100000;
    tmpia[2] = 1;
    check0(bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr)));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    check0(bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr)));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    check0(bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2));
    char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    check0(bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3));
    tmp_str[0] = "LongerStringRequiringBufferReallocation";
    check0(bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3));
    tmp_str[0] = "String1";
    check0(bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3));
    if ( bcf_write1(fp, hdr, rec)!=0 ) error("Failed to write to %s\n", fname);

    // 20     1110696 . A      G,T     67   .   NS=2;DP=10;NEG=-128;AF=0.333,.;AA=T;DB GT 2 1   ./.
    bcf_clear1(rec);
    rec->rid = bcf_hdr_name2id(hdr, "20");
    rec->pos = 1110695;
    check0(bcf_update_alleles_str(hdr, rec, "A,G,T"));
    rec->qual = 67;
    tmpi = 2;
    check0(bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1));
    tmpi = 10;
    check0(bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1));
    tmpi = -128;
    check0(bcf_update_info_int32(hdr, rec, "NEG", &tmpi, 1));
    float *tmpfa = (float*)malloc(2*sizeof(float));
    tmpfa[0] = 0.333;
    bcf_float_set_missing(tmpfa[1]);
    check0(bcf_update_info_float(hdr, rec, "AF", tmpfa, 2));
    check0(bcf_update_info_string(hdr, rec, "AA", "SHORT"));
    check0(bcf_update_info_string(hdr, rec, "AA", "LONGSTRING"));
    check0(bcf_update_info_string(hdr, rec, "AA", "T"));
    check0(bcf_update_info_flag(hdr, rec, "DB", NULL, 1));
    tmpia[0] = bcf_gt_phased(2);
    tmpia[1] = bcf_int32_vector_end;
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_int32_vector_end;
    tmpia[4] = bcf_gt_missing;
    tmpia[5] = bcf_gt_missing;
    check0(bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2));
    if ( bcf_write1(fp, hdr, rec)!=0 ) error("Failed to write to %s\n", fname);

    free(tmpia);
    free(tmpfa);

    // Clean
    free(str.s);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}

void bcf_to_vcf(char *fname)
{
    htsFile *fp    = hts_open(fname,"rb");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
    bcf1_t *rec    = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));

    char *gz_fname = (char*) malloc(strlen(fname)+4);
    if (!gz_fname) error("malloc : %s", strerror(errno));
    snprintf(gz_fname,strlen(fname)+4,"%s.gz",fname);
    htsFile *out   = hts_open(gz_fname,"wg");
    if (!out) error("Couldn't open \"%s\" : %s\n", gz_fname, strerror(errno));

    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    if (!bcf_hdr_get_hrec(hdr_out, BCF_HL_STR,"ID","BB","unused"))
        error("Missing header ##unused=<ID=BB, ...>");
    bcf_hdr_remove(hdr_out,BCF_HL_STR,"BB");
    if (bcf_hdr_get_hrec(hdr_out, BCF_HL_STR,"ID","BB","unused"))
        error("Got pointer to deleted header ##unused=<ID=BB, ...>");

    if (!bcf_hdr_get_hrec(hdr_out,BCF_HL_GEN,"unused","unformatted text 1",NULL))
        error("Missing header ##unused=unformatted text 1");
    bcf_hdr_remove(hdr_out,BCF_HL_GEN,"unused");
    if (bcf_hdr_get_hrec(hdr_out,BCF_HL_GEN,"unused","unformatted text 1",NULL))
        error("Got pointer to deleted header ##unused=unformatted text 1");

    if (!bcf_hdr_get_hrec(hdr_out,BCF_HL_FLT,"ID","Flt",NULL))
        error("Missing header ##FILTER=<ID=Flt, ...>");
    bcf_hdr_remove(hdr_out,BCF_HL_FLT,"Flt");
    if (bcf_hdr_get_hrec(hdr_out,BCF_HL_FLT,"ID","Flt",NULL))
        error("Got pointer to deleted header ##FILTER=<ID=Flt, ...>");

    if (!bcf_hdr_get_hrec(hdr_out,BCF_HL_INFO,"ID","UI",NULL))
        error("Missing header ##INFO=<ID=UI, ...>");
    bcf_hdr_remove(hdr_out,BCF_HL_INFO,"UI");
    if (bcf_hdr_get_hrec(hdr_out,BCF_HL_INFO,"ID","UI",NULL))
        error("Got pointer to deleted header ##INFO=<ID=UI, ...>");

    if (!bcf_hdr_get_hrec(hdr_out,BCF_HL_FMT,"ID","UF",NULL))
        error("Missing header ##INFO=<ID=UF, ...>");
    bcf_hdr_remove(hdr_out,BCF_HL_FMT,"UF");
    if (bcf_hdr_get_hrec(hdr_out,BCF_HL_FMT,"ID","UF",NULL))
        error("Got pointer to deleted header ##INFO=<ID=UF, ...>");

    if (!bcf_hdr_get_hrec(hdr_out,BCF_HL_CTG,"ID","Unused",NULL))
        error("Missing header ##contig=<ID=Unused,length=1>");
    bcf_hdr_remove(hdr_out,BCF_HL_CTG,"Unused");
    if (bcf_hdr_get_hrec(hdr_out,BCF_HL_FMT,"ID","Unused",NULL))
        error("Got pointer to header ##contig=<ID=Unused,length=1>");

    if ( bcf_hdr_write(out, hdr_out)!=0 ) error("Failed to write to %s\n", fname);
    int r;
    while ((r = bcf_read1(fp, hdr, rec)) >= 0)
    {
        if ( bcf_write1(out, hdr_out, rec)!=0 ) error("Failed to write to %s\n", fname);

        // Test problems caused by bcf1_sync: the data block
        // may be realloced, also the unpacked structures must
        // get updated.
        check0(bcf_unpack(rec, BCF_UN_STR));
        check0(bcf_update_id(hdr, rec, 0));
        check0(bcf_update_format_int32(hdr, rec, "GQ", NULL, 0));

        bcf1_t *dup = bcf_dup(rec);     // force bcf1_sync call
        if ( bcf_write1(out, hdr_out, dup)!=0 ) error("Failed to write to %s\n", fname);
        bcf_destroy1(dup);

        check0(bcf_update_alleles_str(hdr_out, rec, "G,A"));
        int32_t tmpi = 99;
        check0(bcf_update_info_int32(hdr_out, rec, "DP", &tmpi, 1));
        int32_t tmpia[] = {9,9,9};
        check0(bcf_update_format_int32(hdr_out, rec, "DP", tmpia, 3));

        if ( bcf_write1(out, hdr_out, rec)!=0 ) error("Failed to write to %s\n", fname);
    }
    if (r < -1) error("bcf_read1");

    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
    if ( (ret=hts_close(out)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",gz_fname,ret);
        exit(ret);
    }


    // read gzip, write stdout
    htsFile *gz_in = hts_open(gz_fname, "r");
    if ( !gz_in )
    {
        fprintf(stderr,"Could not read: %s\n", gz_fname);
        exit(1);
    }

    kstring_t line = {0,0,0};
    while ( hts_getline(gz_in, KS_SEP_LINE, &line)>0 )
    {
        kputc('\n',&line);
        fwrite(line.s,1,line.l,stdout);
    }

    if ( (ret=hts_close(gz_in)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",gz_fname,ret);
        exit(ret);
    }
    free(line.s);
    free(gz_fname);
}

void iterator(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
    hts_idx_t *idx;
    hts_itr_t *iter;

    bcf_index_build(fname, 0);
    idx = bcf_index_load(fname);

    iter = bcf_itr_queryi(idx, bcf_hdr_name2id(hdr, "20"), 1110600, 1110800);
    bcf_itr_destroy(iter);

    iter = bcf_itr_querys(idx, hdr, "20:1110600-1110800");
    bcf_itr_destroy(iter);

    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}

void test_get_info_values(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
    bcf1_t *line = bcf_init();
    if (!line) error("bcf_init : %s", strerror(errno));
    int r;
    while ((r = bcf_read(fp, hdr, line)) == 0)
    {
        float *afs = 0;
        int32_t *negs = NULL;
        int count = 0;
        int ret = bcf_get_info_float(hdr, line, "AF", &afs, &count);

        if (line->pos == 14369)
        {
            if (ret != 1 || afs[0] != 0.5f)
            {
                fprintf(stderr, "AF on position 14370 should be 0.5\n");
                exit(-1);
            }
        }
        else
        {
            if (ret != 2 || afs[0] != 0.333f || !bcf_float_is_missing(afs[1]))
            {
                fprintf(stderr, "AF on position 1110696 should be 0.333, missing\n");
                exit(-1);
            }
        }

        free(afs);

        int32_t expected = (line->pos == 14369)? -127 : -128;
        count = 0;
        ret = bcf_get_info_int32(hdr, line, "NEG", &negs, &count);
        if (ret != 1 || negs[0] != expected)
        {
            if (ret < 0)
                fprintf(stderr, "NEG should be %d, got error ret=%d\n", expected, ret);
            else if (ret == 0)
                fprintf(stderr, "NEG should be %d, got no entries\n", expected);
            else
                fprintf(stderr, "NEG should be %d, got %d entries (first is %d)\n", expected, ret, negs[0]);
            exit(1);
        }
        free(negs);
    }
    if (r < -1) error("bcf_read");

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

void write_format_values(const char *fname)
{
    // Init
    htsFile *fp = hts_open(fname, "wb");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) error("bcf_hdr_init : %s", strerror(errno));
    bcf1_t *rec = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));

    // Create VCF header
    check0(bcf_hdr_append(hdr, "##contig=<ID=1>"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=TF,Number=1,Type=Float,Description=\"Test Float\">"));
    check0(bcf_hdr_add_sample(hdr, "S"));
    check0(bcf_hdr_add_sample(hdr, NULL)); // to update internal structures
    if ( bcf_hdr_write(fp, hdr)!=0 ) error("Failed to write to %s\n", fname);

    // Add a record
    // .. FORMAT
    float test[4];
    bcf_float_set_missing(test[0]);
    test[1] = 47.11f;
    bcf_float_set_vector_end(test[2]);
    test[3] = -1.2e-13;
    check0(bcf_update_format_float(hdr, rec, "TF", test, 4));
    if ( bcf_write1(fp, hdr, rec)!=0 ) error("Failed to write to %s\n", fname);

    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ((ret = hts_close(fp)))
    {
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", fname, ret);
        exit(ret);
    }
}

void check_format_values(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *line = bcf_init();

    while (bcf_read(fp, hdr, line) == 0)
    {
        float *values = 0;
        int count = 0;
        int ret = bcf_get_format_float(hdr, line, "TF", &values, &count);

        // NOTE the return value from bcf_get_format_float is different from
        // bcf_get_info_float in the sense that vector-end markers also count.
        if (ret != 4 ||
            count < ret ||
            !bcf_float_is_missing(values[0]) ||
            values[1] != 47.11f ||
            !bcf_float_is_vector_end(values[2]) ||
            !bcf_float_is_vector_end(values[3]))
        {
            fprintf(stderr, "bcf_get_format_float didn't produce the expected output.\n");
            exit(-1);
        }

        free(values);
    }

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

void test_get_format_values(const char *fname)
{
    write_format_values(fname);
    check_format_values(fname);
}

void test_invalid_end_tag(void)
{
    static const char vcf_data[] = "data:,"
        "##fileformat=VCFv4.1\n"
        "##contig=<ID=X,length=155270560>\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "X\t86470037\trs59780433a\tTTTCA\tTGGTT,T\t.\t.\tEND=85725113\n"
        "X\t86470038\trs59780433b\tT\tTGGTT,T\t.\t.\tEND=86470047\n";

    htsFile *fp;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    int ret;
    int32_t tmpi;
    enum htsLogLevel logging = hts_get_log_level();

    // Silence warning messages
    hts_set_log_level(HTS_LOG_ERROR);

    fp = hts_open(vcf_data, "r");
    if (!fp) error("Failed to open vcf data : %s", strerror(errno));
    rec = bcf_init1();
    if (!rec) error("Failed to allocate BCF record : %s", strerror(errno));

    hdr = bcf_hdr_read(fp);
    if (!hdr) error("Failed to read BCF header : %s", strerror(errno));

    check0(bcf_read(fp, hdr, rec));
    // rec->rlen should ignore the bogus END tag value on the first read
    if (rec->rlen != 5) {
        error("Incorrect rlen - expected 5 got %"PRIhts_pos"\n", rec->rlen);
    }

    check0(bcf_read(fp, hdr, rec));
    // While on the second it should use it
    if (rec->rlen != 10) {
        error("Incorrect rlen - expected 10 got %"PRIhts_pos"\n", rec->rlen);
    }

    // Try to break it - will change rlen
    tmpi = 85725113;
    check0(bcf_update_info_int32(hdr, rec, "END", &tmpi, 1));

    if (rec->rlen != 1) {
        error("Incorrect rlen - expected 1 got %"PRIhts_pos"\n", rec->rlen);
    }

    ret = bcf_read(fp, hdr, rec);
    if (ret != -1) {
        error("Unexpected return code %d from bcf_read at EOF", ret);
    }

    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    ret = hts_close(fp);
    if (ret != 0) {
        error("Unexpected return code %d from hts_close", ret);
    }

    hts_set_log_level(logging);
}

void test_open_format(void) {
    char mode[5];
    int ret;
    strcpy(mode, "r");
    ret = vcf_open_mode(mode+1, "mode1.bcf", NULL);
    if (strncmp(mode, "rb", 2) || ret)
        error("Mode '%s' does not match the expected value '%s'", mode, "rb");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.vcf", NULL);
    if (strncmp(mode, "r", 1) || ret)
        error("Mode '%s' does not match the expected value '%s'", mode, "r");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.vcf.gz", NULL);
    if (strncmp(mode, "rz", 2) || ret)
        error("Mode '%s' does not match the expected value '%s'", mode, "rz");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.vcf.bgz", NULL);
    if (strncmp(mode, "rz", 2) || ret)
        error("Mode '%s' does not match the expected value '%s'", mode, "rz");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.xcf", NULL);
    if (!ret)
        error("Expected failure for wrong extension 'xcf'");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.vcf.gbz", NULL);
    if (!ret)
        error("Expected failure for wrong extension 'vcf.gbz'");
    mode[1] = 0;
    ret = vcf_open_mode(mode+1, "mode1.bvcf.bgz", NULL);
    if (!ret)
        error("Expected failure for wrong extension 'vcf.bvcf.bgz'");
}

//tests on rlen calculation
void test_rlen_values(void)
{
    bcf_hdr_t *hdr;
    bcf1_t *rec, *rec2;
    int i, j;
    int32_t tmpi;
    //data common for all versions, interpreted differently based on version
#define data "##reference=file://tmp\n" \
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" \
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n" \
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n" \
    "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"svlen\">\n" \
    "##INFO=<ID=CN,Number=A,Type=Float,Description=\"Copy number\">\n" \
    "##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"svclaim\">\n" \
    "##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"fmt len\">\n" \
    "##contig=<ID=1,Length=40>\n" \
    "##ALT=<ID=INS,Description=\"INS\">\n" \
    "##ALT=<ID=DEL,Description=\"DEL\">\n" \
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n" \
    "1\t4310\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0|1\n" \
    "1\t4311\t.\tC\tCT\t.\t.\t.\tGT\t0/0\t1/0\n" \
    "1\t4312\t.\tTC\tT\t213.73\t.\t.\tGT\t0/1\t0|0\n" \
    "1\t4314\t.\tG\t<INS>\t213.73\t.\tSVLEN=10;SVCLAIM=J\tGT\t0/1\t0|0\n" \
    "1\t4315\t.\tG\t<DEL>\t213.73\t.\tSVLEN=-10;SVCLAIM=D\tGT\t0/1\t0|0\n" \
    "1\t4326\t.\tG\t<INS>\t213.73\t.\tEND=4326;SVLEN=10;SVCLAIM=J\tGT\t0/1\t0|0\n" \
    "1\t4327\t.\tG\t<DEL>\t213.73\t.\tEND=4337;SVLEN=-10;SVCLAIM=J\tGT\t0/1\t0|0\n" \
    "1\t4338\t.\tG\t<*>\t213.73\t.\tEND=4342;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:8\n" \
    "1\t4353\t.\tG\t<*>\t213.73\t.\tEND=4357;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n" \
    "1\t4363\t.\tG\t<*>\t213.73\t.\tEND=4367;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n" \
    "1\t4370\t.\tG\t<INS>,<*>\t213.73\t.\tEND=4371;SVLEN=.;SVCLAIM=.\tGT:LEN\t0/1:7\t0|0:.\n" \
    "1\t4378\t.\tG\t<DEL>,<INS>,<*>\t213.73\t.\tEND=4379;SVLEN=3,5,.;SVCLAIM=D,J,.\tGT:LEN\t0/1:7\t0|0:.\n" \
    "1\t4385\t.\tG\tT,<DEL>\t213.73\t.\tEND=4387;SVLEN=.,180\tGT\t0/1\t0|0\n" \
    "1\t4585\t.\tG\tT,<DEL:ME>\t213.73\t.\tEND=4587;SVLEN=.,180\tGT\t0/1\t0|0\n" \
    "1\t4685\t.\tG\t<DUP>,<DUP>\t213.73\t.\tEND=4687;SVLEN=10,10\tGT\t0/1\t0|0\n" \
    "1\t4705\t.\tG\t<CNV>\t213.73\t.\tEND=4707;SVLEN=11;CN=2\tGT\t0/1\t0|0\n" \
    "1\t4725\t.\tG\t<CNV:TR>\t213.73\t.\tEND=4727;SVLEN=12;CN=1.5\tGT\t0/1\t0|0\n" \
    "1\t4745\t.\tG\t<INV>\t213.73\t.\tEND=4747;SVLEN=10\tGT\t0/1\t0|0\n" \
    "1\t4885\t.\tG\tT,<*>\t213.73\t.\tEND=4887\tGT:LEN\t0/1:190\t0|0:.\n" \
    "1\t5885\t.\tG\tT\t213.73\t.\tEND=5887;SVLEN=8;SVCLAIM=.\tGT:LEN\t0/1:.\t0|0:10\n"
    // For the last line, SVLEN should be ignored as the allele is not symbolic
    // and LEN should be ignored as there is no <*> allele.  The END tag will
    // be used as it's higher than POS + length of REF.

    //test vcf with different versions
    const int testsz = 3;
    static char d43[] = "data:,"
    "##fileformat=VCFv4.3\n" data;
    static char d44[] = "data:,"
    "##fileformat=VCFv4.4\n" data;
    static char d45[] = "data:,"
    "##fileformat=VCFv4.5\n" data;
    /* ideal expected rlen for tests
    int rlen43[] = {1, 1, 2, 1, 1, 1, 11, 5, 5, 5, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};
    int rlen44[] = {1, 1, 2, 1, 11, 1, 11, 5, 5, 5, 2, 4, 181, 181, 11, 12, 13, 11, 3, 3};
    int rlen45[] = {1, 1, 2, 1, 11, 1, 11, 8, 7, 7, 7, 7, 181, 181, 11, 12, 13, 11, 190, 3};*/
    // but we dont check version and hence resulting rlen should be
    int rlen[] = {1, 1, 2, 1, 11, 1, 11, 8, 7, 7, 7, 7, 181, 181, 11, 12, 13, 11, 190, 3};

    const int vcfsz = sizeof(rlen)/sizeof(rlen[0]);

    char *darr[] = {&d43[0], &d44[0], &d45[0]};     //data array
    int *rarr[] = {&rlen[0], &rlen[0], &rlen[0]};   //result array

    enum htsLogLevel logging = hts_get_log_level();

    // Silence warning messages
    hts_set_log_level(HTS_LOG_ERROR);

    rec = bcf_init1(); rec2 = bcf_init1();
    if (!rec || !rec2) error("Failed to allocate BCF record : %s", strerror(errno));
    //calculating rlen with different vcf versions
    for (i = 0; i < testsz; ++i) {
        htsFile *fp = hts_open(darr[i], "r");
        if (!fp) error("Failed to open vcf data with %d : %s", i + 1, strerror(errno));
        bcf_clear1(rec);
        hdr = bcf_hdr_read(fp);
        if (!hdr) error("Failed to read BCF header with %d : %s", i + 1, strerror(errno));
        for (j = 0; j < vcfsz; ++j) {
            check0(bcf_read(fp, hdr, rec));
            if (rec->rlen != rarr[i][j]) {
                error("Incorrect rlen @ vcf %d on test %d - expected %d got %"PRIhts_pos"\n", j + 1, i + 1, rarr[i][j], rec->rlen);
            }
        }
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }

    //calculating rlen with update to fields
    htsFile *fp = hts_open(d45, "r");
    int id  = 1, val[]= { 1, 15 };
    bcf_clear1(rec);
    hdr = bcf_hdr_read(fp);
    if (!hdr) error("Failed to read BCF header : %s", strerror(errno));
    check0(bcf_read(fp, hdr, rec));
    //"1\t4310\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0|1\n"
    if (rec->rlen != 1)
        error("Incorrect rlen set, expected 1 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    //alt change A->AT
    ++id; check0(bcf_update_alleles_str(hdr, rec, "G,AT"));
    if (rec->rlen != 1)
        error("Incorrect rlen set, expected 1 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    //ref change G->GC, alt AT -> A, "1\t4310\t.\tGC\tA\t.\t.\t.\tGT\t0/0\t0|1\n"
    id++; check0(bcf_update_alleles_str(hdr, rec, "GC,A"));
    if (rec->rlen != 2)
        error("Incorrect rlen set, expected 2 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    ++id; check0(bcf_update_alleles_str(hdr, rec, "G,<*>"));
    if (rec->rlen != 1)
        error("Incorrect rlen set, expected 1 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    tmpi = 4323;
    ++id; check0(bcf_update_info_int32(hdr, rec, "END", &tmpi, 1));  //SVLEN to infer from END and use
    if (rec->rlen != 14)
        error("Incorrect rlen set, expected 14 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    ++id; check0(bcf_update_format_int32(hdr, rec, "LEN", &val, 2));
    if (rec->rlen != 15)
        error("Incorrect rlen set, expected 15 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    ++id; check0(bcf_update_info_int32(hdr, rec, "END", &tmpi, 0));  //removes END
    if (rec->rlen != 15)
        error("Incorrect rlen set, expected 15 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    ++id; check0(bcf_update_format_int32(hdr, rec, "LEN", &tmpi, 0));  //removes LEN
    if (rec->rlen != 1)
        error("Incorrect rlen set, expected 1 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    //ALT -> DEL, "1\t4310\t.\tG\tT,<DEL>\t.\t.\t.\tGT\t0/0\t0|1\n"
    ++id; check0(bcf_update_alleles_str(hdr, rec, "G,T,<DEL>"));
    if (rec->rlen != 1)
        error("Incorrect rlen set, expected 1 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    val[0] = 0; val[1] = -5;
    ++id; check0(bcf_update_info_int32(hdr, rec, "SVLEN", &val, 2));  //Add svlen
    if (rec->rlen != 6)
        error("Incorrect rlen set, expected 6 got %"PRIhts_pos" @ %d\n", rec->rlen, id);
    ++id; bcf_copy(rec2, rec);
    if (rec2->rlen != 6)
        error("Incorrect rlen set, expected 6 got %"PRIhts_pos" @ %d\n", rec->rlen, id);

    //needs update when header version change is handled
    bcf_destroy1(rec);
    bcf_destroy1(rec2);
    bcf_hdr_destroy(hdr);

    hts_close(fp);

    hts_set_log_level(logging);
}

static void test_vl_types(void)
{
    const char *test_vcf = "data:,"
        "##fileformat=VCFv4.5\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=5,length=243199373>\n"
        "##INFO=<ID=FIXED_1_INFO,Number=1,Type=Integer,Description=\"Fixed number 1\">\n"
        "##INFO=<ID=FIXED_4_INFO,Number=4,Type=Float,Description=\"Fixed number 4\">\n"
        "##INFO=<ID=VL_DOT_INFO,Number=.,Type=Integer,Description=\"Variable number\">\n"
        "##INFO=<ID=VL_A_INFO,Number=A,Type=Integer,Description=\"One value for each ALT allele\">\n"
        "##INFO=<ID=VL_G_INFO,Number=G,Type=Integer,Description=\"One value for each possible genotype\">\n"
        "##INFO=<ID=VL_R_INFO,Number=R,Type=Integer,Description=\"One value for each allele including REF\">\n"
        "##FORMAT=<ID=FIXED_1_FMT,Number=1,Type=String,Description=\"Fixed number 1\">\n"
        "##FORMAT=<ID=FIXED_4_FMT,Number=4,Type=String,Description=\"Fixed number 4\">\n"
        "##FORMAT=<ID=VL_DOT_FMT,Number=.,Type=String,Description=\"Variable number\">\n"
        "##FORMAT=<ID=VL_A_FMT,Number=A,Type=Integer,Description=\"One value for each ALT allele\">\n"
        "##FORMAT=<ID=VL_G_FMT,Number=G,Type=Integer,Description=\"One value for each possible genotype\">\n"
        "##FORMAT=<ID=VL_R_FMT,Number=R,Type=Integer,Description=\"One value for each allele including REF\">\n"
        "##FORMAT=<ID=VL_P_FMT,Number=P,Type=String,Description=\"One value for each allele value defined in GT\">\n"
        "##FORMAT=<ID=VL_LA_FMT,Number=LA,Type=Integer,Description=\"One value for each local ALT allele\">\n"
        "##FORMAT=<ID=VL_LG_FMT,Number=LG,Type=Integer,Description=\"One value for each local genotype\">\n"
        "##FORMAT=<ID=VL_LR_FMT,Number=LR,Type=Integer,Description=\"One value for each local allele including REF\">\n"
        "##FORMAT=<ID=VL_M_FMT,Number=M,Type=Integer,Description=\"One value for each posible base modification of the given type\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAAA\n";

    typedef struct {
        const char *id;
        int type;
        int expected_vl_code;
        int expected_number;
    } expected_types_t;

    expected_types_t expected[] = {
        { "FIXED_1_INFO", BCF_HL_INFO, BCF_VL_FIXED, 1 },
        { "FIXED_4_INFO", BCF_HL_INFO, BCF_VL_FIXED, 4 },
        { "VL_DOT_INFO",  BCF_HL_INFO, BCF_VL_VAR,   0xfffff },
        { "VL_A_INFO",    BCF_HL_INFO, BCF_VL_A,     0xfffff },
        { "VL_G_INFO",    BCF_HL_INFO, BCF_VL_G,     0xfffff },
        { "VL_R_INFO",    BCF_HL_INFO, BCF_VL_R,     0xfffff },
        { "FIXED_1_FMT",  BCF_HL_FMT,  BCF_VL_FIXED, 1 },
        { "FIXED_4_FMT",  BCF_HL_FMT,  BCF_VL_FIXED, 4 },
        { "VL_DOT_FMT",   BCF_HL_FMT,  BCF_VL_VAR,   0xfffff },
        { "VL_A_FMT",     BCF_HL_FMT,  BCF_VL_A,     0xfffff },
        { "VL_G_FMT",     BCF_HL_FMT,  BCF_VL_G,     0xfffff },
        { "VL_R_FMT",     BCF_HL_FMT,  BCF_VL_R,     0xfffff },
        { "VL_P_FMT",     BCF_HL_FMT,  BCF_VL_P,     0xfffff },
        { "VL_LA_FMT",    BCF_HL_FMT,  BCF_VL_LA,    0xfffff },
        { "VL_LG_FMT",    BCF_HL_FMT,  BCF_VL_LG,    0xfffff },
        { "VL_LR_FMT",    BCF_HL_FMT,  BCF_VL_LR,    0xfffff },
        { "VL_M_FMT",     BCF_HL_FMT,  BCF_VL_M,     0xfffff },
    };

    htsFile *fp = hts_open(test_vcf, "r");
    if (!fp) error("Failed to open test data\n");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("Failed to read BCF header\n");
    size_t i;
    for (i = 0; i < sizeof(expected)/sizeof(expected[0]); i++)
    {
        int id_num = bcf_hdr_id2int(hdr, BCF_DT_ID, expected[i].id);
        if (id_num < 0)
        {
            error("Couldn't look up VCF header ID %s\n", expected[i].id);
        }
        int vl_code = bcf_hdr_id2length(hdr, expected[i].type, id_num);
        if (vl_code != expected[i].expected_vl_code)
        {
            const char *length_types[] = {
                "BCF_VL_FIXED", "BCF_VL_VAR", "BCF_VL_A", "BCF_VL_G",
                "BCF_VL_R", "BCF_VL_P", "BCF_VL_LA", "BCF_VL_LG",
                "BCF_VL_LR", "BCF_VL_M"
            };
            if (vl_code >= 0 && vl_code < sizeof(length_types)/sizeof(length_types[0]))
            {
                error("Unexpected length code for %s: expected %s got %s\n",
                      expected[i].id,
                      length_types[expected[i].expected_vl_code],
                      length_types[vl_code]);
            }
            else
            {
                error("Unexpected length code for %s: expected %s got %d\n",
                      expected[i].id,
                      length_types[expected[i].expected_vl_code], vl_code);
            }
        }
        int num = bcf_hdr_id2number(hdr, expected[i].type, id_num);
        if (num != expected[i].expected_number)
        {
            error("Unexpected number for %s: expected %d%s got %d%s\n",
                  expected[i].id,
                  expected[i].expected_number,
                  (expected[i].expected_number == 0xfffff
                   ? " (= code for not fixed)" : ""),
                  num,
                  num == 0xfffff ? " (= code for not fixed)" : "");
        }
    }
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}

static int read_vcf_line(const char *line, bcf_hdr_t *hdr, bcf1_t *rec,
                         kstring_t *kstr)
{
    int ret;
    if (kputsn(line, strlen(line), ks_clear(kstr)) < 0)
        return -1;

    ret = vcf_parse(kstr, hdr, rec);
    if (ret != 0) {
        fprintf(stderr, "vcf_parse(\"%s\", hdr, rec) returned %d\n",
                ks_c_str(kstr), ret);
    }
    return ret;
}

static void chomp(kstring_t *kstr)
{
    if (kstr->l < 1)
        return;
    if (kstr->s[kstr->l - 1] == '\n')
        kstr->s[--kstr->l] = '\0';
}

// Test bcf_remove_allele_set()
void test_bcf_remove_allele_set(void)
{
    char header[] = "##fileformat=VCFv4.5\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=5,length=243199373>\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
        "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
        "##INFO=<ID=CN,Number=A,Type=Float,Description=\"Copy number of CNV/breakpoint\">\n"
        "##INFO=<ID=CICN,Number=.,Type=Float,Description=\"Confidence interval around copy number\">\n"
        "##INFO=<ID=CIEND,Number=.,Type=Integer,Description=\"Confidence interval around the inferred END for symbolic structural variants\">\n"
        "##INFO=<ID=CILEN,Number=.,Type=Integer,Description=\"Confidence interval for the SVLEN field\">\n"
        "##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"Confidence interval around POS for symbolic structural variants\">\n"
        "##INFO=<ID=CIRB,Number=.,Type=Integer,Description=\"Confidence interval around RB\">\n"
        "##INFO=<ID=CIRUC,Number=.,Type=Float,Description=\"Confidence interval around RUC\">\n"
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n"
        "##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">\n"
        "##INFO=<ID=METRANS,Number=.,Type=String,Description=\"Mobile element transduction info of the form CHR,START,END,POLARITY\">\n"
        "##INFO=<ID=RB,Number=.,Type=Integer,Description=\"Total number of bases in the corresponding repeat sequence\">\n"
        "##INFO=<ID=RN,Number=A,Type=Integer,Description=\"Total number of repeat sequences in this allele\">\n"
        "##INFO=<ID=RUB,Number=.,Type=Integer,Description=\"Number of bases in each individual repeat unit\">\n"
        "##INFO=<ID=RUC,Number=.,Type=Float,Description=\"Repeat unit count of corresponding repeat sequence\">\n"
        "##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"Repeat unit length of the corresponding repeat sequence\">\n"
        "##INFO=<ID=RUS,Number=.,Type=String,Description=\"Repeat unit sequence of the corresponding repeat sequence\">\n"
        "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Length of structural variant\">\n"
        "##INFO=<ID=VL_A_STR_INFO,Number=A,Type=String,Description=\"INFO string Number=A\">\n"
        "##INFO=<ID=VL_R_STR_INFO,Number=R,Type=String,Description=\"INFO string Number=R\">\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">\n"
        "##FORMAT=<ID=EC,Number=A,Type=Integer,Description=\"Expected allele count\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
        "##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"Local alleles\">\n"
        "##FORMAT=<ID=LAD,Number=LR,Type=Integer,Description=\"Local allele depth\">\n"
        "##FORMAT=<ID=LEC,Number=LA,Type=Integer,Description=\"Local expected allele count\">\n"
        "##FORMAT=<ID=LPL,Number=LG,Type=Integer,Description=\"List of Phred-scaled local genotype likelihoods\">\n"
        "##FORMAT=<ID=VL_A_STR_FMT,Number=A,Type=String,Description=\"FMT string Number=A\">\n"
        "##FORMAT=<ID=VL_G_STR_FMT,Number=G,Type=String,Description=\"FMT string Number=G\">\n"
        "##FORMAT=<ID=VL_LA_STR_FMT,Number=LA,Type=String,Description=\"FMT string Number=LA\">\n"
        "##FORMAT=<ID=VL_LG_STR_FMT,Number=LG,Type=String,Description=\"FMT string Number=LG\">\n"
        "##FORMAT=<ID=VL_LR_STR_FMT,Number=LR,Type=String,Description=\"FMT string Number=LR\">\n"
        "##FORMAT=<ID=VL_R_STR_FMT,Number=R,Type=String,Description=\"FMT string Number=R\">\n"
        "##ALT=<ID=CNV:TR,Description=\"Tandem repeat determined based on DNA abundance\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAAA\tBBB\tCCC\n";
    const char * inputs[] = {
        "5\t110285\t.\tT\tC,<*>\t.\tPASS\tAC=1,0;AD=6,5,0;AF=0.99,0.01;VL_A_STR_INFO=alt_c,alt_nonref;VL_R_STR_INFO=ref,alt_c,alt_nonref\tGT:AD:EC:PL:VL_A_STR_FMT:VL_G_STR_FMT:VL_R_STR_FMT\t.:.:.:.:.:.:.\t0/1:6,5,0:4,0:114,0,15,35,73,113:alt_c,alt_nonref:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_nonref\t.:.:.:.:.:.:.",
        "5\t110290\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
        "5\t110350\t.\tT\t<INS>,<INS>\t.\tPASS\tIMPRECISE;SVLEN=100,200;CIEND=-50,50,-25,25;CIPOS=-10,10,-20,20\tGT\t0/1\t0/1\t0/1",
        "5\t110500\t.\tT\t<CNV>,<CNV>\t.\tPASS\tIMPRECISE;SVLEN=50,100;CILEN=0,25,-25,25;CN=2,4;CICN=-0.5,1,-1.5,1.5\tGT\t0/1\t0/1\t0/1",
        "5\t110700\t.\tA\t<INS:ME>,<INS:ME>\t.\tPASS\tMEINFO=AluY,1,260,+,FLAM_C,1,110,-;METRANS=1,94820,95080,+,1,129678,129788,-\tGT\t0/1\t0/1\t0/1",
        "5\t112000\t.\tC\t<CNV:TR>,<CNV:TR>\t.\tPASS\tRN=2,1;RUS=CAG,TTG,CA;RUL=3,3,2;RB=12,6,6;RUC=4,2,3;RUB=3,3,3,3,3,3,2,2,2;SVLEN=18,6",
        "5\t113000\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
        "5\t114000\t.\tT\tC,A\t.\tPASS\tAC=90,1;AD=6,5,6;AF=0.009,0.0001;VL_A_STR_INFO=alt_c,alt_a;VL_R_STR_INFO=ref,alt_c,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1,2:3,2,0:44,27:114,0,15,35,73,113:alt_c,alt_a:gt_00,gt_01,gt_11,gt_02,gt_12,gt_22:ref,alt_c,alt_a\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
        "5\t115000\t.\tC\t<CNV:TR>,<CNV:TR>\t.\tPASS\tRN=2,1;RUS=CAG,TTG,CA;RUL=3,3,2;RB=12,6,6;RUC=4,2,3;RUB=3,3,3,3,3,3,2,2,2;SVLEN=18,6",
    };
    const char * expected[] = {
        // 2nd allele removed
        "5\t110285\t.\tT\tC\t.\tPASS\tAC=1;AD=6,5;AF=0.99;VL_A_STR_INFO=alt_c;VL_R_STR_INFO=ref,alt_c\tGT:AD:EC:PL:VL_A_STR_FMT:VL_G_STR_FMT:VL_R_STR_FMT\t.:.:.:.:.:.:.\t0/1:6,5:4:114,0,15:alt_c:gt_00,gt_01,gt_11:ref,alt_c\t.:.:.:.:.:.:.",
        "5\t110290\t.\tT\tC\t.\tPASS\tAC=90;AD=6,5;AF=0.009;VL_A_STR_INFO=alt_c;VL_R_STR_INFO=ref,alt_c\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/1:1:3,2:44:114,0,15:alt_c:gt_00,gt_01,gt_11:ref,alt_c\t1/1:1:0,3:46:110,15,0:alt_c:gt_00,gt_01,gt_11:ref,alt_c",
        "5\t110350\t.\tT\t<INS>\t.\tPASS\tIMPRECISE;SVLEN=100;CIEND=-50,50;CIPOS=-10,10\tGT\t0/1\t0/1\t0/1",
        "5\t110500\t.\tT\t<CNV>\t.\tPASS\tIMPRECISE;SVLEN=50;CILEN=0,25;CN=2;CICN=-0.5,1\tGT\t0/1\t0/1\t0/1",
        "5\t110700\t.\tA\t<INS:ME>\t.\tPASS\tMEINFO=AluY,1,260,+;METRANS=1,94820,95080,+\tGT\t0/1\t0/1\t0/1",
        "5\t112000\t.\tC\t<CNV:TR>\t.\tPASS\tRN=2;RUS=CAG,TTG;RUL=3,3;RB=12,6;RUC=4,2;RUB=3,3,3,3,3,3;SVLEN=18",
        // 1st allele removed
        "5\t113000\t.\tT\tA\t.\tPASS\tAC=1;AD=6,6;AF=0.0001;VL_A_STR_INFO=alt_a;VL_R_STR_INFO=ref,alt_a\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/.:1:3,0:27:114,35,113:alt_a:gt_00,gt_02,gt_22:ref,alt_a\t./.:.:0:.:110:.:gt_00:ref",
        // Both alleles removed
        "5\t114000\t.\tT\t.\t.\tPASS\tAD=6;VL_R_STR_INFO=ref\tGT:LAA:LAD:LEC:LPL:VL_LA_STR_FMT:VL_LG_STR_FMT:VL_LR_STR_FMT\t0/0:.:3:.:0:.:gt_00:ref\t0/.:.:3:.:114:.:gt_00:ref\t./.:.:0:.:110:.:gt_00:ref",
        "5\t115000\t.\tC\t.\t.\tPASS\t.",
    };

    kstring_t kstr = KS_INITIALIZE;

    bcf_hdr_t *hdr = bcf_hdr_init("r");
    bcf1_t *rec = bcf_init();
    struct kbitset_t *rm_set = kbs_init(3);
    size_t i;

    if (!hdr)
        error("bcf_hdr_init() failed\n");

    if (!rec)
        error("bcf_init() failed\n");

    if (!rm_set)
        error("kbs_init() failed\n");

    check0(ks_resize(&kstr, 1000));
    check0(bcf_hdr_parse(hdr, header));
    for (i = 0; i < sizeof(inputs)/sizeof(inputs[0]); i++)
    {
        check0(read_vcf_line(inputs[i], hdr, rec, &kstr));
        kbs_clear(rm_set);
        if (rec->pos == 113000 - 1) {
            kbs_insert(rm_set, 1);
        } else if (rec->pos >= 114000 - 1) {
            kbs_insert(rm_set, 1);
            kbs_insert(rm_set, 2);
        } else {
            kbs_insert(rm_set, 2);
        }
        check0(bcf_remove_allele_set(hdr, rec, rm_set));
        check0(vcf_format(hdr, rec, ks_clear(&kstr)));
        chomp(&kstr);
        if (strcmp(expected[i], ks_c_str(&kstr)) != 0)
        {
            error("bcf_remove_allele_set() output differs\nExpected:\n%s\nGot:\n%s\n",
                  expected[i], ks_c_str(&kstr));
        }
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    ks_free(&kstr);
    kbs_destroy(rm_set);
}

int main(int argc, char **argv)
{
    char *fname = argc>1 ? argv[1] : "rmme.bcf";

    // format test. quiet unless there's a failure
    test_get_format_values(fname);

    // main test. writes to stdout
    write_bcf(fname);
    bcf_to_vcf(fname);
    iterator(fname);

    // additional tests. quiet unless there's a failure.
    test_vl_types();
    test_get_info_values(fname);
    test_invalid_end_tag();
    test_open_format();
    test_rlen_values();
    test_bcf_remove_allele_set();
    return 0;
}
