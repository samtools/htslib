/*  test/test-vcf-api.c -- VCF test harness.

    Copyright (C) 2013, 2014, 2017-2019 Genome Research Ltd.

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

#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

void error(const char *format, ...)
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

void write_bcf(char *fname)
{
    // Init
    htsFile *fp    = hts_open(fname,"wb");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    if (!hdr) error("bcf_hdr_init : %s", strerror(errno));
    bcf1_t *rec    = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));

    // Create VCF header
    kstring_t str = {0,0,0};
    check0(bcf_hdr_append(hdr, "##fileDate=20090805"));
    check0(bcf_hdr_append(hdr, "##FORMAT=<ID=UF,Number=1,Type=Integer,Description=\"Unused FORMAT\">"));
    check0(bcf_hdr_append(hdr, "##INFO=<ID=UI,Number=1,Type=Integer,Description=\"Unused INFO\">"));
    check0(bcf_hdr_append(hdr, "##FILTER=<ID=Flt,Description=\"Unused FILTER\">"));
    check0(bcf_hdr_append(hdr, "##unused=<XX=AA,Description=\"Unused generic\">"));
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
    check0(bcf_update_alleles_str(hdr, rec, "G,A"));
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
    bcf_hdr_remove(hdr_out,BCF_HL_STR,"unused");
    bcf_hdr_remove(hdr_out,BCF_HL_GEN,"unused");
    bcf_hdr_remove(hdr_out,BCF_HL_FLT,"Flt");
    bcf_hdr_remove(hdr_out,BCF_HL_INFO,"UI");
    bcf_hdr_remove(hdr_out,BCF_HL_FMT,"UF");
    bcf_hdr_remove(hdr_out,BCF_HL_CTG,"Unused");
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
    test_get_info_values(fname);

    return 0;
}
