#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>

void write_bcf(char *fname)
{
    // Init
    htsFile *fp    = hts_open(fname,"wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf1_t *rec    = bcf_init1();

    // Create VCF header
    kstring_t str = {0,0,0};
    kputs("##contig=<ID=20,length=63025520>", &str);
    bcf_hdr_append(hdr, str.s);
    bcf_hdr_append(hdr, "##fileDate=20090805");
    bcf_hdr_append(hdr, "##source=myImputationProgramV3.1");
    bcf_hdr_append(hdr, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    bcf_hdr_append(hdr, "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    bcf_hdr_append(hdr, "##phasing=partial");
    bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    bcf_hdr_add_sample(hdr, "NA00001");
    bcf_hdr_add_sample(hdr, "NA00002");
    bcf_hdr_add_sample(hdr, "NA00003");

    bcf_hdr_fmt_text(hdr);
    bcf_hdr_write(fp, hdr);


    // Add a record
    // 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    // .. CHROM
    rec->rid = bcf_hdr_name2id(hdr, "20");
    // .. POS
    rec->pos = 14369;
    // .. ID
    bcf_update_id(hdr, rec, "rs6054257");
    // .. REF and ALT
    bcf_update_alleles_str(hdr, rec, "G,A");
    // .. QUAL
    rec->qual = 29;
    // .. FILTER
    int tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = 3;  
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 14; 
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    float tmpf = 0.5;
    bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
    // .. FORMAT
    int *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_phased(0); 
    tmpia[1] = bcf_gt_phased(0);
    tmpia[2] = bcf_gt_phased(1); 
    tmpia[3] = bcf_gt_phased(0);
    tmpia[4] = bcf_gt_unphased(1); 
    tmpia[5] = bcf_gt_unphased(1);
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
    bcf_write1(fp, hdr, rec);

    // 20     1110696 . A      G,T     67   .   NS=2;DP=10;AF=0.333,.;AA=T;DB GT 2 1   ./.
    bcf_clear1(rec);
    rec->rid = bcf_hdr_name2id(hdr, "20");
    rec->pos = 1110695;
    bcf_update_alleles_str(hdr, rec, "A,G,T");
    rec->qual = 67;
    tmpi = 2;
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 10;
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    float *tmpfa = (float*)malloc(2*sizeof(float));
    tmpfa[0] = 0.333; 
    bcf_float_set_missing(tmpfa[1]);
    bcf_update_info_float(hdr, rec, "AF", tmpfa, 2);
    bcf_update_info_string(hdr, rec, "AA", "T");
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    tmpia[0] = bcf_gt_phased(2);
    tmpia[1] = bcf_int32_vector_end;
    tmpia[2] = bcf_gt_phased(1);    
    tmpia[3] = bcf_int32_vector_end;
    tmpia[4] = bcf_gt_missing;
    tmpia[5] = bcf_gt_missing;
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    bcf_write1(fp, hdr, rec);

    free(tmpia);
    free(tmpfa);

    // Clean
    free(str.s);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

void bcf_to_vcf(char *fname)
{
    htsFile *fp    = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init1();
    htsFile *out   = hts_open("-","w");

    bcf_hdr_write(out, hdr);
    while ( bcf_read1(fp, hdr, rec)>=0 )
    {
        bcf_write1(out, hdr, rec);
    }
    
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    hts_close(out);
}

int main(int argc, char **argv)
{
    char *fname = argc>1 ? argv[1] : "rmme.bcf";
    write_bcf(fname);
    bcf_to_vcf(fname);
    return 0;
}

