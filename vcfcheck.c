#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
    int min, max, step, m_vals;
    uint64_t *vals;
}
idist_t;

typedef struct
{
    int n_snps, n_indels, n_mnps, n_others, n_mals;
    int *af_ts, *af_tv, *af_snps, *af_indels;
    #if QUAL_STATS
        int *qual_ts, *qual_tv, *qual_snps, *qual_indels;
    #endif
    int *insertions, *deletions, m_indel;   // maximum indel length
    int in_frame, out_frame;
    int subst[15];
    int *smpl_hets, *smpl_homRR, *smpl_homAA, *smpl_ts, *smpl_tv, *smpl_indels, *smpl_ndp;
    unsigned long int *smpl_dp;
    idist_t dp;
}
stats_t;

typedef struct
{
    int m[3], mm[3];        // number of hom, het and non-ref hom matches and mismatches
    float r2sum;
    int r2n;
}
gtcmp_t;

typedef struct
{
    // stats
    stats_t stats[3];
    int *tmp_iaf, ntmp_iaf, m_af, m_qual;
    int dp_min, dp_max, dp_step;
    gtcmp_t *af_gts_snps, *af_gts_indels, *smpl_gts_snps, *smpl_gts_indels;

    // other
    readers_t *files;
    regions_t regions;
    int prev_reg;
    char **argv, *exons_file, *samples_file;
    int argc, debug;
    int split_by_id, nstats;
}
args_t;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

inline int acgt2int(char c)
{
    if ( (int)c>96 ) c -= 32;
    if ( c=='A' ) return 0;
    if ( c=='C' ) return 1;
    if ( c=='G' ) return 2;
    if ( c=='T' ) return 3;
    return -1;
}
#define int2acgt(i) "ACGT"[i]

void idist_init(idist_t *d, int min, int max, int step)
{
    d->min = min; d->max = max; d->step = step;
    d->m_vals = 4 + (d->max - d->min)/d->step;
    d->vals = (uint64_t*) calloc(d->m_vals,sizeof(uint64_t));
}
void idist_destroy(idist_t *d)
{
    if ( d->vals ) free(d->vals);
}
inline uint64_t *idist(idist_t *d, int val)
{
    if ( val < d->min ) return &d->vals[0];
    if ( val > d->max ) return &d->vals[d->m_vals-1];
    return &d->vals[1 + (val - d->min) / d->step];
}
inline int idist_i2bin(idist_t *d, int i)
{
    if ( i<=0 ) return d->min;
    if ( i>= d->m_vals ) return d->max;
    return i-1+d->min;
}

void init_stats(args_t *args)
{
    int i;
    args->nstats = args->files->nreaders==1 ? 1 : 3;
    if ( args->split_by_id ) args->nstats = 2;

    // AF corresponds to AC but is more robust for mixture of haploid and diploid GTs
    args->m_af = 101;
    for (i=0; i<args->files->nreaders; i++)
        if ( args->files->readers[i].header->n[BCF_DT_SAMPLE] + 1> args->m_af )
            args->m_af = args->files->readers[i].header->n[BCF_DT_SAMPLE] + 1;

    #if QUAL_STATS
        args->m_qual = 500;
    #endif

    if ( args->samples_file )
    {
        if ( !bcf_sr_set_samples(args->files,args->samples_file) ) error("Could not initialize samples: %s\n", args->samples_file);
        args->af_gts_snps     = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
        args->af_gts_indels   = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
        args->smpl_gts_snps   = (gtcmp_t *) calloc(args->files->n_smpl,sizeof(gtcmp_t));
        args->smpl_gts_indels = (gtcmp_t *) calloc(args->files->n_smpl,sizeof(gtcmp_t));
    }
    for (i=0; i<args->nstats; i++)
    {
        stats_t *stats = &args->stats[i];
        stats->m_indel     = 60;
        stats->insertions  = (int*) calloc(stats->m_indel,sizeof(int));
        stats->deletions   = (int*) calloc(stats->m_indel,sizeof(int));
        stats->af_ts       = (int*) calloc(args->m_af,sizeof(int));
        stats->af_tv       = (int*) calloc(args->m_af,sizeof(int));
        stats->af_snps     = (int*) calloc(args->m_af,sizeof(int));
        stats->af_indels   = (int*) calloc(args->m_af,sizeof(int));
        #if QUAL_STATS
            stats->qual_ts     = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_tv     = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_snps   = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_indels = (int*) calloc(args->m_qual,sizeof(int));
        #endif
        if ( args->files->n_smpl )
        {
            stats->smpl_hets   = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homAA  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homRR  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_ts     = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_tv     = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_indels = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_dp     = (unsigned long int *) calloc(args->files->n_smpl,sizeof(unsigned long int));
            stats->smpl_ndp    = (int *) calloc(args->files->n_smpl,sizeof(int));
        }
        idist_init(&stats->dp, args->dp_min,args->dp_max,args->dp_step);
    }

    if ( args->exons_file )
    {
        if ( !init_regions(args->exons_file, &args->regions) )
            error("Error occurred while reading, was the file compressed with bgzip: %s?\n", args->exons_file);
        args->prev_reg = -1;
    }
}
void destroy_stats(args_t *args)
{
    int id;
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        if (stats->af_ts) free(stats->af_ts);
        if (stats->af_tv) free(stats->af_tv);
        if (stats->af_snps) free(stats->af_snps);
        if (stats->af_indels) free(stats->af_indels);
        #if QUAL_STATS
            if (stats->qual_ts) free(stats->qual_ts);
            if (stats->qual_tv) free(stats->qual_tv);
            if (stats->qual_snps) free(stats->qual_snps);
            if (stats->qual_indels) free(stats->qual_indels);
        #endif
        free(stats->insertions);
        free(stats->deletions);
        if (stats->smpl_hets) free(stats->smpl_hets);
        if (stats->smpl_homAA) free(stats->smpl_homAA);
        if (stats->smpl_homRR) free(stats->smpl_homRR);
        if (stats->smpl_ts) free(stats->smpl_ts);
        if (stats->smpl_tv) free(stats->smpl_tv);
        if (stats->smpl_indels) free(stats->smpl_indels);
        if (stats->smpl_dp) free(stats->smpl_dp);
        if (stats->smpl_ndp) free(stats->smpl_ndp);
        idist_destroy(&stats->dp);
    }
    if (args->tmp_iaf) free(args->tmp_iaf);
    if (args->exons_file) destroy_regions(&args->regions);
    if (args->af_gts_snps) free(args->af_gts_snps);
    if (args->af_gts_indels) free(args->af_gts_indels);
    if (args->smpl_gts_snps) free(args->smpl_gts_snps);
    if (args->smpl_gts_indels) free(args->smpl_gts_indels);
}

void init_iaf(args_t *args, reader_t *reader)
{
    bcf1_t *line = reader->buffer[0];
    if ( args->ntmp_iaf < line->n_allele )
    {
        args->tmp_iaf = (int*)realloc(args->tmp_iaf, line->n_allele*sizeof(int));
        args->ntmp_iaf = line->n_allele;
    }
    int ret = calc_ac(reader->header, line, args->tmp_iaf, args->samples_file ? BCF_UN_INFO|BCF_UN_FMT : BCF_UN_INFO);
    if ( ret )
    {
        int i, an=0;
        for (i=0; i<line->n_allele; i++)
            an += args->tmp_iaf[i];
        
        for (i=1; i<line->n_allele; i++)
        {
            if ( args->tmp_iaf[i]==1 ) 
                args->tmp_iaf[i] = 0; // singletons into the first bin
            else if ( !an )
                args->tmp_iaf[i] = 1;   // no genotype at all, put to the AF=0 bin
            else
                args->tmp_iaf[i] = 1 + args->tmp_iaf[i] * (args->m_af-2.0) / an;
        }
    }
            
    // todo: otherwise use AF 
}

inline void do_mnp_stats(args_t *args, stats_t *stats, reader_t *reader)
{
    stats->n_mnps++;
}

inline void do_other_stats(args_t *args, stats_t *stats, reader_t *reader)
{
    stats->n_others++;
}

void do_indel_stats(args_t *args, stats_t *stats, reader_t *reader)
{
    stats->n_indels++;

    bcf1_t *line = reader->buffer[0];

    #if QUAL_STATS
        int iqual = line->qual >= args->m_qual || isnan(line->qual) ? args->m_qual - 1 : line->qual;
        stats->qual_indels[iqual]++;
    #endif

    // Check if the indel is near an exon for the frameshift statistics
    pos_t *reg = NULL, *reg_next = NULL;
    if ( args->regions.nseqs )
    {
        if ( args->files->iseq!=args->prev_reg )
        {
            reset_regions(&args->regions, args->files->seqs[args->files->iseq]);
            args->prev_reg = args->files->iseq;
        }
        reg = is_in_regions(&args->regions, line->pos+1);
        if ( !reg && args->regions.cseq>=0 && args->regions.cpos < args->regions.npos[args->regions.cseq] )
            reg_next = &args->regions.pos[args->regions.cseq][args->regions.cpos];
    }

    int i;
    for (i=1; i<line->n_allele; i++)
    {
        if ( line->d.var[i].type!=VCF_INDEL ) continue;
        stats->af_indels[ args->tmp_iaf[i] ]++;
        int len = line->d.var[i].n;

        // Check the frameshifts
        int tlen = 0;
        if ( reg )
        {
            tlen = abs(len);
            if ( len<0 ) 
            {
                int to = line->pos+1 + tlen;
                if ( to > reg->to ) tlen -= to-reg->to;
            }
        }
        else if ( reg_next && len<0 )
        {
            tlen = abs(len) - reg_next->from + line->pos+1;
            if ( tlen<0 ) tlen = 0;
        }
        if ( tlen )
        {
            if ( tlen%3 ) stats->out_frame++;
            else stats->in_frame++;

            //if ( tlen%3 ) printf("%s\t%d\t%d\t%d\tframeshift (tlen=%d, next=%d)\n", args->files->seqs[args->files->iseq],line->pos+1,reg->from,reg->to,tlen,reg_next);
            //else printf("%s\t%d\t%d\t%d\tin-frame\n", args->files->seqs[args->files->iseq],line->pos+1,reg->from,reg->to);
        }

        // Indel length distribution
        int *ptr = stats->insertions;
        if ( len<0 ) 
        {
            len *= -1;
            ptr = stats->deletions;
        }
        if ( --len >= stats->m_indel ) len = stats->m_indel-1;
        ptr[len]++;
    }
}

void do_snp_stats(args_t *args, stats_t *stats, reader_t *reader)
{
    stats->n_snps++;

    bcf1_t *line = reader->buffer[0];
    int ref = acgt2int(*line->d.allele[0]);
    if ( ref<0 ) return;

    #if QUAL_STATS
        int iqual = line->qual >= args->m_qual || isnan(line->qual) ? args->m_qual - 1 : line->qual;
        stats->qual_snps[iqual]++;
    #endif

    int i;
    for (i=1; i<line->n_allele; i++)
    {
        if ( !(line->d.var[i].type&VCF_SNP) ) continue;
        int alt = acgt2int(*line->d.allele[i]);
        if ( alt<0 ) continue;
        stats->subst[ref<<2|alt]++;
        int iaf = args->tmp_iaf[i];
        stats->af_snps[iaf]++;
        if ( abs(ref-alt)==2 ) 
        {
            stats->af_ts[iaf]++;
            #if QUAL_STATS
                stats->qual_ts[iqual]++;
            #endif
        }
        else 
        {
            stats->af_tv[iaf]++;
            #if QUAL_STATS
                stats->qual_tv[iqual]++;
            #endif
        }
    }
}

void do_sample_stats(args_t *args, stats_t *stats, reader_t *reader, int matched)
{
    readers_t *files = args->files;
    bcf1_t *line = reader->buffer[0];
    bcf_fmt_t *fmt_ptr;

    if ( (fmt_ptr = get_fmt_ptr(reader->header,reader->buffer[0],"GT")) )
    {
        int ref = acgt2int(*line->d.allele[0]);
        int is;
        for (is=0; is<args->files->n_smpl; is++)
        {
            int ial;
            int gt = gt_type(fmt_ptr, reader->samples[is], &ial);
            if ( gt == GT_UNKN ) continue;
            if ( line->d.var_type&VCF_SNP )
            {
                if ( gt == GT_HET_RA || gt == GT_HET_AA ) stats->smpl_hets[is]++;
                else if ( gt == GT_HOM_RR ) stats->smpl_homRR[is]++;
                else if ( gt == GT_HOM_AA ) stats->smpl_homAA[is]++;
                if ( gt != GT_HOM_RR && line->d.var[ial].type&VCF_SNP )
                {
                    int alt = acgt2int(*line->d.allele[ial]);
                    if ( alt<0 ) continue;
                    if ( abs(ref-alt)==2 ) 
                        stats->smpl_ts[is]++;
                    else
                        stats->smpl_tv[is]++;
                }
            }
            if ( line->d.var_type&VCF_INDEL )
            {
                if ( gt != GT_HOM_RR ) stats->smpl_indels[is]++;
            }
        }
    }

    if ( (fmt_ptr = get_fmt_ptr(reader->header,reader->buffer[0],"DP")) )
    {
        #define BRANCH_INT(type_t,missing) { \
            type_t *p = (type_t *) fmt_ptr->p; \
            for (is=0; is<args->files->n_smpl; is++) \
                if (p[is]!=missing) { \
                    (*idist(&stats->dp, p[is]))++; \
                    stats->smpl_ndp[is]++; \
                    stats->smpl_dp[is] += p[is]; \
                } \
            }

        int is;
        if ( fmt_ptr->type==BCF_BT_INT8 ) { BRANCH_INT(uint8_t, 0x80) }
        else if ( fmt_ptr->type==BCF_BT_INT16 ) { BRANCH_INT(uint16_t, 0x8000) } 
        else if ( fmt_ptr->type==BCF_BT_INT32 ) { BRANCH_INT(uint32_t, 0x80000000) }

        #undef BRANCH_INT
    }
   
    if ( matched==3 )
    {
        int is;
        bcf_fmt_t *fmt0, *fmt1;
        fmt0 = get_fmt_ptr(files->readers[0].header,files->readers[0].buffer[0],"GT"); if ( !fmt0 ) return;
        fmt1 = get_fmt_ptr(files->readers[1].header,files->readers[1].buffer[0],"GT"); if ( !fmt1 ) return;

        int iaf = args->tmp_iaf[1]; // only first ALT alelle considered
        gtcmp_t *af_stats = files->readers[0].buffer[0]->d.var_type&VCF_SNP ? args->af_gts_snps : args->af_gts_indels;
        gtcmp_t *smpl_stats = files->readers[0].buffer[0]->d.var_type&VCF_SNP ? args->smpl_gts_snps : args->smpl_gts_indels;

        int r2n = 0;
        float r2a2 = 0, r2a = 0, r2b2 = 0, r2b = 0, r2ab = 0;
        for (is=0; is<files->n_smpl; is++)
        {
            // Simplified comparison: only 0/0, 0/1, 1/1 is looked at as the identity of 
            //  actual alleles can be enforced by running without the -c option.
            int gt = gt_type(fmt0, files->readers[0].samples[is], NULL);
            if ( gt == GT_UNKN ) continue;

            int match = 1;
            int gt2 = gt_type(fmt1, files->readers[1].samples[is], NULL);
            if ( gt2 == GT_UNKN ) match = -1;
            else if ( gt != gt2 ) match = 0;

            if ( match == -1 ) continue;
            if ( gt == GT_HET_AA ) gt = GT_HOM_AA;  // rare case, treat as AA hom
            if ( gt2 == GT_HET_AA ) gt2 = GT_HOM_AA;
            if ( match ) 
            {
                af_stats[iaf].m[gt]++;
                smpl_stats[is].m[gt]++;
            }
            else 
            {
                af_stats[iaf].mm[gt]++;
                smpl_stats[is].mm[gt]++;
            }
            r2a2 += gt*gt;
            r2a  += gt;
            r2b2 += gt2*gt2;
            r2b  += gt2;
            r2ab += gt*gt2;
            r2n++;
        }
        if ( r2n )
        {
            float cov  = r2ab - r2a*r2b/r2n;
            float var2 = (r2a2 - r2a*r2a/r2n) * (r2b2 - r2b*r2b/r2n);
            af_stats[iaf].r2sum += var2==0 ? 1 : cov*cov/var2;
            af_stats[iaf].r2n++;
        }

        if ( args->debug )
        {
            for (is=0; is<files->n_smpl; is++)
            {
                int gt = gt_type(fmt0, files->readers[0].samples[is], NULL);
                if ( gt == GT_UNKN ) continue;
                int gt2 = gt_type(fmt1, files->readers[1].samples[is], NULL);
                if ( gt != gt2 ) 
                {
                    fprintf(stderr,"%s\t%d\t%s\t%d\t%d\n",args->files->seqs[args->files->iseq],files->readers[0].buffer[0]->pos+1,files->samples[is],gt,gt2);
                }
            }
        }
    }
}

void check_vcf(args_t *args)
{
    int ret,i;
    readers_t *files = args->files;
    while ( (ret=bcf_sr_next_line(files)) )
    {
        reader_t *reader = NULL;
        bcf1_t *line = NULL;
        for (i=0; i<files->nreaders; i++)
        {
            if ( !(ret&1<<i) ) continue;
            reader = &files->readers[i];
            line = files->readers[i].buffer[0];
            break;
        }
        set_variant_types(line);
        init_iaf(args, reader);

        stats_t *stats = &args->stats[ret-1];
        if ( args->split_by_id && line->d.id[0]=='.' && !line->d.id[1] )
            stats = &args->stats[1];

        if ( line->d.var_type&VCF_SNP ) 
            do_snp_stats(args, stats, reader);
        if ( line->d.var_type&VCF_INDEL )
            do_indel_stats(args, stats, reader);
        if ( line->d.var_type&VCF_MNP )
            do_mnp_stats(args, stats, reader);
        if ( line->d.var_type&VCF_OTHER )
            do_other_stats(args, stats, reader);

        if ( line->n_allele>2 ) stats->n_mals++;

        if ( files->n_smpl )
            do_sample_stats(args, stats, reader, ret);
    }
}

void print_header(args_t *args)
{
    int i;
    printf("# This file was produced by vcfcheck and can be plotted using plot-vcfcheck.\n");
    printf("# The command line was:\thtscmd %s ", args->argv[0]);
    for (i=1; i<args->argc; i++)
        printf(" %s",args->argv[i]);
    printf("\n#\n");

    printf("# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names\n");
    if ( args->files->nreaders==1 )
    {
        if ( args->split_by_id )
        {
            printf("ID\t0\t%s:known (sites with ID different from \".\")\n", args->files->readers[0].fname);
            printf("ID\t1\t%s:novel (sites where ID column is \".\")\n", args->files->readers[0].fname);
        }
        else
            printf("ID\t0\t%s\n", args->files->readers[0].fname);
    }
    else
    {
        printf("ID\t0\t%s\n", args->files->readers[0].fname);
        printf("ID\t1\t%s\n", args->files->readers[1].fname);
        printf("ID\t2\t%s\t%s\n", args->files->readers[0].fname,args->files->readers[1].fname);
    }
}

void print_stats(args_t *args)
{
    int i, id;
    printf("# SN, Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
    for (id=0; id<args->files->nreaders; id++)
        printf("SN\t%d\tnumber of samples:\t%d\n", id, args->files->readers[id].header->n[BCF_DT_SAMPLE]);
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        printf("SN\t%d\tnumber of SNPs:\t%d\n", id, stats->n_snps);
        printf("SN\t%d\tnumber of MNPs:\t%d\n", id, stats->n_mnps);
        printf("SN\t%d\tnumber of indels:\t%d\n", id, stats->n_indels);
        printf("SN\t%d\tnumber of others:\t%d\n", id, stats->n_others);
        printf("SN\t%d\tnumber of multiallelic sites:\t%d\n", id, stats->n_mals);

        int ts=0,tv=0;
        for (i=0; i<args->m_af; i++) { ts += stats->af_ts[i]; tv += stats->af_tv[i];  }
        printf("SN\t%d\tts/tv:\t%.2f\n", id, tv?(float)ts/tv:0);
    }
    if ( args->exons_file )
    {
        printf("# FS, Indel frameshifts:\n# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]out/(in+out) ratio\n");
        for (id=0; id<args->nstats; id++)
        {
            int in=args->stats[id].in_frame, out=args->stats[id].out_frame;
            printf("FS\t%d\t%d\t%d\t%.2f\n", id, in,out,out?(float)out/(in+out):0);
        }
    }
    printf("# Sis, Singleton stats:\n# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        printf("SiS\t%d\t%d\t%d\t%d\t%d\t%d\n", id,1,stats->af_snps[0],stats->af_ts[0],stats->af_tv[0],stats->af_indels[0]);
        stats->af_snps[1]   += stats->af_snps[0];
        stats->af_ts[1]     += stats->af_ts[0];
        stats->af_tv[1]     += stats->af_tv[0];
        stats->af_indels[1] += stats->af_indels[0];
    }
    printf("# AF, Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=1; i<args->m_af; i++)
        {
            if ( stats->af_snps[i]+stats->af_ts[i]+stats->af_tv[i]+stats->af_indels[i] == 0  ) continue;
            printf("AF\t%d\t%f\t%d\t%d\t%d\t%d\n", id,100.*(i-1)/(args->m_af-2),stats->af_snps[i],stats->af_ts[i],stats->af_tv[i],stats->af_indels[i]);
        }
    }
    #if QUAL_STATS
        printf("# QUAL, Stats by quality:\n# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->m_qual; i++)
            {
                if ( stats->qual_snps[i]+stats->qual_ts[i]+stats->qual_tv[i]+stats->qual_indels[i] == 0  ) continue;
                printf("QUAL\t%d\t%d\t%d\t%d\t%d\t%d\n", id,i,stats->qual_snps[i],stats->qual_ts[i],stats->qual_tv[i],stats->qual_indels[i]);
            }
        }
    #endif
    printf("# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]count\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=stats->m_indel-1; i>=0; i--)
            if ( stats->deletions[i] ) printf("IDD\t%d\t%d\t%d\n", id,-i-1,stats->deletions[i]);
        for (i=0; i<stats->m_indel; i++)
            if ( stats->insertions[i] ) printf("IDD\t%d\t%d\t%d\n", id,i+1,stats->insertions[i]);
    }
    printf("# ST, Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
    for (id=0; id<args->nstats; id++)
    {
        int t;
        for (t=0; t<15; t++)
        {
            if ( t>>2 == (t&3) ) continue;
            printf("ST\t%d\t%c>%c\t%d\n", id, int2acgt(t>>2),int2acgt(t&3),args->stats[id].subst[t]);
        }
    }
    if ( args->files->nreaders>1 && args->files->n_smpl )
    {
        printf("SN\t%d\tNumber of samples:\t%d\n", 2, args->files->n_smpl);

        printf("# GCsAF, Genotype concordance by non-reference allele frequency (SNPs)\n# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of sites\n");
        gtcmp_t *stats = args->af_gts_snps;
        int nrd_m[3] = {0,0,0}, nrd_mm[3] = {0,0,0};
        for (i=1; i<args->m_af; i++)
        {
            int j, n = 0;
            for (j=0; j<3; j++) 
            {
                n += stats[i].m[j] + stats[i].mm[j];
                nrd_m[j]  += stats[i].m[j];
                nrd_mm[j] += stats[i].mm[j];
            }
            if ( !n ) continue;
            printf("GCsAF\t2\t%f", 100.*(i-1)/(args->m_af-2));
            printf("\t%d\t%d\t%d", stats[i].m[GT_HOM_RR],stats[i].m[GT_HET_RA],stats[i].m[GT_HOM_AA]);
            printf("\t%d\t%d\t%d", stats[i].mm[GT_HOM_RR],stats[i].mm[GT_HET_RA],stats[i].mm[GT_HOM_AA]);
            printf("\t%f\t%d\n", stats[i].r2n ? stats[i].r2sum/stats[i].r2n : -1.0, stats[i].r2n);
        }

        int m  = nrd_m[GT_HET_RA] + nrd_m[GT_HOM_AA];
        int mm = nrd_mm[GT_HOM_RR] + nrd_mm[GT_HET_RA] + nrd_mm[GT_HOM_AA];
        printf("# Non-Reference Discordance (NRD)\n# NRD\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n");
        printf("NRD\t2\t%f\t%f\t%f\t%f\n", 
            m+mm ? mm*100.0/(m+mm) : 0, 
            nrd_m[GT_HOM_RR]+nrd_mm[GT_HOM_RR] ? nrd_mm[GT_HOM_RR]*100.0/(nrd_m[GT_HOM_RR]+nrd_mm[GT_HOM_RR]) : 0,
            nrd_m[GT_HET_RA]+nrd_mm[GT_HET_RA] ? nrd_mm[GT_HET_RA]*100.0/(nrd_m[GT_HET_RA]+nrd_mm[GT_HET_RA]) : 0,
            nrd_m[GT_HOM_AA]+nrd_mm[GT_HOM_AA] ? nrd_mm[GT_HOM_AA]*100.0/(nrd_m[GT_HOM_AA]+nrd_mm[GT_HOM_AA]) : 0
            );

        printf("# GCcS, Genotype concordance by sample (SNPs)\n# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\n");
        stats = args->smpl_gts_snps;
        for (i=0; i<args->files->n_smpl; i++)
        {
            int m  = stats[i].m[GT_HET_RA] + stats[i].m[GT_HOM_AA];
            int mm = stats[i].mm[GT_HOM_RR] + stats[i].mm[GT_HET_RA] + stats[i].mm[GT_HOM_AA];
            printf("GCsS\t2\t%s\t%.3f", args->files->samples[i], m+mm ? mm*100.0/(m+mm) : 0);
            printf("\t%d\t%d\t%d", stats[i].m[GT_HOM_RR],stats[i].m[GT_HET_RA],stats[i].m[GT_HOM_AA]);
            printf("\t%d\t%d\t%d\n", stats[i].mm[GT_HOM_RR],stats[i].mm[GT_HET_RA],stats[i].mm[GT_HOM_AA]);
        }
    }

    if ( args->files->n_smpl )
    {
        printf("# PSC, Per-sample counts\n# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->files->n_smpl; i++)
            {
                float dp = stats->smpl_ndp[i] ? stats->smpl_dp[i]/stats->smpl_ndp[i] : 0;
                printf("PSC\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.0f\n", id,args->files->samples[i], 
                    stats->smpl_homRR[i], stats->smpl_homAA[i], stats->smpl_hets[i], stats->smpl_ts[i], stats->smpl_tv[i], stats->smpl_indels[i],dp);
            }
        }

        printf("# DP, Depth distribution\n# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%%)\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            long unsigned int sum = 0;
            for (i=0; i<stats->dp.m_vals; i++) { sum += stats->dp.vals[i]; }
            for (i=0; i<stats->dp.m_vals; i++)
            {
                if ( stats->dp.vals[i]==0 ) continue;
                printf("DP\t%d\t", id);
                if ( i==0 ) printf("<%d", stats->dp.min);
                else if ( i+1==stats->dp.m_vals ) printf(">%d", stats->dp.max);
                else printf("%d", idist_i2bin(&stats->dp,i));
                printf("\t%ld\t%f\n", stats->dp.vals[i], stats->dp.vals[i]*100./sum);
            }
        }
    }
}

static void usage(void)
{
    fprintf(stderr, "\nAbout:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfcheck.\n");
    fprintf(stderr, "         When two files are given, the program generates separate stats for intersection\n");
    fprintf(stderr, "         and the complements.\n");
    fprintf(stderr, "Usage:   vcfcheck [options] <A.vcf.gz> [<B.vcf.gz>]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --collapse <string>           treat sites with differing alleles as same for <snps|indels|both|any>\n");
    fprintf(stderr, "    -d, --depth <int,int,int>         depth distribution: min,max,bin size [0,500,1]\n");
    fprintf(stderr, "        --debug                       produce verbose per-site and per-sample output\n");
    fprintf(stderr, "    -e, --exons <file.gz>             tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed)\n");
    fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
    fprintf(stderr, "    -i, --split-by-ID                 collect stats for sites with ID separately (known vs novel)\n");
    fprintf(stderr, "    -r, --region <chr|chr:from-to>    collect stats in the given region only\n");
    fprintf(stderr, "    -s, --samples <list|file>         produce sample stats, \"-\" to include all samples\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfcheck(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->dp_min = 0; args->dp_max = 500; args->dp_step = 1;

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"collapse",1,0,'c'},
        {"debug",0,0,1},
        {"depth",1,0,'d'},
        {"apply-filters",0,0,'f'},
        {"exons",1,0,'e'},
        {"samples",1,0,'s'},
        {"split-by-ID",0,0,'i'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hc:fr:e:s:d:i1",loptions,NULL)) >= 0) {
        switch (c) {
            case 'c':
                if ( !strcmp(optarg,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) args->files->collapse |= COLLAPSE_ANY;
                break;
            case  1 : args->debug = 1; break;
            case 'd': 
                if ( sscanf(optarg,"%d,%d,%d",&args->dp_min,&args->dp_max,&args->dp_step)!=3 )
                    error("Could not parse --depth %s\n", optarg); 
                if ( args->dp_min<0 || args->dp_min >= args->dp_max || args->dp_step > args->dp_max - args->dp_min + 1 )
                    error("Is this a typo? --depth %s\n", optarg);
                break;
            case 'f': args->files->apply_filters = 1; break;
            case 'r': args->files->region = optarg; break;
            case 'e': args->exons_file = optarg; break;
            case 's': args->samples_file = optarg; break;
            case 'i': args->split_by_id = 1; break;
            case 'h': 
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if (argc == optind) usage();

    if ( argc-optind>2 ) usage();
    if ( args->split_by_id && argc-optind>1 ) error("Only one file can be given with -i.\n");
    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Could not read the file: %s\n", argv[optind]);
        optind++;
    }

    init_stats(args);
    print_header(args);
    check_vcf(args);
    print_stats(args);
    destroy_stats(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

