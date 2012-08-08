#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
	int n_snps, n_indels, n_mals;
	int *af_ts, *af_tv, *af_snps, *af_indels;
	int *insertions, *deletions, m_indel;	// maximum indel length
	int in_frame, out_frame;
	int subst[15];
	int *smpl_hets, *smpl_ts, *smpl_tv;
}
stats_t;

typedef struct
{
	int m[3], mm[3];		// number of hom, het and non-ref hom matches and mismatches
}
gtcmp_t;

typedef struct
{
	// stats
	stats_t stats[3];
	int *tmp_iaf, ntmp_iaf, m_af;
	gtcmp_t *af_gts_snps, *af_gts_indels, *smpl_gts_snps, *smpl_gts_indels;

	// other
	readers_t files;
	regions_t regions;
	int prev_reg;
	char **argv, *exons_file, *samples_file;
	int argc;
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

void init_stats(args_t *args)
{
	int i,nstats = args->files.nreaders==1 ? 1 : 3;

	// AF corresponds to AC but is more robust for mixture of haploid and diploid GTs
	args->m_af = 101;
	for (i=0; i<args->files.nreaders; i++)
		if ( args->files.readers[i].header->n[BCF_DT_SAMPLE] + 1> args->m_af )
			args->m_af = args->files.readers[i].header->n[BCF_DT_SAMPLE] + 1;

	if ( args->samples_file )
	{
		if ( !init_samples(args->samples_file,&args->files) ) error("Could not initialize samples: %s\n", args->samples_file);
		args->af_gts_snps     = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
		args->af_gts_indels   = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
		args->smpl_gts_snps   = (gtcmp_t *) calloc(args->files.n_smpl,sizeof(gtcmp_t));
		args->smpl_gts_indels = (gtcmp_t *) calloc(args->files.n_smpl,sizeof(gtcmp_t));
	}
	for (i=0; i<nstats; i++)
	{
		stats_t *stats = &args->stats[i];
		stats->m_indel    = 60;
		stats->insertions = (int*) calloc(stats->m_indel,sizeof(int));
		stats->deletions  = (int*) calloc(stats->m_indel,sizeof(int));
		stats->af_ts      = (int*)calloc(args->m_af,sizeof(int));
		stats->af_tv      = (int*)calloc(args->m_af,sizeof(int));
		stats->af_snps    = (int*)calloc(args->m_af,sizeof(int));
		stats->af_indels  = (int*)calloc(args->m_af,sizeof(int));
		if ( args->files.n_smpl )
		{
			stats->smpl_hets = (int *) calloc(args->files.n_smpl,sizeof(int));
			stats->smpl_ts   = (int *) calloc(args->files.n_smpl,sizeof(int));
			stats->smpl_tv   = (int *) calloc(args->files.n_smpl,sizeof(int));
		}
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
	int id,nstats = args->files.nreaders==1 ? 1 : 3;
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		if (stats->af_ts) free(stats->af_ts);
		if (stats->af_tv) free(stats->af_tv);
		if (stats->af_snps) free(stats->af_snps);
		if (stats->af_indels) free(stats->af_indels);
		free(stats->insertions);
		free(stats->deletions);
		if (stats->smpl_hets) free(stats->smpl_hets);
		if (stats->smpl_ts) free(stats->smpl_ts);
		if (stats->smpl_tv) free(stats->smpl_tv);
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
	bcf1_t *line = reader->line;
	if ( args->ntmp_iaf < line->n_allele )
	{
		args->tmp_iaf = (int*)realloc(args->tmp_iaf, line->n_allele*sizeof(int));
		args->ntmp_iaf = line->n_allele;
	}
	int ret = calc_ac(reader->header, line, args->tmp_iaf, BCF_UN_FMT);
	if ( ret )
	{
		int i, an=0;
		for (i=0; i<line->n_allele; i++)
			an += args->tmp_iaf[i];
		
		for (i=1; i<line->n_allele; i++)
		{
			if ( args->tmp_iaf[i]==1 ) 
				args->tmp_iaf[i]=0; // singletons into the first bin
			else 
				args->tmp_iaf[i] = 1 + args->tmp_iaf[i] * (args->m_af-2.0) / an;
		}
	}
	// todo: otherwise use AF 
}

void do_indel_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_indels++;

	bcf1_t *line = reader->line;

	// Check if the indel is near an exon for the frameshift statistics
	pos_t *reg = NULL, *reg_next = NULL;
	if ( args->regions.nseqs )
	{
		if ( args->files.iseq!=args->prev_reg )
		{
			reset_regions(&args->regions, args->files.seqs[args->files.iseq]);
			args->prev_reg = args->files.iseq;
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

			//if ( tlen%3 ) printf("%s\t%d\t%d\t%d\tframeshift (tlen=%d, next=%d)\n", args->files.seqs[args->files.iseq],line->pos+1,reg->from,reg->to,tlen,reg_next);
			//else printf("%s\t%d\t%d\t%d\tin-frame\n", args->files.seqs[args->files.iseq],line->pos+1,reg->from,reg->to);
		}

		// Indel length distribution
		int *ptr = stats->insertions;
		if ( len<0 ) 
		{
			len *= -1;
			ptr = stats->deletions;
		}
		if ( --len > stats->m_indel ) len = stats->m_indel-1;
		ptr[len]++;
	}
}
void do_snp_stats(args_t *args, stats_t *stats, reader_t *reader)
{
	stats->n_snps++;

	bcf1_t *line = reader->line;
	int ref = acgt2int(*line->d.allele[0]);
	if ( ref<0 ) return;

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
			stats->af_ts[iaf]++;
		else 
			stats->af_tv[iaf]++;
	}

	if ( args->files.n_smpl )
	{
		if ( !set_fmt_ptr(reader,"GT") ) return;

		int is;
		for (is=0; is<args->files.n_smpl; is++)
		{
			int ial;
			int gt = gt_type(reader->fmt_ptr, reader->samples[is], &ial);
			if ( gt == GT_UNKN ) continue;
			if ( gt == GT_HET_RA || gt == GT_HET_AA ) stats->smpl_hets[is]++;
			if ( gt != GT_HOM_RR )
			{
				if ( !(line->d.var[ial].type&VCF_SNP) ) continue;
				int alt = acgt2int(*line->d.allele[ial]);
				if ( alt<0 ) continue;
				if ( abs(ref-alt)==2 ) 
					stats->smpl_ts[is]++;
				else
					stats->smpl_tv[is]++;
			}
		}
	}
}

void do_sample_stats(args_t *args)
{
	readers_t *files = &args->files;

	if ( files->nreaders>1 )
	{
		int is,ir;
		for (ir=0; ir<files->nreaders; ir++)
			if ( !set_fmt_ptr(&files->readers[ir],"GT") ) return;

		int iaf = args->tmp_iaf[1];	// only first ALT alelle considered
		gtcmp_t *af_stats = files->readers[0].line->d.var_type&VCF_SNP ? args->af_gts_snps : args->af_gts_indels;
		gtcmp_t *smpl_stats = files->readers[0].line->d.var_type&VCF_SNP ? args->smpl_gts_snps : args->smpl_gts_indels;

		for (is=0; is<files->n_smpl; is++)
		{
			// Simplified comparison: only 0/0, 0/1, 1/1 is looked at as the identity of 
			//	actual alleles can be enforced by running without the -c option.
			int gt = gt_type(files->readers[0].fmt_ptr, files->readers[0].samples[is], NULL);
			if ( gt == GT_UNKN ) continue;
			int match = 1;
			for (ir=1; ir<files->nreaders; ir++)
			{
				if ( gt != gt_type(files->readers[ir].fmt_ptr, files->readers[ir].samples[is], NULL) ) { match = 0; break; }
			}
			if ( gt == GT_HET_AA ) gt = GT_HOM_AA;	// rare case, treat as AA hom
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
		}
	}
}

void check_vcf(args_t *args)
{
	int ret,i;
	readers_t *files = &args->files;
	while ( (ret=next_line(files)) )
	{
		reader_t *reader = NULL;
		bcf1_t *line = NULL;
		for (i=0; i<files->nreaders; i++)
		{
			if ( !(ret&1<<i) ) continue;
			reader = &files->readers[i];
			line = files->readers[i].line;
			break;
		}
		set_variant_types(line);
		init_iaf(args, reader);

		if ( line->d.var_type&VCF_SNP ) 
			do_snp_stats(args, &args->stats[ret-1], reader);
		if ( line->d.var_type&VCF_INDEL )
			do_indel_stats(args, &args->stats[ret-1], reader);

		if ( line->n_allele>2 ) args->stats[ret-1].n_mals++;

		if ( files->n_smpl && ret==3 )
			do_sample_stats(args);
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
	if ( args->files.nreaders==1 )
		printf("ID\t0\t%s\n", args->files.readers[0].fname);
	else
	{
		printf("ID\t0\t%s\n", args->files.readers[0].fname);
		printf("ID\t1\t%s\n", args->files.readers[1].fname);
		printf("ID\t2\t%s\t%s\n", args->files.readers[0].fname,args->files.readers[1].fname);
	}
}

void print_stats(args_t *args)
{
	int i, id, nstats = args->files.nreaders==1 ? 1 : 3;
	printf("# Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
	for (id=0; id<args->files.nreaders; id++)
		printf("SN\t%d\tnumber of samples:\t%d\n", id, args->files.readers[id].header->n[BCF_DT_SAMPLE]);
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		printf("SN\t%d\tnumber of SNPs:\t%d\n", id, stats->n_snps);
		printf("SN\t%d\tnumber of indels:\t%d\n", id, stats->n_indels);
		printf("SN\t%d\tnumber of multiallelic sites:\t%d\n", id, stats->n_mals);

		int ts=0,tv=0;
		for (i=0; i<args->m_af; i++) { ts += stats->af_ts[i]; tv += stats->af_tv[i];  }
		printf("SN\t%d\tts/tv:\t%.2f\n", id, tv?(float)ts/tv:0);
	}
	printf("# Indel frameshifts:\n# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]out/(in+out) ratio\n");
	for (id=0; id<nstats; id++)
	{
		int in=args->stats[id].in_frame, out=args->stats[id].out_frame;
		printf("FS\t%d\t%d\t%d\t%.2f\n", id, in,out,out?(float)out/(in+out):0);
	}
	printf("# Singleton stats:\n# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		printf("SiS\t%d\t%d\t%d\t%d\t%d\t%d\n", id,1,stats->af_snps[0],stats->af_ts[0],stats->af_tv[0],stats->af_indels[0]);
		stats->af_snps[1]   += stats->af_snps[0];
		stats->af_ts[1]     += stats->af_ts[0];
		stats->af_tv[1]     += stats->af_tv[0];
		stats->af_indels[1] += stats->af_indels[0];
	}
	printf("# Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=1; i<args->m_af; i++)
		{
			if ( stats->af_snps[i]+stats->af_ts[i]+stats->af_tv[i]+stats->af_indels[i] == 0  ) continue;
			printf("AF\t%d\t%f\t%d\t%d\t%d\t%d\n", id,100.*(i-1)/(args->m_af-2),stats->af_snps[i],stats->af_ts[i],stats->af_tv[i],stats->af_indels[i]);
		}
	}
	printf("# InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]count\n");
	for (id=0; id<nstats; id++)
	{
		stats_t *stats = &args->stats[id];
		for (i=stats->m_indel-1; i>=0; i--)
			if ( stats->deletions[i] ) printf("IDD\t%d\t%d\t%d\n", id,-i-1,stats->deletions[i]);
		for (i=0; i<stats->m_indel; i++)
			if ( stats->insertions[i] ) printf("IDD\t%d\t%d\t%d\n", id,i+1,stats->insertions[i]);
	}
	printf("# Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
	for (id=0; id<nstats; id++)
	{
		int t;
		for (t=0; t<15; t++)
		{
			if ( t>>2 == (t&3) ) continue;
			printf("ST\t%d\t%c>%c\t%d\n", id, int2acgt(t>>2),int2acgt(t&3),args->stats[id].subst[t]);
		}
	}
	if ( args->files.nreaders>1 && args->files.n_smpl )
	{
		printf("# Genotype concordance by non-reference allele frequency (SNPs)\n# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\n");
		gtcmp_t *stats = args->af_gts_snps;
		int ndr_m=0, ndr_mm=0;
		for (i=1; i<args->m_af; i++)
		{
			int j, n = 0;
			for (j=0; j<3; j++) n += stats[i].m[j] + stats[i].mm[j];
			if ( !n ) continue;
			printf("GCsAF\t2\t%f", 100.*(i-1)/(args->m_af-2));
			printf("\t%d\t%d\t%d", stats[i].m[GT_HOM_RR],stats[i].m[GT_HET_RA],stats[i].m[GT_HOM_AA]);
			printf("\t%d\t%d\t%d\n", stats[i].mm[GT_HOM_RR],stats[i].mm[GT_HET_RA],stats[i].mm[GT_HOM_AA]);
			ndr_m  += stats[i].m[GT_HET_RA] + stats[i].m[GT_HOM_AA];
			ndr_mm += stats[i].mm[GT_HOM_RR] + stats[i].mm[GT_HET_RA] + stats[i].mm[GT_HOM_AA];
		}

		printf("SN\t%d\tNon-reference Discordance Rate (NDR):\t%.2f\n", 2, ndr_m+ndr_mm ? ndr_mm*100.0/(ndr_m+ndr_mm) : 0);
		printf("SN\t%d\tNumber of samples:\t%d\n", 2, args->files.n_smpl);

		printf("# Genotype concordance by sample (SNPs)\n# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\n");
		stats = args->smpl_gts_snps;
		for (i=0; i<args->files.n_smpl; i++)
		{
			int m  = stats[i].m[GT_HET_RA] + stats[i].m[GT_HOM_AA];
			int mm = stats[i].mm[GT_HOM_RR] + stats[i].mm[GT_HET_RA] + stats[i].mm[GT_HOM_AA];
			printf("GCsS\t2\t%s\t%.3f", args->files.samples[i], m+mm ? mm*100.0/(m+mm) : 0);
			printf("\t%d\t%d\t%d", stats[i].m[GT_HOM_RR],stats[i].m[GT_HET_RA],stats[i].m[GT_HOM_AA]);
			printf("\t%d\t%d\t%d\n", stats[i].mm[GT_HOM_RR],stats[i].mm[GT_HET_RA],stats[i].mm[GT_HOM_AA]);
		}
	}

	if ( args->files.n_smpl )
	{
		printf("# Per-sample counts\n# PSC\t[2]id\t[3]sample\t[4]nHets\t[5]nTransversions\t[6]nTransitions\n");
		for (id=0; id<nstats; id++)
		{
			stats_t *stats = &args->stats[id];
			for (i=0; i<args->files.n_smpl; i++)
			{
				if ( 0 == stats->smpl_hets[i] + stats->smpl_ts[i] + stats->smpl_tv[i] ) continue;
				printf("PSC\t%d\t%s\t%d\t%d\t%d\n", id,args->files.samples[i], stats->smpl_hets[i], stats->smpl_ts[i], stats->smpl_tv[i]);
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
	fprintf(stderr, "    -e, --exons <file.gz>             tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed)\n");
	fprintf(stderr, "    -f, --apply-filters               skip sites where FILTER is other than PASS\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>    collect statistics in the given region only\n");
	fprintf(stderr, "    -s, --samples <list|file>         create sample stats, \"-\" to include all samples\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfcheck(int argc, char *argv[])
{
	int c;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"collapse",0,0,'c'},
		{"apply-filters",0,0,'f'},
		{"exons",0,0,'e'},
		{"samples",0,0,'s'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hc:fr:e:s:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'c':
				if ( !strcmp(optarg,"snps") ) args->files.collapse |= COLLAPSE_SNPS;
				else if ( !strcmp(optarg,"indels") ) args->files.collapse |= COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"both") ) args->files.collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
				else if ( !strcmp(optarg,"any") ) args->files.collapse |= COLLAPSE_ANY;
				break;
			case 'f': args->files.apply_filters = 1; break;
			case 'r': args->files.region = optarg; break;
			case 'e': args->exons_file = optarg; break;
			case 's': args->samples_file = optarg; break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if (argc == optind) usage();

	if ( argc-optind>2 ) usage();
	while (optind<argc)
	{
		if ( !add_reader(argv[optind], &args->files) ) error("Could not load the index: %s\n", argv[optind]);
		optind++;
	}

	init_stats(args);
	print_header(args);
	check_vcf(args);
	print_stats(args);
	destroy_stats(args);
	destroy_readers(&args->files);
	free(args);
	return 0;
}

