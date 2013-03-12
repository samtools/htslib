#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

typedef struct
{
	bcf_srs_t *files;
    bcf_hdr_t *gt_hdr, *sm_hdr;
    double *lks, *sigs;
    int *cnts, *dps, hom_only;
	char **argv, *gt_fname, *plot, *sample;
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

FILE *open_file(char **fname, const char *mode, const char *fmt, ...);
void py_plot(char *script);
char *msprintf(const char *fmt, ...);

static void do_plot(args_t *args)
{
    char *fname;
    FILE *fp = open_file(&fname, "w", "%s.py", args->plot);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import matplotlib.gridspec as gridspec\n"
            "import csv\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "dat = []\n"
            "with open('%s.tab', 'rb') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        if row[0][0]=='#': continue\n"
            "        tgt = 0\n"
            "        if row[4]=='%s': tgt = 1\n"
            "        dat.append([float(row[0]), float(row[1]), float(row[2]), tgt])\n"
            "\n"
            "dat = sorted(dat, reverse=True)\n"
            "\n"
            "iq = -1; dp = 0\n"
            "for i in range(len(dat)):\n"
            "    if iq==-1 and dat[i][3]==1: iq = i\n"
            "    dp += dat[i][2]\n"
            "dp /= len(dat)\n"
            "\n"
            "fig = plt.figure()\n"
            "gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1])\n"
            "ax1 = plt.subplot(gs[0])\n"
            "ax3 = plt.subplot(gs[1])\n"
            "ax2 = ax1.twinx()\n"
            "ax1.plot([x[0] for x in dat],'g-')\n"
            "ax2.plot([x[1] for x in dat], '^', ms=3, color='r', mec='r')\n"
            "ax3.plot([x[2] for x in dat],'^', ms=3, color='k')\n"
            "if iq!=-1:\n"
            "   ax1.plot([iq],[dat[iq][0]],'ko', ms=4)\n"
            "   ax1.annotate('%s',xy=(iq,dat[iq][0]), xytext=(5,5), textcoords='offset points',fontsize='xx-small',rotation=45,va='bottom',ha='left')\n"
            "   ax2.plot([iq],[dat[iq][1]],'^', ms=4)\n"
            "for tl in ax1.get_yticklabels(): tl.set_color('g')\n"
            "for tl in ax2.get_yticklabels(): tl.set_color('r')\n"
            "ax1.set_title('Average dp = %%.1fx' %% dp)\n"
            "ax1.set_xticks([])\n"
            "ax3.set_xlabel('Sample ID')\n"
            "ax3.set_ylabel('Avg DP')\n"
            "ax1.set_ylabel('Concordance',color='g')\n"
            "ax2.set_ylabel('Uncertainty',color='r')\n"
            "ax1.set_zorder(ax2.get_zorder()+1)\n"
            "ax1.patch.set_visible(False)\n"
            "plt.subplots_adjust(hspace=0.06)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", args->plot, args->sample ? args->sample : "", args->sample ? args->sample : "", args->plot
           );
    fclose(fp);
    py_plot(fname);
    free(fname);
}

static void init_data(args_t *args)
{
    args->gt_hdr = args->files->readers[0].header;
    args->sm_hdr = args->files->readers[1].header;
    if ( !args->gt_hdr->n[BCF_DT_SAMPLE] || !args->sm_hdr->n[BCF_DT_SAMPLE] ) error("No samples?\n");
    args->lks  = (double*) calloc(args->gt_hdr->n[BCF_DT_SAMPLE],sizeof(double));
    args->sigs = (double*) calloc(args->gt_hdr->n[BCF_DT_SAMPLE],sizeof(double));
    args->cnts = (int*) calloc(args->gt_hdr->n[BCF_DT_SAMPLE],sizeof(int));
    args->dps  = (int*) calloc(args->gt_hdr->n[BCF_DT_SAMPLE],sizeof(int));
}

static void destroy_data(args_t *args)
{
    free(args->lks); free(args->sigs); free(args->cnts); free(args->dps);
}

static int allele_to_int(bcf1_t *line, char *allele)
{
    int i;
    for (i=0; i<line->n_allele; i++)
        if ( !strcmp(allele,line->d.allele[i]) ) return i;
    if ( strcmp(line->d.allele[i-1],"X") ) return -1;
    return i-1;
}

static int init_gt2ipl(args_t *args, bcf1_t *gt_line, bcf1_t *sm_line, int *gt2ipl, int n_gt2ipl)
{
    int i, j;
    for (i=0; i<n_gt2ipl; i++) gt2ipl[i] = -1;
    for (i=0; i<gt_line->n_allele; i++)
    {
        // find which of the sm_alleles (k) corresponds to the gt_allele (i)
        int k = allele_to_int(sm_line, gt_line->d.allele[i]);
        if ( k<0 ) return 0;
        for (j=0; j<=i; j++)
        {
            int l = allele_to_int(sm_line, gt_line->d.allele[j]);
            if ( l<0 ) return 0;
            gt2ipl[ bcf_ij2G(j,i) ] = k<=l ? bcf_ij2G(k,l) : bcf_ij2G(l,k);
        }
    }
    //for (i=0; i<n_gt2ipl; i++) printf("%d .. %d\n", i,gt2ipl[i]);
    return 1;
}

static void check_gt(args_t *args)
{
    int i,j, ret, *gt2ipl = NULL, m_gt2ipl = 0, *gt_cnts = NULL;
    int gt_id = bcf_id2int(args->gt_hdr, BCF_DT_ID, "GT");
    int pl_id = bcf_id2int(args->sm_hdr, BCF_DT_ID, "PL");
    int dp_id = bcf_id2int(args->sm_hdr, BCF_DT_ID, "DP");
    if ( gt_id<0 ) error("GT not present in the header?\n");
    if ( pl_id<0 ) error("PL not present in the header?\n");
    if ( dp_id<0 ) error("DP not present in the header?\n");
    int isample = -1;
    if ( args->sample ) 
    {
        isample = bcf_id2int(args->gt_hdr, BCF_DT_SAMPLE, args->sample);
        if ( isample<0 ) error("No such sample: [%s]\n", args->sample);
        if ( args->plot ) isample = -1; // different kind of output with -p
    }
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        if ( ret!=3 ) continue;
        bcf1_t *gt_line = args->files->readers[0].buffer[0];
        bcf1_t *sm_line = args->files->readers[1].buffer[0];
        bcf_unpack(gt_line, BCF_UN_FMT);
        bcf_unpack(sm_line, BCF_UN_ALL);

        // Init the map between genotype and PL fields, assuming diploid genotypes
        int n_gt2ipl = gt_line->n_allele*(gt_line->n_allele + 1)/2;
        if ( n_gt2ipl > m_gt2ipl )
        {
            m_gt2ipl = n_gt2ipl;
            gt2ipl   = (int*) realloc(gt2ipl, sizeof(int)*m_gt2ipl);
            gt_cnts  = (int*) realloc(gt_cnts, sizeof(int)*(m_gt2ipl+1));
        }
        if ( !init_gt2ipl(args, gt_line, sm_line, gt2ipl, n_gt2ipl) ) continue;

        // Get BCF handler for GT
        bcf_fmt_t *gt_fmt = NULL, *pl_fmt = NULL;
        for (i=0; i<(int)gt_line->n_fmt; i++)  
        {
            if ( gt_line->d.fmt[i].id==gt_id )
            {
                gt_fmt = &gt_line->d.fmt[i];
                break;
            }
        }
        if ( !gt_fmt ) error("GT not present at %s:%d?", args->gt_hdr->id[BCF_DT_CTG][gt_line->rid].key, gt_line->pos+1);
        assert( gt_fmt->n==2 && gt_fmt->size==2*sizeof(int8_t) );

        // Get BCF handler for PL
        for (i=0; i<(int)sm_line->n_fmt; i++)  
        {
            if ( sm_line->d.fmt[i].id==pl_id ) 
            {
                pl_fmt = &sm_line->d.fmt[i];
                break;
            }
        }
        if ( !pl_fmt ) error("PL not present at %s:%d?", args->sm_hdr->id[BCF_DT_CTG][sm_line->rid].key, sm_line->pos+1);

        // Set counts of genotypes to test significance
        for (i=0; i<=n_gt2ipl; i++) gt_cnts[i] = 0;
        for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++)
        {
            int8_t *gt_ptr = (int8_t*)(gt_fmt->p + i*gt_fmt->size); /* FIXME: does not work with n_alt >= 64 */
            int a = (gt_ptr[0]>>1) - 1;
            int b = (gt_ptr[1]>>1) - 1;
            if ( a<0 || b<0 ) 
                gt_cnts[n_gt2ipl - 1]++;
            else
            {
                int igt = a<=b ? bcf_ij2G(a,b) : bcf_ij2G(b,a);
                gt_cnts[igt]++;
            }
        }

        // With -s but no -p, print LKs at all sites for debugging
        int igt = -1;
        if ( isample>=0 )
        {
            int8_t *gt_ptr = (int8_t*)(gt_fmt->p + isample*gt_fmt->size);
            int a = (gt_ptr[0]>>1) - 1;
            int b = (gt_ptr[1]>>1) - 1; 
            printf("%s\t%d", args->gt_hdr->id[BCF_DT_CTG][gt_line->rid].key, gt_line->pos+1);
            for (i=0; i<gt_line->n_allele; i++) printf("%c%s", i==0?'\t':',', gt_line->d.allele[i]);
            printf("\t%d\t%s/%s", sm_line->d.info[dp_id].v1.i, a>=0 ? gt_line->d.allele[a] : ".", b>=0 ? gt_line->d.allele[b] : ".");
            if (a>=0 || b>=0) igt = a<=b ? bcf_ij2G(a,b) : bcf_ij2G(b,a);
        }

        // Calculate likelihoods for all samples, assuming diploid genotypes
        #define BRANCH(type_t, pl_is_missing) { \
            type_t *pl_ptr = (type_t*) pl_fmt->p; /* for now, only the first sample */ \
            if ( isample<0 ) \
            { \
                for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++) \
                { \
                    int8_t *gt_ptr = (int8_t*)(gt_fmt->p + i*gt_fmt->size); /* FIXME: does not work with n_alt >= 64 */ \
                    int a = (gt_ptr[0]>>1) - 1; \
                    int b = (gt_ptr[1]>>1) - 1; \
                    if ( a<0 || b<0 ) continue; /* missing genotype */ \
                    if ( args->hom_only && a!=b ) continue; /* heterozygous genotype */ \
                    \
                    igt = a<=b ? bcf_ij2G(a,b) : bcf_ij2G(b,a); \
                    assert( igt<m_gt2ipl ); \
                    args->sigs[i] += -log(1./(gt_cnts[igt]+gt_cnts[n_gt2ipl-1])); \
                    igt = gt2ipl[igt]; \
                    if ( !(pl_is_missing) ) { \
                        double sum = 0; for (j=0; j<pl_fmt->n; j++) sum += pow(10, -0.1*pl_ptr[j]); \
                        args->lks[i] += -log(pow(10, -0.1*pl_ptr[igt])/sum); \
                        args->cnts[i]++; \
                    } \
                    args->dps[i] += sm_line->d.info[dp_id].v1.i; \
                } \
            } \
            else \
            { \
                type_t *pl_ptr = (type_t*) pl_fmt->p; /* for now, only the first sample */ \
                printf("\t%d", igt<0 ? '.' : (int)pl_ptr[ gt2ipl[igt] ]); \
                for (i=0; i<sm_line->n_allele; i++) printf("%c%s", i==0?'\t':',', sm_line->d.allele[i]); \
                for (igt=0; igt<pl_fmt->n && !(pl_is_missing); igt++) printf("\t%d", (int)pl_ptr[igt]); \
                printf("\n"); \
            } \
        }
        switch (pl_fmt->type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  pl_ptr[igt]==INT8_MIN); break;
            case BCF_BT_INT16: BRANCH(int16_t, pl_ptr[igt]==INT16_MIN); break;
            case BCF_BT_INT32: BRANCH(int32_t, pl_ptr[igt]==INT32_MIN); break;
            case BCF_BT_FLOAT: BRANCH(float,  *(uint32_t*)(&pl_ptr[igt])==bcf_missing_float); break;
            default: error("todo: type %d\n", pl_fmt->type); break;
        }
        #undef BRANCH
    }
    if ( gt2ipl ) free(gt2ipl);
    if ( gt_cnts ) free(gt_cnts);

    // Scale LKs and certainties
    double max = args->lks[0];
    for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++) if ( max<args->lks[i] ) max = args->lks[i];
    for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++) args->lks[i] = (max - args->lks[i]) / max;
    max = args->sigs[0];
    for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++) if ( max<args->sigs[i] ) max = args->sigs[i];
    for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++) args->sigs[i] = (max - args->sigs[i]) / max;

    // Output
    if ( args->plot )
    {
        FILE *fp = open_file(NULL, "w", "%s.tab", args->plot);
        fprintf(fp, "# [1]Concordance\t[2]Uncertainty\t[3]Average depth\t[4]Number of sites\t[5]Sample\n");
        for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++)
            fprintf(fp, "%f\t%f\t%.1f\t%d\t%s\n", args->lks[i], args->sigs[i], args->cnts[i]?(double)args->dps[i]/args->cnts[i]:0.0, args->cnts[i], args->gt_hdr->samples[i]);
        fclose(fp);
        do_plot(args);
    }
    else if ( isample<0 )
    {
        printf("# [1]Concordance\t[2]Uncertainty\t[3]Average depth\t[4]Number of sites\t[5]Sample\n");
        for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++)
            printf("%f\t%f\t%.1f\t%d\t%s\n", args->lks[i], args->sigs[i], args->cnts[i]?(double)args->dps[i]/args->cnts[i]:0.0, args->cnts[i], args->gt_hdr->samples[i]);
    }
}

static char *init_prefix(char *prefix)
{
    int len = strlen(prefix);
    if ( prefix[len-1] == '/' || prefix[len-1] == '\\' )
        return msprintf("%sgtcheck", prefix);
    return strdup(prefix);
}

static void usage(void)
{
	fprintf(stderr, "About:   Check sample identity\n");
	fprintf(stderr, "Usage:   vcfgtcheck [options] -g <genotypes.vcf.gz> <file.vcf.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -g, --genotypes <file>     genotypes to compare against in VCF\n");
	fprintf(stderr, "    -H, --homs-only            homozygous genotypes only for very low coverage data\n");
	fprintf(stderr, "    -p, --plot <prefix>        plot\n");
	fprintf(stderr, "    -s, --sample <string>      target sample (used for plotting only)\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfgtcheck(int argc, char *argv[])
{
	int c;
	args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
	args->argc   = argc; args->argv = argv;

	static struct option loptions[] = 
	{
		{"homs-only",0,0,'H'},
		{"help",0,0,'h'},
		{"genotypes",1,0,'g'},
		{"plot",1,0,'p'},
		{"sample",1,0,'s'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hg:p:s:H",loptions,NULL)) >= 0) {
		switch (c) {
			case 'H': args->hom_only = 1; break;
			case 'g': args->gt_fname = optarg; break;
			case 'p': args->plot = optarg; break;
			case 's': args->sample = optarg; break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if ( argc!=optind+1 || !args->gt_fname ) usage();   // none or too many files given
    if ( !bcf_sr_add_reader(args->files, args->gt_fname) ) error("Failed to open or the file not indexed: %s\n", args->gt_fname);
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    args->files->collapse = COLLAPSE_SNPS|COLLAPSE_INDELS;
    if ( args->plot ) args->plot = init_prefix(args->plot);
    init_data(args);
    check_gt(args);
    destroy_data(args);
	bcf_sr_destroy(args->files);
    if (args->plot) free(args->plot);
	free(args);
	return 0;
}

