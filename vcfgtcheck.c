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
    int *cnts, *dps, hom_only, cross_check;
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

static void plot_check(args_t *args)
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

static void plot_cross_check(args_t *args)
{
    char *fname;
    FILE *fp = open_file(&fname, "w", "%s.py", args->plot);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "avg   = []\n"
            "sm2id = {}\n"
            "dat   = None\n"
            "min   = None\n"
            "with open('%s.tab', 'rb') as f:\n"
            "   reader = csv.reader(f, 'tab')\n"
            "   i = 0\n"
            "   for row in reader:\n"
            "       if row[0]=='SM':\n"
            "           sm2id[row[2]] = i\n"
            "           avg.append([i,float(row[1])])\n"
            "           i += 1\n"
            "       elif row[0]=='CN':\n"
            "           val = float(row[1])/int(row[2])\n"
            "           if not dat:\n"
            "               dat = [[0]*len(sm2id) for x in xrange(len(sm2id))]\n"
            "               min = val\n"
            "           id_i = sm2id[row[4]]\n"
            "           id_j = sm2id[row[5]]\n"
            "           if id_i<id_j: dat[id_i][id_j] = val\n"
            "           else: dat[id_j][id_i] = val\n"
            "           if min > val: min = val\n"
            "\n"
            "if len(sm2id)<=1: exit(1)\n"
            "\n"
            "fig = plt.figure()\n"
            "fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,7))\n"
            "\n"
            "ax1.plot([x[0] for x in avg],[x[1] for x in avg],'^-', ms=3, color='k')\n"
            "im = ax2.imshow(dat,clim=(min),interpolation='nearest')\n"
            "cb1  = plt.colorbar(im,shrink=0.8)\n"
            "cb1.set_label('Pairwise discordance')\n"
            "for t in cb1.ax.get_yticklabels(): t.set_fontsize(9)\n"
            "\n"
            "ax1.tick_params(axis='both', which='major', labelsize=9)\n"
            "ax1.tick_params(axis='both', which='minor', labelsize=9)\n"
            "ax2.tick_params(axis='both', which='major', labelsize=9)\n"
            "ax2.tick_params(axis='both', which='minor', labelsize=9)\n"
            "\n"
            "ax1.set_title('Sample Discordance Score')\n"
            "ax2.set_ylabel('Sample ID')\n"
            "ax2.set_xlabel('Sample ID')\n"
            "ax1.set_xlabel('Sample ID')\n"
            "ax1.set_ylabel('Average discordance')\n"
            "\n"
            "plt.subplots_adjust(left=0.15,right=0.95,bottom=0.08,top=0.95,hspace=0.25)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", args->plot,args->plot
           );
    fclose(fp);
    py_plot(fname);
    free(fname);
}

static void init_data(args_t *args)
{
    args->gt_hdr = args->files->readers[0].header;
    int nsamples = args->gt_hdr->n[BCF_DT_SAMPLE];
    if ( !nsamples ) error("No samples?\n");

    if ( !args->cross_check )
    {
        args->sm_hdr = args->files->readers[1].header;
        if ( !args->sm_hdr->n[BCF_DT_SAMPLE] ) error("No samples?\n");
        args->lks  = (double*) calloc(nsamples,sizeof(double));
        args->sigs = (double*) calloc(nsamples,sizeof(double));
        args->cnts = (int*) calloc(nsamples,sizeof(int));
        args->dps  = (int*) calloc(nsamples,sizeof(int));
    }
    else
    {
        int narr = (nsamples-1)*nsamples/2;
        args->lks  = (double*) calloc(narr,sizeof(double));
        args->cnts = (int*) calloc(narr,sizeof(int));
        args->dps  = (int*) calloc(narr,sizeof(int));
}
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
        if ( isample>=0 )
            printf("# [1]Chromosome\t[2]Position\t[3]Alleles in -g file\t[4]Coverage\t[5]Genotype\t[6]Alleles in sample file\t[7-]PL likelihoods\n");
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

        // Set the counts of genotypes to test significance
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
            if ( args->hom_only && a!=b ) continue; /* heterozygous genotype */
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
        plot_check(args);
    }
    else if ( isample<0 )
    {
        printf("# [1]Concordance\t[2]Uncertainty\t[3]Average depth\t[4]Number of sites\t[5]Sample\n");
        for (i=0; i<args->gt_hdr->n[BCF_DT_SAMPLE]; i++)
            printf("%f\t%f\t%.1f\t%d\t%s\n", args->lks[i], args->sigs[i], args->cnts[i]?(double)args->dps[i]/args->cnts[i]:0.0, args->cnts[i], args->gt_hdr->samples[i]);
    }
}

static int cmp_doubleptr(const void *_a, const void *_b)
{
    double *a = *((double**)_a);
    double *b = *((double**)_b);
    if ( *a < *b ) return -1;
    else if ( *a == *b ) return 0;
    return 1;
}


static void cross_check_gts(args_t *args)
{
    int nsamples = args->gt_hdr->n[BCF_DT_SAMPLE], ndp_arr = 0, npl_arr = 0;
    int i,j,k,idx, ret, *dp_arr = NULL, *pl_arr = NULL;
    int pl_id = bcf_id2int(args->gt_hdr, BCF_DT_ID, "PL");
    int dp_id = bcf_id2int(args->gt_hdr, BCF_DT_ID, "DP");
    if ( pl_id<0 ) error("PL not present in the header?\n");
    if ( dp_id<0 ) error("DP not present in the header?\n");
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, BCF_UN_FMT);

        // Get BCF handler for PL and DP
        bcf_fmt_t *dp_fmt = NULL, *pl_fmt = NULL;
        for (i=0; i<(int)line->n_fmt; i++)  
        {
            if ( line->d.fmt[i].id==pl_id ) pl_fmt = &line->d.fmt[i];
            if ( line->d.fmt[i].id==dp_id ) dp_fmt = &line->d.fmt[i];
        }
        if ( !pl_fmt ) error("PL not present at %s:%d?", args->gt_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
        if ( !dp_fmt ) error("DP not present at %s:%d?", args->gt_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);

        dp_arr = bcf_set_iarray(dp_fmt, nsamples, dp_arr, &ndp_arr);
        pl_arr = bcf_set_iarray(pl_fmt, nsamples, pl_arr, &npl_arr);

        idx = 0;
        for (i=1; i<nsamples; i++)
        {
            int *ipl = &pl_arr[i*pl_fmt->n];
            if ( dp_arr[i]==INT_MIN || !dp_arr[i] || *ipl==INT_MIN ) 
            {
                idx += i;
                continue;
            }
            for (j=0; j<i; j++)
            {
                int *jpl = &pl_arr[j*pl_fmt->n];
                if ( dp_arr[j]==INT_MIN || !dp_arr[j] || *jpl==INT_MIN ) 
                {
                    idx++;
                    continue;
                }
                int min_pl = ipl[0] + jpl[0];
                for (k=1; k<pl_fmt->n; k++)
                    if ( min_pl > ipl[k]+jpl[k] ) min_pl = ipl[k]+jpl[k];
                args->lks[idx] += min_pl;
                args->cnts[idx]++;
                args->dps[idx] += dp_arr[i] < dp_arr[j] ? dp_arr[i] : dp_arr[j];
                idx++;
            }
        }
    }
    if ( dp_arr ) free(dp_arr);
    if ( pl_arr ) free(pl_arr);

    FILE *fp = args->plot ? open_file(NULL, "w", "%s.tab", args->plot) : stdout;
    // Output samples sorted by average discordance
    double *score = (double*) calloc(nsamples,sizeof(double));
    idx = 0;
    for (i=1; i<nsamples; i++)
    {
        for (j=0; j<i; j++)
        {
            score[i] += args->cnts[idx] ? (double)args->lks[idx]/args->cnts[idx] : 0;
            score[j] += args->cnts[idx] ? (double)args->lks[idx]/args->cnts[idx] : 0;
            idx++;
        }
    }
    double **p = (double**) malloc(sizeof(double*)*nsamples), avg_score = 0;
    for (i=0; i<nsamples; i++) p[i] = &score[i];
    qsort(p, nsamples, sizeof(int*), cmp_doubleptr);
    fprintf(fp, "# [1]SM\t[2]Average Discordance/Number of sites\t[3]Sample\n");
    for (i=0; i<nsamples; i++)
    {
        idx = p[i] - score;
        double tmp = (double)score[idx]/nsamples;
        avg_score += tmp;
        fprintf(fp, "SM\t%lf\t%s\n", tmp, args->gt_hdr->samples[idx]);
    }

    // Overall score: maximum absolute deviation from the average score
    fprintf(fp, "# [1] MD\t[2]Maximum deviation\t[3]The culprit\n");
    fprintf(fp, "MD\t%f\t%s\n", (double)score[idx]/nsamples - avg_score/nsamples, args->gt_hdr->samples[idx]);    // idx still set
    free(p);
    free(score);

    // Pairwise discordances
    fprintf(fp, "# [1]CN\t[2]Discordance\t[3]Number of sites\t[4]Average minimum depth\t[5]Sample i\t[6]Sample j\n");
    double avg = 0;
    idx = 0;
    for (i=0; i<nsamples; i++)
    {
        for (j=0; j<i; j++)
        {
            avg += args->lks[idx];
            fprintf(fp, "CN\t%.0f\t%d\t%.6f\t%s\t%s\n", args->lks[idx], args->cnts[idx], args->cnts[idx]?(double)args->dps[idx]/args->cnts[idx]:0.0, 
                    args->gt_hdr->samples[i],args->gt_hdr->samples[j]);
            idx++;
        }
    }
    fclose(fp);
    if ( args->plot )
        plot_cross_check(args);
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
	fprintf(stderr, "About:   Check sample identity. With no VCF given, multi-sample cross-check is performed.\n");
	fprintf(stderr, "Usage:   vcfgtcheck [options] -g <genotypes.vcf.gz> [<file.vcf.gz>]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -g, --genotypes <file>             genotypes to compare against in VCF\n");
	fprintf(stderr, "    -H, --homs-only                    homozygous genotypes only for very low coverage data\n");
	fprintf(stderr, "    -p, --plot <prefix>                plot\n");
    fprintf(stderr, "    -r, --region <chr|chr:from-to>     perform the check in the given region only\n");
	fprintf(stderr, "    -s, --sample <string>              target sample (used for plotting only)\n");
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
        {"region",1,0,'r'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hg:p:s:Hr:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'H': args->hom_only = 1; break;
			case 'g': args->gt_fname = optarg; break;
			case 'p': args->plot = optarg; break;
			case 's': args->sample = optarg; break;
            case 'r': args->files->region = optarg; break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
    if ( !args->gt_fname ) usage();
    if ( argc==optind ) args->cross_check = 1;
	else if ( argc>optind+1 ) usage();   // too many files given
    if ( !args->cross_check && !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    if ( !bcf_sr_add_reader(args->files, args->gt_fname) ) error("Failed to open or the file not indexed: %s\n", args->gt_fname);
    args->files->collapse = COLLAPSE_SNPS|COLLAPSE_INDELS;
    if ( args->plot ) args->plot = init_prefix(args->plot);
    init_data(args);
    if ( args->cross_check )
        cross_check_gts(args);
    else
        check_gt(args);
    destroy_data(args);
	bcf_sr_destroy(args->files);
    if (args->plot) free(args->plot);
	free(args);
	return 0;
}

