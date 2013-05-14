#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "vcfutils.h"

#define NFIXED 5
#define MASK_GOOD 2
#define MASK_BAD  4
#define IS_GOOD(mask) ((mask)&MASK_GOOD) 
#define IS_BAD(mask) ((mask)&MASK_BAD)

typedef struct
{
    int nbin, kdim;     // number of bins and dimension
    int nt, t;          // total number of learning cycles and the current cycle
    double *w, *c;      // weights and counts (sum of learning influence)
    double radius, decay, learn, th;     // SOM parameters
}
som_t;

typedef struct
{
    unsigned int ngood, nbad, nall, nmissing;
    unsigned int nbins, *good_data, *all_data, *bad_data;  // used only for plotting
    double good_min, good_max;      // extremes of sites in the good set
    double all_min, all_max;        // overall extremes
    double bad_min, bad_max;        // overall extremes
    double scale_min, scale_max;    // 0.5 and 0.95 percentile
}
dist_t;

#define FLT_LE  2       // less or equal
#define FLT_LT  1       // less than
#define FLT_EQ  0       // equal
#define FLT_BT -1       // bigger than
#define FLT_BE -2       // bigger or equal

typedef struct
{
    int type;       // one of the FLT_* keys above
    double value; 
}
filter_t;

typedef struct
{
    int *nfilt, n;
    filter_t **filt;
}
filters_t;

typedef struct
{   
    char **names;               // annotation names (0..nann-1), list of annotations actually used
    char **colnames;            // all columns' names
    int *col2names;             // mapping from column number to i-th annot. The index of unused annots is set to -1. Includes first NFIXED columns.
    int *ann2cols;              // mapping from i-th annot to j-th column of annot.tab.gz. First annotation has index NFIXED.
    int nann, nann_som;         // number of used annotations (total, included fixed filters) and for SOM
    int ncols;                  // number of columns (total, including first NFIXED columns)
    dist_t *dists;              // annots distributions (all annots, including those not requested)
    int scale;                  // should the annotations be rescaled according to lo_,hi_pctl?
    double lo_pctl, hi_pctl;    // the percentage of sites to zoom out at both ends

    // annots_reader_* data
    htsFile *file;              // reader
    kstring_t str;              // temporary string for annots_reader_*
    int *ignore;                // columns to ignore (ncols, ro)
    int pos, mask, is_ts, type; // filled by annots_reader_next
    char *chr;                  // filled by annots_reader_next
    double *vals;               // filled by annots_reader_next (ncols), scaled if scale set
    double *raw_vals;           // same as vals, unscaled
    int *missing;               // columns with missing values, filled by annots_reader_next (ncols)
    int nset, nset_mask;        // number of filled coordinates and their mask, filled by annots_reader_next

    som_t som, bad_som;
    int filt_type;
    filters_t filt_excl, filt_learn;
    double snp_th, indel_th;

    int rand_seed;
    char *annot_str, *fixed_filters, *learning_filters;
	char **argv, *fname, *out_prefix, *region;
	int argc, do_plot;
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

static void usage(void);
FILE *open_file(char **fname, const char *mode, const char *fmt, ...);
void mkdir_p(const char *fmt, ...);

/**
 *  ks_getline() - Read next line from $fp, appending it to $str.  The newline
 *  is stripped and \0 appended. Returns the number of characters read
 *  excluding the null byte.
 */
size_t ks_getline(FILE *fp, kstring_t *str)
{
    size_t nread=0;
    int c;
    while ((c=getc(fp))!= EOF && c!='\n')
    {
        nread++;
        if ( str->l+nread > str->m )
        {
            str->m += 1024;
            str->s = (char*) realloc(str->s, sizeof(char)*str->m);
        }
        str->s[str->l+nread-1] = c;
    }
    if ( str->l >= str->m )
    {
        str->m += 1024;
        str->s = (char*) realloc(str->s, sizeof(char)*str->m);
    }
    str->l += nread;
    str->s[ str->l ] = 0;
    return nread;
}
char *msprintf(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *str = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(str, n, fmt, ap);
    va_end(ap);

    return str;
}
static char *_strndup(const char *ptr, int len)
{
    char *tmp = (char*) malloc(sizeof(char)*(len+1));
    int i;
    for (i=0; i<len; i++) 
    {
        if ( !ptr[i] ) break;
        tmp[i] = ptr[i];
    }
    tmp[i] = 0;
    return tmp;
}
char **read_list(char *fname, int *_n)
{
    int n = 0;
    char **list = NULL;

    FILE *fp = fopen(fname,"r");
    if ( !fp ) error("%s: %s\n", fname, strerror(errno));

    kstring_t str = {0,0,0};
    while ( ks_getline(fp, &str) )
    {
        list = (char**) realloc(list, sizeof(char*)*(++n));
        list[n-1] = strdup(str.s);
        str.l = 0;
    }
    fclose(fp);
    if ( str.m ) free(str.s);
    *_n = n;
    return list;
}
char **split_list(const char *str, int delim, int *n)
{
    *n = 0;
    char **out = NULL;
    const char *t = str, *p = t;
    while (*t)
    {
        if ( *t == delim ) 
        { 
            (*n)++;
            out = (char**) realloc(out, sizeof(char*)*(*n));
            out[(*n)-1] = _strndup(p,t-p);
            p = t+1;
        }
        t++;
    }
    (*n)++;
    out = (char**) realloc(out, sizeof(char*)*(*n));
    out[(*n)-1] = _strndup(p,t-p);
    return out;
}
void destroy_list(char **list, int n)
{
    int i;
    for (i=0; i<n; i++)
        free(list[i]);
    free(list);
}
void py_plot(char *script)
{
    mkdir_p(script);
    int len = strlen(script);
    char *cmd = !strcmp(".py",script+len-3) ? msprintf("python %s", script) : msprintf("python %s.py", script);
    int ret = system(cmd);
    if ( ret ) fprintf(stderr, "The command returned non-zero status %d: %s\n", ret, cmd);
    free(cmd);
}

/*
 *  char *t, *p = str;
 *  t = column_next(p, '\t'); 
 *  if ( strlen("<something>")==t-p && !strncmp(p,"<something>",t-p) ) printf("found!\n");
 *
 *  char *t;
 *  t = column_next(str, '\t'); if ( !*t ) error("expected field\n", str);
 *  t = column_next(t+1, '\t'); if ( !*t ) error("expected field\n", str);
 */
inline char *column_next(char *start, char delim)
{
    char *end = start;
    while (*end && *end!=delim) end++;
    return end;
}

inline double scale_value(dist_t *dist, double val)
{
    if ( val < dist->scale_min ) val = 0;
    else if ( val > dist->scale_max ) val = 1;
    else val = (val - dist->scale_min) / (dist->scale_max - dist->scale_min);
    assert( val>=0 && val<=1 );
    return val;
}
/**
 *  annots_reader_next() - reads next line from annots.tab.gz and sets the appropriate fields
 *
 *    mask:         2 if a good site, otherwise 1
 *    is_ts:        1 if transition, 0 if transversion (meaningless with indels)
 *    type:         one of VCF_SNP, VCF_INDEL, VCF_OTHER
 *    nset:         number of non-missing values read into the $vals array
 *    nset_mask:    which values are non-missing. Note that this limits the number of maximum annotations
 *    vals:         the parsed values the order of which corresponds to $names (scaled to <0,1> if $scaled is set)
 *    raw_vals:     same as while $vals, but unscaled
 *    missing:      array of flags indicating whether the value is set or not (the same order as above)
 *
 *   In addition also these self-explanatory fields are set: chr, pos
 *
 *   Returns 1 on successful read or 0 if no further record could be read.
 */
int annots_reader_next(args_t *args)
{
    args->str.l = 0;
    if ( hts_getline(args->file,'\n',&args->str)<=0 ) return 0;

    char *t, *line = args->str.s;

    // CHR
    t = column_next(line, '\t'); 
    if ( !*t ) error("Could not parse CHR: [%s]\n", line);
    if ( !args->ignore[0] )
    {
        *t = 0;
        args->chr = line;
    }

    // POS
    t++;
    if ( !args->ignore[1] )
        args->pos = strtol(t, NULL, 10);
    t = column_next(t, '\t');  
    if ( !*t ) error("Could not parse POS: [%s]\n", line);

    // MASK
    if ( !args->ignore[2] )
        args->mask = t[2]=='1' ? 2 : 1;     // good site t="\t11", otherwise t="\t10"
    t = column_next(t+1, '\t');  
    if ( !*t ) error("Could not parse MASK: [%s]\n", line);

    // IS_TS
    if ( !args->ignore[3] )
        args->is_ts = strtol(t, NULL, 10);
    t = column_next(t+1, '\t');  
    if ( !*t ) error("Could not parse IS_TS: [%s]\n", line);

    // TYPE
    if ( !args->ignore[4] )
    {
        if ( !strncmp(t+1,"SNP",3) || !strncmp(t+1,"MNP",3) ) args->type = VCF_SNP;
        else if ( !strncmp(t+1,"INDEL",5) ) args->type = VCF_INDEL;
        else args->type = VCF_OTHER;
    }
    t = column_next(t+1, '\t');  
    if ( !*t ) error("Could not parse TYPE: [%s]\n", line);

    args->nset = 0;
    args->nset_mask = 0;

    int icol, iann = -1;
    for (icol=NFIXED; icol<args->ncols; icol++)
    {
        if ( !*t ) error("Could not parse %d-th data field: is the line truncated?\nThe line was: [%s]\n",icol,line);

        if ( args->ignore[icol] )
        {
            // this column is to be ignored
            t = column_next(t+1,'\t');
            continue;
        }
        iann = args->col2names[icol];

        t++;
        if ( t[0]=='.' && (t[1]=='\t' || t[1]=='\n' || !t[1]) )
        {
            args->missing[iann] = 1;
            t++;
            continue;
        }
        if ( (sscanf(t,"%le", &args->vals[iann])!=1) ) 
            error("Could not parse %d-th data field: [%s]\nThe line was: [%s]\n", icol,t,line);

        if ( args->vals[iann] == HUGE_VAL || args->vals[iann] == -HUGE_VAL || isnan(args->vals[iann]) )
        {
            args->missing[iann] = 1;
            t = column_next(t,'\t');
            continue;
        }

        args->raw_vals[iann] = args->vals[iann];

        if ( args->scale && args->dists )
            args->vals[iann] = scale_value(&args->dists[icol], args->vals[iann]);

        args->nset++;
        args->nset_mask |= 1 << iann; 

        args->missing[iann] = 0;
        t = column_next(t,'\t');
    }
    return 1;
}
void annots_reader_reset(args_t *args)
{
    if ( args->file ) hts_close(args->file);
    if ( !args->fname ) error("annots_reader_reset: no fname\n");
    args->file = hts_open(args->fname, "r", NULL);
    hts_getline(args->file,'\n',&args->str);  // eat the header
}

static void smooth_dist(double *dst, unsigned int *src, int n)
{
    int i, j, win = 5;
    for (i=0; i<n; i++)
    {
        // dst[i] = src[i]; continue;  // no smoothing

        int jfrom = i - win;
        if ( jfrom<0 ) jfrom = 0;
        int jto = i + win;
        if ( jto>n ) jto = n;

        double value = 0, norm = 0;
        for (j=jfrom; j<jto; j++)
        {
            // Epanechnikov
            double tmp = (double)(i - j)/win;
            tmp = 3/4. * (1 - tmp*tmp);
            norm  += tmp;
            value += tmp * src[j];
        }
        dst[i] = norm ? value/norm : 0;
    }
}

static void plot_dists(args_t *args, dist_t *dists, int ndist)
{
    fprintf(stderr,"Plotting distributions...\n");

    char *fname;
    FILE *fp = open_file(&fname,"w","%s.dists.py", args->out_prefix);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "\n"
            "dat = [\n"
           );
    int i, nbins = dists[0].nbins;
    double *good_data = (double*) malloc(sizeof(double)*nbins);
    double *all_data  = (double*) malloc(sizeof(double)*nbins);
    double *bad_data  = (double*) malloc(sizeof(double)*nbins);
    for (i=0; i<ndist; i++)
    {
        smooth_dist(good_data, dists[i].good_data, nbins);
        smooth_dist(all_data, dists[i].all_data, nbins);
        smooth_dist(bad_data, dists[i].bad_data, nbins);

        int j, jfrom, jto;
        for (jfrom=0; jfrom<nbins-1 && all_data[jfrom]==0; jfrom++) if ( all_data[jfrom+1]!=0 ) break;
        for (jto=nbins-1; jto>0 && all_data[jto]==0; jto--) if ( all_data[jto-1]!=0 ) break;
        if ( jto < jfrom ) jto = jfrom = 0;
        fprintf(fp,"[ # %s - distribution of all sites\n", args->colnames[i+NFIXED]);
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", all_data[j]);
        }
        fprintf(fp,"],\n[");
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", (float)j/(nbins-1)*(dists[i].all_max - dists[i].all_min) + dists[i].all_min);
        }
        fprintf(fp,"],\n[ # good sites\n");

        for (jfrom=0; jfrom<nbins-1 && good_data[jfrom]==0; jfrom++) if ( good_data[jfrom+1]!=0 ) break;
        for (jto=nbins-1; jto>0 && good_data[jto]==0; jto--) if ( good_data[jto-1]!=0 ) break;
        if ( jto < jfrom ) jto = jfrom = 0;
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", good_data[j]);
        }
        fprintf(fp,"],\n[");
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", (float)j/(nbins-1)*(dists[i].all_max - dists[i].all_min) + dists[i].all_min);
        }
        fprintf(fp,"],\n[ # bad sites\n");

        for (jfrom=0; jfrom<nbins-1 && bad_data[jfrom]==0; jfrom++) if ( bad_data[jfrom+1]!=0 ) break;
        for (jto=nbins-1; jto>0 && bad_data[jto]==0; jto--) if ( bad_data[jto-1]!=0 ) break;
        if ( jto < jfrom ) jto = jfrom = 0;
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", bad_data[j]);
        }
        fprintf(fp,"],\n[");
        for (j=jfrom; j<=jto; j++)
        {
            if ( j>jfrom ) fprintf(fp,",");
            fprintf(fp,"%e", (float)j/(nbins-1)*(dists[i].all_max - dists[i].all_min) + dists[i].all_min);
        }
        fprintf(fp,"],\n");
    }
    free(good_data);
    free(all_data);
    fprintf(fp, 
        "]\n"
        "\n"
        "for i in range(len(dat)):\n"
        "   if i%%6==1 or i%%6==3 or i%%6==5: continue\n"
        "   m = max(dat[i])\n"
        "   if m==0: continue\n"
        "   for j in range(len(dat[i])):\n"
        "       dat[i][j] = 1.0*dat[i][j]/m\n"
        );
    for (i=0; i<ndist; i++)
    {
        fprintf(fp, 
            "\n"
            "fig = plt.figure()\n"
            "fig, (ax1, ax2) = plt.subplots(2, 1)\n"
            "ax1.plot(dat[%d*6+1],dat[%d*6],   'rv-', mec='r', ms=3.0, label='All')\n"
            "ax1.plot(dat[%d*6+3],dat[%d*6+2], 'g^-', mec='g', ms=3.0, label='Good')\n"
            "ax1.plot(dat[%d*6+5],dat[%d*6+4], 'ko-', mec='k', ms=3.0, label='Bad')\n"
            "ax2.plot(dat[%d*6+1],dat[%d*6],   'rv-', mec='r', ms=3.0)\n"
            "ax2.plot(dat[%d*6+3],dat[%d*6+2], 'g^-', mec='g', ms=3.0)\n"
            "ax2.plot(dat[%d*6+5],dat[%d*6+4], 'ko-', mec='k', ms=3.0)\n"
            "ax1.set_title('%s distribution')\n"
            "ax1.set_ylim([-0.1,1.1])\n"
            "ax2.set_ylim([-0.1,1.1])\n"
            "xmin=%e; xmax=%e; xd=(xmax-xmin)*0.05\n"
            "ax1.set_xlim([xmin-xd,xmax+xd])\n"
            "xmin=%e; xmax=%e; xd=(xmax-xmin)*0.05\n"
            "ax2.set_xlim([xmin-xd,xmax+xd])\n"
            "ax1.legend(loc='best',frameon=False,numpoints=1)\n"
            "ax2.annotate('(Cropped tails)', xy=(0.6,0.8), xycoords='axes fraction')\n"
            "plt.savefig('%s.dists.%s.png')\n"
            "plt.close()\n"
            "\n\n", i,i,i,i,i,i,i,i,i,i,i,i,args->colnames[i+NFIXED], 
                dists[i].all_min,dists[i].all_max, 
                dists[i].good_min,dists[i].good_max, 
                args->out_prefix, args->colnames[i+NFIXED]
           );
    }
    fclose(fp);
    if ( args->do_plot ) py_plot(fname);
    free(fname);

}
static void init_filters(args_t *args, filters_t *filts, char *str, int scale);
static void destroy_filters(filters_t *filt);
static int pass_filters(filters_t *filt, double *vec);

static void create_dists(args_t *args)
{
    fprintf(stderr,"Sorting annotations and creating distribution stats...\n");

    // Create stats for all annotations in args->file. Sort each annotation using unix sort
    //  and determine percentiles.

    // Hackety hack: annots_reader_next will subset the annotations and return them in
    //  requested order. This is the only place where we want to read and return all the values
    //  as they are. Change the scale, ignore, and col2names arrays temporarily.
    //
    int *ignore_ori    = args->ignore;
    int *col2names_ori = args->col2names;
    args->ignore       = (int*) calloc(args->ncols, sizeof(int));
    args->col2names    = (int*) calloc(args->ncols, sizeof(int));
    args->missing      = (int*) calloc(args->ncols, sizeof(int));

    filters_t *tmp_flt = NULL;
    if ( args->fixed_filters ) 
    {
        tmp_flt = (filters_t *) calloc(1, sizeof(filters_t));
        init_filters(args, tmp_flt, args->fixed_filters, 0);
    }

    int i, nann = args->ncols - NFIXED;
    dist_t *dists = (dist_t*) calloc(nann, sizeof(dist_t));
    FILE **fps = (FILE**) malloc(nann*sizeof(FILE*));
    for (i=0; i<nann; i++)
    {
        args->str.l = 0;
        ksprintf(&args->str, "cat | sort -k1,1g > %s.%s", args->out_prefix, args->colnames[i+NFIXED]);
        fps[i] = popen(args->str.s,"w");
        if ( !fps[i] ) error("%s: %s\n", args->str.s, strerror(errno));
        args->col2names[i + NFIXED] = i;
        dists[i].nbins = 200;
        dists[i].good_data = (unsigned int*) calloc(dists[i].nbins,sizeof(unsigned int));
        dists[i].all_data  = (unsigned int*) calloc(dists[i].nbins,sizeof(unsigned int));
        dists[i].bad_data  = (unsigned int*) calloc(dists[i].nbins,sizeof(unsigned int));
    }

    annots_reader_reset(args);
    while ( annots_reader_next(args) )
    {
        int is_bad = tmp_flt && tmp_flt->nfilt && !pass_filters(tmp_flt, args->raw_vals) ? 1 : 0;
        for (i=0; i<nann; i++)
        {
            if ( args->missing[i] )
            {
                dists[i].nmissing++;
                continue;
            }
            if ( IS_GOOD(args->mask) )
            {
                if ( !dists[i].ngood ) dists[i].good_min = dists[i].good_max = args->raw_vals[i];
                if ( args->raw_vals[i] < dists[i].good_min ) dists[i].good_min = args->raw_vals[i];
                if ( args->raw_vals[i] > dists[i].good_max ) dists[i].good_max = args->raw_vals[i];
                dists[i].ngood++;
            }
            if ( is_bad )
            {
                if ( !dists[i].nbad ) dists[i].bad_min = dists[i].bad_max = args->raw_vals[i];
                if ( args->raw_vals[i] < dists[i].bad_min ) dists[i].bad_min = args->raw_vals[i];
                if ( args->raw_vals[i] > dists[i].bad_max ) dists[i].bad_max = args->raw_vals[i];
                dists[i].nbad++;
            }
            if ( !dists[i].nall ) dists[i].all_min = dists[i].all_max = args->raw_vals[i];
            if ( args->raw_vals[i] < dists[i].all_min ) dists[i].all_min = args->raw_vals[i];
            if ( args->raw_vals[i] > dists[i].all_max ) dists[i].all_max = args->raw_vals[i];
            dists[i].nall++;
            fprintf(fps[i], "%le\t%d\t%d\n", args->raw_vals[i], IS_GOOD(args->mask)?1:0, is_bad);
        }
    }
    // Change the arrays back and clean
    free(args->ignore);
    free(args->col2names);
    args->ignore    = ignore_ori;
    args->col2names = col2names_ori;
    for (i=0; i<nann; i++)
        if ( pclose(fps[i]) ) error("An error occured while processing %s.%s\n", args->out_prefix, args->colnames[i+NFIXED]);
    free(fps);
    if ( tmp_flt )
    {
        destroy_filters(tmp_flt);
        free(tmp_flt);
    }

    // Find the extremes
    FILE *fp;
    for (i=0; i<nann; i++)
    {
        char *fname;
        fp = open_file(&fname,"r","%s.%s", args->out_prefix, args->colnames[i+NFIXED]);
        if ( !fp ) error("Cannot not read %s.%s: %s\n", args->out_prefix, args->colnames[i+NFIXED], strerror(errno));
        int count = 0, is_good, is_bad;
        double val;
        dists[i].scale_min = dists[i].scale_max = HUGE_VAL;
        while ( fscanf(fp,"%le\t%d\t%d",&val,&is_good,&is_bad)==3 )
        {
            count++;
            if ( dists[i].scale_min==HUGE_VAL || (double)100.*count/dists[i].nall < args->lo_pctl ) dists[i].scale_min = val;
            if ( dists[i].scale_max==HUGE_VAL && (double)100.*count/dists[i].nall > args->hi_pctl ) dists[i].scale_max = val;

            // For the distribution plots
            int ibin = (float)(dists[i].nbins - 1)*(val - dists[i].all_min)/(dists[i].all_max - dists[i].all_min);
            if ( is_good ) dists[i].good_data[ibin]++;
            if ( is_bad ) dists[i].bad_data[ibin]++;
            dists[i].all_data[ibin]++;
        }
        if ( dists[i].scale_max==HUGE_VAL ) dists[i].scale_max = val;
        if ( fclose(fp) ) error("An error occurred while processing %s.%s\n", args->out_prefix, args->colnames[i+NFIXED]);
        unlink(fname);
        free(fname);
    }
    plot_dists(args, dists, nann);

    fp = open_file(NULL,"w","%s.n", args->out_prefix);
    fprintf(fp, "# [1]nAll\t[2]nGood\t[3]nBad\t[4]nMissing\t[5]minGood\t[6]maxGood\t[7]minBad\t[8]maxBad\t[9]minAll\t[10]maxAll\t[11]%f percentile\t[12]%f percentile\t[13]Annotation\n", args->lo_pctl, args->hi_pctl);
    for (i=0; i<nann; i++)
    {
        fprintf(fp, "%u\t%u\t%u\t%u\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%s\n", dists[i].nall, dists[i].ngood, dists[i].nbad, dists[i].nmissing, 
            dists[i].good_min, dists[i].good_max, dists[i].bad_min, dists[i].bad_max, dists[i].all_min, dists[i].all_max, 
            dists[i].scale_min, dists[i].scale_max, args->colnames[i+NFIXED]);
        free(dists[i].all_data);
        free(dists[i].good_data);
        free(dists[i].bad_data);
    }
    fclose(fp);
    free(dists);
}

static void init_dists(args_t *args)
{
    FILE *fp = open_file(NULL,"r","%s.n", args->out_prefix);
    if ( !fp ) fp = open_file(NULL,"r","%s.n", args->fname);
    if ( !fp ) 
    {
        create_dists(args);
        fp = open_file(NULL,"r","%s.n", args->out_prefix);
    }
    if ( !fp ) error("Could not read %s.n nor %s.n\n", args->out_prefix,args->fname);

    args->dists = (dist_t*) calloc(args->ncols, sizeof(dist_t));
    args->str.l = 0;
    ks_getline(fp, &args->str); // the header
    int i;
    for (i=0; i<args->ncols - NFIXED; i++)
    {
        args->str.l = 0;
        ks_getline(fp, &args->str);
        char *p = args->str.s;
        int j = 0;
        while ( *p && j<12 ) 
        {
            if ( *p=='\t' ) j++;
            p++;
        }
        if ( !*p ) error("Could not parse the line, expected 13 fields: [%s] [%s]\n", args->str.s, p);
        for (j=NFIXED; j<args->ncols; j++)
            if ( !strcmp(args->colnames[j],p) ) break;
        if ( j==args->ncols ) continue;
        sscanf(args->str.s, "%u\t%u\t%u\t%u\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le", &args->dists[j].nall, &args->dists[j].ngood, &args->dists[j].nbad, &args->dists[j].nmissing,
            &args->dists[j].good_min, &args->dists[j].good_max, &args->dists[j].bad_min, &args->dists[j].bad_max, &args->dists[j].all_min, &args->dists[j].all_max, 
            &args->dists[j].scale_min, &args->dists[j].scale_max);
        if ( args->dists[j].scale_min==args->dists[j].scale_max ) error("The annotation %s is not looking good, please leave it out\n", args->colnames[j]);
    }
    for (i=NFIXED; i<args->ncols; i++)
        if ( !args->dists[i].nall && !args->dists[i].nmissing ) error("No extremes found for the annotation: %s\n", args->colnames[i]);
    fclose(fp);
}

static void init_extra_annot(args_t *args, char *annot)
{
    int i;
    for (i=NFIXED; i<args->ncols; i++)
        if ( !strcmp(annot,args->colnames[i]) ) break;
    if ( i==args->ncols ) error("The annotation \"%s\" is not available.\n", annot);
    args->names = (char**) realloc(args->names, sizeof(char*)*(args->nann+1));
    args->names[ args->nann ] = strdup(annot);
    args->ignore[i] = 0;
    args->col2names[i] = args->nann;
    args->ann2cols[ args->nann  ] = i;
    args->nann++;
}
static void init_annots(args_t *args)
{
    args->file = hts_open(args->fname, "r", NULL);
    if ( !args->file ) error("Could not read %s\n", args->fname);
    hts_getline(args->file,'\n',&args->str);
    const char *exp = "# [1]CHROM\t[2]POS\t[3]MASK\t[4]IS_TS\t[5]TYPE";
    if ( args->str.s[0] != '#' ) error("Missing header line in %s, was vcf query run with -H?\n", args->fname);
    if ( strncmp(args->str.s,exp,strlen(exp)) ) error("Version mismatch? %s: %s\n", args->fname, args->str.s);

    int i, j;
    char **colnames = split_list(args->str.s, '\t', &args->ncols);
    if ( args->ncols >= 8*sizeof(int)-1 ) error("Fixme: Too many columns (%d), currently limited by max %d\n", args->ncols, 8*sizeof(int)-1);

    args->colnames  = (char**) malloc(sizeof(char*)*args->ncols);
    for (i=0; i<args->ncols; i++)
    {
        char *p = colnames[i];
        while (*p && *p!=']') p++;
        assert(*p);
        args->colnames[i] = strdup(p+1);
        free(colnames[i]);
    }
    free(colnames);
    // check that column names are unique
    for (i=0; i<args->ncols; i++)
        for (j=0; j<i; j++)
            if ( !strcmp(args->colnames[i],args->colnames[j]) ) error("Error: duplicate column names in %s [%s]\n", args->fname, args->colnames[i]);
    args->col2names = (int*) malloc(sizeof(int)*args->ncols);
    args->missing   = (int*) malloc(sizeof(int)*args->ncols);
    args->ignore    = (int*) malloc(sizeof(int)*args->ncols);
    args->vals      = (double*) malloc(sizeof(double)*args->ncols);
    args->raw_vals  = (double*) malloc(sizeof(double)*args->ncols);
    args->ann2cols  = (int*) malloc(sizeof(int)*args->ncols);
    for (i=0; i<args->ncols; i++)
    {
        args->col2names[i] = -1;
        args->ignore[i] = 1;
    }
    for (i=0; i<NFIXED; i++)
        args->ignore[i] = 0;

    if ( !args->annot_str )
    {
        // all annotations
        args->nann     = args->ncols - NFIXED;
        args->nann_som = args->nann;
        args->names    = (char**) malloc(sizeof(char*)*args->nann_som);
        for (i=NFIXED; i<args->ncols; i++)
        {
            args->col2names[i] = i-NFIXED;
            args->ignore[i] = 0;
            args->names[i-NFIXED] = strdup(args->colnames[i]);
            args->ann2cols[i-NFIXED] = i;
        }
        init_dists(args);
        return;
    }

    args->names     = split_list(args->annot_str, ',', &args->nann);
    args->nann_som  = args->nann;
    for (i=0; i<args->nann; i++)
    {
        for (j=NFIXED; j<args->ncols; j++)
            if ( !strcmp(args->colnames[j], args->names[i]) ) break;
        if ( j==args->ncols ) 
            error("The requested annotation \"%s\" not in %s\n", args->names[i], args->fname);
        if ( args->col2names[j]!=-1 )
            error("The annotation \"%s\" given multiple times?\n", args->names[i]);
        args->col2names[j] = i;
        args->ann2cols[i] = j;
        args->ignore[j] = 0;
    }
    init_dists(args);
}
static void destroy_annots(args_t *args)
{
    free(args->dists);
    if ( args->file ) hts_close(args->file);
    if ( args->str.m ) free(args->str.s);
    if ( args->nann ) 
    {
        destroy_list(args->names, args->nann);
        free(args->col2names);
        free(args->ann2cols);
    }
    int i;
    for (i=0; i<args->ncols; i++)
        free(args->colnames[i]);
    free(args->colnames);
    if ( args->missing ) free(args->missing);
    if ( args->ignore ) free(args->ignore);
    if ( args->vals ) free(args->vals);
    if ( args->raw_vals ) free(args->raw_vals);
}

static int pass_filters(filters_t *filt, double *vec)
{
    int i, j;
    for (i=0; i<filt->n; i++)
    {
        if ( !filt->nfilt[i] ) continue;
        for (j=0; j<filt->nfilt[i]; j++)
        {
            switch (filt->filt[i][j].type) 
            {
                case FLT_BE: if ( vec[i] < filt->filt[i][j].value ) return 0; break;
                case FLT_BT: if ( vec[i] <= filt->filt[i][j].value ) return 0; break;
                case FLT_EQ: if ( vec[i] != filt->filt[i][j].value ) return 0; break;
                case FLT_LT: if ( vec[i] >= filt->filt[i][j].value ) return 0; break;
                case FLT_LE: if ( vec[i] > filt->filt[i][j].value ) return 0; break;
            }
        }
    }
    return 1;
}
static void init_filters(args_t *args, filters_t *filts, char *str, int scale)
{
    filts->n     = args->nann;
    filts->filt  = (filter_t**) calloc(args->ncols-NFIXED,sizeof(filter_t*));
    filts->nfilt = (int*) calloc(args->ncols-NFIXED,sizeof(int));

    char *s, *e = str;
    args->str.l = 0;
    while ( *e )
    {
        if ( !isspace(*e) ) kputc(*e, &args->str);
        e++;
    }

    kstring_t left = {0,0,0}, right = {0,0,0};
    int type = 0;

    e = s = args->str.s;
    while ( *e )
    {
        if ( *e=='&' ) { s = ++e; continue; }
        if ( *e=='<' || *e=='>' || *e=='=' )
        {
            left.l = 0;
            kputsn(s, e-s, &left);
            s = e;
            while ( *e && (*e=='<' || *e=='>' || *e=='=') ) e++;
            if ( !*e ) error("Could not parse filter expression: %s\n", str);
            if ( e-s==2 )
            {
                if ( !strncmp(s,"==",2) ) type = FLT_EQ;
                else if ( !strncmp(s,"<=",2) ) type = FLT_LE;
                else if ( !strncmp(s,">=",2) ) type = FLT_BE;
                else error("Could not parse filter expression: %s\n", str);
            }
            else if ( e-s==1 )
            {
                switch (*s)
                {
                    case '>': type = FLT_BT; break;
                    case '<': type = FLT_LT; break;
                    case '=': type = FLT_EQ; break;
                    default: error("Could not parse filter expression: %s\n", str); break;
                }
            }
            else error("Could not parse filter expression: %s\n", str);
            s = e;
            while ( *e && *e!='&' ) e++;
            right.l = 0;
            kputsn(s, e-s, &right);

            char *ann = NULL;
            int i;
            for (i=NFIXED; i<args->ncols; i++)
            {
                if ( !strcmp(args->colnames[i],left.s) ) { s = right.s; ann = left.s; break; }
                if ( !strcmp(args->colnames[i],right.s) ) { s = left.s; ann = right.s; type *= -1; break; }
            }
            if ( i==args->ncols ) error("No such annotation is available: %s\n", str);
            if ( args->col2names[i]==-1 )
            {
                init_extra_annot(args, ann);
                filts->n = args->nann;
            }
            i = args->col2names[i];
            filts->nfilt[i]++;
            filts->filt[i] = (filter_t*) realloc(filts->filt[i],sizeof(filter_t)*filts->nfilt[i]);
            filter_t *filt = &filts->filt[i][filts->nfilt[i]-1];
            filt->type = type;
            if ( sscanf(s,"%le",&filt->value)!=1 ) error("Could not parse filter expression: %s\n", args->fixed_filters);
            if ( scale )
                filt->value = scale_value(&args->dists[args->ann2cols[i]], filt->value);
            s = e;
            continue;
        }
        e++;
    }
    if ( left.m ) free(left.s);
    if ( right.m ) free(right.s);
}
static void destroy_filters(filters_t *filt)
{
    if ( !filt->nfilt ) return;
    int i;
    for (i=0; i<filt->n; i++)
        if ( filt->filt[i] ) free(filt->filt[i]);
    free(filt->filt);
    free(filt->nfilt);
}

static void som_plot(som_t *som, char *bname, int do_plot)
{
    char *fname;
    FILE *fp = open_file(&fname,"w","%s.py", bname);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "\n"
            "dat = [\n"
           );
    int i,j;
    double *val = som->c;
    for (i=0; i<som->nbin; i++)
    {
        fprintf(fp,"[");
        for (j=0; j<som->nbin; j++)
        {
            if ( j>0 ) fprintf(fp,",");
            fprintf(fp,"%e", *val);
            val++;
        }
        fprintf(fp,"],\n");
    }
    fprintf(fp,
            "]\n"
            "fig = plt.figure()\n"
            "ax1 = plt.subplot(111)\n"
            "im1 = ax1.imshow(dat, interpolation='nearest')\n"
            "fig.colorbar(im1)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", bname
           );
    fclose(fp);
    if ( do_plot ) py_plot(fname);
    free(fname);
}
static void som_train(som_t *som, double *vec)
{
    // find the best matching unit: a node with minimum distance from the input vector
    double *ptr = som->w, *cnt = som->c;
    double min_dist = HUGE_VAL;
    int i, j, k, imin = 0, jmin = 0;
    for (i=0; i<som->nbin; i++)
    {
        for (j=0; j<som->nbin; j++)
        {
            double dist = 0;
            for (k=0; k<som->kdim; k++)
                dist += (vec[k] - ptr[k]) * (vec[k] - ptr[k]);
            if ( dist < min_dist )
            {
                min_dist = dist;
                imin = i;
                jmin = j;
            }
            ptr += som->kdim;
        }
    }

    // calculate the radius
    som->t++;
    double radius = som->radius * exp(-som->t/som->decay);
    radius *= radius;

    double learning_rate = som->learn * exp(-som->t/som->decay);

    // update the weights
    ptr = som->w;
    for (i=0; i<som->nbin; i++)
    {
        for (j=0; j<som->nbin; j++)
        {
            double dist = (i-imin)*(i-imin) + (j-jmin)*(j-jmin);
            if ( dist <= radius ) 
            {
                double influence = exp(-dist*dist*0.5/radius) * learning_rate;
                for (k=0; k<som->kdim; k++)
                    ptr[k] += influence * (vec[k] - ptr[k]);
                *cnt += influence;
            }
            ptr += som->kdim;
            cnt++;
        }
    }
}
static double som_calc_dist(som_t *som, double *vec)
{
    double min_dist = HUGE_VAL, *cnt = som->c, *ptr = som->w;
    int i, j, k;
    for (i=0; i<som->nbin; i++)
    {
        for (j=0; j<som->nbin; j++)
        {
            if ( *cnt >= som->th )
            {
                double dist = 0;
                for (k=0; k<som->kdim; k++)
                    dist += (vec[k] - ptr[k]) * (vec[k] - ptr[k]);

                if ( dist < min_dist )
                    min_dist = dist;
            }
            cnt++;
            ptr += som->kdim;
        }
    }
    return min_dist;
}
static void som_norm(som_t *som)
{
    int i;
    double max = 0;
    for (i=0; i<som->nbin*som->nbin; i++) 
        if ( max < som->c[i] ) max = som->c[i];
    for (i=0; i<som->nbin*som->nbin; i++) som->c[i] /= max;
}
static void som_destroy(som_t *som)
{
    free(som->w);
    free(som->c);
}
static void som_init(args_t *args)
{
    int i, j;
    som_t *som  = &args->som;
    som->kdim   = args->nann_som;
    som->radius = som->nbin / 2;
    som->w      = (double*) malloc(sizeof(double) * som->kdim * som->nbin * som->nbin);
    som->c      = (double*) calloc(som->nbin*som->nbin, sizeof(double));
    int n = INT_MAX;
    for (i=0; i<args->nann_som; i++)
        if ( args->dists[args->ann2cols[i]].ngood < n ) n = args->dists[args->ann2cols[i]].ngood;
    if ( !som->nt || som->nt > n ) som->nt = n;
    som->decay  = som->nt / log(som->radius);
    srandom(args->rand_seed);
    int n3 = som->nbin * som->nbin * som->kdim;
    for (i=0; i<n3; i++)
        som->w[i] = (double)random()/RAND_MAX;
    int nvals = 0, nbad_vals = 0;
    double *vals = (double*) malloc(sizeof(double)*som->nt*args->nann_som);
    double *bad_vals = (double*) malloc(sizeof(double)*som->nt*args->nann_som);


    // bad sites' SOM
    som_t *bad_som  = &args->bad_som;
    *bad_som = *som;
    bad_som->w      = (double*) malloc(sizeof(double) * som->kdim * som->nbin * som->nbin);
    bad_som->c      = (double*) calloc(som->nbin*som->nbin, sizeof(double));
    for (i=0; i<n3; i++)
        bad_som->w[i] = som->w[i];


    srandom(args->rand_seed);
    annots_reader_reset(args);
    while ( annots_reader_next(args) )
    {
        if ( args->type != args->filt_type ) continue;

        // Note that filtering and training requires that all annotations used in
        //  hard-filtering and learning are present, even when not used in SOM filtering. 
        //  So far the decision was that it is the responsibility of the caller
        //  to supply meaningful default values if the annotation is not
        //  applicable for a site (for example some tests cannot be performed
        //  on a pure HOM site). This may not be an optimal solution and it
        //  might be better to treat the missing values here, possibly with a
        //  user switch to decide whether the missing values by default pass or
        //  fail hard filters.
        //
        if ( args->nset != args->nann ) continue;
        if ( args->filt_excl.nfilt && !pass_filters(&args->filt_excl, args->raw_vals) ) 
        {
            i = nbad_vals < bad_som->nt ? nbad_vals : (double)(bad_som->nt-1)*random()/RAND_MAX;
            assert( i>=0 && i<bad_som->nt );

            for (j=0; j<args->nann_som; j++)
                bad_vals[i*args->nann_som + j] = args->vals[j];

            nbad_vals++;
            continue;
        }
        if ( !IS_GOOD(args->mask) )
        {
            if ( !args->filt_learn.nfilt ) continue;    // this is supervised learning, ignore non-training sites
            if ( !pass_filters(&args->filt_learn, args->vals) ) continue;
        }
        i = nvals < som->nt ? nvals : (double)(som->nt-1)*random()/RAND_MAX;
        assert( i>=0 && i<som->nt );

        for (j=0; j<args->nann_som; j++)
            vals[i*args->nann_som + j] = args->vals[j];

        nvals++;
    }
    if ( nvals < som->nt ) som->nt = nvals;
    if ( nbad_vals < bad_som->nt ) bad_som->nt = nbad_vals;
    fprintf(stderr,"Selected %d training values for good and %d values for bad SOM\n", som->nt, bad_som->nt);

    for (i=0; i<som->nt; i++)
        som_train(som, &vals[i*args->nann_som]);
    free(vals);

    for (i=0; i<bad_som->nt; i++)
        som_train(bad_som, &bad_vals[i*args->nann_som]);
    free(bad_vals);

    som_norm(som);
    som_norm(bad_som);

    char *fname = msprintf("%s.som", args->out_prefix);
    som_plot(som, fname, args->do_plot);
    free(fname);
    fname = msprintf("%s.bad_som", args->out_prefix);
    som_plot(bad_som, fname, args->do_plot);
    free(fname);
}

static void init_data(args_t *args)
{
    fprintf(stderr,"Initializing and training...\n");
    if ( !args->out_prefix ) args->out_prefix = strdup(args->fname);
    else
    {
        int len = strlen(args->out_prefix);
        if ( args->out_prefix[len-1] == '/' || args->out_prefix[len-1] == '\\' )
            args->out_prefix = msprintf("%sannots", args->out_prefix);
        else
            args->out_prefix = strdup(args->out_prefix);
    }
    if ( !args->out_prefix ) error("Missing the -o option.\n");
    init_annots(args);
    if ( args->fixed_filters ) init_filters(args, &args->filt_excl, args->fixed_filters, 0);
    if ( args->learning_filters ) init_filters(args, &args->filt_learn, args->learning_filters, 0);
    som_init(args);
}
static void destroy_data(args_t *args)
{
    destroy_filters(&args->filt_excl);
    destroy_filters(&args->filt_learn);
    destroy_annots(args);
    free(args->out_prefix);
    som_destroy(&args->som);
    som_destroy(&args->bad_som);
}

static void eval_filters(args_t *args)
{
    init_data(args);
    char *type = args->filt_type==VCF_SNP ? "SNP" : "INDEL";

    char *fname; 
    open_file(&fname,NULL,"%s.%s.sites.gz", args->out_prefix, type); // creates also directory
    htsFile *file = hts_open(fname, "wb", NULL);

    // Calculate scores
    kstring_t str = {0,0,0};
    kputs("# [1]score\t[2]good SOM\t[3]bad SOM\t[4]is transition\t[5]mask, good(&2) or bad(&4)\t[6]chromsome\t[7]position\n", &str);
    bgzf_write((BGZF *)file->fp, str.s, str.l);
    str.l = 0;

    fprintf(stderr,"Classifying...\n");
    annots_reader_reset(args);
    int ngood = 0, nall = 0;
    double max_dist = args->som.kdim;
    while ( annots_reader_next(args) )
    {
        if ( args->type != args->filt_type ) continue;
        if ( args->nset != args->nann ) continue;

        double dist = max_dist, bad_dist = 0;
        if ( !args->filt_excl.nfilt || pass_filters(&args->filt_excl, args->raw_vals) )
        {
            dist = som_calc_dist(&args->som, args->vals);
            bad_dist = som_calc_dist(&args->bad_som, args->vals);
            assert( dist>=0 && dist<=max_dist );
        }
        else
            args->mask |= MASK_BAD;
        nall++;

        if ( IS_GOOD(args->mask) ) ngood++;

        str.l = 0;
        double score = dist/max_dist + (1-bad_dist/max_dist);        // 2.154    835770  99.48   2.07    1.103583e+00 (mpileup-test-2013-05-10/filter-sf-sum)
        //* double score = dist/max_dist;                            // 2.153    837182  99.48   2.07 (mpileup-test-2013-05-10/filter-hf)
        // double score = bad_dist>0 ? dist / bad_dist : 1e+10;      // 2.153    823529  99.35   2.07
        // double score = dist/max_dist*(1-bad_dist/max_dist);       // 2.153    837196  99.48   2.07
        // double score = 1-bad_dist/max_dist;                       // 2.154    822909  99.35   2.07
        // double score = dist/max_dist + (bad_dist/max_dist<0.1 ? 0.1 : 0);    // 2.153    837198  99.48   2.07    2.037545e-01
        // double score = dist/max_dist * (bad_dist/max_dist<1e-2 ? 2 : 1);     // 2.153    837199  99.48   2.07    2.075090e-01
        // double score = dist/max_dist + 0.5*(1-bad_dist/max_dist);            // 2.154    835747  99.48   2.07    6.035832e-01
        ksprintf(&str, "%le\t%le\t%le\t%d\t%d\t%s\t%d\n", score, dist/max_dist, bad_dist/max_dist, args->is_ts, args->mask, args->chr, args->pos);
        bgzf_write((BGZF *)file->fp, str.s, str.l);
    }
    hts_close(file);
    tbx_conf_t conf = { 0, 6, 7, 0, '#', 0 };
    tbx_index_build(fname,0,&conf);
    free(fname);


    // Evaluate
    fprintf(stderr,"Evaluating...\n");
    int ngood_read = 0, nall_read = 0, ntstv[2] = {0,0}, ntstv_novel[2] = {0,0};
    double prev_tstv = 0;
    args->str.l = 0;
    ksprintf(&args->str, "gunzip -c %s.%s.sites.gz | sort -k1,1g", args->out_prefix, type);
    FILE *fp = popen(args->str.s, "r");
    if ( !fp ) error("Could not run \"%s\": %s\n", args->str.s, strerror(errno));
    FILE *out = open_file(NULL,"w","%s.%s.tab", args->out_prefix, type);
    fprintf(out,"# [1]ts/tv (all)\t[2]nAll\t[3]sensitivity\t[4]ts/tv novel\t[5]threshold\n");
    while ( ks_getline(fp, &args->str) )
    {
        if ( args->str.s[0]=='#' ) 
        {
            args->str.l = 0;
            continue;
        }

        // SCORE
        double dist = atof(args->str.s);
        char *t = column_next(args->str.s, '\t'); 
        if ( !*t ) error("Could not parse score: [%s]\n", args->str.s);
        // GOOD-SOM DIST
        t = column_next(t+1, '\t'); 
        if ( !*t ) error("Could not parse bad score: [%s]\n", args->str.s);
        // BAD-SOM DIST
        t = column_next(t+1, '\t'); 
        if ( !*t ) error("Could not parse bad score: [%s]\n", args->str.s);
        // IS_TS
        int is_ts = strtol(t, NULL, 10);
        t = column_next(t+1, '\t');  
        if ( !*t ) error("Could not parse POS: [%s]\n", args->str.s);
        // IS_GOOD
        int mask = strtol(t, NULL, 10);
        args->str.l = 0;

        nall_read++;
        ntstv[is_ts]++;
        if ( IS_GOOD(mask) ) ngood_read++;
        else if ( ngood ) ntstv_novel[is_ts]++;
        if ( (double)nall_read/nall < 0.1 ) continue;
        double tstv = (double)ntstv[1]/ntstv[0];
        if ( prev_tstv==0 || fabs(tstv-prev_tstv)>0.005 )
        {
            double ntstv2 = ntstv_novel[0] ? (double)ntstv_novel[1]/ntstv_novel[0] : 0;
            fprintf(out, "%.3f\t%d\t%.2f\t%.2f\t%le\n", tstv, nall_read, ngood ? (double)100.*ngood_read/ngood : 0, ntstv2, dist);
            prev_tstv = tstv;
        }
    }
    if ( str.m ) free(str.s);
    fclose(out);
    fclose(fp);
    destroy_data(args);
}

typedef struct
{
    htsFile *file;
    hts_itr_t *itr;
    tbx_t *tbx;
    kstring_t str;
    double score;
    int rid, pos;
}
site_t;
static void destroy_site(site_t *site)
{
    if ( site->str.m ) free(site->str.s); 
    if ( site->file ) hts_close(site->file);
    if ( site->itr ) tbx_itr_destroy(site->itr);
    if ( site->tbx ) tbx_destroy(site->tbx);
    free(site);
}
static site_t *init_site(char *prefix, const char *type, char *region)
{
    site_t *site = (site_t*) calloc(1,sizeof(site_t));
    char *fname = msprintf("%s.%s.sites.gz", prefix, type);
    site->file = hts_open(fname, region ? "rb" : "r", NULL);
    if ( !site->file ) error("Error: could not read %s\n", fname);
    if ( region ) 
    {
        site->tbx = tbx_index_load(fname);
        if ( !site->tbx ) error("Error: could not load index of %s\n", fname);
        site->itr = tbx_itr_querys(site->tbx,region);
        if ( !site->itr ) error("Error: could not init itr of %s\n", fname);
    }
    free(fname);
    return site;
}
static int sync_site(bcf_hdr_t *hdr, bcf1_t *line, site_t *site, int type)
{
    while (1)
    {
        if ( !site->str.l )
        {
            // no site in the buffer
            if ( site->itr )
            {
                if ( tbx_itr_next((BGZF*)site->file->fp, site->tbx, site->itr, &site->str) < 0 ) return 0;
            }
            else
                if ( hts_getline(site->file, '\n', &site->str) <= 0 ) return 0;
            
            if ( site->str.s[0]=='#' ) 
            {
                site->str.l = 0;
                continue;
            }

            // SCORE
            site->score = atof(site->str.s);
            char *t = column_next(site->str.s, '\t'); 
            if ( !*t ) error("Could not parse SCORE: [%s]\n", site->str.s);
            // IS_TS
            t = column_next(t+1, '\t');  
            if ( !*t ) error("Could not parse IS_TS: [%s]\n", site->str.s);
            // IS_GOOD mask
            t = column_next(t+1, '\t');  
            if ( !*t ) error("Could not parse IS_GOOD: [%s]\n", site->str.s);
            // CHR
            char *p = t+1;
            t = column_next(t+1, '\t');  
            *t = 0;
            site->rid = bcf_name2id(hdr, p);
            if ( site->rid<0 ) error("The chrom \"%s\" not in the header?\n", p);
            // POS
            site->pos = strtol(t+1, NULL, 10);
        }
        if ( !(line->d.var_type & type) ) return 0;
        if ( line->pos+1 < site->pos || line->rid != site->rid ) return 0;
        if ( line->pos+1 == site->pos && line->rid == site->rid )
        {
            site->str.l = 0;
            return 1;
        }
        assert( line->pos+1 <= site->pos || line->rid != site->rid );
        break;
    }
    return 0;
}

static void apply_filters(args_t *args)
{
    // Init files
    site_t *snp = NULL, *indel = NULL;
    if ( args->snp_th >= 0 ) snp = init_site(args->out_prefix, "SNP", args->region);
    if ( args->indel_th >= 0 ) indel = init_site(args->out_prefix, "INDEL", args->region);

    bcf_srs_t *sr = bcf_sr_init();
    sr->region = args->region;
    if ( !bcf_sr_add_reader(sr, args->fname) ) error("Failed to open or the file not indexed: %s\n", args->fname);
    bcf_hdr_t *hdr = sr->readers[0].header;

    // Add header lines
    kstring_t str = {0,0,0};
    kputs("##FILTER=<ID=FailSOM,Description=\"Failed SOM filter, lower is better:", &str);
    if ( snp )
    {
        ksprintf(&str, " SNP cutoff %e", args->snp_th);
        if ( indel ) kputc(';', &str);
    }
    if ( indel )
        ksprintf(&str, " INDEL cutoff %e", args->indel_th);
    kputs(".\">", &str);
    bcf_hdr_append(hdr, str.s); free(str.s); str.m = str.l = 0;
    bcf_hdr_append(hdr, "##INFO=<ID=FiltScore,Number=1,Type=Float,Description=\"SOM Filtering Score\">");
    bcf_hdr_fmt_text(hdr);
    int flt_pass = bcf_id2int(hdr, BCF_DT_ID, "PASS");
    int flt_fail = bcf_id2int(hdr, BCF_DT_ID, "FailSOM");

    htsFile *out = hts_open("-","w",0);
    vcf_hdr_write(out, hdr);

    while ( bcf_sr_next_line(sr) ) 
    {
        bcf1_t *line = sr->readers[0].buffer[0];
        bcf_unpack(line, BCF_UN_INFO);
        bcf_set_variant_types(line);
        assert( line->d.n_flt>=1 );
        line->d.n_flt = 1;
        if ( snp && sync_site(hdr, line, snp, VCF_SNP) )
        {
            float tmp = snp->score;
            bcf1_update_info_float(hdr, line, "FiltScore", &tmp, 1);
            line->d.flt[0] = snp->score <= args->snp_th ? flt_pass : flt_fail;
            vcf_write1(out, hdr, line);
            continue;
        }
        if ( indel && sync_site(hdr, line, indel, VCF_INDEL) )
        {
            float tmp = snp->score;
            bcf1_update_info_float(hdr, line, "FiltScore", &tmp, 1);
            line->d.flt[0] = indel->score <= args->indel_th ? flt_pass : flt_fail;
            vcf_write1(out, hdr, line);
            continue;
        }
        line->d.n_flt = 0;
        vcf_write1(out, hdr, line);
    }
    if ( snp ) destroy_site(snp);
    if ( indel ) destroy_site(indel);
    hts_close(out);
    bcf_sr_destroy(sr);
}

static void usage(void)
{
	fprintf(stderr, "About:   SOM (Self-Organizing Map) filtering.\n");
	fprintf(stderr, "Usage:   vcffilter [options] <annots.tab.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -a, --annots <list>                        list of annotations (default: use all annotations)\n");
	fprintf(stderr, "    -f, --fixed-filter <expr>                  list of fixed threshold filters to apply (absolute values, e.g. 'QUAL>4')\n");
	fprintf(stderr, "    -i, --indel-threshold <float>              filter INDELs\n");
	fprintf(stderr, "    -l, --learning-filters <expr>              filters for selecting training sites (values scaled to interval [0-1])\n");
	fprintf(stderr, "    -m, --map-parameters <int,float,float>     number of bins, learning constant, unit threshold [20,0.1,0.2]\n");
	fprintf(stderr, "    -n, --ntrain-sites <int>                   number of training sites to use\n");
	fprintf(stderr, "    -o, --output-prefix <string>               prefix of output files\n");
	fprintf(stderr, "    -p, --create-plots                         create plots\n");
	fprintf(stderr, "    -r, --region <chr|chr:from-to>             apply filtering in this region only\n");
	fprintf(stderr, "    -R, --random-seed <int>                    random seed, 0 for time() [1]\n");
	fprintf(stderr, "    -s, --snp-threshold <float>                filter SNPs\n");
	fprintf(stderr, "    -t, --type <SNP|INDEL>                     variant type to filter [SNP]\n");
	fprintf(stderr, "Example:\n");
	fprintf(stderr, "   # 1) Extract annotations from the VCF. This is because several passes through the data are\n");
    fprintf(stderr, "   # required and VCF parsing is slow. The second VCF is required only for supervised learning.\n");
	fprintf(stderr, "   vcf query -Ha QUAL,Annot1,Annot2,... target.vcf.gz | bgzip -c > annots.tab.gz\n");
	fprintf(stderr, "   vcf query -Ha QUAL,Annot1,Annot2,... target.vcf.gz training.vcf.gz | bgzip -c > annots.tab.gz\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "   # 2) Test which annotations and parameters give the best result. SNPs and INDELs are done\n");
    fprintf(stderr, "   # separately. Notice the use of -l for unsupervised learning in the second example.\n");
    fprintf(stderr, "   # Without the -a option, all annotations from annots.tab.gz are used.\n");
	fprintf(stderr, "   vcf filter annots.tab.gz -o prefix -p -f'QUAL>4' -a Annot2,Annot3\n");
	fprintf(stderr, "   vcf filter annots.tab.gz -o prefix -p -f'QUAL>4' -l'QUAL>0.6'\n");
	fprintf(stderr, "   vcf filter annots.tab.gz -o prefix -p -f'QUAL>4' -l'QUAL>0.6' -t INDEL\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "   # 3) Choose threshold in prefix.SNP.tab and prefix.INDEL.tab and apply with -i and -s options\n");
	fprintf(stderr, "   vcf filter target.vcf.gz -o prefix -s 1.054277e-02 | bgzip -c > filtered.vcf.gz\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcffilter(int argc, char *argv[])
{
	int c;
	args_t *args    = (args_t*) calloc(1,sizeof(args_t));
	args->argc      = argc; args->argv = argv;
    args->lo_pctl   = 0.1;
    args->hi_pctl   = 99.9;
    args->som.nbin  = 20;
    args->som.learn = 0.1;
    args->som.th    = 0.2;
    args->scale     = 1;
    args->filt_type = VCF_SNP;
    args->snp_th    = -1;
    args->indel_th  = -1;
    args->rand_seed = 1;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"annots",1,0,'a'},
		{"output-prefix",1,0,'o'},
		{"map-params",1,0,'m'},
		{"create-plots",0,0,'p'},
		{"fixed-filter",1,0,'f'},
		{"do-not-scale",0,0,'S'},
		{"type",1,0,'t'},
		{"ntrain-sites",1,0,'n'},
		{"learning-filters",1,0,'l'},
		{"snp-threshold",1,0,'s'},
		{"indel-threshold",1,0,'i'},
		{"random-seed",1,0,'R'},
		{"region",1,0,'r'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "ho:a:s:i:pf:St:n:l:m:r:R:",loptions,NULL)) >= 0) {
		switch (c) {
			case 'S': args->scale = 0; break;
			case 't': 
                {
                    if ( !strcasecmp("SNP",optarg) ) args->filt_type = VCF_SNP;  
                    else if ( !strcasecmp("INDEL",optarg) ) args->filt_type = VCF_INDEL;
                    else error("The variant type \"%s\" not recognised.\n", optarg);
                    break;
                }
			case 'a': args->annot_str = optarg; break;
			case 'n': args->som.nt = atoi(optarg); break;
			case 's': args->snp_th = atof(optarg); break;
			case 'i': args->indel_th = atof(optarg); break;
			case 'o': args->out_prefix = optarg; break;
			case 'p': args->do_plot = 1; break;
			case 'R': args->rand_seed = atoi(optarg); break;
			case 'r': args->region = optarg; break;
            case 'f': args->fixed_filters = optarg; break;
            case 'l': args->learning_filters = optarg; break;
			case 'm': if (sscanf(optarg,"%d,%le,%le",&args->som.nbin,&args->som.learn,&args->som.th)!=3) error("Could not parse --SOM-params %s\n", optarg); break;
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
    if ( !args->rand_seed ) args->rand_seed = time(NULL);
    if ( argc!=optind+1 ) usage();
    args->fname = argv[optind];
    fprintf(stderr,"Random seed %d\n", args->rand_seed);
    if ( args->snp_th<0 && args->indel_th<0 )
    {
        if ( args->region ) error("The -r option is to be used with -s or -i only.\n");
        eval_filters(args);
    }
    else
        apply_filters(args);
	free(args);
	return 0;
}

