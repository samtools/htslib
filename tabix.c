/*  tabix.c -- Generic indexer for TAB-delimited genome position files.

    Copyright (C) 2009-2011 Broad Institute.
    Copyright (C) 2010-2012, 2014-2020, 2024 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"
#include "htslib/hts_defs.h"
#include "htslib/hts_log.h"
#include "htslib/thread_pool.h"

//for easy coding
#define RELEASE_TPOOL(X) { hts_tpool *ptr = (hts_tpool*)(X); if (ptr) { hts_tpool_destroy(ptr); } }
#define bam_index_build3(fn, min_shift, nthreads) (sam_index_build3((fn), NULL, (min_shift), (nthreads)))

typedef struct
{
    char *regions_fname, *targets_fname;
    int print_header, header_only, cache_megs, download_index, separate_regs, threads;
}
args_t;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error(const char *format, ...)
{
    va_list ap;
    fflush(stdout);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error_errno(const char *format, ...)
{
    va_list ap;
    int eno = errno;
    fflush(stdout);
    if (format) {
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    if (eno) {
        fprintf(stderr, "%s%s\n", format ? ": " : "", strerror(eno));
    } else {
        fprintf(stderr, "\n");
    }
    fflush(stderr);
    exit(EXIT_FAILURE);
}


#define IS_GFF  (1<<0)
#define IS_BED  (1<<1)
#define IS_SAM  (1<<2)
#define IS_VCF  (1<<3)
#define IS_BCF  (1<<4)
#define IS_BAM  (1<<5)
#define IS_CRAM (1<<6)
#define IS_GAF  (1<<7)
#define IS_TXT  (IS_GFF|IS_BED|IS_SAM|IS_VCF)

int file_type(const char *fname)
{
    int l = strlen(fname);
    if (l>=7 && strcasecmp(fname+l-7, ".gff.gz") == 0) return IS_GFF;
    else if (l>=7 && strcasecmp(fname+l-7, ".bed.gz") == 0) return IS_BED;
    else if (l>=7 && strcasecmp(fname+l-7, ".sam.gz") == 0) return IS_SAM;
    else if (l>=7 && strcasecmp(fname+l-7, ".vcf.gz") == 0) return IS_VCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bcf") == 0) return IS_BCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bam") == 0) return IS_BAM;
    else if (l>=4 && strcasecmp(fname+l-5, ".cram") == 0) return IS_CRAM;
    else if (l>=7 && strcasecmp(fname+l-7, ".gaf.gz") == 0) return IS_GAF;

    htsFile *fp = hts_open(fname,"r");
    if (!fp) {
        if (errno == ENOEXEC) {
            // hts_open() uses this to report that it didn't understand the
            // file format.
            error("Couldn't understand format of \"%s\"\n", fname);
        } else {
            error_errno("Couldn't open \"%s\"", fname);
        }
    }
    enum htsExactFormat format = hts_get_format(fp)->format;
    hts_close(fp);
    if ( format == bcf ) return IS_BCF;
    if ( format == bam ) return IS_BAM;
    if ( format == cram ) return IS_CRAM;
    if ( format == vcf ) return IS_VCF;

    return 0;
}

static char **parse_regions(char *regions_fname, char **argv, int argc, int *nregs)
{
    kstring_t str = {0,0,0};
    int iseq = 0, ireg = 0;
    char **regs = NULL;
    *nregs = argc;

    if ( regions_fname )
    {
        // improve me: this is a too heavy machinery for parsing regions...

        regidx_t *idx = regidx_init(regions_fname, NULL, NULL, 0, NULL);
        if ( !idx ) {
            error_errno("Could not build region list for \"%s\"", regions_fname);
        }
        regitr_t *itr = regitr_init(idx);
        if ( !itr ) {
            error_errno("Could not initialize an iterator over \"%s\"",
                        regions_fname);
        }

        (*nregs) += regidx_nregs(idx);
        regs = (char**) malloc(sizeof(char*)*(*nregs));
        if (!regs) error_errno(NULL);

        int nseq;
        char **seqs = regidx_seq_names(idx, &nseq);
        for (iseq=0; iseq<nseq; iseq++)
        {
            if (regidx_overlap(idx, seqs[iseq], 0, HTS_POS_MAX, itr) < 0)
                error_errno("Failed to build overlapping regions list");

            while ( regitr_overlap(itr) )
            {
                str.l = 0;
                if (ksprintf(&str, "%s:%"PRIhts_pos"-%"PRIhts_pos, seqs[iseq], itr->beg+1, itr->end+1) < 0) {
                    error_errno(NULL);
                }
                regs[ireg] = strdup(str.s);
                if (!regs[ireg]) error_errno(NULL);
                ireg++;
            }
        }
        regidx_destroy(idx);
        regitr_destroy(itr);
    }
    free(str.s);

    if ( !ireg )
    {
        if ( argc )
        {
            regs = (char**) malloc(sizeof(char*)*argc);
            if (!regs) error_errno(NULL);
        }
        else
        {
            regs = (char**) malloc(sizeof(char*));
            if (!regs) error_errno(NULL);
            regs[0] = strdup(".");
            if (!regs[0]) error_errno(NULL);
            *nregs = 1;
        }
    }

    for (iseq=0; iseq<argc; iseq++, ireg++) {
        regs[ireg] = strdup(argv[iseq]);
        if (!regs[ireg]) error_errno(NULL);
    }
    return regs;
}
static int query_regions(args_t *args, tbx_conf_t *conf, char *fname, char **regs, int nregs)
{
    int i;
    htsThreadPool tpool = {NULL, 0};
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) error_errno("Could not open \"%s\"", fname);
    enum htsExactFormat format = hts_get_format(fp)->format;
    if (args->cache_megs)
        hts_set_cache_size(fp, args->cache_megs * 1048576);

    //set threads if needed, errors are logged and ignored
    if (args->threads >= 1) {
        if (!(tpool.pool = hts_tpool_init(args->threads))) {
            hts_log_info("Could not initialize thread pool!");
        }
        if (hts_set_thread_pool(fp, &tpool) < 0) {
            hts_log_info("Could not set thread pool!");
        }
    }

    regidx_t *reg_idx = NULL;
    if ( args->targets_fname )
    {
        reg_idx = regidx_init(args->targets_fname, NULL, NULL, 0, NULL);
        if (!reg_idx) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not build region list for \"%s\"",
                        args->targets_fname);
        }
    }

    if ( format == bcf )
    {
        htsFile *out = hts_open("-","w");
        if ( !out ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not open stdout");
        }
        if (hts_set_thread_pool(out, &tpool) < 0) {
            hts_log_info("Could not set thread pool to output file!");
        }
        hts_idx_t *idx = bcf_index_load3(fname, NULL, args->download_index ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !idx ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not load .csi index of \"%s\"", fname);
        }

        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not read the header from \"%s\"", fname);
        }

        if ( args->print_header ) {
            if ( bcf_hdr_write(out,hdr)!=0 ) {
                RELEASE_TPOOL(tpool.pool);
                error_errno("Failed to write to stdout");
            }
        }
        if ( !args->header_only )
        {
            assert(regs != NULL);
            bcf1_t *rec = bcf_init();
            if (!rec) {
                RELEASE_TPOOL(tpool.pool);
                error_errno(NULL);
            }
            for (i=0; i<nregs; i++)
            {
                int ret, found = 0;
                hts_itr_t *itr = bcf_itr_querys(idx,hdr,regs[i]);
                if (!itr) continue;
                while ((ret = bcf_itr_next(fp, itr, rec)) >=0 )
                {
                    if ( reg_idx )
                    {
                        const char *chr = bcf_seqname(hdr,rec);
                        if (!chr) {
                            RELEASE_TPOOL(tpool.pool);
                            error("Bad BCF record in \"%s\" : "
                                  "Invalid CONTIG id %d\n",
                                  fname, rec->rid);
                        }
                        if ( !regidx_overlap(reg_idx,chr,rec->pos,rec->pos+rec->rlen-1, NULL) ) continue;
                    }
                    if (!found) {
                        if (args->separate_regs) printf("%c%s\n", conf->meta_char, regs[i]);
                        found = 1;
                    }
                    if ( bcf_write(out,hdr,rec)!=0 ) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno("Failed to write to stdout");
                    }
                }

                if (ret < -1) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Reading \"%s\" failed", fname);
                }
                bcf_itr_destroy(itr);
            }
            bcf_destroy(rec);
        }
        if ( hts_close(out) ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("hts_close returned non-zero status for stdout");
        }

        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( format==vcf || format==sam || format==bed || format==text_format || format==unknown_format )
    {
        tbx_t *tbx = tbx_index_load3(fname, NULL, args->download_index ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !tbx ) {
            RELEASE_TPOOL(tpool.pool);
            error_errno("Could not load .tbi/.csi index of %s", fname);
        }
        kstring_t str = {0,0,0};
        if ( args->print_header )
        {
            int ret;
            while ((ret = hts_getline(fp, KS_SEP_LINE, &str)) >= 0)
            {
                if ( !str.l || str.s[0]!=tbx->conf.meta_char ) break;
                if (puts(str.s) < 0) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Error writing to stdout");
                }
            }
            if (ret < -1) {
                RELEASE_TPOOL(tpool.pool);
                error_errno("Reading \"%s\" failed", fname);
            }
        }
        if ( !args->header_only )
        {
            int nseq;
            const char **seq = NULL;
            if ( reg_idx ) {
                seq = tbx_seqnames(tbx, &nseq);
                if (!seq) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Failed to get sequence names list");
                }
            }
            for (i=0; i<nregs; i++)
            {
                int ret, found = 0;
                hts_itr_t *itr = tbx_itr_querys(tbx, regs[i]);
                if ( !itr ) continue;
                while ((ret = tbx_itr_next(fp, tbx, itr, &str)) >= 0)
                {
                    if ( reg_idx && !regidx_overlap(reg_idx,seq[itr->curr_tid],itr->curr_beg,itr->curr_end-1, NULL) ) continue;
                    if (!found) {
                        if (args->separate_regs) printf("%c%s\n", conf->meta_char, regs[i]);
                        found = 1;
                    }
                    if (puts(str.s) < 0) {
                        RELEASE_TPOOL(tpool.pool);
                        error_errno("Failed to write to stdout");
                    }
                }
                if (ret < -1) {
                    RELEASE_TPOOL(tpool.pool);
                    error_errno("Reading \"%s\" failed", fname);
                }
                tbx_itr_destroy(itr);
            }
            free(seq);
        }
        free(str.s);
        tbx_destroy(tbx);
    }
    else if ( format==bam ) {
        RELEASE_TPOOL(tpool.pool);
        error("Please use \"samtools view\" for querying BAM files.\n");
    }

    if ( reg_idx ) regidx_destroy(reg_idx);
    if ( hts_close(fp) ) {
        RELEASE_TPOOL(tpool.pool);
        error_errno("hts_close returned non-zero status: %s", fname);
    }

    for (i=0; i<nregs; i++) free(regs[i]);
    free(regs);
    RELEASE_TPOOL(tpool.pool);
    return 0;
}
static int query_chroms(char *fname, int download)
{
    const char **seq;
    int i, nseq, ftype = file_type(fname);
    if ( ftype & IS_TXT || !ftype )
    {
        tbx_t *tbx = tbx_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !tbx ) error_errno("Could not load .tbi index of %s", fname);
        seq = tbx_seqnames(tbx, &nseq);
        if (!seq) error_errno("Couldn't get list of sequence names");
        for (i=0; i<nseq; i++) {
            if (printf("%s\n", seq[i]) < 0)
                error_errno("Couldn't write to stdout");
        }
        free(seq);
        tbx_destroy(tbx);
    }
    else if ( ftype==IS_BCF )
    {
        htsFile *fp = hts_open(fname,"r");
        if ( !fp ) error_errno("Could not open \"%s\"", fname);
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) error_errno("Could not read the header: \"%s\"", fname);
        hts_close(fp);
        hts_idx_t *idx = bcf_index_load3(fname, NULL, download ? HTS_IDX_SAVE_REMOTE : 0);
        if ( !idx ) error_errno("Could not load .csi index of \"%s\"", fname);
        seq = bcf_index_seqnames(idx, hdr, &nseq);
        if (!seq) error_errno("Couldn't get list of sequence names");
        for (i=0; i<nseq; i++) {
            if (printf("%s\n", seq[i]) < 0)
                error_errno("Couldn't write to stdout");
        }
        free(seq);
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( ftype==IS_BAM )   // todo: BAM
        error("BAM: todo\n");
    return 0;
}

int reheader_file(const char *fname, const char *header, int ftype, tbx_conf_t *conf, int threads)
{
    hts_tpool *tpool = NULL;
    if (threads >= 1) {
        if (!(tpool = hts_tpool_init(threads))) {
            hts_log_info("Could not initialize thread pool!");
        }
    }
    if ( ftype & IS_TXT || !ftype )
    {
        BGZF *fp = bgzf_open(fname,"r");
        if (!fp) {
            RELEASE_TPOOL(tpool);
            return -1;
        }
        if (bgzf_thread_pool(fp, tpool, 0) < 0) {
            hts_log_info("Could not set thread pool!");
        }
        if (bgzf_read_block(fp) != 0 || !fp->block_length ) {
            RELEASE_TPOOL(tpool);
            return -1;
        }

        char *buffer = fp->uncompressed_block;
        int skip_until = 0;

        // Skip the header: find out the position of the data block
        if ( buffer[0]==conf->meta_char )
        {
            skip_until = 1;
            while (1)
            {
                if ( buffer[skip_until]=='\n' )
                {
                    skip_until++;
                    if ( skip_until>=fp->block_length )
                    {
                        if ( bgzf_read_block(fp) != 0 || !fp->block_length ) {
                            RELEASE_TPOOL(tpool);
                            error("FIXME: No body in the file: %s\n", fname);
                        }
                        skip_until = 0;
                    }
                    // The header has finished
                    if ( buffer[skip_until]!=conf->meta_char ) break;
                }
                skip_until++;
                if ( skip_until>=fp->block_length )
                {
                    if (bgzf_read_block(fp) != 0 || !fp->block_length) {
                        RELEASE_TPOOL(tpool);
                        error("FIXME: No body in the file: %s\n", fname);
                    }
                    skip_until = 0;
                }
            }
        }

        // Output the new header
        FILE *hdr  = fopen(header,"r");
        if ( !hdr ) {
            RELEASE_TPOOL(tpool);
            error("%s: %s", header,strerror(errno));
        }
        const size_t page_size = 32768;
        char *buf = malloc(page_size);
        BGZF *bgzf_out = bgzf_open("-", "w");
        ssize_t nread;

        if (!buf) {
            RELEASE_TPOOL(tpool);
            error("%s\n", strerror(errno));
        }
        if (!bgzf_out) {
            RELEASE_TPOOL(tpool);
            error_errno("Couldn't open output stream");
        }
        if (bgzf_thread_pool(bgzf_out, tpool, 0) < 0) {
            hts_log_info("Could not set thread pool to output file!");
        }
        while ( (nread=fread(buf,1,page_size-1,hdr))>0 )
        {
            if ( nread<page_size-1 && buf[nread-1]!='\n' ) buf[nread++] = '\n';
            if (bgzf_write(bgzf_out, buf, nread) < 0) {
                RELEASE_TPOOL(tpool);
                error_errno("Write error %d", bgzf_out->errcode);
            }
        }
        if ( ferror(hdr) ) {
            RELEASE_TPOOL(tpool);
            error_errno("Failed to read \"%s\"", header);
        }
        if ( fclose(hdr) ) {
            RELEASE_TPOOL(tpool);
            error_errno("Closing \"%s\" failed", header);
        }

        // Output all remaining data read with the header block
        if ( fp->block_length - skip_until > 0 )
        {
            if (bgzf_write(bgzf_out, buffer+skip_until, fp->block_length-skip_until) < 0) {
                RELEASE_TPOOL(tpool);
                error_errno("Write error %d",fp->errcode);
            }
        }
        if (bgzf_flush(bgzf_out) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Write error %d", bgzf_out->errcode);
        }

        while (1)
        {
            nread = bgzf_raw_read(fp, buf, page_size);
            if ( nread<=0 ) break;

            int count = bgzf_raw_write(bgzf_out, buf, nread);
            if (count != nread) {
                RELEASE_TPOOL(tpool);
                error_errno("Write failed, wrote %d instead of %d bytes", count,(int)nread);
            }
        }
        if (nread < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error reading \"%s\"", fname);
        }
        if (bgzf_close(bgzf_out) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error %d closing output", bgzf_out->errcode);
        }
        if (bgzf_close(fp) < 0) {
            RELEASE_TPOOL(tpool);
            error_errno("Error %d closing \"%s\"", bgzf_out->errcode, fname);
        }
        free(buf);
    }
    else {
        RELEASE_TPOOL(tpool);
        error("todo: reheader BCF, BAM\n");  // BCF is difficult, records contain pointers to the header.
    }
    RELEASE_TPOOL(tpool);
    return 0;
}

static int usage(FILE *fp, int status)
{
    fprintf(fp, "\n");
    fprintf(fp, "Version: %s\n", hts_version());
    fprintf(fp, "Usage:   tabix [OPTIONS] [FILE] [REGION [...]]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Indexing Options:\n");
    fprintf(fp, "   -0, --zero-based           coordinates are zero-based\n");
    fprintf(fp, "   -b, --begin INT            column number for region start [4]\n");
    fprintf(fp, "   -c, --comment CHAR         skip comment lines starting with CHAR [null]\n");
    fprintf(fp, "   -C, --csi                  generate CSI index for VCF (default is TBI)\n");
    fprintf(fp, "   -e, --end INT              column number for region end (if no end, set INT to -b) [5]\n");
    fprintf(fp, "   -f, --force                overwrite existing index without asking\n");
    fprintf(fp, "   -m, --min-shift INT        set minimal interval size for CSI indices to 2^INT [14]\n");
    fprintf(fp, "   -p, --preset STR           gff, bed, sam, vcf, gaf\n");
    fprintf(fp, "   -s, --sequence INT         column number for sequence names (suppressed by -p) [1]\n");
    fprintf(fp, "   -S, --skip-lines INT       skip first INT lines [0]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Querying and other options:\n");
    fprintf(fp, "   -h, --print-header         print also the header lines\n");
    fprintf(fp, "   -H, --only-header          print only the header lines\n");
    fprintf(fp, "   -l, --list-chroms          list chromosome names\n");
    fprintf(fp, "   -r, --reheader FILE        replace the header with the content of FILE\n");
    fprintf(fp, "   -R, --regions FILE         restrict to regions listed in the file\n");
    fprintf(fp, "   -T, --targets FILE         similar to -R but streams rather than index-jumps\n");
    fprintf(fp, "   -D                         do not download the index file\n");
    fprintf(fp, "       --cache INT            set cache size to INT megabytes (0 disables) [10]\n");
    fprintf(fp, "       --separate-regions     separate the output by corresponding regions\n");
    fprintf(fp, "       --verbosity INT        set verbosity [3]\n");
    fprintf(fp, "   -@, --threads INT          number of additional threads to use [0]\n");
    fprintf(fp, "\n");
    return status;
}

int main(int argc, char *argv[])
{
    int c, detect = 1, min_shift = 0, is_force = 0, list_chroms = 0, do_csi = 0;
    tbx_conf_t conf = tbx_conf_gff;
    char *reheader = NULL;
    args_t args;
    memset(&args,0,sizeof(args_t));
    args.cache_megs = 10;
    args.download_index = 1;
    int32_t new_line_skip = -1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 2},
        {"regions", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 'T'},
        {"csi", no_argument, NULL, 'C'},
        {"zero-based", no_argument, NULL, '0'},
        {"print-header", no_argument, NULL, 'h'},
        {"only-header", no_argument, NULL, 'H'},
        {"begin", required_argument, NULL, 'b'},
        {"comment", required_argument, NULL, 'c'},
        {"end", required_argument, NULL, 'e'},
        {"force", no_argument, NULL, 'f'},
        {"min-shift", required_argument, NULL, 'm'},
        {"preset", required_argument, NULL, 'p'},
        {"sequence", required_argument, NULL, 's'},
        {"skip-lines", required_argument, NULL, 'S'},
        {"list-chroms", no_argument, NULL, 'l'},
        {"reheader", required_argument, NULL, 'r'},
        {"version", no_argument, NULL, 1},
        {"verbosity", required_argument, NULL, 3},
        {"cache", required_argument, NULL, 4},
        {"separate-regions", no_argument, NULL, 5},
        {"threads", required_argument, NULL, '@'},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    while ((c = getopt_long(argc, argv, "hH?0b:c:e:fm:p:s:S:lr:CR:T:D@:", loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'R': args.regions_fname = optarg; break;
            case 'T': args.targets_fname = optarg; break;
            case 'C': do_csi = 1; break;
            case 'r': reheader = optarg; break;
            case 'h': args.print_header = 1; break;
            case 'H': args.print_header = 1; args.header_only = 1; break;
            case 'l': list_chroms = 1; break;
            case '0': conf.preset |= TBX_UCSC; detect = 0; break;
            case 'b':
                conf.bc = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -b %s\n", optarg);
                detect = 0;
                break;
            case 'e':
                conf.ec = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -e %s\n", optarg);
                detect = 0;
                break;
            case 'c': conf.meta_char = *optarg; detect = 0; break;
            case 'f': is_force = 1; break;
            case 'm':
                min_shift = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -m %s\n", optarg);
                break;
            case 'p':
                detect = 0;
                if (strcmp(optarg, "gff") == 0) conf = tbx_conf_gff;
                else if (strcmp(optarg, "bed") == 0) conf = tbx_conf_bed;
                else if (strcmp(optarg, "sam") == 0) conf = tbx_conf_sam;
                else if (strcmp(optarg, "vcf") == 0) conf = tbx_conf_vcf;
                else if (strcmp(optarg, "gaf") == 0) conf = tbx_conf_gaf;
                else if (strcmp(optarg, "bcf") == 0) detect = 1; // bcf is autodetected, preset is not needed
                else if (strcmp(optarg, "bam") == 0) detect = 1; // same as bcf
                else error("The preset string not recognised: '%s'\n", optarg);
                break;
            case 's':
                conf.sc = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -s %s\n", optarg);
                detect = 0;
                break;
            case 'S':
                new_line_skip = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: -S %s\n", optarg);
                detect = 0;
                break;
            case 'D':
                args.download_index = 0;
                break;
            case 1:
                printf(
"tabix (htslib) %s\n"
"Copyright (C) 2024 Genome Research Ltd.\n", hts_version());
                return EXIT_SUCCESS;
            case 2:
                return usage(stdout, EXIT_SUCCESS);
            case 3: {
                int v = atoi(optarg);
                if (v < 0) v = 0;
                hts_set_log_level(v);
                break;
            }
            case 4:
                args.cache_megs = atoi(optarg);
                if (args.cache_megs < 0) {
                    args.cache_megs = 0;
                } else if (args.cache_megs >= INT_MAX / 1048576) {
                    args.cache_megs = INT_MAX / 1048576;
                }
                break;
            case 5:
                args.separate_regs = 1;
                break;
            case '@':   //thread count
                args.threads = atoi(optarg);
                break;
            default: return usage(stderr, EXIT_FAILURE);
        }
    }

    if (new_line_skip >= 0)
        conf.line_skip = new_line_skip;

    if ( optind==argc ) return usage(stderr, EXIT_FAILURE);

    if ( list_chroms )
        return query_chroms(argv[optind], args.download_index);

    char *fname = argv[optind];
    int ftype = file_type(fname);
    if ( detect )  // no preset given
    {
        if ( ftype==IS_GFF ) conf = tbx_conf_gff;
        else if ( ftype==IS_BED ) conf = tbx_conf_bed;
        else if ( ftype==IS_GAF ) conf = tbx_conf_gaf;
        else if ( ftype==IS_SAM ) conf = tbx_conf_sam;
        else if ( ftype==IS_VCF )
        {
            conf = tbx_conf_vcf;
            if ( !min_shift && do_csi ) min_shift = 14;
        }
        else if ( ftype==IS_BCF )
        {
            if ( !min_shift ) min_shift = 14;
        }
        else if ( ftype==IS_BAM )
        {
            if ( !min_shift ) min_shift = 14;
        }
    }
    if ( argc > optind+1 || args.header_only || args.regions_fname || args.targets_fname )
    {
        int nregs = 0;
        char **regs = NULL;
        if ( !args.header_only )
            regs = parse_regions(args.regions_fname, argv+optind+1, argc-optind-1, &nregs);
        return query_regions(&args, &conf, fname, regs, nregs);
    }
    if ( do_csi )
    {
        if ( !min_shift ) min_shift = 14;
        min_shift *= do_csi;  // positive for CSIv2, negative for CSIv1
    }
    if ( min_shift!=0 && !do_csi ) do_csi = 1;

    if ( reheader )
        return reheader_file(fname, reheader, ftype, &conf, args.threads);

    char *suffix = ".tbi";
    if ( do_csi ) suffix = ".csi";
    else if ( ftype==IS_BAM ) suffix = ".bai";
    else if ( ftype==IS_CRAM ) suffix = ".crai";

    char *idx_fname = calloc(strlen(fname) + 6, 1);
    if (!idx_fname) error("%s\n", strerror(errno));
    strcat(strcpy(idx_fname, fname), suffix);

    struct stat stat_tbi, stat_file;
    if ( !is_force && stat(idx_fname, &stat_tbi)==0 )
    {
        // Before complaining about existing index, check if the VCF file isn't
        // newer. This is a common source of errors, people tend not to notice
        // that tabix failed
        stat(fname, &stat_file);
        if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            error("[tabix] the index file exists. Please use '-f' to overwrite.\n");
    }
    free(idx_fname);

    int ret;
    if ( ftype==IS_CRAM )
    {
        if ( bam_index_build3(fname, min_shift, args.threads)!=0 ) error("bam_index_build failed: %s\n", fname);
        return 0;
    }
    else if ( do_csi )
    {
        if ( ftype==IS_BCF )
        {
            if ( bcf_index_build3(fname, NULL, min_shift, args.threads)!=0 ) error("bcf_index_build failed: %s\n", fname);
            return 0;
        }
        if ( ftype==IS_BAM )
        {
            if ( bam_index_build3(fname, min_shift, args.threads)!=0 ) error("bam_index_build failed: %s\n", fname);
            return 0;
        }

        switch (ret = tbx_index_build3(fname, NULL, min_shift, args.threads, &conf))
        {
            case 0:
                return 0;
            case -2:
                error("[tabix] the compression of '%s' is not BGZF\n", fname);
            default:
                error("tbx_index_build3 failed: %s\n", fname);
        }
    }
    else    // TBI index
    {
        switch (ret = tbx_index_build3(fname, NULL, min_shift, args.threads, &conf))
        {
            case 0:
                return 0;
            case -2:
                error("[tabix] the compression of '%s' is not BGZF\n", fname);
            default:
                error("tbx_index_build3 failed: %s\n", fname);
        }
    }

    return 0;
}
