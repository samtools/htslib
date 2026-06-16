/*  synced_bcf_reader.c -- stream through multiple VCF files.

    Copyright (C) 2012-2023, 2025-2026 Genome Research Ltd.

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <inttypes.h>
#include <errno.h>
#include <sys/stat.h>
#include "htslib/synced_bcf_reader.h"
#include "htslib/kseq.h"
#include "htslib/khash_str2int.h"
#include "htslib/bgzf.h"
#include "htslib/hts_alloc.h"
#include "htslib/thread_pool.h"
#include "bcf_sr_sort.h"

#define REQUIRE_IDX_      1
#define ALLOW_NO_IDX_     2

// Maximum indexable coordinate of .csi, for default min_shift of 14.
// This comes out to about 17 Tbp.  Limiting factor is the bin number,
// which is a uint32_t in CSI.  The highest number of levels compatible
// with this is 10 (needs 31 bits).
#define MAX_CSI_COOR ((1LL << (14 + 30)) - 1)

typedef struct
{
    hts_pos_t start, end;   // records are marked for skipping have start>end
}
region1_t;

typedef struct bcf_sr_region_t
{
    region1_t *regs;            // regions will sorted and merged, redundant records marked for skipping have start>end
    int nregs, mregs, creg;     // creg: the current active region
}
region_t;

#define BCF_SR_AUX(x) ((aux_t*)((x)->aux))
typedef struct
{
    sr_sort_t sort;
    int regions_overlap, targets_overlap;
    int *closefile;             // close htsfile with sync reader close or not
}
aux_t;

static bcf_sr_regions_t *bcf_sr_regions_alloc(void);
static int _regions_add(bcf_sr_regions_t *reg, const char *chr, hts_pos_t start, hts_pos_t end);
static bcf_sr_regions_t *_regions_init_string(const char *str);
static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec, bcf_sr_error *errnum);
static void _regions_sort_and_merge(bcf_sr_regions_t *reg);
static int _bcf_sr_regions_overlap(bcf_sr_regions_t *reg, const char *seq, hts_pos_t start, hts_pos_t end, int missed_reg_handler);
static void bcf_sr_seek_start(bcf_srs_t *readers);

char *bcf_sr_strerror(int errnum)
{
    switch (errnum)
    {
        case bcf_sr_open_failed:
            return strerror(errno);
        case not_bgzf:
            return "not compressed with bgzip";
        case idx_load_failed:
            return "could not load index";
        case file_type_error:
            return "unknown file type";
        case api_usage_error:
            return "API usage error";
        case header_error:
            return "could not parse header";
        case no_eof:
            return "no BGZF EOF marker; file may be truncated";
        case no_memory:
            return "Out of memory";
        case vcf_parse_error:
            return "VCF parse error";
        case bcf_read_error:
            return "BCF read error";
        case noidx_error:
            return "merge of unindexed files failed";
        case bcf_sr_seek_error:
            return "failed to seek to next location";
        case bcf_sr_regions_error:
            return "failed to read regions/targets list";
        case bcf_sr_samples_error:
            return "failed to read samples list";
        default: return "";
    }
}

int bcf_sr_set_opt(bcf_srs_t *readers, bcf_sr_opt_t opt, ...)
{
    va_list args;
    switch (opt)
    {
        case BCF_SR_REQUIRE_IDX:
            readers->require_index = REQUIRE_IDX_;
            return 0;

        case BCF_SR_ALLOW_NO_IDX:
            readers->require_index = ALLOW_NO_IDX_;
            return 0;

        case BCF_SR_PAIR_LOGIC:
            va_start(args, opt);
            BCF_SR_AUX(readers)->sort.pair = va_arg(args, int);
            return 0;

        case BCF_SR_REGIONS_OVERLAP:
            va_start(args, opt);
            BCF_SR_AUX(readers)->regions_overlap = va_arg(args, int);
            if ( readers->regions ) readers->regions->overlap = BCF_SR_AUX(readers)->regions_overlap;
            return 0;

        case BCF_SR_TARGETS_OVERLAP:
            va_start(args, opt);
            BCF_SR_AUX(readers)->targets_overlap = va_arg(args, int);
            if ( readers->targets ) readers->targets->overlap = BCF_SR_AUX(readers)->targets_overlap;
            return 0;

        default:
            break;
    }
    return 1;
}

static int *init_filters(bcf_hdr_t *hdr, const char *filters, int *nfilters)
{
    kstring_t str = {0,0,0};
    const char *tmp = filters, *prev = filters;
    int nout = 0, *out = NULL;
    while ( 1 )
    {
        if ( *tmp==',' || !*tmp )
        {
            int *otmp = hts_realloc_ps(out, sizeof(*otmp), nout, 1);
            if (!otmp)
                goto err;
            out = otmp;
            if ( tmp-prev==1 && *prev=='.' )
            {
                out[nout] = -1;
                nout++;
            }
            else
            {
                str.l = 0;
                if (kputsn(prev, tmp-prev, &str) < 0)
                    goto err;
                out[nout] = bcf_hdr_id2int(hdr, BCF_DT_ID, str.s);
                if ( out[nout]>=0 ) nout++;
            }
            if ( !*tmp ) break;
            prev = tmp+1;
        }
        tmp++;
    }
    if ( str.m ) free(str.s);
    *nfilters = nout;
    return out;

 err:
    if (str.m) free(str.s);
    free(out);
    return NULL;
}

int bcf_sr_set_regions(bcf_srs_t *readers, const char *regions, int is_file)
{
    if ( readers->nreaders || readers->regions )
    {
        if ( readers->regions ) bcf_sr_regions_destroy(readers->regions);
        readers->regions = bcf_sr_regions_init(regions,is_file,0,1,-2);
        if (!readers->regions)
        {
            readers->errnum = bcf_sr_regions_error;
            return -1;
        }
        bcf_sr_seek_start(readers);
        return 0;
    }

    readers->regions = bcf_sr_regions_init(regions,is_file,0,1,-2);
    if ( !readers->regions )
    {
        readers->errnum = bcf_sr_regions_error;
        return -1;
    }
    readers->explicit_regs = 1;
    readers->require_index = REQUIRE_IDX_;
    readers->regions->overlap = BCF_SR_AUX(readers)->regions_overlap;
    return 0;
}

int bcf_sr_set_targets(bcf_srs_t *readers, const char *targets, int is_file, int alleles)
{
    if ( readers->nreaders || readers->targets )
    {
        hts_log_error("Must call bcf_sr_set_targets() before bcf_sr_add_reader()");
        readers->errnum = api_usage_error;
        return -1;
    }
    if ( targets[0]=='^' )
    {
        readers->targets_exclude = 1;
        targets++;
    }
    readers->targets = bcf_sr_regions_init(targets,is_file,0,1,-2);
    if ( !readers->targets )
    {
        readers->errnum = bcf_sr_regions_error;
        return -1;
    }
    readers->targets_als = alleles;
    readers->targets->overlap = BCF_SR_AUX(readers)->targets_overlap;
    return 0;
}

int bcf_sr_set_threads(bcf_srs_t *files, int n_threads)
{
    if (!(files->n_threads = n_threads))
        return 0;

    files->p = calloc(1, sizeof(*files->p));
    if (!files->p)
    {
        files->errnum = no_memory;
        return -1;
    }
    if (!(files->p->pool = hts_tpool_init(n_threads)))
    {
        files->errnum = no_memory;
        free(files->p);
        files->p = NULL;
        return -1;
    }

    return 0;
}

void bcf_sr_destroy_threads(bcf_srs_t *files) {
    if (!files->p)
        return;

    if (files->p->pool)
        hts_tpool_destroy(files->p->pool);
    free(files->p);
    files->p = NULL;
}

int bcf_sr_add_reader(bcf_srs_t *files, const char *fname)
{
    char fmode[5];
    int ret = 0;
    const char *idxname = NULL;

    strcpy(fmode, "r");
    vcf_open_mode(fmode+1, fname, NULL);
    htsFile* file_ptr = hts_open(fname, fmode);
    if ( ! file_ptr ) {
        files->errnum = bcf_sr_open_failed;
        return 0;
    }
    //get idx name and pass to add_hreader
    idxname = strstr(fname, HTS_IDX_DELIM);
    if (idxname) idxname += strlen(HTS_IDX_DELIM);
    if (!(ret = bcf_sr_add_hreader(files, file_ptr, 1, idxname))) {
        // bcf_sr_add_hreader() will have set files->errnum
        hts_close(file_ptr);    //failed, close the file
    }
    return ret;
}

int bcf_sr_add_hreader(bcf_srs_t *files, htsFile *file_ptr, int autoclose, const char *idxname)
{
    aux_t *auxdata = BCF_SR_AUX(files);
    assert(auxdata != NULL);
    if ( ! file_ptr ) {
        files->errnum = bcf_sr_open_failed;
        return 0;
    }

    int *new_has_line = hts_realloc_ps(files->has_line, sizeof(*files->has_line),
                                       files->nreaders, 1);
    if (!new_has_line)
        goto nomem;
    files->has_line = new_has_line;
    files->has_line[files->nreaders] = 0;
    bcf_sr_t *new_readers = hts_realloc_ps(files->readers, sizeof(*files->readers),
                                           files->nreaders, 1);
    if (!new_readers)
        goto nomem;
    files->readers = new_readers;
    int *new_closefile = hts_realloc_ps(auxdata->closefile, sizeof(int),
                                        files->nreaders, 1);
    if (!new_closefile)
        goto nomem;
    auxdata->closefile = new_closefile;
    auxdata->closefile[files->nreaders] = 0;  // Will be updated later
    bcf_sr_t *reader = &files->readers[files->nreaders++];
    memset(reader,0,sizeof(bcf_sr_t));

    reader->file = file_ptr;

    files->errnum = 0;

    if ( reader->file->format.compression==bgzf )
    {
        BGZF *bgzf = hts_get_bgzfp(reader->file);
        if ( bgzf && bgzf_check_EOF(bgzf) == 0 ) {
            files->errnum = no_eof;
            hts_log_warning("No BGZF EOF marker; file '%s' may be truncated", file_ptr->fn);
        }
        if (files->p)
            bgzf_thread_pool(bgzf, files->p->pool, files->p->qsize);
    }

    if ( files->require_index==REQUIRE_IDX_ )
    {
        if ( reader->file->format.format==vcf )
        {
            if ( reader->file->format.compression!=bgzf )
            {
                files->errnum = not_bgzf;
                return 0;
            }

            reader->tbx_idx = tbx_index_load2(file_ptr->fn, idxname);
            if ( !reader->tbx_idx )
            {
                files->errnum = idx_load_failed;
                return 0;
            }

            reader->header = bcf_hdr_read(reader->file);
        }
        else if ( reader->file->format.format==bcf )
        {
            if ( reader->file->format.compression!=bgzf )
            {
                files->errnum = not_bgzf;
                return 0;
            }

            reader->header = bcf_hdr_read(reader->file);

            reader->bcf_idx = bcf_index_load2(file_ptr->fn, idxname);
            if ( !reader->bcf_idx )
            {
                files->errnum = idx_load_failed;
                return 0;
            }
        }
        else
        {
            files->errnum = file_type_error;
            return 0;
        }
    }
    else
    {
        if ( reader->file->format.format==bcf || reader->file->format.format==vcf )
        {
            reader->header = bcf_hdr_read(reader->file);
        }
        else
        {
            files->errnum = file_type_error;
            return 0;
        }
        files->streaming = 1;
    }
    if ( files->streaming && files->nreaders>1 )
    {
        static int no_index_warned = 0;
        if ( files->require_index==ALLOW_NO_IDX_ && !no_index_warned )
        {
            hts_log_warning("Using multiple unindexed files may produce errors, make sure chromosomes are in the same order!");
            no_index_warned = 1;
        }
        if ( files->require_index!=ALLOW_NO_IDX_ )
        {
            files->errnum = api_usage_error;
            hts_log_error("Must set require_index when the number of readers is greater than one");
            return 0;
        }
    }
    if ( files->streaming && files->regions )
    {
        files->errnum = api_usage_error;
        hts_log_error("Cannot tabix-jump in streaming mode");
        return 0;
    }
    if ( !reader->header )
    {
        files->errnum = header_error;
        return 0;
    }

    reader->fname = strdup(file_ptr->fn);
    if (!reader->fname)
        goto nomem;

    if ( files->apply_filters )
    {
        reader->filter_ids = init_filters(reader->header, files->apply_filters, &reader->nfilter_ids);
        if (!reader->filter_ids)
            goto nomem;
    }
    // Update list of chromosomes
    if ( !files->explicit_regs && !files->streaming )
    {
        int n = 0, i;
        const char **names;

        if ( !files->regions )
        {
            files->regions = bcf_sr_regions_alloc();
            if ( !files->regions )
                goto nomem;
        }

        names = reader->tbx_idx ? tbx_seqnames(reader->tbx_idx, &n) : bcf_hdr_seqnames(reader->header, &n);
        for (i=0; i<n; i++)
        {
            if (_regions_add(files->regions, names[i], -1, -1) < 0)
            {
                free(names);
                goto nomem;
            }
        }
        free(names);
        _regions_sort_and_merge(files->regions);
    }

    if ( files->require_index==ALLOW_NO_IDX_ && files->nreaders > 1 )
    {
        bcf_hdr_t *hdr0 = files->readers[0].header;
        bcf_hdr_t *hdr1 = reader->header;
        if ( hdr0->n[BCF_DT_CTG]!=hdr1->n[BCF_DT_CTG] )
        {
            files->errnum = noidx_error;
            hts_log_error("Different number of sequences in the header, refusing to stream multiple unindexed files");
            return 0;
        }
        int i;
        for (i=0; i<hdr0->n[BCF_DT_CTG]; i++)
        {
            if ( strcmp(bcf_hdr_id2name(hdr0,i),bcf_hdr_id2name(hdr1,i)) )
            {
                files->errnum = noidx_error;
                hts_log_error("Sequences in the header appear in different order, refusing to stream multiple unindexed files");
                return 0;
            }
        }
    }


    // Take ownership of the input file
    auxdata->closefile[files->nreaders - 1] = autoclose;

    return 1;

 nomem:
    files->errnum = no_memory;
    return 0;
}

bcf_srs_t *bcf_sr_init(void)
{
    bcf_srs_t *files = (bcf_srs_t*) calloc(1,sizeof(bcf_srs_t));
    if (!files)
        return NULL;
    files->aux = (aux_t*) calloc(1,sizeof(aux_t));
    if (!files->aux)
    {
        free(files);
        return NULL;
    }

    // sr_sort_t struct already allocated, so this cannot fail
    bcf_sr_sort_init(&BCF_SR_AUX(files)->sort);

    bcf_sr_set_opt(files,BCF_SR_REGIONS_OVERLAP,1);
    bcf_sr_set_opt(files,BCF_SR_TARGETS_OVERLAP,0);
    return files;
}

static void bcf_sr_destroy1(bcf_sr_t *reader, int closefile)
{
    if (!reader)
        return;
    free(reader->fname);
    if ( reader->tbx_idx ) tbx_destroy(reader->tbx_idx);
    if ( reader->bcf_idx ) hts_idx_destroy(reader->bcf_idx);
    bcf_hdr_destroy(reader->header);
    if (closefile) {
        if (hts_close(reader->file) < 0)
            hts_log_error("Error on closing %s", reader->fname);
    }
    if ( reader->itr ) tbx_itr_destroy(reader->itr);
    int j;
    for (j=0; j<reader->mbuffer; j++)
        bcf_destroy1(reader->buffer[j]);
    free(reader->buffer);
    free(reader->samples);
    free(reader->filter_ids);
}

void bcf_sr_destroy(bcf_srs_t *files)
{
    if (!files)
        return;

    int i;
    int *autoclose = BCF_SR_AUX(files)->closefile;

    for (i=0; i<files->nreaders; i++)
        bcf_sr_destroy1(&files->readers[i], autoclose[i]);
    free(files->has_line);
    free(files->readers);
    for (i=0; i<files->n_smpl; i++) free(files->samples[i]);
    free(files->samples);
    if (files->targets) bcf_sr_regions_destroy(files->targets);
    if (files->regions) bcf_sr_regions_destroy(files->regions);
    if (files->tmps.m) free(files->tmps.s);
    if (files->n_threads) bcf_sr_destroy_threads(files);
    bcf_sr_sort_destroy(&BCF_SR_AUX(files)->sort);
    free(autoclose);
    free(files->aux);
    free(files);
}

void bcf_sr_remove_reader(bcf_srs_t *files, int i)
{
    assert( !files->samples );  // not ready for this yet
    int *autoclose = BCF_SR_AUX(files)->closefile;

    if (i < 0 || i >= files->nreaders)
        return;

    bcf_sr_sort_remove_reader(files, &BCF_SR_AUX(files)->sort, i);
    bcf_sr_destroy1(&files->readers[i], autoclose[i]);
    if ( i+1 < files->nreaders )
    {
        memmove(&files->readers[i], &files->readers[i+1], (files->nreaders-i-1)*sizeof(bcf_sr_t));
        memmove(&files->has_line[i], &files->has_line[i+1], (files->nreaders-i-1)*sizeof(int));
        memmove(&autoclose[i], &autoclose[i+1], (files->nreaders-i-1)*sizeof(int));
    }
    files->nreaders--;
}

#if DEBUG_SYNCED_READER
void debug_buffer(FILE *fp, bcf_sr_t *reader)
{
    int j;
    for (j=0; j<=reader->nbuffer; j++)
    {
        bcf1_t *line = reader->buffer[j];
        fprintf(fp,"\t%p\t%s%s\t%s:%"PRIhts_pos"\t%s ", (void*)line,reader->fname,j==0?"*":" ",reader->header->id[BCF_DT_CTG][line->rid].key,line->pos+1,line->n_allele?line->d.allele[0]:"");
        int k;
        for (k=1; k<line->n_allele; k++) fprintf(fp," %s", line->d.allele[k]);
        fprintf(fp,"\n");
    }
}

void debug_buffers(FILE *fp, bcf_srs_t *files)
{
    int i;
    for (i=0; i<files->nreaders; i++)
    {
        fprintf(fp, "has_line: %d\t%s\n", bcf_sr_has_line(files,i),files->readers[i].fname);
        debug_buffer(fp, &files->readers[i]);
    }
    fprintf(fp,"\n");
}
#endif

static inline int has_filter(bcf_sr_t *reader, bcf1_t *line)
{
    int i, j;
    if ( !line->d.n_flt )
    {
        for (j=0; j<reader->nfilter_ids; j++)
            if ( reader->filter_ids[j]<0 ) return 1;
        return 0;
    }
    for (i=0; i<line->d.n_flt; i++)
    {
        for (j=0; j<reader->nfilter_ids; j++)
            if ( line->d.flt[i]==reader->filter_ids[j] ) return 1;
    }
    return 0;
}

// Reset a reader and make it iterate over the range seq:start-end
// Returns  0 on success
//         -1 if seq isn't in the index
//         -2 if it couldn't make a new iterator
static int _reader_seek(bcf_sr_t *reader, const char *seq, hts_pos_t start, hts_pos_t end)
{
    if ( end>=MAX_CSI_COOR )
    {
        hts_log_error("The coordinate is out of csi index limit: %"PRIhts_pos, end+1);
        return -2;
    }
    if ( reader->itr )
    {
        hts_itr_destroy(reader->itr);
        reader->itr = NULL;
    }
    reader->nbuffer = 0;
    if ( reader->tbx_idx )
    {
        int tid = tbx_name2id(reader->tbx_idx, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        reader->itr = tbx_itr_queryi(reader->tbx_idx,tid,start,end+1);
    }
    else
    {
        int tid = bcf_hdr_name2id(reader->header, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        reader->itr = bcf_itr_queryi(reader->bcf_idx,tid,start,end+1);
    }
    if (!reader->itr) {
        hts_log_error("Could not seek: %s:%"PRIhts_pos"-%"PRIhts_pos, seq, start + 1, end + 1);
        return -2;
    }
    return 0;
}

/*
 *  _readers_next_region() - jumps to next region if necessary
 *  Returns  0 on success
 *          -1 when there are no more regions left
 *          -2 on error (also sets files->errnum)
 */
static int _readers_next_region(bcf_srs_t *files)
{
    // Need to open new chromosome? Check number of lines in all readers' buffers
    int i, eos = 0;
    for (i=0; i<files->nreaders; i++)
        if ( !files->readers[i].itr && !files->readers[i].nbuffer ) eos++;

    if ( eos!=files->nreaders )
    {
        // Some of the readers still has buffered lines
        return 0;
    }

    // No lines in the buffer, need to open new region or quit.
    int prev_iseq = files->regions->iseq;
    hts_pos_t prev_end = files->regions->end;
    int res = bcf_sr_regions_next(files->regions);
    if ( res<0 )
    {
        if (res < -1) files->errnum = bcf_sr_seek_error;
        return res;
    }
    files->regions->prev_end = prev_iseq==files->regions->iseq ? prev_end : -1;

    for (i=0; i<files->nreaders; i++)
    {
        res = _reader_seek(&files->readers[i],
                           files->regions->seq_names[files->regions->iseq],
                           files->regions->start,files->regions->end);
        if (res < -1)
        {
            files->errnum = bcf_sr_seek_error;
            return res;
        }
    }

    return 0;
}

static void _set_variant_boundaries(bcf1_t *rec, hts_pos_t *beg, hts_pos_t *end)
{
    hts_pos_t off;
    if ( rec->n_allele )
    {
        off = rec->rlen;
        bcf_unpack(rec, BCF_UN_STR);
        int i;
        for (i=1; i<rec->n_allele; i++)
        {
            // Make symbolic alleles start at POS, although this is not strictly true for
            // <DEL>,<INS> where POS should be the position BEFORE the deletion/insertion.
            // However, since arbitrary symbolic alleles can be defined by the user, we
            // will simplify the interpretation of --targets-overlap and --region-overlap.
            int j = 0;
            char *ref = rec->d.allele[0];
            char *alt = rec->d.allele[i];
            while ( ref[j] && alt[j] && ref[j]==alt[j] ) j++;
            if ( off > j ) off = j;
            if ( !off ) break;
        }
    }
    else
        off = 0;

    *beg = rec->pos + off;
    *end = rec->pos + rec->rlen - 1;
}

/*
 *  _reader_fill_buffer() - buffers all records with the same coordinate
 *  returns 0 on success, -1 on failure (also sets files->errnum)
 */
static int _reader_fill_buffer(bcf_srs_t *files, bcf_sr_t *reader)
{
    // Return if the buffer is full: the coordinate of the last buffered record differs
    if ( reader->nbuffer && reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) return 0;

    // No iterator (sequence not present in this file) and not streaming
    if ( !reader->itr && !files->streaming ) return 0;

    // Fill the buffer with records starting at the same position
    int i, ret = 0;
    bcf_sr_error err = bcf_sr_ok;
    while (1)
    {
        if ( reader->mbuffer - reader->nbuffer <= 1 )
        {
            // Increase buffer size
            bcf1_t **new_buf = hts_realloc_ps(reader->buffer, sizeof(bcf1_t*),
                                              reader->mbuffer, 8);
            if (!new_buf)
            {
                err = files->errnum = no_memory;
                break;
            }
            reader->buffer = new_buf;
            reader->mbuffer += 8;
            for (i=8; i>0; i--)     // initialize
            {
                reader->buffer[reader->mbuffer-i] = bcf_init1();
                if (!reader->buffer[reader->mbuffer-i])
                    break;
                reader->buffer[reader->mbuffer-i]->max_unpack = files->max_unpack;
                reader->buffer[reader->mbuffer-i]->pos = -1;    // for rare cases when VCF starts from 1
            }
            if (i > 0)
            {
                reader->mbuffer -= i; // Adjust to number that got initialized
                err = files->errnum = no_memory;
                break;
            }
        }
        if ( files->streaming )
        {
            if ( reader->file->format.format==vcf )
            {
                ret = hts_getline(reader->file, KS_SEP_LINE, &files->tmps);
                if ( ret < -1 ) err = files->errnum = bcf_read_error;
                if ( ret < 0 ) break; // no more lines or an error
                ret = vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret<0 ) { err = files->errnum = vcf_parse_error; break; }
            }
            else if ( reader->file->format.format==bcf )
            {
                ret = bcf_read1(reader->file, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret < -1 ) err = files->errnum = bcf_read_error;
                if ( ret < 0 ) break; // no more lines or an error
            }
            else
            {
                hts_log_error("Unsupported file type for streaming");
                err = files->errnum = file_type_error;
                break;
            }
        }
        else if ( reader->tbx_idx )
        {
            ret = tbx_itr_next(reader->file, reader->tbx_idx, reader->itr, &files->tmps);
            if ( ret < -1 ) err = files->errnum = bcf_read_error;
            if ( ret < 0 ) break; // no more lines or an error
            ret = vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
            if ( ret<0 ) { err = files->errnum = vcf_parse_error; break; }
        }
        else
        {
            ret = bcf_itr_next(reader->file, reader->itr, reader->buffer[reader->nbuffer+1]);
            if ( ret < -1 ) err = files->errnum = bcf_read_error;
            if ( ret < 0 ) break; // no more lines or an error
            bcf_subset_format(reader->header,reader->buffer[reader->nbuffer+1]);
        }

        // Prevent creation of duplicates from records overlapping multiple regions
        // and recognize true variant overlaps vs record overlaps (e.g. TA>T vs A>-)
        if ( files->regions )
        {
            hts_pos_t beg, end;
            if ( BCF_SR_AUX(files)->regions_overlap==0 )
                beg = end = reader->buffer[reader->nbuffer+1]->pos;
            else if ( BCF_SR_AUX(files)->regions_overlap==1 )
            {
                beg = reader->buffer[reader->nbuffer+1]->pos;
                end = reader->buffer[reader->nbuffer+1]->pos + reader->buffer[reader->nbuffer+1]->rlen - 1;
            }
            else if ( BCF_SR_AUX(files)->regions_overlap==2 )
                _set_variant_boundaries(reader->buffer[reader->nbuffer+1], &beg,&end);
            else
            {
                hts_log_error("Invalid BCF_SR_REGIONS_OVERLAP value: %d",
                              BCF_SR_AUX(files)->regions_overlap);
                err = files->errnum = api_usage_error;
                break;
            }
            if ( beg <= files->regions->prev_end || end < files->regions->start || beg > files->regions->end ) continue;
        }

        // apply filter
        if ( !reader->nfilter_ids )
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR);
        else
        {
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR|BCF_UN_FLT);
            if ( !has_filter(reader, reader->buffer[reader->nbuffer+1]) ) continue;
        }
        reader->nbuffer++;

        if ( reader->buffer[reader->nbuffer]->rid != reader->buffer[1]->rid ) break;
        if ( reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) break;    // the buffer is full
    }
    if ( ret<0 )
    {
        // done for this region
        tbx_itr_destroy(reader->itr);
        reader->itr = NULL;
    }
    if ( files->require_index==ALLOW_NO_IDX_ && reader->buffer[reader->nbuffer]->rid < reader->buffer[1]->rid )
    {
         hts_log_error("Sequences out of order, cannot stream multiple unindexed files: %s", reader->fname);
         err = files->errnum = noidx_error;
    }
    return err == bcf_sr_ok ? 0 : -1;
}

/*
 *  _readers_shift_buffer() - removes the first line
 */
static void _reader_shift_buffer(bcf_sr_t *reader)
{
    if ( !reader->nbuffer ) return;
    int i;
    bcf1_t *tmp = reader->buffer[1];
    for (i=2; i<=reader->nbuffer; i++)
        reader->buffer[i-1] = reader->buffer[i];
    if ( reader->nbuffer > 1 )
        reader->buffer[reader->nbuffer] = tmp;
    reader->nbuffer--;
}

// Returns number of active readers (may be 0)
// On error, returns 0 but readers->errnum will have been set

static int next_line(bcf_srs_t *files)
{
    const char *chr = NULL;
    hts_pos_t min_pos = HTS_POS_MAX;

    // Loop until next suitable line is found or all readers have finished
    while ( 1 )
    {
        // Get all readers ready for the next region.
        if ( files->regions )
        {
            int r = _readers_next_region(files);
            if ( r < -1 ) return 0; // Will have set files->errnum
            if ( r <  0 ) break;
        }
        // Fill buffers and find the minimum chromosome
        int i, min_rid = INT32_MAX;
        for (i=0; i<files->nreaders; i++)
        {
            if (_reader_fill_buffer(files, &files->readers[i]) < 0)
                return 0; // Will have set files->errnum
            if ( files->require_index==ALLOW_NO_IDX_ )
            {
                if ( !files->readers[i].nbuffer ) continue;
                if ( min_rid > files->readers[i].buffer[1]->rid ) min_rid = files->readers[i].buffer[1]->rid;
            }
        }

        for (i=0; i<files->nreaders; i++)
        {
            if ( !files->readers[i].nbuffer ) continue;
            if ( files->require_index==ALLOW_NO_IDX_ && min_rid != files->readers[i].buffer[1]->rid ) continue;

            // Update the minimum coordinate
            if ( min_pos > files->readers[i].buffer[1]->pos )
            {
                min_pos = files->readers[i].buffer[1]->pos;
                chr = bcf_seqname(files->readers[i].header, files->readers[i].buffer[1]);
                if (!chr)
                {
                    files->errnum = bcf_read_error;
                    return 0;
                }
                if (bcf_sr_sort_set_active(&BCF_SR_AUX(files)->sort, i) < 0)
                {
                    files->errnum = no_memory;
                    return 0;
                }
            }
            else if ( min_pos==files->readers[i].buffer[1]->pos )
            {
                if (bcf_sr_sort_add_active(&BCF_SR_AUX(files)->sort, i) < 0)
                {
                    files->errnum = no_memory;
                    return 0;
                }
            }
        }
        if ( min_pos==HTS_POS_MAX )
        {
            if ( !files->regions ) break;
            continue;
        }

        // Skip this position if not present in targets
        if ( files->targets )
        {
            int match = 0;
            for (i=0; i<files->nreaders; i++)
            {
                if ( !files->readers[i].nbuffer || files->readers[i].buffer[1]->pos!=min_pos ) continue;
                hts_pos_t beg, end;
                if ( BCF_SR_AUX(files)->targets_overlap==0 )
                    beg = end = min_pos;
                else if ( BCF_SR_AUX(files)->targets_overlap==1 )
                {
                    beg = min_pos;
                    end = min_pos + files->readers[i].buffer[1]->rlen - 1;
                }
                else if ( BCF_SR_AUX(files)->targets_overlap==2 )
                    _set_variant_boundaries(files->readers[i].buffer[1], &beg,&end);
                else
                {
                    hts_log_error("Invalid BCF_SR_TARGETS_OVERLAP value: %d",
                                  BCF_SR_AUX(files)->targets_overlap);
                    files->errnum = api_usage_error;
                    return 0;
                }
                int overlap = bcf_sr_regions_overlap(files->targets, chr, beg, end)==0 ? 1 : 0;
                if (overlap < -1)
                {
                    files->errnum = no_memory;
                    return 0;
                }
                if ( (!files->targets_exclude && !overlap) || (files->targets_exclude && overlap) )
                    _reader_shift_buffer(&files->readers[i]);
                else
                    match = 1;
            }
            if ( !match )
            {
                min_pos = HTS_POS_MAX;
                chr = NULL;
                continue;
            }
        }
        break;  // done: chr and min_pos are set
    }
    if ( !chr ) return 0;

    return bcf_sr_sort_next(files, &BCF_SR_AUX(files)->sort, chr, min_pos);
}

int bcf_sr_next_line(bcf_srs_t *files)
{
    if ( !files->targets_als )
        return next_line(files);

    while (1)
    {
        int i, ret = next_line(files);
        if ( !ret ) return ret;

        for (i=0; i<files->nreaders; i++)
            if ( files->has_line[i] ) break;

        int r = _regions_match_alleles(files->targets, files->targets_als-1, files->readers[i].buffer[0], &files->errnum);
        if ( r ) return r > 0 ? ret : 0;

        // Check if there are more duplicate lines in the buffers. If not, return this line as if it
        // matched the targets, even if there is a type mismatch
        for (i=0; i<files->nreaders; i++)
        {
            if ( !files->has_line[i] ) continue;
            if ( files->readers[i].nbuffer==0 || files->readers[i].buffer[1]->pos!=files->readers[i].buffer[0]->pos ) continue;
            break;
        }
        if ( i==files->nreaders ) return ret;   // no more lines left, output even if target alleles are not of the same type
    }
}

static void bcf_sr_seek_start(bcf_srs_t *readers)
{
    bcf_sr_regions_t *reg = readers->regions;
    int i;
    for (i=0; i<reg->nseqs; i++)
        reg->regs[i].creg = -1;
    reg->iseq = 0;
    reg->start = -1;
    reg->end   = -1;
    reg->prev_seq = -1;
    reg->prev_start = -1;
    reg->prev_end   = -1;
}


int bcf_sr_seek(bcf_srs_t *readers, const char *seq, hts_pos_t pos)
{
    if ( !readers->regions ) return 0;
    bcf_sr_sort_reset(&BCF_SR_AUX(readers)->sort);
    if ( !seq && !pos )
    {
        // seek to start
        bcf_sr_seek_start(readers);
        return 0;
    }

    int i, nret = 0;

    // Need to position both the readers and the regions. The latter is a bit of a mess
    // because we can have in memory or external regions. The safe way is:
    //  - reset all regions as if they were not read from at all (bcf_sr_seek_start)
    //  - find the requested iseq (stored in the seq_hash)
    //  - position regions to the requested position (bcf_sr_regions_overlap)
    bcf_sr_seek_start(readers);
    if ( khash_str2int_get(readers->regions->seq_hash, seq, &i)>=0 ) readers->regions->iseq = i;
    int r = _bcf_sr_regions_overlap(readers->regions, seq, pos, pos, 0);
    if (r < -2)
    {
        readers->errnum = bcf_sr_seek_error;
        return 0;
    }

    for (i=0; i<readers->nreaders; i++)
    {
        r = _reader_seek(&readers->readers[i],seq,pos,MAX_CSI_COOR-1);
        if (r < -1)
        {
            readers->errnum = bcf_sr_seek_error;
            return 0;
        }
        nret += r >= 0 ? r : 0;
    }
    return nret;
}

int bcf_sr_set_samples(bcf_srs_t *files, const char *fname, int is_file)
{
    int i, j, nsmpl, free_smpl = 0;
    char **smpl = NULL;

    void *exclude = NULL;

    if (fname[0]=='^')
    {
        exclude = khash_str2int_init();
        if (!exclude)
        {
            files->errnum = no_memory;
            return 0;
        }
    }
    if ( exclude || strcmp("-",fname) ) // "-" stands for all samples
    {
        smpl = hts_readlist(fname, is_file, &nsmpl);
        if ( !smpl )
        {
            hts_log_error("Could not read the file: \"%s\"", fname);
            khash_str2int_destroy(exclude);
            files->errnum = bcf_sr_samples_error;
            return 0;
        }
        free_smpl = 1;
        if ( exclude )
        {
            for (i=0; i<nsmpl; i++)
            {
                if ( khash_str2int_inc(exclude, smpl[i]) < 0 )
                    goto nomem;
            }
        }
    }
    if ( !smpl )
    {
        smpl  = files->readers[0].header->samples;   // intersection of all samples
        nsmpl = bcf_hdr_nsamples(files->readers[0].header);
    }

    files->samples = NULL;
    files->n_smpl  = 0;
    for (i=0; i<nsmpl; i++)
    {
        if ( exclude && khash_str2int_has_key(exclude,smpl[i])  ) continue;

        int n_isec = 0;
        for (j=0; j<files->nreaders; j++)
        {
            if ( bcf_hdr_id2int(files->readers[j].header, BCF_DT_SAMPLE, smpl[i])<0 ) break;
            n_isec++;
        }
        if ( n_isec!=files->nreaders )
        {
            hts_log_warning("The sample \"%s\" was not found in %s, skipping",
                smpl[i], files->readers[n_isec].fname);
            continue;
        }

        char **new_samples = hts_realloc_ps(files->samples, sizeof(char *),
                                            files->n_smpl, 1);
        if (!new_samples)
            goto nomem;
        files->samples = new_samples;
        if ( (files->samples[files->n_smpl] = strdup(smpl[i])) == NULL )
            goto nomem;
        files->n_smpl++;
    }

    if ( exclude ) khash_str2int_destroy(exclude);
    if ( free_smpl )
    {
        for (i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
        free_smpl = 0;
    }

    if ( !files->n_smpl )
    {
        if ( files->nreaders>1 )
            hts_log_warning("No samples in common");
        return 0;
    }
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        reader->samples  = hts_malloc_p(sizeof(int), files->n_smpl);
        if (!reader->samples)
        {
            files->errnum = no_memory;
            return 0;
        }
        reader->n_smpl   = files->n_smpl;
        for (j=0; j<files->n_smpl; j++)
            reader->samples[j] = bcf_hdr_id2int(reader->header, BCF_DT_SAMPLE, files->samples[j]);
    }
    return 1;

 nomem:
    if ( exclude ) khash_str2int_destroy(exclude);
    if ( free_smpl )
    {
        for (i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
    }
    files->errnum = no_memory;
    return 0;
}

// Allocate a new region list structure.
static bcf_sr_regions_t *bcf_sr_regions_alloc(void)
{
    bcf_sr_regions_t *reg = (bcf_sr_regions_t *) calloc(1, sizeof(bcf_sr_regions_t));
    if ( !reg ) return NULL;

    reg->start = reg->end = -1;
    reg->prev_start = reg->prev_end = reg->prev_seq = -1;
    return reg;
}

// Add a new region into a list. On input the coordinates are 1-based, inclusive, then stored 0-based,
// inclusive. Sorting and merging step needed afterwards: qsort(..,cmp_regions) and merge_regions().
// Returns 0 on success, -1 on failure
static int _regions_add(bcf_sr_regions_t *reg, const char *chr, hts_pos_t start, hts_pos_t end)
{
    if ( start==-1 && end==-1 )
    {
        start = 0; end = MAX_CSI_COOR-1;
    }
    else
    {
        start--; end--; // store 0-based coordinates
    }

    if ( !reg->seq_hash )
    {
         reg->seq_hash = khash_str2int_init();
         if ( !reg->seq_hash )
             goto nomem;
    }

    int iseq;
    if ( khash_str2int_get(reg->seq_hash, chr, &iseq)<0 )
    {
        // the chromosome block does not exist
        char **new_seq_names = hts_realloc_ps(reg->seq_names,
                                              sizeof(char*), reg->nseqs, 1);
        if (!new_seq_names)
            goto nomem;
        reg->seq_names = new_seq_names;
        region_t *new_regs = hts_realloc_ps(reg->regs,sizeof(region_t),
                                            reg->nseqs, 1);
        if (!new_regs)
            goto nomem;
        reg->regs = new_regs;

        iseq = reg->nseqs;
        memset(&reg->regs[iseq],0,sizeof(region_t));
        reg->seq_names[iseq] = strdup(chr);
        if (!reg->seq_names[iseq])
            goto nomem;
        reg->regs[iseq].creg = -1;
        if (khash_str2int_set(reg->seq_hash,reg->seq_names[iseq],iseq) < 0)
        {
            free(reg->seq_names[iseq]);
            goto nomem;
        }
        reg->nseqs++;
    }

    region_t *creg = &reg->regs[iseq];
    if (hts_resize(region1_t,creg->nregs+1,&creg->mregs,&creg->regs,0) < 0)
        goto nomem;
    creg->regs[creg->nregs].start = start;
    creg->regs[creg->nregs].end   = end;
    creg->nregs++;

    return 0;

 nomem:
    return -1;
}

static int regions_cmp(const void *aptr, const void *bptr)
{
    region1_t *a = (region1_t*)aptr;
    region1_t *b = (region1_t*)bptr;
    if ( a->start < b->start ) return -1;
    if ( a->start > b->start ) return 1;
    if ( a->end < b->end ) return -1;
    if ( a->end > b->end ) return 1;
    return 0;
}
static void regions_merge(region_t *reg)
{
    int i = 0, j;
    while ( i<reg->nregs )
    {
        j = i + 1;
        while ( j<reg->nregs && reg->regs[i].end >= reg->regs[j].start )
        {
            if ( reg->regs[i].end < reg->regs[j].end ) reg->regs[i].end = reg->regs[j].end;
            reg->regs[j].start = 1;  reg->regs[j].end = 0;  // if beg>end, this region marked for skipping
            j++;
        }
        i = j;
    }
}
void _regions_sort_and_merge(bcf_sr_regions_t *reg)
{
    if ( !reg ) return;

    int i;
    for (i=0; i<reg->nseqs; i++)
    {
        qsort(reg->regs[i].regs, reg->regs[i].nregs, sizeof(*reg->regs[i].regs), regions_cmp);
        regions_merge(&reg->regs[i]);
    }
}

// File name or a list of genomic locations. If file name, NULL is returned.
// Recognises regions in the form chr, chr:pos, chr:beg-end, chr:beg-, {weird-chr-name}:pos.
// Cannot use hts_parse_region() as that requires the header and if header is not present,
// wouldn't learn the chromosome name.
// Returns populated bcf_sr_regions_t struct on success; NULL on failure
static bcf_sr_regions_t *_regions_init_string(const char *str)
{
    bcf_sr_regions_t *reg = bcf_sr_regions_alloc();
    if ( !reg )
    {
        hts_log_error("Out of memory");
        return NULL;
    }

    kstring_t tmp = {0,0,0};
    const char *sp = str, *ep = str;
    hts_pos_t from, to;
    while ( 1 )
    {
        tmp.l = 0;
        if ( *ep=='{' )
        {
            while ( *ep && *ep!='}' ) ep++;
            if ( !*ep )
            {
                hts_log_error("Could not parse the region, mismatching braces in: \"%s\"", str);
                goto exit_nicely;
            }
            ep++;
            if (kputsn(sp+1,ep-sp-2,&tmp) < 0)
                goto exit_nicely_nomem;
        }
        else
        {
            while ( *ep && *ep!=',' && *ep!=':' ) ep++;
            if (kputsn(sp,ep-sp,&tmp) < 0)
                goto exit_nicely_nomem;
        }
        if ( *ep==':' )
        {
            sp = ep+1;
            from = hts_parse_decimal(sp,(char**)&ep,0);
            if ( sp==ep )
            {
                hts_log_error("Could not parse the region(s): %s", str);
                goto exit_nicely;
            }
            if ( !*ep || *ep==',' )
            {
                if (_regions_add(reg, tmp.s, from, from) < 0)
                    goto exit_nicely_nomem;
                sp = ep;
                continue;
            }
            if ( *ep!='-' )
            {
                hts_log_error("Could not parse the region(s): %s", str);
                goto exit_nicely;
            }
            ep++;
            sp = ep;
            to = hts_parse_decimal(sp,(char**)&ep,0);
            if ( *ep && *ep!=',' )
            {
                hts_log_error("Could not parse the region(s): %s", str);
                goto exit_nicely;
            }
            if ( sp==ep ) to = MAX_CSI_COOR-1;
            if (_regions_add(reg, tmp.s, from, to) < 0)
                goto exit_nicely_nomem;
            if ( !*ep ) break;
            sp = ep;
        }
        else if ( !*ep || *ep==',' )
        {
            if ( tmp.l )
            {
                if (_regions_add(reg, tmp.s, -1, -1) < 0)
                    goto exit_nicely_nomem;
            }
            if ( !*ep ) break;
            sp = ++ep;
        }
        else
        {
            hts_log_error("Could not parse the region(s): %s", str);
            goto exit_nicely;
        }
    }
    free(tmp.s);
    return reg;

exit_nicely_nomem:
    hts_log_error("Out of memory");

exit_nicely:
    bcf_sr_regions_destroy(reg);
    free(tmp.s);
    return NULL;
}

// ichr,ifrom,ito are 0-based;
// returns -1 on error, 0 if the line is a comment line, 1 on success
static int _regions_parse_line(char *line, int ichr, int ifrom, int ito, char **chr, char **chr_end, hts_pos_t *from, hts_pos_t *to)
{
    if (ifrom < 0 || ito < 0) return -1;
    *chr_end = NULL;

    if ( line[0]=='#' ) return 0;

    int k,l;    // index of the start and end column of the tab-delimited file
    if ( ifrom <= ito )
        k = ifrom, l = ito;
    else
        l = ifrom, k = ito;

    int i;
    char *se = line, *ss = NULL; // start and end
    char *tmp;
    for (i=0; i<=k && *se; i++)
    {
        ss = i==0 ? se++ : ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=k ) return -1;
    if ( k==l )
    {
        *from = *to = hts_parse_decimal(ss, &tmp, 0);
        if ( tmp==ss || (*tmp && *tmp!='\t') ) return -1;
    }
    else
    {
        if ( k==ifrom )
            *from = hts_parse_decimal(ss, &tmp, 0);
        else
            *to = hts_parse_decimal(ss, &tmp, 0);
        if ( ss==tmp || (*tmp && *tmp!='\t') ) return -1;

        for (i=k; i<l && *se; i++)
        {
            ss = ++se;
            while (*se && *se!='\t') se++;
        }
        if ( i<l ) return -1;
        if ( k==ifrom )
            *to = hts_parse_decimal(ss, &tmp, 0);
        else
            *from = hts_parse_decimal(ss, &tmp, 0);
        if ( ss==tmp || (*tmp && *tmp!='\t') ) return -1;
    }

    ss = se = line;
    for (i=0; i<=ichr && *se; i++)
    {
        if ( i>0 ) ss = ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=ichr ) return -1;
    *chr_end = se;
    *chr = ss;
    return 1;
}

bcf_sr_regions_t *bcf_sr_regions_init(const char *regions, int is_file, int ichr, int ifrom, int ito)
{
    bcf_sr_regions_t *reg;

    if (!regions)
        return NULL;

    if ( !is_file )
    {
        reg = _regions_init_string(regions);
        if (!reg)
            return NULL;
        _regions_sort_and_merge(reg);
        return reg;
    }

    reg = bcf_sr_regions_alloc();
    if ( !reg ) return NULL;

    reg->file = hts_open(regions, "rb");
    if ( !reg->file )
    {
        hts_log_error("Could not open file: %s", regions);
        free(reg);
        return NULL;
    }

    reg->tbx = tbx_index_load3(regions, NULL, HTS_IDX_SAVE_REMOTE|HTS_IDX_SILENT_FAIL);
    if ( !reg->tbx )
    {
        size_t iline = 0;
        int len = strlen(regions);
        int is_bed  = strcasecmp(".bed",regions+len-4) ? 0 : 1;
        int r;
        if ( !is_bed && !strcasecmp(".bed.gz",regions+len-7) ) is_bed = 1;

        if ( reg->file->format.format==vcf ) ito = 1;

        // read the whole file, tabix index is not present
        while ( ( r = hts_getline(reg->file, KS_SEP_LINE, &reg->line) ) > 0 )
        {
            iline++;
            char *chr, *chr_end;
            hts_pos_t from, to;
            int ret;
            ret = _regions_parse_line(reg->line.s, ichr,ifrom,abs(ito), &chr,&chr_end,&from,&to);
            if ( ret < 0 )
            {
                if ( ito<0 )
                    ret = _regions_parse_line(reg->line.s, ichr,ifrom,ifrom, &chr,&chr_end,&from,&to);
                if ( ret<0 )
                {
                    hts_log_error("Could not parse %zu-th line of file %s, using the columns %d,%d[,%d]",
                        iline, regions,ichr+1,ifrom+1,ito+1);
                    bcf_sr_regions_destroy(reg);
                    return NULL;
                }
                ito = ifrom;
            }
            else if ( ito<0 )
                ito = abs(ito);
            if ( !ret ) continue;
            if ( is_bed ) from++;
            *chr_end = 0;
            if (_regions_add(reg, chr, from, to) < 0)
                goto nomem;
            *chr_end = '\t';
        }
        if ( r < -1 )
        {
            hts_log_error("Error reading file: \"%s\" : %s",
                          regions, strerror(errno));
            bcf_sr_regions_destroy(reg);
            return NULL;
        }
        hts_close(reg->file); reg->file = NULL;
        if ( !reg->nseqs ) { free(reg); return NULL; }
        _regions_sort_and_merge(reg);
        return reg;
    }

    reg->seq_names = (char**) tbx_seqnames(reg->tbx, &reg->nseqs);
    if ( !reg->seq_hash )
    {
        reg->seq_hash = khash_str2int_init();
        if (!reg->seq_hash)
            goto nomem;
    }
    int i;
    for (i=0; i<reg->nseqs; i++)
    {
        if (khash_str2int_set(reg->seq_hash,reg->seq_names[i],i) < 0)
            goto nomem;
    }
    reg->fname  = strdup(regions);
    if (!reg->fname)
        goto nomem;
    reg->is_bin = 1;
    return reg;

 nomem:
    hts_log_error("Out of memory");
    bcf_sr_regions_destroy(reg);
    return NULL;
}

void bcf_sr_regions_destroy(bcf_sr_regions_t *reg)
{
    int i;
    free(reg->fname);
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    if ( reg->tbx ) tbx_destroy(reg->tbx);
    if ( reg->file ) hts_close(reg->file);
    if ( reg->als ) free(reg->als);
    if ( reg->als_str.s ) free(reg->als_str.s);
    free(reg->line.s);
    if ( reg->regs )
    {
         // free only in-memory names, tbx names are const
        for (i=0; i<reg->nseqs; i++)
        {
            free(reg->seq_names[i]);
            free(reg->regs[i].regs);
        }
    }
    free(reg->regs);
    free(reg->seq_names);
    khash_str2int_destroy(reg->seq_hash);
    free(reg);
}

int bcf_sr_regions_seek(bcf_sr_regions_t *reg, const char *seq)
{
    reg->iseq = reg->start = reg->end = -1;
    if ( khash_str2int_get(reg->seq_hash, seq, &reg->iseq) < 0 ) return -1;  // sequence seq not in regions

    // using in-memory regions
    if ( reg->regs )
    {
        reg->regs[reg->iseq].creg = -1;
        return 0;
    }

    // reading regions from tabix
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    reg->itr = tbx_itr_querys(reg->tbx, seq);
    if ( reg->itr ) return 0;

    return -1;
}

// Returns 0 on success, -1 when done
static int advance_creg(region_t *reg)
{
    int i = reg->creg + 1;
    while ( i<reg->nregs && reg->regs[i].start > reg->regs[i].end ) i++;    // regions with start>end are marked to skip by merge_regions()
    reg->creg = i;
    if ( i>=reg->nregs ) return -1;
    return 0;
}

int bcf_sr_regions_next(bcf_sr_regions_t *reg)
{
    if ( reg->iseq<0 ) return -1;
    reg->start = reg->end = -1;
    reg->nals = 0;

    // using in-memory regions
    if ( reg->regs )
    {
        while ( reg->iseq < reg->nseqs )
        {
            if ( advance_creg(&reg->regs[reg->iseq])==0 ) break;    // a valid record was found
            reg->iseq++;
        }
        if ( reg->iseq >= reg->nseqs ) { reg->iseq = -1; return -1; } // no more regions left
        region1_t *creg = &reg->regs[reg->iseq].regs[reg->regs[reg->iseq].creg];
        reg->start = creg->start;
        reg->end   = creg->end;
        return 0;
    }

    // reading from tabix
    char *chr, *chr_end;
    int ichr = 0, ifrom = 1, ito = 2, is_bed = 0;
    hts_pos_t from, to;
    if ( reg->tbx )
    {
        ichr   = reg->tbx->conf.sc-1;
        ifrom  = reg->tbx->conf.bc-1;
        ito    = reg->tbx->conf.ec-1;
        if ( ito<0 ) ito = ifrom;
        is_bed = reg->tbx->conf.preset==TBX_UCSC ? 1 : 0;
    }

    int ret = 0;
    while ( !ret )
    {
        if ( reg->itr )
        {
            // tabix index present, reading a chromosome block
            ret = tbx_itr_next(reg->file, reg->tbx, reg->itr, &reg->line);
            if ( ret<0 ) { reg->iseq = -1; return ret; }
        }
        else
        {
            if ( reg->is_bin )
            {
                // Waited for seek which never came. Reopen in text mode and stream
                // through the regions, otherwise hts_getline would fail
                hts_close(reg->file);
                reg->file = hts_open(reg->fname, "r");
                if ( !reg->file )
                {
                    hts_log_error("Could not open file: %s", reg->fname);
                    reg->file = NULL;
                    bcf_sr_regions_destroy(reg);
                    return -2;
                }
                reg->is_bin = 0;
            }

            // tabix index absent, reading the whole file
            ret = reg->file ? hts_getline(reg->file, KS_SEP_LINE, &reg->line) : -1;
            if ( ret<-1)
            {
                hts_log_error("Error reading \"%s\": %s\n",
                              reg->fname, strerror(errno));
                return -2;
            }
            if ( ret<0 ) { reg->iseq = -1; return ret; }
        }
        ret = _regions_parse_line(reg->line.s, ichr,ifrom,ito, &chr,&chr_end,&from,&to);
        if ( ret<0 )
        {
            hts_log_error("Could not parse the file %s, using the columns %d,%d,%d",
                reg->fname,ichr+1,ifrom+1,ito+1);
            return -2;
        }
    }
    if ( is_bed ) from++;

    *chr_end = 0;
    if ( khash_str2int_get(reg->seq_hash, chr, &reg->iseq)<0 )
    {
        hts_log_error("Broken tabix index? The sequence \"%s\" not in dictionary [%s]",
            chr, reg->line.s);
        return -2;
    }
    *chr_end = '\t';

    reg->start = from - 1;
    reg->end   = to - 1;
    return 0;
}

static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec, bcf_sr_error *errnum)
{
    if ( reg->regs )
    {
        // payload is not supported for in-memory regions, switch to regidx instead in future
        hts_log_error("Compressed and indexed targets file is required");
        *errnum = api_usage_error;
        return -1;
    }

    int i = 0, max_len = 0;
    if ( !reg->nals )
    {
        char *ss = reg->line.s;
        while ( i<als_idx && *ss )
        {
            if ( *ss=='\t' ) i++;
            ss++;
        }
        char *se = ss;
        reg->nals = 1;
        while ( *se && *se!='\t' )
        {
            if ( *se==',' ) reg->nals++;
            se++;
        }
        if (ks_resize(&reg->als_str, se-ss+1+reg->nals) < 0)
            goto nomem;
        reg->als_str.l = 0;
        if (hts_resize(char*,reg->nals,&reg->mals,&reg->als, 0) < 0)
            goto nomem;
        reg->nals = 0;

        se = ss;
        while ( *(++se) )
        {
            if ( *se=='\t' ) break;
            if ( *se!=',' ) continue;
            reg->als[reg->nals] = &reg->als_str.s[reg->als_str.l];
            if (kputsn(ss,se-ss,&reg->als_str) < 0)
                goto nomem;
            if ( &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals] > max_len ) max_len = &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals];
            reg->als_str.l++;
            reg->nals++;
            ss = ++se;
        }
        reg->als[reg->nals] = &reg->als_str.s[reg->als_str.l];
        if (kputsn(ss,se-ss,&reg->als_str) < 0)
            goto nomem;
        if ( &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals] > max_len ) max_len = &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals];
        reg->nals++;
        reg->als_type = max_len > 1 ? VCF_INDEL : VCF_SNP;  // this is a simplified check, see vcf.c:bcf_set_variant_types
    }
    int type = bcf_get_variant_types(rec);
    if ( reg->als_type & VCF_INDEL )
        return type & VCF_INDEL ? 1 : 0;
    return !(type & VCF_INDEL) ? 1 : 0;

 nomem:
    *errnum = no_memory;
    return -1;
}

int bcf_sr_regions_overlap(bcf_sr_regions_t *reg, const char *seq, hts_pos_t start, hts_pos_t end)
{
    return _bcf_sr_regions_overlap(reg,seq,start,end,1);
}

static int _bcf_sr_regions_overlap(bcf_sr_regions_t *reg, const char *seq, hts_pos_t start, hts_pos_t end, int missed_reg_handler)
{
    int iseq;
    if ( khash_str2int_get(reg->seq_hash, seq, &iseq)<0 ) return -1;    // no such sequence
    if ( missed_reg_handler && !reg->missed_reg_handler ) missed_reg_handler = 0;

    if ( reg->prev_seq==-1 || iseq!=reg->prev_seq || reg->prev_start > start ) // new chromosome or after a seek
    {
        // flush regions left on previous chromosome
        if ( missed_reg_handler && reg->prev_seq!=-1 && reg->iseq!=-1 )
        {
            if ( bcf_sr_regions_flush(reg) < 0 )
                return -3;
        }
        if ( bcf_sr_regions_seek(reg, seq) < 0 )
            return -3;
        reg->start = reg->end = -1;
    }
    if ( reg->prev_seq==iseq && reg->iseq!=iseq ) return -2;    // no more regions on this chromosome
    reg->prev_seq = reg->iseq;
    reg->prev_start = start;

    while ( iseq==reg->iseq && reg->end < start )
    {
        int r = bcf_sr_regions_next(reg);
        if ( r < -1 ) return -3; // Something went wrong
        if ( r <  0 ) return -2;  // no more regions left
        if ( reg->iseq != iseq ) return -1; // does not overlap any regions
        if ( missed_reg_handler && reg->end < start ) reg->missed_reg_handler(reg, reg->missed_reg_data);
    }
    if ( reg->start <= end ) return 0;    // region overlap
    return -1;  // no overlap
}

int bcf_sr_regions_flush(bcf_sr_regions_t *reg)
{
    int r;
    if ( !reg->missed_reg_handler || reg->prev_seq==-1 ) return 0;
    while ( (r = bcf_sr_regions_next(reg)) == 0 )
        reg->missed_reg_handler(reg, reg->missed_reg_data);
    return r > -2 ? 0 : -1;
}

