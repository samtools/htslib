#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>
#include "htslib/synced_bcf_reader.h"
#include "htslib/kseq.h"

typedef struct _region_t 
{
    char *chr;
    int start, end;
}
region_t;

static void _regions_add(bcf_sr_regions_t *reg, const char *chr, int start, int end);
static bcf_sr_regions_t *_regions_init_string(const char *str);
static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec);

static int *init_filters(bcf_hdr_t *hdr, const char *filters, int *nfilters)
{
    kstring_t str = {0,0,0};
    const char *tmp = filters, *prev = filters;
    int nout = 0, *out = NULL;
    while ( 1 )
    {
        if ( *tmp==',' || !*tmp )
        {
            out = (int*) realloc(out, sizeof(int));
            if ( tmp-prev==1 && *prev=='.' )
                out[nout] = -1;
            else
            {
                str.l = 0;
                kputsn(prev, tmp-prev, &str);
                out[nout] = bcf_hdr_id2int(hdr, BCF_DT_ID, str.s);
            }
            nout++;
            if ( !*tmp ) break;
            prev = tmp+1;
        }
        tmp++;
    }
    if ( str.m ) free(str.s);
    *nfilters = nout;
    return out;
}

int bcf_sr_set_regions(bcf_srs_t *readers, const char *regions)
{
    assert( !readers->regions );
    if ( readers->nreaders ) 
    {
        fprintf(stderr,"[%s:%d %s] Error: bcf_sr_set_regions() must be called before bcf_sr_add_reader()\n", __FILE__,__LINE__,__FUNCTION__);
        return -1;
    }
    readers->regions = bcf_sr_regions_init(regions);
    if ( !readers->regions ) return -1;
    readers->explicit_regs = 1;
    readers->require_index = 1;
    return 0;
}
int bcf_sr_set_targets(bcf_srs_t *readers, const char *targets, int alleles)
{
    assert( !readers->targets );
    readers->targets = bcf_sr_regions_init(targets);
    if ( !readers->targets ) return -1;
    readers->targets_als = alleles;
    return 0;
}

int bcf_sr_add_reader(bcf_srs_t *files, const char *fname)
{
    files->has_line = (int*) realloc(files->has_line, sizeof(int)*(files->nreaders+1));
    files->readers  = (bcf_sr_t*) realloc(files->readers, sizeof(bcf_sr_t)*(files->nreaders+1));
    bcf_sr_t *reader = &files->readers[files->nreaders++];
    memset(reader,0,sizeof(bcf_sr_t));

    reader->file = hts_open(fname, "r");
    if ( !reader->file ) return 0;

    reader->type = reader->file->is_bin? FT_BCF : FT_VCF;
    if (reader->file->is_compressed) reader->type |= FT_GZ;

    if ( files->require_index )
    {
        if ( reader->type==FT_VCF_GZ ) 
        {
            reader->tbx_idx = tbx_index_load(fname);
            if ( !reader->tbx_idx )
            {
                fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
                return 0;
            }

            reader->header = bcf_hdr_read(reader->file);
        }
        else if ( reader->type==FT_BCF_GZ ) 
        {
            reader->header = bcf_hdr_read(reader->file);

            reader->bcf_idx = bcf_index_load(fname);
            if ( !reader->bcf_idx ) 
            {
                fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
                return 0;   // not indexed..?
            }
        }
        else
        {
            fprintf(stderr,"Index required, expected .vcf.gz or .bcf file: %s\n", fname);
            return 0;
        }
    }
    else 
    {
        if ( reader->type & FT_BCF )
        {
            reader->header = bcf_hdr_read(reader->file);
        }
        else if ( reader->type & FT_VCF )
        {
            reader->header = bcf_hdr_read(reader->file);
        }
        else
        {
            fprintf(stderr,"File type not recognised: %s\n", fname);
            return 0;
        }
        files->streaming = 1;
    }
    if ( files->streaming && files->nreaders>1 )
    {
        fprintf(stderr,"[%s:%d %s] Error: %d readers, yet require_index not set\n", __FILE__,__LINE__,__FUNCTION__,files->nreaders);
        return 0;
    }
    if ( files->streaming && files->regions )
    {
        fprintf(stderr,"[%s:%d %s] Error: cannot tabix-jump in streaming mode\n", __FILE__,__LINE__,__FUNCTION__);
        return 0;
    }
    if ( !reader->header ) return 0;

    reader->fname = fname;
    if ( files->apply_filters )
        reader->filter_ids = init_filters(reader->header, files->apply_filters, &reader->nfilter_ids);

    // Update list of chromosomes
    if ( !files->explicit_regs && !files->streaming )
    {
        int n,i;
        const char **names = reader->tbx_idx ? tbx_seqnames(reader->tbx_idx, &n) : bcf_hdr_seqnames(reader->header, &n);
        for (i=0; i<n; i++)
        {
            if ( !files->regions )
                files->regions = _regions_init_string(names[i]);
            else
                _regions_add(files->regions, names[i], -1, -1);
        }
        free(names);
    }
    return 1;
}

bcf_srs_t *bcf_sr_init(void)
{
    bcf_srs_t *files = (bcf_srs_t*) calloc(1,sizeof(bcf_srs_t));
    files->crid = -1;
    return files;
}

void bcf_sr_destroy(bcf_srs_t *files)
{
    if ( !files->nreaders ) 
    {
        free(files);
        return;
    }
    int i;
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        if ( reader->tbx_idx ) tbx_destroy(reader->tbx_idx);
        if ( reader->bcf_idx ) hts_idx_destroy(reader->bcf_idx);
        bcf_hdr_destroy(reader->header);
        hts_close(reader->file);
        if ( reader->itr ) tbx_itr_destroy(reader->itr);
        int j;
        for (j=0; j<reader->mbuffer; j++)
            bcf_destroy1(reader->buffer[j]);
        free(reader->buffer);
        free(reader->samples);
        free(reader->filter_ids);
    }
    free(files->has_line);
    free(files->readers);
    for (i=0; i<files->n_smpl; i++) free(files->samples[i]);
    free(files->samples);
    if (files->targets) bcf_sr_regions_destroy(files->targets);
    if (files->regions) bcf_sr_regions_destroy(files->regions);
    if ( files->tmps.m ) free(files->tmps.s);
    free(files);
}

/*
   Removes duplicate records from the buffer. The meaning of "duplicate" is
   controlled by the $collapse variable, which can cause that from multiple
   <indel|snp|any> lines only the first is considered and the rest is ignored.
   The removal is done by setting the redundant lines' positions to -1 and
   moving these lines at the end of the buffer.
 */
static void collapse_buffer(bcf_srs_t *files, bcf_sr_t *reader)
{
    int irec,jrec, has_snp=0, has_indel=0, has_any=0;
    for (irec=1; irec<=reader->nbuffer; irec++)
    {
        bcf1_t *line = reader->buffer[irec];
        if ( line->pos != reader->buffer[1]->pos ) break;
        if ( files->collapse&COLLAPSE_ANY )
        {
            if ( !has_any ) has_any = 1;
            else line->pos = -1;
        }
        int line_type = bcf_get_variant_types(line);
        if ( files->collapse&COLLAPSE_SNPS && line_type&(VCF_SNP|VCF_MNP) )
        {
            if ( !has_snp ) has_snp = 1;
            else line->pos = -1;
        }
        if ( files->collapse&COLLAPSE_INDELS && line_type&VCF_INDEL )
        {
            if ( !has_indel ) has_indel = 1;
            else line->pos = -1;
        }
    }
    bcf1_t *tmp;
    irec = jrec = 1;
    while ( irec<=reader->nbuffer && jrec<=reader->nbuffer )
    {
        if ( reader->buffer[irec]->pos != -1 ) { irec++; continue; }
        if ( jrec<=irec ) jrec = irec+1;
        while ( jrec<=reader->nbuffer && reader->buffer[jrec]->pos==-1 ) jrec++;
        if ( jrec<=reader->nbuffer )
        {
            tmp = reader->buffer[irec]; reader->buffer[irec] = reader->buffer[jrec]; reader->buffer[jrec] = tmp;
        }
    }
    reader->nbuffer = irec - 1;
}

void debug_buffer(FILE *fp, bcf_sr_t *reader)
{
    int j;
    for (j=0; j<=reader->nbuffer; j++)
    {
        bcf1_t *line = reader->buffer[j];
        fprintf(fp,"%s%s\t%s:%d\t%s ", reader->fname,j==0?"*":"",reader->header->id[BCF_DT_CTG][line->rid].key,line->pos+1,line->n_allele?line->d.allele[0]:"");
        int k;
        for (k=1; k<line->n_allele; k++) fprintf(fp," %s", line->d.allele[k]);
        fprintf(fp,"\n");
    }
}

void debug_buffers(FILE *fp, bcf_srs_t *files)
{
    int i;
    for (i=0; i<files->nreaders; i++)
        debug_buffer(fp, &files->readers[i]);
    fprintf(fp,"\n");
}

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

/*
 *  _readers_next_region() - jumps to next region if necessary
 *  Returns 0 on success or -1 when there are no more regions left
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

    // No lines in the buffer, need to open new region or quit
    if ( bcf_sr_regions_next(files->regions)<0 ) return -1;

    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        if ( reader->tbx_idx )
        {
            int tid = tbx_name2id(reader->tbx_idx, files->regions->seq);
            if ( tid==-1 ) continue;    // the sequence not present in this file
            reader->itr = tbx_itr_queryi(reader->tbx_idx,tid,files->regions->start,files->regions->end+1);
        }
        else
        {
            int tid = bcf_hdr_name2id(reader->header,files->regions->seq);
            if ( tid==-1 ) continue;    // the sequence not present in this file
            reader->itr = bcf_itr_queryi(reader->bcf_idx,tid,files->regions->start,files->regions->end+1);
        }
        assert(reader->itr);
    }
    return 0;
}

/*
 *  _reader_fill_buffer() - buffers all records with the same coordinate
 */
static void _reader_fill_buffer(bcf_srs_t *files, bcf_sr_t *reader)
{
    // Return if the buffer is full: the coordinate of the last buffered record differs
    if ( reader->nbuffer && reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) return;

    // No iterator (sequence not present in this file) and not streaming
    if ( !reader->itr && !files->streaming ) return;

    // Fill the buffer with records starting at the same position
    int i, ret = 0;
    while (1)
    {
        if ( reader->nbuffer+1 >= reader->mbuffer ) 
        {
            // Increase buffer size
            reader->mbuffer += 8;
            reader->buffer = (bcf1_t**) realloc(reader->buffer, sizeof(bcf1_t*)*reader->mbuffer);
            for (i=8; i>0; i--)     // initialize
            {
                reader->buffer[reader->mbuffer-i] = bcf_init1();
                reader->buffer[reader->mbuffer-i]->max_unpack = files->max_unpack;
            }
        }
        if ( files->streaming )
        {
            if ( reader->type & FT_VCF )
            {
                if ( (ret=hts_getline(reader->file, KS_SEP_LINE, &files->tmps)) < 0 ) break;   // no more lines
                int ret = vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret<0 ) break;
            }
            else if ( reader->type & FT_BCF )
            {
                if ( (ret=bcf_read1(reader->file, reader->header, reader->buffer[reader->nbuffer+1])) < 0 ) break; // no more lines
            }
            else
            {
                fprintf(stderr,"[%s:%d %s] fixme: not ready for this\n", __FILE__,__LINE__,__FUNCTION__);
                exit(1);
            }
        }
        else if ( reader->tbx_idx )
        {
            if ( (ret=tbx_itr_next(reader->file, reader->tbx_idx, reader->itr, &files->tmps)) < 0 ) break;  // no more lines
            vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
        }
        else
            if ( (ret=bcf_itr_next(reader->file, reader->itr, reader->buffer[reader->nbuffer+1])) < 0 ) break; // no more lines

        // apply filter
        if ( !reader->nfilter_ids )
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR);
        else
        {
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR|BCF_UN_FLT);
            if ( !has_filter(reader, reader->buffer[reader->nbuffer+1]) ) continue;
        }
        reader->nbuffer++;

        if ( reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) break;    // the buffer is full
    }
    if ( ret<0 ) 
    { 
        // done for this region
        tbx_itr_destroy(reader->itr);
        reader->itr = NULL; 
    }
    if ( files->collapse && reader->nbuffer>=2 && reader->buffer[1]->pos==reader->buffer[2]->pos )
        collapse_buffer(files, reader);
}

/*
 *  _readers_shift_buffer() - removes the first line and all subsequent lines with the same position
 */
static void _reader_shift_buffer(bcf_sr_t *reader)
{
    int i;
    for (i=2; i<=reader->nbuffer; i++)
        if ( reader->buffer[i]->pos!=reader->buffer[1]->pos ) break;
    if ( i<=reader->nbuffer )
    {
        // A record with a different position follows, swap it. Because of the reader's logic,
        // only one such line can be present.
        bcf1_t *tmp = reader->buffer[1]; reader->buffer[1] = reader->buffer[i]; reader->buffer[i] = tmp;
        reader->nbuffer = 1;
    }
    else 
        reader->nbuffer = 0;    // no other line
}

/*
 *  _reader_match_alleles() - from multiple buffered lines selects the one which
 *  corresponds best to the template line. The logic is controlled by COLLAPSE_*
 *  Returns 0 on success or -1 when no good matching line is found.
 */
static int _reader_match_alleles(bcf_srs_t *files, bcf_sr_t *reader, bcf1_t *tmpl)
{
    int i, irec = -1;

    // if no template given, use the first available record
    if ( !tmpl )
        irec = 1;
    else
    {
        int tmpl_type = bcf_get_variant_types(tmpl);
        for (i=1; i<=reader->nbuffer; i++)
        {
            bcf1_t *line = reader->buffer[i];
            if ( line->pos != reader->buffer[1]->pos ) break;  // done with this reader

            // Easiest case: matching by position only
            if ( files->collapse&COLLAPSE_ANY ) { irec=i; break; }

            int line_type = bcf_get_variant_types(line);

            // No matter what the alleles are, as long as they are both SNPs
            if ( files->collapse&COLLAPSE_SNPS && tmpl_type&VCF_SNP && line_type&VCF_SNP ) { irec=i; break; }
            // ... or indels
            if ( files->collapse&COLLAPSE_INDELS && tmpl_type&VCF_INDEL && line_type&VCF_INDEL ) { irec=i; break; }

            // More thorough checking: REFs must match
            if ( tmpl->rlen != line->rlen ) continue;  // different length
            if ( strcmp(tmpl->d.allele[0], line->d.allele[0]) ) continue; // the strings do not match

            int ial,jal;
            if ( files->collapse==COLLAPSE_NONE )
            {
                // Exact match, all alleles must be identical
                if ( tmpl->n_allele!=line->n_allele ) continue;   // different number of alleles, skip

                int nmatch = 1; // REF has been already checked
                for (ial=1; ial<tmpl->n_allele; ial++)
                {
                    for (jal=1; jal<line->n_allele; jal++)
                        if ( !strcmp(tmpl->d.allele[ial], line->d.allele[jal]) ) { nmatch++; break; }
                }
                if ( nmatch==tmpl->n_allele ) { irec=i; break; }    // found: exact match
                continue;
            }
            
            // COLLAPSE_SOME: at least some ALTs must match
            for (ial=1; ial<tmpl->n_allele; ial++)
            {
                for (jal=1; jal<line->n_allele; jal++)
                    if ( !strcmp(tmpl->d.allele[ial], line->d.allele[jal]) ) { irec=i; break; }
                if ( irec>=1 ) break;
            }
            if ( irec>=1 ) break;
        }
        if ( irec==-1 ) return -1;  // no matching line was found
    }

    // Set the selected line (irec) as active: set it to buffer[0], move the remaining lines forward
    // and put the old bcf1_t record at the end.
    bcf1_t *tmp = reader->buffer[0];
    reader->buffer[0] = reader->buffer[irec];
    for (i=irec+1; i<=reader->nbuffer; i++) reader->buffer[i-1] = reader->buffer[i];
    reader->buffer[ reader->nbuffer ] = tmp;
    reader->nbuffer--;

    return 0;
}

int _reader_next_line(bcf_srs_t *files)
{
    int i, min_pos = INT_MAX;

    // Loop until next suitable line is found or all readers have finished
    while ( 1 )
    {
        // Get all readers ready for the next region.
        if ( files->regions && _readers_next_region(files)<0 ) break;

        // Fill buffers
        for (i=0; i<files->nreaders; i++)
        {
            _reader_fill_buffer(files, &files->readers[i]);

            // Update the minimum coordinate
            if ( !files->readers[i].nbuffer ) continue;
            if ( min_pos > files->readers[i].buffer[1]->pos ) min_pos = files->readers[i].buffer[1]->pos; 
        }
        if ( min_pos==INT_MAX ) 
        {
            if ( !files->regions ) break;
            continue;
        }

        // Skip this position if not present in targets
        if ( files->targets )
        {
            if  ( files->regions )
            {
                if ( !files->cseq || files->cseq!=files->regions->seq )
                {
                    bcf_sr_regions_seek(files->targets, files->regions->seq);
                    files->cseq = files->regions->seq;  // set the current sequence
                }
            }
            else
            {
                // If here, we must be streaming a single VCF, otherwise either explicit
                // or implicit regions would be set. We can safely use rid as a unique sequence
                // identifier
                int rid = files->readers[0].buffer[1]->rid;
                if ( files->crid<0 || files->crid!=rid )
                {
                    bcf_sr_regions_seek(files->targets, files->readers[0].header->id[BCF_DT_CTG][rid].key);
                    files->crid = rid;
                }
            }
            if ( bcf_sr_regions_query(files->targets, min_pos, min_pos)<0 ) 
            {
                // Remove all lines with this position from the buffer
                for (i=0; i<files->nreaders; i++)
                    if ( files->readers[i].nbuffer && files->readers[i].buffer[1]->pos==min_pos ) 
                        _reader_shift_buffer(&files->readers[i]);
                min_pos = INT_MAX;
                continue;
            }
        }
        
        break;  // done: min_pos is set 
    }

    // There can be records with duplicate positions. Set the active line intelligently so that
    // the alleles match.
    int nret = 0;   // number of readers sharing the position
    bcf1_t *first = NULL;   // record which will be used for allele matching
    for (i=0; i<files->nreaders; i++)
    {
        files->has_line[i] = 0;
        
        // Skip readers with no records at this position
        if ( !files->readers[i].nbuffer || files->readers[i].buffer[1]->pos!=min_pos ) continue;

        // Until now buffer[0] of all reader was empty and the lines started at buffer[1].
        // No lines which are ready to output will be moved to buffer[0].
        if ( _reader_match_alleles(files, &files->readers[i], first) < 0 ) continue;
        if ( !first ) first = files->readers[i].buffer[0];

        nret++;
        files->has_line[i] = 1;
    }
    return nret;
}

int bcf_sr_next_line(bcf_srs_t *files)
{
    if ( !files->targets_als ) 
        return _reader_next_line(files);

    while (1)
    {
        int i, ret = _reader_next_line(files);
        if ( !ret ) return ret;
        
        for (i=0; i<files->nreaders; i++)
            if ( files->has_line[i] ) break;

        if ( _regions_match_alleles(files->targets, files->targets_als-1, files->readers[i].buffer[0]) ) return ret;
        
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

size_t mygetline(char **line, size_t *n, FILE *fp)
{
    if (line == NULL || n == NULL || fp == NULL)
    {
        errno = EINVAL;
        return -1;
    }
    if (*n==0 || !*line)
    {
        *line = NULL;
        *n = 0;
    }

    size_t nread=0;
    int c;
    while ((c=getc(fp))!= EOF && c!='\n')
    {
        if ( ++nread>=*n )
        {
            *n += 255;
            *line = (char*) realloc(*line, sizeof(char)*(*n));
        }
        (*line)[nread-1] = c;
    }
    if ( nread>=*n )
    {
        *n += 255;
        *line = (char*) realloc(*line, sizeof(char)*(*n));
    }
    (*line)[nread] = 0;
    return nread>0 ? nread : -1;

}

int bcf_sr_set_samples(bcf_srs_t *files, const char *fname)
{
    int i;
    struct stat sbuf;
    files->samples = NULL;
    files->n_smpl  = 0;

    int *exclude = NULL;
    if ( fname[0]=='!' )    // exclude samples
    { 
        fname++; 
        exclude = (int*) calloc(bcf_hdr_nsamples(files->readers[0].header),sizeof(int));
        if ( stat(fname, &sbuf)==0 && S_ISREG(sbuf.st_mode) ) // reading from file
        {
            FILE *fp = fopen(fname,"r");
            if ( !fp ) { fprintf(stderr,"%s: %s\n", fname,strerror(errno)); return 0; }
            char *line = NULL;
            size_t len = 0;
            ssize_t nread;
            while ((nread = mygetline(&line, &len, fp)) != -1)
            {
                int id = bcf_hdr_id2int(files->readers[0].header, BCF_DT_SAMPLE, line);
                // sanity check our assumptions for the things to follow
                assert( id < bcf_hdr_nsamples(files->readers[0].header) );
                assert( !strcmp(files->readers[0].header->samples[id],line) );
                if ( id>=0 ) exclude[id] = 1;
            }
        }
        else
        {
            kstring_t str = {0,0,0};
            const char *b = fname;
            while (b)
            {
                str.l = 0;
                const char *e = index(b,','); 
                if ( !(e-b) ) break;
                if ( e ) { kputsn(b, e-b, &str); e++; }
                else kputs(b, &str);
                b = e;
                for (i=0; i<files->nreaders; i++)
                {
                    int id = bcf_hdr_id2int(files->readers[0].header, BCF_DT_SAMPLE, str.s);
                    // sanity check our assumptions for the things to follow
                    assert( id < bcf_hdr_nsamples(files->readers[0].header) );
                    assert( !strcmp(files->readers[0].header->samples[id],str.s) );
                    if ( id>=0 ) exclude[id] = 1;
                }
            }
            free(str.s);
        }
        fname = "-";
    }
    if ( !strcmp(fname,"-") )   // Intersection of all samples across all readers
    {
        int n = files->readers[0].header->n[BCF_DT_SAMPLE];
        char **smpl = files->readers[0].header->samples;
        int ism;
        for (ism=0; ism<n; ism++)
        {
            if ( exclude && exclude[ism] ) continue;
            int n_isec = 1;
            for (i=1; i<files->nreaders; i++)
            {
                if ( bcf_hdr_id2int(files->readers[i].header, BCF_DT_SAMPLE, smpl[ism])==-1 ) break;
                n_isec++;
            }
            if ( n_isec<files->nreaders ) continue;
            files->samples = (char**) realloc(files->samples, (files->n_smpl+1)*sizeof(const char*));
            files->samples[files->n_smpl++] = strdup(smpl[ism]);
        }
    }
    else if ( stat(fname, &sbuf)==0 && S_ISREG(sbuf.st_mode) ) // read samples from file
    {
        FILE *fp = fopen(fname,"r");
        if ( !fp ) { free(exclude); fprintf(stderr,"%s: %s\n", fname,strerror(errno)); return 0; }
        char *line = NULL;
        size_t len = 0;
        ssize_t nread;
        while ((nread = mygetline(&line, &len, fp)) != -1) 
        {
            int n_isec = 0;
            for (i=0; i<files->nreaders; i++)
            {
                if ( bcf_hdr_id2int(files->readers[i].header, BCF_DT_SAMPLE, line)==-1 ) break;
                n_isec++;
            }
            if ( n_isec<files->nreaders ) 
            {
                fprintf(stderr,"[init_samples] sample not found, skipping: [%s]\n", line);
                continue;
            }
            files->samples = (char**) realloc(files->samples, (files->n_smpl+1)*sizeof(const char*));
            files->samples[files->n_smpl++] = strdup(line);
        }
        if (line) free(line);
        fclose(fp);
    }
    else    // samples given as a comma-separated list
    {
        kstring_t str = {0,0,0};
        const char *b = fname;
        while (b)
        {
            str.l = 0;
            const char *e = index(b,','); 
            if ( !(e-b) ) break;
            if ( e ) { kputsn(b, e-b, &str); e++; }
            else kputs(b, &str);
            b = e;

            int n_isec = 0;
            for (i=0; i<files->nreaders; i++)
            {
                if ( bcf_hdr_id2int(files->readers[i].header, BCF_DT_SAMPLE, str.s)==-1 ) break;
                n_isec++;
            }
            if ( n_isec<files->nreaders ) 
            {
                fprintf(stderr,"[init_samples] sample not found, skipping: %s\n", str.s);
                continue;
            }
            files->samples = (char**) realloc(files->samples, (files->n_smpl+1)*sizeof(const char*));
            files->samples[files->n_smpl++] = strdup(str.s);
        }
        if ( str.s ) free(str.s);
    }
    free(exclude);
    if ( !files->n_smpl ) 
    {
        if ( files->nreaders>1 ) fprintf(stderr,"[init_samples] No samples in common.\n");
        return 0;
    }
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        reader->samples  = (int*) malloc(sizeof(int)*files->n_smpl);
        reader->n_smpl   = files->n_smpl;
        int ism;
        for (ism=0; ism<files->n_smpl; ism++)
            reader->samples[ism] = bcf_hdr_id2int(reader->header, BCF_DT_SAMPLE, files->samples[ism]);
    }
    return 1;
}


// Add a new region into a list sorted by start,end (1-based coordinates)
static void _regions_add(bcf_sr_regions_t *reg, const char *chr, int start, int end)
{
    if ( start==-1 && end==-1 )
    {
        start = 0; end = (1<<29) - 1;
    }
    else
    {
        start--; end--; // store 0-based coordinates
    }

    int i;
    for (i=0; i<reg->nregs; i++)
        if ( !strcmp(reg->regs[i].chr,chr) ) break;
    if ( i<reg->nregs && !strcmp(reg->regs[i].chr,chr) ) // the chromosome block already exists
    {
        for (; i<reg->nregs; i++)
            if ( strcmp(reg->regs[i].chr,chr) || reg->regs[i].start >= start ) break;

        // return if the region already exists
        if ( i<reg->nregs && !strcmp(reg->regs[i].chr,chr) && reg->regs[i].start==start && reg->regs[i].end==end ) return;

        reg->regs = (region_t*) realloc(reg->regs,sizeof(region_t)*(reg->nregs+1));
        if ( i<reg->nregs )
            memmove(&reg->regs[i+1],&reg->regs[i],(reg->nregs - i)*sizeof(region_t));
    }
    else
        reg->regs = (region_t*) realloc(reg->regs,sizeof(region_t)*(reg->nregs+1));

    // Check if a new sequence name has to be added
    int j;
    for (j=0; j<reg->nseqs; j++)
        if ( !strcmp(chr,reg->snames[j]) ) break;
    if ( j==reg->nseqs )
    {
        reg->nseqs++;
        reg->snames = (char**) realloc(reg->snames,sizeof(char*)*reg->nseqs);
        reg->snames[j] = strdup(chr);
    }

    reg->nregs++;
    reg->regs[i].chr   = reg->snames[j];
    reg->regs[i].start = start;
    reg->regs[i].end   = end;
}

// File name or a list of genomic locations
static bcf_sr_regions_t *_regions_init_string(const char *str)
{
    struct stat sbuf;
    if ( stat(str, &sbuf)==0 ) return NULL; // it's a file

    bcf_sr_regions_t *reg = (bcf_sr_regions_t *) calloc(1, sizeof(bcf_sr_regions_t));

    kstring_t tmp = {0,0,0};
    const char *sp = str, *ep = str;
    int from, to;
    while ( 1 )
    {
        while ( *ep && *ep!=',' && *ep!=':' ) ep++;
        tmp.l = 0;
        kputsn(sp,ep-sp,&tmp);
        if ( *ep==':' )
        {
            sp = ep+1;
            from = strtol(sp,(char**)&ep,10);
            if ( sp==ep )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            if ( !*ep || *ep==',' )
            {
                _regions_add(reg, tmp.s, from, from);
                sp = ep;
                continue;
            }
            if ( *ep!='-' )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            ep++;
            sp = ep;
            to = strtol(sp,(char**)&ep,10);
            if ( *ep && *ep!=',' )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            if ( sp==ep ) to = 1<<29;
            _regions_add(reg, tmp.s, from, to);
            if ( !*ep ) break;
            sp = ep;
        }
        else
        {
            if ( tmp.l ) _regions_add(reg, tmp.s, -1, -1);
            if ( !*ep ) break;
            sp = ++ep;
        }
    }
    free(tmp.s);
    reg->ireg = -1;
    return reg;
}

// ichr,ifrom,ito are 0-based; line will be modified so that the *chr pointer is 0-terminated
// returns -1 on error, 0 if the line is a comment line, 1 on success
static int _regions_parse_line(char *line, int ichr,int ifrom,int ito, char **chr,int *from,int *to)
{
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
        *from = *to = strtol(ss, &tmp, 10);
        if ( tmp==ss ) return -1;
    }
    else
    {
        if ( k==ifrom ) 
            *from = strtol(ss, &tmp, 10);
        else
            *to = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;

        for (i=k; i<l && *se; i++)
        {
            ss = ++se;
            while (*se && *se!='\t') se++;
        }
        if ( i<l ) return -1;
        if ( k==ifrom ) 
            *to = strtol(ss, &tmp, 10);
        else
            *from = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;
    }

    ss = se = line;
    for (i=0; i<=ichr && *se; i++)
    {
        if ( i>0 ) ss = ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=ichr ) return -1;
    *se = 0;
    *chr = ss;
    return 1;
}

bcf_sr_regions_t *bcf_sr_regions_init(const char *regions)
{
    bcf_sr_regions_t *reg = _regions_init_string(regions);   // file name or a genomic region?
    if ( reg ) return reg;

    reg = (bcf_sr_regions_t *) calloc(1, sizeof(bcf_sr_regions_t));
    reg->ireg = -1;

    reg->file = hts_open(regions, "rb");
    if ( !reg->file )
    {
        fprintf(stderr,"[%s:%d %s] Could not open file: %s\n", __FILE__,__LINE__,__FUNCTION__,regions);
        free(reg);
        return NULL;
    }

    reg->tbx = tbx_index_load(regions);
    if ( !reg->tbx ) 
    {
        int ichr = 0, ifrom = 1, ito = 2;
        int len = strlen(regions);
        int is_bed  = strcasecmp(".bed",regions+len-4) ? 0 : 1;
        if ( !is_bed && !strcasecmp(".bed.gz",regions+len-7) ) is_bed = 1;
        int ft_type = hts_file_type(regions);
        if ( ft_type & FT_VCF ) ito = 1;

        // read the whole file, tabix index is not present
        while ( hts_getline(reg->file, KS_SEP_LINE, &reg->line) > 0 )
        {
            char *chr;
            int from, to, ret;
            ret = _regions_parse_line(reg->line.s, ichr,ifrom,ito, &chr,&from,&to);
            if ( ret < 0 ) 
            {
                fprintf(stderr,"[%s:%d] Could not parse the file %s, using the columns %d,%d,%d\n", __FILE__,__LINE__,regions,ichr+1,ifrom+1,ito+1);
                hts_close(reg->file); reg->file = NULL; free(reg); 
                return NULL;
            }
            if ( !ret ) continue;
            if ( is_bed ) from++;
            _regions_add(reg, chr, from, to);
        }
        hts_close(reg->file); reg->file = NULL;
        if ( !reg->nregs ) { free(reg); return NULL; }
        return reg;
    }
    reg->snames = (char**) tbx_seqnames(reg->tbx, &reg->nseqs);
    reg->fname  = strdup(regions);
    reg->is_bin = 1;
    return reg;
}

void bcf_sr_regions_destroy(bcf_sr_regions_t *reg)
{
    int i;
    free(reg->regs);
    free(reg->fname);
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    if ( reg->tbx ) tbx_destroy(reg->tbx);
    if ( reg->file ) hts_close(reg->file);
    if ( reg->als ) free(reg->als);
    free(reg->line.s);
    if (reg->regs) 
        for (i=0; i<reg->nseqs; i++) free(reg->snames[i]);  // free only in-memory names, tbx names are const
    free(reg->snames);
    free(reg);
}

int bcf_sr_regions_seek(bcf_sr_regions_t *reg, const char *seq)
{
    reg->done = 1;
    reg->start = reg->end = -1;

    int i;
    if ( reg->regs )    // using in-memory regions
    {
        for (i=0; i<reg->nregs; i++)
            if ( !strcmp(seq,reg->regs[i].chr) ) break;
        reg->ireg = i-1;
        if ( i==reg->nregs ) return -1;
        reg->seq  = reg->snames[i]; 
        reg->done = 0;
        return 0;
    }

    // reading regions from tabix
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    reg->itr = tbx_itr_querys(reg->tbx, seq);
    if ( reg->itr )
    {
        for (i=0; i<reg->nseqs; i++) 
            if (!strcmp(seq,reg->snames[i]) ) break;
        if ( i==reg->nseqs ) return -1;
        reg->seq  = reg->snames[i];
        reg->done = 0;
        return 0;
    }
    return -1;
}

int bcf_sr_regions_next(bcf_sr_regions_t *reg)
{
    if ( reg->done ) return -1;
    reg->seq  = NULL; reg->start = reg->end = -1;
    reg->nals = 0;

    if ( reg->regs )    // using in-memory regions
    {
        reg->ireg++;
        if ( reg->ireg>=reg->nregs ) { reg->done = 1; return -1; } // no more regions left
        reg->seq   = reg->regs[reg->ireg].chr;
        reg->start = reg->regs[reg->ireg].start;
        reg->end   = reg->regs[reg->ireg].end;
        return 0;
    }

    char *chr;
    int ichr = 0, ifrom = 1, ito = 2, is_bed = 0, from, to;
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
            if ( ret<0 ) { reg->done = 1; return -1; }
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
                    fprintf(stderr,"[%s:%d %s] Could not open file: %s\n", __FILE__,__LINE__,__FUNCTION__,reg->fname);
                    reg->file = NULL;
                    bcf_sr_regions_destroy(reg);
                    return -1;
                }
                reg->is_bin = 0;
            }

            // tabix index absent, reading the whole file
            ret = hts_getline(reg->file, KS_SEP_LINE, &reg->line);
            if ( ret<0 ) { reg->done = 1; return -1; }
        }
        ret = _regions_parse_line(reg->line.s, ichr,ifrom,ito, &chr,&from,&to);
        if ( ret<0 ) 
        {
            fprintf(stderr,"[%s:%d] Could not parse the file %s, using the columns %d,%d,%d\n", __FILE__,__LINE__,reg->fname,ichr+1,ifrom+1,ito+1);
            return -1;
        }
    }
    if ( is_bed ) from++;

    // find the chromosome name: assuming small number of chromosomes!
    int i;
    for (i=0; i<reg->nseqs; i++) 
        if (!strcmp(chr,reg->snames[i]) ) break;
    if ( !(i<reg->nseqs) ) fprintf(stderr,"i=%d nseq=%d  [%s][%s] [%s]\n", i,reg->nseqs,chr,reg->line.s,reg->snames[0]);
    assert( i<reg->nseqs );

    // This is a bit hacky: unset the chr-terminating 0 set by _regions_parse_line, or
    //  otherwise _regions_match_alleles will be confused.
    int len = strlen(chr);
    if ( len < reg->line.l ) chr[len] = '\t';

    reg->seq    = reg->snames[i];
    reg->start  = from - 1;
    reg->end    = to - 1;
    return 0;
}

static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec)
{
    int i = 0, max_len = 0;
    if ( !reg->nals )
    {
        char *ss = reg->line.s;
        while ( i<als_idx && *ss )
        {
            if ( *ss=='\t' ) i++;
            ss++;
        }
        reg->nals = 1;
        hts_expand(char*,reg->nals,reg->mals,reg->als);
        reg->als[0] = ss;
        while ( *(++ss) )
        {
            if ( *ss=='\t' ) break;
            if ( *ss!=',' ) continue;
            *ss = 0;
            reg->nals++;
            hts_expand(char*,reg->nals,reg->mals,reg->als);
            reg->als[reg->nals-1] = ss+1;
            if ( ss - reg->als[reg->nals-2] > max_len ) max_len = ss - reg->als[reg->nals-2];
        }
        if ( ss - reg->als[reg->nals-1] > max_len ) max_len = ss - reg->als[reg->nals-1];
        reg->als_type = max_len > 1 ? VCF_INDEL : VCF_SNP;  // this is a too-simplified check, see vcf.c:bcf_set_variant_types
    }
    int type = bcf_get_variant_types(rec);
    if ( reg->als_type & VCF_INDEL )
        return type & VCF_INDEL ? 1 : 0;
    return !(type & VCF_INDEL) ? 1 : 0;
}

int bcf_sr_regions_query(bcf_sr_regions_t *reg, int start, int end)
{
    if ( reg->done ) return -2;     // no more regions left

    // init regions if it was not done already
    if ( reg->start==-1 )
        if ( bcf_sr_regions_next(reg) < 0 ) return -2;  // no more regions left

    char *seq = reg->seq;
    while ( reg->seq==seq && reg->end < start )
    {
        if ( bcf_sr_regions_next(reg) < 0 ) return -2;  // no more regions left
    }
    if ( reg->seq != seq ) return -2;   // different chromosome, new seek is necessary
    if ( reg->start <= end ) return 0;  // is in a region
    return -1;  // does not overlap any region 
}

