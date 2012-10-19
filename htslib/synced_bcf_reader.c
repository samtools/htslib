#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>
#include "synced_bcf_reader.h"

int bcf_sr_add_reader(readers_t *files, const char *fname)
{
    files->readers = (reader_t*) realloc(files->readers, sizeof(reader_t)*(files->nreaders+1));
    reader_t *reader = &files->readers[files->nreaders++];
    memset(reader,0,sizeof(reader_t));

    int type = file_type(fname);
    if ( type==IS_VCF_GZ ) 
    {
        reader->tbx = tbx_index_load(fname);
        if ( !reader->tbx )
        {
            fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
            return 0;
        }

        // This is just to read the header
        htsFile *file = hts_open(fname, "r", NULL);
        if ( !file ) return 0;
        reader->header = vcf_hdr_read(file);
        hts_close(file);

        // The VCF opened in binary tabix mode
        reader->file = hts_open(fname, "rb", NULL);
        if ( !reader->file ) return 0;
    }
    else if ( type==IS_BCF ) 
    {
        reader->file = hts_open(fname, "rb", NULL);
        if ( !reader->file ) return 0;
        reader->header = vcf_hdr_read(reader->file);

        reader->bcf = bcf_index_load(fname);
        if ( !reader->bcf ) 
        {
            fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
            return 0;   // not indexed..?
        }
    }
    else
    {
        fprintf(stderr,"Expected .vcf.gz or .bcf file\n");
        return 0;
    }

    reader->fname   = fname;
    reader->filter_id = -1;
    if ( files->apply_filters )
        reader->filter_id = bcf_id2int(reader->header, BCF_DT_ID, "PASS");

    // Update list of chromosomes
    if ( files->region )
    {
        if ( !files->seqs )
        {
            files->mseqs = files->nseqs = 1;
            files->seqs = (const char**) malloc(sizeof(const char*));
            files->seqs[0] = files->region;
        }
    }
    else
    {
        int n,i,j;
        const char **names = bcf_seqnames(reader->header, &n);
        for (i=0; i<n; i++)
        {
            for (j=0; j<files->nseqs; j++)
                if ( !strcmp(names[i],files->seqs[j]) ) break;
            if ( j<files->nseqs ) continue;     // already have this chr
            files->mseqs += 30;
            files->seqs = (const char**) realloc(files->seqs, sizeof(const char*)*files->mseqs);
            files->seqs[files->nseqs++] = names[i];
        }
        free(names);
    }
    files->iseq = -1;
    return 1;
}

readers_t *bcf_sr_init(void)
{
    readers_t *files = (readers_t*) calloc(1,sizeof(readers_t));
    return files;
}

void bcf_sr_destroy(readers_t *files)
{
    if ( !files->nreaders ) return;
    int i;
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        if ( reader->tbx ) tbx_destroy(reader->tbx);
        if ( reader->bcf ) hts_idx_destroy(reader->bcf);
        bcf_hdr_destroy(reader->header);
        hts_close(reader->file);
        if ( reader->itr ) tbx_itr_destroy(reader->itr);
        int j;
        for (j=0; j<reader->mbuffer; j++)
            bcf_destroy1(reader->buffer[j]);
        free(reader->buffer);
        if ( reader->samples ) free(reader->samples);
    }
    free(files->readers);
    free(files->seqs);
    for (i=0; i<files->n_smpl; i++) free(files->samples[i]);
    free(files->samples);
    if (files->targets)
    {
        if (files->targets->itr) tbx_itr_destroy(files->targets->itr);
        tbx_destroy(files->targets->tbx);
        if (files->targets->line.m) free(files->targets->line.s);
        hts_close(files->targets->file);
        free(files->targets->seq_names);
        free(files->targets);
    }
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
static void collapse_buffer(readers_t *files, reader_t *reader)
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
        if ( files->collapse&COLLAPSE_SNPS && line->d.var_type&(VCF_SNP|VCF_MNP) )
        {
            if ( !has_snp ) has_snp = 1;
            else line->pos = -1;
        }
        if ( files->collapse&COLLAPSE_INDELS && line->d.var_type&VCF_INDEL )
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

void debug_buffer(FILE *fp, reader_t *reader)
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

void debug_buffers(FILE *fp, readers_t *files)
{
    int i;
    for (i=0; i<files->nreaders; i++)
        debug_buffer(fp, &files->readers[i]);
    fprintf(fp,"\n");
}

int bcf_sr_set_targets(readers_t *files, const char *fname)
{
    regions_t *tgts = (regions_t *) calloc(1,sizeof(regions_t));
    tgts->file = hts_open(fname, "rb", NULL);
    if ( !tgts->file ) return 0;
    tgts->tbx = tbx_index_load(fname);
    tgts->seq_names = (char**) tbx_seqnames(tgts->tbx, &tgts->nseqs);
    tgts->cseq = -1;
    files->targets = tgts;
    return 1;
}

static char *tgt_next_seq(regions_t *tgt)
{
    if ( ++tgt->cseq >= tgt->nseqs ) return NULL;
    if ( tgt->itr ) tbx_itr_destroy(tgt->itr);
    tgt->itr = tbx_itr_querys(tgt->tbx,tgt->seq_names[tgt->cseq]);
    tgt->tpos.to = -1;
    return tgt->seq_names[tgt->cseq];
}

// 1 if position is present, 0 if not, -1 if not and done
static int tgt_has_position(regions_t *tgt, int32_t pos)
{
    while ( tgt->tpos.to < pos )
    {
        int ret = tbx_itr_next((BGZF*)tgt->file->fp, tgt->tbx, tgt->itr, &tgt->line);
        if ( ret<0 ) return -1;

        // parse line
        int k,l;
        if ( tgt->tbx->conf.bc <= tgt->tbx->conf.ec ) 
            k = tgt->tbx->conf.bc, l = tgt->tbx->conf.ec;
        else 
            l = tgt->tbx->conf.bc, k = tgt->tbx->conf.ec;

        int i = 0;
        char *end = tgt->line.s, *start = NULL;
        for (i=0; i<k; i++)
        {
            start = i==0 ? end++ : ++end;
            while (*end && *end!='\t') end++;
        }
        if ( k==l )
            tgt->tpos.from = tgt->tpos.to = strtol(start, NULL, 10);
        else
        {
            if ( k==tgt->tbx->conf.bc ) 
                tgt->tpos.from = strtol(start, NULL, 10);
            else
                tgt->tpos.to = strtol(start, NULL, 10);

            for (i=k; i<l; i++)
            {
                start = ++end;
                while (*end && *end!='\t') end++;
            }
            if ( k==tgt->tbx->conf.bc ) 
                tgt->tpos.to = strtol(start, NULL, 10);
            else
                tgt->tpos.from = strtol(start, NULL, 10);
        }
        tgt->tpos.from--;
        tgt->tpos.to--;
    }

    if ( pos >= tgt->tpos.from && pos <= tgt->tpos.to ) return 1;
    return 0;
}

int bcf_sr_next_line(readers_t *files)
{
    int32_t min_pos = INT_MAX;
    int ret,i,j;
    kstring_t *str = &files->tmps;

    while ( min_pos==INT_MAX )
    {
        // Need to open new chromosome?
        int eos = 0;
        for (i=0; i<files->nreaders; i++)
            if ( !files->readers[i].itr && !files->readers[i].nbuffer ) eos++;
        if ( eos==files->nreaders )
        {
            const char *seq;
            if ( files->targets )
            {
                seq = tgt_next_seq(files->targets);
                if ( !seq ) return 0;   // all chroms scanned
            }
            else
            {
                if ( ++files->iseq >= files->nseqs ) return 0;  // all chroms scanned
                seq = files->seqs[files->iseq];
            }
            for (i=0; i<files->nreaders; i++)
            {
                reader_t *reader = &files->readers[i];
                if ( reader->tbx )
                    reader->itr = tbx_itr_querys(reader->tbx,seq);
                else
                    reader->itr = bcf_itr_querys(reader->bcf,reader->header,seq);
            }
        }

        // Find the smallest coordinate
        for (i=0; i<files->nreaders; i++)
        {
            reader_t *reader = &files->readers[i];
            int buffer_full = ( reader->nbuffer && reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) ? 1 : 0;
            if ( reader->itr && !buffer_full )
            {
                // Fill the buffer with records starting at the same position
                while (1)
                {
                    if ( reader->nbuffer+1 >= reader->mbuffer ) 
                    {
                        reader->mbuffer += 8;
                        reader->buffer = (bcf1_t**) realloc(reader->buffer, sizeof(bcf1_t*)*reader->mbuffer);
                        for (j=8; j>0; j--)
                            reader->buffer[reader->mbuffer-j] = bcf_init1();
                    }
                    if ( reader->tbx )
                    {
                        ret = tbx_itr_next((BGZF*)reader->file->fp, reader->tbx, reader->itr, str);
                        if ( ret<0 ) break;
                        vcf_parse1(str, reader->header, reader->buffer[reader->nbuffer+1]);
                    }
                    else
                    {
                        ret = bcf_itr_next((BGZF*)reader->file->fp, reader->itr, reader->buffer[reader->nbuffer+1]);
                        if ( ret<0 ) break;
                    }
                    bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR|BCF_UN_FLT);
                    // apply filter
                    if ( reader->filter_id!=-1 && reader->buffer[reader->nbuffer+1]->d.n_flt && reader->filter_id!=reader->buffer[reader->nbuffer+1]->d.flt[0] ) continue;
                    set_variant_types(reader->buffer[reader->nbuffer+1]);
                    reader->nbuffer++;
                    if ( reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) break;
                }
                if ( ret<0 ) { tbx_itr_destroy(reader->itr); reader->itr = NULL; } // done for this chromosome
            }
            if ( reader->nbuffer )
            {
                if ( min_pos > reader->buffer[1]->pos ) min_pos = reader->buffer[1]->pos; 
            }
            // The buffer is full - either there is nothing else to read or the last record has a different coordinate
            if ( files->collapse && reader->nbuffer>2 && reader->buffer[1]->pos==reader->buffer[2]->pos )
            {
                collapse_buffer(files, reader);
            }
        }
        if ( files->targets && min_pos!=INT_MAX )
        {
            int ret = tgt_has_position(files->targets, min_pos);
            if ( ret==1 ) continue;

            // The position must be skipped
            if ( ret==-1 )
            {
                // done for this chromosome, don't read the rest
                for (i=0; i<files->nreaders; i++) 
                {
                    files->readers[i].nbuffer = 0;
                    if ( files->readers[i].itr )
                    {
                        tbx_itr_destroy(files->readers[i].itr);
                        files->readers[i].itr = NULL;
                    }
                }
                min_pos = INT_MAX;
                continue;
            }

            // remove the active line, save the buffer line
            for (i=0; i<files->nreaders; i++)
            {
                reader_t *reader = &files->readers[i];
                for (j=1; j<=reader->nbuffer; j++)
                    if ( reader->buffer[j]->pos!=min_pos ) break;
                if ( j==1 ) continue;
                if ( j<=reader->nbuffer )
                {
                    bcf1_t *tmp = reader->buffer[1]; reader->buffer[1] = reader->buffer[j]; reader->buffer[j] = tmp;
                    reader->nbuffer = 1;
                }
                else 
                    reader->nbuffer = 0;
            }
            min_pos = INT_MAX;
        }
    }

    //printf("[next_line] min_pos=%d\n", min_pos+1);
    //debug_buffers(files);

    // Set the current line
    ret = 0;
    bcf1_t *first = NULL;
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        if ( !reader->nbuffer || reader->buffer[1]->pos!=min_pos ) continue;

        // Match the records by REF and ALT
        int j, irec = -1;
        if ( first )
        {
            for (j=1; j<=reader->nbuffer; j++)
            {
                bcf1_t *line = reader->buffer[j];
                if ( min_pos != line->pos ) break;  // done with this buffer

                if ( files->collapse&COLLAPSE_ANY ) { irec=j; break; }  // checking position only
                if ( files->collapse&COLLAPSE_SNPS && first->d.var_type&VCF_SNP && line->d.var_type&VCF_SNP ) { irec=j; break; }
                if ( files->collapse&COLLAPSE_INDELS && first->d.var_type&VCF_INDEL && line->d.var_type&VCF_INDEL ) { irec=j; break; }

                if ( first->rlen != line->rlen ) continue;  // REFs do not match
                if ( strcmp(first->d.allele[0], line->d.allele[0]) ) continue; // REFs do not match
                int ial,jal;
                if ( files->collapse==COLLAPSE_NONE )
                {
                    // require exact match, all alleles must be identical
                    if ( first->n_allele!=line->n_allele ) continue;   // different number of alleles
                    int nmatch = 1; // REF has been already checked
                    for (ial=1; ial<first->n_allele; ial++)
                    {
                        for (jal=1; jal<line->n_allele; jal++)
                            if ( !strcmp(first->d.allele[ial], line->d.allele[jal]) ) { nmatch++; break; }
                    }
                    if ( nmatch>=first->n_allele ) { irec=j; break; }
                }
                else
                {
                    // thorough check: the REFs and some of the alleles have to be shared
                    // (neglecting different representations of the same indel for now)
                    for (ial=1; ial<first->n_allele; ial++)
                    {
                        for (jal=1; jal<line->n_allele; jal++)
                            if ( !strcmp(first->d.allele[ial], line->d.allele[jal]) ) { irec=j; break; }
                        if ( irec>=1 ) break;
                    }
                }
                if ( irec>=1 ) break;
            }
            if ( irec==-1 ) continue;
        }
        else 
        {
            first = reader->buffer[1];
            irec  = 1;
        }
        bcf1_t *tmp = reader->buffer[0];
        reader->buffer[0] = reader->buffer[irec];
        for (j=irec+1; j<=reader->nbuffer; j++)
            reader->buffer[j-1] = reader->buffer[j];
        reader->buffer[ reader->nbuffer ] = tmp;
        reader->nbuffer--;
        ret |= 1<<i;
    }
    // fprintf(stdout,"[next_line] min_pos=%d mask=%d\n", min_pos+1, ret);
    // debug_buffers(stdout,files);

    return ret;
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

int bcf_sr_set_samples(readers_t *files, const char *fname)
{
    int i;
    struct stat sbuf;
    files->samples = NULL;
    files->n_smpl  = 0;
    if ( !strcmp(fname,"-") )   // Intersection of all samples across all readers
    {
        int n = files->readers[0].header->n[BCF_DT_SAMPLE];
        char **smpl = files->readers[0].header->samples;
        int ism;
        for (ism=0; ism<n; ism++)
        {
            int n_isec = 1;
            for (i=1; i<files->nreaders; i++)
            {
                if ( bcf_id2int(files->readers[i].header, BCF_DT_SAMPLE, smpl[ism])==-1 ) break;
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
        if ( !fp ) { fprintf(stderr,"%s: %s\n", fname,strerror(errno)); return 0; }
        char *line = NULL;
        size_t len = 0;
        ssize_t nread;
        while ((nread = mygetline(&line, &len, fp)) != -1) 
        {
            int n_isec = 0;
            for (i=0; i<files->nreaders; i++)
            {
                if ( bcf_id2int(files->readers[i].header, BCF_DT_SAMPLE, line)==-1 ) break;
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
                if ( bcf_id2int(files->readers[i].header, BCF_DT_SAMPLE, str.s)==-1 ) break;
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
    if ( !files->n_smpl ) 
    {
        if ( files->nreaders>1 ) fprintf(stderr,"[init_samples] No samples in common.\n");
        return 0;
    }
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        reader->samples  = (int*) malloc(sizeof(int)*files->n_smpl);
        reader->n_smpl   = files->n_smpl;
        int ism;
        for (ism=0; ism<files->n_smpl; ism++)
            reader->samples[ism] = bcf_id2int(reader->header, BCF_DT_SAMPLE, files->samples[ism]);
    }
    return 1;
}

int init_regions(const char *fname, regions_t *reg)
{
    int bgzf_getline(BGZF *fp, int delim, kstring_t *str);

    BGZF *zfp = bgzf_open(fname, "r");
    if ( !zfp ) 
    {
        fprintf(stderr,"%s: %s\n",fname,strerror(errno));
        return 0;
    }

    int i, mseqs = 10, mpos = 0;
    reg->nseqs = 0;
    reg->pos   = (pos_t **)calloc(mseqs,sizeof(pos_t*));
    reg->npos  = (int*) calloc(mseqs,sizeof(int));
    reg->seq_names = (char **) calloc(mseqs,sizeof(char*));

    kstring_t str = {0,0,0};
    ssize_t nread;
    while ((nread = bgzf_getline(zfp, '\n', &str)) > 0) 
    {
        char *line = str.s;
        if ( line[0] == '#' ) continue;

        int i = 0;
        while ( i<nread && !isspace(line[i]) ) i++;
        if ( i>=nread ) 
        { 
            fprintf(stderr,"Could not parse the file: %s [%s]\n", fname,line); 
            return 0; 
        }
        line[i] = 0;

        if ( reg->nseqs==0 || strcmp(line,reg->seq_names[reg->nseqs-1]) )
        {
            // New sequence
            reg->nseqs++;
            if ( reg->nseqs >= mseqs )
            {
                mseqs++;
                reg->pos  = (pos_t **) realloc(reg->pos,sizeof(pos_t*)*mseqs); reg->pos[mseqs-1] = NULL;
                reg->npos = (int *) realloc(reg->npos,sizeof(int)*mseqs); reg->npos[mseqs-1] = 0;
                reg->seq_names = (char**) realloc(reg->seq_names,sizeof(char*)*mseqs);
            }
            reg->seq_names[reg->nseqs-1] = strdup(line);
            mpos = 0;
        }

        int iseq = reg->nseqs-1;
        if ( reg->npos[iseq] >= mpos )
        {
            mpos += 100;
            reg->pos[iseq] = (pos_t*) realloc(reg->pos[iseq],sizeof(pos_t)*mpos);
        }
        int ipos = reg->npos[iseq];
        pos_t *pos = reg->pos[iseq];
        reg->npos[iseq]++;
        if ( (sscanf(line+i+1,"%d %d",&pos[ipos].from,&pos[ipos].to))!=2 ) 
        {
            if ( (sscanf(line+i+1,"%d",&pos[ipos].from))!=1 )
            {
                fprintf(stderr,"Could not parse the region [%s]\n",line+i+1);
                return 0;
            }
            pos[ipos].to = pos[ipos].from;
        }

        // Check that the file is sorted
        if ( ipos>0 && (pos[ipos].from < pos[ipos-1].from || (pos[ipos].from==pos[ipos-1].from && pos[ipos].to<pos[ipos-1].to)) )
        {
            fprintf(stderr,"The file is not sorted: %s\n", fname);
            return 0;
        }
    }

    // Check that chromosomes come in blocks
    int j;
    for (i=0; i<reg->nseqs; i++)
    {
        for (j=0; j<i; j++)
        {
            if ( !strcmp(reg->seq_names[i],reg->seq_names[j]) ) 
            {
                fprintf(stderr,"The file is not sorted: %s\n", fname);
                return 0;
            }
        }
    }

    if (str.m) free(str.s);
    else return 0;

    bgzf_close(zfp);
    return 1;
}

void destroy_regions(regions_t *reg)
{
    int i;
    for (i=0; i<reg->nseqs; i++)
    {
        free(reg->pos[i]);
        free(reg->seq_names[i]);
    }
    free(reg->seq_names);
    free(reg->pos);
    free(reg->npos);
}

int reset_regions(regions_t *reg, const char *seq)
{
    reg->cpos = 0;
    reg->cseq = -1;
    int i;
    for (i=0; i<reg->nseqs; i++)
    {
        int n = strlen(reg->seq_names[i]);
        if ( strncmp(reg->seq_names[i],seq,n) ) continue;
        if ( seq[n] && seq[n]!=':' ) continue;
        reg->cseq = i;
        return 1;
    }
    return 0;
}

pos_t *is_in_regions(regions_t *reg, int32_t pos)
{
    if ( reg->cseq==-1 ) return NULL;

    int ipos = reg->cpos;
    int npos = reg->npos[reg->cseq];
    if ( ipos==npos ) return NULL;  // done for this chr

    pos_t *p = reg->pos[reg->cseq];

    // Find a matching interval
    while ( ipos < npos && pos > p[ipos].to ) ipos++;
    reg->cpos = ipos;
    if ( ipos >= npos || pos < p[ipos].from ) return NULL;

    return &p[ipos];
}



