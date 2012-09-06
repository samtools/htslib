#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>
#include "synced_bcf_reader.h"

#define VCF    1
#define VCF_GZ 2
#define BCF    3

int file_type(const char *fname)
{
    int len = strlen(fname);
    if ( !strcasecmp(".vcf.gz",fname+len-7) ) return VCF_GZ;
    if ( !strcasecmp(".vcf",fname+len-4) ) return VCF;
    if ( !strcasecmp(".bcf",fname+len-4) ) return BCF;
    return 0;
}

int add_reader(const char *fname, readers_t *files)
{
    files->readers = (reader_t*) realloc(files->readers, sizeof(reader_t)*(files->nreaders+1));
    reader_t *reader = &files->readers[files->nreaders++];
    memset(reader,0,sizeof(reader_t));

    int type = file_type(fname);
    if ( type==VCF_GZ ) 
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
    else if ( type==BCF ) 
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

    reader->line  = bcf_init1();
    reader->fname = fname;
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
        const char **names = bcf_seq_names(reader->header, &n);
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

void destroy_readers(readers_t *files)
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
        bcf_destroy1(reader->line);
        free(reader->buffer);
        if ( reader->samples ) free(reader->samples);
    }
    free(files->readers);
    free(files->seqs);
    for (i=0; i<files->n_smpl; i++) free(files->samples[i]);
    free(files->samples);
    memset(files,0,sizeof(readers_t));
}


void collapse_buffer(readers_t *files, reader_t *reader)
{
    int irec,jrec, has_snp=0, has_indel=0, has_any=0;
    for (irec=0; irec<reader->nbuffer; irec++)
    {
        bcf1_t *line = reader->buffer[irec];
        if ( line->pos != reader->buffer[0]->pos ) break;
        set_variant_types(line);
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
    irec = jrec = 0;
    while ( irec<reader->nbuffer && jrec<reader->nbuffer )
    {
        if ( reader->buffer[irec]->pos != -1 ) { irec++; continue; }
        if ( jrec<=irec ) jrec = irec+1;
        while ( jrec<reader->nbuffer && reader->buffer[jrec]->pos==-1 ) jrec++;
        if ( jrec<reader->nbuffer )
        {
            tmp = reader->buffer[irec]; reader->buffer[irec] = reader->buffer[jrec]; reader->buffer[jrec] = tmp;
        }
    }
    reader->nbuffer = irec;
}

void debug_buffer(readers_t *files)
{
    int i,j;
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        for (j=0; j<reader->nbuffer; j++)
        {
            bcf1_t *line = reader->buffer[j];
            printf("%s\t%s:%d\t%s %s\n", reader->fname,files->seqs[files->iseq],line->pos+1,line->d.allele[0],line->d.allele[1]);
        }
    }
    printf("\n\n");
}

int next_line(readers_t *files)
{
    int32_t min_pos = INT_MAX;
    int ret,i,j;
    kstring_t s = {0,0,0};

    while ( min_pos==INT_MAX )
    {
        // Need to open new chromosome?
        int eos = 0;
        for (i=0; i<files->nreaders; i++)
            if ( !files->readers[i].itr && !files->readers[i].nbuffer ) eos++;
        if ( eos==files->nreaders )
        {
            if ( ++files->iseq >= files->nseqs ) return 0;  // all chroms scanned
            for (i=0; i<files->nreaders; i++)
            {
                reader_t *reader = &files->readers[i];
                if ( reader->tbx )
                    reader->itr = tbx_itr_querys(reader->tbx,files->seqs[files->iseq]);
                else
                    reader->itr = bcf_itr_querys(reader->bcf,reader->header,files->seqs[files->iseq]);
            }
        }

        // Find the smallest coordinate
        for (i=0; i<files->nreaders; i++)
        {
            reader_t *reader = &files->readers[i];
            int buffer_full = ( reader->nbuffer && reader->buffer[reader->nbuffer-1]->pos != reader->buffer[0]->pos ) ? 1 : 0;
            if ( reader->itr && !buffer_full )
            {
                // Fill the buffer with records starting at the same position
                while (1)
                {
                    if ( reader->nbuffer >= reader->mbuffer ) 
                    {
                        reader->mbuffer += 5;
                        reader->buffer = (bcf1_t**) realloc(reader->buffer, sizeof(bcf1_t*)*reader->mbuffer);
                        for (j=5; j>0; j--)
                            reader->buffer[reader->mbuffer-j] = bcf_init1();
                    }
                    if ( reader->tbx )
                    {
                        ret = tbx_itr_next((BGZF*)reader->file->fp, reader->tbx, reader->itr, &s);
                        if ( ret<0 ) break;
                        vcf_parse1(&s, reader->header, reader->buffer[reader->nbuffer]);
                    }
                    else
                    {
                        ret = bcf_itr_next((BGZF*)reader->file->fp, reader->itr, reader->buffer[reader->nbuffer]);
                        if ( ret<0 ) break;
                    }
                    bcf_unpack(reader->buffer[reader->nbuffer], BCF_UN_STR|BCF_UN_FLT);
                    if ( reader->filter_id!=-1 && reader->buffer[reader->nbuffer]->d.n_flt && reader->filter_id!=reader->buffer[reader->nbuffer]->d.flt[0] ) continue;
                    reader->nbuffer++;
                    if ( reader->buffer[reader->nbuffer-1]->pos != reader->buffer[0]->pos ) break;
                }
                if ( ret<0 ) { tbx_itr_destroy(reader->itr); reader->itr = NULL; } // done for this chromosome
            }
            if ( reader->nbuffer )
            {
                if ( min_pos > reader->buffer[0]->pos ) min_pos = reader->buffer[0]->pos; 
            }
            if ( files->collapse && reader->nbuffer>1 && reader->buffer[0]->pos==reader->buffer[1]->pos )
                collapse_buffer(files, reader);
        }
    }

    //debug_buffer(files);

    // Set the current line
    ret = 0;
    bcf1_t *first = NULL;
    for (i=0; i<files->nreaders; i++)
    {
        reader_t *reader = &files->readers[i];
        if ( !reader->nbuffer || reader->buffer[0]->pos!=min_pos ) continue;

        // Match the records by REF and ALT
        int j, irec = -1;
        if ( first )
        {
            for (j=0; j<reader->nbuffer; j++)
            {
                bcf1_t *line = reader->buffer[j];
                if ( min_pos != line->pos ) break;  // done with this buffer

                if ( files->collapse&COLLAPSE_ANY ) { irec=j; break; }  // checking position only
                if ( files->collapse&COLLAPSE_SNPS && first->d.var_type&VCF_SNP && line->d.var_type&VCF_SNP ) { irec=j; break; }
                if ( files->collapse&COLLAPSE_INDELS && first->d.var_type&VCF_INDEL && line->d.var_type&VCF_INDEL ) { irec=j; break; }

                // thorough check: the REFs and some of the alleles have to be shared
                // (neglecting different representations of the same indel for now)
                if ( first->rlen != line->rlen ) continue;  // REFs do not match
                if ( strcmp(first->d.allele[0], line->d.allele[0]) ) continue; // REFs do not match
                int ial,jal;
                for (ial=1; ial<first->n_allele; ial++)
                {
                    for (jal=1; jal<line->n_allele; jal++)
                        if ( !strcmp(first->d.allele[ial], line->d.allele[jal]) ) { irec=j; break; }
                    if ( irec>=0 ) break;
                }
                if ( irec>=0 ) break;
            }
            if ( irec==-1 ) continue;
        }
        else 
        {
            first = reader->buffer[0];
            irec  = 0;
        }
        bcf1_t *tmp = reader->line;
        reader->line = reader->buffer[irec];
        for (j=irec+1; j<reader->nbuffer; j++)
            reader->buffer[j-1] = reader->buffer[j];
        reader->nbuffer--;
        reader->buffer[ reader->nbuffer ] = tmp;
        ret |= 1<<i;
    }
    if ( s.m ) free(s.s);

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

int init_samples(const char *fname, readers_t *files)
{
    int i;
    struct stat sbuf;
    files->samples = NULL;
    files->n_smpl  = 0;
    if ( !strcmp(fname,"-") )   // Intersection of all samples across all readers
    {
        int n;
        const char **smpl = (const char**) bcf_sample_names(files->readers[0].header,&n);
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
        free(smpl);
    }
    else if ( stat(fname, &sbuf) == 0 ) // read samples from file
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

int set_fmt_ptr(reader_t *reader, char *fmt)
{
    int i, gt_id = bcf_id2int(reader->header,BCF_DT_ID,fmt);
    if ( gt_id<0 ) return 0;
    bcf_unpack(reader->line, BCF_UN_FMT);
    reader->fmt_ptr = NULL;
    for (i=0; i<(int)reader->line->n_fmt; i++) 
        if ( reader->line->d.fmt[i].id==gt_id ) { reader->fmt_ptr = &reader->line->d.fmt[i]; break; }
    return reader->fmt_ptr ? 1 : 0;
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



