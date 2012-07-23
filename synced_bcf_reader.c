#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include "synced_bcf_reader.h"

int add_reader(const char *fname, readers_t *files)
{
	files->readers = (reader_t*) realloc(files->readers, sizeof(reader_t)*(files->nreaders+1));
	reader_t *reader = &files->readers[files->nreaders++];
	memset(reader,0,sizeof(reader_t));

	reader->idx = tbx_index_load(fname);
	if ( !reader->idx ) return 0;	// not indexed..?

	// Isn't there a better way then this: need to read VCF.gz header ("r") but also random access with tabix ("rb")
	reader->file = hts_open(fname, "r", NULL);
	if ( !reader->file ) return 0;
	reader->header = vcf_hdr_read(reader->file);
	hts_close(reader->file);

	reader->file = hts_open(fname, "rb", NULL);
	if ( !reader->file ) return 0;
	reader->line  = bcf_init1();
	reader->fname = fname;

	// Update list of chromosomes
	int n,i,j;
	const char **names = bcf_seqnames(reader->header, &n);
	for (i=0; i<n; i++)
	{
		for (j=0; j<files->nseqs; j++)
			if ( !strcmp(names[i],files->seqs[j]) ) break;
		if ( j<files->nseqs ) continue;		// already have this chr
		files->mseqs += 30;
		files->seqs = (const char**) realloc(files->seqs, sizeof(const char*)*files->mseqs);
		files->seqs[files->nseqs++] = names[i];
	}
	free(names);
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
		tbx_destroy(reader->idx);
		bcf_hdr_destroy(reader->header);
		hts_close(reader->file);
		if ( reader->itr ) tbx_itr_destroy(reader->itr);
		// destroy buffer, line
	}
	free(files->readers);
	free(files->seqs);
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
			if ( !files->readers[i].itr ) eos++;
		if ( eos==files->nreaders )
		{
			if ( ++files->iseq >= files->nseqs ) return 0;	// all chroms scanned

			char buf[100]; snprintf(buf,100,"%s:50000", files->seqs[files->iseq]);	// hack to temporarily work around a bug on Mac
			for (i=0; i<files->nreaders; i++)
			{
				reader_t *reader = &files->readers[i];
				reader->itr = tbx_itr_querys(reader->idx,buf);
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
				while ( (ret=tbx_itr_next((BGZF*)reader->file->fp, reader->idx,reader->itr, &s))>=0 )
				{
					if ( reader->nbuffer >= reader->mbuffer ) 
					{
						reader->mbuffer += 5;
						reader->buffer = (bcf1_t**) realloc(reader->buffer, sizeof(bcf1_t*)*reader->mbuffer);
						for (j=5; j>0; j--)
							reader->buffer[reader->mbuffer-j] = bcf_init1();
					}
					vcf_parse1(&s, reader->header, reader->buffer[reader->nbuffer]);
					bcf_unpack(reader->buffer[reader->nbuffer], BCF_UN_STR);
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
				if ( min_pos != line->pos ) break;	// done with this buffer

				if ( files->collapse&COLLAPSE_ANY ) { irec=j; break; }	// checking position only
				if ( files->collapse&COLLAPSE_SNPS && first->d.var_type&VCF_SNP && line->d.var_type&VCF_SNP ) { irec=j; break; }
				if ( files->collapse&COLLAPSE_INDELS && first->d.var_type&VCF_INDEL && line->d.var_type&VCF_INDEL ) { irec=j; break; }

				// thorough check: the REFs and some of the alleles have to be shared
				// (neglecting different representations of the same indel for now)
				if ( first->rlen != line->rlen ) continue;	// REFs do not match
				if ( strcmp(first->d.allele[0], line->d.allele[0]) ) continue; // REFs do not match
				int ial,jal;
				for (ial=1; ial<first->n_allele; ial++)
				{
					for (jal=1; jal<first->n_allele; jal++)
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


