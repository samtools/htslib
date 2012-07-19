#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include "synced_bcf_reader.h"

int add_reader(char *fname, readers_t *files)
{
	files->readers = (reader_t*) realloc(files->readers, sizeof(reader_t)*(files->nreaders+1));
	reader_t *reader = &files->readers[files->nreaders++];
	memset(reader,0,sizeof(reader_t));

	reader->idx = tbx_index_load(fname);
	if ( !reader->idx ) return 0;	// not indexed..?

	// There must be a better way then this: we need to be able to read VCF.gz header ("r") 
	//	but also acces with tabix ("rb")
	reader->file = hts_open(fname, "r", NULL);
	if ( !reader->file ) return 0;
	reader->header = vcf_hdr_read(reader->file);
	hts_close(reader->file);

	reader->file = hts_open(fname, "rb", NULL);
	if ( !reader->file ) return 0;

	// Update list of chromosomes
	int n,i,j;
	const char **names = tbx_seqnames(reader->idx, &n);
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

int next_line(readers_t *files)
{
	// Need to open new chromosome?
	int ret, i,j, eos = 0;
	for (i=0; i<files->nreaders; i++)
		if ( !files->readers[i].itr ) eos++;
	if ( eos==files->nreaders )
	{
		if ( ++files->iseq >= files->nseqs ) return 0;	// all chroms scanned

		char buf[100]; snprintf(buf,100,"%s:50000", files->seqs[files->iseq]);
		for (i=0; i<files->nreaders; i++)
		{
			reader_t *reader = &files->readers[i];
			reader->itr = tbx_itr_querys(reader->idx,buf);
	
	kstring_t s = {0,0,0};
	tbx_itr_next((BGZF*)reader->file->fp, reader->idx,reader->itr, &s);
	fprintf(stderr,"opening: %d -> %s, %lu .. %s\n", i,buf,(long unsigned int)reader->itr,s.s);
		}
	}

	// Find the smallest coordinate
	int32_t min_pos = INT_MAX;
	kstring_t s = {0,0,0};
	for (i=0; i<files->nreaders; i++)
	{
		reader_t *reader = &files->readers[i];
		if ( reader->itr )
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
				reader->nbuffer++;
				if ( reader->buffer[reader->nbuffer-1]->pos != reader->buffer[0]->pos ) break;
			}
			if ( ret<0 ) { tbx_itr_destroy(reader->itr); reader->itr = NULL; } // done for this chromosome
		}
		if ( reader->nbuffer )
		{
			if ( min_pos>reader->buffer[0]->pos ) min_pos = reader->buffer[0]->pos; 
			continue;
		}
	}
	if ( min_pos==INT_MAX ) return next_line(files); // No record available for this chromosome

	// Set the current line
	ret = 0;
	for (i=0; i<files->nreaders; i++)
	{
		reader_t *reader = &files->readers[i];
		if ( !reader->nbuffer || reader->buffer[0]->pos!=min_pos ) continue;
// here must match snps and indels together
		reader->line = reader->buffer[0];
		for (j=1; j<reader->nbuffer; j++)
			reader->buffer[j-1] = reader->buffer[j];
		reader->nbuffer--;
		ret |= 1<<i;
	}
	if ( s.m ) free(s.s);

	return ret;
}


