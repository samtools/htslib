#ifndef SYNCED_BCF_READER_H
#define SYNCED_BCF_READER_H

#include "hts.h"
#include "vcf.h"
#include "tbx.h"

typedef struct
{
	htsFile *file;
	tbx_t *idx;
	bcf_hdr_t *header;
	hts_itr_t *itr;
	const char *fname;
	bcf1_t **buffer, *line;
	int nbuffer, mbuffer;
}
reader_t;

typedef struct
{
	reader_t *readers;
	int nreaders;
	const char **seqs;
	int iseq,nseqs,mseqs;
}
readers_t;

int add_reader(const char *fname, readers_t *readers);
void destroy_readers(readers_t *readers);
int next_line(readers_t *readers);

#endif
