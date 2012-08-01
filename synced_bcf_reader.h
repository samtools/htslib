#ifndef SYNCED_BCF_READER_H
#define SYNCED_BCF_READER_H

#include "hts.h"
#include "vcf.h"
#include "tbx.h"

// How should be treated sites with the same position but different alleles
#define COLLAPSE_NONE   0
#define COLLAPSE_SNPS   1
#define COLLAPSE_INDELS 2
#define COLLAPSE_ANY    4

typedef struct { int32_t from, to; } pos_t;
typedef struct
{
	int *npos,nseqs,cpos,cseq;
	pos_t **pos;
	char **seq_names;
}
regions_t;

typedef struct
{
	htsFile *file;
	tbx_t *idx;
	bcf_hdr_t *header;
	hts_itr_t *itr;
	const char *fname;
	bcf1_t **buffer, *line;
	int nbuffer, mbuffer;
	int filter_id;
}
reader_t;

typedef struct
{
	reader_t *readers;
	int nreaders;
	const char **seqs, *region;
	int iseq,nseqs,mseqs;
	int collapse;
	int apply_filters;
}
readers_t;

int add_reader(const char *fname, readers_t *readers);
void destroy_readers(readers_t *readers);
int next_line(readers_t *readers);

int init_regions(const char *fname, regions_t *regions);
int reset_regions(regions_t *regions, const char *seq);
pos_t *is_in_regions(regions_t *regions, int32_t pos);
void destroy_regions(regions_t *regions);

#endif
