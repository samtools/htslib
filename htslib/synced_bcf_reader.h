/*
	The synced_bcf_reader allows to keep multiple VCFs open and stream them
	using the next_line iterator in a seamless matter without worrying about
	chromosomes and synchronizing the sites. This is used by vcfcheck to
	compare multiple VCFs simultaneously and will be used also for merging,
	creating intersections, etc.

	The *_regions class of functions provide a convenient way for checking if a
	coordinate is inside one of the regions.
*/

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
	pos_t **pos, tpos;  // **pos and npos will be deprecated with *_regions functions
	char **seq_names;
    tbx_t *tbx;
    htsFile *file;
    hts_itr_t *itr;
    kstring_t line;
}
regions_t;

typedef struct
{
	htsFile *file;
    tbx_t *tbx;
    hts_idx_t *bcf;
	bcf_hdr_t *header;
	hts_itr_t *itr;
	const char *fname;
	bcf1_t **buffer;
	int nbuffer, mbuffer;
	int filter_id;
	int *samples, n_smpl;	// list of columns in the order consistent with readers_t.samples
//	bcf_fmt_t *fmt_ptr;	    // set by set_fmt_ptr
}
reader_t;

typedef struct
{
	// Parameters controlling the logic
	int collapse;  // How should the duplicate sites be treated
	int apply_filters;	// If set to non-zero, sites with FILTER value other
						// than "." or "PASS" will be skipped. Active only at
						// the time of initialization, that is during the
						// add_reader() calls.
	// Auxiliary data
	reader_t *readers;
	int nreaders;
	const char **seqs, *region;
	int iseq,nseqs,mseqs;
	char **samples;	// List of samples 
    regions_t *targets;
    kstring_t tmps;
	int n_smpl;
}
readers_t;

/** Init readers_t struct */
readers_t *bcf_sr_init();

/** Destroy  readers_t struct */
void bcf_sr_destroy(readers_t *readers);

/**
 *  bcf_sr_add_reader() - open new reader
 *  @readers: holder of the open readers
 *  @fname:   the VCF file
 *
 *  Returns 1 if the call succeeded, or 0 on error.
 *
 *  See also the readers_t data structure for parameters controlling
 *  the reader's logic.
 */
int bcf_sr_add_reader(readers_t *readers, const char *fname);

/** 
 * bcf_sr_next_line() - the iterator
 * @readers:    holder of the open readers
 *
 * Returns 0 when all lines from all files have been read or a bit mask
 * indicating which of the readers have a reader_t.line set at this position. 
 */
int bcf_sr_next_line(readers_t *readers);

/**
 * bcf_sr_set_samples() - sets active samples
 * @readers: holder of the open readers
 * @samples: this can be one of: file name with one sample per line;
 *           or column-separated list of samples; or '-' for a list of 
 *           samples shared by all files
 *
 * Returns 1 if the call succeeded, or 0 on error.
 */
int bcf_sr_set_samples(readers_t *readers, const char *samples);

/**
 *  bcf_sr_set_targets() - init positions 
 *  @readers: holder of the open readers
 *  @fname:   tabix indexed tab delimited file <chr,pos> or <chr,from,to>,
 *            coordinates 1-based and inclusive
 *
 *  Returns 1 if the call succeeded, or 0 on error.
 *
 *  Sets target regions; positions not listed will be skipped by
 *  bcf_sr_next_line(). Note that the files are streamed anyway, the index is
 *  required only for reading whole chromosomes blocks.
 */
int bcf_sr_set_targets(readers_t *readers, const char *fname);



/**
 * init_regions() - initialize regions_t structure
 * @fname:   bgzip compressed tab delimited file with chr,from,to.
 *           Coordinates are one-based and inclusive
 * @regions: regions_t structure
 *
 * Returns 1 if the call succeeded, or 0 on error.
 *
 * In the current implementation the regions list is assumed to be
 * relatively small as the regions are kept in memory. (It is
 * OK for human exons for example.) If this assumption is no longer
 * valid for some applications, it is straightforward to require
 * tabix indexed file and read the regions from the file on fly.
 */
int init_regions(const char *fname, regions_t *regions);

/**
 * reset_regions() - position regions to next chromosome
 * @regions: regions_t structure
 * @seq:     chromosome name or region
 *
 * Returns 1 if the call succeeded, or 0 when $seq not found.
 */
int reset_regions(regions_t *regions, const char *seq);

/**
 * is_in_regions() - looks up the region containing $pos
 * @regions: regions_t structure
 * @pos:     query coordinate
 *
 * Returns pointer to the region containing the position or
 * NULL if no such position was found.
 */
pos_t *is_in_regions(regions_t *regions, int32_t pos);

/**
 * destroy_regions() - free memory occupied by regions
 * @regions: regions_t structure
 */
void destroy_regions(regions_t *regions);

#endif
