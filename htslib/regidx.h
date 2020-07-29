/// @file htslib/regidx.h
/// Region indexing.
/*
    Copyright (C) 2014-2019 Genome Research Ltd.

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Region indexing with an optional payload.

    Example of usage:

        // Init the parser and print regions. In this example the payload is a
        // pointer to a string. For the description of parse_custom and
        // free_custom functions, see regidx_parse_f and regidx_free_f below,
        // and for working example see test/test-regidx.c.
        regidx_t *idx = regidx_init(in_fname,parse_custom,free_custom,sizeof(char*),NULL);

        // Query overlap with chr:beg-end (beg,end are 1-based coordinates)
        regitr_t *itr = regitr_init(idx);
        if ( regidx_overlap(idx, chr,beg-1,end-1, itr) ) printf("There is an overlap!\n");

        while ( regitr_overlap(itr) )
        {
            printf("[%"PRIhts_pos",%"PRIhts_pos"] overlaps with [%"PRIhts_pos",%"PRIhts_pos"], payload=%s\n",
                   beg, end, itr->beg+1, itr->end+1, regitr_payload(itr,char*));
        }

        regidx_destroy(idx);
        regitr_destroy(itr);


    Another example, loop over all regions:

        regidx_t *idx = regidx_init(in_fname,NULL,NULL,0,NULL);
        regitr_t *itr = regitr_init(idx);

        while ( regitr_loop(itr) )
            printf("chr=%s  beg=%d  end=%d\n", itr->seq, itr->beg+1, itr->end+1);

        regidx_destroy(idx);
        regitr_destroy(itr);
*/

#ifndef HTSLIB_REGIDX_H
#define HTSLIB_REGIDX_H

#include "hts.h"

#ifdef __cplusplus
extern "C" {
#endif

// maximum regidx position (0-based).  Used to represent the end point of
// regions which do not explicitly set one.  regidx_push() also limits
// positions passed to it to be no bigger than this.

// Limit is set to ensure some internal values used by regidx keep within 32
// bits and to stop the index from getting too big.

#define REGIDX_MAX (1ULL << 35)

typedef struct regidx_t regidx_t;
typedef struct regitr_t
{
    hts_pos_t beg,end;
    void *payload;
    char *seq;
    void *itr;
}
regitr_t;

#define regitr_payload(itr,type_t) (*((type_t*)(itr)->payload))

// Old API for backwards compatibility
#define REGITR_START(itr) (itr).beg
#define REGITR_END(itr)   (itr).end
#define REGITR_PAYLOAD(itr,type_t) ((type_t*)(itr).payload)
#define REGITR_OVERLAP(itr,from,to) regidx_overlap((itr));

/*
 *  regidx_parse_f - Function to parse one input line, such as regidx_parse_bed
 *  or regidx_parse_tab below. The function is expected to set `chr_from` and
 *  `chr_to` to point to first and last character of chromosome name and set
 *  coordinates `beg` and `end` (0-based, inclusive). If regidx_init() was
 *  called with non-zero payload_size, the `payload` points to a memory
 *  location of the payload_size and `usr` is the data passed to regidx_init().
 *  Any memory allocated by the function will be freed by regidx_free_f called
 *  by regidx_destroy().
 *
 *  Return value: 0 on success, -1 to skip a record, -2 on fatal error.
 */
typedef int  (*regidx_parse_f)(const char *line, char **chr_beg, char **chr_end, hts_pos_t *beg, hts_pos_t *end, void *payload, void *usr);
typedef void (*regidx_free_f)(void *payload);

/*
 *  A note about the parsers:
 *      - leading spaces are ignored
 *      - lines starting with "#" are ignored
 */
HTSLIB_EXPORT
int regidx_parse_bed(const char*,char**,char**,hts_pos_t*,hts_pos_t*,void*,void*);   // CHROM or whitespace-sepatated CHROM,FROM,TO (0-based,right-open)
HTSLIB_EXPORT
int regidx_parse_tab(const char*,char**,char**,hts_pos_t*,hts_pos_t*,void*,void*);   // CHROM or whitespace-separated CHROM,POS (1-based, inclusive)
HTSLIB_EXPORT
int regidx_parse_reg(const char*,char**,char**,hts_pos_t*,hts_pos_t*,void*,void*);   // CHROM, CHROM:POS, CHROM:FROM-TO, CHROM:FROM- (1-based, inclusive)
HTSLIB_EXPORT
int regidx_parse_vcf(const char*,char**,char**,hts_pos_t*,hts_pos_t*,void*,void*);

/*
 *  regidx_init() - creates new index
 *  regidx_init_string() - creates new index, from a string rather than from a file
 *
 *  @param fname:  input file name or NULL if regions will be added one-by-one via regidx_insert()
 *  @param parsef: regidx_parse_bed, regidx_parse_tab or see description of regidx_parse_f. If NULL,
 *                 the format will be autodected, currently either regidx_parse_tab (the default) or
 *                 regidx_parse_bed (file must be named 'bed' or 'bed.gz') will be used. Note that
 *                 the exact autodetection algorithm will change.
 *  @param freef:  NULL or see description of regidx_parse_f
 *  @param payload_size: 0 with regidx_parse_bed, regidx_parse_tab or see regidx_parse_f
 *  @param usr:    optional user data passed to regidx_parse_f
 *
 *  Returns index on success or NULL on error.
 *
 *  The regidx_t index struct returned by a successful call should be freed
 *  via regidx_destroy() when it is no longer needed.
 */
HTSLIB_EXPORT
regidx_t *regidx_init(const char *fname, regidx_parse_f parsef, regidx_free_f freef, size_t payload_size, void *usr);
HTSLIB_EXPORT
regidx_t *regidx_init_string(const char *string, regidx_parse_f parsef, regidx_free_f freef, size_t payload_size, void *usr);

/*
 *  regidx_destroy() - free memory allocated by regidx_init
 */
HTSLIB_EXPORT
void regidx_destroy(regidx_t *idx);

/*
 *  regidx_overlap() - check overlap of the location chr:from-to with regions
 *  @param beg,end:     0-based start, end coordinate (inclusive)
 *  @param itr:         pointer to iterator, can be NULL if regidx_loop not needed
 *
 *  Returns 0 if there is no overlap or 1 if overlap is found. The overlapping
 *  regions can be iterated as shown in the example above.
 */
HTSLIB_EXPORT
int regidx_overlap(regidx_t *idx, const char *chr, hts_pos_t beg, hts_pos_t end, regitr_t *itr);

/*
 *  regidx_insert() - add a new region.
 *  regidx_insert_list() - add new regions from a list
 *  regidx_push() - low level insertion of a new region
 *
 *  Returns 0 on success or -1 on error.
 */
HTSLIB_EXPORT
int regidx_insert(regidx_t *idx, char *line);
HTSLIB_EXPORT
int regidx_insert_list(regidx_t *idx, char *line, char delim);
HTSLIB_EXPORT
int regidx_push(regidx_t *idx, char *chr_beg, char *chr_end, hts_pos_t beg, hts_pos_t end, void *payload);

/*
 *  regidx_seq_names() - return list of all sequence names
 */
HTSLIB_EXPORT
char **regidx_seq_names(regidx_t *idx, int *n);

/*
 *  regidx_seq_nregs() - number of regions
 *  regidx_nregs()  - total number of regions
 */
HTSLIB_EXPORT
int regidx_seq_nregs(regidx_t *idx, const char *seq);

HTSLIB_EXPORT
int regidx_nregs(regidx_t *idx);

/*
 *  regitr_init() - initialize an iterator. The idx parameter is required only
 *                  with regitr_loop. If only regitr_overlap is called, NULL
 *                  can be given.
 *
 *                  The regitr_t struct returned by a successful regitr_init()
 *                  call should be freed via regitr_destroy() when it is no
 *                  longer needed.
 *
 *  regitr_reset() - initialize an iterator for a repeated regitr_loop cycle.
 *                  Not required with regitr_overlap.
 */
HTSLIB_EXPORT
regitr_t *regitr_init(regidx_t *idx);
HTSLIB_EXPORT
void regitr_destroy(regitr_t *itr);
HTSLIB_EXPORT
void regitr_reset(regidx_t *idx, regitr_t *itr);

/*
 *  regitr_overlap() - next overlapping region
 *  Returns 0 when done or 1 when itr is set to next region
 */
HTSLIB_EXPORT
int regitr_overlap(regitr_t *itr);

/*
 *  regitr_loop() - loop over all regions
 *  Returns 0 when done or 1 when itr is set to next region
 */
HTSLIB_EXPORT
int regitr_loop(regitr_t *itr);

/*
 *  regitr_copy() - create a copy of an iterator for a repeated iteration with regitr_loop
 */
HTSLIB_EXPORT
void regitr_copy(regitr_t *dst, regitr_t *src);

#ifdef __cplusplus
}
#endif

#endif
