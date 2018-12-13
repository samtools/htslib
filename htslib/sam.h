/// @file htslib/sam.h
/// High-level SAM/BAM/CRAM sequence file operations.
/*
    Copyright (C) 2008, 2009, 2013-2017 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#ifndef HTSLIB_SAM_H
#define HTSLIB_SAM_H

#include <stdint.h>
#include "hts.h"

#ifdef __cplusplus
extern "C" {
#endif

/// Highest SAM format version supported by this library
#define SAM_FORMAT_VERSION "1.6"

/**********************
 *** SAM/BAM header ***
 **********************/

typedef struct sam_hdr bam_hrecs_t;

/*! @typedef
 @abstract Structure for the alignment header.
 @field n_targets   number of reference sequences
 @field l_text      length of the plain text in the header
 @field target_len  lengths of the reference sequences
 @field target_name names of the reference sequences
 @field text        plain text
 @field sdict       header dictionary
 */

typedef struct bam_hdr_t {
    int32_t n_targets, ignore_sam_err;
    uint32_t l_text;
    uint32_t *target_len;
    int8_t *cigar_tab;
    char **target_name;
    char *text;
    void *sdict;
    bam_hrecs_t *hrecs;
    uint32_t ref_count;
} bam_hdr_t;

/****************************
 *** CIGAR related macros ***
 ****************************/

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
// Note that BAM_CIGAR_STR is padded to length 16 bytes below so that
// the array look-up will not fall off the end.  '?' is chosen as the
// padding character so it's easy to spot if one is emitted, and will
// result in a parsing failure (in sam_parse1(), at least) if read.
#define bam_cigar_opchr(c) (BAM_CIGAR_STR "??????" [bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

/*************************
 *** Alignment records ***
 *************************/

/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 */
typedef struct {
    int32_t tid;
    int32_t pos;
    uint16_t bin;
    uint8_t qual;
    uint8_t l_qname;
    uint16_t flag;
    uint8_t unused1;
    uint8_t l_extranul;
    uint32_t n_cigar;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

 @discussion Notes:

 1. The data blob should be accessed using bam_get_qname, bam_get_cigar,
    bam_get_seq, bam_get_qual and bam_get_aux macros.  These returns pointers
    to the start of each type of data.
 2. qname is terminated by one to four NULs, so that the following
    cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
    while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
 3. Cigar data is encoded 4 bytes per CIGAR operation.
    See the bam_cigar_* macros for manipulation.
 4. seq is nibble-encoded according to bam_nt16_table.
    See the bam_seqi macro for retrieving individual bases.
 5. Per base qualilties are stored in the Phred scale with no +33 offset.
    Ie as per the BAM specification and not the SAM ASCII printable method.
 */
typedef struct {
    bam1_core_t core;
    int l_data;
    uint32_t m_data;
    uint8_t *data;
#ifndef BAM_NO_ID
    uint64_t id;
#endif
} bam1_t;

/*! @function
 @abstract  Get whether the query is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
#define bam_is_rev(b) (((b)->core.flag&BAM_FREVERSE) != 0)
/*! @function
 @abstract  Get whether the query's mate is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query's mate on the reverse strand
 */
#define bam_is_mrev(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
/*! @function
 @abstract  Get the name of the query
 @param  b  pointer to an alignment
 @return    pointer to the name string, null terminated
 */
#define bam_get_qname(b) ((char*)(b)->data)
/*! @function
 @abstract  Get the CIGAR array
 @param  b  pointer to an alignment
 @return    pointer to the CIGAR array

 @discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
 */
#define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
/*! @function
 @abstract  Get query sequence
 @param  b  pointer to an alignment
 @return    pointer to sequence

 @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
 8 for T and 15 for N. Two bases are packed in one byte with the base
 at the higher 4 bits having smaller coordinate on the read. It is
 recommended to use bam_seqi() macro to get the base.
 */
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
/*! @function
 @abstract  Get query quality
 @param  b  pointer to an alignment
 @return    pointer to quality string
 */
#define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
/*! @function
 @abstract  Get auxiliary data
 @param  b  pointer to an alignment
 @return    pointer to the concatenated auxiliary data
 */
#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
/*! @function
 @abstract  Get length of auxiliary data
 @param  b  pointer to an alignment
 @return    length of the concatenated auxiliary data
 */
#define bam_get_l_aux(b) ((b)->l_data - ((b)->core.n_cigar<<2) - (b)->core.l_qname - (b)->core.l_qseq - (((b)->core.l_qseq + 1)>>1))
/*! @function
 @abstract  Get a base on read
 @param  s  Query sequence returned by bam_get_seq()
 @param  i  The i-th position, 0-based
 @return    4-bit integer representing the base.
 */
#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

/**************************
 *** Exported functions ***
 **************************/

/***************
 *** BAM I/O ***
 ***************/

/* Header */

bam_hdr_t *bam_hdr_init(void);
bam_hdr_t *bam_hdr_read(BGZF *fp);
int bam_hdr_write(BGZF *fp, bam_hdr_t *h) HTS_RESULT_USED;
void bam_hdr_destroy(bam_hdr_t *h);
int bam_name2id(bam_hdr_t *h, const char *ref);
const char *sam_parse_region(bam_hdr_t *h, const char *s, int *tid, int64_t *beg, int64_t *end, int flags);
bam_hdr_t* bam_hdr_dup(const bam_hdr_t *h0);

typedef htsFile samFile;
bam_hdr_t *sam_hdr_parse(int l_text, const char *text);
bam_hdr_t *sam_hdr_read(samFile *fp);
int sam_hdr_write(samFile *fp, bam_hdr_t *h) HTS_RESULT_USED;

/*!
 * Returns the current length of the bam_hdr_t in text form.
 */
int sam_hdr_length2(bam_hdr_t *bh);

/*!
 * Returns the string form of the bam_hdr_t.
 */
const char *sam_hdr_str2(bam_hdr_t *bh);

/*! Reconstructs the kstring from the header hash table.
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_rebuild2(bam_hdr_t *bh);

/* ==== Line level methods ==== */

/*! Appends a formatted line to an existing SAM header.
 *
 * Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 * optional new-line. If it contains more than 1 line then multiple lines
 * will be added in order.
 *
 * Input text is of maximum length len or as terminated earlier by a NUL.
 * Len may be 0 if unknown, in which case lines must be NUL-terminated.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add_lines2(bam_hdr_t *bh, const char *lines, int len);

/*! Adds a single line to a SAM header.
 *
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. sam_hdr_add_line(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add_line2(bam_hdr_t *bh, const char *type, ...);

/*!
 *
 * Returns a complete line of formatted text
 * for a specific head type/ID combination. If ID is NULL then it returns
 * the first line of the specified type.
 *
 * The returned string is malloced and should be freed by the calling
 * function with free().
 *
 * @return
 * Returns NULL if no type/ID is found.
 */
char *sam_hdr_find_line2(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value);

/*!
 *
 * Remove a line from the header by specifying a tag:value that uniquely
 * identifies a line, i.e. the @SQ line containing "SN:ref1".
 * @SQ line is uniquely identified by SN tag.
 * @RG line is uniquely identified by ID tag.
 * @PG line is uniquely identified by ID tag.
 *
 * @return
 * Returns 0 on success and -1 on error.
 */
int sam_hdr_remove_line_key2(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value);

/*!
 *
 * Remove a line from the header by specifying the position in the type
 * group, i.e. 3rd @SQ line.
 */
int sam_hdr_remove_line_pos2(bam_hdr_t *bh, const char *type, int position);

/*! Adds or updates tag key,value pairs in a header line.
 *
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line.
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_update_line(bam_hdr_t *bh, const char *type,
        const char *ID_key, const char *ID_value, ...);

int sam_hdr_update_HD(bam_hdr_t *bh, ...);

/*!
 *
 * Remove all lines of type <type> from the header, except the one specified by tag:value, i.e. the @SQ line containing "SN:ref1".
 * @param type     - the selected type, e.g. SQ, RG.
 * @param ID_key   - the key that uniquely identifies the respective type, e.g. SN for SQ, ID for RG.
 * @param ID_value - the value of the key, e.g. SN:ref1 for SQ, ID:rg1 for RG.
 *
 * @return
 * Returns 0 on success and -1 on error.
 */
int sam_hdr_keep_line(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value);

/* ==== Key:val level methods ==== */

/*! Looks for a specific key in a single sam header line and returns the
 * associated value.
 *
 * @return
 * Returns a char pointer on success;
 *         NULL on failure
 */
const char *sam_hdr_find_tag2(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value, const char *key);

/* ! Remove the key from the line identified by type, ID_key and ID_value
 *
 */
int sam_hdr_remove_tag2(bam_hdr_t *bh, const char *type, const char *ID_key, const char *ID_value, const char *key);

/*! Looks up a reference sequence by name and returns the numerical ID.
 * @return
 * Returns -1 if unknown reference.
 */
int sam_hdr_name2ref2(bam_hdr_t *bh, const char *ref);

/*! Fixes any PP links in @PG headers.
 *
 * If the entries are in order then this doesn't need doing, but incase
 * our header is out of order this goes through the hrecs->pg[] array
 * setting the prev_id field.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure (indicating broken PG/PP records)
 */
int sam_hdr_link_pg2(bam_hdr_t *bh);

/*
 * Returns a unique ID from a base name.
 *
 * The value returned is valid until the next call to
 * this function.
 */
const char *sam_hdr_pg_id(bam_hdr_t *bh, const char *name);

/*! Add an @PG line.
 *
 * If we wish complete control over this use sam_hdr_add() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple @PG records if we have multiple PG chains.
 *
 * Call it as per sam_hdr_add() with a series of key,value pairs ending
 * in NULL.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add_pg(bam_hdr_t *bh, const char *name, ...);

/* Alignment */

bam1_t *bam_init1(void);
void bam_destroy1(bam1_t *b);
int bam_read1(BGZF *fp, bam1_t *b) HTS_RESULT_USED;
int bam_write1(BGZF *fp, const bam1_t *b) HTS_RESULT_USED;
bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc);
bam1_t *bam_dup1(const bam1_t *bsrc);

int bam_cigar2qlen(int n_cigar, const uint32_t *cigar);
int bam_cigar2rlen(int n_cigar, const uint32_t *cigar);

/*!
      @abstract Calculate the rightmost base position of an alignment on the
      reference genome.

      @param  b  pointer to an alignment
      @return    the coordinate of the first base after the alignment, 0-based

      @discussion For a mapped read, this is just b->core.pos + bam_cigar2rlen.
      For an unmapped read (either according to its flags or if it has no cigar
      string), we return b->core.pos + 1 by convention.
 */
int32_t bam_endpos(const bam1_t *b);

int   bam_str2flag(const char *str);    /** returns negative value on error */
char *bam_flag2str(int flag);   /** The string must be freed by the user */

/*************************
 *** BAM/CRAM indexing ***
 *************************/

// These BAM iterator functions work only on BAM files.  To work with either
// BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
#define bam_itr_destroy(iter) hts_itr_destroy(iter)
#define bam_itr_queryi(idx, tid, beg, end) sam_itr_queryi(idx, tid, beg, end)
#define bam_itr_querys(idx, hdr, region) sam_itr_querys(idx, hdr, region)
#define bam_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), 0)

// Load/build .csi or .bai BAM index file.  Does not work with CRAM.
// It is recommended to use the sam_index_* functions below instead.
#define bam_index_load(fn) hts_idx_load((fn), HTS_FMT_BAI)
#define bam_index_build(fn, min_shift) (sam_index_build((fn), (min_shift)))

/// Initialise fp->idx for the current format type for SAM, BAM and CRAM types .
/** @param fp        File handle for the data file being written.
    @param h         Bam header structured (needed for BAI and CSI).
    @param min_shift 0 for BAI, or larger for CSI (CSI defaults to 14).
    @param fnidx     Filename to write index to.  This pointer must remain valid
                     until after sam_idx_save is called.
    @return          0 on success, <0 on failure.

    @note This must be called after the header has been written, but before
          any other data.
*/
int sam_idx_init(htsFile *fp, bam_hdr_t *h, int min_shift, const char *fnidx);

/// Writes the index initialised with sam_idx_init to disk.
/** @param fp        File handle for the data file being written.
    @return          0 on success, <0 on filaure.
*/
int sam_idx_save(htsFile *fp) HTS_RESULT_USED;

/// Load a BAM (.csi or .bai) or CRAM (.crai) index file
/** @param fp  File handle of the data file whose index is being opened
    @param fn  BAM/CRAM/etc filename to search alongside for the index file
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *sam_index_load(htsFile *fp, const char *fn);

/// Load a specific BAM (.csi or .bai) or CRAM (.crai) index file
/** @param fp     File handle of the data file whose index is being opened
    @param fn     BAM/CRAM/etc data file filename
    @param fnidx  Index filename, or NULL to search alongside @a fn
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx);

/// Generate and save an index file
/** @param fn        Input BAM/etc filename, to which .csi/etc will be added
    @param min_shift Positive to generate CSI, or 0 to generate BAI
    @return  0 if successful, or negative if an error occurred (usually -1; or
             -2: opening fn failed; -3: format not indexable; -4:
             failed to create and/or save the index)
*/
int sam_index_build(const char *fn, int min_shift) HTS_RESULT_USED;

/// Generate and save an index to a specific file
/** @param fn        Input BAM/CRAM/etc filename
    @param fnidx     Output filename, or NULL to add .bai/.csi/etc to @a fn
    @param min_shift Positive to generate CSI, or 0 to generate BAI
    @return  0 if successful, or negative if an error occurred (see
             sam_index_build for error codes)
*/
int sam_index_build2(const char *fn, const char *fnidx, int min_shift) HTS_RESULT_USED;

/// Generate and save an index to a specific file
/** @param fn        Input BAM/CRAM/etc filename
    @param fnidx     Output filename, or NULL to add .bai/.csi/etc to @a fn
    @param min_shift Positive to generate CSI, or 0 to generate BAI
    @param nthreads  Number of threads to use when building the index
    @return  0 if successful, or negative if an error occurred (see
             sam_index_build for error codes)
*/
int sam_index_build3(const char *fn, const char *fnidx, int min_shift, int nthreads) HTS_RESULT_USED;

/// Free a SAM iterator
/// @param iter     Iterator to free
#define sam_itr_destroy(iter) hts_itr_destroy(iter)

/// Create a BAM/CRAM iterator
/** @param idx     Index
    @param tid     Target id
    @param beg     Start position in target
    @param end     End position in target
    @return An iterator on success; NULL on failure

The following special values (defined in htslib/hts.h)can be used for @p tid.
When using one of these values, @p beg and @p end are ignored.

  HTS_IDX_NOCOOR iterates over unmapped reads sorted at the end of the file
  HTS_IDX_START  iterates over the entire file
  HTS_IDX_REST   iterates from the current position to the end of the file
  HTS_IDX_NONE   always returns "no more alignment records"

When using HTS_IDX_REST or HTS_IDX_NONE, NULL can be passed in to @p idx.
 */
hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);

/// Create a SAM/BAM/CRAM iterator
/** @param idx     Index
    @param hdr     Header
    @param region  Region specification
    @return An iterator on success; NULL on failure

Regions are parsed by hts_parse_reg(), and take one of the following forms:

region          | Outputs
--------------- | -------------
REF             | All reads with RNAME REF
REF:            | All reads with RNAME REF
REF:START       | Reads with RNAME REF overlapping START to end of REF
REF:-END        | Reads with RNAME REF overlapping start of REF to END
REF:START-END   | Reads with RNAME REF overlapping START to END
.               | All reads from the start of the file
*               | Unmapped reads at the end of the file (RNAME '*' in SAM)

The form `REF:` should be used when the reference name itself contains a colon.

Note that SAM files must be bgzf-compressed for iterators to work.
 */
hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region);

/// Create a multi-region iterator
/** @param idx       Index
    @param hdr       Header
    @param reglist   Array of regions to iterate over
    @param regcount  Number of items in reglist

Each @p reglist entry should have the reference name in the `reg` field, an
array of regions for that reference in `intervals` and the number of items
in `intervals` should be stored in `count`.  No other fields need to be filled
in.

The iterator will return all reads overlapping the given regions.  If a read
overlaps more than one region, it will only be returned once.
 */
hts_itr_t *sam_itr_regions(const hts_idx_t *idx, bam_hdr_t *hdr, hts_reglist_t *reglist, unsigned int regcount);

/// Create a multi-region iterator
/** @param idx       Index
    @param hdr       Header
    @param regarray  Array of ref:interval region specifiers
    @param regcount  Number of items in regarray

Each @p regarray entry is parsed by hts_parse_reg(), and takes one of the
following forms:

region          | Outputs
--------------- | -------------
REF             | All reads with RNAME REF
REF:            | All reads with RNAME REF
REF:START       | Reads with RNAME REF overlapping START to end of REF
REF:-END        | Reads with RNAME REF overlapping start of REF to END
REF:START-END   | Reads with RNAME REF overlapping START to END
.               | All reads from the start of the file
*               | Unmapped reads at the end of the file (RNAME '*' in SAM)

The form `REF:` should be used when the reference name itself contains a colon.

The iterator will return all reads overlapping the given regions.  If a read
overlaps more than one region, it will only be returned once.
 */
hts_itr_t *sam_itr_regarray(const hts_idx_t *idx, bam_hdr_t *hdr, char **regarray, unsigned int regcount);

/// Get the next read from a SAM/BAM/CRAM iterator
/** @param htsfp       Htsfile pointer for the input file
    @param itr         Iterator
    @param r           Pointer to a bam1_t struct
    @return >= 0 on success; -1 when there is no more data; < -1 on error
 */
static inline int sam_itr_next(htsFile *htsfp, hts_itr_t *itr, bam1_t *r) {
    if (!htsfp->is_bgzf && !htsfp->is_cram) {
        hts_log_error("%s not BGZF compressed", htsfp->fn ? htsfp->fn : "File");
        return -2;
    }
    if (!itr) {
        hts_log_error("Null iterator");
        return -2;
    }

    if (itr->multi)
        return hts_itr_multi_next(htsfp, itr, r);
    else
        return hts_itr_next(htsfp->is_bgzf ? htsfp->fp.bgzf : NULL, itr, r, htsfp);
}

/// Get the next read from a BAM/CRAM multi-iterator
/** @param htsfp       Htsfile pointer for the input file
    @param itr         Iterator
    @param r           Pointer to a bam1_t struct
    @return >= 0 on success; -1 when there is no more data; < -1 on error
 */
#define sam_itr_multi_next(htsfp, itr, r) sam_itr_next(htsfp, itr, r)

    /***************
     *** SAM I/O ***
     ***************/

    #define sam_open(fn, mode) (hts_open((fn), (mode)))
    #define sam_open_format(fn, mode, fmt) (hts_open_format((fn), (mode), (fmt)))
    #define sam_close(fp) hts_close(fp)

    int sam_open_mode(char *mode, const char *fn, const char *format);

    // A version of sam_open_mode that can handle ,key=value options.
    // The format string is allocated and returned, to be freed by the caller.
    // Prefix should be "r" or "w",
    char *sam_open_mode_opts(const char *fn,
                             const char *mode,
                             const char *format);

    int sam_hdr_change_HD(bam_hdr_t *h, const char *key, const char *val);

    int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b) HTS_RESULT_USED;
    int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str) HTS_RESULT_USED;

/// sam_read1 - Read a record from a file
/** @param fp   Pointer to the source file
 *  @param h    Pointer to the header previously read (fully or partially)
 *  @param b    Pointer to the record placeholder
 *  @return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error
 */
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b) HTS_RESULT_USED;
/// sam_write1 - Write a record to a file
/** @param fp    Pointer to the destination file
 *  @param h     Pointer to the header structure previously read
 *  @param b     Pointer to the record to be written
 *  @return >= 0 on successfully writing the record, -1 on error
 */
    int sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b) HTS_RESULT_USED;

    /*************************************
     *** Manipulating auxiliary fields ***
     *************************************/

/// Return a pointer to an aux record
/** @param b   Pointer to the bam record
    @param tag Desired aux tag
    @return Pointer to the tag data, or NULL if tag is not present or on error
    If the tag is not present, this function returns NULL and sets errno to
    ENOENT.  If the bam record's aux data is corrupt (either a tag has an
    invalid type, or the last record is incomplete) then errno is set to
    EINVAL and NULL is returned.
 */
uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);

/// Get an integer aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get()
    @return The value, or 0 if the tag was not an integer type
    If the tag is not an integer type, errno is set to EINVAL.  This function
    will not return the value of floating-point tags.
*/
int64_t bam_aux2i(const uint8_t *s);

/// Get an integer aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get()
    @return The value, or 0 if the tag was not an integer type
    If the tag is not an numeric type, errno is set to EINVAL.  The value of
    integer flags will be returned cast to a double.
*/
double bam_aux2f(const uint8_t *s);

/// Get a character aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return The value, or 0 if the tag was not a character ('A') type
    If the tag is not a character type, errno is set to EINVAL.
*/
char bam_aux2A(const uint8_t *s);

/// Get a string aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return Pointer to the string, or NULL if the tag was not a string type
    If the tag is not a string type ('Z' or 'H'), errno is set to EINVAL.
*/
char *bam_aux2Z(const uint8_t *s);

/// Get the length of an array-type ('B') tag
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return The length of the array, or 0 if the tag is not an array type.
    If the tag is not an array type, errno is set to EINVAL.
 */
uint32_t bam_auxB_len(const uint8_t *s);

/// Get an integer value from an array-type tag
/** @param s   Pointer to the tag data, as returned by bam_aux_get().
    @param idx 0-based Index into the array
    @return The idx'th value, or 0 on error.
    If the array is not an integer type, errno is set to EINVAL.  If idx
    is greater than or equal to  the value returned by bam_auxB_len(s),
    errno is set to ERANGE.  In both cases, 0 will be returned.
 */
int64_t bam_auxB2i(const uint8_t *s, uint32_t idx);

/// Get a floating-point value from an array-type tag
/** @param s   Pointer to the tag data, as returned by bam_aux_get().
    @param idx 0-based Index into the array
    @return The idx'th value, or 0.0 on error.
    If the array is not a numeric type, errno is set to EINVAL.  This can
    only actually happen if the input record has an invalid type field.  If
    idx is greater than or equal to  the value returned by bam_auxB_len(s),
    errno is set to ERANGE.  In both cases, 0.0 will be returned.
 */
double bam_auxB2f(const uint8_t *s, uint32_t idx);

/// Append tag data to a bam record
/* @param b    The bam record to append to.
   @param tag  Tag identifier
   @param type Tag data type
   @param len  Length of the data in bytes
   @param data The data to append
   @return 0 on success; -1 on failure.
If there is not enough space to store the additional tag, errno is set to
ENOMEM.  If the type is invalid, errno may be set to EINVAL.  errno is
also set to EINVAL if the bam record's aux data is corrupt.
*/
int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data);

/// Delete tag data from a bam record
/* @param b The bam record to update
   @param s Pointer to the tag to delete, as returned by bam_aux_get().
   @return 0 on success; -1 on failure
   If the bam record's aux data is corrupt, errno is set to EINVAL and this
   function returns -1;
*/
int bam_aux_del(bam1_t *b, uint8_t *s);

/// Update or add a string-type tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param len  The length of the new string
   @param data The new string
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of type 'Z'.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data);

/// Update or add an integer tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param val  The new value
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of an integer type (c, C, s, S, i or I).

   EOVERFLOW (or ERANGE on systems that do not have EOVERFLOW): val is
   outside the range that can be stored in an integer bam tag (-2147483647
   to 4294967295).

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_int(bam1_t *b, const char tag[2], int64_t val);

/// Update or add a floating-point tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param val  The new value
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of a float type.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_float(bam1_t *b, const char tag[2], float val);

/// Update or add an array tag
/* @param b     The bam record to update
   @param tag   Tag identifier
   @param type  Data type (one of c, C, s, S, i, I or f)
   @param items Number of items
   @param data  Pointer to data
   @return 0 on success, -1 on failure
   The type parameter indicates the how the data is interpreted:

   Letter code | Data type | Item Size (bytes)
   ----------- | --------- | -----------------
   c           | int8_t    | 1
   C           | uint8_t   | 1
   s           | int16_t   | 2
   S           | uint16_t  | 2
   i           | int32_t   | 4
   I           | uint32_t  | 4
   f           | float     | 4

   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.  The bam record
   will grow or shrink in order to accomodate the new data.

   The data parameter must not point to any data in the bam record itself or
   undefined behaviour may result.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt, an existing tag with the
   given ID is not of an array type or the type parameter is not one of
   the values listed above.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_array(bam1_t *b, const char tag[2],
                         uint8_t type, uint32_t items, void *data);

/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

/*! @typedef
 @abstract Generic pileup 'client data'.

 @discussion The pileup iterator allows setting a constructor and
 destructor function, which will be called every time a sequence is
 fetched and discarded.  This permits caching of per-sequence data in
 a tidy manner during the pileup process.  This union is the cached
 data to be manipulated by the "client" (the caller of pileup).
*/
typedef union {
    void *p;
    int64_t i;
    double f;
} bam_pileup_cd;

/*! @typedef
 @abstract Structure for one alignment covering the pileup position.
 @field  b          pointer to the alignment
 @field  qpos       position of the read base at the pileup site, 0-based
 @field  indel      indel length; 0 for no indel, positive for ins and negative for del
 @field  level      the level of the read in the "viewer" mode
 @field  is_del     1 iff the base on the padded read is a deletion
 @field  is_head    1 iff this is the first base in the query sequence
 @field  is_tail    1 iff this is the last base in the query sequence
 @field  is_refskip 1 iff the base on the padded read is part of CIGAR N op
 @field  aux        (used by bcf_call_gap_prep())
 @field  cigar_ind  index of the CIGAR operator that has just been processed

 @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
 difference between the two functions is that the former does not
 set bam_pileup1_t::level, while the later does. Level helps the
 implementation of alignment viewers, but calculating this has some
 overhead.
 */
typedef struct {
    bam1_t *b;
    int32_t qpos;
    int indel, level;
    uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, /* reserved */ :1, aux:27;
    bam_pileup_cd cd; // generic per-struct data, owned by caller.
    int cigar_ind;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

struct __bam_plp_t;
typedef struct __bam_plp_t *bam_plp_t;

struct __bam_mplp_t;
typedef struct __bam_mplp_t *bam_mplp_t;

    /**
     *  bam_plp_init() - sets an iterator over multiple
     *  @func:      see mplp_func in bam_plcmd.c in samtools for an example. Expected return
     *              status: 0 on success, -1 on end, < -1 on non-recoverable errors
     *  @data:      user data to pass to @func
     */
    bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
    void bam_plp_destroy(bam_plp_t iter);
    int bam_plp_push(bam_plp_t iter, const bam1_t *b);
    const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
    void bam_plp_reset(bam_plp_t iter);

    /**
     *  bam_plp_constructor() - sets a callback to initialise any per-pileup1_t fields.
     *  @plp:       The bam_plp_t initialised using bam_plp_init.
     *  @func:      The callback function itself.  When called, it is given the
     *              data argument (specified in bam_plp_init), the bam structure and
     *              a pointer to a locally allocated bam_pileup_cd union.  This union
     *              will also be present in each bam_pileup1_t created.
     */
    void bam_plp_constructor(bam_plp_t plp,
                             int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
    void bam_plp_destructor(bam_plp_t plp,
                            int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));

    /// Get pileup padded insertion sequence
    /**
     * @param p       pileup data
     * @param ins     the kstring where the insertion sequence will be written
     * @param del_len location for deletion length
     * @return the length of insertion string on success; -1 on failure.
     *
     * Fills out the kstring with the padded insertion sequence for the current
     * location in 'p'.  If this is not an insertion site, the string is blank.
     *
     * If del_len is not NULL, the location pointed to is set to the length of
     * any deletion immediately following the insertion, or zero if none.
     */
    int bam_plp_insertion(const bam_pileup1_t *p, kstring_t *ins, int *del_len) HTS_RESULT_USED;

    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
    /// Set up mpileup overlap detection
    /**
     * @param iter    mpileup iterator
     * @return 0 on success; a negative value on error
     *
     *  If called, mpileup will detect overlapping
     *  read pairs and for each base pair set the base quality of the
     *  lower-quality base to zero, thus effectively discarding it from
     *  calling. If the two bases are identical, the quality of the other base
     *  is increased to the sum of their qualities (capped at 200), otherwise
     *  it is multiplied by 0.8.
     */
    int bam_mplp_init_overlaps(bam_mplp_t iter);
    void bam_mplp_destroy(bam_mplp_t iter);
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);
    void bam_mplp_reset(bam_mplp_t iter);
    void bam_mplp_constructor(bam_mplp_t iter,
                              int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
    void bam_mplp_destructor(bam_mplp_t iter,
                             int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));

#endif // ~!defined(BAM_NO_PILEUP)


/***********************************
 * BAQ calculation and realignment *
 ***********************************/

int sam_cap_mapq(bam1_t *b, const char *ref, int ref_len, int thres);

/// Calculate BAQ scores
/** @param b   BAM record
    @param ref     Reference sequence
    @param ref_len Reference sequence length
    @param flag    Flags, see description
    @return 0 on success \n
           -1 if the read was unmapped, zero length, had no quality values, did not have at least one M, X or = CIGAR operator, or included a reference skip. \n
           -3 if BAQ alignment has already been done and does not need to be applied, or has already been applied. \n
           -4 if alignment failed (most likely due to running out of memory)

This function calculates base alignment quality (BAQ) values using the method
described in "Improving SNP discovery by base alignment quality", Heng Li,
Bioinformatics, Volume 27, Issue 8 (https://doi.org/10.1093/bioinformatics/btr076).

The following @param flag bits can be used:

Bit 0: Adjust the quality values using the BAQ values

 If set, the data in the BQ:Z tag is used to adjust the quality values, and
 the BQ:Z tag is renamed to ZQ:Z.

 If clear, and a ZQ:Z tag is present, the quality values are reverted using
 the data in the tag, and the tag is renamed to BQ:Z.

Bit 1: Use "extended" BAQ.

 Changes the BAQ calculation to increase sensitivity at the expense of
 reduced specificity.

Bit 2: Recalculate BAQ, even if a BQ tag is present.

 Force BAQ to be recalculated.  Note that a ZQ:Z tag will always disable
 recalculation.

@bug
If the input read has both BQ:Z and ZQ:Z tags, the ZQ:Z one will be removed.
Depending on what previous processing happened, this may or may not be the
correct thing to do.  It would be wise to avoid this situation if possible.
*/

int sam_prob_realn(bam1_t *b, const char *ref, int ref_len, int flag);

#ifdef __cplusplus
}
#endif

#endif
