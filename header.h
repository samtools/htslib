/*
Copyright (c) 2013-2018 Genome Research Ltd.
Authors: James Bonfield <jkb@sanger.ac.uk>, Valeriu Ohan <vo2@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! \file
 * SAM header parsing.
 *
 * These functions can be shared between SAM, BAM and CRAM file
 * formats as all three internally use the same string encoding for
 * header fields.
 */


#ifndef HEADER_H_
#define HEADER_H_

#include <stdarg.h>

#include "cram/string_alloc.h"
#include "cram/pooled_alloc.h"

#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

// For structure assignment. Eg kstring_t s = KS_INITIALIZER;
#define KS_INITIALIZER {0,0,0}

// For initialisation elsewhere. Eg KS_INIT(x->str);
#define KS_INIT(ks) ((ks)->l = 0, (ks)->m = 0, (ks)->s = NULL)

// Frees the string subfield only. Assumes 's' itself is static.
#define KS_FREE(ks) do { if ((ks)->s) {free((ks)->s); (ks)->s = NULL;} } while(0)

#define K(a) (((a)[0]<<8)|((a)[1]))

#define SAM_HDR_LINES 32

/*
 * Proposed new SAM header parsing

1 @SQ ID:foo LN:100
2 @SQ ID:bar LN:200
3 @SQ ID:ram LN:300 UR:xyz
4 @RG ID:r ...
5 @RG ID:s ...

Hash table for 2-char @keys without dup entries.
If dup lines, we form a circular linked list. Ie hash keys = {RG, SQ}.

HASH("SQ")--\
            |
    (3) <-> 1 <-> 2 <-> 3 <-> (1)

HASH("RG")--\
            |
    (5) <-> 4 <-> 5 <-> (4)

Items stored in the hash values also form their own linked lists:
Ie SQ->ID(foo)->LN(100)
   SQ->ID(bar)->LN(200)
   SQ->ID(ram)->LN(300)->UR(xyz)
   RG->ID(r)
 */

/*! A single key:value pair on a header line
 *
 * These form a linked list and hold strings. The strings are
 * allocated from a string_alloc_t pool referenced in the master
 * bam_hrecs_t structure. Do not attempt to free, malloc or manipulate
 * these strings directly.
 */
typedef struct bam_hrec_tag_s {
    struct bam_hrec_tag_s *next;
    char *str;
    int   len;
} bam_hrec_tag_t;

/*! The parsed version of the SAM header string.
 *
 * Each header type (SQ, RG, HD, etc) points to its own sam_hdr_type
 * struct via the main hash table h in the bam_hrecs_t struct.
 *
 * These in turn consist of circular bi-directional linked lists (ie
 * rings) to hold the multiple instances of the same header type
 * code. For example if we have 5 \@SQ lines the primary hash table
 * will key on \@SQ pointing to the first sam_hdr_type and that in turn
 * will be part of a ring of 5 elements.
 *
 * For each sam_hdr_type structure we also point to a sam_hdr_tag
 * structure which holds the tokenised attributes; the tab separated
 * key:value pairs per line.
 */
typedef struct bam_hrec_type_s {
    struct bam_hrec_type_s *next; // circular list of this type
    struct bam_hrec_type_s *prev; // circular list of this type
    struct bam_hrec_type_s *global_next; // circular list of all lines
    struct bam_hrec_type_s *global_prev; // circular list of all lines
    bam_hrec_tag_t *tag;          // first tag
    khint32_t type;               // Two-letter type code as an int
} bam_hrec_type_t;

/*! Parsed \@SQ lines */
typedef struct {
    char *name;
    uint32_t len;
    bam_hrec_type_t *ty;
} bam_hrec_sq_t;

/*! Parsed \@RG lines */
typedef struct {
    char *name;
    bam_hrec_type_t *ty;
    int name_len;
    int id;           // numerical ID
} bam_hrec_rg_t;

/*! Parsed \@PG lines */
typedef struct {
    char *name;
    bam_hrec_type_t *ty;
    int name_len;
    int id;           // numerical ID
    int prev_id;      // -1 if none
} bam_hrec_pg_t;


/*! Sort order parsed from @HD line */
enum sam_sort_order {
    ORDER_UNKNOWN  =-1,
    ORDER_UNSORTED = 0,
    ORDER_NAME     = 1,
    ORDER_COORD    = 2
  //ORDER_COLLATE  = 3 // maybe one day!
};

enum sam_group_order {
    ORDER_NONE      =-1,
    ORDER_QUERY     = 0,
    ORDER_REFERENCE = 1
};

KHASH_MAP_INIT_INT(bam_hrecs_t, bam_hrec_type_t*)
KHASH_MAP_INIT_STR(m_s2i, int)

/*! Primary structure for header manipulation
 *
 * The initial header text is held in the text kstring_t, but is also
 * parsed out into SQ, RG and PG arrays. These have a hash table
 * associated with each to allow lookup by ID or SN fields instead of
 * their numeric array indices. Additionally PG has an array to hold
 * the linked list start points (the last in a PP chain).
 *
 * Use the appropriate sam_hdr_* functions to edit the header, and
 * call sam_hdr_rebuild() any time the textual form needs to be
 * updated again.
 */
struct sam_hdr {
    khash_t(bam_hrecs_t) *h;
    bam_hrec_type_t *first_line; //!< First line (usually @HD)
    string_alloc_t *str_pool; //!< Pool of sam_hdr_tag->str strings
    pool_alloc_t   *type_pool;//!< Pool of sam_hdr_type structs
    pool_alloc_t   *tag_pool; //!< Pool of sam_hdr_tag structs

    // @SQ lines / references
    int nref;                  //!< Number of \@SQ lines
    int ref_sz;                //!< Number of entries available in ref[]
    bam_hrec_sq_t *ref;        //!< Array of parsed \@SQ lines
    khash_t(m_s2i) *ref_hash;  //!< Maps SQ SN field to ref[] index

    // @RG lines / read-groups
    int nrg;                   //!< Number of \@RG lines
    int rg_sz;                 //!< number of entries available in rg[]
    bam_hrec_rg_t *rg;         //!< Array of parsed \@RG lines
    khash_t(m_s2i) *rg_hash;   //!< Maps RG ID field to rg[] index

    // @PG lines / programs
    int npg;                   //!< Number of \@PG lines
    int pg_sz;                //!< Number of entries available in pg[]
    int npg_end;               //!< Number of terminating \@PG lines
    int npg_end_alloc;         //!< Size of pg_end field
    bam_hrec_pg_t *pg;         //!< Array of parsed \@PG lines
    khash_t(m_s2i) *pg_hash;   //!< Maps PG ID field to pg[] index
    int *pg_end;               //!< \@PG chain termination IDs

    // @cond internal
    char *ID_buf;             // temporary buffer for bam_hdr_pg_id
    uint32_t ID_buf_sz;
    int ID_cnt;
    // @endcond

    int dirty;                // marks the header as modified, so it can be rebuilt
    int refs_changed;         // Index of first changed ref (-1 if unchanged)
    int pgs_changed;          // New PG line added
    int type_count;
    char (*type_order)[3];
};

/*!
 * Method for parsing the header text and populating the
 * internal hash tables. After calling this method, the
 * parsed representation becomes the single source of truth.
 *
 * @param bh    Header structure, previously initialised by a
 *              bam_hdr_init call
 * @return      0 on success, -1 on failure
 */
int bam_hdr_parse(bam_hdr_t *bh);

/*!
 * Reconstructs the text representation of the header from
 * the hash table data after a change has been performed on
 * the header.
 *
 * @return  0 on success, -1 on failure
 */
int bam_hdr_rebuild(bam_hdr_t *bh);

/*! Creates an empty SAM header, ready to be populated.
 *
 * @return
 * Returns a bam_hrecs_t struct on success (free with bam_hrecs_free())
 *         NULL on failure
 */
bam_hrecs_t *bam_hrecs_new(void);

/*! Produces a duplicate copy of hrecs and returns it.
 * @return
 * Returns NULL on failure
 */
bam_hrecs_t *bam_hrecs_dup(bam_hrecs_t *hrecs);

/*! Update bam_hdr_t target_name and target_len arrays
 *
 *  bam_hdr_t and bam_hrecs_t are specified separately so that bam_hdr_dup
 *  can use it to construct target arrays from the source header.
 *
 *  @return 0 on success; -1 on failure
 */
int update_target_arrays(bam_hdr_t *bh, const bam_hrecs_t *hrecs,
                         int refs_changed);

/*! Reconstructs a kstring from the header hash table.
 *
 * @return
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hrecs_rebuild_text(const bam_hrecs_t *hrecs, kstring_t *ks);

/*! Deallocates all storage used by a bam_hrecs_t struct.
 *
 * This also decrements the header reference count. If after decrementing
 * it is still non-zero then the header is assumed to be in use by another
 * caller and the free is not done.
 *
 * This is a synonym for sam_hdr_dec_ref().
 */
void bam_hrecs_free(bam_hrecs_t *hrecs);

/*!
 * @return
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
bam_hrec_type_t *bam_hrecs_find_type_id(bam_hrecs_t *hrecs, const char *type,
                                     const char *ID_key, const char *ID_value);

/*
 * Adds or updates tag key,value pairs in a header line.
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line.
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bam_hrecs_update(bam_hrecs_t *hrecs, bam_hrec_type_t *type, va_list ap);

bam_hrec_tag_t *bam_hrecs_find_key(bam_hrec_type_t *type,
                                   const char *key,
                                   bam_hrec_tag_t **prev);

int bam_hrecs_remove_key(bam_hrecs_t *hrecs,
                         bam_hrec_type_t *type,
                         const char *key);

/*! Looks up a read-group by name and returns a pointer to the start of the
 * associated tag list.
 *
 * @return
 * Returns NULL on failure
 */
bam_hrec_rg_t *bam_hrecs_find_rg(bam_hrecs_t *hrecs, const char *rg);

/*! Returns the sort order from the @HD SO: field */
enum sam_sort_order bam_hrecs_sort_order(bam_hrecs_t *hrecs);

/*! Returns the group order from the @HD SO: field */
enum sam_group_order bam_hrecs_group_order(bam_hrecs_t *hrecs);

#ifdef __cplusplus
}
#endif

#endif /* HEADER_H_ */
