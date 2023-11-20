/*
Copyright (c) 2012-2016, 2018-2020, 2023 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

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

#ifndef HTSLIB_CRAM_STRUCTS_H
#define HTSLIB_CRAM_STRUCTS_H

/*
 * Defines in-memory structs for the basic file-format objects in the
 * CRAM format.
 *
 * The basic file format is:
 *     File-def SAM-hdr Container Container ...
 *
 * Container:
 *     Service-block data-block data-block ...
 *
 * Multiple blocks in a container are grouped together as slices,
 * also sometimes referred to as landmarks in the spec.
 */


#include <pthread.h>
#include <stdint.h>
#include <sys/types.h>

#include "../htslib/thread_pool.h"
#include "../htslib/cram.h"
#include "string_alloc.h"
#include "mFILE.h"
#include "../htslib/khash.h"

#ifdef __cplusplus
extern "C" {
#endif

// Generic hash-map integer -> integer
KHASH_MAP_INIT_INT64(m_i2i, int)

// Generic hash-set integer -> (existence)
KHASH_SET_INIT_INT(s_i2i)

// For brevity
typedef unsigned char uc;

/*
 * A union for the preservation map. Required for khash.
 */
typedef union {
    int i;
    char *p;
} pmap_t;

// Generates static functions here which isn't ideal, but we have no way
// currently to declare the kh_map_t structure here without also declaring a
// duplicate in the .c files due to the nature of the KHASH macros.
KHASH_MAP_INIT_STR(map, pmap_t)

struct hFILE;

#define SEQS_PER_SLICE 10000
#define BASES_PER_SLICE (SEQS_PER_SLICE*500)
#define SLICE_PER_CNT  1

#define CRAM_SUBST_MATRIX "CGTNGTANCATNGCANACGT"

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct cram_stats {
    int freqs[MAX_STAT_VAL];
    khash_t(m_i2i) *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
    int64_t min_val, max_val;
} cram_stats;

/* NB: matches java impl, not the spec */
enum cram_encoding {
    E_NULL               = 0,
    E_EXTERNAL           = 1,  // Only for BYTE type in CRAM 4
    E_GOLOMB             = 2,  // Not in CRAM 4
    E_HUFFMAN            = 3,  // Not in CRAM 4
    E_BYTE_ARRAY_LEN     = 4,
    E_BYTE_ARRAY_STOP    = 5,
    E_BETA               = 6,  // Not in CRAM 4
    E_SUBEXP             = 7,  // Not in CRAM 4
    E_GOLOMB_RICE        = 8,  // Not in CRAM 4
    E_GAMMA              = 9,  // Not in CRAM 4

    // CRAM 4 specific codecs
    E_VARINT_UNSIGNED    = 41, // Specialisation of EXTERNAL
    E_VARINT_SIGNED      = 42, // Specialisation of EXTERNAL
    E_CONST_BYTE         = 43, // Alternative to HUFFMAN with 1 symbol
    E_CONST_INT          = 44, // Alternative to HUFFMAN with 1 symbol

    // More experimental ideas, not documented in spec yet
    E_XHUFFMAN           = 50, // To external block
    E_XPACK              = 51, // Transform to sub-codec
    E_XRLE               = 52, // Transform to sub-codec
    E_XDELTA             = 53, // Transform to sub-codec

    // Total number of codecs, not a real one.
    E_NUM_CODECS,
};

enum cram_external_type {
    E_INT                = 1,
    E_LONG               = 2,
    E_BYTE               = 3,
    E_BYTE_ARRAY         = 4,
    E_BYTE_ARRAY_BLOCK   = 5,
    E_SINT               = 6, // signed INT
    E_SLONG              = 7, // signed LONG
};

/* External IDs used by this implementation (only assumed during writing) */
enum cram_DS_ID {
    DS_CORE   = 0,
    DS_aux    = 1, // aux_blk
    DS_aux_OQ = 2,
    DS_aux_BQ = 3,
    DS_aux_BD = 4,
    DS_aux_BI = 5,
    DS_aux_FZ = 6, // also ZM:B
    DS_aux_oq = 7, // other qualities
    DS_aux_os = 8, // other sequences
    DS_aux_oz = 9, // other strings
    DS_ref,
    DS_RN, // name_blk
    DS_QS, // qual_blk
    DS_IN, // base_blk
    DS_SC, // soft_blk

    DS_BF, // start loop
    DS_CF,
    DS_AP,
    DS_RG,
    DS_MQ,
    DS_NS,
    DS_MF,
    DS_TS,
    DS_NP,
    DS_NF,
    DS_RL,
    DS_FN,
    DS_FC,
    DS_FP,
    DS_DL,
    DS_BA,
    DS_BS,
    DS_TL,
    DS_RI,
    DS_RS,
    DS_PD,
    DS_HC,
    DS_BB,
    DS_QQ,

    DS_TN, // end loop

    DS_RN_len,
    DS_SC_len,
    DS_BB_len,
    DS_QQ_len,

    DS_TC, // CRAM v1.0 tags
    DS_TM, // test
    DS_TV, // test

    DS_END,
};

/* "File Definition Structure" */
struct cram_file_def {
    char    magic[4];
    uint8_t major_version;
    uint8_t minor_version;
    char    file_id[20] HTS_NONSTRING; // Filename or SHA1 checksum
};

#define CRAM_MAJOR_VERS(v) ((v) >> 8)
#define CRAM_MINOR_VERS(v) ((v) & 0xff)

struct cram_slice;

// Internal version of htslib/cram.h enum.
// Note these have to match the laout of methmap and methcost in
// cram_io.c:cram_compress_block2
enum cram_block_method_int {
    // Public methods as defined in the CRAM spec.
    BM_ERROR = -1,

    // CRAM 2.x and 3.0
    RAW      = 0,
    GZIP     = 1,
    BZIP2    = 2,
    LZMA     = 3,
    RANS     = 4, RANS0 = RANS,

    // CRAM 3.1 onwards
    RANSPR   = 5, RANS_PR0  = RANSPR,
    ARITH    = 6, ARITH_PR0 = ARITH,
    FQZ      = 7,
    TOK3     = 8,
    // BSC = 9, ZSTD = 10

    // Methods not externalised, but used in metrics.
    // Externally they become one of the above methods.
    GZIP_RLE = 11,
    GZIP_1,      // Z_DEFAULT_STRATEGY level 1, NB: not externalised in CRAM

    FQZ_b, FQZ_c, FQZ_d, // Various preset FQZ methods

  //RANS0,       // Order 0
    RANS1,

  //RANS_PR0,    // Order 0
    RANS_PR1,    // Order 1
    RANS_PR64,   // O0 + RLE
    RANS_PR9,    // O1 + X4
    RANS_PR128,  // O0 + Pack
    RANS_PR129,  // O1 + Pack
    RANS_PR192,  // O0 + RLE + pack
    RANS_PR193,  // O1 + RLE + pack

  //TOK3,   // tok+rans
    TOKA,   // tok+arith

  //ARITH_PR0,   // Order 0
    ARITH_PR1,   // Order 1
    ARITH_PR64,  // O0 + RLE
    ARITH_PR9,   // O1 + X4
    ARITH_PR128, // O0 + Pack
    ARITH_PR129, // O1 + Pack
    ARITH_PR192, // O0 + RLE + pack
    ARITH_PR193, // O1 + RLE + pack

    // NB: must end on no more than 31 unless we change to a
    // 64-bit method type.
};

/* Now in htslib/cram.h
enum cram_content_type {
    CT_ERROR           = -1,
    FILE_HEADER        = 0,
    COMPRESSION_HEADER = 1,
    MAPPED_SLICE       = 2,
    UNMAPPED_SLICE     = 3, // CRAM V1.0 only
    EXTERNAL           = 4,
    CORE               = 5,
};
*/

/* Maximum simultaneous codecs allowed, 1 per bit */
#define CRAM_MAX_METHOD 32

/* Compression metrics */
struct cram_metrics {
    // number of trials and time to next trial
    int trial;
    int next_trial;
    int consistency;

    // aggregate sizes during trials
    int sz[CRAM_MAX_METHOD];
    int input_avg_sz, input_avg_delta;

    // resultant method from trials
    int method, revised_method;
    int strat;

    // Revisions of method, to allow culling of continually failing ones.
    int cnt[CRAM_MAX_METHOD];

    double extra[CRAM_MAX_METHOD];

    // Not amenable to rANS bit-packing techniques; cardinality > 16
    int unpackable;
};

// Hash aux key (XX:i) to cram_metrics
KHASH_MAP_INIT_INT(m_metrics, cram_metrics*)


/* Block */
struct cram_block {
    enum cram_block_method_int  method, orig_method;
    enum cram_content_type  content_type;
    int32_t  content_id;
    int32_t  comp_size;
    int32_t  uncomp_size;
    uint32_t crc32;
    int32_t  idx; /* offset into data */
    unsigned char    *data;

    // For bit I/O
    size_t alloc;
    size_t byte;
    int bit;

    // To aid compression
    cram_metrics *m; // used to track aux block compression only

    int crc32_checked;
    uint32_t crc_part;
};

struct cram_codec; /* defined in cram_codecs.h */
struct cram_map;

#define CRAM_MAP_HASH 32
#define CRAM_MAP(a,b) (((a)*3+(b))&(CRAM_MAP_HASH-1))

/* Compression header block */
struct cram_block_compression_hdr {
    int32_t ref_seq_id;
    int64_t ref_seq_start;
    int64_t ref_seq_span;
    int32_t num_records;
    int32_t num_landmarks;
    int32_t *landmark;

    /* Flags from preservation map */
    int read_names_included;
    int AP_delta;
    // indexed by ref-base and subst. code
    char substitution_matrix[5][4];
    int no_ref;
    int qs_seq_orient; // 1 => same as seq. 0 => original orientation

    // TD Dictionary as a concatenated block
    cram_block *TD_blk;          // Tag Dictionary
    int nTL;                     // number of TL entries in TD
    unsigned char **TL;          // array of size nTL, pointer into TD_blk.
    khash_t(m_s2i) *TD_hash;     // Keyed on TD strings, map to TL[] indices
    string_alloc_t *TD_keys;     // Pooled keys for TD hash.

    khash_t(map) *preservation_map;
    struct cram_map *rec_encoding_map[CRAM_MAP_HASH];
    struct cram_map *tag_encoding_map[CRAM_MAP_HASH];

    struct cram_codec *codecs[DS_END];

    char *uncomp; // A single block of uncompressed data
    size_t uncomp_size, uncomp_alloc;

    // Total codec count, used for index to block_by_id for transforms
    int ncodecs;
};

typedef struct cram_map {
    int key;    /* 0xe0 + 3 bytes */
    enum cram_encoding encoding;
    int offset; /* Offset into a single block of memory */
    int size;   /* Size */
    struct cram_codec *codec;
    struct cram_map *next; // for noddy internal hash
} cram_map;

typedef struct cram_tag_map {
    struct cram_codec *codec;
    cram_block *blk;
    cram_block *blk2;
    cram_metrics *m;
} cram_tag_map;

// Hash aux key (XX:i) to cram_tag_map
KHASH_MAP_INIT_INT(m_tagmap, cram_tag_map*)

/* Mapped or unmapped slice header block */
struct cram_block_slice_hdr {
    enum cram_content_type content_type;
    int32_t ref_seq_id;     /* if content_type == MAPPED_SLICE */
    int64_t ref_seq_start;  /* if content_type == MAPPED_SLICE */
    int64_t ref_seq_span;   /* if content_type == MAPPED_SLICE */
    int32_t num_records;
    int64_t record_counter;
    int32_t num_blocks;
    int32_t num_content_ids;
    int32_t *block_content_ids;
    int32_t ref_base_id;    /* if content_type == MAPPED_SLICE */
    unsigned char md5[16];
};

struct ref_entry;

/*
 * Container.
 *
 * Conceptually a container is split into slices, and slices into blocks.
 * However on disk it's just a list of blocks and we need to query the
 * block types to identify the start/end points of the slices.
 *
 * OR... are landmarks the start/end points of slices?
 */
struct cram_container {
    int32_t  length;
    int32_t  ref_seq_id;
    int64_t  ref_seq_start;
    int64_t  ref_seq_span;
    int64_t  record_counter;
    int64_t  num_bases;
    int32_t  num_records;
    int32_t  num_blocks;
    int32_t  num_landmarks;
    int32_t *landmark;

    /* Size of container header above */
    size_t   offset;

    /* Compression header is always the first block? */
    cram_block_compression_hdr *comp_hdr;
    cram_block *comp_hdr_block;

    /* For construction purposes */
    int max_slice, curr_slice;   // maximum number of slices
    int curr_slice_mt;           // Curr_slice when reading ahead (via threads)
    int max_rec, curr_rec;       // current and max recs per slice
    int max_c_rec, curr_c_rec;   // current and max recs per container
    int slice_rec;               // rec no. for start of this slice
    int curr_ref;                // current ref ID. -2 for no previous
    int64_t last_pos;                // last record position
    struct cram_slice **slices, *slice;
    int pos_sorted;              // boolean, 1=>position sorted data
    int64_t max_apos;                // maximum position, used if pos_sorted==0
    int last_slice;              // number of reads in last slice (0 for 1st)
    int multi_seq;               // true if packing multi seqs per cont/slice
    int unsorted;                // true is AP_delta is 0.
    int qs_seq_orient;           // 1 => same as seq. 0 => original orientation

    /* Copied from fd before encoding, to allow multi-threading */
    int ref_id;
    hts_pos_t ref_start, first_base, last_base, ref_end;
    char *ref;
    int embed_ref;               // 1 if embedding ref, 2 if embedding cons
    int no_ref;                  // true if referenceless
    //struct ref_entry *ref;

    /* For multi-threading */
    bam_seq_t **bams;

    /* Statistics for encoding */
    cram_stats *stats[DS_END];

    khash_t(m_tagmap) *tags_used; // set of tag types in use, for tag encoding map
    int *refs_used;       // array of frequency of ref seq IDs

    uint32_t crc32;       // CRC32

    uint64_t s_num_bases; // number of bases in this slice
    uint64_t s_aux_bytes; // number of bytes of aux in BAM

    uint32_t n_mapped;    // Number of mapped reads
    int ref_free;         // whether 'ref' is owned by us and must be freed.
};

/*
 * A single cram record
 */
typedef struct cram_record {
    struct cram_slice *s; // Filled out by cram_decode only

    int32_t ref_id;       // fixed for all recs in slice?
    int32_t flags;        // BF
    int32_t cram_flags;   // CF
    int32_t len;          // RL
    int64_t apos;         // AP
    int32_t rg;           // RG
    int32_t name;         // RN; idx to s->names_blk
    int32_t name_len;
    int32_t mate_line;    // index to another cram_record
    int32_t mate_ref_id;
    int64_t mate_pos;     // NP
    int64_t tlen;         // TS
    int64_t explicit_tlen;// TS, but PNEXT/RNEXT still need auto-computing

    // Auxiliary data
    int32_t ntags;        // TC
    uint32_t aux;         // idx to s->aux_blk
    uint32_t aux_size;    // total size of packed ntags in aux_blk
#ifndef TN_external
    int32_t TN_idx;       // TN; idx to s->TN;
#else
    int32_t tn;           // idx to s->tn_blk
#endif
    int     TL;

    uint32_t seq;         // idx to s->seqs_blk
    uint32_t qual;        // idx to s->qual_blk
    uint32_t cigar;       // idx to s->cigar
    int32_t ncigar;
    int64_t aend;         // alignment end
    int32_t mqual;        // MQ

    uint32_t feature;     // idx to s->feature
    uint32_t nfeature;    // number of features
    int32_t mate_flags;   // MF
} cram_record;

// Accessor macros as an analogue of the bam ones
#define cram_qname(c)    (&(c)->s->name_blk->data[(c)->name])
#define cram_seq(c)      (&(c)->s->seqs_blk->data[(c)->seq])
#define cram_qual(c)     (&(c)->s->qual_blk->data[(c)->qual])
#define cram_aux(c)      (&(c)->s->aux_blk->data[(c)->aux])
#define cram_seqi(c,i)   (cram_seq((c))[(i)])
#define cram_name_len(c) ((c)->name_len)
#define cram_strand(c)   (((c)->flags & BAM_FREVERSE) != 0)
#define cram_mstrand(c)  (((c)->flags & BAM_FMREVERSE) != 0)
#define cram_cigar(c)    (&((cr)->s->cigar)[(c)->cigar])

/*
 * A feature is a base difference, used for the sequence reference encoding.
 * (We generate these internally when writing CRAM.)
 */
typedef union cram_feature {
    struct {
        int pos;
        int code;
        int base;    // substitution code
    } X;
    struct {
        int pos;
        int code;
        int base;    // actual base & qual
        int qual;
    } B;
    struct {
        int pos;
        int code;
        int seq_idx; // index to s->seqs_blk
        int len;
    } b;
    struct {
        int pos;
        int code;
        int qual;
    } Q;
    struct {
        int pos;
        int code;
        int len;
        int seq_idx; // soft-clip multiple bases
    } S;
    struct {
        int pos;
        int code;
        int len;
        int seq_idx; // insertion multiple bases
    } I;
    struct {
        int pos;
        int code;
        int base; // insertion single base
    } i;
    struct {
        int pos;
        int code;
        int len;
    } D;
    struct {
        int pos;
        int code;
        int len;
    } N;
    struct {
        int pos;
        int code;
        int len;
    } P;
    struct {
        int pos;
        int code;
        int len;
    } H;
} cram_feature;

/*
 * A slice is really just a set of blocks, but it
 * is the logical unit for decoding a number of
 * sequences.
 */
struct cram_slice {
    cram_block_slice_hdr *hdr;
    cram_block *hdr_block;
    cram_block **block;
    cram_block **block_by_id;

    /* State used during encoding/decoding */
    int64_t last_apos, max_apos;

    /* Array of decoded cram records */
    cram_record *crecs;

    /* An dynamically growing buffers for data pointed
     * to by crecs[] array.
     */
    uint32_t  *cigar;
    uint32_t   cigar_alloc;
    uint32_t   ncigar;

    cram_feature *features;
    uint32_t      nfeatures;
    uint32_t      afeatures; // allocated size of features

#ifndef TN_external
    // TN field (Tag Name)
    uint32_t      *TN;
    int           nTN, aTN;  // used and allocated size for TN[]
#else
    cram_block *tn_blk;
    int tn_id;
#endif

    // For variable sized elements which are always external blocks.
    cram_block *name_blk;
    cram_block *seqs_blk;
    cram_block *qual_blk;
    cram_block *base_blk;
    cram_block *soft_blk;
    cram_block *aux_blk;       // BAM aux block, created while decoding CRAM

    string_alloc_t *pair_keys; // Pooled keys for pair hash.
    khash_t(m_s2i) *pair[2];   // for identifying read-pairs in this slice.

    char *ref;                 // slice of current reference
    hts_pos_t ref_start;       // start position of current reference;
    hts_pos_t ref_end;         // end position of current reference;
    int ref_id;

    // For going from BAM to CRAM; an array of auxiliary blocks per type
    int naux_block;
    cram_block **aux_block;

    unsigned int data_series; // See cram_fields enum
    int decode_md;

    int max_rec, curr_rec;       // current and max recs per slice
    int slice_num;               // To be copied into c->curr_slice in decode
};

/*-----------------------------------------------------------------------------
 * Consider moving reference handling to cram_refs.[ch]
 */
// from fa.fai / samtools faidx files
typedef struct ref_entry {
    char *name;
    char *fn;
    int64_t length;
    int64_t offset;
    int bases_per_line;
    int line_length;
    int64_t count;         // for shared references so we know to dealloc seq
    char *seq;
    mFILE *mf;
    int is_md5;            // Reference comes from a raw seq found by MD5
    int validated_md5;
} ref_entry;

KHASH_MAP_INIT_STR(refs, ref_entry*)

// References structure.
struct refs_t {
    string_alloc_t *pool;  // String pool for holding filenames and SN vals

    khash_t(refs) *h_meta; // ref_entry*, index by name
    ref_entry **ref_id;    // ref_entry*, index by ID
    int nref;              // number of ref_entry

    char *fn;              // current file opened
    BGZF *fp;              // and the hFILE* to go with it.

    int count;             // how many cram_fd sharing this refs struct

    pthread_mutex_t lock;  // Mutex for multi-threaded updating
    ref_entry *last;       // Last queried sequence
    int last_id;           // Used in cram_ref_decr_locked to delay free
};

/*-----------------------------------------------------------------------------
 * CRAM index
 *
 * Detect format by number of entries per line.
 * 5 => 1.0 (refid, start, nseq, C offset, slice)
 * 6 => 1.1 (refid, start, span, C offset, S offset, S size)
 *
 * Indices are stored in a nested containment list, which is trivial to set
 * up as the indices are on sorted data so we're appending to the nclist
 * in sorted order. Basically if a slice entirely fits within a previous
 * slice then we append to that slices list. This is done recursively.
 *
 * Lists are sorted on two dimensions: ref id + slice coords.
 */
typedef struct cram_index {
    int nslice, nalloc;   // total number of slices
    struct cram_index *e; // array of size nslice

    int     refid;  // 1.0                 1.1
    int     start;  // 1.0                 1.1
    int     end;    //                     1.1
    int     nseq;   // 1.0 - undocumented
    int     slice;  // 1.0 landmark index, 1.1 landmark value
    int     len;    //                     1.1 - size of slice in bytes
    int64_t offset; // 1.0                 1.1

    // Linked list of cram_index entries. Used to convert recursive
    // NCList back to a linear list.
    struct cram_index *e_next;
} cram_index;

typedef struct {
    int refid;
    int64_t start;
    int64_t end;
} cram_range;

/*-----------------------------------------------------------------------------
 */
/* CRAM File handle */

typedef struct spare_bams {
    bam_seq_t **bams;
    struct spare_bams *next;
} spare_bams;

struct cram_fd;
typedef struct varint_vec {
    // Returns number of bytes decoded from fd, 0 on error
    int (*varint_decode32_crc)(struct cram_fd *fd, int32_t *val_p, uint32_t *crc);
    int (*varint_decode32s_crc)(struct cram_fd *fd, int32_t *val_p, uint32_t *crc);
    int (*varint_decode64_crc)(struct cram_fd *fd, int64_t *val_p, uint32_t *crc);

    // Returns the value and increments *cp.  Sets err to 1 iff an error occurs.
    // NOTE: Does not set err to 0 on success.
    int64_t (*varint_get32) (char **cp, const char *endp, int *err);
    int64_t (*varint_get32s)(char **cp, const char *endp, int *err);
    int64_t (*varint_get64) (char **cp, const char *endp, int *err);
    int64_t (*varint_get64s)(char **cp, const char *endp, int *err);

    // Returns the number of bytes written, <= 0 on error.
    int (*varint_put32) (char *cp, char *endp, int32_t val_p);
    int (*varint_put32s)(char *cp, char *endp, int32_t val_p);
    int (*varint_put64) (char *cp, char *endp, int64_t val_p);
    int (*varint_put64s)(char *cp, char *endp, int64_t val_p);

    // Returns the number of bytes written, <= 0 on error.
    int (*varint_put32_blk) (cram_block *blk, int32_t val_p);
    int (*varint_put32s_blk)(cram_block *blk, int32_t val_p);
    int (*varint_put64_blk) (cram_block *blk, int64_t val_p);
    int (*varint_put64s_blk)(cram_block *blk, int64_t val_p);

    // Returns number of bytes needed to encode 'val'
    int (*varint_size)(int64_t val);
} varint_vec;

struct cram_fd {
    struct hFILE  *fp;
    int            mode;     // 'r' or 'w'
    int            version;
    cram_file_def *file_def;
    sam_hdr_t     *header;

    char          *prefix;
    int64_t        record_counter;
    int            err;

    // Most recent compression header decoded
    //cram_block_compression_hdr *comp_hdr;
    //cram_block_slice_hdr       *slice_hdr;

    // Current container being processed
    cram_container *ctr;

    // Current container used for decoder threads
    cram_container *ctr_mt;

    // positions for encoding or decoding
    int first_base, last_base; // copied to container

    // cached reference portion
    refs_t   *refs;                // ref meta-data structure
    char     *ref, *ref_free;      // current portion held in memory
    int       ref_id;              // copied to container
    hts_pos_t ref_start;           // copied to container
    hts_pos_t ref_end;             // copied to container
    char     *ref_fn;              // reference fasta filename

    // compression level and metrics
    int level;
    cram_metrics *m[DS_END];
    khash_t(m_metrics) *tags_used; // cram_metrics[], per tag types in use.

    // options
    int decode_md; // Whether to export MD and NM tags
    int seqs_per_slice;
    int bases_per_slice;
    int slices_per_container;
    int embed_ref; // copied to container
    int no_ref;    // copied to container
    int no_ref_counter; // decide if permanent switch
    int ignore_md5;
    int use_bz2;
    int use_rans;
    int use_lzma;
    int use_fqz;
    int use_tok;
    int use_arith;
    int shared_ref;
    unsigned int required_fields;
    int store_md;
    int store_nm;
    cram_range range;

    // lookup tables, stored here so we can be trivially multi-threaded
    unsigned int bam_flag_swap[0x1000]; // cram -> bam flags
    unsigned int cram_flag_swap[0x1000];// bam -> cram flags
    unsigned char L1[256];              // ACGT{*} ->0123{4}
    unsigned char L2[256];              // ACGTN{*}->01234{5}
    char cram_sub_matrix[32][32];       // base substitution codes

    int         index_sz;
    cram_index *index;                  // array, sizeof index_sz
    off_t first_container;
    off_t curr_position;
    int eof;
    int last_slice;                     // number of recs encoded in last slice
    int last_RI_count;                  // number of references encoded in last container
    int multi_seq;                      // -1 is auto, 0 is one ref per container, 1 is multi...
    int multi_seq_user;                 // Original user setting (CRAM_OPT_MULTI_SEQ_PER_SLICE)
    int unsorted;
    int last_mapped;                    // number of mapped reads in last container
    int empty_container;                // Marker for EOF block

    // thread pool
    int own_pool;
    hts_tpool *pool;
    hts_tpool_process *rqueue;
    pthread_mutex_t metrics_lock;
    pthread_mutex_t ref_lock;
    pthread_mutex_t range_lock;
    spare_bams *bl;
    pthread_mutex_t bam_list_lock;
    void *job_pending;
    int ooc;                            // out of containers.

    int lossy_read_names;               // boolean
    int tlen_approx;                    // max TLEN calculation offset.
    int tlen_zero;                      // If true, permit tlen 0 (=> tlen calculated)

    BGZF *idxfp;                        // File pointer for on-the-fly index creation

    // variable integer decoding callbacks.
    // This changed in CRAM4.0 to a data-size agnostic encoding.
    varint_vec vv;

    // Force AP delta even on non positional sorted data.
    // This can be beneficial for pairs where pairs are nearby each other.
    // We suffer with delta to unrelated things (previous pair), but gain
    // in delta between them.  (Ideal would be a per read setting.)
    int ap_delta;
};

// Translation of required fields to cram data series
enum cram_fields {
    CRAM_BF = 0x00000001,
    CRAM_AP = 0x00000002,
    CRAM_FP = 0x00000004,
    CRAM_RL = 0x00000008,
    CRAM_DL = 0x00000010,
    CRAM_NF = 0x00000020,
    CRAM_BA = 0x00000040,
    CRAM_QS = 0x00000080,
    CRAM_FC = 0x00000100,
    CRAM_FN = 0x00000200,
    CRAM_BS = 0x00000400,
    CRAM_IN = 0x00000800,
    CRAM_RG = 0x00001000,
    CRAM_MQ = 0x00002000,
    CRAM_TL = 0x00004000,
    CRAM_RN = 0x00008000,
    CRAM_NS = 0x00010000,
    CRAM_NP = 0x00020000,
    CRAM_TS = 0x00040000,
    CRAM_MF = 0x00080000,
    CRAM_CF = 0x00100000,
    CRAM_RI = 0x00200000,
    CRAM_RS = 0x00400000,
    CRAM_PD = 0x00800000,
    CRAM_HC = 0x01000000,
    CRAM_SC = 0x02000000,
    CRAM_BB = 0x04000000,
    CRAM_BB_len = 0x08000000,
    CRAM_QQ = 0x10000000,
    CRAM_QQ_len = 0x20000000,
    CRAM_aux= 0x40000000,
    CRAM_ALL= 0x7fffffff,
};

// A CIGAR opcode, but not necessarily the implications of it. Eg FC/FP may
// encode a base difference, but we don't need to know what it is for CIGAR.
// If we have a soft-clip or insertion, we do need SC/IN though to know how
// long that array is.
#define CRAM_CIGAR (CRAM_FN | CRAM_FP | CRAM_FC | CRAM_DL | CRAM_IN | \
                    CRAM_SC | CRAM_HC | CRAM_PD | CRAM_RS | CRAM_RL | CRAM_BF)

#define CRAM_SEQ (CRAM_CIGAR | CRAM_BA | CRAM_BS | \
                  CRAM_RL    | CRAM_AP | CRAM_BB)

#define CRAM_QUAL (CRAM_CIGAR | CRAM_RL | CRAM_AP | CRAM_QS | CRAM_QQ)

/* BF bitfields */
/* Corrected in 1.1. Use bam_flag_swap[bf] and BAM_* macros for 1.0 & 1.1 */
#define CRAM_FPAIRED      256
#define CRAM_FPROPER_PAIR 128
#define CRAM_FUNMAP        64
#define CRAM_FREVERSE      32
#define CRAM_FREAD1        16
#define CRAM_FREAD2         8
#define CRAM_FSECONDARY     4
#define CRAM_FQCFAIL        2
#define CRAM_FDUP           1

#define DS_aux_S "\001"
#define DS_aux_OQ_S "\002"
#define DS_aux_BQ_S "\003"
#define DS_aux_BD_S "\004"
#define DS_aux_BI_S "\005"
#define DS_aux_FZ_S "\006"
#define DS_aux_oq_S "\007"
#define DS_aux_os_S "\010"
#define DS_aux_oz_S "\011"

#define CRAM_M_REVERSE  1
#define CRAM_M_UNMAP    2


/* CF bitfields */
#define CRAM_FLAG_PRESERVE_QUAL_SCORES (1<<0)
#define CRAM_FLAG_DETACHED             (1<<1)
#define CRAM_FLAG_MATE_DOWNSTREAM      (1<<2)
#define CRAM_FLAG_NO_SEQ               (1<<3)
#define CRAM_FLAG_EXPLICIT_TLEN        (1<<4)
#define CRAM_FLAG_MASK                 ((1<<5)-1)

/* Internal only */
#define CRAM_FLAG_STATS_ADDED          (1<<30)
#define CRAM_FLAG_DISCARD_NAME         (1U<<31)

#ifdef __cplusplus
}
#endif

#endif /* HTSLIB_CRAM_STRUCTS_H */
