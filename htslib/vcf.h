/// @file htslib/vcf.h
/// High-level VCF/BCF variant calling file operations.
/*
    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2020, 2022-2023 Genome Research Ltd.

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

/*
    todo:
        - make the function names consistent
        - provide calls to abstract away structs as much as possible
 */

#ifndef HTSLIB_VCF_H
#define HTSLIB_VCF_H

#include <stdint.h>
#include <limits.h>
#include <errno.h>
#include "hts.h"
#include "kstring.h"
#include "hts_defs.h"
#include "hts_endian.h"

/* Included only for backwards compatibility with e.g. bcftools 1.10 */
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/*****************
 * Header struct *
 *****************/

#define BCF_HL_FLT  0 // header line
#define BCF_HL_INFO 1
#define BCF_HL_FMT  2
#define BCF_HL_CTG  3
#define BCF_HL_STR  4 // structured header line TAG=<A=..,B=..>
#define BCF_HL_GEN  5 // generic header line

#define BCF_HT_FLAG 0 // header type
#define BCF_HT_INT  1
#define BCF_HT_REAL 2
#define BCF_HT_STR  3
#define BCF_HT_LONG (BCF_HT_INT | 0x100) // BCF_HT_INT, but for int64_t values; VCF only!

#define BCF_VL_FIXED 0 // variable length
#define BCF_VL_VAR   1
#define BCF_VL_A     2
#define BCF_VL_G     3
#define BCF_VL_R     4

/* === Dictionary ===

   The header keeps three dictionaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
   in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/

#define BCF_DT_ID       0 // dictionary type
#define BCF_DT_CTG      1
#define BCF_DT_SAMPLE   2

// Complete textual representation of a header line
typedef struct bcf_hrec_t {
    int type;       // One of the BCF_HL_* type
    char *key;      // The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
    char *value;    // Set only for generic lines, NULL for FILTER/INFO, etc.
    int nkeys;              // Number of structured fields
    char **keys, **vals;    // The key=value pairs
} bcf_hrec_t;

typedef struct bcf_idinfo_t {
    uint64_t info[3];  // stores Number:20, var:4, Type:4, ColType:4 in info[0..2]
                       // for BCF_HL_FLT,INFO,FMT and contig length in info[0] for BCF_HL_CTG
    bcf_hrec_t *hrec[3];
    int id;
} bcf_idinfo_t;

typedef struct bcf_idpair_t {
    const char *key;
    const bcf_idinfo_t *val;
} bcf_idpair_t;

// Note that bcf_hdr_t structs must always be created via bcf_hdr_init()
typedef struct bcf_hdr_t {
    int32_t n[3];           // n:the size of the dictionary block in use, (allocated size, m, is below to preserve ABI)
    bcf_idpair_t *id[3];
    void *dict[3];          // ID dictionary, contig dict and sample dict
    char **samples;
    bcf_hrec_t **hrec;
    int nhrec, dirty;
    int ntransl, *transl[2];    // for bcf_translate()
    int nsamples_ori;           // for bcf_hdr_set_samples()
    uint8_t *keep_samples;
    kstring_t mem;
    int32_t m[3];          // m: allocated size of the dictionary block in use (see n above)
} bcf_hdr_t;

HTSLIB_EXPORT
extern uint8_t bcf_type_shift[];

/**************
 * VCF record *
 **************/

#define BCF_BT_NULL     0
#define BCF_BT_INT8     1
#define BCF_BT_INT16    2
#define BCF_BT_INT32    3
#define BCF_BT_INT64    4  // Unofficial, for internal use only.
#define BCF_BT_FLOAT    5
#define BCF_BT_CHAR     7

#define VCF_REF         0
#define VCF_SNP     (1<<0)
#define VCF_MNP     (1<<1)
#define VCF_INDEL   (1<<2)
#define VCF_OTHER   (1<<3)
#define VCF_BND     (1<<4)      // breakend
#define VCF_OVERLAP (1<<5)      // overlapping deletion, ALT=*
#define VCF_INS     (1<<6)      // implies VCF_INDEL
#define VCF_DEL     (1<<7)      // implies VCF_INDEL
#define VCF_ANY     (VCF_SNP|VCF_MNP|VCF_INDEL|VCF_OTHER|VCF_BND|VCF_OVERLAP|VCF_INS|VCF_DEL)       // any variant type (but not VCF_REF)

typedef struct bcf_variant_t {
    int type, n;    // variant type and the number of bases affected, negative for deletions
} bcf_variant_t;

typedef struct bcf_fmt_t {
    int id;             // id: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$id].key
    int n, size, type;  // n: number of values per-sample; size: number of bytes per-sample; type: one of BCF_BT_* types
    uint8_t *p;         // same as vptr and vptr_* in bcf_info_t below
    uint32_t p_len;
    uint32_t p_off:31, p_free:1;
} bcf_fmt_t;

typedef struct bcf_info_t {
    int key;        // key: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$key].key
    int type;  // type: one of BCF_BT_* types
    union {
        int64_t i; // integer value
        float f;   // float value
    } v1; // only set if $len==1; for easier access
    uint8_t *vptr;          // pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
    uint32_t vptr_len;      // length of the vptr block or, when set, of the vptr_mod block, excluding offset
    uint32_t vptr_off:31,   // vptr offset, i.e., the size of the INFO key plus size+type bytes
            vptr_free:1;    // indicates that vptr-vptr_off must be freed; set only when modified and the new
                            //    data block is bigger than the original
    int len;                // vector length, 1 for scalars
} bcf_info_t;


#define BCF1_DIRTY_ID  1
#define BCF1_DIRTY_ALS 2
#define BCF1_DIRTY_FLT 4
#define BCF1_DIRTY_INF 8

typedef struct bcf_dec_t {
    int m_fmt, m_info, m_id, m_als, m_allele, m_flt; // allocated size (high-water mark); do not change
    int n_flt;  // Number of FILTER fields
    int *flt;   // FILTER keys in the dictionary
    char *id, *als;     // ID and REF+ALT block (\0-separated)
    char **allele;      // allele[0] is the REF (allele[] pointers to the als block); all null terminated
    bcf_info_t *info;   // INFO
    bcf_fmt_t *fmt;     // FORMAT and individual sample
    bcf_variant_t *var; // $var and $var_type set only when set_variant_types called
    int n_var, var_type;
    int shared_dirty;   // if set, shared.s must be recreated on BCF output
    int indiv_dirty;    // if set, indiv.s must be recreated on BCF output
} bcf_dec_t;


#define BCF_ERR_CTG_UNDEF 1
#define BCF_ERR_TAG_UNDEF 2
#define BCF_ERR_NCOLS     4
#define BCF_ERR_LIMITS    8
#define BCF_ERR_CHAR     16
#define BCF_ERR_CTG_INVALID   32
#define BCF_ERR_TAG_INVALID   64

/// Get error description for bcf error code
/** @param errorcode  The error code which is to be described
    @param buffer     The buffer in which description to be added
    @param maxbuffer  The size of buffer passed
    @return NULL on invalid buffer; buffer on other cases

The buffer will be an empty string when @p errorcode is 0.
Description of errors present in code will be appended to @p buffer with ',' separation.
The buffer has to be at least 4 characters long. NULL will be returned if it is smaller or when buffer is NULL.

'...' will be appended if the description doesn't fit in the given buffer.
 */

HTSLIB_EXPORT
const char *bcf_strerror(int errorcode, char *buffer, size_t maxbuffer);

/*
    The bcf1_t structure corresponds to one VCF/BCF line. Reading from VCF file
    is slower because the string is first to be parsed, packed into BCF line
    (done in vcf_parse), then unpacked into internal bcf1_t structure. If it
    is known in advance that some of the fields will not be required (notably
    the sample columns), parsing of these can be skipped by setting max_unpack
    appropriately.
    Similarly, it is fast to output a BCF line because the columns (kept in
    shared.s, indiv.s, etc.) are written directly by bcf_write, whereas a VCF
    line must be formatted in vcf_format.
 */
typedef struct bcf1_t {
    hts_pos_t pos;  // POS
    hts_pos_t rlen; // length of REF
    int32_t rid;  // CHROM
    float qual;   // QUAL
    uint32_t n_info:16, n_allele:16;
    uint32_t n_fmt:8, n_sample:24;
    kstring_t shared, indiv;
    bcf_dec_t d; // lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
    int max_unpack;         // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields won't be needed
    int unpacked;           // remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
    int unpack_size[3];     // the original block size of ID, REF+ALT and FILTER
    int errcode;    // one of BCF_ERR_* codes
} bcf1_t;

/*******
 * API *
 *******/

    /***********************************************************************
     *  BCF and VCF I/O
     *
     *  A note about naming conventions: htslib internally represents VCF
     *  records as bcf1_t data structures, therefore most functions are
     *  prefixed with bcf_. There are a few exceptions where the functions must
     *  be aware of both BCF and VCF worlds, such as bcf_parse vs vcf_parse. In
     *  these cases, functions prefixed with bcf_ are more general and work
     *  with both BCF and VCF.
     *
     ***********************************************************************/

    /** These macros are defined only for consistency with other parts of htslib */
    #define bcf_init1()         bcf_init()
    #define bcf_read1(fp,h,v)   bcf_read((fp),(h),(v))
    #define vcf_read1(fp,h,v)   vcf_read((fp),(h),(v))
    #define bcf_write1(fp,h,v)  bcf_write((fp),(h),(v))
    #define vcf_write1(fp,h,v)  vcf_write((fp),(h),(v))
    #define bcf_destroy1(v)     bcf_destroy(v)
    #define bcf_empty1(v)       bcf_empty(v)
    #define vcf_parse1(s,h,v)   vcf_parse((s),(h),(v))
    #define bcf_clear1(v)       bcf_clear(v)
    #define vcf_format1(h,v,s)  vcf_format((h),(v),(s))

    /**
     *  bcf_hdr_init() - create an empty BCF header.
     *  @param mode    "r" or "w"
     *
     *  When opened for writing, the mandatory fileFormat and
     *  FILTER=PASS lines are added automatically.
     *
     * The bcf_hdr_t struct returned by a successful call should be freed
     * via bcf_hdr_destroy() when it is no longer needed.
     */
    HTSLIB_EXPORT
    bcf_hdr_t *bcf_hdr_init(const char *mode);

    /** Destroy a BCF header struct */
    HTSLIB_EXPORT
    void bcf_hdr_destroy(bcf_hdr_t *h);

    /** Allocate and initialize a bcf1_t object.
     *
     * The bcf1_t struct returned by a successful call should be freed
     * via bcf_destroy() when it is no longer needed.
     */
    HTSLIB_EXPORT
    bcf1_t *bcf_init(void);

    /** Deallocate a bcf1_t object */
    HTSLIB_EXPORT
    void bcf_destroy(bcf1_t *v);

    /**
     *  Same as bcf_destroy() but frees only the memory allocated by bcf1_t,
     *  not the bcf1_t object itself.
     */
    HTSLIB_EXPORT
    void bcf_empty(bcf1_t *v);

    /**
     *  Make the bcf1_t object ready for next read. Intended mostly for
     *  internal use, the user should rarely need to call this function
     *  directly.
     */
    HTSLIB_EXPORT
    void bcf_clear(bcf1_t *v);


    /** bcf_open and vcf_open mode: please see hts_open() in hts.h */
    typedef htsFile vcfFile;
    #define bcf_open(fn, mode) hts_open((fn), (mode))
    #define vcf_open(fn, mode) hts_open((fn), (mode))
    #define bcf_flush(fp) hts_flush((fp))
    #define bcf_close(fp) hts_close(fp)
    #define vcf_close(fp) hts_close(fp)

    /// Read a VCF or BCF header
    /** @param  fp  The file to read the header from
        @return Pointer to a populated header structure on success;
                NULL on failure

        The bcf_hdr_t struct returned by a successful call should be freed
        via bcf_hdr_destroy() when it is no longer needed.
    */
    HTSLIB_EXPORT
    bcf_hdr_t *bcf_hdr_read(htsFile *fp) HTS_RESULT_USED;

    /**
     *  bcf_hdr_set_samples() - for more efficient VCF parsing when only one/few samples are needed
     *  @param samples  samples to include or exclude from file or as a comma-separated string.
     *              LIST|FILE   .. select samples in list/file
     *              ^LIST|FILE  .. exclude samples from list/file
     *              -           .. include all samples
     *              NULL        .. exclude all samples
     *  @param is_file  @p samples is a file (1) or a comma-separated list (0)
     *
     *  The bottleneck of VCF reading is parsing of genotype fields. If the
     *  reader knows in advance that only subset of samples is needed (possibly
     *  no samples at all), the performance of bcf_read() can be significantly
     *  improved by calling bcf_hdr_set_samples after bcf_hdr_read().
     *  The function bcf_read() will subset the VCF/BCF records automatically
     *  with the notable exception when reading records via bcf_itr_next().
     *  In this case, bcf_subset_format() must be called explicitly, because
     *  bcf_readrec() does not see the header.
     *
     *  Returns 0 on success, -1 on error or a positive integer if the list
     *  contains samples not present in the VCF header. In such a case, the
     *  return value is the index of the offending sample.
     */
    HTSLIB_EXPORT
    int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file) HTS_RESULT_USED;

    HTSLIB_EXPORT
    int bcf_subset_format(const bcf_hdr_t *hdr, bcf1_t *rec);

    /// Write a VCF or BCF header
    /** @param  fp  Output file
        @param  h   The header to write
        @return 0 on success; -1 on failure
     */
    HTSLIB_EXPORT
    int bcf_hdr_write(htsFile *fp, bcf_hdr_t *h) HTS_RESULT_USED;

    /**
     * Parse VCF line contained in kstring and populate the bcf1_t struct
     * The line must not end with \n or \r characters.
     */
    HTSLIB_EXPORT
    int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v);

    /**
     * Complete the file opening mode, according to its extension.
     * @param mode      Preallocated mode string to be completed.
     * @param fn        File name to be opened.
     * @param format    Format string (vcf|bcf|vcf.gz)
     * @return          0 on success; -1 on failure
     */
    HTSLIB_EXPORT
    int vcf_open_mode(char *mode, const char *fn, const char *format);

    /** The opposite of vcf_parse. It should rarely be called directly, see vcf_write */
    HTSLIB_EXPORT
    int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

    /// Read next VCF or BCF record
    /** @param fp  The file to read the record from
        @param h   The header for the vcf/bcf file
        @param v   The bcf1_t structure to populate
        @return 0 on success; -1 on end of file; < -1 on critical error

On errors which are not critical for reading, such as missing header
definitions in vcf files, zero will be returned but v->errcode will have been
set to one of BCF_ERR* codes and must be checked before calling bcf_write().
     */
    HTSLIB_EXPORT
    int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v) HTS_RESULT_USED;

    /**
     *  bcf_unpack() - unpack/decode a BCF record (fills the bcf1_t::d field)
     *
     *  Note that bcf_unpack() must be called even when reading VCF. It is safe
     *  to call the function repeatedly, it will not unpack the same field
     *  twice.
     */
    #define BCF_UN_STR  1       // up to ALT inclusive
    #define BCF_UN_FLT  2       // up to FILTER
    #define BCF_UN_INFO 4       // up to INFO
    #define BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
    #define BCF_UN_FMT  8                           // unpack format and each sample
    #define BCF_UN_IND  BCF_UN_FMT                  // a synonym of BCF_UN_FMT
    #define BCF_UN_ALL  (BCF_UN_SHR|BCF_UN_FMT)     // everything
    HTSLIB_EXPORT
    int bcf_unpack(bcf1_t *b, int which);

    /*
     *  bcf_dup() - create a copy of BCF record.
     *
     *  Note that bcf_unpack() must be called on the returned copy as if it was
     *  obtained from bcf_read(). Also note that bcf_dup() calls bcf_sync1(src)
     *  internally to reflect any changes made by bcf_update_* functions.
     *
     *  The bcf1_t struct returned by a successful call should be freed
     *  via bcf_destroy() when it is no longer needed.
     */
    HTSLIB_EXPORT
    bcf1_t *bcf_dup(bcf1_t *src);

    HTSLIB_EXPORT
    bcf1_t *bcf_copy(bcf1_t *dst, bcf1_t *src);

    /// Write one VCF or BCF record. The type is determined at the open() call.
    /** @param  fp  The file to write to
        @param  h   The header for the vcf/bcf file
        @param  v   The bcf1_t structure to write
        @return 0 on success; -1 on error
     */
    HTSLIB_EXPORT
    int bcf_write(htsFile *fp, bcf_hdr_t *h, bcf1_t *v) HTS_RESULT_USED;

    /**
     *  The following functions work only with VCFs and should rarely be called
     *  directly. Usually one wants to use their bcf_* alternatives, which work
     *  transparently with both VCFs and BCFs.
     */
    /// Read a VCF format header
    /** @param  fp  The file to read the header from
        @return Pointer to a populated header structure on success;
                NULL on failure

        Use bcf_hdr_read() instead.

        The bcf_hdr_t struct returned by a successful call should be freed
        via bcf_hdr_destroy() when it is no longer needed.
    */
    HTSLIB_EXPORT
    bcf_hdr_t *vcf_hdr_read(htsFile *fp) HTS_RESULT_USED;

    /// Write a VCF format header
    /** @param  fp  Output file
        @param  h   The header to write
        @return 0 on success; -1 on failure

        Use bcf_hdr_write() instead
    */
    HTSLIB_EXPORT
    int vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h) HTS_RESULT_USED;

    /// Read a record from a VCF file
    /** @param fp  The file to read the record from
        @param h   The header for the vcf file
        @param v   The bcf1_t structure to populate
        @return 0 on success; -1 on end of file; < -1 on error

        Use bcf_read() instead
    */
    HTSLIB_EXPORT
    int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v) HTS_RESULT_USED;

    /// Write a record to a VCF file
    /** @param  fp  The file to write to
        @param h   The header for the vcf file
        @param v   The bcf1_t structure to write
        @return 0 on success; -1 on error

        Use bcf_write() instead
    */
    HTSLIB_EXPORT
    int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v) HTS_RESULT_USED;

    /** Helper function for the bcf_itr_next() macro; internal use, ignore it */
    HTSLIB_EXPORT
    int bcf_readrec(BGZF *fp, void *null, void *v, int *tid, hts_pos_t *beg, hts_pos_t *end);

    /// Write a line to a VCF file
    /** @param line   Line to write
        @param fp     File to write it to
        @return 0 on success; -1 on failure

        @note No checks are done on the line being added, apart from
              ensuring that it ends with a newline.  This function
              should therefore be used with care.
    */
    HTSLIB_EXPORT
    int vcf_write_line(htsFile *fp, kstring_t *line);

    /**************************************************************************
     *  Header querying and manipulation routines
     **************************************************************************/

    /** Create a new header using the supplied template
     *
     *  The bcf_hdr_t struct returned by a successful call should be freed
     *  via bcf_hdr_destroy() when it is no longer needed.
     *  @return NULL on failure, header otherwise
     */
    HTSLIB_EXPORT
    bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr);

    /**
     *  Copy header lines from src to dst if not already present in dst. See also bcf_translate().
     *  Returns 0 on success or sets a bit on error:
     *      1 .. conflicting definitions of tag length
     *      // todo
     */
    HTSLIB_EXPORT
    int bcf_hdr_combine(bcf_hdr_t *dst, const bcf_hdr_t *src) HTS_DEPRECATED("Please use bcf_hdr_merge instead");

    /**
     *  bcf_hdr_merge() - copy header lines from src to dst, see also bcf_translate()
     *  @param dst: the destination header to be merged into, NULL on the first pass
     *  @param src: the source header
     *  @return NULL on failure, header otherwise
     *
     *  Notes:
     *      - use as:
     *          bcf_hdr_t *dst = NULL;
     *          for (i=0; i<nsrc; i++) dst = bcf_hdr_merge(dst,src[i]);
     *
     *      - bcf_hdr_merge() replaces bcf_hdr_combine() which had a problem when
     *      combining multiple BCF headers. The current bcf_hdr_combine()
     *      does not have this problem, but became slow when used for many files.
     */
    HTSLIB_EXPORT
    bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const bcf_hdr_t *src);

    /**
     *  bcf_hdr_add_sample() - add a new sample.
     *  @param sample:  sample name to be added
     *
     *  Note:
     *      After all samples have been added, the internal header structure must be updated
     *      by calling bcf_hdr_sync(). This is normally done automatically by the first bcf_hdr_write()
     *      or bcf_write() call. Otherwise, the caller must force the update by calling bcf_hdr_sync()
     *      explicitly.
     */
    HTSLIB_EXPORT
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample);

    /** Read VCF header from a file and update the header */
    HTSLIB_EXPORT
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname);

    /// Appends formatted header text to _str_.
    /** If _is_bcf_ is zero, `IDX` fields are discarded.
     *  @return 0 if successful, or negative if an error occurred
     *  @since 1.4
     */
    HTSLIB_EXPORT
    int bcf_hdr_format(const bcf_hdr_t *hdr, int is_bcf, kstring_t *str);

    /** Returns formatted header (newly allocated string) and its length,
     *  excluding the terminating \0. If is_bcf parameter is unset, IDX
     *  fields are discarded.
     *  @deprecated Use bcf_hdr_format() instead as it can handle huge headers.
     */
    HTSLIB_EXPORT
    char *bcf_hdr_fmt_text(const bcf_hdr_t *hdr, int is_bcf, int *len)
        HTS_DEPRECATED("use bcf_hdr_format() instead");

    /** Append new VCF header line, returns 0 on success */
    HTSLIB_EXPORT
    int bcf_hdr_append(bcf_hdr_t *h, const char *line);

    HTSLIB_EXPORT
    int bcf_hdr_printf(bcf_hdr_t *h, const char *format, ...)
    HTS_FORMAT(HTS_PRINTF_FMT, 2, 3);

    /** VCF version, e.g. VCFv4.2 */
    HTSLIB_EXPORT
    const char *bcf_hdr_get_version(const bcf_hdr_t *hdr);

    /// Set version in bcf header
    /**
       @param hdr     BCF header struct
       @param version Version to set, e.g. "VCFv4.3"
       @return 0 on success; < 0 on error
     */
    HTSLIB_EXPORT
    int bcf_hdr_set_version(bcf_hdr_t *hdr, const char *version);

    /**
     *  bcf_hdr_remove() - remove VCF header tag
     *  @param type:      one of BCF_HL_*
     *  @param key:       tag name or NULL to remove all tags of the given type
     */
    HTSLIB_EXPORT
    void bcf_hdr_remove(bcf_hdr_t *h, int type, const char *key);

    /**
     *  bcf_hdr_subset() - creates a new copy of the header removing unwanted samples
     *  @param n:        number of samples to keep
     *  @param samples:  names of the samples to keep
     *  @param imap:     mapping from index in @samples to the sample index in the original file
     *  @return NULL on failure, header otherwise
     *
     *  Sample names not present in h0 are ignored. The number of unmatched samples can be checked
     *  by comparing n and bcf_hdr_nsamples(out_hdr).
     *  This function can be used to reorder samples.
     *  See also bcf_subset() which subsets individual records.
     *  The bcf_hdr_t struct returned by a successful call should be freed
     *  via bcf_hdr_destroy() when it is no longer needed.
     */
    HTSLIB_EXPORT
    bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);

    /**
     *  Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names).
     *  NB: sequence name indexes returned by bcf_hdr_seqnames() may not correspond to bcf1_t.rid, use
     *  bcf_hdr_id2name() or bcf_seqname() instead.
     */
    HTSLIB_EXPORT
    const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *nseqs);

    /** Get number of samples */
    #define bcf_hdr_nsamples(hdr) (hdr)->n[BCF_DT_SAMPLE]


    /** The following functions are for internal use and should rarely be called directly */
    HTSLIB_EXPORT
    int bcf_hdr_parse(bcf_hdr_t *hdr, char *htxt);

    /// Synchronize internal header structures
    /** @param h  Header
        @return 0 on success, -1 on failure

        This function updates the id, sample and contig arrays in the
        bcf_hdr_t structure so that they point to the same locations as
        the id, sample and contig dictionaries.
    */
    HTSLIB_EXPORT
    int bcf_hdr_sync(bcf_hdr_t *h) HTS_RESULT_USED;

    /**
     * bcf_hdr_parse_line() - parse a single line of VCF textual header
     * @param h     BCF header struct
     * @param line  One or more lines of header text
     * @param len   Filled out with length data parsed from 'line'.
     * @return bcf_hrec_t* on success;
     *         NULL on error or on end of header text.
     *         NB: to distinguish error from end-of-header, check *len:
     *           *len == 0 indicates @p line did not start with "##"
     *           *len == -1 indicates failure, likely due to out of memory
     *           *len > 0 indicates a malformed header line
     *
     * If *len > 0 on exit, it will contain the full length of the line
     * including any trailing newline (this includes cases where NULL was
     * returned due to a malformed line).  Callers can use this to skip to
     * the next header line.
     */
    HTSLIB_EXPORT
    bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len);
    /// Convert a bcf header record to string form
    /**
     * @param hrec    Header record
     * @param str     Destination kstring
     * @return 0 on success; < 0 on error
     */
    HTSLIB_EXPORT
    int bcf_hrec_format(const bcf_hrec_t *hrec, kstring_t *str);

    /// Add a header record into a header
    /**
     *  @param hdr  Destination header
     *  @param hrec Header record
     *  @return 0 on success, -1 on failure
     *
     *  If this function returns success, ownership of @p hrec will have
     *  been transferred to the header structure.  It may also have been
     *  freed if it was a duplicate of a record already in the header.
     *  Therefore the @p hrec pointer should not be used after a successful
     *  return from this function.
     *
     *  If this function returns failure, ownership will not have been taken
     *  and the caller is responsible for cleaning up @p hrec.
     */

    HTSLIB_EXPORT
    int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec);

    /**
     *  bcf_hdr_get_hrec() - get header line info
     *  @param type:  one of the BCF_HL_* types: FLT,INFO,FMT,CTG,STR,GEN
     *  @param key:   the header key for generic lines (e.g. "fileformat"), any field
     *                  for structured lines, typically "ID".
     *  @param value: the value which pairs with key. Can be be NULL for BCF_HL_GEN
     *  @param str_class: the class of BCF_HL_STR line (e.g. "ALT" or "SAMPLE"), otherwise NULL
     */
    HTSLIB_EXPORT
    bcf_hrec_t *bcf_hdr_get_hrec(const bcf_hdr_t *hdr, int type, const char *key, const char *value, const char *str_class);

    /// Duplicate a header record
    /** @param hrec   Header record to copy
        @return A new header record on success; NULL on failure

        The bcf_hrec_t struct returned by a successful call should be freed
        via bcf_hrec_destroy() when it is no longer needed.
    */
    HTSLIB_EXPORT
    bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec);

    /// Add a new header record key
    /** @param hrec  Header record
        @param str   Key name
        @param len   Length of @p str
        @return 0 on success; -1 on failure
    */
    HTSLIB_EXPORT
    int bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, size_t len) HTS_RESULT_USED;

    /// Set a header record value
    /** @param hrec      Header record
        @param i         Index of value
        @param str       Value to set
        @param len       Length of @p str
        @param is_quoted Value should be quoted
        @return 0 on success; -1 on failure
    */
    HTSLIB_EXPORT
    int bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, size_t len, int is_quoted) HTS_RESULT_USED;

    HTSLIB_EXPORT
    int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key);


    /// Add an IDX header record
    /** @param hrec   Header record
        @param idx    IDX value to add
        @return 0 on success; -1 on failure
    */
    HTSLIB_EXPORT
    int hrec_add_idx(bcf_hrec_t *hrec, int idx) HTS_RESULT_USED;

    /// Free up a header record and associated structures
    /** @param hrec  Header record
     */
    HTSLIB_EXPORT
    void bcf_hrec_destroy(bcf_hrec_t *hrec);



    /**************************************************************************
     *  Individual record querying and manipulation routines
     **************************************************************************/

    /** See the description of bcf_hdr_subset() */
    HTSLIB_EXPORT
    int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap);

    /**
     *  bcf_translate() - translate tags ids to be consistent with different header. This function
     *                    is useful when lines from multiple VCF need to be combined.
     *  @dst_hdr:   the destination header, to be used in bcf_write(), see also bcf_hdr_combine()
     *  @src_hdr:   the source header, used in bcf_read()
     *  @src_line:  line obtained by bcf_read()
     */
    HTSLIB_EXPORT
    int bcf_translate(const bcf_hdr_t *dst_hdr, bcf_hdr_t *src_hdr, bcf1_t *src_line);

    /// Get variant types in a BCF record
    /**
     *  @param rec   BCF/VCF record
     *  @return Types of variant present
     *
     *  The return value will be a bitwise-or of VCF_SNP, VCF_MNP,
     *  VCF_INDEL, VCF_OTHER, VCF_BND or VCF_OVERLAP.  If will return
     *  VCF_REF (i.e. 0) if none of the other types is present.
     *  @deprecated Please use bcf_has_variant_types() instead
     */
    HTSLIB_EXPORT
    int bcf_get_variant_types(bcf1_t *rec);

    /// Get variant type in a BCF record, for a given allele
    /**
     *  @param  rec        BCF/VCF record
     *  @param  ith_allele Allele to check
     *  @return Type of variant present
     *
     *  The return value will be one of VCF_REF, VCF_SNP, VCF_MNP,
     *  VCF_INDEL, VCF_OTHER, VCF_BND or VCF_OVERLAP.
     *  @deprecated Please use bcf_has_variant_type() instead
     */
    HTSLIB_EXPORT
    int bcf_get_variant_type(bcf1_t *rec, int ith_allele);

    /// Match mode for bcf_has_variant_types()
    enum bcf_variant_match {
        bcf_match_exact,   ///< Types present exactly match tested for
        bcf_match_overlap, ///< At least one variant type in common
        bcf_match_subset,  ///< Test set is a subset of types present
    };

    /// Check for presence of variant types in a BCF record
    /**
     *  @param rec      BCF/VCF record
     *  @param bitmask  Set of variant types to test for
     *  @param mode     Match mode
     *  @return >0 if the variant types are present,
     *           0 if not present,
     *          -1 on error
     *
     *  @p bitmask should be the bitwise-or of the variant types (VCF_SNP,
     *     VCF_MNP, etc.) to test for.
     *
     *  The return value is the bitwise-and of the set of types present
     *  and @p bitmask.  Callers that want to check for the presence of more
     *  than one type can avoid function call overhead by passing all the
     *  types to be checked for in a single call to this function, in
     *  bcf_match_overlap mode, and then check for them individually in the
     *  returned value.
     *
     *  As VCF_REF is represented by 0 (i.e. the absence of other variants)
     *  it should be tested for using
     *    bcf_has_variant_types(rec, VCF_REF, bcf_match_exact)
     *  which will return 1 if no other variant type is present, otherwise 0.
     */
    HTSLIB_EXPORT
    int bcf_has_variant_types(bcf1_t *rec, uint32_t bitmask, enum bcf_variant_match mode);

    /// Check for presence of variant types in a BCF record, for a given allele
    /**
     *  @param rec         BCF/VCF record
     *  @param ith_allele  Allele to check
     *  @param bitmask     Set of variant types to test for
     *  @return >0 if one of the variant types is present,
     *           0 if not present,
     *          -1 on error
     *
     *  @p bitmask should be the bitwise-or of the variant types (VCF_SNP,
     *     VCF_MNP, etc.) to test for, or VCF_REF on its own.
     *
     *  The return value is the bitwise-and of the set of types present
     *  and @p bitmask.  Callers that want to check for the presence of more
     *  than one type can avoid function call overhead by passing all the
     *  types to be checked for in a single call to this function, and then
     *  check for them individually in the returned value.
     *
     *  As a special case, if @p bitmask is VCF_REF (i.e. 0), the function
     *  tests for an exact match.  The return value will be 1 if the
     *  variant type calculated for the allele is VCF_REF, otherwise if
     *  any other type is present it will be 0.
     */
    HTSLIB_EXPORT
    int bcf_has_variant_type(bcf1_t *rec, int ith_allele, uint32_t bitmask);

    /// Return the number of bases affected by a variant, for a given allele
    /**
     *  @param rec         BCF/VCF record
     *  @param ith_allele  Allele index
     *  @return The number of bases affected (negative for deletions),
     *          or bcf_int32_missing on error.
     */
    HTSLIB_EXPORT
    int bcf_variant_length(bcf1_t *rec, int ith_allele);

    HTSLIB_EXPORT
    int bcf_is_snp(bcf1_t *v);

    /**
     *  bcf_update_filter() - sets the FILTER column
     *  @flt_ids:  The filter IDs to set, numeric IDs returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
     *  @n:        Number of filters. If n==0, all filters are removed
     */
    HTSLIB_EXPORT
    int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n);
    /**
     *  bcf_add_filter() - adds to the FILTER column
     *  @flt_id:   filter ID to add, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
     *
     *  If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed.
     */
    HTSLIB_EXPORT
    int bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id);
    /**
     *  bcf_remove_filter() - removes from the FILTER column
     *  @flt_id:   filter ID to remove, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
     *  @pass:     when set to 1 and no filters are present, set to PASS
     */
    HTSLIB_EXPORT
    int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int pass);
    /**
     *  Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
     */
    HTSLIB_EXPORT
    int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter);
    /**
     *  bcf_update_alleles() and bcf_update_alleles_str() - update REF and ALT column
     *  @alleles:           Array of alleles
     *  @nals:              Number of alleles
     *  @alleles_string:    Comma-separated alleles, starting with the REF allele
     */
    HTSLIB_EXPORT
    int bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals);

    HTSLIB_EXPORT
    int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string);

    /**
      *  bcf_update_id() - sets new ID string
      *  bcf_add_id() - adds to the ID string checking for duplicates
      */
    HTSLIB_EXPORT
    int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);

    HTSLIB_EXPORT
    int bcf_add_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);

    /**
     *  bcf_update_info_*() - functions for updating INFO fields
     *  @param hdr:       the BCF header
     *  @param line:      VCF line to be edited
     *  @param key:       the INFO tag to be updated
     *  @param values:    pointer to the array of values. Pass NULL to remove the tag.
     *  @param n:         number of values in the array. When set to 0, the INFO tag is removed
     *  @return 0 on success or negative value on error.
     *
     *  The @p string in bcf_update_info_flag() is optional,
     *  @p n indicates whether the flag is set or removed.
     *
     *  Note that updating an END info tag will cause line->rlen to be
     *  updated as a side-effect (removing the tag will set it to the
     *  string length of the REF allele). If line->pos is being changed as
     *  well, it is important that this is done before calling
     *  bcf_update_info_int32() to update the END tag, otherwise rlen will be
     *  set incorrectly.  If the new END value is less than or equal to
     *  line->pos, a warning will be printed and line->rlen will be set to
     *  the length of the REF allele.
     */
    #define bcf_update_info_int32(hdr,line,key,values,n)   bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf_update_info_float(hdr,line,key,values,n)   bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf_update_info_flag(hdr,line,key,string,n)    bcf_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
    #define bcf_update_info_string(hdr,line,key,string)    bcf_update_info((hdr),(line),(key),(string),1,BCF_HT_STR)
    HTSLIB_EXPORT
    int bcf_update_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);

    /// Set or update 64-bit integer INFO values
    /**
     *  @param hdr:       the BCF header
     *  @param line:      VCF line to be edited
     *  @param key:       the INFO tag to be updated
     *  @param values:    pointer to the array of values. Pass NULL to remove the tag.
     *  @param n:         number of values in the array. When set to 0, the INFO tag is removed
     *  @return 0 on success or negative value on error.
     *
     *  This function takes an int64_t values array as input.  The data
     *  actually stored will be shrunk to the minimum size that can
     *  accept all of the values.
     *
     *  INFO values outside of the range BCF_MIN_BT_INT32 to BCF_MAX_BT_INT32
     *  can only be written to VCF files.
     */
    static inline int bcf_update_info_int64(const bcf_hdr_t *hdr, bcf1_t *line,
                                            const char *key,
                                            const int64_t *values, int n)
    {
        return bcf_update_info(hdr, line, key, values, n, BCF_HT_LONG);
    }

    /*
     *  bcf_update_format_*() - functions for updating FORMAT fields
     *  @values:    pointer to the array of values, the same number of elements
     *              is expected for each sample. Missing values must be padded
     *              with bcf_*_missing or bcf_*_vector_end values.
     *  @n:         number of values in the array. If n==0, existing tag is removed.
     *
     *  The function bcf_update_format_string() is a higher-level (slower) variant of
     *  bcf_update_format_char(). The former accepts array of \0-terminated strings
     *  whereas the latter requires that the strings are collapsed into a single array
     *  of fixed-length strings. In case of strings with variable length, shorter strings
     *  can be \0-padded. Note that the collapsed strings passed to bcf_update_format_char()
     *  are not \0-terminated.
     *
     *  Returns 0 on success or negative value on error.
     */
    #define bcf_update_format_int32(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf_update_format_float(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf_update_format_char(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_STR)
    #define bcf_update_genotypes(hdr,line,gts,n) bcf_update_format((hdr),(line),"GT",(gts),(n),BCF_HT_INT)     // See bcf_gt_ macros below

    HTSLIB_EXPORT
    int bcf_update_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char **values, int n);

    HTSLIB_EXPORT
    int bcf_update_format(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);

    // Macros for setting genotypes correctly, for use with bcf_update_genotypes only; idx corresponds
    // to VCF's GT (1-based index to ALT or 0 for the reference allele) and val is the opposite, obtained
    // from bcf_get_genotypes() below.
    #define bcf_gt_phased(idx)      (((idx)+1)<<1|1)
    #define bcf_gt_unphased(idx)    (((idx)+1)<<1)
    #define bcf_gt_missing          0
    #define bcf_gt_is_missing(val)  ((val)>>1 ? 0 : 1)
    #define bcf_gt_is_phased(idx)   ((idx)&1)
    #define bcf_gt_allele(val)      (((val)>>1)-1)

    /** Conversion between alleles indexes to Number=G genotype index (assuming diploid, all 0-based) */
    #define bcf_alleles2gt(a,b) ((a)>(b)?((a)*((a)+1)/2+(b)):((b)*((b)+1)/2+(a)))
    static inline void bcf_gt2alleles(int igt, int *a, int *b)
    {
        int k = 0, dk = 1;
        while ( k<igt ) { dk++; k += dk; }
        *b = dk - 1; *a = igt - k + *b;
    }

    /**
     * bcf_get_fmt() - returns pointer to FORMAT's field data
     * @header: for access to BCF_DT_ID dictionary
     * @line:   VCF line obtained from vcf_parse1
     * @fmt:    one of GT,PL,...
     *
     * Returns bcf_fmt_t* if the call succeeded, or returns NULL when the field
     * is not available.
     */
    HTSLIB_EXPORT
    bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);

    HTSLIB_EXPORT
    bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);

    /**
     * bcf_get_*_id() - returns pointer to FORMAT/INFO field data given the header index instead of the string ID
     * @line: VCF line obtained from vcf_parse1
     * @id:  The header index for the tag, obtained from bcf_hdr_id2int()
     *
     * Returns bcf_fmt_t* / bcf_info_t*. These functions do not check if the index is valid
     * as their goal is to avoid the header lookup.
     */
    HTSLIB_EXPORT
    bcf_fmt_t *bcf_get_fmt_id(bcf1_t *line, const int id);

    HTSLIB_EXPORT
    bcf_info_t *bcf_get_info_id(bcf1_t *line, const int id);

    /**
     *  bcf_get_info_*() - get INFO values, integers or floats
     *  @param hdr:    BCF header
     *  @param line:   BCF record
     *  @param tag:    INFO tag to retrieve
     *  @param dst:    *dst is pointer to a memory location, can point to NULL
     *  @param ndst:   pointer to the size of allocated memory
     *  @return  >=0 on success
     *          -1 .. no such INFO tag defined in the header
     *          -2 .. clash between types defined in the header and encountered in the VCF record
     *          -3 .. tag is not present in the VCF record
     *          -4 .. the operation could not be completed (e.g. out of memory)
     *
     *  Returns negative value on error or the number of values (including
     *  missing values) put in *dst on success. bcf_get_info_string() returns
     *  on success the number of characters stored excluding the nul-
     *  terminating byte. bcf_get_info_flag() does not store anything in *dst
     *  but returns 1 if the flag is set or 0 if not.
     *
     *  *dst will be reallocated if it is not big enough (i.e. *ndst is too
     *  small) or NULL on entry.  The new size will be stored in *ndst.
     */
    #define bcf_get_info_int32(hdr,line,tag,dst,ndst)  bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_INT)
    #define bcf_get_info_float(hdr,line,tag,dst,ndst)  bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_REAL)
    #define bcf_get_info_string(hdr,line,tag,dst,ndst) bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_STR)
    #define bcf_get_info_flag(hdr,line,tag,dst,ndst)   bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_FLAG)

    HTSLIB_EXPORT
    int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);

    /// Put integer INFO values into an int64_t array
    /**
     *  @param hdr:    BCF header
     *  @param line:   BCF record
     *  @param tag:    INFO tag to retrieve
     *  @param dst:    *dst is pointer to a memory location, can point to NULL
     *  @param ndst:   pointer to the size of allocated memory
     *  @return  >=0 on success
     *          -1 .. no such INFO tag defined in the header
     *          -2 .. clash between types defined in the header and encountered in the VCF record
     *          -3 .. tag is not present in the VCF record
     *          -4 .. the operation could not be completed (e.g. out of memory)
     *
     *  Returns negative value on error or the number of values (including
     *  missing values) put in *dst on success.
     *
     *  *dst will be reallocated if it is not big enough (i.e. *ndst is too
     *  small) or NULL on entry.  The new size will be stored in *ndst.
     */
    static inline int bcf_get_info_int64(const bcf_hdr_t *hdr, bcf1_t *line,
                                         const char *tag, int64_t **dst,
                                         int *ndst)
    {
        return bcf_get_info_values(hdr, line, tag,
                                   (void **) dst, ndst, BCF_HT_LONG);
    }

    /**
     *  bcf_get_format_*() - same as bcf_get_info*() above
     *
     *  The function bcf_get_format_string() is a higher-level (slower) variant of bcf_get_format_char().
     *  see the description of bcf_update_format_string() and bcf_update_format_char() above.
     *  Unlike other bcf_get_format__*() functions, bcf_get_format_string() allocates two arrays:
     *  a single block of \0-terminated strings collapsed into a single array and an array of pointers
     *  to these strings. Both arrays must be cleaned by the user.
     *
     *  Returns negative value on error or the number of written values on success.
     *
     *  Use the returned number of written values for accessing valid entries of dst, as ndst is only a
     *  watermark that can be higher than the returned value, i.e. the end of dst can contain carry-over
     *  values from previous calls to bcf_get_format_*() on lines with more values per sample.
     *
     *  Example:
     *      int ndst = 0; char **dst = NULL;
     *      if ( bcf_get_format_string(hdr, line, "XX", &dst, &ndst) > 0 )
     *          for (i=0; i<bcf_hdr_nsamples(hdr); i++) printf("%s\n", dst[i]);
     *      free(dst[0]); free(dst);
     *
     *  Example:
     *      int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
     *      int32_t *gt_arr = NULL, ngt_arr = 0;
     *
     *      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
     *      if ( ngt<=0 ) return; // GT not present
     *
     *      int max_ploidy = ngt/nsmpl;
     *      for (i=0; i<nsmpl; i++)
     *      {
     *        int32_t *ptr = gt_arr + i*max_ploidy;
     *        for (j=0; j<max_ploidy; j++)
     *        {
     *           // if true, the sample has smaller ploidy
     *           if ( ptr[j]==bcf_int32_vector_end ) break;
     *
     *           // missing allele
     *           if ( bcf_gt_is_missing(ptr[j]) ) continue;
     *
     *           // the VCF 0-based allele index
     *           int allele_index = bcf_gt_allele(ptr[j]);
     *
     *           // is phased?
     *           int is_phased = bcf_gt_is_phased(ptr[j]);
     *
     *           // .. do something ..
     *         }
     *      }
     *      free(gt_arr);
     *
     */
    #define bcf_get_format_int32(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_INT)
    #define bcf_get_format_float(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_REAL)
    #define bcf_get_format_char(hdr,line,tag,dst,ndst)   bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_STR)
    #define bcf_get_genotypes(hdr,line,dst,ndst)         bcf_get_format_values(hdr,line,"GT",(void**)(dst),ndst,BCF_HT_INT)

    HTSLIB_EXPORT
    int bcf_get_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst);

    HTSLIB_EXPORT
    int bcf_get_format_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);



    /**************************************************************************
     *  Helper functions
     **************************************************************************/

    /**
     *  bcf_hdr_id2int() - Translates string into numeric ID
     *  bcf_hdr_int2id() - Translates numeric ID into string
     *  @type:     one of BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE
     *  @id:       tag name, such as: PL, DP, GT, etc.
     *
     *  Returns -1 if string is not in dictionary, otherwise numeric ID which identifies
     *  fields in BCF records.
     */
    HTSLIB_EXPORT
    int bcf_hdr_id2int(const bcf_hdr_t *hdr, int type, const char *id);
    #define bcf_hdr_int2id(hdr,type,int_id) ((hdr)->id[type][int_id].key)

    /**
     *  bcf_hdr_name2id() - Translates sequence names (chromosomes) into numeric ID
     *  bcf_hdr_id2name() - Translates numeric ID to sequence name
     */
    static inline int bcf_hdr_name2id(const bcf_hdr_t *hdr, const char *id) { return bcf_hdr_id2int(hdr, BCF_DT_CTG, id); }
    static inline const char *bcf_hdr_id2name(const bcf_hdr_t *hdr, int rid)
    {
        if ( !hdr || rid<0 || rid>=hdr->n[BCF_DT_CTG] ) return NULL;
        return hdr->id[BCF_DT_CTG][rid].key;
    }
    static inline const char *bcf_seqname(const bcf_hdr_t *hdr, const bcf1_t *rec) {
        return bcf_hdr_id2name(hdr, rec ? rec->rid : -1);
    }

    /** Return CONTIG name, or "(unknown)"

        Like bcf_seqname(), but this function will never return NULL.  If
        the contig name cannot be found (either because @p hdr was not
        supplied or rec->rid was out of range) it returns the string
        "(unknown)".
    */
    static inline const char *bcf_seqname_safe(const bcf_hdr_t *hdr, const bcf1_t *rec) {
        const char *name = bcf_seqname(hdr, rec);
        return name ? name : "(unknown)";
    }

    /**
     *  bcf_hdr_id2*() - Macros for accessing bcf_idinfo_t
     *  @type:      one of BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT
     *  @int_id:    return value of bcf_hdr_id2int, must be >=0
     *
     *  The returned values are:
     *     bcf_hdr_id2length   ..  whether the number of values is fixed or variable, one of BCF_VL_*
     *     bcf_hdr_id2number   ..  the number of values, 0xfffff for variable length fields
     *     bcf_hdr_id2type     ..  the field type, one of BCF_HT_*
     *     bcf_hdr_id2coltype  ..  the column type, one of BCF_HL_*
     *
     *  Notes: Prior to using the macros, the presence of the info should be
     *  tested with bcf_hdr_idinfo_exists().
     */
    #define bcf_hdr_id2length(hdr,type,int_id)  ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>8 & 0xf)
    #define bcf_hdr_id2number(hdr,type,int_id)  ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>12)
    #define bcf_hdr_id2type(hdr,type,int_id)    (uint32_t)((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>4 & 0xf)
    #define bcf_hdr_id2coltype(hdr,type,int_id) (uint32_t)((hdr)->id[BCF_DT_ID][int_id].val->info[type] & 0xf)
    #define bcf_hdr_idinfo_exists(hdr,type,int_id)  ((int_id)>=0 && (int_id)<(hdr)->n[BCF_DT_ID] && (hdr)->id[BCF_DT_ID][int_id].val && bcf_hdr_id2coltype((hdr),(type),(int_id))!=0xf)
    #define bcf_hdr_id2hrec(hdr,dict_type,col_type,int_id)    ((hdr)->id[(dict_type)==BCF_DT_CTG?BCF_DT_CTG:BCF_DT_ID][int_id].val->hrec[(dict_type)==BCF_DT_CTG?0:(col_type)])
    /// Convert BCF FORMAT data to string form
    /**
     * @param s    kstring to write into
     * @param n    number of items in @p data
     * @param type type of items in @p data
     * @param data BCF format data
     * @return  0 on success
     *         -1 if out of memory
     */
    HTSLIB_EXPORT
    int bcf_fmt_array(kstring_t *s, int n, int type, void *data);

    HTSLIB_EXPORT
    uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr);

    /// Encode a variable-length char array in BCF format
    /**
     * @param s    kstring to write into
     * @param l    length of input
     * @param a    input data to encode
     * @return 0 on success; < 0 on error
     */
    HTSLIB_EXPORT
    int bcf_enc_vchar(kstring_t *s, int l, const char *a);

    /// Encode a variable-length integer array in BCF format
    /**
     * @param s      kstring to write into
     * @param n      total number of items in @p a (<= 0 to encode BCF_BT_NULL)
     * @param a      input data to encode
     * @param wsize  vector length (<= 0 is equivalent to @p n)
     * @return 0 on success; < 0 on error
     * @note @p n should be an exact multiple of @p wsize
     */
    HTSLIB_EXPORT
    int bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize);

    /// Encode a variable-length float array in BCF format
    /**
     * @param s      kstring to write into
     * @param n      total number of items in @p a (<= 0 to encode BCF_BT_NULL)
     * @param a      input data to encode
     * @return 0 on success; < 0 on error
     */
    HTSLIB_EXPORT
    int bcf_enc_vfloat(kstring_t *s, int n, float *a);


    /**************************************************************************
     *  BCF index
     *
     *  Note that these functions work with BCFs only. See synced_bcf_reader.h
     *  which provides (amongst other things) an API to work transparently with
     *  both indexed BCFs and VCFs.
     **************************************************************************/

    #define bcf_itr_destroy(iter) hts_itr_destroy(iter)
    #define bcf_itr_queryi(idx, tid, beg, end) hts_itr_query((idx), (tid), (beg), (end), bcf_readrec)
    #define bcf_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bcf_hdr_name2id), (hdr), hts_itr_query, bcf_readrec)

    static inline int bcf_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r) {
        if (htsfp->is_bgzf)
            return hts_itr_next(htsfp->fp.bgzf, itr, r, 0);

        hts_log_error("Only bgzf compressed files can be used with iterators");
        errno = EINVAL;
        return -2;
    }
/// Load a BCF index
/** @param fn   BCF file name
    @return The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.
*/
    #define bcf_index_load(fn) hts_idx_load(fn, HTS_FMT_CSI)
    #define bcf_index_seqnames(idx, hdr, nptr) hts_idx_seqnames((idx),(nptr),(hts_id2name_f)(bcf_hdr_id2name),(hdr))

/// Load a BCF index from a given index file name
/**  @param fn     Input BAM/BCF/etc filename
     @param fnidx  The input index filename
     @return  The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.
*/
    HTSLIB_EXPORT
    hts_idx_t *bcf_index_load2(const char *fn, const char *fnidx);

/// Load a BCF index from a given index file name
/**  @param fn     Input BAM/BCF/etc filename
     @param fnidx  The input index filename
     @param flags  Flags to alter behaviour (see description)
     @return  The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.

     The @p flags parameter can be set to a combination of the following
     values:

        HTS_IDX_SAVE_REMOTE   Save a local copy of any remote indexes
        HTS_IDX_SILENT_FAIL   Fail silently if the index is not present

     Equivalent to hts_idx_load3(fn, fnidx, HTS_FMT_CSI, flags);
*/
    HTSLIB_EXPORT
    hts_idx_t *bcf_index_load3(const char *fn, const char *fnidx, int flags);

    /**
     *  bcf_index_build() - Generate and save an index file
     *  @fn:         Input VCF(compressed)/BCF filename
     *  @min_shift:  log2(width of the smallest bin), e.g. a value of 14
     *  imposes a 16k base lower limit on the width of index bins.
     *  Positive to generate CSI, or 0 to generate TBI. However, a small
     *  value of min_shift would create a large index, which would lead to
     *  reduced performance when using the index. A recommended value is 14.
     *  For BCF files, only the CSI index can be generated.
     *
     *  Returns 0 if successful, or negative if an error occurred.
     *
     *  List of error codes:
     *      -1 .. indexing failed
     *      -2 .. opening @fn failed
     *      -3 .. format not indexable
     *      -4 .. failed to create and/or save the index
     */
    HTSLIB_EXPORT
    int bcf_index_build(const char *fn, int min_shift);

    /**
     *  bcf_index_build2() - Generate and save an index to a specific file
     *  @fn:         Input VCF/BCF filename
     *  @fnidx:      Output filename, or NULL to add .csi/.tbi to @fn
     *  @min_shift:  Positive to generate CSI, or 0 to generate TBI
     *
     *  Returns 0 if successful, or negative if an error occurred.
     *
     *  List of error codes:
     *      -1 .. indexing failed
     *      -2 .. opening @fn failed
     *      -3 .. format not indexable
     *      -4 .. failed to create and/or save the index
     */
    HTSLIB_EXPORT
    int bcf_index_build2(const char *fn, const char *fnidx, int min_shift);

    /**
     *  bcf_index_build3() - Generate and save an index to a specific file
     *  @fn:         Input VCF/BCF filename
     *  @fnidx:      Output filename, or NULL to add .csi/.tbi to @fn
     *  @min_shift:  Positive to generate CSI, or 0 to generate TBI
     *  @n_threads:  Number of VCF/BCF decoder threads
     *
     *  Returns 0 if successful, or negative if an error occurred.
     *
     *  List of error codes:
     *      -1 .. indexing failed
     *      -2 .. opening @fn failed
     *      -3 .. format not indexable
     *      -4 .. failed to create and/or save the index
     */
     HTSLIB_EXPORT
     int bcf_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads);

     /// Initialise fp->idx for the current format type, for VCF and BCF files.
     /** @param fp        File handle for the data file being written.
         @param h         BCF header structured (needed for BAI and CSI).
         @param min_shift CSI bin size (CSI default is 14).
         @param fnidx     Filename to write index to.  This pointer must remain valid
                          until after bcf_idx_save is called.
         @return          0 on success, <0 on failure.
         @note This must be called after the header has been written, but before
               any other data.
     */
     HTSLIB_EXPORT
     int bcf_idx_init(htsFile *fp, bcf_hdr_t *h, int min_shift, const char *fnidx);

     /// Writes the index initialised with bcf_idx_init to disk.
     /** @param fp        File handle for the data file being written.
         @return          0 on success, <0 on failure.
     */
     HTSLIB_EXPORT
     int bcf_idx_save(htsFile *fp);

/*******************
 * Typed value I/O *
 *******************/

/*
    Note that in contrast with BCFv2.1 specification, HTSlib implementation
    allows missing values in vectors. For integer types, the values 0x80,
    0x8000, 0x80000000 are interpreted as missing values and 0x81, 0x8001,
    0x80000001 as end-of-vector indicators.  Similarly for floats, the value of
    0x7F800001 is interpreted as a missing value and 0x7F800002 as an
    end-of-vector indicator.
    Note that the end-of-vector byte is not part of the vector.

    This trial BCF version (v2.2) is compatible with the VCF specification and
    enables to handle correctly vectors with different ploidy in presence of
    missing values.
 */
#define bcf_int8_vector_end  (-127)         /* INT8_MIN  + 1 */
#define bcf_int16_vector_end (-32767)       /* INT16_MIN + 1 */
#define bcf_int32_vector_end (-2147483647)  /* INT32_MIN + 1 */
#define bcf_int64_vector_end (-9223372036854775807LL)  /* INT64_MIN + 1 */
#define bcf_str_vector_end   0
#define bcf_int8_missing     (-128)          /* INT8_MIN  */
#define bcf_int16_missing    (-32767-1)      /* INT16_MIN */
#define bcf_int32_missing    (-2147483647-1) /* INT32_MIN */
#define bcf_int64_missing    (-9223372036854775807LL - 1LL)  /* INT64_MIN */

// All of the above are values, which may occur multiple times in lists of
// integers or lists of floating point.  Strings in VCF don't have
// lists - a list of strings is just another (comma-separated) string.
//
// Hence bcf_str_missing is the whole string being missing rather than
// an element of a list.  Ie a string of length zero: (0<<4)|BCF_BT_CHAR.
#define bcf_str_missing      BCF_BT_CHAR

// Limits on BCF values stored in given types.  Max values are the same
// as for the underlying type.  Min values are slightly different as
// the last 8 values for each type were reserved by BCFv2.2.
#define BCF_MAX_BT_INT8  (0x7f)        /* INT8_MAX  */
#define BCF_MAX_BT_INT16 (0x7fff)      /* INT16_MAX */
#define BCF_MAX_BT_INT32 (0x7fffffff)  /* INT32_MAX */
#define BCF_MIN_BT_INT8  (-120)        /* INT8_MIN  + 8 */
#define BCF_MIN_BT_INT16 (-32760)      /* INT16_MIN + 8 */
#define BCF_MIN_BT_INT32 (-2147483640) /* INT32_MIN + 8 */

HTSLIB_EXPORT
extern uint32_t bcf_float_vector_end;
HTSLIB_EXPORT
extern uint32_t bcf_float_missing;
static inline void bcf_float_set(float *ptr, uint32_t value)
{
    union { uint32_t i; float f; } u;
    u.i = value;
    *ptr = u.f;
}
#define bcf_float_set_vector_end(x) bcf_float_set(&(x),bcf_float_vector_end)
#define bcf_float_set_missing(x)    bcf_float_set(&(x),bcf_float_missing)
static inline int bcf_float_is_missing(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==bcf_float_missing ? 1 : 0;
}
static inline int bcf_float_is_vector_end(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==bcf_float_vector_end ? 1 : 0;
}

static inline int bcf_format_gt(bcf_fmt_t *fmt, int isample, kstring_t *str)
{
    uint32_t e = 0;
    #define BRANCH(type_t, convert, missing, vector_end) { \
        uint8_t *ptr = fmt->p + isample*fmt->size; \
        int i; \
        for (i=0; i<fmt->n; i++, ptr += sizeof(type_t)) \
        { \
            type_t val = convert(ptr); \
            if ( val == vector_end ) break; \
            if ( i ) e |= kputc("/|"[val&1], str) < 0; \
            if ( !(val>>1) ) e |= kputc('.', str) < 0; \
            else e |= kputw((val>>1) - 1, str) < 0; \
        } \
        if (i == 0) e |= kputc('.', str) < 0; \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, bcf_int32_missing, bcf_int32_vector_end); break;
        case BCF_BT_NULL:  e |= kputc('.', str) < 0; break;
        default: hts_log_error("Unexpected type %d", fmt->type); return -2;
    }
    #undef BRANCH
    return e == 0 ? 0 : -1;
}

static inline int bcf_enc_size(kstring_t *s, int size, int type)
{
    // Most common case is first
    if (size < 15) {
        if (ks_resize(s, s->l + 1) < 0)
            return -1;
        uint8_t *p = (uint8_t *)s->s + s->l;
        *p++ = (size<<4) | type;
        s->l++;
        return 0;
    }

    if (ks_resize(s, s->l + 6) < 0)
        return -1;
    uint8_t *p = (uint8_t *)s->s + s->l;
    *p++ = 15<<4|type;

    if (size < 128) {
        *p++ = 1<<4|BCF_BT_INT8;
        *p++ = size;
        s->l += 3;
    } else {
        if (size < 32768) {
            *p++ = 1<<4|BCF_BT_INT16;
            i16_to_le(size, p);
            s->l += 4;
        } else {
            *p++ = 1<<4|BCF_BT_INT32;
            i32_to_le(size, p);
            s->l += 6;
        }
    }
    return 0;
}

static inline int bcf_enc_inttype(long x)
{
    if (x <= BCF_MAX_BT_INT8 && x >= BCF_MIN_BT_INT8) return BCF_BT_INT8;
    if (x <= BCF_MAX_BT_INT16 && x >= BCF_MIN_BT_INT16) return BCF_BT_INT16;
    return BCF_BT_INT32;
}

static inline int bcf_enc_int1(kstring_t *s, int32_t x)
{
    if (ks_resize(s, s->l + 5) < 0)
        return -1;
    uint8_t *p = (uint8_t *)s->s + s->l;

    if (x == bcf_int32_vector_end) {
        // An inline implementation of bcf_enc_size with size==1 and
        // memory allocation already accounted for.
        *p = (1<<4) | BCF_BT_INT8;
        p[1] = bcf_int8_vector_end;
        s->l+=2;
    } else if (x == bcf_int32_missing) {
        *p = (1<<4) | BCF_BT_INT8;
        p[1] = bcf_int8_missing;
        s->l+=2;
    } else if (x <= BCF_MAX_BT_INT8 && x >= BCF_MIN_BT_INT8) {
        *p = (1<<4) | BCF_BT_INT8;
        p[1] = x;
        s->l+=2;
    } else if (x <= BCF_MAX_BT_INT16 && x >= BCF_MIN_BT_INT16) {
        *p = (1<<4) | BCF_BT_INT16;
        i16_to_le(x, p+1);
        s->l+=3;
    } else {
        *p = (1<<4) | BCF_BT_INT32;
        i32_to_le(x, p+1);
        s->l+=5;
    }

    return 0;
}

/// Return the value of a single typed integer.
/** @param      p    Pointer to input data block.
    @param      type One of the BCF_BT_INT* type codes
    @param[out] q    Location to store an updated value for p
    @return The integer value, or zero if @p type is not valid.

If @p type is not one of BCF_BT_INT8, BCF_BT_INT16, BCF_BT_INT32 or
BCF_BT_INT64, zero will be returned and @p *q will not be updated.
Otherwise, the integer value will be returned and @p *q will be set
to the memory location immediately following the integer value.

Cautious callers can detect invalid type codes by checking that *q has
actually been updated.
*/

static inline int64_t bcf_dec_int1(const uint8_t *p, int type, uint8_t **q)
{
    if (type == BCF_BT_INT8) {
        *q = (uint8_t*)p + 1;
        return le_to_i8(p);
    } else if (type == BCF_BT_INT16) {
        *q = (uint8_t*)p + 2;
        return le_to_i16(p);
    } else if (type == BCF_BT_INT32) {
        *q = (uint8_t*)p + 4;
        return le_to_i32(p);
    } else if (type == BCF_BT_INT64) {
        *q = (uint8_t*)p + 8;
        return le_to_i64(p);
    } else { // Invalid type.
        return 0;
    }
}

/// Return the value of a single typed integer from a byte stream.
/** @param      p    Pointer to input data block.
    @param[out] q    Location to store an updated value for p
    @return The integer value, or zero if the type code was not valid.

Reads a one-byte type code from @p p, and uses it to decode an integer
value from the following bytes in @p p.

If the type is not one of BCF_BT_INT8, BCF_BT_INT16 or BCF_BT_INT32, zero
will be returned and @p *q will unchanged.  Otherwise, the integer value will
be returned and @p *q will be set to the memory location immediately following
the integer value.

Cautious callers can detect invalid type codes by checking that *q has
actually been updated.
*/
static inline int64_t bcf_dec_typed_int1(const uint8_t *p, uint8_t **q)
{
    return bcf_dec_int1(p + 1, *p&0xf, q);
}

static inline int32_t bcf_dec_size(const uint8_t *p, uint8_t **q, int *type)
{
    *type = *p & 0xf;
    if (*p>>4 != 15) {
        *q = (uint8_t*)p + 1;
        return *p>>4;
    } else return bcf_dec_typed_int1(p + 1, q);
}

#ifdef __cplusplus
}
#endif

#endif
