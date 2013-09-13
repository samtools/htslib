/*
    todo: 
        - make the function names consistent
        - provide calls to abstract away structs as much as possible
 */

#ifndef BCF_H
#define BCF_H

#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "hts.h"
#include "kstring.h"


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

#define BCF_VL_FIXED 0 // variable length
#define BCF_VL_VAR   1
#define BCF_VL_A     2
#define BCF_VL_G     3

/* === Dictionary ===

   The header keeps three dictonaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
   in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/

#define BCF_DT_ID		0 // dictionary type
#define BCF_DT_CTG		1
#define BCF_DT_SAMPLE	2

// Complete textual representation of a header line
typedef struct {
    int type;       // One of the BCF_HL_* type
    char *key;      // The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
    char *value;    // Set only for generic lines, NULL for FILTER/INFO, etc.
    int nkeys;              // Number of structured fields
    char **keys, **vals;    // The key=value pairs
} bcf_hrec_t;

typedef struct {
	uint32_t info[3];  // stores Number:20, var:4, Type:4, ColType:4 for BCF_HL_FLT,INFO,FMT
    bcf_hrec_t *hrec[3];
	int id;
} bcf_idinfo_t;

typedef struct {
	const char *key;
	const bcf_idinfo_t *val;
} bcf_idpair_t;

typedef struct {
	int32_t l_text, n[3];
	bcf_idpair_t *id[3];
	void *dict[3]; // ID dictionary, contig dict and sample dict
	char *text, **samples;
    bcf_hrec_t **hrec;
    int nhrec;
	kstring_t mem;
} bcf_hdr_t;

extern uint8_t bcf_type_shift[];

/**************
 * VCF record *
 **************/

#define BCF_BT_NULL		0
#define BCF_BT_INT8		1
#define BCF_BT_INT16	2
#define BCF_BT_INT32	3
#define BCF_BT_FLOAT	5
#define BCF_BT_CHAR		7

#define VCF_REF   0
#define VCF_SNP   1
#define VCF_MNP   2
#define VCF_INDEL 4
#define VCF_OTHER 8

typedef struct {
	int type, n;	// variant type and the number of bases affected, negative for deletions
} variant_t;

typedef struct {
	int id, n, type, size; // bcf_hdr_t::id[BCF_DT_ID][$id].key is the key in string; $n is the number of values per-sample; $size is the per-sample size in bytes; $type is one of BCF_BT_* types
	uint8_t *p;     // same as vptr and vptr_* in bcf_info_t below
    uint32_t p_len;
    uint32_t p_off:31, p_free:1;
} bcf_fmt_t;

typedef struct {
	int key, type, len; // bcf_hdr_t::id[BCF_DT_ID][$key].key is the key in string; $len: the length of the vector
	union {
		int32_t i; // integer value
		float f;   // float value
	} v1; // only set if $len==1; for easier access
	uint8_t *vptr;          // pointer to data array in bcf1_t->shared.s, excluding sized bytes
    uint32_t vptr_len;      // length of the vptr block or, when set, of the vptr_mod block, excluding offset
    uint32_t vptr_off:31,   // vptr offset, i.e., the size of the INFO key plus sized bytes
            vptr_free:1;    // indicates that vptr-vptr_off must be freed; set only when modified and the new 
                            //    data block is bigger than the original
} bcf_info_t;


#define BCF1_DIRTY_ID  1
#define BCF1_DIRTY_ALS 2
#define BCF1_DIRTY_FLT 4
#define BCF1_DIRTY_INF 8

typedef struct {
	int m_fmt, m_info, m_id, m_als, m_allele, m_flt; // allocated size (high-water mark); do not change
	int n_flt;  // Number of FILTER fields
	int *flt;   // FILTER keys in the dictionary
    char *id, *als;     // ID and REF+ALT block (\0-seperated)
	char **allele;      // allele[0] is the REF (allele[] pointers to the als block); all null terminated
	bcf_info_t *info;   // INFO
	bcf_fmt_t *fmt;     // FORMAT and individual sample
	variant_t *var;	    // $var and $var_type set only when set_variant_types called
	int n_var, var_type;
    int shared_dirty;   // if set, shared.s must be recreated on BCF output
    int indiv_dirty;    // if set, indiv.s must be recreated on BCF output
} bcf_dec_t;


#define BCF_ERR_CTG_UNDEF 1
#define BCF_ERR_TAG_UNDEF 2

/*
    The bcf1_t structure corresponds to one VCF/BCF line. Reading from VCF file
    is slower because the string is first to be parsed, packed into BCF line
    (done in vcf_parse1), then unpacked into internal bcf1_t structure. If it
    is known in advance that some of the fields will not be required (notably
    the sample columns), parsing of these can be skipped by setting max_unpack
    appropriately.
    Similarly, it is fast to output a BCF line because the columns (kept in
    shared.s, indiv.s, etc.) are written directly by bcf_write1, whereas a VCF
    line must be formatted in vcf_format1. 
 */
typedef struct {
	int32_t rid;  // CHROM
	int32_t pos;  // POS
	int32_t rlen; // length of REF
	float qual;   // QUAL
	uint32_t n_info:16, n_allele:16;
	uint32_t n_fmt:8, n_sample:24;
	kstring_t shared, indiv;
	bcf_dec_t d; // lazy evaluation: $d is not generated by bcf_read1(), but by explicitly calling bcf_unpack()
    int max_unpack;         // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse1 when some of the fields won't be needed
	int unpacked;           // remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
	uint8_t *unpack_ptr;    // position of the last unpack call
    int errcode;    // one of BCF_ERR_* codes
} bcf1_t;

/*******
 * API *
 *******/

#ifdef __cplusplus
extern "C" {
#endif

	/***************
	 *** BCF I/O ***
	 ***************/

	bcf_hdr_t *bcf_hdr_init(const char *mode);  // "r" or "w"
	int bcf_hdr_parse(bcf_hdr_t *h);

	/**
	 * Read BCF header
	 *
	 * @param fp     BGZF file pointer; file offset must be placed at the beginning
	 * 
	 * @return BCF header struct
	 */
	bcf_hdr_t *bcf_hdr_read(BGZF *fp);

	/**
	 * Write BCF header to BCF
	 *
	 * @param fp    BGZF file pointer; file offset placed at the beginning
	 * @param h     BCF header
	 */
	void bcf_hdr_write(BGZF *fp, const bcf_hdr_t *h);

	/** Destroy a BCF header struct */
	void bcf_hdr_destroy(bcf_hdr_t *h);

	/** Initialize a bcf1_t object; equivalent to calloc(1, sizeof(bcf1_t)) */
	bcf1_t *bcf_init1();
	
	/** Deallocate a bcf1_t object */
	void bcf_destroy1(bcf1_t *v);
	void bcf_clear1(bcf1_t *v);

	/**
	 * Read one BCF record
	 *
	 * @param fp     BGZF file pointer
	 * @param v      BCF record read from $fp
	 *
	 * @return  0 on success; -1 on normal file end; <-1 on error
	 */
	int bcf_read1(BGZF *fp, bcf1_t *v);

	/**
	 * Write one BCF record
	 *
	 * @param fp     BGZF file pointer
	 * @param v      BCF record to write
	 *
	 */
	int bcf_write1(BGZF *fp, const bcf1_t *v);

	/** Helper function for the bcf_itr_next() macro; ignore it */
	int bcf_readrec(BGZF *fp, void *null, bcf1_t *v, int *tid, int *beg, int *end);

	#define BCF_UN_STR  1 // up to ALT inclusive
	#define BCF_UN_FLT  2 // up to FILTER
	#define BCF_UN_INFO 4 // up to INFO
	#define BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
	#define BCF_UN_FMT  8 // unpack format and each sample
	#define BCF_UN_IND  BCF_UN_FMT // a synonymous of BCF_UN_FMT
	#define BCF_UN_ALL  (BCF_UN_SHR|BCF_UN_FMT) // everything

	/**
	 * Unpack/decode a BCF record (fill the bcf1_t::d field
	 */
	int bcf_unpack(bcf1_t *b, int which); // to unpack everything, set $which to BCF_UN_ALL

    #define bcf_nsamples(hdr) (hdr)->n[BCF_DT_SAMPLE]

    /**
      *  bcf_id2int() - Translates string into numeric ID
      *  @type:     one of BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE
      *  @id:       tag name, such as: PL, DP, GT, etc.
      *
      *  Returns -1 if string is not in dictionary, otherwise numeric ID which identifies
      *  fields in BCF records.
      */
	int bcf_id2int(const bcf_hdr_t *hdr, int type, const char *id);

    /**
      *  bcf_name2id() - Translates sequence names (chromosomes) into numeric ID
      */
	int bcf_name2id(const bcf_hdr_t *hdr, const char *id);


    /**
      *  bcf_id2*() - Macros for accessing bcf_idinfo_t 
      *  @type:      one of BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT
      *  @int_id:    return value of bcf_id2int, must be >=0
      *
      *  The returned values are:
      *     bcf_id2length   ..  whether the number of values is fixed or variable, one of BCF_VL_* 
      *     bcf_id2number   ..  the number of values, 0xfffff for variable length fields
      *     bcf_id2type     ..  the field type, one of BCF_HT_*
      *     bcf_id2coltype  ..  the column type, one of BCF_HL_*
      *
      *  Notes: Prior to using the macros, the presence of the info should be
      *  tested with bcf_idinfo_exists().
      */
    #define bcf_id2length(hdr,type,int_id)  ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>8 & 0xf)
    #define bcf_id2number(hdr,type,int_id)  ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>12)
    #define bcf_id2type(hdr,type,int_id)    ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>4 & 0xf)
    #define bcf_id2coltype(hdr,type,int_id) ((hdr)->id[BCF_DT_ID][int_id].val->info[type] & 0xf)
    #define bcf_idinfo_exists(hdr,type,int_id)  (bcf_id2coltype(hdr,type,int_id)==0xf ? 0 : 1)
    
	void bcf_fmt_array(kstring_t *s, int n, int type, void *data);
	uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr);

	void bcf_enc_vchar(kstring_t *s, int l, const char *a);
	void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize);
	void bcf_enc_vfloat(kstring_t *s, int n, float *a);

	/*****************
	 *** BCF index ***
	 *****************/

	#define bcf_itr_destroy(iter) hts_itr_destroy(iter)
	#define bcf_itr_queryi(idx, tid, beg, end) hts_itr_query((idx), (tid), (beg), (end))
	#define bcf_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bcf_name2id), (hdr))
	#define bcf_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), (hts_readrec_f)(bcf_readrec), 0)
	#define bcf_index_load(fn) hts_idx_load(fn, HTS_FMT_CSI)

	int bcf_index_build(const char *fn, int min_shift);

	/***************
	 *** VCF I/O ***
	 ***************/

	typedef htsFile vcfFile;
	#define vcf_open(fn, mode, fn_ref) hts_open((fn), (mode), (fn_ref)) // strchr(mode, 'b')!=0 for BCF; otherwise VCF
	#define vcf_close(fp) hts_close(fp)

	bcf_hdr_t *vcf_hdr_read(htsFile *fp);
	void vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h);

	int vcf_parse1(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v);
	int vcf_format1(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);
    /**
      *  vcf_read1() - read next VCF or BCF record. 
      *
      *  Returns -1 on critical errors, 0 otherwise. On errors which are not
      *  critical for reading, such as missing header definitions, v->errcode is
      *  set and must be checked before calling vcf_write1().
      */
	int vcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);
	int vcf_write1(htsFile *fp, const bcf_hdr_t *h, const bcf1_t *v);


	/****************************************
	 *** VCF header manipulation routines ***
	 ****************************************/

    bcf_hdr_t *bcf_hdr_init(const char *mode);  // "w" or "r"
    bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr);
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample);
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname);
    void bcf_hdr_fmt_text(bcf_hdr_t *hdr);
    int bcf_hdr_append(bcf_hdr_t *h, const char *line);
    bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len);
    int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec);
    bcf_hrec_t *bcf_hdr_get_hrec(bcf_hdr_t *hdr, int type, const char *id);   // type is one of BCF_HL_FLT,..,BCF_HL_CTG
    bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec);
    void bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, int len);
    void bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, int len, int is_quoted);
    int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key);
    void bcf_hrec_destroy(bcf_hrec_t *hrec);


	/*************************
	 *** VCF/BCF utilities ***
	 *************************/

	int bcf_hdr_append(bcf_hdr_t *h, const char *line);

    /**
      *  @n:        number of samples to keep
      *  @samples:  names of the samples to keep
      *  @imap:     mapping from index in @samples to the sample index in the original file
      *
      *  Sample names not present in @h0 are ignored. The number of unmatched samples can be checked
      *  by comparing @n and out_hdr->n[BCF_DT_SAMPLE].
      */
	bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);
	int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap);
	const char **bcf_seqnames(const bcf_hdr_t *h, int *nseqs);
	int bcf_is_snp(bcf1_t *v);
	void bcf_set_variant_types(bcf1_t *v);

    /**
      *  bcf1_sync() - update BCF record
      *
      *  It is necessary to call iff one or more of the bcf1_update_* functions were used
      *  and output is BCF.
      */
    int bcf1_sync(bcf1_t *line);

    /**
     *  bcf1_update_filter() - sets the FILTER column
     *  @flt_ids:           Set of filters to set, numeric IDs returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
     *  @n:                 Number of filters. If n==0, all filters are removed
     */
    int bcf1_update_filter(bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n);
    /**
     *  bcf1_add_filter() - adds to the FILTER column
     *  @flt_id:            filter ID to set, numeric IDs returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
     *  
     *  If PASS, all existing filters are moved first. If other than PASS and PASS present in $line, it
     *  will be removed.
     */
    int bcf1_add_filter(bcf_hdr_t *hdr, bcf1_t *line, int flt_id);
    /**
     *  bcf1_update_alleles() and bcf1_update_alleles_str() - update REF and ALLT column
     *  @alleles:           Array of alleles
     *  @nals:              Number of alleles
     *  @alleles_string:    Comma-separated alleles, starting with the REF allele  
     */
    int bcf1_update_alleles(bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals);
    int bcf1_update_alleles_str(bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string);
    int bcf1_update_id(bcf_hdr_t *hdr, bcf1_t *line, const char *id);

    // If n==0, existing tag is removed. Otherwise it is updated or appended. With *_flag, $string is optional.
    // With *_string, existing tag is removed when $string set to NULL.
    #define bcf1_update_info_int32(hdr,line,key,values,n)   bcf1_update_info((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf1_update_info_float(hdr,line,key,values,n)   bcf1_update_info((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf1_update_info_flag(hdr,line,key,string,n)    bcf1_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
    #define bcf1_update_info_string(hdr,line,key,string)    bcf1_update_info((hdr),(line),(key),(string),1,BCF_HT_STR)
    int bcf1_update_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key, void *values, int n, int type);

    // If n==0, existing tag is removed. Otherwise it is updated or appended.
    #define bcf1_update_format_int32(hdr,line,key,values,n) bcf1_update_format((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf1_update_format_float(hdr,line,key,values,n) bcf1_update_format((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf1_update_format_char(hdr,line,key,values,n) bcf1_update_format((hdr),(line),(key),(values),(n),BCF_HT_STR)
    #define bcf1_update_genotypes(hdr,line,gts,n) bcf1_update_format((hdr),(line),"GT",(gts),(n),BCF_HT_INT)
    int bcf1_update_format(bcf_hdr_t *hdr, bcf1_t *line, const char *key, void *values, int n, int type);

    // For use with bcf1_update_genotypes only
    #define bcf_gt_phased(id)       ((id+1)<<1|1)
    #define bcf_gt_unphased(id)     ((id+1)<<1)
    #define bcf_gt_missing          0
    #define bcf_gt_is_phased(id)    ((id)&1)
    #define bcf_gt_allele(id)       (((id)>>1)-1)

    bcf_fmt_t *bcf_get_fmt(bcf_hdr_t *hdr, bcf1_t *line, const char *key);
    bcf_info_t *bcf_get_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key);

    /**
     *  bcf_get_info*() - get INFO values, integers or floats
     *  @hdr:       BCF header
     *  @line:      BCF record
     *  @tag:       INFO tag to retrieve
     *  @dst:       *dst is pointer to a memory location, can point to NULL
     *  @ndst:      the number of elements in the dst array
     *
     *  Returns the number of written values or a negative value on error.
     */
    #define bcf_get_info_int(hdr,line,tag,dst,ndst)    bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_INT)
    #define bcf_get_info_float(hdr,line,tag,dst,ndst)  bcf_get_info_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_REAL)
    int bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);

    /**
     *  bcf_get_format*() - same as bcf_get_info*() above
     */
    #define bcf_get_format_int(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_INT)
    #define bcf_get_format_float(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_REAL)
    int bcf_get_format_values(bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
	

#ifdef __cplusplus
}
#endif

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
#define bcf_int8_vector_end  (INT8_MIN+1)
#define bcf_int16_vector_end (INT16_MIN+1)
#define bcf_int32_vector_end (INT32_MIN+1)
#define bcf_str_vector_end   0
#define bcf_int8_missing     INT8_MIN
#define bcf_int16_missing    INT16_MIN
#define bcf_int32_missing    INT32_MIN
#define bcf_str_missing      0x07
extern uint32_t bcf_float_vector_end;
extern uint32_t bcf_float_missing;
#define bcf_float_set_vector_end(x) (*(uint32_t*)(&(x)) = bcf_float_vector_end)
#define bcf_float_set_missing(x) (*(uint32_t*)(&(x)) = bcf_float_missing)
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

static inline void bcf_format_gt(bcf_fmt_t *fmt, int isample, kstring_t *str)
{
    #define BRANCH(type_t, missing, vector_end) { \
        type_t *ptr = (type_t*) (fmt->p + isample*fmt->size); \
        int i; \
        for (i=0; i<fmt->n && ptr[i]!=vector_end; i++) \
        { \
            if ( i ) kputc("/|"[ptr[i]&1], str); \
            if ( !(ptr[i]>>1) ) kputc('.', str); \
            else kputw((ptr[i]>>1) - 1, str); \
        } \
        if (i == 0) kputc('.', str); \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr,"FIXME: type %d in bcf_format_gt?\n", fmt->type); abort(); break;
    }
    #undef BRANCH
}

static inline void bcf_enc_size(kstring_t *s, int size, int type)
{
	if (size >= 15) {
		kputc(15<<4|type, s);
		if (size >= 128) {
			if (size >= 32768) {
				int32_t x = size;
				kputc(1<<4|BCF_BT_INT32, s);
				kputsn((char*)&x, 4, s);
			} else {
				int16_t x = size;
				kputc(1<<4|BCF_BT_INT16, s);
				kputsn((char*)&x, 2, s);
			}
		} else {
			kputc(1<<4|BCF_BT_INT8, s);
			kputc(size, s);
		}
	} else kputc(size<<4|type, s);
}

static inline int bcf_enc_inttype(long x)
{
	if (x <= INT8_MAX && x > bcf_int8_missing) return BCF_BT_INT8;
	if (x <= INT16_MAX && x > bcf_int16_missing) return BCF_BT_INT16;
	return BCF_BT_INT32;
}

static inline void bcf_enc_int1(kstring_t *s, int32_t x)
{
	if (x == bcf_int32_vector_end) {
		bcf_enc_size(s, 1, BCF_BT_INT8);
		kputc(bcf_int8_vector_end, s);
	} else if (x == bcf_int32_missing) {
		bcf_enc_size(s, 1, BCF_BT_INT8);
		kputc(bcf_int8_missing, s);
	} else if (x <= INT8_MAX && x > bcf_int8_missing) {
		bcf_enc_size(s, 1, BCF_BT_INT8);
		kputc(x, s);
	} else if (x <= INT16_MAX && x > bcf_int16_missing) {
		int16_t z = x;
		bcf_enc_size(s, 1, BCF_BT_INT16);
		kputsn((char*)&z, 2, s);
	} else {
		int32_t z = x;
		bcf_enc_size(s, 1, BCF_BT_INT32);
		kputsn((char*)&z, 4, s);
	}
}

static inline int32_t bcf_dec_int1(const uint8_t *p, int type, uint8_t **q)
{
	if (type == BCF_BT_INT8) {
		*q = (uint8_t*)p + 1;
		return *(int8_t*)p;
	} else if (type == BCF_BT_INT16) {
		*q = (uint8_t*)p + 2;
		return *(int16_t*)p;
	} else {
		*q = (uint8_t*)p + 4;
		return *(int32_t*)p;
	}
}

static inline int32_t bcf_dec_typed_int1(const uint8_t *p, uint8_t **q)
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

#endif
