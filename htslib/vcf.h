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
    (done in vcf_parse), then unpacked into internal bcf1_t structure. If it
    is known in advance that some of the fields will not be required (notably
    the sample columns), parsing of these can be skipped by setting max_unpack
    appropriately.
    Similarly, it is fast to output a BCF line because the columns (kept in
    shared.s, indiv.s, etc.) are written directly by bcf_write, whereas a VCF
    line must be formatted in vcf_format. 
 */
typedef struct {
	int32_t rid;  // CHROM
	int32_t pos;  // POS
	int32_t rlen; // length of REF
	float qual;   // QUAL
	uint32_t n_info:16, n_allele:16;
	uint32_t n_fmt:8, n_sample:24;
	kstring_t shared, indiv;
	bcf_dec_t d; // lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
    int max_unpack;         // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields won't be needed
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
    #define vcf_parse1(s,h,v)   vcf_parse((s),(h),(v))
    #define bcf_clear1(v)       bcf_clear(v)
    #define vcf_format1(h,v,s)  vcf_format((h),(v),(s))

    /**
     *  bcf_hdr_init() - create an empty BCF header.
     *  @param mode    "r" or "w"
     *
     *  When opened for writing, the mandatory fileFormat and
     *  FILTER=PASS lines are added automatically.
     */
	bcf_hdr_t *bcf_hdr_init(const char *mode);

   	/** Destroy a BCF header struct */
	void bcf_hdr_destroy(bcf_hdr_t *h);

	/** Initialize a bcf1_t object; equivalent to calloc(1, sizeof(bcf1_t)) */
	bcf1_t *bcf_init();
	
	/** Deallocate a bcf1_t object */
	void bcf_destroy(bcf1_t *v);

    /** 
     *  Same as bcf_destroy() but frees only the memory allocated by bcf1_t,
     *  not the bcf1_t object itself. 
     */
	void bcf_empty(bcf1_t *v);

	/** 
     *  Make the bcf1_t object ready for next read. Intended mostly for
     *  internal use, the user should rarely need to call this function
     *  directly.
     */
	void bcf_clear(bcf1_t *v);


    /** bcf_open and vcf_open mode: please see hts_open() in hts.h */
	typedef htsFile vcfFile;
	#define bcf_open(fn, mode) hts_open((fn), (mode))
	#define vcf_open(fn, mode) hts_open((fn), (mode))
	#define bcf_close(fp) hts_close(fp)
	#define vcf_close(fp) hts_close(fp)

    /** Reads VCF or BCF header */
	bcf_hdr_t *bcf_hdr_read(htsFile *fp);

    /** Writes VCF or BCF header */
	int bcf_hdr_write(htsFile *fp, const bcf_hdr_t *h);

    /** Parse VCF line contained in kstring and populate the bcf1_t struct */
	int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v);

    /** The opposite of vcf_parse. It should rarely be called directly, see vcf_write */
	int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

    /**
     *  bcf_read() - read next VCF or BCF record
     *
     *  Returns -1 on critical errors, 0 otherwise. On errors which are not
     *  critical for reading, such as missing header definitions, v->errcode is
     *  set to one of BCF_ERR* code and must be checked before calling
     *  vcf_write().
     */
	int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);

	/**
	 *  bcf_unpack() - unpack/decode a BCF record (fills the bcf1_t::d field)
     *  
     *  Note that bcf_unpack() must be called even when reading VCF. It is save
     *  to call the function repeatedly, it will not unpack the same field
     *  twice.
	 */
	#define BCF_UN_STR  1       // up to ALT inclusive
	#define BCF_UN_FLT  2       // up to FILTER
	#define BCF_UN_INFO 4       // up to INFO
	#define BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
	#define BCF_UN_FMT  8                           // unpack format and each sample
	#define BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
	#define BCF_UN_ALL  (BCF_UN_SHR|BCF_UN_FMT)     // everything
	int bcf_unpack(bcf1_t *b, int which);

    /**
     *  bcf_write() - write one VCF or BCF record. The type is determined at the open() call.
     */
	int bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);

    /**
     *  The following functions work only with VCFs and should rarely be called
     *  directly. Usually one wants to use their bcf_* alternatives, which work
     *  transparently with both VCFs and BCFs.
     */
	bcf_hdr_t *vcf_hdr_read(htsFile *fp);
	int vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h);
	int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);
	int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);

	/** Helper function for the bcf_itr_next() macro; internal use, ignore it */
	int bcf_readrec(BGZF *fp, void *null, bcf1_t *v, int *tid, int *beg, int *end);



	/**************************************************************************
	 *  Header querying and manipulation routines
	 **************************************************************************/

    /** Create a new header using the supplied template */
    bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr);
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample);

    /** Read VCF header from a file and update the header */
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname);

    /** Append new VCF header line, returns 0 on success */
    int bcf_hdr_append(bcf_hdr_t *h, const char *line);
    int bcf_hdr_printf(bcf_hdr_t *h, const char *format, ...);

    /** Make the header ready for output, required if bcf_hdr_append was called */
    void bcf_hdr_fmt_text(bcf_hdr_t *hdr);

    /**
     *  bcf_hdr_subset() - creates a new copy of the header removing unwanted samples 
     *  @param n:        number of samples to keep
     *  @param samples:  names of the samples to keep
     *  @param imap:     mapping from index in @samples to the sample index in the original file
     *
     *  Sample names not present in h0 are ignored. The number of unmatched samples can be checked
     *  by comparing n and bcf_hdr_nsamples(out_hdr).
     *  This function can be used to reorder samples.
     *  See also bcf_subset() which subsets individual records.
     */
	bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);

    /** Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names) */
	const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *nseqs);

    /** Get number of samples */
    #define bcf_hdr_nsamples(hdr) (hdr)->n[BCF_DT_SAMPLE]


    /** The following functions are for internal use and should rarely be called directly */
    bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len);
    int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec);
    bcf_hrec_t *bcf_hdr_get_hrec(bcf_hdr_t *hdr, int type, const char *id);   // type is one of BCF_HL_FLT,..,BCF_HL_CTG
    bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec);
    void bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, int len);
    void bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, int len, int is_quoted);
    int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key);
    void bcf_hrec_destroy(bcf_hrec_t *hrec);



	/**************************************************************************
	 *  Individual record querying and manipulation routines
	 **************************************************************************/

    /** See the description of bcf_hdr_subset() */
	int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap);

    /**
     *  bcf_get_variant_type[s]()  - returns one of VCF_REF, VCF_SNP, etc
     */
    int bcf_get_variant_types(bcf1_t *rec);
    int bcf_get_variant_type(bcf1_t *rec, int ith_allele);
	int bcf_is_snp(bcf1_t *v);

    /**
     *  bcf_update_filter() - sets the FILTER column
     *  @flt_ids:  The filter IDs to set, numeric IDs returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
     *  @n:        Number of filters. If n==0, all filters are removed
     */
    int bcf_update_filter(bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n);
    /**
     *  bcf_add_filter() - adds to the FILTER column
     *  @flt_id:   filter ID to add, numeric ID returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
     *  
     *  If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed.
     */
    int bcf_add_filter(bcf_hdr_t *hdr, bcf1_t *line, int flt_id);
    /**
     *  bcf_update_alleles() and bcf_update_alleles_str() - update REF and ALLT column
     *  @alleles:           Array of alleles
     *  @nals:              Number of alleles
     *  @alleles_string:    Comma-separated alleles, starting with the REF allele  
     */
    int bcf_update_alleles(bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals);
    int bcf_update_alleles_str(bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string);
    int bcf_update_id(bcf_hdr_t *hdr, bcf1_t *line, const char *id);

    // If n==0, existing tag is removed. Otherwise it is updated or appended. With *_flag, $string is optional.
    // With *_string, existing tag is removed when $string set to NULL.
    #define bcf_update_info_int32(hdr,line,key,values,n)   bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf_update_info_float(hdr,line,key,values,n)   bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf_update_info_flag(hdr,line,key,string,n)    bcf_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
    #define bcf_update_info_string(hdr,line,key,string)    bcf_update_info((hdr),(line),(key),(string),1,BCF_HT_STR)
    int bcf_update_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);

    // If n==0, existing tag is removed. Otherwise it is updated or appended.
    #define bcf_update_format_int32(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_INT)
    #define bcf_update_format_float(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf_update_format_char(hdr,line,key,values,n) bcf_update_format((hdr),(line),(key),(values),(n),BCF_HT_STR)
    #define bcf_update_genotypes(hdr,line,gts,n) bcf_update_format((hdr),(line),"GT",(gts),(n),BCF_HT_INT)     // See bcf_gt_ macros below
    int bcf_update_format(bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);

    // Macros for setting genotypes correctly, for use with bcf_update_genotypes only; idx corresponds
    // to VCF's GT (1-based index to ALT or 0 for the reference allele) and val is the opposite, obtained
    // from bcf_get_genotypes() below.
    #define bcf_gt_phased(idx)      ((idx+1)<<1|1)
    #define bcf_gt_unphased(idx)    ((idx+1)<<1)
    #define bcf_gt_missing          0
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
    bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);
    bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);

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
     *  Example: 
     *      int ngt, *gt_arr = NULL, ngt_arr = 0;
     *      ngt = bcf_get_format_int(hdr, line, "GT", &gt_arr, &ngt_arr);
     */
    #define bcf_get_format_int(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_INT)
    #define bcf_get_format_float(hdr,line,tag,dst,ndst)  bcf_get_format_values(hdr,line,tag,(void**)(dst),ndst,BCF_HT_REAL)
    #define bcf_get_genotypes(hdr,line,dst,ndst)  bcf_get_format_values(hdr,line,"GT",(void**)(dst),ndst,BCF_HT_STR)
    int bcf_get_format_values(bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);


	
	/**************************************************************************
	 *  Helper functions
	 **************************************************************************/
 
    /**
     *  bcf_hdr_id2int() - Translates string into numeric ID
     *  @type:     one of BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE
     *  @id:       tag name, such as: PL, DP, GT, etc.
     *
     *  Returns -1 if string is not in dictionary, otherwise numeric ID which identifies
     *  fields in BCF records.
     */
	int bcf_hdr_id2int(const bcf_hdr_t *hdr, int type, const char *id);

    /**
     *  bcf_hdr_name2id() - Translates sequence names (chromosomes) into numeric ID
     *  bcf_hdr_id2name() - Translates numeric ID to sequence name
     */
	static inline int bcf_hdr_name2id(const bcf_hdr_t *hdr, const char *id) { return bcf_hdr_id2int(hdr, BCF_DT_CTG, id); }
    static inline const char *bcf_hdr_id2name(const bcf_hdr_t *hdr, int rid) { return hdr->id[BCF_DT_CTG][rid].key; }
    static inline const char *bcf_seqname(const bcf_hdr_t *hdr, bcf1_t *rec) { return hdr->id[BCF_DT_CTG][rec->rid].key; }

    /**
     *  bcf_hdr_id2*() - Macros for accessing bcf_idinfo_t 
     *  @type:      one of BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT
     *  @int_id:    return value of bcf_id2int, must be >=0
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
    #define bcf_hdr_id2type(hdr,type,int_id)    ((hdr)->id[BCF_DT_ID][int_id].val->info[type]>>4 & 0xf)
    #define bcf_hdr_id2coltype(hdr,type,int_id) ((hdr)->id[BCF_DT_ID][int_id].val->info[type] & 0xf)
    #define bcf_hdr_idinfo_exists(hdr,type,int_id)  ((int_id<0 || bcf_hdr_id2coltype(hdr,type,int_id)==0xf) ? 0 : 1)
    
	void bcf_fmt_array(kstring_t *s, int n, int type, void *data);
	uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr);

	void bcf_enc_vchar(kstring_t *s, int l, const char *a);
	void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize);
	void bcf_enc_vfloat(kstring_t *s, int n, float *a);

 
	/**************************************************************************
	 *  BCF index
     *
     *  Note that these functions work with BCFs only. See synced_bcf_reader.h
     *  which provides (amongst other things) an API to work transparently with
     *  both indexed BCFs and VCFs.
	 **************************************************************************/
 
	#define bcf_itr_destroy(iter) hts_itr_destroy(iter)
	#define bcf_itr_queryi(idx, tid, beg, end) hts_itr_query((idx), (tid), (beg), (end))
	#define bcf_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bcf_name2id), (hdr))
	#define bcf_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), (hts_readrec_f)(bcf_readrec), 0)
	#define bcf_index_load(fn) hts_idx_load(fn, HTS_FMT_CSI)

	int bcf_index_build(const char *fn, int min_shift);

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
