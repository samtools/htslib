/*
    Time will show if this module will be merged into others 
    or perhaps removed completely.
*/
#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include "vcf.h"

/**
 * get_fmt_ptr() - returns pointer to a FORMAT field
 * @header: for access to BCF_DT_ID dictionary
 * @line:   VCF line obtained from vcf_parse1
 * @fmt:    one of GT,PL,...
 *
 * Returns bcf_fmt_t* if the call succeeded, or returns NULL when the field
 * is not available.
 */
bcf_fmt_t *get_fmt_ptr(const bcf_hdr_t *header, bcf1_t *line, char *tag);


/**
 *  trim_alleles() - remove ALT alleles unused in genotype fields
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtain from vcf_parse1
 *
 *  Returns the number of removed alleles.
 */
int trim_alleles(const bcf_hdr_t *header, bcf1_t *line);


/**
 *  trim_alleles() - remove ALT alleles according to bitmask @mask
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @mask:    alleles to remove
 */
void remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int mask);


/**
 *  calc_ac() - calculate the number of REF and ALT alleles
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtained from vcf_parse1
 *  @ac:      array of length line->n_allele
 *  @which:   determine if INFO/AN,AC and indv fields be used
 *
 *  Returns 1 if the call succeeded, or 0 if the value could not 
 *  be determined.
 *
 *  The value of @which determines if existing INFO/AC,AN can be 
 *  used (BCF_UN_INFO) and and if indv fields can be splitted 
 *  (BCF_UN_FMT). 
 */
int calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which);


/**
 * gt_type() - determines type of the genotype
 * @fmt_ptr:  the GT format field as set for example by set_fmt_ptr
 * @isample:  sample index (starting from 0)
 * @ial:      index of the non-reference allele (starting from 1)
 *
 * Returns the type of the genotype (one of GT_HOM_RR, GT_HET_RA,
 * GT_HOM_AA, GT_HET_AA, or GT_UNKN). If $ial is not NULL and the
 * genotype has one or more non-reference alleles, $ial will be set.
 * In case of GT_HET_AA, the allele which appeared first in ALT is
 * used.
 */
#define GT_HOM_RR 0 // note: the actual value of GT_* matters, used in dosage r2 calculation
#define GT_HOM_AA 1
#define GT_HET_RA 2
#define GT_HET_AA 3
#define GT_UNKN   4
inline int gt_type(bcf_fmt_t *fmt_ptr, int isample, int *ial);

#endif
