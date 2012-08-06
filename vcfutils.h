#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include "vcf.h"

/**
 *  calc_ac() - calculate the number of REF and ALT alleles
 *  @header:  for access to BCF_DT_ID dictionary
 *  @line:    VCF line obtain from vcf_parse1
 *  @ac:      array of length line->n_allele
 *  @which:   determined if INFO/AN,AC and indv fields be used
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
 * gt_type() - determines type of genotype
 * @fmt_ptr:  the GT format field as set for example by set_fmt_ptr
 * @i:        sample
 *
 * Returns one of GT_HOM_RR, GT_HET_RA, GT_HOM_AA.
 *
 * In the current version does not recognise unknown genotypes (".")
 * or non-reference hets.
 */
#define GT_HOM_RR 0
#define GT_HET_RA 1
#define GT_HOM_AA 2
inline int gt_type(bcf_fmt_t *fmt_ptr, int i);

#endif
