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

#endif
