/*
Copyright (c) 2010-2013, 2017-2019 Genome Research Ltd.
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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "cram.h"
#include "../htslib/sam.h"
#include "../sam_internal.h"

/*---------------------------------------------------------------------------
 * Samtools compatibility portion
 */
int bam_construct_seq(bam_seq_t **bp, size_t extra_len,
                      const char *qname, size_t qname_len,
                      int flag,
                      int rname,      // Ref ID
                      int64_t pos,
                      int64_t end,        // aligned start/end coords
                      int mapq,
                      uint32_t ncigar, const uint32_t *cigar,
                      int mrnm,       // Mate Ref ID
                      int64_t mpos,
                      int64_t isize,
                      int len,
                      const char *seq,
                      const char *qual)
{
    int r = bam_construct((bam1_t*)*bp,
                          qname_len, qname,
                          flag, rname, pos - 1, mapq,
                          ncigar, cigar,
                          mrnm, mpos - 1, isize,
                          len, seq, qual,
                          extra_len);
    if (r < 0) {
        return r;
    }

    // mark the buffer space allocated for aux data as containing actual aux data
    ((bam1_t*)*bp)->l_data += extra_len;
    return r + extra_len;
}
