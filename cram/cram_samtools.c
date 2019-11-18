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

#include "cram/cram.h"
#include "htslib/sam.h"
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
                      const char *qual) {
    static const char L[256] = {
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };
    bam1_t *b = (bam1_t *)*bp;
    uint8_t *cp;
    int i, qname_nuls, bam_len;

    //b->l_aux = extra_len; // we fill this out later

    qname_nuls = 4 - qname_len%4;
    bam_len = qname_len + qname_nuls + ncigar*4 + (len+1)/2 + len + extra_len;
    if (realloc_bam_data(b, bam_len) < 0)
        return -1;
    b->l_data = bam_len;

    b->core.tid     = rname;
    b->core.pos     = pos-1;
    b->core.bin     = bam_reg2bin(pos-1, end);
    b->core.qual    = mapq;
    b->core.l_qname = qname_len+qname_nuls;
    b->core.l_extranul = qname_nuls-1;
    b->core.flag    = flag;
    b->core.n_cigar = ncigar;
    b->core.l_qseq  = len;
    b->core.mtid    = mrnm;
    b->core.mpos    = mpos-1;
    b->core.isize   = isize;

    cp = b->data;

    strncpy((char *)cp, qname, qname_len);
    for (i = 0; i < qname_nuls; i++)
        cp[qname_len+i] = '\0';
    cp += qname_len+qname_nuls;
    if (ncigar > 0) memcpy(cp, cigar, ncigar*4);
    cp += ncigar*4;

    for (i = 0; i+1 < len; i+=2) {
        *cp++ = (L[(uc)seq[i]]<<4) + L[(uc)seq[i+1]];
    }
    if (i < len)
        *cp++ = L[(uc)seq[i]]<<4;

    if (qual)
        memcpy(cp, qual, len);
    else
        memset(cp, '\xff', len);

    return bam_len;
}
