#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "cram/cram.h"
#include "sam.h"

/*---------------------------------------------------------------------------
 * Samtools compatibility portion
 */
int bam_construct_seq(bam_seq_t **bp, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end,        // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
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
    int i, bam_len;

    //b->l_aux = extra_len; // we fill this out later

    bam_len = qname_len + 1 + ncigar*4 + (len+1)/2 + len + extra_len;
    if (b->m_data < bam_len) {
	b->m_data = bam_len;
	kroundup32(b->m_data);
	b->data = (uint8_t*)realloc(b->data, b->m_data);
	if (!b->data)
	    return -1;
    }
    b->l_data = bam_len;

    b->core.tid     = rname;
    b->core.pos     = pos-1;
    b->core.bin     = bam_reg2bin(pos, end);
    b->core.qual    = mapq;
    b->core.l_qname = qname_len+1;
    b->core.flag    = flag;
    b->core.n_cigar = ncigar;
    b->core.l_qseq  = len;
    b->core.mtid    = mrnm;
    b->core.mpos    = mpos-1;
    b->core.isize   = isize;

    cp = b->data;

    strncpy((char *)cp, qname, qname_len);
    cp[qname_len] = 0;
    cp += qname_len+1;
    memcpy(cp, cigar, ncigar*4);
    cp += ncigar*4;

    for (i = 0; i+1 < len; i+=2) {
	*cp++ = (L[(uc)seq[i]]<<4) + L[(uc)seq[i+1]];
    }
    if (i < len)
	*cp++ = L[(uc)seq[i]]<<4;

    memcpy(cp, qual, len);
    cp += len;

    return 0;
}

bam_hdr_t *cram_header_to_bam(SAM_hdr *h) {
    int i;
    bam_hdr_t *header = bam_hdr_init();

    header->l_text = ks_len(&h->text);
    header->text = malloc(header->l_text+1);
    memcpy(header->text, ks_str(&h->text), header->l_text);
    header->text[header->l_text] = 0;

    header->n_targets = h->nref;
    header->target_name = (char **)calloc(header->n_targets,
					  sizeof(char *));
    header->target_len = (uint32_t *)calloc(header->n_targets, 4);

    for (i = 0; i < h->nref; i++) {
	header->target_name[i] = strdup(h->ref[i].name);
	header->target_len[i] = h->ref[i].len;
    }

    return header;
}

SAM_hdr *bam_header_to_cram(bam_hdr_t *h) {
    return sam_hdr_parse_(h->text, h->l_text);
}
