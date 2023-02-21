/*  sam.c -- SAM and BAM file I/O and manipulation.

    Copyright (C) 2008-2010, 2012-2023 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include <signal.h>
#include <inttypes.h>
#include <unistd.h>

// Suppress deprecation message for cigar_tab, which we initialise
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hts_internal.h"
#include "sam_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/hts_expr.h"
#include "header.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_SET_INIT_INT(tag)

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif

/**********************
 *** BAM header I/O ***
 **********************/

HTSLIB_EXPORT
const int8_t bam_cigar_table[256] = {
    // 0 .. 47
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 48 .. 63  (including =)
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, BAM_CEQUAL, -1, -1,

    // 64 .. 79  (including MIDNHB)
    -1, -1, BAM_CBACK, -1,  BAM_CDEL, -1, -1, -1,
        BAM_CHARD_CLIP, BAM_CINS, -1, -1,  -1, BAM_CMATCH, BAM_CREF_SKIP, -1,

    // 80 .. 95  (including SPX)
    BAM_CPAD, -1, -1, BAM_CSOFT_CLIP,  -1, -1, -1, -1,
        BAM_CDIFF, -1, -1, -1,  -1, -1, -1, -1,

    // 96 .. 127
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,

    // 128 .. 255
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
    -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
};

sam_hdr_t *sam_hdr_init()
{
    sam_hdr_t *bh = (sam_hdr_t*)calloc(1, sizeof(sam_hdr_t));
    if (bh == NULL) return NULL;

    bh->cigar_tab = bam_cigar_table;
    return bh;
}

void sam_hdr_destroy(sam_hdr_t *bh)
{
    int32_t i;

    if (bh == NULL) return;

    if (bh->ref_count > 0) {
        --bh->ref_count;
        return;
    }

    if (bh->target_name) {
        for (i = 0; i < bh->n_targets; ++i)
            free(bh->target_name[i]);
        free(bh->target_name);
        free(bh->target_len);
    }
    free(bh->text);
    if (bh->hrecs)
        sam_hrecs_free(bh->hrecs);
    if (bh->sdict)
        kh_destroy(s2i, (khash_t(s2i) *) bh->sdict);
    free(bh);
}

// Copy the sam_hdr_t::sdict hash, used to store the real lengths of long
// references before sam_hdr_t::hrecs is populated
int sam_hdr_dup_sdict(const sam_hdr_t *h0, sam_hdr_t *h)
{
    const khash_t(s2i) *src_long_refs = (khash_t(s2i) *) h0->sdict;
    khash_t(s2i) *dest_long_refs = kh_init(s2i);
    int i;
    if (!dest_long_refs) return -1;

    for (i = 0; i < h->n_targets; i++) {
        int ret;
        khiter_t ksrc, kdest;
        if (h->target_len[i] < UINT32_MAX) continue;
        ksrc = kh_get(s2i, src_long_refs, h->target_name[i]);
        if (ksrc == kh_end(src_long_refs)) continue;
        kdest = kh_put(s2i, dest_long_refs, h->target_name[i], &ret);
        if (ret < 0) {
            kh_destroy(s2i, dest_long_refs);
            return -1;
        }
        kh_val(dest_long_refs, kdest) = kh_val(src_long_refs, ksrc);
    }

    h->sdict = dest_long_refs;
    return 0;
}

sam_hdr_t *sam_hdr_dup(const sam_hdr_t *h0)
{
    if (h0 == NULL) return NULL;
    sam_hdr_t *h;
    if ((h = sam_hdr_init()) == NULL) return NULL;
    // copy the simple data
    h->n_targets = 0;
    h->ignore_sam_err = h0->ignore_sam_err;
    h->l_text = 0;

    // Then the pointery stuff

    if (!h0->hrecs) {
        h->target_len = (uint32_t*)calloc(h0->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto fail;
        h->target_name = (char**)calloc(h0->n_targets, sizeof(char*));
        if (!h->target_name) goto fail;

        int i;
        for (i = 0; i < h0->n_targets; ++i) {
            h->target_len[i] = h0->target_len[i];
            h->target_name[i] = strdup(h0->target_name[i]);
            if (!h->target_name[i]) break;
        }
        h->n_targets = i;
        if (i < h0->n_targets) goto fail;

        if (h0->sdict) {
            if (sam_hdr_dup_sdict(h0, h) < 0) goto fail;
        }
    }

    if (h0->hrecs) {
        kstring_t tmp = { 0, 0, NULL };
        if (sam_hrecs_rebuild_text(h0->hrecs, &tmp) != 0) {
            free(ks_release(&tmp));
            goto fail;
        }

        h->l_text = tmp.l;
        h->text   = ks_release(&tmp);

        if (sam_hdr_update_target_arrays(h, h0->hrecs, 0) != 0)
            goto fail;
    } else {
        h->l_text = h0->l_text;
        h->text = malloc(h->l_text + 1);
        if (!h->text) goto fail;
        memcpy(h->text, h0->text, h->l_text);
        h->text[h->l_text] = '\0';
    }

    return h;

 fail:
    sam_hdr_destroy(h);
    return NULL;
}

sam_hdr_t *bam_hdr_read(BGZF *fp)
{
    sam_hdr_t *h;
    uint8_t buf[4];
    int magic_len, has_EOF;
    int32_t i, name_len, num_names = 0;
    size_t bufsize;
    ssize_t bytes;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) {
        perror("[W::bam_hdr_read] bgzf_check_EOF");
    } else if (has_EOF == 0) {
        hts_log_warning("EOF marker is absent. The input is probably truncated");
    }
    // read "BAM1"
    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || memcmp(buf, "BAM\1", 4)) {
        hts_log_error("Invalid BAM binary header");
        return 0;
    }
    h = sam_hdr_init();
    if (!h) goto nomem;

    // read plain text and the number of reference sequences
    bytes = bgzf_read(fp, buf, 4);
    if (bytes != 4) goto read_err;
    h->l_text = le_to_u32(buf);

    bufsize = h->l_text + 1;
    if (bufsize < h->l_text) goto nomem; // so large that adding 1 overflowed
    h->text = (char*)malloc(bufsize);
    if (!h->text) goto nomem;
    h->text[h->l_text] = 0; // make sure it is NULL terminated
    bytes = bgzf_read(fp, h->text, h->l_text);
    if (bytes != h->l_text) goto read_err;

    bytes = bgzf_read(fp, &h->n_targets, 4);
    if (bytes != 4) goto read_err;
    if (fp->is_be) ed_swap_4p(&h->n_targets);

    if (h->n_targets < 0) goto invalid;

    // read reference sequence names and lengths
    if (h->n_targets > 0) {
        h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
        if (!h->target_name) goto nomem;
        h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto nomem;
    }
    else {
        h->target_name = NULL;
        h->target_len = NULL;
    }

    for (i = 0; i != h->n_targets; ++i) {
        bytes = bgzf_read(fp, &name_len, 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&name_len);
        if (name_len <= 0) goto invalid;

        h->target_name[i] = (char*)malloc(name_len);
        if (!h->target_name[i]) goto nomem;
        num_names++;

        bytes = bgzf_read(fp, h->target_name[i], name_len);
        if (bytes != name_len) goto read_err;

        if (h->target_name[i][name_len - 1] != '\0') {
            /* Fix missing NUL-termination.  Is this being too nice?
               We could alternatively bail out with an error. */
            char *new_name;
            if (name_len == INT32_MAX) goto invalid;
            new_name = realloc(h->target_name[i], name_len + 1);
            if (new_name == NULL) goto nomem;
            h->target_name[i] = new_name;
            h->target_name[i][name_len] = '\0';
        }

        bytes = bgzf_read(fp, &h->target_len[i], 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&h->target_len[i]);
    }
    return h;

 nomem:
    hts_log_error("Out of memory");
    goto clean;

 read_err:
    if (bytes < 0) {
        hts_log_error("Error reading BGZF stream");
    } else {
        hts_log_error("Truncated BAM header");
    }
    goto clean;

 invalid:
    hts_log_error("Invalid BAM binary header");

 clean:
    if (h != NULL) {
        h->n_targets = num_names; // ensure we free only allocated target_names
        sam_hdr_destroy(h);
    }
    return NULL;
}

int bam_hdr_write(BGZF *fp, const sam_hdr_t *h)
{
    int32_t i, name_len, x;
    kstring_t hdr_ks = { 0, 0, NULL };
    char *text;
    uint32_t l_text;

    if (!h) return -1;

    if (h->hrecs) {
        if (sam_hrecs_rebuild_text(h->hrecs, &hdr_ks) != 0) return -1;
        if (hdr_ks.l > UINT32_MAX) {
            hts_log_error("Header too long for BAM format");
            free(hdr_ks.s);
            return -1;
        } else if (hdr_ks.l > INT32_MAX) {
            hts_log_warning("Header too long for BAM specification (>2GB)");
            hts_log_warning("Output file may not be portable");
        }
        text = hdr_ks.s;
        l_text = hdr_ks.l;
    } else {
        if (h->l_text > UINT32_MAX) {
            hts_log_error("Header too long for BAM format");
            return -1;
        } else if (h->l_text > INT32_MAX) {
            hts_log_warning("Header too long for BAM specification (>2GB)");
            hts_log_warning("Output file may not be portable");
        }
        text = h->text;
        l_text = h->l_text;
    }
    // write "BAM1"
    if (bgzf_write(fp, "BAM\1", 4) < 0) { free(hdr_ks.s); return -1; }
    // write plain text and the number of reference sequences
    if (fp->is_be) {
        x = ed_swap_4(l_text);
        if (bgzf_write(fp, &x, 4) < 0) { free(hdr_ks.s); return -1; }
        if (l_text) {
            if (bgzf_write(fp, text, l_text) < 0) { free(hdr_ks.s); return -1; }
        }
        x = ed_swap_4(h->n_targets);
        if (bgzf_write(fp, &x, 4) < 0) { free(hdr_ks.s); return -1; }
    } else {
        if (bgzf_write(fp, &l_text, 4) < 0) { free(hdr_ks.s); return -1; }
        if (l_text) {
            if (bgzf_write(fp, text, l_text) < 0) { free(hdr_ks.s); return -1; }
        }
        if (bgzf_write(fp, &h->n_targets, 4) < 0) { free(hdr_ks.s); return -1; }
    }
    free(hdr_ks.s);
    // write sequence names and lengths
    for (i = 0; i != h->n_targets; ++i) {
        char *p = h->target_name[i];
        name_len = strlen(p) + 1;
        if (fp->is_be) {
            x = ed_swap_4(name_len);
            if (bgzf_write(fp, &x, 4) < 0) return -1;
        } else {
            if (bgzf_write(fp, &name_len, 4) < 0) return -1;
        }
        if (bgzf_write(fp, p, name_len) < 0) return -1;
        if (fp->is_be) {
            x = ed_swap_4(h->target_len[i]);
            if (bgzf_write(fp, &x, 4) < 0) return -1;
        } else {
            if (bgzf_write(fp, &h->target_len[i], 4) < 0) return -1;
        }
    }
    if (bgzf_flush(fp) < 0) return -1;
    return 0;
}

const char *sam_parse_region(sam_hdr_t *h, const char *s, int *tid,
                             hts_pos_t *beg, hts_pos_t *end, int flags) {
    return hts_parse_region(s, tid, beg, end, (hts_name2id_f)bam_name2id, h, flags);
}

/*************************
 *** BAM alignment I/O ***
 *************************/

bam1_t *bam_init1()
{
    return (bam1_t*)calloc(1, sizeof(bam1_t));
}

int sam_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

void bam_destroy1(bam1_t *b)
{
    if (b == 0) return;
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        free(b->data);
        if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) != 0) {
            // In case of reuse
            b->data = NULL;
            b->m_data = 0;
            b->l_data = 0;
        }
    }

    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_STRUCT) == 0)
        free(b);
}

bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
    if (realloc_bam_data(bdst, bsrc->l_data) < 0) return NULL;
    memcpy(bdst->data, bsrc->data, bsrc->l_data); // copy var-len data
    memcpy(&bdst->core, &bsrc->core, sizeof(bsrc->core)); // copy the rest
    bdst->l_data = bsrc->l_data;
    bdst->id = bsrc->id;
    return bdst;
}

bam1_t *bam_dup1(const bam1_t *bsrc)
{
    if (bsrc == NULL) return NULL;
    bam1_t *bdst = bam_init1();
    if (bdst == NULL) return NULL;
    if (bam_copy1(bdst, bsrc) == NULL) {
        bam_destroy1(bdst);
        return NULL;
    }
    return bdst;
}

static void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                             hts_pos_t *rlen, hts_pos_t *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;
    }
}

static int subtract_check_underflow(size_t length, size_t *limit)
{
    if (length <= *limit) {
        *limit -= length;
        return 0;
    }

    return -1;
}

int bam_set1(bam1_t *bam,
             size_t l_qname, const char *qname,
             uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
             size_t n_cigar, const uint32_t *cigar,
             int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
             size_t l_seq, const char *seq, const char *qual,
             size_t l_aux)
{
    // use a default qname "*" if none is provided
    if (l_qname == 0) {
        l_qname = 1;
        qname = "*";
    }

    // note: the qname is stored nul terminated and padded as described in the
    // documentation for the bam1_t struct.
    size_t qname_nuls = 4 - l_qname % 4;

    // the aligment length, needed for bam_reg2bin(), is calculated as in bam_endpos().
    // can't use bam_endpos() directly as some fields not yet set up.
    hts_pos_t rlen = 0, qlen = 0;
    if (!(flag & BAM_FUNMAP)) {
        bam_cigar2rqlens((int)n_cigar, cigar, &rlen, &qlen);
    }
    if (rlen == 0) {
        rlen = 1;
    }

    // validate parameters
    if (l_qname > 254) {
        hts_log_error("Query name too long");
        errno = EINVAL;
        return -1;
    }
    if (HTS_POS_MAX - rlen <= pos) {
        hts_log_error("Read ends beyond highest supported position");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && n_cigar == 0) {
        hts_log_error("Mapped query must have a CIGAR");
        errno = EINVAL;
        return -1;
    }
    if (!(flag & BAM_FUNMAP) && l_seq > 0 && l_seq != qlen) {
        hts_log_error("CIGAR and query sequence are of different length");
        errno = EINVAL;
        return -1;
    }

    size_t limit = INT32_MAX;
    int u = subtract_check_underflow(l_qname + qname_nuls, &limit);
    u    += subtract_check_underflow(n_cigar * 4, &limit);
    u    += subtract_check_underflow((l_seq + 1) / 2, &limit);
    u    += subtract_check_underflow(l_seq, &limit);
    u    += subtract_check_underflow(l_aux, &limit);
    if (u != 0) {
        hts_log_error("Size overflow");
        errno = EINVAL;
        return -1;
    }

    // re-allocate the data buffer as needed.
    size_t data_len = l_qname + qname_nuls + n_cigar * 4 + (l_seq + 1) / 2 + l_seq;
    if (realloc_bam_data(bam, data_len + l_aux) < 0) {
        return -1;
    }

    bam->l_data = (int)data_len;
    bam->core.pos = pos;
    bam->core.tid = tid;
    bam->core.bin = bam_reg2bin(pos, pos + rlen);
    bam->core.qual = mapq;
    bam->core.l_extranul = (uint8_t)(qname_nuls - 1);
    bam->core.flag = flag;
    bam->core.l_qname = (uint16_t)(l_qname + qname_nuls);
    bam->core.n_cigar = (uint32_t)n_cigar;
    bam->core.l_qseq = (int32_t)l_seq;
    bam->core.mtid = mtid;
    bam->core.mpos = mpos;
    bam->core.isize = isize;

    uint8_t *cp = bam->data;
    strncpy((char *)cp, qname, l_qname);
    int i;
    for (i = 0; i < qname_nuls; i++) {
        cp[l_qname + i] = '\0';
    }
    cp += l_qname + qname_nuls;

    if (n_cigar > 0) {
        memcpy(cp, cigar, n_cigar * 4);
    }
    cp += n_cigar * 4;

    for (i = 0; i + 1 < l_seq; i += 2) {
        *cp++ = (seq_nt16_table[(unsigned char)seq[i]] << 4) | seq_nt16_table[(unsigned char)seq[i + 1]];
    }
    for (; i < l_seq; i++) {
        *cp++ = seq_nt16_table[(unsigned char)seq[i]] << 4;
    }

    if (qual) {
        memcpy(cp, qual, l_seq);
    }
    else {
        memset(cp, '\xff', l_seq);
    }

    return (int)data_len;
}

hts_pos_t bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}

hts_pos_t bam_endpos(const bam1_t *b)
{
    hts_pos_t rlen = (b->core.flag & BAM_FUNMAP)? 0 : bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    if (rlen == 0) rlen = 1;
    return b->core.pos + rlen;
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data, *cigar0, CG_len, fake_bytes;
    uint8_t *CG;

    // test where there is a real CIGAR in the CG tag to move
    if (c->n_cigar == 0 || c->tid < 0 || c->pos < 0) return 0;
    cigar0 = bam_get_cigar(b);
    if (bam_cigar_op(cigar0[0]) != BAM_CSOFT_CLIP || bam_cigar_oplen(cigar0[0]) != c->l_qseq) return 0;
    fake_bytes = c->n_cigar * 4;
    int saved_errno = errno;
    CG = bam_aux_get(b, "CG");
    if (!CG) {
        if (errno != ENOENT) return -1;  // Bad aux data
        errno = saved_errno; // restore errno on expected no-CG-tag case
        return 0;
    }
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0; // don't move if the real CIGAR length is shorter than the fake cigar length

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    b->l_data = b->l_data - fake_bytes + n_cigar4; // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes)); // insert c->n_cigar-fake_bytes empty space to make room
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4); // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_error("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

// Fix bad records where qname is not terminated correctly.
static int fixup_missing_qname_nul(bam1_t *b) {
    bam1_core_t *c = &b->core;

    // Note this is called before c->l_extranul is added to c->l_qname
    if (c->l_extranul > 0) {
        b->data[c->l_qname++] = '\0';
        c->l_extranul--;
    } else {
        if (b->l_data > INT_MAX - 4) return -1;
        if (realloc_bam_data(b, b->l_data + 4) < 0) return -1;
        b->l_data += 4;
        b->data[c->l_qname++] = '\0';
        c->l_extranul = 3;
    }
    return 0;
}

/*
 * Note a second interface that returns a bam pointer instead would avoid bam_copy1
 * in multi-threaded handling.  This may be worth considering for htslib2.
 */
int bam_read1(BGZF *fp, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t x[8], new_l_data;

    b->l_data = 0;

    if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
    if (fp->is_be)
        ed_swap_4p(&block_len);
    if (block_len < 32) return -4;  // block_len includes core data
    if (bgzf_read(fp, x, 32) != 32) return -3;
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
    }
    c->tid = x[0]; c->pos = (int32_t)x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = (int32_t)x[6]; c->isize = (int32_t)x[7];

    new_l_data = block_len - 32 + c->l_extranul;
    if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
        return -4;
    if (realloc_bam_data(b, new_l_data) < 0) return -4;
    b->l_data = new_l_data;

    if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
    if (b->data[c->l_qname - 1] != '\0') { // Try to fix missing NUL termination
        if (fixup_missing_qname_nul(b) < 0) return -4;
    }
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;
    if (b->l_data < c->l_qname ||
        bgzf_read(fp, b->data + c->l_qname, b->l_data - c->l_qname) != b->l_data - c->l_qname)
        return -4;
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    if (bam_tag2cigar(b, 0, 0) < 0)
        return -4;

    if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
        hts_pos_t rlen, qlen;
        bam_cigar2rqlens(c->n_cigar, bam_get_cigar(b), &rlen, &qlen);
        if ((b->core.flag & BAM_FUNMAP) || rlen == 0) rlen = 1;
        b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                    bam_get_qname(b));
            return -4;
        }
    }

    return 4 + block_len;
}

int bam_write1(BGZF *fp, const bam1_t *b)
{
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (bgzf_flush_try(fp, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (bgzf_write(fp, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) ok = (bgzf_write(fp, &block_len, 4) >= 0);
    }
    if (ok) ok = (bgzf_write(fp, x, 32) >= 0);
    if (ok) ok = (bgzf_write(fp, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (bgzf_write(fp, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1<<28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"PRIhts_pos
                          " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                          bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t*)bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t)cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) ok = (bgzf_write(fp, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        if (ok) ok = (bgzf_write(fp, &b->data[cigar_en], b->l_data - cigar_en) >= 0); // write data after CIGAR
        if (ok) ok = (bgzf_write(fp, "CGBI", 4) >= 0); // write CG:B,I
        u32_to_le(c->n_cigar, buf);
        if (ok) ok = (bgzf_write(fp, buf, 4) >= 0); // write the true CIGAR length
        if (ok) ok = (bgzf_write(fp, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
}

/*
 * Write a BAM file and append to the in-memory index simultaneously.
 */
static int bam_write_idx1(htsFile *fp, const sam_hdr_t *h, const bam1_t *b) {
    BGZF *bfp = fp->fp.bgzf;

    if (!fp->idx)
        return bam_write1(bfp, b);

    uint32_t block_len = b->l_data - b->core.l_extranul + 32;
    if (bgzf_flush_try(bfp, 4 + block_len) < 0)
        return -1;
    if (!bfp->mt)
        hts_idx_amend_last(fp->idx, bgzf_tell(bfp));
    else
        bgzf_idx_amend_last(bfp, fp->idx, bgzf_tell(bfp));

    int ret = bam_write1(bfp, b);
    if (ret < 0)
        return -1;

    if (bgzf_idx_push(bfp, fp->idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(bfp), !(b->core.flag&BAM_FUNMAP)) < 0) {
        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
        ret = -1;
    }

    return ret;
}

/*
 * Set the qname in a BAM record
 */
int bam_set_qname(bam1_t *rec, const char *qname)
{
    if (!rec) return -1;
    if (!qname || !*qname) return -1;

    size_t old_len = rec->core.l_qname;
    size_t new_len = strlen(qname) + 1;
    if (new_len < 1 || new_len > 255) return -1;

    int extranul = (new_len%4 != 0) ? (4 - new_len%4) : 0;

    size_t new_data_len = rec->l_data - old_len + new_len + extranul;
    if (realloc_bam_data(rec, new_data_len) < 0) return -1;

    // Make room
    if (new_len + extranul != rec->core.l_qname)
        memmove(rec->data + new_len + extranul, rec->data + rec->core.l_qname, rec->l_data - rec->core.l_qname);
    // Copy in new name and pad if needed
    memcpy(rec->data, qname, new_len);
    int n;
    for (n = 0; n < extranul; n++) rec->data[new_len + n] = '\0';

    rec->l_data = new_data_len;
    rec->core.l_qname = new_len + extranul;
    rec->core.l_extranul = extranul;

    return 0;
}

/********************
 *** BAM indexing ***
 ********************/

static hts_idx_t *sam_index(htsFile *fp, int min_shift)
{
    int n_lvls, i, fmt, ret;
    bam1_t *b;
    hts_idx_t *idx;
    sam_hdr_t *h;
    h = sam_hdr_read(fp);
    if (h == NULL) return NULL;
    if (min_shift > 0) {
        hts_pos_t max_len = 0, s;
        for (i = 0; i < h->n_targets; ++i) {
            hts_pos_t len = sam_hdr_tid2len(h, i);
            if (max_len < len) max_len = len;
        }
        max_len += 256;
        for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);
        fmt = HTS_FMT_CSI;
    } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;
    idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
    b = bam_init1();
    while ((ret = sam_read1(fp, h, b)) >= 0) {
        ret = hts_idx_push(idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP));
        if (ret < 0) { // unsorted or doesn't fit
            hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed", bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
            goto err;
        }
    }
    if (ret < -1) goto err; // corrupted BAM file

    hts_idx_finish(idx, bgzf_tell(fp->fp.bgzf));
    sam_hdr_destroy(h);
    bam_destroy1(b);
    return idx;

err:
    bam_destroy1(b);
    hts_idx_destroy(idx);
    return NULL;
}

int sam_index_build3(const char *fn, const char *fnidx, int min_shift, int nthreads)
{
    hts_idx_t *idx;
    htsFile *fp;
    int ret = 0;

    if ((fp = hts_open(fn, "r")) == 0) return -2;
    if (nthreads)
        hts_set_threads(fp, nthreads);

    switch (fp->format.format) {
    case cram:

        ret = cram_index_build(fp->fp.cram, fn, fnidx);
        break;

    case bam:
    case sam:
        if (fp->format.compression != bgzf) {
            hts_log_error("%s file \"%s\" not BGZF compressed",
                          fp->format.format == bam ? "BAM" : "SAM", fn);
            ret = -1;
            break;
        }
        idx = sam_index(fp, min_shift);
        if (idx) {
            ret = hts_idx_save_as(idx, fn, fnidx, (min_shift > 0)? HTS_FMT_CSI : HTS_FMT_BAI);
            if (ret < 0) ret = -4;
            hts_idx_destroy(idx);
        }
        else ret = -1;
        break;

    default:
        ret = -3;
        break;
    }
    hts_close(fp);

    return ret;
}

int sam_index_build2(const char *fn, const char *fnidx, int min_shift)
{
    return sam_index_build3(fn, fnidx, min_shift, 0);
}

int sam_index_build(const char *fn, int min_shift)
{
    return sam_index_build3(fn, NULL, min_shift, 0);
}

// Provide bam_index_build() symbol for binary compatibility with earlier HTSlib
#undef bam_index_build
int bam_index_build(const char *fn, int min_shift)
{
    return sam_index_build2(fn, NULL, min_shift);
}

// Initialise fp->idx for the current format type.
// This must be called after the header has been written but no other data.
int sam_idx_init(htsFile *fp, sam_hdr_t *h, int min_shift, const char *fnidx) {
    fp->fnidx = fnidx;
    if (fp->format.format == bam || fp->format.format == bcf ||
        (fp->format.format == sam && fp->format.compression == bgzf)) {
        int n_lvls, fmt = HTS_FMT_CSI;
        if (min_shift > 0) {
            int64_t max_len = 0, s;
            int i;
            for (i = 0; i < h->n_targets; ++i)
                if (max_len < h->target_len[i]) max_len = h->target_len[i];
            max_len += 256;
            for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);

        } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;

        fp->idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
        return fp->idx ? 0 : -1;
    }

    if (fp->format.format == cram) {
        fp->fp.cram->idxfp = bgzf_open(fnidx, "wg");
        return fp->fp.cram->idxfp ? 0 : -1;
    }

    return -1;
}

// Finishes an index. Call after the last record has been written.
// Returns 0 on success, <0 on failure.
int sam_idx_save(htsFile *fp) {
    if (fp->format.format == bam || fp->format.format == bcf ||
        fp->format.format == vcf || fp->format.format == sam) {
        int ret;
        if ((ret = sam_state_destroy(fp)) < 0) {
            errno = -ret;
            return -1;
        }
        if (bgzf_flush(fp->fp.bgzf) < 0)
            return -1;
        hts_idx_amend_last(fp->idx, bgzf_tell(fp->fp.bgzf));

        if (hts_idx_finish(fp->idx, bgzf_tell(fp->fp.bgzf)) < 0)
            return -1;

        return hts_idx_save_as(fp->idx, NULL, fp->fnidx, hts_idx_fmt(fp->idx));

    } else if (fp->format.format == cram) {
        // flushed and closed by cram_close
    }

    return 0;
}

static int sam_readrec(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = (htsFile *)fpv;
    bam1_t *b = bv;
    fp->line.l = 0;
    int ret = sam_read1(fp, fp->bam_header, b);
    if (ret >= 0) {
        *tid = b->core.tid;
        *beg = b->core.pos;
        *end = bam_endpos(b);
    }
    return ret;
}

// This is used only with read_rest=1 iterators, so need not set tid/beg/end.
static int sam_readrec_rest(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = (htsFile *)fpv;
    bam1_t *b = bv;
    fp->line.l = 0;
    int ret = sam_read1(fp, fp->bam_header, b);
    return ret;
}

// Internal (for now) func used by bam_sym_lookup.  This is copied from
// samtools/bam.c.
static const char *bam_get_library(const bam_hdr_t *h, const bam1_t *b)
{
    const char *rg;
    kstring_t lib = { 0, 0, NULL };
    rg = (char *)bam_aux_get(b, "RG");

    if (!rg)
        return NULL;
    else
        rg++;

    if (sam_hdr_find_tag_id((bam_hdr_t *)h, "RG", "ID", rg, "LB", &lib)  < 0)
        return NULL;

    static char LB_text[1024];
    int len = lib.l < sizeof(LB_text) - 1 ? lib.l : sizeof(LB_text) - 1;

    memcpy(LB_text, lib.s, len);
    LB_text[len] = 0;

    free(lib.s);

    return LB_text;
}


// Bam record pointer and SAM header combined
typedef struct {
    const sam_hdr_t *h;
    const bam1_t *b;
} hb_pair;

// Looks up variable names in str and replaces them with their value.
// Also supports aux tags.
//
// Note the expression parser deliberately overallocates str size so it
// is safe to use memcmp over strcmp.
static int bam_sym_lookup(void *data, char *str, char **end,
                          hts_expr_val_t *res) {
    hb_pair *hb = (hb_pair *)data;
    const bam1_t *b = hb->b;

    res->is_str = 0;
    switch(*str) {
    case 'c':
        if (memcmp(str, "cigar", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            ks_clear(&res->s);
            uint32_t *cigar = bam_get_cigar(b);
            int i, n = b->core.n_cigar, r = 0;
            if (n) {
                for (i = 0; i < n; i++) {
                    r |= kputw (bam_cigar_oplen(cigar[i]), &res->s) < 0;
                    r |= kputc_(bam_cigar_opchr(cigar[i]), &res->s) < 0;
                }
                r |= kputs("", &res->s) < 0;
            } else {
                r |= kputs("*", &res->s) < 0;
            }
            return r ? -1 : 0;
        }
        break;

    case 'e':
        if (memcmp(str, "endpos", 6) == 0) {
            *end = str+6;
            res->d = bam_endpos(b);
            return 0;
        }
        break;

    case 'f':
        if (memcmp(str, "flag", 4) == 0) {
            str = *end = str+4;
            if (*str != '.') {
                res->d = b->core.flag;
                return 0;
            } else {
                str++;
                if (!memcmp(str, "paired", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FPAIRED;
                    return 0;
                } else if (!memcmp(str, "proper_pair", 11)) {
                    *end = str+11;
                    res->d = b->core.flag & BAM_FPROPER_PAIR;
                    return 0;
                } else if (!memcmp(str, "unmap", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FUNMAP;
                    return 0;
                } else if (!memcmp(str, "munmap", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FMUNMAP;
                    return 0;
                } else if (!memcmp(str, "reverse", 7)) {
                    *end = str+7;
                    res->d = b->core.flag & BAM_FREVERSE;
                    return 0;
                } else if (!memcmp(str, "mreverse", 8)) {
                    *end = str+8;
                    res->d = b->core.flag & BAM_FMREVERSE;
                    return 0;
                } else if (!memcmp(str, "read1", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FREAD1;
                    return 0;
                } else if (!memcmp(str, "read2", 5)) {
                    *end = str+5;
                    res->d = b->core.flag & BAM_FREAD2;
                    return 0;
                } else if (!memcmp(str, "secondary", 9)) {
                    *end = str+9;
                    res->d = b->core.flag & BAM_FSECONDARY;
                    return 0;
                } else if (!memcmp(str, "qcfail", 6)) {
                    *end = str+6;
                    res->d = b->core.flag & BAM_FQCFAIL;
                    return 0;
                } else if (!memcmp(str, "dup", 3)) {
                    *end = str+3;
                    res->d = b->core.flag & BAM_FDUP;
                    return 0;
                } else if (!memcmp(str, "supplementary", 13)) {
                    *end = str+13;
                    res->d = b->core.flag & BAM_FSUPPLEMENTARY;
                    return 0;
                } else {
                    hts_log_error("Unrecognised flag string");
                    return -1;
                }
            }
        }
        break;

    case 'l':
        if (memcmp(str, "library", 7) == 0) {
            *end = str+7;
            res->is_str = 1;
            const char *lib = bam_get_library(hb->h, b);
            kputs(lib ? lib : "", ks_clear(&res->s));
            return 0;
        }
        break;

    case 'm':
        if (memcmp(str, "mapq", 4) == 0) {
            *end = str+4;
            res->d = b->core.qual;
            return 0;
        } else if (memcmp(str, "mpos", 4) == 0) {
            *end = str+4;
            res->d = b->core.mpos+1;
            return 0;
        } else if (memcmp(str, "mrname", 6) == 0) {
            *end = str+6;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.mtid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "mrefid", 6) == 0) {
            *end = str+6;
            res->d = b->core.mtid;
            return 0;
        }
        break;

    case 'n':
        if (memcmp(str, "ncigar", 6) == 0) {
            *end = str+6;
            res->d = b->core.n_cigar;
            return 0;
        }
        break;

    case 'p':
        if (memcmp(str, "pos", 3) == 0) {
            *end = str+3;
            res->d = b->core.pos+1;
            return 0;
        } else if (memcmp(str, "pnext", 5) == 0) {
            *end = str+5;
            res->d = b->core.mpos+1;
            return 0;
        }
        break;

    case 'q':
        if (memcmp(str, "qlen", 4) == 0) {
            *end = str+4;
            res->d = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
            return 0;
        } else if (memcmp(str, "qname", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            kputs(bam_get_qname(b), ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "qual", 4) == 0) {
            *end = str+4;
            ks_clear(&res->s);
            if (ks_resize(&res->s, b->core.l_qseq+1) < 0)
                return -1;
            memcpy(res->s.s, bam_get_qual(b), b->core.l_qseq);
            res->s.l = b->core.l_qseq;
            res->is_str = 1;
            return 0;
        }
        break;

    case 'r':
        if (memcmp(str, "rlen", 4) == 0) {
            *end = str+4;
            res->d = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
            return 0;
        } else if (memcmp(str, "rname", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.tid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "rnext", 5) == 0) {
            *end = str+5;
            res->is_str = 1;
            const char *rn = sam_hdr_tid2name(hb->h, b->core.mtid);
            kputs(rn ? rn : "*", ks_clear(&res->s));
            return 0;
        } else if (memcmp(str, "refid", 5) == 0) {
            *end = str+5;
            res->d = b->core.tid;
            return 0;
        }
        break;

    case 's':
        if (memcmp(str, "seq", 3) == 0) {
            *end = str+3;
            ks_clear(&res->s);
            if (ks_resize(&res->s, b->core.l_qseq+1) < 0)
                return -1;
            nibble2base(bam_get_seq(b), res->s.s, b->core.l_qseq);
            res->s.s[b->core.l_qseq] = 0;
            res->s.l = b->core.l_qseq;
            res->is_str = 1;
            return 0;
        } else if (memcmp(str, "sclen", 5) == 0) {
            int sclen = 0;
            uint32_t *cigar = bam_get_cigar(b);
            int ncigar = b->core.n_cigar;
            int left = 0;

            // left
            if (ncigar > 0
                && bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)
                left = 0, sclen += bam_cigar_oplen(cigar[0]);
            else if (ncigar > 1
                     && bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP
                     && bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP)
                left = 1, sclen += bam_cigar_oplen(cigar[1]);

            // right
            if (ncigar-1 > left
                && bam_cigar_op(cigar[ncigar-1]) == BAM_CSOFT_CLIP)
                sclen += bam_cigar_oplen(cigar[ncigar-1]);
            else if (ncigar-2 > left
                     && bam_cigar_op(cigar[ncigar-1]) == BAM_CHARD_CLIP
                     && bam_cigar_op(cigar[ncigar-2]) == BAM_CSOFT_CLIP)
                sclen += bam_cigar_oplen(cigar[ncigar-2]);

            *end = str+5;
            res->d = sclen;
            return 0;
        }
        break;

    case 't':
        if (memcmp(str, "tlen", 4) == 0) {
            *end = str+4;
            res->d = b->core.isize;
            return 0;
        }
        break;

    case '[':
        if (*str == '[' && str[1] && str[2] && str[3] == ']') {
            /* aux tags */
            *end = str+4;

            uint8_t *aux = bam_aux_get(b, str+1);
            if (aux) {
                // we define the truth of a tag to be its presence, even if 0.
                res->is_true = 1;
                switch (*aux) {
                case 'Z':
                case 'H':
                    res->is_str = 1;
                    kputs((char *)aux+1, ks_clear(&res->s));
                    break;

                case 'A':
                    res->is_str = 1;
                    kputsn((char *)aux+1, 1, ks_clear(&res->s));
                    break;

                case 'i': case 'I':
                case 's': case 'S':
                case 'c': case 'C':
                    res->is_str = 0;
                    res->d = bam_aux2i(aux);
                    break;

                case 'f':
                case 'd':
                    res->is_str = 0;
                    res->d = bam_aux2f(aux);
                    break;

                default:
                    hts_log_error("Aux type '%c not yet supported by filters",
                                  *aux);
                    return -1;
                }
                return 0;

            } else {
                // hence absent tags are always false (and strings)
                res->is_str = 1;
                res->s.l = 0;
                res->d = 0;
                res->is_true = 0;
                return 0;
            }
        }
        break;
    }

    // All successful matches in switch should return 0.
    // So if we didn't match, it's a parse error.
    return -1;
}

// Returns 1 when accepted by the filter, 0 if not, -1 on error.
int sam_passes_filter(const sam_hdr_t *h, const bam1_t *b, hts_filter_t *filt)
{
    hb_pair hb = {h, b};
    hts_expr_val_t res = HTS_EXPR_VAL_INIT;
    if (hts_filter_eval2(filt, &hb, bam_sym_lookup, &res)) {
        hts_log_error("Couldn't process filter expression");
        hts_expr_val_free(&res);
        return -1;
    }

    int t = res.is_true;
    hts_expr_val_free(&res);

    return t;
}

static int cram_readrec(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = fpv;
    bam1_t *b = bv;
    int pass_filter, ret;

    do {
        ret = cram_get_bam_seq(fp->fp.cram, &b);
        if (ret < 0)
            return cram_eof(fp->fp.cram) ? -1 : -2;

        if (bam_tag2cigar(b, 1, 1) < 0)
            return -2;

        *tid = b->core.tid;
        *beg = b->core.pos;
        *end = bam_endpos(b);

        if (fp->filter) {
            pass_filter = sam_passes_filter(fp->bam_header, b, fp->filter);
            if (pass_filter < 0)
                return -2;
        } else {
            pass_filter = 1;
        }
    } while (pass_filter == 0);

    return ret;
}

static int cram_pseek(void *fp, int64_t offset, int whence)
{
    cram_fd *fd =  (cram_fd *)fp;

    if ((0 != cram_seek(fd, offset, SEEK_SET))
     && (0 != cram_seek(fd, offset - fd->first_container, SEEK_CUR)))
        return -1;

    fd->curr_position = offset;

    if (fd->ctr) {
        cram_free_container(fd->ctr);
        if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
            cram_free_container(fd->ctr_mt);

        fd->ctr = NULL;
        fd->ctr_mt = NULL;
        fd->ooc = 0;
    }

    return 0;
}

/*
 * cram_ptell is a pseudo-tell function, because it matches the position of the disk cursor only
 *   after a fresh seek call. Otherwise it indicates that the read takes place inside the buffered
 *   container previously fetched. It was designed like this to integrate with the functionality
 *   of the iterator stepping logic.
 */

static int64_t cram_ptell(void *fp)
{
    cram_fd *fd = (cram_fd *)fp;
    cram_container *c;
    cram_slice *s;
    int64_t ret = -1L;

    if (fd) {
        if ((c = fd->ctr) != NULL) {
            if ((s = c->slice) != NULL && s->max_rec) {
                if ((c->curr_slice + s->curr_rec/s->max_rec) >= (c->max_slice + 1))
                    fd->curr_position += c->offset + c->length;
            }
        }
        ret = fd->curr_position;
    }

    return ret;
}

static int bam_pseek(void *fp, int64_t offset, int whence)
{
    BGZF *fd = (BGZF *)fp;

    return bgzf_seek(fd, offset, whence);
}

static int64_t bam_ptell(void *fp)
{
    BGZF *fd = (BGZF *)fp;
    if (!fd)
        return -1L;

    return bgzf_tell(fd);
}



static hts_idx_t *index_load(htsFile *fp, const char *fn, const char *fnidx, int flags)
{
    switch (fp->format.format) {
    case bam:
    case sam:
        return hts_idx_load3(fn, fnidx, HTS_FMT_BAI, flags);

    case cram: {
        if (cram_index_load(fp->fp.cram, fn, fnidx) < 0) return NULL;

        // Cons up a fake "index" just pointing at the associated cram_fd:
        hts_cram_idx_t *idx = malloc(sizeof (hts_cram_idx_t));
        if (idx == NULL) return NULL;
        idx->fmt = HTS_FMT_CRAI;
        idx->cram = fp->fp.cram;
        return (hts_idx_t *) idx;
        }

    default:
        return NULL; // TODO Would use tbx_index_load if it returned hts_idx_t
    }
}

hts_idx_t *sam_index_load3(htsFile *fp, const char *fn, const char *fnidx, int flags)
{
    return index_load(fp, fn, fnidx, flags);
}

hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx) {
    return index_load(fp, fn, fnidx, HTS_IDX_SAVE_REMOTE);
}

hts_idx_t *sam_index_load(htsFile *fp, const char *fn)
{
    return index_load(fp, fn, NULL, HTS_IDX_SAVE_REMOTE);
}

static hts_itr_t *cram_itr_query(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end, hts_readrec_func *readrec)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    hts_itr_t *iter = (hts_itr_t *) calloc(1, sizeof(hts_itr_t));
    if (iter == NULL) return NULL;

    // Cons up a dummy iterator for which hts_itr_next() will simply invoke
    // the readrec function:
    iter->is_cram = 1;
    iter->read_rest = 1;
    iter->off = NULL;
    iter->bins.a = NULL;
    iter->readrec = readrec;

    if (tid >= 0 || tid == HTS_IDX_NOCOOR || tid == HTS_IDX_START) {
        cram_range r = { tid, beg+1, end };
        int ret = cram_set_option(cidx->cram, CRAM_OPT_RANGE, &r);

        iter->curr_off = 0;
        // The following fields are not required by hts_itr_next(), but are
        // filled in in case user code wants to look at them.
        iter->tid = tid;
        iter->beg = beg;
        iter->end = end;

        switch (ret) {
        case 0:
            break;

        case -2:
            // No data vs this ref, so mark iterator as completed.
            // Same as HTS_IDX_NONE.
            iter->finished = 1;
            break;

        default:
            free(iter);
            return NULL;
        }
    }
    else switch (tid) {
    case HTS_IDX_REST:
        iter->curr_off = 0;
        break;
    case HTS_IDX_NONE:
        iter->curr_off = 0;
        iter->finished = 1;
        break;
    default:
        hts_log_error("Query with tid=%d not implemented for CRAM files", tid);
        abort();
        break;
    }

    return iter;
}

hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, hts_pos_t beg, hts_pos_t end)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    if (idx == NULL)
        return hts_itr_query(NULL, tid, beg, end, sam_readrec_rest);
    else if (cidx->fmt == HTS_FMT_CRAI)
        return cram_itr_query(idx, tid, beg, end, sam_readrec);
    else
        return hts_itr_query(idx, tid, beg, end, sam_readrec);
}

static int cram_name2id(void *fdv, const char *ref)
{
    cram_fd *fd = (cram_fd *) fdv;
    return sam_hdr_name2tid(fd->header, ref);
}

hts_itr_t *sam_itr_querys(const hts_idx_t *idx, sam_hdr_t *hdr, const char *region)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    return hts_itr_querys(idx, region, (hts_name2id_f)(bam_name2id), hdr,
                          cidx->fmt == HTS_FMT_CRAI ? cram_itr_query : hts_itr_query,
                          sam_readrec);
}

hts_itr_t *sam_itr_regarray(const hts_idx_t *idx, sam_hdr_t *hdr, char **regarray, unsigned int regcount)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    hts_reglist_t *r_list = NULL;
    int r_count = 0;

    if (!cidx || !hdr)
        return NULL;

    hts_itr_t *itr = NULL;
    if (cidx->fmt == HTS_FMT_CRAI) {
        r_list = hts_reglist_create(regarray, regcount, &r_count, cidx->cram, cram_name2id);
        if (!r_list)
            return NULL;
        itr = hts_itr_regions(idx, r_list, r_count, cram_name2id, cidx->cram,
                   hts_itr_multi_cram, cram_readrec, cram_pseek, cram_ptell);
    } else {
        r_list = hts_reglist_create(regarray, regcount, &r_count, hdr, (hts_name2id_f)(bam_name2id));
        if (!r_list)
            return NULL;
        itr = hts_itr_regions(idx, r_list, r_count, (hts_name2id_f)(bam_name2id), hdr,
                   hts_itr_multi_bam, sam_readrec, bam_pseek, bam_ptell);
    }

    if (!itr)
        hts_reglist_free(r_list, r_count);

    return itr;
}

hts_itr_t *sam_itr_regions(const hts_idx_t *idx, sam_hdr_t *hdr, hts_reglist_t *reglist, unsigned int regcount)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;

    if(!cidx || !hdr || !reglist)
        return NULL;

    if (cidx->fmt == HTS_FMT_CRAI)
        return hts_itr_regions(idx, reglist, regcount, cram_name2id, cidx->cram,
                   hts_itr_multi_cram, cram_readrec, cram_pseek, cram_ptell);
    else
        return hts_itr_regions(idx, reglist, regcount, (hts_name2id_f)(bam_name2id), hdr,
                   hts_itr_multi_bam, sam_readrec, bam_pseek, bam_ptell);
}

/**********************
 *** SAM header I/O ***
 **********************/

#include "htslib/kseq.h"
#include "htslib/kstring.h"

sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text)
{
    sam_hdr_t *bh = sam_hdr_init();
    if (!bh) return NULL;

    if (sam_hdr_add_lines(bh, text, l_text) != 0) {
        sam_hdr_destroy(bh);
        return NULL;
    }

    return bh;
}

static int valid_sam_header_type(const char *s) {
    if (s[0] != '@') return 0;
    switch (s[1]) {
    case 'H':
        return s[2] == 'D' && s[3] == '\t';
    case 'S':
        return s[2] == 'Q' && s[3] == '\t';
    case 'R':
    case 'P':
        return s[2] == 'G' && s[3] == '\t';
    case 'C':
        return s[2] == 'O';
    }
    return 0;
}

// Minimal sanitisation of a header to ensure.
// - null terminated string.
// - all lines start with @ (also implies no blank lines).
//
// Much more could be done, but currently is not, including:
// - checking header types are known (HD, SQ, etc).
// - syntax (eg checking tab separated fields).
// - validating n_targets matches @SQ records.
// - validating target lengths against @SQ records.
static sam_hdr_t *sam_hdr_sanitise(sam_hdr_t *h) {
    if (!h)
        return NULL;

    // Special case for empty headers.
    if (h->l_text == 0)
        return h;

    size_t i;
    unsigned int lnum = 0;
    char *cp = h->text, last = '\n';
    for (i = 0; i < h->l_text; i++) {
        // NB: l_text excludes terminating nul.  This finds early ones.
        if (cp[i] == 0)
            break;

        // Error on \n[^@], including duplicate newlines
        if (last == '\n') {
            lnum++;
            if (cp[i] != '@') {
                hts_log_error("Malformed SAM header at line %u", lnum);
                sam_hdr_destroy(h);
                return NULL;
            }
        }

        last = cp[i];
    }

    if (i < h->l_text) { // Early nul found.  Complain if not just padding.
        size_t j = i;
        while (j < h->l_text && cp[j] == '\0') j++;
        if (j < h->l_text)
            hts_log_warning("Unexpected NUL character in header. Possibly truncated");
    }

    // Add trailing newline and/or trailing nul if required.
    if (last != '\n') {
        hts_log_warning("Missing trailing newline on SAM header. Possibly truncated");

        if (h->l_text < 2 || i >= h->l_text - 2) {
            if (h->l_text >= SIZE_MAX - 2) {
                hts_log_error("No room for extra newline");
                sam_hdr_destroy(h);
                return NULL;
            }

            cp = realloc(h->text, (size_t) h->l_text+2);
            if (!cp) {
                sam_hdr_destroy(h);
                return NULL;
            }
            h->text = cp;
        }
        cp[i++] = '\n';

        // l_text may be larger already due to multiple nul padding
        if (h->l_text < i)
            h->l_text = i;
        cp[h->l_text] = '\0';
    }

    return h;
}

static void known_stderr(const char *tool, const char *advice) {
    hts_log_warning("SAM file corrupted by embedded %s error/log message", tool);
    hts_log_warning("%s", advice);
}

static void warn_if_known_stderr(const char *line) {
    if (strstr(line, "M::bwa_idx_load_from_disk") != NULL)
        known_stderr("bwa", "Use `bwa mem -o file.sam ...` or `bwa sampe -f file.sam ...` instead of `bwa ... > file.sam`");
    else if (strstr(line, "M::mem_pestat") != NULL)
        known_stderr("bwa", "Use `bwa mem -o file.sam ...` instead of `bwa mem ... > file.sam`");
    else if (strstr(line, "loaded/built the index") != NULL)
        known_stderr("minimap2", "Use `minimap2 -o file.sam ...` instead of `minimap2 ... > file.sam`");
}

static sam_hdr_t *sam_hdr_create(htsFile* fp) {
    kstring_t str = { 0, 0, NULL };
    khint_t k;
    sam_hdr_t* h = sam_hdr_init();
    const char *q, *r;
    char* sn = NULL;
    khash_t(s2i) *d = kh_init(s2i);
    khash_t(s2i) *long_refs = NULL;
    if (!h || !d)
        goto error;

    int ret, has_SQ = 0;
    int next_c = '@';
    while (next_c == '@' && (ret = hts_getline(fp, KS_SEP_LINE, &fp->line)) >= 0) {
        if (fp->line.s[0] != '@')
            break;

        if (fp->line.l > 3 && strncmp(fp->line.s, "@SQ", 3) == 0) {
            has_SQ = 1;
            hts_pos_t ln = -1;
            for (q = fp->line.s + 4;; ++q) {
                if (strncmp(q, "SN:", 3) == 0) {
                    q += 3;
                    for (r = q;*r != '\t' && *r != '\n' && *r != '\0';++r);

                    if (sn) {
                        hts_log_warning("SQ header line has more than one SN: tag");
                        free(sn);
                    }
                    sn = (char*)calloc(r - q + 1, 1);
                    if (!sn)
                        goto error;

                    strncpy(sn, q, r - q);
                    q = r;
                } else {
                    if (strncmp(q, "LN:", 3) == 0)
                        ln = strtoll(q + 3, (char**)&q, 10);
                }

                while (*q != '\t' && *q != '\n' && *q != '\0')
                    ++q;
                if (*q == '\0' || *q == '\n')
                    break;
            }
            if (sn) {
                if (ln >= 0) {
                    int absent;
                    k = kh_put(s2i, d, sn, &absent);
                    if (absent < 0)
                        goto error;

                    if (!absent) {
                        hts_log_warning("Duplicated sequence \"%s\" in file \"%s\"", sn, fp->fn);
                        free(sn);
                    } else {
                        sn = NULL;
                        if (ln >= UINT32_MAX) {
                            // Stash away ref length that
                            // doesn't fit in target_len array
                            int k2;
                            if (!long_refs) {
                                long_refs = kh_init(s2i);
                                if (!long_refs)
                                    goto error;
                            }
                            k2 = kh_put(s2i, long_refs, kh_key(d, k), &absent);
                            if (absent < 0)
                                goto error;
                            kh_val(long_refs, k2) = ln;
                            kh_val(d, k) = ((int64_t) (kh_size(d) - 1) << 32
                                            | UINT32_MAX);
                        } else {
                            kh_val(d, k) = (int64_t) (kh_size(d) - 1) << 32 | ln;
                        }
                    }
                } else {
                    hts_log_warning("Ignored @SQ SN:%s : bad or missing LN tag", sn);
                    warn_if_known_stderr(fp->line.s);
                    free(sn);
                }
            } else {
                hts_log_warning("Ignored @SQ line with missing SN: tag");
                warn_if_known_stderr(fp->line.s);
            }
            sn = NULL;
        }
        else if (!valid_sam_header_type(fp->line.s)) {
            hts_log_error("Invalid header line: must start with @HD/@SQ/@RG/@PG/@CO");
            warn_if_known_stderr(fp->line.s);
            goto error;
        }

        if (kputsn(fp->line.s, fp->line.l, &str) < 0)
            goto error;

        if (kputc('\n', &str) < 0)
            goto error;

        if (fp->is_bgzf) {
            next_c = bgzf_peek(fp->fp.bgzf);
        } else {
            unsigned char nc;
            ssize_t pret = hpeek(fp->fp.hfile, &nc, 1);
            next_c = pret > 0 ? nc : pret - 1;
        }
        if (next_c < -1)
            goto error;
    }
    if (next_c != '@')
        fp->line.l = 0;

    if (ret < -1)
        goto error;

    if (!has_SQ && fp->fn_aux) {
        kstring_t line = { 0, 0, NULL };

        /* The reference index (.fai) is actually needed here */
        char *fai_fn = fp->fn_aux;
        char *fn_delim = strstr(fp->fn_aux, HTS_IDX_DELIM);
        if (fn_delim)
            fai_fn = fn_delim + strlen(HTS_IDX_DELIM);

        hFILE* f = hopen(fai_fn, "r");
        int e = 0, absent;
        if (f == NULL)
            goto error;

        while (line.l = 0, kgetline(&line, (kgets_func*) hgets, f) >= 0) {
            char* tab = strchr(line.s, '\t');
            hts_pos_t ln;

            if (tab == NULL)
                continue;

            sn = (char*)calloc(tab-line.s+1, 1);
            if (!sn) {
                e = 1;
                break;
            }
            memcpy(sn, line.s, tab-line.s);
            k = kh_put(s2i, d, sn, &absent);
            if (absent < 0) {
                e = 1;
                break;
            }

            ln = strtoll(tab, NULL, 10);

            if (!absent) {
                hts_log_warning("Duplicated sequence \"%s\" in the file \"%s\"", sn, fai_fn);
                free(sn);
                sn = NULL;
            } else {
                sn = NULL;
                if (ln >= UINT32_MAX) {
                    // Stash away ref length that
                    // doesn't fit in target_len array
                    khint_t k2;
                    int absent = -1;
                    if (!long_refs) {
                        long_refs = kh_init(s2i);
                        if (!long_refs) {
                            e = 1;
                            break;
                        }
                    }
                    k2 = kh_put(s2i, long_refs, kh_key(d, k), &absent);
                    if (absent < 0) {
                         e = 1;
                         break;
                    }
                    kh_val(long_refs, k2) = ln;
                    kh_val(d, k) = ((int64_t) (kh_size(d) - 1) << 32
                                    | UINT32_MAX);
                } else {
                    kh_val(d, k) = (int64_t) (kh_size(d) - 1) << 32 | ln;
                }
                has_SQ = 1;
            }

            e |= kputs("@SQ\tSN:", &str) < 0;
            e |= kputsn(line.s, tab - line.s, &str) < 0;
            e |= kputs("\tLN:", &str) < 0;
            e |= kputll(ln, &str) < 0;
            e |= kputc('\n', &str) < 0;
            if (e)
                break;
        }

        ks_free(&line);
        if (hclose(f) != 0) {
            hts_log_error("Error on closing %s", fai_fn);
            e = 1;
        }
        if (e)
            goto error;
    }

    if (has_SQ) {
        // Populate the targets array
        h->n_targets = kh_size(d);

        h->target_name = (char**) malloc(sizeof(char*) * h->n_targets);
        if (!h->target_name) {
            h->n_targets = 0;
            goto error;
        }

        h->target_len = (uint32_t*) malloc(sizeof(uint32_t) * h->n_targets);
        if (!h->target_len) {
            h->n_targets = 0;
            goto error;
        }

        for (k = kh_begin(d); k != kh_end(d); ++k) {
            if (!kh_exist(d, k))
                continue;

            h->target_name[kh_val(d, k) >> 32] = (char*) kh_key(d, k);
            h->target_len[kh_val(d, k) >> 32] = kh_val(d, k) & 0xffffffffUL;
            kh_val(d, k) >>= 32;
        }
    }

    // Repurpose sdict to hold any references longer than UINT32_MAX
    h->sdict = long_refs;

    kh_destroy(s2i, d);

    if (str.l == 0)
        kputsn("", 0, &str);
    h->l_text = str.l;
    h->text = ks_release(&str);
    fp->bam_header = sam_hdr_sanitise(h);
    fp->bam_header->ref_count = 1;

    return fp->bam_header;

 error:
    if (h && d && (!h->target_name || !h->target_len)) {
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((void *)kh_key(d, k));
    }
    sam_hdr_destroy(h);
    ks_free(&str);
    kh_destroy(s2i, d);
    kh_destroy(s2i, long_refs);
    if (sn) free(sn);
    return NULL;
}

sam_hdr_t *sam_hdr_read(htsFile *fp)
{
    if (!fp) {
        errno = EINVAL;
        return NULL;
    }

    switch (fp->format.format) {
    case bam:
        return sam_hdr_sanitise(bam_hdr_read(fp->fp.bgzf));

    case cram:
        return sam_hdr_sanitise(sam_hdr_dup(fp->fp.cram->header));

    case sam:
        return sam_hdr_create(fp);

    case fastq_format:
    case fasta_format:
        return sam_hdr_init();

    case empty_format:
        errno = EPIPE;
        return NULL;

    default:
        errno = EFTYPE;
        return NULL;
    }
}

int sam_hdr_write(htsFile *fp, const sam_hdr_t *h)
{
    if (!fp || !h) {
        errno = EINVAL;
        return -1;
    }

    switch (fp->format.format) {
    case binary_format:
        fp->format.category = sequence_data;
        fp->format.format = bam;
        /* fall-through */
    case bam:
        if (bam_hdr_write(fp->fp.bgzf, h) < 0) return -1;
        break;

    case cram: {
        cram_fd *fd = fp->fp.cram;
        if (cram_set_header2(fd, h) < 0) return -1;
        if (fp->fn_aux)
            cram_load_reference(fd, fp->fn_aux);
        if (cram_write_SAM_hdr(fd, fd->header) < 0) return -1;
        }
        break;

    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam: {
        if (!h->hrecs && !h->text)
            return 0;
        char *text;
        kstring_t hdr_ks = { 0, 0, NULL };
        size_t l_text;
        ssize_t bytes;
        int r = 0, no_sq = 0;

        if (h->hrecs) {
            if (sam_hrecs_rebuild_text(h->hrecs, &hdr_ks) != 0)
                return -1;
            text = hdr_ks.s;
            l_text = hdr_ks.l;
        } else {
            const char *p = NULL;
            do {
                const char *q = p == NULL ? h->text : p + 4;
                p = strstr(q, "@SQ\t");
            } while (!(p == NULL || p == h->text || *(p - 1) == '\n'));
            no_sq = p == NULL;
            text = h->text;
            l_text = h->l_text;
        }

        if (fp->is_bgzf) {
            bytes = bgzf_write(fp->fp.bgzf, text, l_text);
        } else {
            bytes = hwrite(fp->fp.hfile, text, l_text);
        }
        free(hdr_ks.s);
        if (bytes != l_text)
            return -1;

        if (no_sq) {
            int i;
            for (i = 0; i < h->n_targets; ++i) {
                fp->line.l = 0;
                r |= kputsn("@SQ\tSN:", 7, &fp->line) < 0;
                r |= kputs(h->target_name[i], &fp->line) < 0;
                r |= kputsn("\tLN:", 4, &fp->line) < 0;
                r |= kputw(h->target_len[i], &fp->line) < 0;
                r |= kputc('\n', &fp->line) < 0;
                if (r != 0)
                    return -1;

                if (fp->is_bgzf) {
                    bytes = bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l);
                } else {
                    bytes = hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
                }
                if (bytes != fp->line.l)
                    return -1;
            }
        }
        if (fp->is_bgzf) {
            if (bgzf_flush(fp->fp.bgzf) != 0) return -1;
        } else {
            if (hflush(fp->fp.hfile) != 0) return -1;
        }
        }
        break;

    case fastq_format:
    case fasta_format:
        // Nothing to output; FASTQ has no file headers.
        break;

    default:
        errno = EBADF;
        return -1;
    }
    return 0;
}

static int old_sam_hdr_change_HD(sam_hdr_t *h, const char *key, const char *val)
{
    char *p, *q, *beg = NULL, *end = NULL, *newtext;
    size_t new_l_text;
    if (!h || !key)
        return -1;

    if (h->l_text > 3) {
        if (strncmp(h->text, "@HD", 3) == 0) { //@HD line exists
            if ((p = strchr(h->text, '\n')) == 0) return -1;
            *p = '\0'; // for strstr call

            char tmp[5] = { '\t', key[0], key[0] ? key[1] : '\0', ':', '\0' };

            if ((q = strstr(h->text, tmp)) != 0) { // key exists
                *p = '\n'; // change back

                // mark the key:val
                beg = q;
                for (q += 4; *q != '\n' && *q != '\t'; ++q);
                end = q;

                if (val && (strncmp(beg + 4, val, end - beg - 4) == 0)
                    && strlen(val) == end - beg - 4)
                     return 0; // val is the same, no need to change

            } else {
                beg = end = p;
                *p = '\n';
            }
        }
    }
    if (beg == NULL) { // no @HD
        new_l_text = h->l_text;
        if (new_l_text > SIZE_MAX - strlen(SAM_FORMAT_VERSION) - 9)
            return -1;
        new_l_text += strlen(SAM_FORMAT_VERSION) + 8;
        if (val) {
            if (new_l_text > SIZE_MAX - strlen(val) - 5)
                return -1;
            new_l_text += strlen(val) + 4;
        }
        newtext = (char*)malloc(new_l_text + 1);
        if (!newtext) return -1;

        if (val)
            snprintf(newtext, new_l_text + 1,
                    "@HD\tVN:%s\t%s:%s\n%s", SAM_FORMAT_VERSION, key, val, h->text);
        else
            snprintf(newtext, new_l_text + 1,
                    "@HD\tVN:%s\n%s", SAM_FORMAT_VERSION, h->text);
    } else { // has @HD but different or no key
        new_l_text = (beg - h->text) + (h->text + h->l_text - end);
        if (val) {
            if (new_l_text > SIZE_MAX - strlen(val) - 5)
                return -1;
            new_l_text += strlen(val) + 4;
        }
        newtext = (char*)malloc(new_l_text + 1);
        if (!newtext) return -1;

        if (val) {
            snprintf(newtext, new_l_text + 1, "%.*s\t%s:%s%s",
                    (int) (beg - h->text), h->text, key, val, end);
        } else { //delete key
            snprintf(newtext, new_l_text + 1, "%.*s%s",
                    (int) (beg - h->text), h->text, end);
        }
    }
    free(h->text);
    h->text = newtext;
    h->l_text = new_l_text;
    return 0;
}


int sam_hdr_change_HD(sam_hdr_t *h, const char *key, const char *val)
{
    if (!h || !key)
        return -1;

    if (!h->hrecs)
        return old_sam_hdr_change_HD(h, key, val);

    if (val) {
        if (sam_hdr_update_line(h, "HD", NULL, NULL, key, val, NULL) != 0)
            return -1;
    } else {
        if (sam_hdr_remove_tag_id(h, "HD", NULL, NULL, key) != 0)
            return -1;
    }
    return sam_hdr_rebuild(h);
}
/**********************
 *** SAM record I/O ***
 **********************/

static int sam_parse_B_vals(char type, uint32_t n, char *in, char **end,
                            char *r, bam1_t *b)
{
    int orig_l = b->l_data;
    char *q = in;
    int32_t size;
    size_t bytes;
    int overflow = 0;

    size = aux_type2size(type);
    if (size <= 0 || size > 4) {
        hts_log_error("Unrecognized type B:%c", type);
        return -1;
    }

    // Ensure space for type + values
    bytes = (size_t) n * (size_t) size;
    if (bytes / size != n
        || possibly_expand_bam_data(b, bytes + 2 + sizeof(uint32_t))) {
        hts_log_error("Out of memory");
        return -1;
    }

    b->data[b->l_data++] = 'B';
    b->data[b->l_data++] = type;
    i32_to_le(n, b->data + b->l_data);
    b->l_data += sizeof(uint32_t);
    // This ensures that q always ends up at the next comma after
    // reading a number even if it's followed by junk.  It
    // prevents the possibility of trying to read more than n items.
#define skip_to_comma_(q) do { while (*(q) > '\t' && *(q) != ',') (q)++; } while (0)
    if (type == 'c') {
        while (q < r) {
            *(b->data + b->l_data) = hts_str2int(q + 1, &q, 8, &overflow);
            b->l_data++;
            skip_to_comma_(q);
        }
    } else if (type == 'C') {
        while (q < r) {
            if (*q != '-') {
                *(b->data + b->l_data) = hts_str2uint(q + 1, &q, 8, &overflow);
                b->l_data++;
            } else {
                overflow = 1;
            }
            skip_to_comma_(q);
        }
    } else if (type == 's') {
        while (q < r) {
            i16_to_le(hts_str2int(q + 1, &q, 16, &overflow), b->data + b->l_data);
            b->l_data += 2;
            skip_to_comma_(q);
        }
    } else if (type == 'S') {
        while (q < r) {
            if (*q != '-') {
                u16_to_le(hts_str2uint(q + 1, &q, 16, &overflow), b->data + b->l_data);
                b->l_data += 2;
            } else {
                overflow = 1;
            }
            skip_to_comma_(q);
        }
    } else if (type == 'i') {
        while (q < r) {
            i32_to_le(hts_str2int(q + 1, &q, 32, &overflow), b->data + b->l_data);
            b->l_data += 4;
            skip_to_comma_(q);
        }
    } else if (type == 'I') {
        while (q < r) {
            if (*q != '-') {
                u32_to_le(hts_str2uint(q + 1, &q, 32, &overflow), b->data + b->l_data);
                b->l_data += 4;
            } else {
                overflow = 1;
            }
            skip_to_comma_(q);
        }
    } else if (type == 'f') {
        while (q < r) {
            float_to_le(strtod(q + 1, &q), b->data + b->l_data);
            b->l_data += 4;
            skip_to_comma_(q);
        }
    } else {
        hts_log_error("Unrecognized type B:%c", type);
        return -1;
    }

    if (!overflow) {
        *end = q;
        return 0;
    } else {
        int64_t max = 0, min = 0, val;
        // Given type was incorrect.  Try to rescue the situation.
        q = in;
        overflow = 0;
        b->l_data = orig_l;
        // Find out what range of values is present
        while (q < r) {
            val = hts_str2int(q + 1, &q, 64, &overflow);
            if (max < val) max = val;
            if (min > val) min = val;
            skip_to_comma_(q);
        }
        // Retry with appropriate type
        if (!overflow) {
            if (min < 0) {
                if (min >= INT8_MIN && max <= INT8_MAX) {
                    return sam_parse_B_vals('c', n, in, end, r, b);
                } else if (min >= INT16_MIN && max <= INT16_MAX) {
                    return sam_parse_B_vals('s', n, in, end, r, b);
                } else if (min >= INT32_MIN && max <= INT32_MAX) {
                    return sam_parse_B_vals('i', n, in, end, r, b);
                }
            } else {
                if (max < UINT8_MAX) {
                    return sam_parse_B_vals('C', n, in, end, r, b);
                } else if (max <= UINT16_MAX) {
                    return sam_parse_B_vals('S', n, in, end, r, b);
                } else if (max <= UINT32_MAX) {
                    return sam_parse_B_vals('I', n, in, end, r, b);
                }
            }
        }
        // If here then at least one of the values is too big to store
        hts_log_error("Numeric value in B array out of allowed range");
        return -1;
    }
#undef skip_to_comma_
}

static inline unsigned int parse_sam_flag(char *v, char **rv, int *overflow) {
    if (*v >= '1' && *v <= '9') {
        return hts_str2uint(v, rv, 16, overflow);
    }
    else if (*v == '0') {
        // handle single-digit "0" directly; otherwise it's hex or octal
        if (v[1] == '\t') { *rv = v+1; return 0; }
        else {
            unsigned long val = strtoul(v, rv, 0);
            if (val > 65535) { *overflow = 1; return 65535; }
            return val;
        }
    }
    else {
        // TODO implement symbolic flag letters
        *rv = v;
        return 0;
    }
}

// Parse tag line and append to bam object b.
// Shared by both SAM and FASTQ parsers.
//
// The difference between the two is how lenient we are to recognising
// non-compliant strings.  The FASTQ parser glosses over arbitrary
// non-SAM looking strings.
static inline int aux_parse(char *start, char *end, bam1_t *b, int lenient,
                            khash_t(tag) *tag_whitelist) {
    int overflow = 0;
    int checkpoint;
    char logbuf[40];
    char *q = start, *p = end;

#define _parse_err(cond, ...)                   \
    do {                                        \
        if (cond) {                             \
            if (lenient) {                      \
                while (q < p && !isspace_c(*q))   \
                    q++;                        \
                while (q < p && isspace_c(*q))    \
                    q++;                        \
                b->l_data = checkpoint;         \
                goto loop;                      \
            } else {                            \
                hts_log_error(__VA_ARGS__);     \
                goto err_ret;                   \
            }                                   \
        }                                       \
    } while (0)

    while (q < p) loop: {
        char type;
        checkpoint = b->l_data;
        if (p - q < 5) {
            if (lenient) {
                break;
            } else {
                hts_log_error("Incomplete aux field");
                goto err_ret;
            }
        }
        _parse_err(q[0] < '!' || q[1] < '!', "invalid aux tag id");

        if (lenient && (q[2] | q[4]) != ':') {
            while (q < p && !isspace_c(*q))
                q++;
            while (q < p && isspace_c(*q))
                q++;
            continue;
        }

        if (tag_whitelist) {
            int tt = q[0]*256 + q[1];
            if (kh_get(tag, tag_whitelist, tt) == kh_end(tag_whitelist)) {
                while (q < p && *q != '\t')
                    q++;
                continue;
            }
        }

        // Copy over id
        if (possibly_expand_bam_data(b, 2) < 0) goto err_ret;
        memcpy(b->data + b->l_data, q, 2); b->l_data += 2;
        q += 3; type = *q++; ++q; // q points to value
        if (type != 'Z' && type != 'H') // the only zero length acceptable fields
            _parse_err(*q <= '\t', "incomplete aux field");

        // Ensure enough space for a double + type allocated.
        if (possibly_expand_bam_data(b, 16) < 0) goto err_ret;

        if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
            b->data[b->l_data++] = 'A';
            b->data[b->l_data++] = *q++;
        } else if (type == 'i' || type == 'I') {
            if (*q == '-') {
                int32_t x = hts_str2int(q, &q, 32, &overflow);
                if (x >= INT8_MIN) {
                    b->data[b->l_data++] = 'c';
                    b->data[b->l_data++] = x;
                } else if (x >= INT16_MIN) {
                    b->data[b->l_data++] = 's';
                    i16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'i';
                    i32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            } else {
                uint32_t x = hts_str2uint(q, &q, 32, &overflow);
                if (x <= UINT8_MAX) {
                    b->data[b->l_data++] = 'C';
                    b->data[b->l_data++] = x;
                } else if (x <= UINT16_MAX) {
                    b->data[b->l_data++] = 'S';
                    u16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'I';
                    u32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            }
        } else if (type == 'f') {
            b->data[b->l_data++] = 'f';
            float_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(float);
        } else if (type == 'd') {
            b->data[b->l_data++] = 'd';
            double_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(double);
        } else if (type == 'Z' || type == 'H') {
            char *end = strchr(q, '\t');
            if (!end) end = q + strlen(q);
            _parse_err(type == 'H' && ((end-q)&1) != 0,
                       "hex field does not have an even number of digits");
            b->data[b->l_data++] = type;
            if (possibly_expand_bam_data(b, end - q + 1) < 0) goto err_ret;
            memcpy(b->data + b->l_data, q, end - q);
            b->l_data += end - q;
            b->data[b->l_data++] = '\0';
            q = end;
        } else if (type == 'B') {
            uint32_t n;
            char *r;
            type = *q++; // q points to the first ',' following the typing byte
            _parse_err(*q && *q != ',' && *q != '\t',
                       "B aux field type not followed by ','");

            for (r = q, n = 0; *r > '\t'; ++r)
                if (*r == ',') ++n;

            if (sam_parse_B_vals(type, n, q, &q, r, b) < 0)
                goto err_ret;
        } else _parse_err(1, "unrecognized type %s", hts_strprint(logbuf, sizeof logbuf, '\'', &type, 1));

        while (*q > '\t') { q++; } // Skip any junk to next tab
        q++;
    }

    _parse_err(!lenient && overflow != 0, "numeric value out of allowed range");
#undef _parse_err

    return 0;

err_ret:
    return -2;
}

int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret; *tab = '\0'; (_p) = tab + 1; } while (0)

#if HTS_ALLOW_UNALIGNED != 0 && ULONG_MAX == 0xffffffffffffffff

// Macro that operates on 64-bits at a time.
#define COPY_MINUS_N(to,from,n,l,failed)                        \
    do {                                                        \
        uint64_u *from8 = (uint64_u *)(from);                   \
        uint64_u *to8 = (uint64_u *)(to);                       \
        uint64_t uflow = 0;                                     \
        size_t l8 = (l)>>3, i;                                  \
        for (i = 0; i < l8; i++) {                              \
            to8[i] = from8[i] - (n)*0x0101010101010101UL;       \
            uflow |= to8[i];                                    \
        }                                                       \
        for (i<<=3; i < (l); ++i) {                             \
            to[i] = from[i] - (n);                              \
            uflow |= to[i];                                     \
        }                                                       \
        failed = (uflow & 0x8080808080808080UL) > 0;            \
    } while (0)

#else

// Basic version which operates a byte at a time
#define COPY_MINUS_N(to,from,n,l,failed) do {                \
        uint8_t uflow = 0;                                   \
        for (i = 0; i < (l); ++i) {                          \
            (to)[i] = (from)[i] - (n);                       \
            uflow |= (uint8_t) (to)[i];                      \
        }                                                    \
        failed = (uflow & 0x80) > 0;                         \
    } while (0)

#endif

#define _get_mem(type_t, x, b, l) if (possibly_expand_bam_data((b), (l)) < 0) goto err_ret; *(x) = (type_t*)((b)->data + (b)->l_data); (b)->l_data += (l)
#define _parse_err(cond, ...) do { if (cond) { hts_log_error(__VA_ARGS__); goto err_ret; } } while (0)
#define _parse_warn(cond, ...) do { if (cond) { hts_log_warning(__VA_ARGS__); } } while (0)

    uint8_t *t;

    char *p = s->s, *q;
    int i, overflow = 0;
    char logbuf[40];
    hts_pos_t cigreflen;
    bam1_core_t *c = &b->core;

    b->l_data = 0;
    memset(c, 0, 32);

    // qname
    q = _read_token(p);

    _parse_warn(p - q <= 1, "empty query name");
    _parse_err(p - q > 255, "query name too long");
    // resize large enough for name + extranul
    if (possibly_expand_bam_data(b, (p - q) + 4) < 0) goto err_ret;
    memcpy(b->data + b->l_data, q, p-q); b->l_data += p-q;

    c->l_extranul = (4 - (b->l_data & 3)) & 3;
    memcpy(b->data + b->l_data, "\0\0\0\0", c->l_extranul);
    b->l_data += c->l_extranul;

    c->l_qname = p - q + c->l_extranul;

    // flag
    c->flag = parse_sam_flag(p, &p, &overflow);
    if (*p++ != '\t') goto err_ret; // malformated flag

    // chr
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(h->n_targets == 0, "no SQ lines present in the header");
        c->tid = bam_name2id(h, q);
        _parse_err(c->tid < -1, "failed to parse header");
        _parse_warn(c->tid < 0, "unrecognized reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    } else c->tid = -1;

    // pos
    c->pos = hts_str2uint(p, &p, 63, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->pos < 0 && c->tid >= 0) {
        _parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
        c->tid = -1;
    }
    if (c->tid < 0) c->flag |= BAM_FUNMAP;

    // mapq
    c->qual = hts_str2uint(p, &p, 8, &overflow);
    if (*p++ != '\t') goto err_ret;
    // cigar
    if (*p != '*') {
        uint32_t *cigar = NULL;
        int old_l_data = b->l_data;
        int n_cigar = bam_parse_cigar(p, &p, b);
        if (n_cigar < 1 || *p++ != '\t') goto err_ret;
        cigar = (uint32_t *)(b->data + old_l_data);
        c->n_cigar = n_cigar;

        // can't use bam_endpos() directly as some fields not yet set up
        cigreflen = (!(c->flag&BAM_FUNMAP))? bam_cigar2rlen(c->n_cigar, cigar) : 1;
        if (cigreflen == 0) cigreflen = 1;
    } else {
        _parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
        c->flag |= BAM_FUNMAP;
        q = _read_token(p);
        cigreflen = 1;
    }
    _parse_err(HTS_POS_MAX - cigreflen <= c->pos,
               "read ends beyond highest supported position");
    c->bin = hts_reg2bin(c->pos, c->pos + cigreflen, 14, 5);
    // mate chr
    q = _read_token(p);
    if (strcmp(q, "=") == 0) {
        c->mtid = c->tid;
    } else if (strcmp(q, "*") == 0) {
        c->mtid = -1;
    } else {
        c->mtid = bam_name2id(h, q);
        _parse_err(c->mtid < -1, "failed to parse header");
        _parse_warn(c->mtid < 0, "unrecognized mate reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    }
    // mpos
    c->mpos = hts_str2uint(p, &p, 63, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->mpos < 0 && c->mtid >= 0) {
        _parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
        c->mtid = -1;
    }
    // tlen
    c->isize = hts_str2int(p, &p, 64, &overflow);
    if (*p++ != '\t') goto err_ret;
    // seq
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(p - q - 1 > INT32_MAX, "read sequence is too long");
        c->l_qseq = p - q - 1;
        hts_pos_t ql = bam_cigar2qlen(c->n_cigar, (uint32_t*)(b->data + c->l_qname));
        _parse_err(c->n_cigar && ql != c->l_qseq, "CIGAR and query sequence are of different length");
        i = (c->l_qseq + 1) >> 1;
        _get_mem(uint8_t, &t, b, i);

        unsigned int lqs2 = c->l_qseq&~1, i;
        for (i = 0; i < lqs2; i+=2)
            t[i>>1] = (seq_nt16_table[(unsigned char)q[i]] << 4) | seq_nt16_table[(unsigned char)q[i+1]];
        for (; i < c->l_qseq; ++i)
            t[i>>1] = seq_nt16_table[(unsigned char)q[i]] << ((~i&1)<<2);
    } else c->l_qseq = 0;
    // qual
    _get_mem(uint8_t, &t, b, c->l_qseq);
    if (p[0] == '*' && (p[1] == '\t' || p[1] == '\0')) {
        memset(t, 0xff, c->l_qseq);
        p += 2;
    } else {
        int failed = 0;
        _parse_err(s->l - (p - s->s) < c->l_qseq
                   || (p[c->l_qseq] != '\t' && p[c->l_qseq] != '\0'),
                   "SEQ and QUAL are of different length");
        COPY_MINUS_N(t, p, 33, c->l_qseq, failed);
        _parse_err(failed, "invalid QUAL character");
        p += c->l_qseq + 1;
    }

    // aux
    if (aux_parse(p, s->s + s->l, b, 0, NULL) < 0)
        goto err_ret;

    if (bam_tag2cigar(b, 1, 1) < 0)
        return -2;
    return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token
err_ret:
    return -2;
}

static uint32_t read_ncigar(const char *q) {
    uint32_t n_cigar = 0;
    for (; *q && *q != '\t'; ++q)
        if (!isdigit_c(*q)) ++n_cigar;
    if (!n_cigar) {
        hts_log_error("No CIGAR operations");
        return 0;
    }
    if (n_cigar >= 2147483647) {
        hts_log_error("Too many CIGAR operations");
        return 0;
    }

    return n_cigar;
}

/*! @function
 @abstract  Parse a CIGAR string into preallocated a uint32_t array
 @param  in      [in]  pointer to the source string
 @param  a_cigar [out]  address of the destination uint32_t buffer
 @return         number of processed input characters; 0 on error
 */
static int parse_cigar(const char *in, uint32_t *a_cigar, uint32_t n_cigar) {
    int i, overflow = 0;
    const char *p = in;
    for (i = 0; i < n_cigar; i++) {
        uint32_t len;
        int op;
        char *q;
        len = hts_str2uint(p, &q, 28, &overflow)<<BAM_CIGAR_SHIFT;
        if (q == p) {
            hts_log_error("CIGAR length invalid at position %d (%s)", (int)(i+1), p);
            return 0;
        }
        if (overflow) {
            hts_log_error("CIGAR length too long at position %d (%.*s)", (int)(i+1), (int)(q-p+1), p);
            return 0;
        }
        p = q;
        op = bam_cigar_table[(unsigned char)*p++];
        if (op < 0) {
            hts_log_error("Unrecognized CIGAR operator");
            return 0;
        }
        a_cigar[i] = len;
        a_cigar[i] |= op;
    }

    return p-in;
}

ssize_t sam_parse_cigar(const char *in, char **end, uint32_t **a_cigar, size_t *a_mem) {
    size_t n_cigar = 0;
    int diff;

    if (!in || !a_cigar || !a_mem) {
        hts_log_error("NULL pointer arguments");
        return -1;
    }
    if (end) *end = (char *)in;

    if (*in == '*') {
        if (end) (*end)++;
        return 0;
    }
    n_cigar = read_ncigar(in);
    if (!n_cigar) return 0;
    if (n_cigar > *a_mem) {
        uint32_t *a_tmp = realloc(*a_cigar, n_cigar*sizeof(**a_cigar));
        if (a_tmp) {
            *a_cigar = a_tmp;
            *a_mem = n_cigar;
        } else {
            hts_log_error("Memory allocation error");
            return -1;
        }
    }

    if (!(diff = parse_cigar(in, *a_cigar, n_cigar))) return -1;
    if (end) *end = (char *)in+diff;

    return n_cigar;
}

ssize_t bam_parse_cigar(const char *in, char **end, bam1_t *b) {
    size_t n_cigar = 0;
    int diff;

    if (!in || !b) {
        hts_log_error("NULL pointer arguments");
        return -1;
    }
    if (end) *end = (char *)in;

    if (*in == '*') {
        if (end) (*end)++;
        return 0;
    }
    n_cigar = read_ncigar(in);
    if (!n_cigar) return 0;
    if (possibly_expand_bam_data(b, n_cigar * sizeof(uint32_t)) < 0) {
        hts_log_error("Memory allocation error");
        return -1;
    }

    if (!(diff = parse_cigar(in, (uint32_t *)(b->data + b->l_data), n_cigar))) return -1;
    b->l_data += (n_cigar * sizeof(uint32_t));
    if (end) *end = (char *)in+diff;

    return n_cigar;
}

/*
 * -----------------------------------------------------------------------------
 * SAM threading
 */
// Size of SAM text block (reading)
#define SAM_NBYTES 240000

// Number of BAM records (writing, up to NB_mem in size)
#define SAM_NBAM 1000

struct SAM_state;

// Output job - a block of BAM records
typedef struct sp_bams {
    struct sp_bams *next;
    int serial;

    bam1_t *bams;
    int nbams, abams; // used and alloc for bams[] array
    size_t bam_mem;   // very approximate total size

    struct SAM_state *fd;
} sp_bams;

// Input job - a block of SAM text
typedef struct sp_lines {
    struct sp_lines *next;
    int serial;

    char *data;
    int data_size;
    int alloc;

    struct SAM_state *fd;
    sp_bams *bams;
} sp_lines;

enum sam_cmd {
    SAM_NONE = 0,
    SAM_CLOSE,
    SAM_CLOSE_DONE,
};

typedef struct SAM_state {
    sam_hdr_t *h;

    hts_tpool *p;
    int own_pool;
    pthread_mutex_t lines_m;
    hts_tpool_process *q;
    pthread_t dispatcher;
    int dispatcher_set;

    sp_lines *lines;
    sp_bams *bams;

    sp_bams *curr_bam;
    int curr_idx;
    int serial;

    // Be warned: moving these mutexes around in this struct can reduce
    // threading performance by up to 70%!
    pthread_mutex_t command_m;
    pthread_cond_t command_c;
    enum sam_cmd command;

    // One of the E* errno codes
    int errcode;

    htsFile *fp;
} SAM_state;

// Returns a SAM_state struct from a generic hFILE.
//
// Returns NULL on failure.
static SAM_state *sam_state_create(htsFile *fp) {
    // Ideally sam_open wouldn't be a #define to hts_open but instead would
    // be a redirect call with an additional 'S' mode.  This in turn would
    // correctly set the designed format to sam instead of a generic
    // text_format.
    if (fp->format.format != sam && fp->format.format != text_format)
        return NULL;

    SAM_state *fd = calloc(1, sizeof(*fd));
    if (!fd)
        return NULL;

    fp->state = fd;
    fd->fp = fp;

    return fd;
}

static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str);
static void *sam_format_worker(void *arg);

static void sam_state_err(SAM_state *fd, int errcode) {
    pthread_mutex_lock(&fd->command_m);
    if (!fd->errcode)
        fd->errcode = errcode;
    pthread_mutex_unlock(&fd->command_m);
}

static void sam_free_sp_bams(sp_bams *b) {
    if (!b)
        return;

    if (b->bams) {
        int i;
        for (i = 0; i < b->abams; i++) {
            if (b->bams[i].data)
                free(b->bams[i].data);
        }
        free(b->bams);
    }
    free(b);
}

// Destroys the state produce by sam_state_create.
int sam_state_destroy(htsFile *fp) {
    int ret = 0;

    if (!fp->state)
        return 0;

    SAM_state *fd = fp->state;
    if (fd->p) {
        if (fd->h) {
            // Notify sam_dispatcher we're closing
            pthread_mutex_lock(&fd->command_m);
            if (fd->command != SAM_CLOSE_DONE)
                fd->command = SAM_CLOSE;
            pthread_cond_signal(&fd->command_c);
            ret = -fd->errcode;
            if (fd->q)
                hts_tpool_wake_dispatch(fd->q); // unstick the reader

            if (!fp->is_write && fd->q && fd->dispatcher_set) {
                for (;;) {
                    // Avoid deadlocks with dispatcher
                    if (fd->command == SAM_CLOSE_DONE)
                        break;
                    hts_tpool_wake_dispatch(fd->q);
                    pthread_mutex_unlock(&fd->command_m);
                    usleep(10000);
                    pthread_mutex_lock(&fd->command_m);
                }
            }
            pthread_mutex_unlock(&fd->command_m);

            if (fp->is_write) {
                // Dispatch the last partial block.
                sp_bams *gb = fd->curr_bam;
                if (!ret && gb && gb->nbams > 0 && fd->q)
                    ret = hts_tpool_dispatch(fd->p, fd->q, sam_format_worker, gb);

                // Flush and drain output
                if (fd->q)
                    hts_tpool_process_flush(fd->q);
                pthread_mutex_lock(&fd->command_m);
                if (!ret) ret = -fd->errcode;
                pthread_mutex_unlock(&fd->command_m);

                while (!ret && fd->q && !hts_tpool_process_empty(fd->q)) {
                    usleep(10000);
                    pthread_mutex_lock(&fd->command_m);
                    ret = -fd->errcode;
                    // not empty but shutdown implies error
                    if (hts_tpool_process_is_shutdown(fd->q) && !ret)
                        ret = EIO;
                    pthread_mutex_unlock(&fd->command_m);
                }
                if (fd->q)
                    hts_tpool_process_shutdown(fd->q);
            }

            // Wait for it to acknowledge
            if (fd->dispatcher_set)
                pthread_join(fd->dispatcher, NULL);
            if (!ret) ret = -fd->errcode;
        }

        // Tidy up memory
        if (fd->q)
            hts_tpool_process_destroy(fd->q);

        if (fd->own_pool && fp->format.compression == no_compression) {
            hts_tpool_destroy(fd->p);
            fd->p = NULL;
        }
        pthread_mutex_destroy(&fd->lines_m);
        pthread_mutex_destroy(&fd->command_m);
        pthread_cond_destroy(&fd->command_c);

        sp_lines *l = fd->lines;
        while (l) {
            sp_lines *n = l->next;
            free(l->data);
            free(l);
            l = n;
        }

        sp_bams *b = fd->bams;
        while (b) {
            if (fd->curr_bam == b)
                fd->curr_bam = NULL;
            sp_bams *n = b->next;
            sam_free_sp_bams(b);
            b = n;
        }

        if (fd->curr_bam)
            sam_free_sp_bams(fd->curr_bam);

        // Decrement counter by one, maybe destroying too.
        // This is to permit the caller using bam_hdr_destroy
        // before sam_close without triggering decode errors
        // in the background threads.
        bam_hdr_destroy(fd->h);
    }

    free(fp->state);
    fp->state = NULL;
    return ret;
}

// Cleanup function - job for sam_parse_worker; result for sam_format_worker
static void cleanup_sp_lines(void *arg) {
    sp_lines *gl = (sp_lines *)arg;
    if (!gl) return;

    // Should always be true for lines passed to / from thread workers.
    assert(gl->next == NULL);

    free(gl->data);
    sam_free_sp_bams(gl->bams);
    free(gl);
}

// Run from one of the worker threads.
// Convert a passed in array of lines to array of BAMs, returning
// the result back to the thread queue.
static void *sam_parse_worker(void *arg) {
    sp_lines *gl = (sp_lines *)arg;
    sp_bams *gb = NULL;
    char *lines = gl->data;
    int i;
    bam1_t *b;
    SAM_state *fd = gl->fd;

    // Use a block of BAM structs we had earlier if available.
    pthread_mutex_lock(&fd->lines_m);
    if (fd->bams) {
        gb = fd->bams;
        fd->bams = gb->next;
    }
    pthread_mutex_unlock(&fd->lines_m);

    if (gb == NULL) {
        gb = calloc(1, sizeof(*gb));
        if (!gb) {
            return NULL;
        }
        gb->abams = 100;
        gb->bams = b = calloc(gb->abams, sizeof(*b));
        if (!gb->bams) {
            sam_state_err(fd, ENOMEM);
            goto err;
        }
        gb->nbams = 0;
        gb->bam_mem = 0;
    }
    gb->serial = gl->serial;
    gb->next = NULL;

    b = (bam1_t *)gb->bams;
    if (!b) {
        sam_state_err(fd, ENOMEM);
        goto err;
    }

    i = 0;
    char *cp = lines, *cp_end = lines + gl->data_size;
    while (cp < cp_end) {
        if (i >= gb->abams) {
            int old_abams = gb->abams;
            gb->abams *= 2;
            b = (bam1_t *)realloc(gb->bams, gb->abams*sizeof(bam1_t));
            if (!b) {
                gb->abams /= 2;
                sam_state_err(fd, ENOMEM);
                goto err;
            }
            memset(&b[old_abams], 0, (gb->abams - old_abams)*sizeof(*b));
            gb->bams = b;
        }

        // Ideally we'd get sam_parse1 to return the number of
        // bytes decoded and to be able to stop on newline as
        // well as \0.
        //
        // We can then avoid the additional strchr loop.
        // It's around 6% of our CPU cost, albeit threadable.
        //
        // However this is an API change so for now we copy.

        char *nl = strchr(cp, '\n');
        char *line_end;
        if (nl) {
            line_end = nl;
            if (line_end > cp && *(line_end - 1) == '\r')
                line_end--;
            nl++;
        } else {
            nl = line_end = cp_end;
        }
        *line_end = '\0';
        kstring_t ks = { line_end - cp, gl->alloc, cp };
        if (sam_parse1(&ks, fd->h, &b[i]) < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            cleanup_sp_lines(gl);
            goto err;
        }

        cp = nl;
        i++;
    }
    gb->nbams = i;

    pthread_mutex_lock(&fd->lines_m);
    gl->next = fd->lines;
    fd->lines = gl;
    pthread_mutex_unlock(&fd->lines_m);
    return gb;

 err:
    sam_free_sp_bams(gb);
    return NULL;
}

static void *sam_parse_eof(void *arg) {
    return NULL;
}

// Cleanup function - result for sam_parse_worker; job for sam_format_worker
static void cleanup_sp_bams(void *arg) {
    sam_free_sp_bams((sp_bams *) arg);
}

// Runs in its own thread.
// Reads a block of text (SAM) and sends a new job to the thread queue to
// translate this to BAM.
static void *sam_dispatcher_read(void *vp) {
    htsFile *fp = vp;
    kstring_t line = {0};
    int line_frag = 0;
    SAM_state *fd = fp->state;
    sp_lines *l = NULL;

    // Pre-allocate buffer for left-over bits of line (exact size doesn't
    // matter as it will grow if necessary).
    if (ks_resize(&line, 1000) < 0)
        goto err;

    for (;;) {
        // Check for command
        pthread_mutex_lock(&fd->command_m);
        switch (fd->command) {

        case SAM_CLOSE:
            pthread_cond_signal(&fd->command_c);
            pthread_mutex_unlock(&fd->command_m);
            hts_tpool_process_shutdown(fd->q);
            goto tidyup;

        default:
            break;
        }
        pthread_mutex_unlock(&fd->command_m);

        pthread_mutex_lock(&fd->lines_m);
        if (fd->lines) {
            // reuse existing line buffer
            l = fd->lines;
            fd->lines = l->next;
        }
        pthread_mutex_unlock(&fd->lines_m);

        if (l == NULL) {
            // none to reuse, to create a new one
            l = calloc(1, sizeof(*l));
            if (!l)
                goto err;
            l->alloc = SAM_NBYTES;
            l->data = malloc(l->alloc+8); // +8 for optimisation in sam_parse1
            if (!l->data) {
                free(l);
                l = NULL;
                goto err;
            }
            l->fd = fd;
        }
        l->next = NULL;

        if (l->alloc < line_frag+SAM_NBYTES/2) {
            char *rp = realloc(l->data, line_frag+SAM_NBYTES/2 +8);
            if (!rp)
                goto err;
            l->alloc = line_frag+SAM_NBYTES/2;
            l->data = rp;
        }
        memcpy(l->data, line.s, line_frag);

        l->data_size = line_frag;
        ssize_t nbytes;
    longer_line:
        if (fp->is_bgzf)
            nbytes = bgzf_read(fp->fp.bgzf, l->data + line_frag, l->alloc - line_frag);
        else
            nbytes = hread(fp->fp.hfile, l->data + line_frag, l->alloc - line_frag);
        if (nbytes < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            goto err;
        } else if (nbytes == 0)
            break; // EOF
        l->data_size += nbytes;

        // trim to last \n. Maybe \r\n, but that's still fine
        if (nbytes == l->alloc - line_frag) {
            char *cp_end = l->data + l->data_size;
            char *cp = cp_end-1;

            while (cp > (char *)l->data && *cp != '\n')
                cp--;

            // entire buffer is part of a single line
            if (cp == l->data) {
                line_frag = l->data_size;
                char *rp = realloc(l->data, l->alloc * 2 + 8);
                if (!rp)
                    goto err;
                l->alloc *= 2;
                l->data = rp;
                assert(l->alloc >= l->data_size);
                assert(l->alloc >= line_frag);
                assert(l->alloc >= l->alloc - line_frag);
                goto longer_line;
            }
            cp++;

            // line holds the remainder of our line.
            if (ks_resize(&line, cp_end - cp) < 0)
                goto err;
            memcpy(line.s, cp, cp_end - cp);
            line_frag = cp_end - cp;
            l->data_size = l->alloc - line_frag;
        } else {
            // out of buffer
            line_frag = 0;
        }

        l->serial = fd->serial++;
        //fprintf(stderr, "Dispatching %p, %d bytes, serial %d\n", l, l->data_size, l->serial);
        if (hts_tpool_dispatch3(fd->p, fd->q, sam_parse_worker, l,
                                cleanup_sp_lines, cleanup_sp_bams, 0) < 0)
            goto err;
        pthread_mutex_lock(&fd->command_m);
        if (fd->command == SAM_CLOSE) {
            pthread_mutex_unlock(&fd->command_m);
            l = NULL;
            goto tidyup;
        }
        l = NULL;  // Now "owned" by sam_parse_worker()
        pthread_mutex_unlock(&fd->command_m);
    }

    if (hts_tpool_dispatch(fd->p, fd->q, sam_parse_eof, NULL) < 0)
        goto err;

    // At EOF, wait for close request.
    // (In future if we add support for seek, this is where we need to catch it.)
    for (;;) {
        pthread_mutex_lock(&fd->command_m);
        if (fd->command == SAM_NONE)
            pthread_cond_wait(&fd->command_c, &fd->command_m);
        switch (fd->command) {
        case SAM_CLOSE:
            pthread_cond_signal(&fd->command_c);
            pthread_mutex_unlock(&fd->command_m);
            hts_tpool_process_shutdown(fd->q);
            goto tidyup;

        default:
            pthread_mutex_unlock(&fd->command_m);
            break;
        }
    }

 tidyup:
    pthread_mutex_lock(&fd->command_m);
    fd->command = SAM_CLOSE_DONE;
    pthread_cond_signal(&fd->command_c);
    pthread_mutex_unlock(&fd->command_m);

    if (l) {
        pthread_mutex_lock(&fd->lines_m);
        l->next = fd->lines;
        fd->lines = l;
        pthread_mutex_unlock(&fd->lines_m);
    }
    free(line.s);

    return NULL;

 err:
    sam_state_err(fd, errno ? errno : ENOMEM);
    hts_tpool_process_shutdown(fd->q);
    goto tidyup;
}

// Runs in its own thread.
// Takes encoded blocks of SAM off the thread results queue and writes them
// to our output stream.
static void *sam_dispatcher_write(void *vp) {
    htsFile *fp = vp;
    SAM_state *fd = fp->state;
    hts_tpool_result *r;

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(fd->q))) {
        sp_lines *gl = (sp_lines *)hts_tpool_result_data(r);
        if (!gl) {
            sam_state_err(fd, ENOMEM);
            goto err;
        }

        if (fp->idx) {
            sp_bams *gb = gl->bams;
            int i = 0, count = 0;
            while (i < gl->data_size) {
                int j = i;
                while (i < gl->data_size && gl->data[i] != '\n')
                    i++;
                if (i < gl->data_size)
                    i++;

                if (fp->is_bgzf) {
                    if (bgzf_flush_try(fp->fp.bgzf, i-j) < 0)
                        goto err;
                    if (bgzf_write(fp->fp.bgzf, &gl->data[j], i-j) != i-j)
                        goto err;
                } else {
                    if (hwrite(fp->fp.hfile, &gl->data[j], i-j) != i-j)
                        goto err;
                }

                bam1_t *b = &gb->bams[count++];
                if (fp->format.compression == bgzf) {
                    if (bgzf_idx_push(fp->fp.bgzf, fp->idx,
                                      b->core.tid, b->core.pos, bam_endpos(b),
                                      bgzf_tell(fp->fp.bgzf),
                                      !(b->core.flag&BAM_FUNMAP)) < 0) {
                        sam_state_err(fd, errno ? errno : ENOMEM);
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(fd->h, b->core.tid), sam_hdr_tid2len(fd->h, b->core.tid), b->core.flag, b->core.pos+1);
                        goto err;
                    }
                } else {
                    if (hts_idx_push(fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                     bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        sam_state_err(fd, errno ? errno : ENOMEM);
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(fd->h, b->core.tid), sam_hdr_tid2len(fd->h, b->core.tid), b->core.flag, b->core.pos+1);
                        goto err;
                    }
                }
            }

            assert(count == gb->nbams);

            // Add bam array to free-list
            pthread_mutex_lock(&fd->lines_m);
            gb->next = fd->bams;
            fd->bams = gl->bams;
            gl->bams = NULL;
            pthread_mutex_unlock(&fd->lines_m);
        } else {
            if (fp->is_bgzf) {
                // We keep track of how much in the current block we have
                // remaining => R.  We look for the last newline in input
                // [i] to [i+R], backwards => position N.
                //
                // If we find a newline, we write out bytes i to N.
                // We know we cannot fit the next record in this bgzf block,
                // so we flush what we have and copy input N to i+R into
                // the start of a new block, and recompute a new R for that.
                //
                // If we don't find a newline (i==N) then we cannot extend
                // the current block at all, so flush whatever is in it now
                // if it ends on a newline.
                // We still copy i(==N) to i+R to the next block and
                // continue as before with a new R.
                //
                // The only exception on the flush is when we run out of
                // data in the input.  In that case we skip it as we don't
                // yet know if the next record will fit.
                //
                // Both conditions share the same code here:
                // - Look for newline (pos N)
                // - Write i to N (which maybe 0)
                // - Flush if block ends on newline and not end of input
                // - write N to i+R

                int i = 0;
                BGZF *fb = fp->fp.bgzf;
                while (i < gl->data_size) {
                    // remaining space in block
                    int R = BGZF_BLOCK_SIZE - fb->block_offset;
                    int eod = 0;
                    if (R > gl->data_size-i)
                        R = gl->data_size-i, eod = 1;

                    // Find last newline in input data
                    int N = i + R;
                    while (--N > i) {
                        if (gl->data[N] == '\n')
                            break;
                    }

                    if (N != i) {
                        // Found a newline
                        N++;
                        if (bgzf_write(fb, &gl->data[i], N-i) != N-i)
                            goto err;
                    }

                    // Flush bgzf block
                    int b_off = fb->block_offset;
                    if (!eod && b_off &&
                        ((char *)fb->uncompressed_block)[b_off-1] == '\n')
                        if (bgzf_flush_try(fb, BGZF_BLOCK_SIZE) < 0)
                            goto err;

                    // Copy from N onwards into next block
                    if (i+R > N)
                        if (bgzf_write(fb, &gl->data[N], i+R - N)
                            != i+R - N)
                            goto err;

                    i = i+R;
                }
            } else {
                if (hwrite(fp->fp.hfile, gl->data, gl->data_size) != gl->data_size)
                    goto err;
            }
        }

        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&fd->lines_m);
        gl->next = fd->lines;
        fd->lines = gl;
        pthread_mutex_unlock(&fd->lines_m);
    }

    sam_state_err(fd, 0); // success
    hts_tpool_process_shutdown(fd->q);
    return NULL;

 err:
    sam_state_err(fd, errno ? errno : EIO);
    return (void *)-1;
}

// Run from one of the worker threads.
// Convert a passed in array of BAMs (sp_bams) and converts to a block
// of text SAM records (sp_lines).
static void *sam_format_worker(void *arg) {
    sp_bams *gb = (sp_bams *)arg;
    sp_lines *gl = NULL;
    int i;
    SAM_state *fd = gb->fd;
    htsFile *fp = fd->fp;

    // Use a block of SAM strings we had earlier if available.
    pthread_mutex_lock(&fd->lines_m);
    if (fd->lines) {
        gl = fd->lines;
        fd->lines = gl->next;
    }
    pthread_mutex_unlock(&fd->lines_m);

    if (gl == NULL) {
        gl = calloc(1, sizeof(*gl));
        if (!gl) {
            sam_state_err(fd, ENOMEM);
            return NULL;
        }
        gl->alloc = gl->data_size = 0;
        gl->data = NULL;
    }
    gl->serial = gb->serial;
    gl->next = NULL;

    kstring_t ks = {0, gl->alloc, gl->data};

    for (i = 0; i < gb->nbams; i++) {
        if (sam_format1_append(fd->h, &gb->bams[i], &ks) < 0) {
            sam_state_err(fd, errno ? errno : EIO);
            goto err;
        }
        kputc('\n', &ks);
    }

    pthread_mutex_lock(&fd->lines_m);
    gl->data_size = ks.l;
    gl->alloc = ks.m;
    gl->data = ks.s;

    if (fp->idx) {
        // Keep hold of the bam array a little longer as
        // sam_dispatcher_write needs to use them for building the index.
        gl->bams = gb;
    } else {
        // Add bam array to free-list
        gb->next = fd->bams;
        fd->bams = gb;
    }
    pthread_mutex_unlock(&fd->lines_m);

    return gl;

 err:
    // Possible race between this and fd->curr_bam.
    // Easier to not free and leave it on the input list so it
    // gets freed there instead?
    // sam_free_sp_bams(gb);
    if (gl) {
        free(gl->data);
        free(gl);
    }
    return NULL;
}

int sam_set_thread_pool(htsFile *fp, htsThreadPool *p) {
    if (fp->state)
        return 0;

    if (!(fp->state = sam_state_create(fp)))
        return -1;
    SAM_state *fd = (SAM_state *)fp->state;

    pthread_mutex_init(&fd->lines_m, NULL);
    pthread_mutex_init(&fd->command_m, NULL);
    pthread_cond_init(&fd->command_c, NULL);
    fd->p = p->pool;
    int qsize = p->qsize;
    if (!qsize)
        qsize = 2*hts_tpool_size(fd->p);
    fd->q = hts_tpool_process_init(fd->p, qsize, 0);
    if (!fd->q) {
        sam_state_destroy(fp);
        return -1;
    }

    if (fp->format.compression == bgzf)
        return bgzf_thread_pool(fp->fp.bgzf, p->pool, p->qsize);

    return 0;
}

int sam_set_threads(htsFile *fp, int nthreads) {
    if (nthreads <= 0)
        return 0;

    htsThreadPool p;
    p.pool = hts_tpool_init(nthreads);
    p.qsize = nthreads*2;

    int ret = sam_set_thread_pool(fp, &p);
    if (ret < 0)
        return ret;

    SAM_state *fd = (SAM_state *)fp->state;
    fd->own_pool = 1;

    return 0;
}

typedef struct {
    kstring_t name;
    kstring_t comment; // NB: pointer into name, do not free
    kstring_t seq;
    kstring_t qual;
    int casava;
    int aux;
    int rnum;
    char BC[3];         // aux tag ID for barcode
    khash_t(tag) *tags; // which aux tags to use (if empty, use all).
    char nprefix;
    int sra_names;
} fastq_state;

// Initialise fastq state.
// Name char of '@' or '>' distinguishes fastq vs fasta variant
static fastq_state *fastq_state_init(int name_char) {
    fastq_state *x = (fastq_state *)calloc(1, sizeof(*x));
    if (!x)
        return NULL;
    strcpy(x->BC, "BC");
    x->nprefix = name_char;

    return x;
}

void fastq_state_destroy(htsFile *fp) {
    if (fp->state) {
        fastq_state *x = (fastq_state *)fp->state;
        if (x->tags)
            kh_destroy(tag, x->tags);
        ks_free(&x->name);
        ks_free(&x->seq);
        ks_free(&x->qual);
        free(fp->state);
    }
}

int fastq_state_set(samFile *fp, enum hts_fmt_option opt, ...) {
    va_list args;

    if (!fp)
        return -1;
    if (!fp->state)
        if (!(fp->state = fastq_state_init(fp->format.format == fastq_format
                                           ? '@' : '>')))
            return -1;

    fastq_state *x = (fastq_state *)fp->state;

    switch (opt) {
    case FASTQ_OPT_CASAVA:
        x->casava = 1;
        break;

    case FASTQ_OPT_NAME2:
        x->sra_names = 1;
        break;

    case FASTQ_OPT_AUX: {
        va_start(args, opt);
        x->aux = 1;
        char *tag = va_arg(args, char *);
        va_end(args);
        if (tag && strcmp(tag, "1") != 0) {
            if (!x->tags)
                if (!(x->tags = kh_init(tag)))
                    return -1;

            size_t i, tlen = strlen(tag);
            for (i = 0; i+3 <= tlen+1; i += 3) {
                if (tag[i+0] == ',' || tag[i+1] == ',' ||
                    !(tag[i+2] == ',' || tag[i+2] == '\0')) {
                    hts_log_warning("Bad tag format '%.3s'; skipping option", tag+i);
                    break;
                }
                int ret, tcode = tag[i+0]*256 + tag[i+1];
                kh_put(tag, x->tags, tcode, &ret);
                if (ret < 0)
                    return -1;
            }
        }
        break;
    }

    case FASTQ_OPT_BARCODE: {
        va_start(args, opt);
        char *bc = va_arg(args, char *);
        va_end(args);
        strncpy(x->BC, bc, 2);
        x->BC[2] = 0;
        break;
    }

    case FASTQ_OPT_RNUM:
        x->rnum = 1;
        break;

    default:
        break;
    }
    return 0;
}

static int fastq_parse1(htsFile *fp, bam1_t *b) {
    fastq_state *x = (fastq_state *)fp->state;
    size_t i, l;
    int ret = 0;

    if (fp->format.format == fasta_format && fp->line.s) {
        // For FASTA we've already read the >name line; steal it
        // Not the most efficient, but we don't optimise for fasta reading.
        if (fp->line.l == 0)
            return -1; // EOF

        free(x->name.s);
        x->name = fp->line;
        fp->line.l = fp->line.m = 0;
        fp->line.s = NULL;
    } else {
        // Read a FASTQ format entry.
        ret = hts_getline(fp, KS_SEP_LINE, &x->name);
        if (ret == -1)
            return -1;  // EOF
        else if (ret < -1)
            return ret; // ERR
    }

    // Name
    if (*x->name.s != x->nprefix)
        return -2;

    // Reverse the SRA strangeness of putting the run_name.number before
    // the read name.
    i = 0;
    char *name = x->name.s+1;
    if (x->sra_names) {
        char *cp = strpbrk(x->name.s, " \t");
        if (cp) {
            while (*cp == ' ' || *cp == '\t')
                cp++;
            *--cp = '@';
            i = cp - x->name.s;
            name = cp+1;
        }
    }

    l = x->name.l;
    char *s = x->name.s;
    while (i < l && !isspace_c(s[i]))
        i++;
    if (i < l) {
        s[i] = 0;
        x->name.l = i++;
    }

    // Comment; a kstring struct, but pointer into name line.  (Do not free)
    while (i < l && isspace_c(s[i]))
        i++;
    x->comment.s = s+i;
    x->comment.l = l - i;

    // Seq
    x->seq.l = 0;
    for (;;) {
        if ((ret = hts_getline(fp, KS_SEP_LINE, &fp->line)) < 0)
            if (fp->format.format == fastq_format || ret < -1)
                return -2;
        if (ret == -1 ||
            *fp->line.s == (fp->format.format == fastq_format ? '+' : '>'))
            break;
        if (kputsn(fp->line.s, fp->line.l, &x->seq) < 0)
            return -2;
    }

    // Qual
    if (fp->format.format == fastq_format) {
        size_t remainder = x->seq.l;
        x->qual.l = 0;
        do {
            if (hts_getline(fp, KS_SEP_LINE, &fp->line) < 0)
                return -2;
            if (fp->line.l > remainder)
                return -2;
            if (kputsn(fp->line.s, fp->line.l, &x->qual) < 0)
                return -2;
            remainder -= fp->line.l;
        } while (remainder > 0);

        // Decr qual
        for (i = 0; i < x->qual.l; i++)
            x->qual.s[i] -= '!';
    }

    int flag = BAM_FUNMAP; int pflag = BAM_FMUNMAP | BAM_FPAIRED;
    if (x->name.l > 2 &&
        x->name.s[x->name.l-2] == '/' &&
        isdigit_c(x->name.s[x->name.l-1])) {
        switch(x->name.s[x->name.l-1]) {
        case '1': flag |= BAM_FREAD1 | pflag; break;
        case '2': flag |= BAM_FREAD2 | pflag; break;
        default : flag |= BAM_FREAD1 | BAM_FREAD2 | pflag; break;
        }
        x->name.s[x->name.l-=2] = 0;
    }

    // Convert to BAM
    ret = bam_set1(b,
                   x->name.s + x->name.l - name, name,
                   flag,
                   -1, -1, 0, // ref '*', pos, mapq,
                   0, NULL,     // no cigar,
                   -1, -1, 0,    // mate
                   x->seq.l, x->seq.s, x->qual.s,
                   0);

    // Identify Illumina CASAVA strings.
    // <read>:<is_filtered>:<control_bits>:<barcode_sequence>
    char *barcode = NULL;
    int barcode_len = 0;
    kstring_t *kc = &x->comment;
    char *endptr;
    if (x->casava &&
        // \d:[YN]:\d+:[ACGTN]+
        kc->l > 6 && (kc->s[1] | kc->s[3]) == ':' && isdigit_c(kc->s[0]) &&
        strtol(kc->s+4, &endptr, 10) >= 0 && endptr != kc->s+4
        && *endptr == ':') {

        // read num
        switch(kc->s[0]) {
        case '1': b->core.flag |= BAM_FREAD1 | pflag; break;
        case '2': b->core.flag |= BAM_FREAD2 | pflag; break;
        default : b->core.flag |= BAM_FREAD1 | BAM_FREAD2 | pflag; break;
        }

        if (kc->s[2] == 'Y')
            b->core.flag |= BAM_FQCFAIL;

        // Barcode, maybe numeric in which case we skip it
        if (!isdigit_c(endptr[1])) {
            barcode = endptr+1;
            for (i = barcode - kc->s; i < kc->l; i++)
                if (isspace_c(kc->s[i]))
                    break;

            kc->s[i] = 0;
            barcode_len = i+1-(barcode - kc->s);
        }
    }

    if (ret >= 0 && barcode_len)
        if (bam_aux_append(b, x->BC, 'Z', barcode_len, (uint8_t *)barcode) < 0)
            ret = -2;

    if (!x->aux)
        return ret;

    // Identify any SAM style aux tags in comments too.
    if (aux_parse(&kc->s[barcode_len], kc->s + kc->l, b, 1, x->tags) < 0)
        ret = -2;

    return ret;
}

// Internal component of sam_read1 below
static inline int sam_read1_bam(htsFile *fp, sam_hdr_t *h, bam1_t *b) {
    int ret = bam_read1(fp->fp.bgzf, b);
    if (h && ret >= 0) {
        if (b->core.tid  >= h->n_targets || b->core.tid  < -1 ||
            b->core.mtid >= h->n_targets || b->core.mtid < -1) {
            errno = ERANGE;
            return -3;
        }
    }
    return ret;
}

// Internal component of sam_read1 below
static inline int sam_read1_cram(htsFile *fp, sam_hdr_t *h, bam1_t **b) {
    int ret = cram_get_bam_seq(fp->fp.cram, b);
    if (ret < 0)
        return cram_eof(fp->fp.cram) ? -1 : -2;

    if (bam_tag2cigar(*b, 1, 1) < 0)
        return -2;

    return ret;
}

// Internal component of sam_read1 below
static inline int sam_read1_sam(htsFile *fp, sam_hdr_t *h, bam1_t *b) {
    int ret;

    // Consume 1st line after header parsing as it wasn't using peek
    if (fp->line.l != 0) {
        ret = sam_parse1(&fp->line, h, b);
        fp->line.l = 0;
        return ret;
    }

    if (fp->state) {
        SAM_state *fd = (SAM_state *)fp->state;

        if (fp->format.compression == bgzf && fp->fp.bgzf->seeked) {
            // We don't support multi-threaded SAM parsing with seeks yet.
            int ret;
            if ((ret = sam_state_destroy(fp)) < 0) {
                errno = -ret;
                return -2;
            }
            if (bgzf_seek(fp->fp.bgzf, fp->fp.bgzf->seeked, SEEK_SET) < 0)
                return -1;
            fp->fp.bgzf->seeked = 0;
            goto err_recover;
        }

        if (!fd->h) {
            fd->h = h;
            fd->h->ref_count++;
            // Ensure hrecs is initialised now as we don't want multiple
            // threads trying to do this simultaneously.
            if (!fd->h->hrecs && sam_hdr_fill_hrecs(fd->h) < 0)
                return -2;

            // We can only do this once we've got a header
            if (pthread_create(&fd->dispatcher, NULL, sam_dispatcher_read,
                               fp) != 0)
                return -2;
            fd->dispatcher_set = 1;
        }

        if (fd->h != h) {
            hts_log_error("SAM multi-threaded decoding does not support changing header");
            return -1;
        }

        sp_bams *gb = fd->curr_bam;
        if (!gb) {
            if (fd->errcode) {
                // In case reader failed
                errno = fd->errcode;
                return -2;
            }
            hts_tpool_result *r = hts_tpool_next_result_wait(fd->q);
            if (!r)
                return -2;
            fd->curr_bam = gb = (sp_bams *)hts_tpool_result_data(r);
            hts_tpool_delete_result(r, 0);
        }
        if (!gb)
            return fd->errcode ? -2 : -1;
        bam1_t *b_array = (bam1_t *)gb->bams;
        if (fd->curr_idx < gb->nbams)
            if (!bam_copy1(b, &b_array[fd->curr_idx++]))
                return -2;
        if (fd->curr_idx == gb->nbams) {
            pthread_mutex_lock(&fd->lines_m);
            gb->next = fd->bams;
            fd->bams = gb;
            pthread_mutex_unlock(&fd->lines_m);

            fd->curr_bam = NULL;
            fd->curr_idx = 0;
        }

        ret = 0;

    } else  {
    err_recover:
        ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
        if (ret < 0) return ret;

        ret = sam_parse1(&fp->line, h, b);
        fp->line.l = 0;
        if (ret < 0) {
            hts_log_warning("Parse error at line %lld", (long long)fp->lineno);
            if (h && h->ignore_sam_err) goto err_recover;
        }
    }

    return ret;
}

// Returns 0 on success,
//        -1 on EOF,
//       <-1 on error
int sam_read1(htsFile *fp, sam_hdr_t *h, bam1_t *b)
{
    int ret, pass_filter;

    do {
        switch (fp->format.format) {
        case bam:
            ret = sam_read1_bam(fp, h, b);
            break;

        case cram:
            ret = sam_read1_cram(fp, h, &b);
            break;

        case sam:
            ret = sam_read1_sam(fp, h, b);
            break;

        case fasta_format:
        case fastq_format: {
            fastq_state *x = (fastq_state *)fp->state;
            if (!x) {
                if (!(fp->state = fastq_state_init(fp->format.format
                                                   == fastq_format ? '@' : '>')))
                    return -2;
            }

            return fastq_parse1(fp, b);
        }

        case empty_format:
            errno = EPIPE;
            return -3;

        default:
            errno = EFTYPE;
            return -3;
        }

        pass_filter = (ret >= 0 && fp->filter)
            ? sam_passes_filter(h, b, fp->filter)
            : 1;
    } while (pass_filter == 0);

    return pass_filter < 0 ? -2 : ret;
}


static int sam_format1_append(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    int i, r = 0;
    uint8_t *s, *end;
    const bam1_core_t *c = &b->core;

    if (c->l_qname == 0)
        return -1;
    r |= kputsn_(bam_get_qname(b), c->l_qname-1-c->l_extranul, str);
    r |= kputc_('\t', str); // query name
    r |= kputw(c->flag, str); r |= kputc_('\t', str); // flag
    if (c->tid >= 0) { // chr
        r |= kputs(h->target_name[c->tid] , str);
        r |= kputc_('\t', str);
    } else r |= kputsn_("*\t", 2, str);
    r |= kputll(c->pos + 1, str); r |= kputc_('\t', str); // pos
    r |= kputw(c->qual, str); r |= kputc_('\t', str); // qual
    if (c->n_cigar) { // cigar
        uint32_t *cigar = bam_get_cigar(b);
        for (i = 0; i < c->n_cigar; ++i) {
            r |= kputw(bam_cigar_oplen(cigar[i]), str);
            r |= kputc_(bam_cigar_opchr(cigar[i]), str);
        }
    } else r |= kputc_('*', str);
    r |= kputc_('\t', str);
    if (c->mtid < 0) r |= kputsn_("*\t", 2, str); // mate chr
    else if (c->mtid == c->tid) r |= kputsn_("=\t", 2, str);
    else {
        r |= kputs(h->target_name[c->mtid], str);
        r |= kputc_('\t', str);
    }
    r |= kputll(c->mpos + 1, str); r |= kputc_('\t', str); // mate pos
    r |= kputll(c->isize, str); r |= kputc_('\t', str); // template len
    if (c->l_qseq) { // seq and qual
        uint8_t *s = bam_get_seq(b);
        if (ks_resize(str, str->l+2+2*c->l_qseq) < 0) goto mem_err;
        char *cp = str->s + str->l;

        // Sequence, 2 bases at a time
        nibble2base(s, cp, c->l_qseq);
        cp[c->l_qseq] = '\t';
        cp += c->l_qseq+1;

        // Quality
        s = bam_get_qual(b);
        i = 0;
        if (s[0] == 0xff) {
            cp[i++] = '*';
        } else {
            // local copy of c->l_qseq to aid unrolling
            uint32_t lqseq = c->l_qseq;
            for (i = 0; i < lqseq; ++i)
                cp[i]=s[i]+33;
        }
        cp[i] = 0;
        cp += i;
        str->l = cp - str->s;
    } else r |= kputsn_("*\t*", 3, str);

    s = bam_get_aux(b); // aux
    end = b->data + b->l_data;

    while (end - s >= 4) {
        r |= kputc_('\t', str);
        if ((s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)) == NULL)
            goto bad_aux;
    }
    r |= kputsn("", 0, str); // nul terminate
    if (r < 0) goto mem_err;

    return str->l;

 bad_aux:
    hts_log_error("Corrupted aux data for read %.*s",
                  b->core.l_qname, bam_get_qname(b));
    errno = EINVAL;
    return -1;

 mem_err:
    hts_log_error("Out of memory");
    errno = ENOMEM;
    return -1;
}

int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
    str->l = 0;
    return sam_format1_append(h, b, str);
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end);
int fastq_format1(fastq_state *x, const bam1_t *b, kstring_t *str)
{
    unsigned flag = b->core.flag;
    int i, e = 0, len = b->core.l_qseq;
    uint8_t *seq, *qual;

    str->l = 0;

    if (len == 0) return 0;

    // Name
    if (kputc(x->nprefix, str) == EOF || kputs(bam_get_qname(b), str) == EOF)
        return -1;

    // /1 or /2 suffix
    if (x && x->rnum && (flag & BAM_FPAIRED)) {
        int r12 = flag & (BAM_FREAD1 | BAM_FREAD2);
        if (r12 == BAM_FREAD1) {
            if (kputs("/1", str) == EOF)
                return -1;
        } else if (r12 == BAM_FREAD2) {
            if (kputs("/2", str) == EOF)
                return -1;
        }
    }

    // Illumina CASAVA tag.
    // This is <rnum>:<Y/N qcfail>:<control-bits>:<barcode-or-zero>
    if (x && x->casava) {
        int rnum = (flag & BAM_FREAD1)? 1 : (flag & BAM_FREAD2)? 2 : 0;
        char filtered = (flag & BAM_FQCFAIL)? 'Y' : 'N';
        uint8_t *bc = bam_aux_get(b, x->BC);
        if (ksprintf(str, " %d:%c:0:%s", rnum, filtered,
                     bc ? (char *)bc+1 : "0") < 0)
            return -1;

        if (bc && (*bc != 'Z' || (!isupper_c(bc[1]) && !islower_c(bc[1])))) {
            hts_log_warning("BC tag starts with non-sequence base; using '0'");
            str->l -= strlen((char *)bc)-2; // limit to 1 char
            str->s[str->l-1] = '0';
            str->s[str->l] = 0;
            bc = NULL;
        }

        // Replace any non-alpha with '+'.  Ie seq-seq to seq+seq
        if (bc) {
            int l = strlen((char *)bc+1);
            char *c = (char *)str->s + str->l - l;
            for (i = 0; i < l; i++) {
                if (!isalpha_c(c[i]))
                    c[i] = '+';
                else if (islower_c(c[i]))
                    c[i] = toupper_c(c[i]);
            }
        }
    }

    // Aux tags
    if (x && x->aux) {
        uint8_t *s = bam_get_aux(b), *end = b->data + b->l_data;
        while (s && end - s >= 4) {
            int tt = s[0]*256 + s[1];
            if (x->tags == NULL ||
                kh_get(tag, x->tags, tt) != kh_end(x->tags)) {
                e |= kputc_('\t', str) < 0;
                if (!(s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)))
                    return -1;
            } else {
                s = skip_aux(s+2, end);
            }
        }
        e |= kputsn("", 0, str) < 0; // nul terminate
    }

    if (ks_resize(str, str->l + 1 + len+1 + 2 + len+1 + 1) < 0) return -1;
    e |= kputc_('\n', str) < 0;

    // Seq line
    seq = bam_get_seq(b);
    if (flag & BAM_FREVERSE)
        for (i = len-1; i >= 0; i--)
            e |= kputc_("!TGKCYSBAWRDMHVN"[bam_seqi(seq, i)], str) < 0;
    else
        for (i = 0; i < len; i++)
            e |= kputc_(seq_nt16_str[bam_seqi(seq, i)], str) < 0;


    // Qual line
    if (x->nprefix == '@') {
        kputsn("\n+\n", 3, str);
        qual = bam_get_qual(b);
        if (qual[0] == 0xff)
            for (i = 0; i < len; i++)
                e |= kputc_('B', str) < 0;
        else if (flag & BAM_FREVERSE)
            for (i = len-1; i >= 0; i--)
                e |= kputc_(33 + qual[i], str) < 0;
        else
            for (i = 0; i < len; i++)
                e |= kputc_(33 + qual[i], str) < 0;

    }
    e |= kputc('\n', str) < 0;

    return e ? -1 : str->l;
}

// Sadly we need to be able to modify the bam_hdr here so we can
// reference count the structure.
int sam_write1(htsFile *fp, const sam_hdr_t *h, const bam1_t *b)
{
    switch (fp->format.format) {
    case binary_format:
        fp->format.category = sequence_data;
        fp->format.format = bam;
        /* fall-through */
    case bam:
        return bam_write_idx1(fp, h, b);

    case cram:
        return cram_put_bam_seq(fp->fp.cram, (bam1_t *)b);

    case text_format:
        fp->format.category = sequence_data;
        fp->format.format = sam;
        /* fall-through */
    case sam:
        if (fp->state) {
            SAM_state *fd = (SAM_state *)fp->state;

            // Threaded output
            if (!fd->h) {
                // NB: discard const.  We don't actually modify sam_hdr_t here,
                // just data pointed to by it (which is a bit weasely still),
                // but out cached pointer must be non-const as we want to
                // destroy it later on and sam_hdr_destroy takes non-const.
                //
                // We do this because some tools do sam_hdr_destroy; sam_close
                // while others do sam_close; sam_hdr_destroy.  The former is
                // an issue as we need the header still when flushing.
                fd->h = (sam_hdr_t *)h;
                fd->h->ref_count++;

                if (pthread_create(&fd->dispatcher, NULL, sam_dispatcher_write,
                                   fp) != 0)
                    return -2;
                fd->dispatcher_set = 1;
            }

            if (fd->h != h) {
                hts_log_error("SAM multi-threaded decoding does not support changing header");
                return -2;
            }

            // Find a suitable BAM array to copy to
            sp_bams *gb = fd->curr_bam;
            if (!gb) {
                pthread_mutex_lock(&fd->lines_m);
                if (fd->bams) {
                    fd->curr_bam = gb = fd->bams;
                    fd->bams = gb->next;
                    gb->next = NULL;
                    gb->nbams = 0;
                    gb->bam_mem = 0;
                    pthread_mutex_unlock(&fd->lines_m);
                } else {
                    pthread_mutex_unlock(&fd->lines_m);
                    if (!(gb = calloc(1, sizeof(*gb)))) return -1;
                    if (!(gb->bams = calloc(SAM_NBAM, sizeof(*gb->bams)))) {
                        free(gb);
                        return -1;
                    }
                    gb->nbams = 0;
                    gb->abams = SAM_NBAM;
                    gb->bam_mem = 0;
                    gb->fd = fd;
                    fd->curr_idx = 0;
                    fd->curr_bam = gb;
                }
            }

            if (!bam_copy1(&gb->bams[gb->nbams++], b))
                return -2;
            gb->bam_mem += b->l_data + sizeof(*b);

            // Dispatch if full
            if (gb->nbams == SAM_NBAM || gb->bam_mem > SAM_NBYTES*0.8) {
                gb->serial = fd->serial++;
                pthread_mutex_lock(&fd->command_m);
                if (fd->errcode != 0) {
                    pthread_mutex_unlock(&fd->command_m);
                    return -fd->errcode;
                }
                if (hts_tpool_dispatch3(fd->p, fd->q, sam_format_worker, gb,
                                        cleanup_sp_bams,
                                        cleanup_sp_lines, 0) < 0) {
                    pthread_mutex_unlock(&fd->command_m);
                    return -1;
                }
                pthread_mutex_unlock(&fd->command_m);
                fd->curr_bam = NULL;
            }

            // Dummy value as we don't know how long it really is.
            // We could track file sizes via a SAM_state field, but I don't think
            // it is necessary.
            return 1;
        } else {
            if (sam_format1(h, b, &fp->line) < 0) return -1;
            kputc('\n', &fp->line);
            if (fp->is_bgzf) {
                if (bgzf_flush_try(fp->fp.bgzf, fp->line.l) < 0)
                    return -1;
                if ( bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l) != fp->line.l ) return -1;
            } else {
                if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;
            }

            if (fp->idx) {
                if (fp->format.compression == bgzf) {
                    if (bgzf_idx_push(fp->fp.bgzf, fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                      bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
                        return -1;
                    }
                } else {
                    if (hts_idx_push(fp->idx, b->core.tid, b->core.pos, bam_endpos(b),
                                     bgzf_tell(fp->fp.bgzf), !(b->core.flag&BAM_FUNMAP)) < 0) {
                        hts_log_error("Read '%s' with ref_name='%s', ref_length=%"PRIhts_pos", flags=%d, pos=%"PRIhts_pos" cannot be indexed",
                                bam_get_qname(b), sam_hdr_tid2name(h, b->core.tid), sam_hdr_tid2len(h, b->core.tid), b->core.flag, b->core.pos+1);
                        return -1;
                    }
                }
            }

            return fp->line.l;
        }


    case fasta_format:
    case fastq_format: {
        fastq_state *x = (fastq_state *)fp->state;
        if (!x) {
            if (!(fp->state = fastq_state_init(fp->format.format
                                               == fastq_format ? '@' : '>')))
                return -2;
        }

        if (fastq_format1(fp->state, b, &fp->line) < 0)
            return -1;
        if (fp->is_bgzf) {
            if (bgzf_flush_try(fp->fp.bgzf, fp->line.l) < 0)
                return -1;
            if (bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l) != fp->line.l)
                return -1;
        } else {
            if (hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l)
                return -1;
        }
        return fp->line.l;
    }

    default:
        errno = EBADF;
        return -1;
    }
}

/************************
 *** Auxiliary fields ***
 ************************/
#ifndef HTS_LITTLE_ENDIAN
static int aux_to_le(char type, uint8_t *out, const uint8_t *in, size_t len) {
    int tsz = aux_type2size(type);

    if (tsz >= 2 && tsz <= 8 && (len & (tsz - 1)) != 0) return -1;

    switch (tsz) {
        case 'H': case 'Z': case 1:  // Trivial
            memcpy(out, in, len);
            break;

#define aux_val_to_le(type_t, store_le) do {                            \
        type_t v;                                                       \
        size_t i;                                                       \
        for (i = 0; i < len; i += sizeof(type_t), out += sizeof(type_t)) { \
            memcpy(&v, in + i, sizeof(type_t));                         \
            store_le(v, out);                                           \
        }                                                               \
    } while (0)

        case 2: aux_val_to_le(uint16_t, u16_to_le); break;
        case 4: aux_val_to_le(uint32_t, u32_to_le); break;
        case 8: aux_val_to_le(uint64_t, u64_to_le); break;

#undef aux_val_to_le

        case 'B': { // Recurse!
            uint32_t n;
            if (len < 5) return -1;
            memcpy(&n, in + 1, 4);
            out[0] = in[0];
            u32_to_le(n, out + 1);
            return aux_to_le(in[0], out + 5, in + 5, len - 5);
        }

        default: // Unknown type code
            return -1;
    }



    return 0;
}
#endif

int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data)
{
    uint32_t new_len;

    assert(b->l_data >= 0);
    new_len = b->l_data + 3 + len;
    if (new_len > INT32_MAX || new_len < b->l_data) goto nomem;

    if (realloc_bam_data(b, new_len) < 0) return -1;

    b->data[b->l_data] = tag[0];
    b->data[b->l_data + 1] = tag[1];
    b->data[b->l_data + 2] = type;

#ifdef HTS_LITTLE_ENDIAN
    memcpy(b->data + b->l_data + 3, data, len);
#else
    if (aux_to_le(type, b->data + b->l_data + 3, data, len) != 0) {
        errno = EINVAL;
        return -1;
    }
#endif

    b->l_data = new_len;

    return 0;

 nomem:
    errno = ENOMEM;
    return -1;
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

uint8_t *bam_aux_first(const bam1_t *b)
{
    uint8_t *s = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    if (s >= end) { errno = ENOENT; return NULL; }
    return s+2;
}

uint8_t *bam_aux_next(const bam1_t *b, const uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = s? skip_aux((uint8_t *) s, end) : end;
    if (next == NULL) goto bad_aux;
    if (next >= end) { errno = ENOENT; return NULL; }
    return next+2;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s", bam_get_qname(b));
    errno = EINVAL;
    return NULL;
}

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    for (s = bam_aux_first(b); s; s = bam_aux_next(b, s))
        if (s[-2] == tag[0] && s[-1] == tag[1]) {
            // Check the tag value is valid and complete
            uint8_t *e = skip_aux(s, b->data + b->l_data);
            if (e == NULL) goto bad_aux;
            if ((*s == 'Z' || *s == 'H') && *(e - 1) != '\0') goto bad_aux;

            return s;
        }

    // errno now as set by bam_aux_first()/bam_aux_next()
    return NULL;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s", bam_get_qname(b));
    errno = EINVAL;
    return NULL;
}

int bam_aux_del(bam1_t *b, uint8_t *s)
{
    s = bam_aux_remove(b, s);
    return (s || errno == ENOENT)? 0 : -1;
}

uint8_t *bam_aux_remove(bam1_t *b, uint8_t *s)
{
    uint8_t *end = b->data + b->l_data;
    uint8_t *next = skip_aux(s, end);
    if (next == NULL) goto bad_aux;

    b->l_data -= next - (s-2);
    if (next >= end) { errno = ENOENT; return NULL; }

    memmove(s-2, next, end - next);
    return s;

 bad_aux:
    hts_log_error("Corrupted aux data for read %s", bam_get_qname(b));
    errno = EINVAL;
    return NULL;
}

int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data)
{
    // FIXME: This is not at all efficient!
    size_t ln = len >= 0 ? len : strlen(data) + 1;
    size_t old_ln = 0;
    int need_nul = ln == 0 || data[ln - 1] != '\0';
    int save_errno = errno;
    int new_tag = 0;
    uint8_t *s = bam_aux_get(b,tag), *e;

    if (s) {  // Replacing existing tag
        char type = *s;
        if (type != 'Z') {
            hts_log_error("Called bam_aux_update_str for type '%c' instead of 'Z'", type);
            errno = EINVAL;
            return -1;
        }
        s++;
        e = memchr(s, '\0', b->data + b->l_data - s);
        old_ln = (e ? e - s : b->data + b->l_data - s) + 1;
        s -= 3;
    } else {
        if (errno != ENOENT) { // Invalid aux data, give up
            return -1;
        } else { // Tag doesn't exist - put it on the end
            errno = save_errno;
            s = b->data + b->l_data;
            new_tag = 3;
        }
    }

    if (old_ln < ln + need_nul + new_tag) {
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, ln + need_nul + new_tag - old_ln) < 0)
            return -1;
        s = b->data + s_offset;
    }
    if (!new_tag) {
        memmove(s + 3 + ln + need_nul,
                s + 3 + old_ln,
                b->l_data - (s + 3 - b->data) - old_ln);
    }
    b->l_data += new_tag + ln + need_nul - old_ln;

    s[0] = tag[0];
    s[1] = tag[1];
    s[2] = 'Z';
    memmove(s+3,data,ln);
    if (need_nul) s[3 + ln] = '\0';
    return 0;
}

int bam_aux_update_int(bam1_t *b, const char tag[2], int64_t val)
{
    uint32_t sz, old_sz = 0, new = 0;
    uint8_t *s, type;

    if (val < INT32_MIN || val > UINT32_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    if (val < INT16_MIN)       { type = 'i'; sz = 4; }
    else if (val < INT8_MIN)   { type = 's'; sz = 2; }
    else if (val < 0)          { type = 'c'; sz = 1; }
    else if (val < UINT8_MAX)  { type = 'C'; sz = 1; }
    else if (val < UINT16_MAX) { type = 'S'; sz = 2; }
    else                       { type = 'I'; sz = 4; }

    s = bam_aux_get(b, tag);
    if (s) {  // Tag present - how big was the old one?
        switch (*s) {
            case 'c': case 'C': old_sz = 1; break;
            case 's': case 'S': old_sz = 2; break;
            case 'i': case 'I': old_sz = 4; break;
            default: errno = EINVAL; return -1;  // Not an integer
        }
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            s = b->data + b->l_data;
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    if (new || old_sz < sz) {
        // Make room for new tag
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, (new ? 3 : 0) + sz - old_sz) < 0)
            return -1;
        s =  b->data + s_offset;
        if (new) { // Add tag id
            *s++ = tag[0];
            *s++ = tag[1];
        } else {   // Shift following data so we have space
            memmove(s + sz, s + old_sz, b->l_data - s_offset - old_sz);
        }
    } else {
        // Reuse old space.  Data value may be bigger than necessary but
        // we avoid having to move everything else
        sz = old_sz;
        type = (val < 0 ? "\0cs\0i" : "\0CS\0I")[old_sz];
        assert(type > 0);
    }
    *s++ = type;
#ifdef HTS_LITTLE_ENDIAN
    memcpy(s, &val, sz);
#else
    switch (sz) {
        case 4:  u32_to_le(val, s); break;
        case 2:  u16_to_le(val, s); break;
        default: *s = val; break;
    }
#endif
    b->l_data += (new ? 3 : 0) + sz - old_sz;
    return 0;
}

int bam_aux_update_float(bam1_t *b, const char tag[2], float val)
{
    uint8_t *s = bam_aux_get(b, tag);
    int shrink = 0, new = 0;

    if (s) { // Tag present - what was it?
        switch (*s) {
            case 'f': break;
            case 'd': shrink = 1; break;
            default: errno = EINVAL; return -1;  // Not a float
        }
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    if (new) { // Ensure there's room
        if (possibly_expand_bam_data(b, 3 + 4) < 0)
            return -1;
        s = b->data + b->l_data;
        *s++ = tag[0];
        *s++ = tag[1];
    } else if (shrink) { // Convert non-standard double tag to float
        memmove(s + 5, s + 9, b->l_data - ((s + 9) - b->data));
        b->l_data -= 4;
    }
    *s++ = 'f';
    float_to_le(val, s);
    if (new) b->l_data += 7;

    return 0;
}

int bam_aux_update_array(bam1_t *b, const char tag[2],
                         uint8_t type, uint32_t items, void *data)
{
    uint8_t *s = bam_aux_get(b, tag);
    size_t old_sz = 0, new_sz;
    int new = 0;

    if (s) { // Tag present
        if (*s != 'B') { errno = EINVAL; return -1; }
        old_sz = aux_type2size(s[1]);
        if (old_sz < 1 || old_sz > 4) { errno = EINVAL; return -1; }
        old_sz *= le_to_u32(s + 2);
    } else {
        if (errno == ENOENT) {  // Tag doesn't exist - add a new one
            s = b->data + b->l_data;
            new = 1;
        }  else { // Invalid aux data, give up.
            return -1;
        }
    }

    new_sz = aux_type2size(type);
    if (new_sz < 1 || new_sz > 4) { errno = EINVAL; return -1; }
    if (items > INT32_MAX / new_sz) { errno = ENOMEM; return -1; }
    new_sz *= items;

    if (new || old_sz < new_sz) {
        // Make room for new tag
        ptrdiff_t s_offset = s - b->data;
        if (possibly_expand_bam_data(b, (new ? 8 : 0) + new_sz - old_sz) < 0)
            return -1;
        s =  b->data + s_offset;
    }
    if (new) { // Add tag id and type
        *s++ = tag[0];
        *s++ = tag[1];
        *s = 'B';
        b->l_data += 8 + new_sz;
    } else if (old_sz != new_sz) { // shift following data if necessary
        memmove(s + 6 + new_sz, s + 6 + old_sz,
                b->l_data - ((s + 6 + old_sz) - b->data));
        b->l_data -= old_sz;
        b->l_data += new_sz;
    }

    s[1] = type;
    u32_to_le(items, s + 2);
#ifdef HTS_LITTLE_ENDIAN
    memcpy(s + 6, data, new_sz);
    return 0;
#else
    return aux_to_le(type, s + 6, data, new_sz);
#endif
}

static inline int64_t get_int_aux_val(uint8_t type, const uint8_t *s,
                                      uint32_t idx)
{
    switch (type) {
        case 'c': return le_to_i8(s + idx);
        case 'C': return s[idx];
        case 's': return le_to_i16(s + 2 * idx);
        case 'S': return le_to_u16(s + 2 * idx);
        case 'i': return le_to_i32(s + 4 * idx);
        case 'I': return le_to_u32(s + 4 * idx);
        default:
            errno = EINVAL;
            return 0;
    }
}

int64_t bam_aux2i(const uint8_t *s)
{
    int type;
    type = *s++;
    return get_int_aux_val(type, s, 0);
}

double bam_aux2f(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'd') return le_to_double(s);
    else if (type == 'f') return le_to_float(s);
    else return get_int_aux_val(type, s, 0);
}

char bam_aux2A(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'A') return *(char*)s;
    errno = EINVAL;
    return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
    int type;
    type = *s++;
    if (type == 'Z' || type == 'H') return (char*)s;
    errno = EINVAL;
    return 0;
}

uint32_t bam_auxB_len(const uint8_t *s)
{
    if (s[0] != 'B') {
        errno = EINVAL;
        return 0;
    }
    return le_to_u32(s + 2);
}

int64_t bam_auxB2i(const uint8_t *s, uint32_t idx)
{
    uint32_t len = bam_auxB_len(s);
    if (idx >= len) {
        errno = ERANGE;
        return 0;
    }
    return get_int_aux_val(s[1], s + 6, idx);
}

double bam_auxB2f(const uint8_t *s, uint32_t idx)
{
    uint32_t len = bam_auxB_len(s);
    if (idx >= len) {
        errno = ERANGE;
        return 0.0;
    }
    if (s[1] == 'f') return le_to_float(s + 6 + 4 * idx);
    else return get_int_aux_val(s[1], s + 6, idx);
}

int sam_open_mode(char *mode, const char *fn, const char *format)
{
    // TODO Parse "bam5" etc for compression level
    if (format == NULL) {
        // Try to pick a format based on the filename extension
        char extension[HTS_MAX_EXT_LEN];
        if (find_file_extension(fn, extension) < 0) return -1;
        return sam_open_mode(mode, fn, extension);
    }
    else if (strcasecmp(format, "bam") == 0) strcpy(mode, "b");
    else if (strcasecmp(format, "cram") == 0) strcpy(mode, "c");
    else if (strcasecmp(format, "sam") == 0) strcpy(mode, "");
    else if (strcasecmp(format, "sam.gz") == 0) strcpy(mode, "z");
    else if (strcasecmp(format, "fastq") == 0 ||
             strcasecmp(format, "fq") == 0) strcpy(mode, "f");
    else if (strcasecmp(format, "fastq.gz") == 0 ||
             strcasecmp(format, "fq.gz") == 0) strcpy(mode, "fz");
    else if (strcasecmp(format, "fasta") == 0 ||
             strcasecmp(format, "fa") == 0) strcpy(mode, "F");
    else if (strcasecmp(format, "fasta.gz") == 0 ||
             strcasecmp(format, "fa.gz") == 0) strcpy(mode, "Fz");
    else return -1;

    return 0;
}

// A version of sam_open_mode that can handle ,key=value options.
// The format string is allocated and returned, to be freed by the caller.
// Prefix should be "r" or "w",
char *sam_open_mode_opts(const char *fn,
                         const char *mode,
                         const char *format)
{
    char *mode_opts = malloc((format ? strlen(format) : 1) +
                             (mode   ? strlen(mode)   : 1) + 12);
    char *opts, *cp;
    int format_len;

    if (!mode_opts)
        return NULL;

    strcpy(mode_opts, mode ? mode : "r");
    cp = mode_opts + strlen(mode_opts);

    if (format == NULL) {
        // Try to pick a format based on the filename extension
        char extension[HTS_MAX_EXT_LEN];
        if (find_file_extension(fn, extension) < 0) {
            free(mode_opts);
            return NULL;
        }
        if (sam_open_mode(cp, fn, extension) == 0) {
            return mode_opts;
        } else {
            free(mode_opts);
            return NULL;
        }
    }

    if ((opts = strchr(format, ','))) {
        format_len = opts-format;
    } else {
        opts="";
        format_len = strlen(format);
    }

    if (strncmp(format, "bam", format_len) == 0) {
        *cp++ = 'b';
    } else if (strncmp(format, "cram", format_len) == 0) {
        *cp++ = 'c';
    } else if (strncmp(format, "cram2", format_len) == 0) {
        *cp++ = 'c';
        strcpy(cp, ",VERSION=2.1");
        cp += 12;
    } else if (strncmp(format, "cram3", format_len) == 0) {
        *cp++ = 'c';
        strcpy(cp, ",VERSION=3.0");
        cp += 12;
    } else if (strncmp(format, "sam", format_len) == 0) {
        ; // format mode=""
    } else if (strncmp(format, "sam.gz", format_len) == 0) {
        *cp++ = 'z';
    } else if (strncmp(format, "fastq", format_len) == 0 ||
               strncmp(format, "fq", format_len) == 0) {
        *cp++ = 'f';
    } else if (strncmp(format, "fastq.gz", format_len) == 0 ||
               strncmp(format, "fq.gz", format_len) == 0) {
        *cp++ = 'f';
        *cp++ = 'z';
    } else if (strncmp(format, "fasta", format_len) == 0 ||
               strncmp(format, "fa", format_len) == 0) {
        *cp++ = 'F';
    } else if (strncmp(format, "fasta.gz", format_len) == 0 ||
               strncmp(format, "fa", format_len) == 0) {
        *cp++ = 'F';
        *cp++ = 'z';
    } else {
        free(mode_opts);
        return NULL;
    }

    strcpy(cp, opts);

    return mode_opts;
}

#define STRNCMP(a,b,n) (strncasecmp((a),(b),(n)) || strlen(a)!=(n))
int bam_str2flag(const char *str)
{
    char *end, *beg = (char*) str;
    long int flag = strtol(str, &end, 0);
    if ( end!=str ) return flag;    // the conversion was successful
    flag = 0;
    while ( *str )
    {
        end = beg;
        while ( *end && *end!=',' ) end++;
        if ( !STRNCMP("PAIRED",beg,end-beg) ) flag |= BAM_FPAIRED;
        else if ( !STRNCMP("PROPER_PAIR",beg,end-beg) ) flag |= BAM_FPROPER_PAIR;
        else if ( !STRNCMP("UNMAP",beg,end-beg) ) flag |= BAM_FUNMAP;
        else if ( !STRNCMP("MUNMAP",beg,end-beg) ) flag |= BAM_FMUNMAP;
        else if ( !STRNCMP("REVERSE",beg,end-beg) ) flag |= BAM_FREVERSE;
        else if ( !STRNCMP("MREVERSE",beg,end-beg) ) flag |= BAM_FMREVERSE;
        else if ( !STRNCMP("READ1",beg,end-beg) ) flag |= BAM_FREAD1;
        else if ( !STRNCMP("READ2",beg,end-beg) ) flag |= BAM_FREAD2;
        else if ( !STRNCMP("SECONDARY",beg,end-beg) ) flag |= BAM_FSECONDARY;
        else if ( !STRNCMP("QCFAIL",beg,end-beg) ) flag |= BAM_FQCFAIL;
        else if ( !STRNCMP("DUP",beg,end-beg) ) flag |= BAM_FDUP;
        else if ( !STRNCMP("SUPPLEMENTARY",beg,end-beg) ) flag |= BAM_FSUPPLEMENTARY;
        else return -1;
        if ( !*end ) break;
        beg = end + 1;
    }
    return flag;
}

char *bam_flag2str(int flag)
{
    kstring_t str = {0,0,0};
    if ( flag&BAM_FPAIRED ) ksprintf(&str,"%s%s", str.l?",":"","PAIRED");
    if ( flag&BAM_FPROPER_PAIR ) ksprintf(&str,"%s%s", str.l?",":"","PROPER_PAIR");
    if ( flag&BAM_FUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","UNMAP");
    if ( flag&BAM_FMUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","MUNMAP");
    if ( flag&BAM_FREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","REVERSE");
    if ( flag&BAM_FMREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","MREVERSE");
    if ( flag&BAM_FREAD1 ) ksprintf(&str,"%s%s", str.l?",":"","READ1");
    if ( flag&BAM_FREAD2 ) ksprintf(&str,"%s%s", str.l?",":"","READ2");
    if ( flag&BAM_FSECONDARY ) ksprintf(&str,"%s%s", str.l?",":"","SECONDARY");
    if ( flag&BAM_FQCFAIL ) ksprintf(&str,"%s%s", str.l?",":"","QCFAIL");
    if ( flag&BAM_FDUP ) ksprintf(&str,"%s%s", str.l?",":"","DUP");
    if ( flag&BAM_FSUPPLEMENTARY ) ksprintf(&str,"%s%s", str.l?",":"","SUPPLEMENTARY");
    if ( str.l == 0 ) kputsn("", 0, &str);
    return str.s;
}


/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

#include <assert.h>

/*******************
 *** Memory pool ***
 *******************/

typedef struct {
    int k, y;
    hts_pos_t x, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
    bam1_t b;
    hts_pos_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
    bam_pileup_cd cd;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init(void)
{
    mempool_t *mp;
    mp = (mempool_t*)calloc(1, sizeof(mempool_t));
    return mp;
}
static void mp_destroy(mempool_t *mp)
{
    int k;
    for (k = 0; k < mp->n; ++k) {
        free(mp->buf[k]->b.data);
        free(mp->buf[k]);
    }
    free(mp->buf);
    free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
    else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
    --mp->cnt; p->next = 0; // clear lbnode_t::next here
    if (mp->n == mp->max) {
        mp->max = mp->max? mp->max<<1 : 256;
        mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
    }
    mp->buf[mp->n++] = p;
}

/**********************
 *** CIGAR resolver ***
 **********************/

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordinate of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, hts_pos_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

    bam1_t *b = p->b;
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    int k;
    // determine the current CIGAR operation
    //fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam_get_qname(b), pos, s->end, s->k, s->x, s->y);
    if (s->k == -1) { // never processed
        p->qpos = 0;
        if (c->n_cigar == 1) { // just one operation, save a loop
          if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
        } else { // find the first match or deletion
            for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
                int op = _cop(cigar[k]);
                int l = _cln(cigar[k]);
                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                    op == BAM_CEQUAL || op == BAM_CDIFF) break;
                else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
            }
            assert(k < c->n_cigar);
            s->k = k;
        }
    } else { // the read has been processed before
        int op, l = _cln(cigar[s->k]);
        if (pos - s->x >= l) { // jump to the next operation
            assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
            op = _cop(cigar[s->k+1]);
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                ++s->k;
            } else { // find the next M/D/N/=/X
              if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
                s->x += l;
                for (k = s->k + 1; k < c->n_cigar; ++k) {
                    op = _cop(cigar[k]), l = _cln(cigar[k]);
                    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
                    else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
                }
                s->k = k;
            }
            assert(s->k < c->n_cigar); // otherwise a bug
        } // else, do nothing
    }
    { // collect pileup information
        int op, l;
        op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
        p->is_del = p->indel = p->is_refskip = 0;
        if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
            int op2 = _cop(cigar[s->k+1]);
            int l2 = _cln(cigar[s->k+1]);
            if (op2 == BAM_CDEL && op != BAM_CDEL) {
                // At start of a new deletion, merge e.g. 1D2D to 3D.
                // Within a deletion (the 2D in 1D2D) we keep p->indel=0
                // and rely on is_del=1 as we would for 3D.
                p->indel = -(int)l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CDEL) p->indel -= l2;
                    else break;
                }
            } else if (op2 == BAM_CINS) {
                p->indel = l2;
                for (k = s->k+2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) p->indel += l2;
                    else if (op2 != BAM_CPAD) break;
                }
            } else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) {
                int l3 = 0;
                for (k = s->k + 2; k < c->n_cigar; ++k) {
                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
                    if (op2 == BAM_CINS) l3 += l2;
                    else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
                }
                if (l3 > 0) p->indel = l3;
            }
        }
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            p->qpos = s->y + (pos - s->x);
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
            p->is_refskip = (op == BAM_CREF_SKIP);
        } // cannot be other operations; otherwise a bug
        p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
    }
    p->cigar_ind = s->k;
    return 1;
}

/*******************************
 *** Expansion of insertions ***
 *******************************/

/*
 * Fills out the kstring with the padded insertion sequence for the current
 * location in 'p'.  If this is not an insertion site, the string is blank.
 *
 * This variant handles base modifications, but only when "m" is non-NULL.
 *
 * Returns the number of inserted base on success, with string length being
 *        accessable via ins->l;
 *        -1 on failure.
 */
int bam_plp_insertion_mod(const bam_pileup1_t *p,
                          hts_base_mod_state *m,
                          kstring_t *ins, int *del_len) {
    int j, k, indel, nb = 0;
    uint32_t *cigar;

    if (p->indel <= 0) {
        if (ks_resize(ins, 1) < 0)
            return -1;
        ins->l = 0;
        ins->s[0] = '\0';
        return 0;
    }

    if (del_len)
        *del_len = 0;

    // Measure indel length including pads
    indel = 0;
    k = p->cigar_ind+1;
    cigar = bam_get_cigar(p->b);
    while (k < p->b->core.n_cigar) {
        switch (cigar[k] & BAM_CIGAR_MASK) {
        case BAM_CPAD:
        case BAM_CINS:
            indel += (cigar[k] >> BAM_CIGAR_SHIFT);
            break;
        default:
            k = p->b->core.n_cigar;
            break;
        }
        k++;
    }
    nb = ins->l = indel;

    // Produce sequence
    if (ks_resize(ins, indel+1) < 0)
        return -1;
    indel = 0;
    k = p->cigar_ind+1;
    j = 1;
    while (k < p->b->core.n_cigar) {
        int l, c;
        switch (cigar[k] & BAM_CIGAR_MASK) {
        case BAM_CPAD:
            for (l = 0; l < (cigar[k]>>BAM_CIGAR_SHIFT); l++)
                ins->s[indel++] = '*';
            break;
        case BAM_CINS:
            for (l = 0; l < (cigar[k]>>BAM_CIGAR_SHIFT); l++, j++) {
                c = p->qpos + j - p->is_del < p->b->core.l_qseq
                    ? seq_nt16_str[bam_seqi(bam_get_seq(p->b),
                                            p->qpos + j - p->is_del)]
                    : 'N';
                ins->s[indel++] = c;
                int nm;
                hts_base_mod mod[256];
                if (m && (nm = bam_mods_at_qpos(p->b, p->qpos + j - p->is_del,
                                                m, mod, 256)) > 0) {
                    int o_indel = indel;
                    if (ks_resize(ins, ins->l + nm*16+3) < 0)
                        return -1;
                    ins->s[indel++] = '[';
                    int j;
                    for (j = 0; j < nm; j++) {
                        char qual[20];
                        if (mod[j].qual >= 0)
                            sprintf(qual, "%d", mod[j].qual);
                        else
                            *qual=0;
                        if (mod[j].modified_base < 0)
                            // ChEBI
                            indel += sprintf(&ins->s[indel], "%c(%d)%s",
                                             "+-"[mod[j].strand],
                                             -mod[j].modified_base,
                                             qual);
                        else
                            indel += sprintf(&ins->s[indel], "%c%c%s",
                                             "+-"[mod[j].strand],
                                             mod[j].modified_base,
                                             qual);
                    }
                    ins->s[indel++] = ']';
                    ins->l += indel - o_indel; // grow by amount we used
                }
            }
            break;
        case BAM_CDEL:
            // eg cigar 1M2I1D gives mpileup output in T+2AA-1C style
            if (del_len)
                *del_len = cigar[k]>>BAM_CIGAR_SHIFT;
            // fall through
        default:
            k = p->b->core.n_cigar;
            break;
        }
        k++;
    }
    ins->s[indel] = '\0';
    ins->l = indel; // string length

    return nb;      // base length
}

/*
 * Fills out the kstring with the padded insertion sequence for the current
 * location in 'p'.  If this is not an insertion site, the string is blank.
 *
 * This is the original interface with no capability for reporting base
 * modifications.
 *
 * Returns the length of insertion string on success;
 *        -1 on failure.
 */
int bam_plp_insertion(const bam_pileup1_t *p, kstring_t *ins, int *del_len) {
    return bam_plp_insertion_mod(p, NULL, ins, del_len);
}

/***********************
 *** Pileup iterator ***
 ***********************/

// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct bam_plp_s {
    mempool_t *mp;
    lbnode_t *head, *tail;
    int32_t tid, max_tid;
    hts_pos_t pos, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;

    // For notification of creation and destruction events
    // and associated client-owned pointer.
    int (*plp_construct)(void *data, const bam1_t *b, bam_pileup_cd *cd);
    int (*plp_destruct )(void *data, const bam1_t *b, bam_pileup_cd *cd);
};

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
{
    bam_plp_t iter;
    iter = (bam_plp_t)calloc(1, sizeof(struct bam_plp_s));
    iter->mp = mp_init();
    iter->head = iter->tail = mp_alloc(iter->mp);
    iter->max_tid = iter->max_pos = -1;
    iter->maxcnt = 8000;
    if (func) {
        iter->func = func;
        iter->data = data;
        iter->b = bam_init1();
    }
    return iter;
}

int bam_plp_init_overlaps(bam_plp_t iter)
{
    iter->overlaps = kh_init(olap_hash);  // hash for tweaking quality of bases in overlapping reads
    return iter->overlaps ? 0 : -1;
}

void bam_plp_destroy(bam_plp_t iter)
{
    lbnode_t *p, *pnext;
    if ( iter->overlaps ) kh_destroy(olap_hash, iter->overlaps);
    for (p = iter->head; p != NULL; p = pnext) {
        pnext = p->next;
        mp_free(iter->mp, p);
    }
    mp_destroy(iter->mp);
    if (iter->b) bam_destroy1(iter->b);
    free(iter->plp);
    free(iter);
}

void bam_plp_constructor(bam_plp_t plp,
                         int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd)) {
    plp->plp_construct = func;
}

void bam_plp_destructor(bam_plp_t plp,
                        int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd)) {
    plp->plp_destruct = func;
}

//---------------------------------
//---  Tweak overlapping reads
//---------------------------------

/**
 *  cigar_iref2iseq_set()  - find the first CMATCH setting the ref and the read index
 *  cigar_iref2iseq_next() - get the next CMATCH base
 *  @cigar:       pointer to current cigar block (rw)
 *  @cigar_max:   pointer just beyond the last cigar block
 *  @icig:        position within the current cigar block (rw)
 *  @iseq:        position in the sequence (rw)
 *  @iref:        position with respect to the beginning of the read (iref_pos - b->core.pos) (rw)
 *
 *  Returns BAM_CMATCH, -1 when there is no more cigar to process or the requested position is not covered,
 *  or -2 on error.
 */
static inline int cigar_iref2iseq_set(const uint32_t **cigar,
                                      const uint32_t *cigar_max,
                                      hts_pos_t *icig,
                                      hts_pos_t *iseq,
                                      hts_pos_t *iref)
{
    hts_pos_t pos = *iref;
    if ( pos < 0 ) return -1;
    *icig = 0;
    *iseq = 0;
    *iref = 0;
    while ( *cigar<cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            pos -= ncig;
            if ( pos < 0 ) { *icig = ncig + pos; *iseq += *icig; *iref += *icig; return BAM_CMATCH; }
            (*cigar)++; *iseq += ncig; *icig = 0; *iref += ncig;
            continue;
        }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP )
        {
            pos -= ncig;
            if ( pos<0 ) pos = 0;
            (*cigar)++; *icig = 0; *iref += ncig;
            continue;
        }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    return -1;
}
static inline int cigar_iref2iseq_next(const uint32_t **cigar,
                                       const uint32_t *cigar_max,
                                       hts_pos_t *icig,
                                       hts_pos_t *iseq,
                                       hts_pos_t *iref)
{
    while ( *cigar < cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            if ( *icig >= ncig - 1 ) { *icig = -1;  (*cigar)++; continue; }
            (*iseq)++; (*icig)++; (*iref)++;
            return BAM_CMATCH;
        }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP ) { (*cigar)++; (*iref) += ncig; *icig = -1; continue; }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = -1; continue; }
        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = -1; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = -1; continue; }
        hts_log_error("Unexpected cigar %d", cig);
        return -2;
    }
    *iseq = -1;
    *iref = -1;
    return -1;
}

// Given overlapping read 'a' (left) and 'b' (right) on the same
// template, adjust quality values to zero for either a or b.
// Note versions 1.12 and earlier always removed quality from 'b' for
// matching bases.  Now we select a or b semi-randomly based on name hash.
// Returns 0 on success,
//        -1 on failure
static int tweak_overlap_quality(bam1_t *a, bam1_t *b)
{
    const uint32_t *a_cigar = bam_get_cigar(a),
        *a_cigar_max = a_cigar + a->core.n_cigar;
    const uint32_t *b_cigar = bam_get_cigar(b),
        *b_cigar_max = b_cigar + b->core.n_cigar;
    hts_pos_t a_icig = 0, a_iseq = 0;
    hts_pos_t b_icig = 0, b_iseq = 0;
    uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
    uint8_t *a_seq  = bam_get_seq(a), *b_seq = bam_get_seq(b);

    hts_pos_t iref   = b->core.pos;
    hts_pos_t a_iref = iref - a->core.pos;
    hts_pos_t b_iref = iref - b->core.pos;

    int a_ret = cigar_iref2iseq_set(&a_cigar, a_cigar_max,
                                    &a_icig, &a_iseq, &a_iref);
    if ( a_ret<0 )
        // no overlap or error
        return a_ret<-1 ? -1:0;

    int b_ret = cigar_iref2iseq_set(&b_cigar, b_cigar_max,
                                    &b_icig, &b_iseq, &b_iref);
    if ( b_ret<0 )
        // no overlap or error
        return b_ret<-1 ? -1:0;

    // Determine which seq is the one getting modified qualities.
    uint8_t amul, bmul;
    if (__ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(a))) & 1) {
        amul = 1;
        bmul = 0;
    } else {
        amul = 0;
        bmul = 1;
    }

    // Loop over the overlapping region nulling qualities in either
    // seq a or b.
    int err = 0;
    while ( 1 )
    {
        // Step to next matching reference position in a and b
        while ( a_ret >= 0 && a_iref>=0 && a_iref < iref - a->core.pos )
            a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max,
                                         &a_icig, &a_iseq, &a_iref);
        if ( a_ret<0 ) { // done
            err = a_ret<-1?-1:0;
            break;
        }
        if ( iref < a_iref + a->core.pos )
            iref = a_iref + a->core.pos;

        while ( b_ret >= 0 && b_iref>=0 && b_iref < iref - b->core.pos )
            b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig,
                                         &b_iseq, &b_iref);
        if ( b_ret<0 ) { // done
            err = b_ret<-1?-1:0;
            break;
        }
        if ( iref < b_iref + b->core.pos )
            iref = b_iref + b->core.pos;

        iref++;

        if ( a_iref+a->core.pos != b_iref+b->core.pos )
            // only CMATCH positions, don't know what to do with indels
            continue;

        if (a_iseq > a->core.l_qseq || b_iseq > b->core.l_qseq)
            // Fell off end of sequence, bad CIGAR?
            return -1;

        // We're finally at the same ref base in both a and b.
        // Check if the bases match (confident) or mismatch
        // (not so confident).
        if ( bam_seqi(a_seq,a_iseq) == bam_seqi(b_seq,b_iseq) ) {
            // We are very confident about this base.  Use sum of quals
            int qual = a_qual[a_iseq] + b_qual[b_iseq];
            a_qual[a_iseq] = amul * (qual>200 ? 200 : qual);
            b_qual[b_iseq] = bmul * (qual>200 ? 200 : qual);;
        } else {
            // Not so confident about anymore given the mismatch.
            // Reduce qual for lowest quality base.
            if ( a_qual[a_iseq] > b_qual[b_iseq] ) {
                // A highest qual base; keep
                a_qual[a_iseq] = 0.8 * a_qual[a_iseq];
                b_qual[b_iseq] = 0;
            } else if (a_qual[a_iseq] < b_qual[b_iseq] ) {
                // B highest qual base; keep
                b_qual[b_iseq] = 0.8 * b_qual[b_iseq];
                a_qual[a_iseq] = 0;
            } else {
                // Both equal, so pick randomly
                a_qual[a_iseq] = amul * 0.8 * a_qual[a_iseq];
                b_qual[b_iseq] = bmul * 0.8 * b_qual[b_iseq];
            }
        }
    }

    return err;
}

// Fix overlapping reads. Simple soft-clipping did not give good results.
// Lowering qualities of unwanted bases is more selective and works better.
//
// Returns 0 on success, -1 on failure
static int overlap_push(bam_plp_t iter, lbnode_t *node)
{
    if ( !iter->overlaps ) return 0;

    // mapped mates and paired reads only
    if ( node->b.core.flag&BAM_FMUNMAP || !(node->b.core.flag&BAM_FPROPER_PAIR) ) return 0;

    // no overlap possible, unless some wild cigar
    if ( (node->b.core.mtid >= 0 && node->b.core.tid != node->b.core.mtid)
         || (llabs(node->b.core.isize) >= 2*node->b.core.l_qseq
         && node->b.core.mpos >= node->end) // for those wild cigars
       ) return 0;

    khiter_t kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(&node->b));
    if ( kitr==kh_end(iter->overlaps) )
    {
        // Only add reads where the mate is still to arrive
        if (node->b.core.mpos >= node->b.core.pos ||
            ((node->b.core.flag & BAM_FPAIRED) && node->b.core.mpos == -1)) {
            int ret;
            kitr = kh_put(olap_hash, iter->overlaps, bam_get_qname(&node->b), &ret);
            if (ret < 0) return -1;
            kh_value(iter->overlaps, kitr) = node;
        }
    }
    else
    {
        lbnode_t *a = kh_value(iter->overlaps, kitr);
        int err = tweak_overlap_quality(&a->b, &node->b);
        kh_del(olap_hash, iter->overlaps, kitr);
        assert(a->end-1 == a->s.end);
        return err;
    }
    return 0;
}

static void overlap_remove(bam_plp_t iter, const bam1_t *b)
{
    if ( !iter->overlaps ) return;

    khiter_t kitr;
    if ( b )
    {
        kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(b));
        if ( kitr!=kh_end(iter->overlaps) )
            kh_del(olap_hash, iter->overlaps, kitr);
    }
    else
    {
        // remove all
        for (kitr = kh_begin(iter->overlaps); kitr<kh_end(iter->overlaps); kitr++)
            if ( kh_exist(iter->overlaps, kitr) ) kh_del(olap_hash, iter->overlaps, kitr);
    }
}



// Prepares next pileup position in bam records collected by bam_plp_auto -> user func -> bam_plp_push. Returns
// pointer to the piled records if next position is ready or NULL if there is not enough records in the
// buffer yet (the current position is still the maximum position across all buffered reads).
const bam_pileup1_t *bam_plp64_next(bam_plp_t iter, int *_tid, hts_pos_t *_pos, int *_n_plp)
{
    if (iter->error) { *_n_plp = -1; return NULL; }
    *_n_plp = 0;
    if (iter->is_eof && iter->head == iter->tail) return NULL;
    while (iter->is_eof || iter->max_tid > iter->tid || (iter->max_tid == iter->tid && iter->max_pos > iter->pos)) {
        int n_plp = 0;
        // write iter->plp at iter->pos
        lbnode_t **pptr = &iter->head;
        while (*pptr != iter->tail) {
            lbnode_t *p = *pptr;
            if (p->b.core.tid < iter->tid || (p->b.core.tid == iter->tid && p->end <= iter->pos)) { // then remove
                overlap_remove(iter, &p->b);
                if (iter->plp_destruct)
                    iter->plp_destruct(iter->data, &p->b, &p->cd);
                *pptr = p->next; mp_free(iter->mp, p);
            }
            else {
                if (p->b.core.tid == iter->tid && p->beg <= iter->pos) { // here: p->end > pos; then add to pileup
                    if (n_plp == iter->max_plp) { // then double the capacity
                        iter->max_plp = iter->max_plp? iter->max_plp<<1 : 256;
                        iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
                    }
                    iter->plp[n_plp].b = &p->b;
                    iter->plp[n_plp].cd = p->cd;
                    if (resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s)) ++n_plp; // actually always true...
                }
                pptr = &(*pptr)->next;
            }
        }
        *_n_plp = n_plp; *_tid = iter->tid; *_pos = iter->pos;
        // update iter->tid and iter->pos
        if (iter->head != iter->tail) {
            if (iter->tid > iter->head->b.core.tid) {
                hts_log_error("Unsorted input. Pileup aborts");
                iter->error = 1;
                *_n_plp = -1;
                return NULL;
            }
        }
        if (iter->tid < iter->head->b.core.tid) { // come to a new reference sequence
            iter->tid = iter->head->b.core.tid; iter->pos = iter->head->beg; // jump to the next reference
        } else if (iter->pos < iter->head->beg) { // here: tid == head->b.core.tid
            iter->pos = iter->head->beg; // jump to the next position
        } else ++iter->pos; // scan contiguously
        // return
        if (n_plp) return iter->plp;
        if (iter->is_eof && iter->head == iter->tail) break;
    }
    return NULL;
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
    hts_pos_t pos64 = 0;
    const bam_pileup1_t *p = bam_plp64_next(iter, _tid, &pos64, _n_plp);
    if (pos64 < INT_MAX) {
        *_pos = pos64;
    } else {
        hts_log_error("Position %"PRId64" too large", pos64);
        *_pos = INT_MAX;
        iter->error = 1;
        *_n_plp = -1;
        return NULL;
    }
    return p;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b)
{
    if (iter->error) return -1;
    if (b) {
        if (b->core.tid < 0) { overlap_remove(iter, b); return 0; }
        // Skip only unmapped reads here, any additional filtering must be done in iter->func
        if (b->core.flag & BAM_FUNMAP) { overlap_remove(iter, b); return 0; }
        if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt)
        {
            overlap_remove(iter, b);
            return 0;
        }
        if (bam_copy1(&iter->tail->b, b) == NULL)
            return -1;
        iter->tail->b.id = iter->id++;
        iter->tail->beg = b->core.pos;
        // Use raw rlen rather than bam_endpos() which adjusts rlen=0 to rlen=1
        iter->tail->end = b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        iter->tail->s = g_cstate_null; iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
        if (b->core.tid < iter->max_tid) {
            hts_log_error("The input is not sorted (chromosomes out of order)");
            iter->error = 1;
            return -1;
        }
        if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
            hts_log_error("The input is not sorted (reads out of order)");
            iter->error = 1;
            return -1;
        }
        iter->max_tid = b->core.tid; iter->max_pos = iter->tail->beg;
        if (iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
            lbnode_t *next = mp_alloc(iter->mp);
            if (!next) {
                iter->error = 1;
                return -1;
            }
            if (iter->plp_construct) {
                if (iter->plp_construct(iter->data, &iter->tail->b,
                                        &iter->tail->cd) < 0) {
                    mp_free(iter->mp, next);
                    iter->error = 1;
                    return -1;
                }
            }
            if (overlap_push(iter, iter->tail) < 0) {
                mp_free(iter->mp, next);
                iter->error = 1;
                return -1;
            }
            iter->tail->next = next;
            iter->tail = iter->tail->next;
        }
    } else iter->is_eof = 1;
    return 0;
}

const bam_pileup1_t *bam_plp64_auto(bam_plp_t iter, int *_tid, hts_pos_t *_pos, int *_n_plp)
{
    const bam_pileup1_t *plp;
    if (iter->func == 0 || iter->error) { *_n_plp = -1; return 0; }
    if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
    else { // no pileup line can be obtained; read alignments
        *_n_plp = 0;
        if (iter->is_eof) return 0;
        int ret;
        while ( (ret=iter->func(iter->data, iter->b)) >= 0) {
            if (bam_plp_push(iter, iter->b) < 0) {
                *_n_plp = -1;
                return 0;
            }
            if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
            // otherwise no pileup line can be returned; read the next alignment.
        }
        if ( ret < -1 ) { iter->error = ret; *_n_plp = -1; return 0; }
        if (bam_plp_push(iter, 0) < 0) {
            *_n_plp = -1;
            return 0;
        }
        if ((plp = bam_plp64_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
        return 0;
    }
}

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
    hts_pos_t pos64 = 0;
    const bam_pileup1_t *p = bam_plp64_auto(iter, _tid, &pos64, _n_plp);
    if (pos64 < INT_MAX) {
        *_pos = pos64;
    } else {
        hts_log_error("Position %"PRId64" too large", pos64);
        *_pos = INT_MAX;
        iter->error = 1;
        *_n_plp = -1;
        return NULL;
    }
    return p;
}

void bam_plp_reset(bam_plp_t iter)
{
    overlap_remove(iter, NULL);
    iter->max_tid = iter->max_pos = -1;
    iter->tid = iter->pos = 0;
    iter->is_eof = 0;
    while (iter->head != iter->tail) {
        lbnode_t *p = iter->head;
        iter->head = p->next;
        mp_free(iter->mp, p);
    }
}

void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
{
    iter->maxcnt = maxcnt;
}

/************************
 *** Mpileup iterator ***
 ************************/

struct bam_mplp_s {
    int n;
    int32_t min_tid, *tid;
    hts_pos_t min_pos, *pos;
    bam_plp_t *iter;
    int *n_plp;
    const bam_pileup1_t **plp;
};

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
{
    int i;
    bam_mplp_t iter;
    iter = (bam_mplp_t)calloc(1, sizeof(struct bam_mplp_s));
    iter->pos = (hts_pos_t*)calloc(n, sizeof(hts_pos_t));
    iter->tid = (int32_t*)calloc(n, sizeof(int32_t));
    iter->n_plp = (int*)calloc(n, sizeof(int));
    iter->plp = (const bam_pileup1_t**)calloc(n, sizeof(bam_pileup1_t*));
    iter->iter = (bam_plp_t*)calloc(n, sizeof(bam_plp_t));
    iter->n = n;
    iter->min_pos = HTS_POS_MAX;
    iter->min_tid = (uint32_t)-1;
    for (i = 0; i < n; ++i) {
        iter->iter[i] = bam_plp_init(func, data[i]);
        iter->pos[i] = iter->min_pos;
        iter->tid[i] = iter->min_tid;
    }
    return iter;
}

int bam_mplp_init_overlaps(bam_mplp_t iter)
{
    int i, r = 0;
    for (i = 0; i < iter->n; ++i)
        r |= bam_plp_init_overlaps(iter->iter[i]);
    return r == 0 ? 0 : -1;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
{
    int i;
    for (i = 0; i < iter->n; ++i)
        iter->iter[i]->maxcnt = maxcnt;
}

void bam_mplp_destroy(bam_mplp_t iter)
{
    int i;
    for (i = 0; i < iter->n; ++i) bam_plp_destroy(iter->iter[i]);
    free(iter->iter); free(iter->pos); free(iter->tid);
    free(iter->n_plp); free(iter->plp);
    free(iter);
}

int bam_mplp64_auto(bam_mplp_t iter, int *_tid, hts_pos_t *_pos, int *n_plp, const bam_pileup1_t **plp)
{
    int i, ret = 0;
    hts_pos_t new_min_pos = HTS_POS_MAX;
    uint32_t new_min_tid = (uint32_t)-1;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
            int tid;
            hts_pos_t pos;
            iter->plp[i] = bam_plp64_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
            if ( iter->iter[i]->error ) return -1;
            if (iter->plp[i]) {
                iter->tid[i] = tid;
                iter->pos[i] = pos;
            } else {
                iter->tid[i] = 0;
                iter->pos[i] = 0;
            }
        }
        if (iter->plp[i]) {
            if (iter->tid[i] < new_min_tid) {
                new_min_tid = iter->tid[i];
                new_min_pos = iter->pos[i];
            } else if (iter->tid[i] == new_min_tid && iter->pos[i] < new_min_pos) {
                new_min_pos = iter->pos[i];
            }
        }
    }
    iter->min_pos = new_min_pos;
    iter->min_tid = new_min_tid;
    if (new_min_pos == HTS_POS_MAX) return 0;
    *_tid = new_min_tid; *_pos = new_min_pos;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
            n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
            ++ret;
        } else n_plp[i] = 0, plp[i] = 0;
    }
    return ret;
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
    hts_pos_t pos64 = 0;
    int ret = bam_mplp64_auto(iter, _tid, &pos64, n_plp, plp);
    if (ret >= 0) {
        if (pos64 < INT_MAX) {
            *_pos = pos64;
        } else {
            hts_log_error("Position %"PRId64" too large", pos64);
            *_pos = INT_MAX;
            return -1;
        }
    }
    return ret;
}

void bam_mplp_reset(bam_mplp_t iter)
{
    int i;
    iter->min_pos = HTS_POS_MAX;
    iter->min_tid = (uint32_t)-1;
    for (i = 0; i < iter->n; ++i) {
        bam_plp_reset(iter->iter[i]);
        iter->pos[i] = HTS_POS_MAX;
        iter->tid[i] = (uint32_t)-1;
        iter->n_plp[i] = 0;
        iter->plp[i] = NULL;
    }
}

void bam_mplp_constructor(bam_mplp_t iter,
                          int (*func)(void *arg, const bam1_t *b, bam_pileup_cd *cd)) {
    int i;
    for (i = 0; i < iter->n; ++i)
        bam_plp_constructor(iter->iter[i], func);
}

void bam_mplp_destructor(bam_mplp_t iter,
                         int (*func)(void *arg, const bam1_t *b, bam_pileup_cd *cd)) {
    int i;
    for (i = 0; i < iter->n; ++i)
        bam_plp_destructor(iter->iter[i], func);
}

#endif // ~!defined(BAM_NO_PILEUP)

// ---------------------------
// Base Modification retrieval
//
// These operate by recording state in an opaque type, allocated and freed
// via the functions below.
//
// Initially we call bam_parse_basemod to process the tags and record the
// modifications in the state structure, and then functions such as
// bam_next_basemod can iterate over this cached state.

/*
 * Base modification are stored in MM/Mm tags as <mod_list> defined as
 *
 * <mod_list>        ::= <mod_chain><mod_list> | ""
 * <mod_chain>       ::= <canonical_base><strand><mod-list><delta-list>
 *
 * <canonical_base>  ::= "A" | "C" | "G" | "T" | "N".
 *
 * <strand>          ::= "+" | "-".
 *
 * <mod-list>        ::= <simple-mod-list> | <ChEBI-code>
 * <simple-mod-list> ::= <simple-mod><simple-mod-list> | <simple-mod>
 * <ChEBI-code>      ::= <integer>
 * <simple-mod>      ::= <letter>
 *
 * <delta-list>      ::= "," <integer> <delta-list> | ";"
 *
 * We do not allocate additional memory other than the fixed size
 * state, thus we track up to 256 pointers to different locations
 * within the MM and ML tags.  Each pointer is for a distinct
 * modification code (simple or ChEBI), meaning some may point to the
 * same delta-list when multiple codes are combined together
 * (e.g. "C+mh,1,5,18,3;").  This is the MM[] array.
 *
 * Each numeric in the delta-list is tracked in MMcount[], counted
 * down until it hits zero in which case the next delta is fetched.
 *
 * ML array similarly holds the locations in the quality (ML) tag per
 * type, but these are interleaved so C+mhfc,10,15 will have 4 types
 * all pointing to the same delta position, but in ML we store
 * Q(m0)Q(h0)Q(f0)Q(c0) followed by Q(m1)Q(h1)Q(f1)Q(c1).  This ML
 * also has MLstride indicating how many positions along ML to jump
 * each time we consume a base. (4 in our above example, but usually 1
 * for the simple case).
 *
 * One complexity of the base modification system is that mods are
 * always stored in the original DNA orientation.  This is so that
 * tools that may reverse-complement a sequence (eg "samtools fastq -T
 * MM,ML") can pass through these modification tags irrespective of
 * whether they have any knowledge of their internal workings.
 *
 * Because we don't wish to allocate extra memory, we cannot simply
 * reverse the MM and ML tags.  Sadly this means we have to manage the
 * reverse complementing ourselves on-the-fly.
 * For reversed reads we start at the right end of MM and no longer
 * stop at the semicolon.  Instead we use MMend[] array to mark the
 * termination point.
 */
#define MAX_BASE_MOD 256
struct hts_base_mod_state {
    int type[MAX_BASE_MOD];     // char or minus-CHEBI
    int canonical[MAX_BASE_MOD];// canonical base, as seqi (1,2,4,8,15)
    char strand[MAX_BASE_MOD];  // strand of modification; + or -
    int MMcount[MAX_BASE_MOD];  // no. canonical bases left until next mod
    char *MM[MAX_BASE_MOD];     // next pos delta (string)
    char *MMend[MAX_BASE_MOD];  // end of pos-delta string
    uint8_t *ML[MAX_BASE_MOD];  // next qual
    int MLstride[MAX_BASE_MOD]; // bytes between quals for this type
    int implicit[MAX_BASE_MOD]; // treat unlisted positions as non-modified?
    int seq_pos;                // current position along sequence
    int nmods;                  // used array size (0 to MAX_BASE_MOD-1).
};

hts_base_mod_state *hts_base_mod_state_alloc(void) {
    return calloc(1, sizeof(hts_base_mod_state));
}

void hts_base_mod_state_free(hts_base_mod_state *state) {
    free(state);
}

/*
 * Count frequency of A, C, G, T and N canonical bases in the sequence
 */
static void seq_freq(const bam1_t *b, int freq[16]) {
    int i;

    memset(freq, 0, 16*sizeof(*freq));
    uint8_t *seq = bam_get_seq(b);
    for (i = 0; i < b->core.l_qseq; i++)
        freq[bam_seqi(seq, i)]++;
    freq[15] = b->core.l_qseq; // all bases count as N for base mods
}

//0123456789ABCDEF
//=ACMGRSVTWYHKDBN  aka seq_nt16_str[]
//=TGKCYSBAWRDMHVN  comp1ement of seq_nt16_str
//084C2A6E195D3B7F
static int seqi_rc[] = { 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 };

/*
 * Parse the MM and ML tags to populate the base mod state.
 * This structure will have been previously allocated via
 * hts_base_mod_state_alloc, but it does not need to be repeatedly
 * freed and allocated for each new bam record. (Although obviously
 * it requires a new call to this function.)
 *
 */
int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state) {
    // Read MM and ML tags
    uint8_t *mm = bam_aux_get(b, "MM");
    if (!mm) mm = bam_aux_get(b, "Mm");
    if (!mm)
        return 0;
    if (mm[0] != 'Z') {
        hts_log_error("MM tag is not of type Z");
        return -1;
    }

    uint8_t *ml = bam_aux_get(b, "ML");
    if (!ml) ml = bam_aux_get(b, "Ml");
    if (ml && (ml[0] != 'B' || ml[1] != 'C')) {
        hts_log_error("ML tag is not of type B,C");
        return -1;
    }
    uint8_t *ml_end = ml ? ml+6 + le_to_u32(ml+2) : NULL;
    if (ml) ml += 6;

    state->seq_pos = 0;

    // Aggregate freqs of ACGTN if reversed, to get final-delta (later)
    int freq[16];
    if (b->core.flag & BAM_FREVERSE)
        seq_freq(b, freq);

    char *cp = (char *)mm+1;
    int mod_num = 0;
    int implicit = 1;
    while (*cp) {
        for (; *cp; cp++) {
            // cp should be [ACGTNU][+-]([a-zA-Z]+|[0-9]+)[.?]?(,\d+)*;
            unsigned char btype = *cp++;

            if (btype != 'A' && btype != 'C' &&
                btype != 'G' && btype != 'T' &&
                btype != 'U' && btype != 'N')
                return -1;
            if (btype == 'U') btype = 'T';

            btype = seq_nt16_table[btype];

            // Strand
            if (*cp != '+' && *cp != '-')
                return -1; // malformed
            char strand = *cp++;

            // List of modification types
            char *ms = cp, *me; // mod code start and end
            char *cp_end = NULL;
            int chebi = 0;
            if (isdigit_c(*cp)) {
                chebi = strtol(cp, &cp_end, 10);
                cp = cp_end;
                ms = cp-1;
            } else {
                while (*cp && isalpha_c(*cp))
                    cp++;
                if (*cp == '\0')
                    return -1;
            }

            me = cp;

            // Optional explicit vs implicit marker
            if (*cp == '.') {
                // default is implicit = 1;
                cp++;
            } else if (*cp == '?') {
                implicit = 0;
                cp++;
            } else if (*cp != ',' && *cp != ';') {
                // parse error
                return -1;
            }

            long delta;
            int n = 0; // nth symbol in a multi-mod string
            int stride = me-ms;
            int ndelta = 0;

            if (b->core.flag & BAM_FREVERSE) {
                // We process the sequence in left to right order,
                // but delta is successive count of bases to skip
                // counting right to left.  This also means the number
                // of bases to skip at left edge is unrecorded (as it's
                // the remainder).
                //
                // To output mods in left to right, we step through the
                // MM list in reverse and need to identify the left-end
                // "remainder" delta.
                int total_seq = 0;
                for (;;) {
                    cp += (*cp == ',');
                    if (*cp == 0 || *cp == ';')
                        break;

                    delta = strtol(cp, &cp_end, 10);
                    if (cp_end == cp) {
                        hts_log_error("Hit end of MM tag. Missing semicolon?");
                        return -1;
                    }

                    cp = cp_end;
                    total_seq += delta+1;
                    ndelta++;
                }
                delta = freq[seqi_rc[btype]] - total_seq; // remainder
            } else {
                delta = *cp == ','
                    ? strtol(cp+1, &cp_end, 10)
                    : 0;
                if (!cp_end) {
                    // empty list
                    delta = INT_MAX;
                    cp_end = cp+1;
                }
            }
            // Now delta is first in list or computed remainder,
            // and cp_end is either start or end of the MM list.
            while (ms < me) {
                state->type     [mod_num] = chebi ? -chebi : *ms;
                state->strand   [mod_num] = (strand == '-');
                state->canonical[mod_num] = btype;
                state->MLstride [mod_num] = stride;
                state->implicit [mod_num] = implicit;

                if (delta < 0) {
                    hts_log_error("MM tag refers to bases beyond sequence "
                                  "length");
                    return -1;
                }
                state->MMcount  [mod_num] = delta;
                if (b->core.flag & BAM_FREVERSE) {
                    state->MM   [mod_num] = cp+1;
                    state->MMend[mod_num] = cp_end;
                    state->ML   [mod_num] = ml ? ml+n +(ndelta-1)*stride: NULL;
                } else {
                    state->MM   [mod_num] = cp_end;
                    state->MMend[mod_num] = NULL;
                    state->ML   [mod_num] = ml ? ml+n : NULL;
                }

                if (++mod_num >= MAX_BASE_MOD) {
                    hts_log_error("Too many base modification types");
                    return -1;
                }
                ms++; n++;
            }

            // Skip modification deltas
            if (ml) {
                if (b->core.flag & BAM_FREVERSE) {
                    ml += ndelta*stride;
                } else {
                    while (*cp && *cp != ';') {
                        if (*cp == ',')
                            ml+=stride;
                        cp++;
                    }
                }
                if (ml > ml_end) {
                    hts_log_error("Insufficient number of entries in ML tag");
                    return -1;
                }
            } else {
                // cp_end already known if FREVERSE
                if (cp_end && (b->core.flag & BAM_FREVERSE))
                    cp = cp_end;
                else
                    while (*cp && *cp != ';')
                        cp++;
            }
            if (!*cp) {
                hts_log_error("Hit end of MM tag. Missing semicolon?");
                return -1;
            }
        }
    }

    state->nmods = mod_num;

    return 0;
}

/*
 * Fills out mods[] with the base modifications found.
 * Returns the number found (0 if none), which may be more than
 * the size of n_mods if more were found than reported.
 * Returns <= -1 on error.
 *
 * This always marches left to right along sequence, irrespective of
 * reverse flag or modification strand.
 */
int bam_mods_at_next_pos(const bam1_t *b, hts_base_mod_state *state,
                         hts_base_mod *mods, int n_mods) {
    if (b->core.flag & BAM_FREVERSE) {
        if (state->seq_pos < 0)
            return -1;
    } else {
        if (state->seq_pos >= b->core.l_qseq)
            return -1;
    }

    int i, j, n = 0;
    unsigned char base = bam_seqi(bam_get_seq(b), state->seq_pos);
    state->seq_pos++;
    if (b->core.flag & BAM_FREVERSE)
        base = seqi_rc[base];

    for (i = 0; i < state->nmods; i++) {
        if (state->canonical[i] != base && state->canonical[i] != 15/*N*/)
            continue;

        if (state->MMcount[i]-- > 0)
            continue;

        char *MMptr = state->MM[i];
        if (n < n_mods) {
            mods[n].modified_base = state->type[i];
            mods[n].canonical_base = seq_nt16_str[state->canonical[i]];
            mods[n].strand = state->strand[i];
            mods[n].qual = state->ML[i] ? *state->ML[i] : -1;
        }
        n++;
        if (state->ML[i])
            state->ML[i] += (b->core.flag & BAM_FREVERSE)
                ? -state->MLstride[i]
                : +state->MLstride[i];

        if (b->core.flag & BAM_FREVERSE) {
            // process MM list backwards
            char *cp;
            for (cp = state->MMend[i]-1; cp != state->MM[i]; cp--)
                if (*cp == ',')
                    break;
            state->MMend[i] = cp;
            if (cp != state->MM[i])
                state->MMcount[i] = strtol(cp+1, NULL, 10);
            else
                state->MMcount[i] = INT_MAX;
        } else {
            if (*state->MM[i] == ',')
                state->MMcount[i] = strtol(state->MM[i]+1, &state->MM[i], 10);
            else
                state->MMcount[i] = INT_MAX;
        }

        // Multiple mods at the same coords.
        for (j=i+1; j < state->nmods && state->MM[j] == MMptr; j++) {
            if (n < n_mods) {
                mods[n].modified_base = state->type[j];
                mods[n].canonical_base = seq_nt16_str[state->canonical[j]];
                mods[n].strand = state->strand[j];
                mods[n].qual = state->ML[j] ? *state->ML[j] : -1;
            }
            n++;
            state->MMcount[j] = state->MMcount[i];
            state->MM[j]      = state->MM[i];
            if (state->ML[j])
                state->ML[j] += (b->core.flag & BAM_FREVERSE)
                    ? -state->MLstride[j]
                    : +state->MLstride[j];
        }
        i = j-1;
    }

    return n;
}

/*
 * Looks for the next location with a base modification.
 */
int bam_next_basemod(const bam1_t *b, hts_base_mod_state *state,
                     hts_base_mod *mods, int n_mods, int *pos) {
    if (state->seq_pos >= b->core.l_qseq)
        return 0;

    // Look through state->MMcount arrays to see when the next lowest is
    // per base type;
    int next[16], freq[16] = {0}, i;
    memset(next, 0x7f, 16*sizeof(*next));
    if (b->core.flag & BAM_FREVERSE) {
        for (i = 0; i < state->nmods; i++) {
            if (next[seqi_rc[state->canonical[i]]] > state->MMcount[i])
                next[seqi_rc[state->canonical[i]]] = state->MMcount[i];
        }
    } else {
        for (i = 0; i < state->nmods; i++) {
            if (next[state->canonical[i]] > state->MMcount[i])
                next[state->canonical[i]] = state->MMcount[i];
        }
    }

    // Now step through the sequence counting off base types.
    for (i = state->seq_pos; i < b->core.l_qseq; i++) {
        unsigned char bc = bam_seqi(bam_get_seq(b), i);
        if (next[bc] <= freq[bc] || next[15] <= freq[15])
            break;
        freq[bc]++;
        if (bc != 15) // N
            freq[15]++;
    }
    *pos = state->seq_pos = i;

    if (i >= b->core.l_qseq) {
        // Check for more MM elements than bases present.
        for (i = 0; i < state->nmods; i++) {
            if (!(b->core.flag & BAM_FREVERSE) &&
                state->MMcount[i] < 0x7f000000) {
                hts_log_warning("MM tag refers to bases beyond sequence length");
                return -1;
            }
        }
        return 0;
    }

    if (b->core.flag & BAM_FREVERSE) {
        for (i = 0; i < state->nmods; i++)
            state->MMcount[i] -= freq[seqi_rc[state->canonical[i]]];
    } else {
        for (i = 0; i < state->nmods; i++)
            state->MMcount[i] -= freq[state->canonical[i]];
    }

    int r = bam_mods_at_next_pos(b, state, mods, n_mods);
    return r > 0 ? r : 0;
}

/*
 * As per bam_mods_at_next_pos, but at a specific qpos >= the previous qpos.
 * This can only march forwards along the read, but can do so by more than
 * one base-pair.
 *
 * This makes it useful for calling from pileup iterators where qpos may
 * start part way through a read for the first occurrence of that record.
 */
int bam_mods_at_qpos(const bam1_t *b, int qpos, hts_base_mod_state *state,
                    hts_base_mod *mods, int n_mods) {
    // FIXME: for now this is inefficient in implementation.
    int r = 0;
    while (state->seq_pos <= qpos)
        if ((r = bam_mods_at_next_pos(b, state, mods, n_mods)) < 0)
            break;

    return r;
}

/*
 * Returns the list of base modification codes provided for this
 * alignment record as an array of character codes (+ve) or ChEBI numbers
 * (negative).
 *
 * Returns the array, with *ntype filled out with the size.
 *         The array returned should not be freed.
 *         It is a valid pointer until the state is freed using
 *         hts_base_mod_free().
 */
int *bam_mods_recorded(hts_base_mod_state *state, int *ntype) {
    *ntype = state->nmods;
    return state->type;
}

/*
 * Returns data about a specific modification type for the alignment record.
 * Code is either positive (eg 'm') or negative for ChEBI numbers.
 *
 * Return 0 on success or -1 if not found.  The strand, implicit and canonical
 * fields are filled out if passed in as non-NULL pointers.
 */
int bam_mods_query_type(hts_base_mod_state *state, int code,
                        int *strand, int *implicit, char *canonical) {
    // Find code entry
    int i;
    for (i = 0; i < state->nmods; i++) {
        if (state->type[i] == code)
            break;
    }
    if (i == state->nmods)
        return -1;

    // Return data
    if (strand)    *strand    = state->strand[i];
    if (implicit)  *implicit  = state->implicit[i];
    if (canonical) *canonical = "?AC?G???T??????N"[state->canonical[i]];

    return 0;
}
