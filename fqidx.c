/*  fqidx.c -- FASTQ random access.

    Copyright (C) 2008, 2009, 2013-2017, 2018 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

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

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>

#include "htslib/bgzf.h"
#include "htslib/fqidx.h"
#include "htslib/hfile.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "hts_internal.h"

typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} fqidx1_t;
KHASH_MAP_INIT_STR(s, fqidx1_t)

struct __fqidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int fqi_insert_index(fqidx_t *idx, const char* name, int64_t len, int line_len, int line_blen, uint64_t seq_offset, uint64_t qual_offset)
{
    if (!name) {
        hts_log_error("Malformed line");
        return -1;
    }

    char *name_key = strdup(name);
    int absent;
    khint_t k = kh_put(s, idx->hash, name_key, &absent);
    fqidx1_t *v = &kh_value(idx->hash, k);

    if (! absent) {
        hts_log_warning("Ignoring duplicate sequence \"%s\" at byte offset %"PRIu64"", name, seq_offset);
        free(name_key);
        return 0;
    }

    if (idx->n == idx->m) {
        char **tmp;
        idx->m = idx->m? idx->m<<1 : 16;
        if (!(tmp = (char**)realloc(idx->name, sizeof(char*) * idx->m))) {
            hts_log_error("Out of memory");
            return -1;
        }
        idx->name = tmp;
    }
    idx->name[idx->n++] = name_key;
    v->len = len;
    v->line_len = line_len;
    v->line_blen = line_blen;
    v->seq_offset = seq_offset;
    v->qual_offset = qual_offset;

    return 0;
}

fqidx_t *fqi_build_core(BGZF *bgzf)
{
    kstring_t name = { 0, 0, NULL };
    int c;
    int line_len, line_blen;
     enum state_e {GOTSEQ, GOTNAME, GOTLASTSEQ, BEYONDLASTSEQ, GOTQUALHEAD, GOTQUAL, GOTLASTQUAL, BEYONDLASTQUAL} state;
    int l1, l2;
    fqidx_t *idx;
    uint64_t offset, qual_offset;
    int64_t len,q_len;

    idx = (fqidx_t*)calloc(1, sizeof(fqidx_t));
    idx->hash = kh_init(s);
    len = q_len = line_len = line_blen = -1; state = GOTSEQ; l1 = l2 = -1; offset = 0, qual_offset = 0;
    while ( (c=bgzf_getc(bgzf))>=0 ) {
        if (c == '\n') { // an empty line
            if (state == GOTNAME) {
                offset = bgzf_utell(bgzf);
                continue;
            } else if ((state == GOTSEQ && len < 0) || state == GOTLASTSEQ) continue;
            else if (state == GOTSEQ) { state = GOTLASTSEQ; continue; }
        }
        else if (c == '@') { // fastq header
            if (len >= 0) {
                if (fqi_insert_index(idx, name.s, len, line_len, line_blen, offset, qual_offset) != 0)
                    goto fail;
            }

            // Read name in until we hit a newline or space
            name.l = 0;
            while ((c = bgzf_getc(bgzf)) >= 0)
                if (! isspace(c)) kputc_(c, &name);
                else if (name.l > 0 || c == '\n') break;
            kputsn("", 0, &name); // null terminate name

            if ( c<0 ) {
                hts_log_error("The last entry has no sequence");
                goto fail;
            }
            // If there is data after the name separated by a space then eat it
            if (c != '\n') while ( (c=bgzf_getc(bgzf))>=0 && c != '\n');
            state = GOTNAME, len = 0, q_len = 0;
            offset = bgzf_utell(bgzf);
            continue;
        }
        else if (c == '+') {
            state = GOTQUALHEAD;
            if (c != '\n') while ( (c=bgzf_getc(bgzf))>=0 && c != '\n');
            qual_offset = bgzf_utell(bgzf);
        } else {
            switch (state) {
                case BEYONDLASTSEQ: {
                hts_log_error("Inlined empty line is not allowed in sequence '%s'", name.s);
                goto fail;
                }
                case GOTLASTSEQ: state = BEYONDLASTSEQ;
                case GOTSEQ:
                case GOTNAME:
                    {
                        l1 = l2 = 0;
                        do {
                            ++l1; // count all length
                            if (isgraph(c)) ++l2; // count only seq length
                        } while ( (c=bgzf_getc(bgzf))>=0 && c != '\n'); // Read seq until
                        if (state == BEYONDLASTSEQ && l2) {
                            hts_log_error("Different line length in sequence '%s'", name.s);
                            goto fail;
                        }
                        ++l1; len += l2;
                        if (state == GOTNAME) { line_len = l1, line_blen = l2, state = GOTSEQ; }
                        else if (state == GOTSEQ) {
                            if (l1 != line_len || l2 != line_blen) state = GOTLASTSEQ;
                        }
                    }
                    break;
                case BEYONDLASTQUAL:{
                    hts_log_error("Inlined empty line is not allowed in quality of sequence '%s'", name.s);
                    goto fail;
                }
                case GOTLASTQUAL: state = BEYONDLASTQUAL;
                case GOTQUALHEAD:
                case GOTQUAL:
                    {
                        l1 = l2 = 0;
                        do {
                            ++l1; // count all length
                            if (isgraph(c)) ++l2; // count only seq length
                        } while ( (c=bgzf_getc(bgzf))>=0 && c != '\n'); // Read seq until
                        if (state == BEYONDLASTQUAL && l2) {
                            hts_log_error("Different line length in sequence '%s'", name.s);
                            goto fail;
                        }
                        ++l1; q_len += l2;
                        if (state == GOTQUALHEAD) { state = GOTQUAL; }
                        if (state == GOTQUAL) {
                            if (l1 != line_len || l2 != line_blen) state = GOTLASTQUAL;
                        }
                    }
                    break;
                default:
                    hts_log_error("Invalid state last known sequence '%s'", name.s);
                    goto fail;
            }
        }
    }

    if (len >= 0) {
        // Insert last entry into index
        if (fqi_insert_index(idx, name.s, len, line_len, line_blen, offset, qual_offset) != 0)
            goto fail;
    } else {
        goto fail;
    }

    free(name.s);
    return idx;

fail:
    free(name.s);
    fqi_destroy(idx);
    return NULL;
}

static int fqi_save(const fqidx_t *fqi, hFILE *fp) {
    khint_t k;
    int i;
    char buf[96]; // Must be big enough for format below.

    for (i = 0; i < fqi->n; ++i) {
        fqidx1_t x;
        k = kh_get(s, fqi->hash, fqi->name[i]);
        assert(k < kh_end(fqi->hash));
        x = kh_value(fqi->hash, k);
        snprintf(buf, sizeof(buf),
                 "\t%"PRId64"\t%"PRIu64"\t%"PRIu64"\t%"PRId32"\t%"PRId32"\n",
                 x.len, x.seq_offset, x.qual_offset, x.line_blen, x.line_len);
        if (hputs(fqi->name[i], fp) != 0) return -1;
        if (hputs(buf, fp) != 0) return -1;
    }
    return 0;
}

static fqidx_t *fqi_read(hFILE *fp, const char *fname)
{
    fqidx_t *fqi;
    char *buf = NULL, *p;
    int line_len, line_blen, n;
    int64_t len;
    uint64_t offset, qual_offset;
    ssize_t l, lnum = 1;

    fqi = (fqidx_t*)calloc(1, sizeof(fqidx_t));
    if (!fqi) return NULL;

    fqi->hash = kh_init(s);
    if (!fqi->hash) goto fail;

    buf = (char*)calloc(0x10000, 1);
    if (!buf) goto fail;

    while ((l = hgetln(buf, 0x10000, fp)) > 0) {
        for (p = buf; *p && !isspace_c(*p); ++p);
        if (p - buf < l) {
            *p = 0; ++p;
        }
        n = sscanf(p, "%"SCNd64"%"SCNu64"%"SCNu64"%d%d", &len, &offset, &qual_offset, &line_blen, &line_len);
        if (n != 5) {
            hts_log_error("Could not understand FAI %s line %zd", fname, lnum);
            goto fail;
        }
        if (fqi_insert_index(fqi, buf, len, line_len, line_blen, offset, qual_offset) != 0) {
            goto fail;
        }
        if (buf[l - 1] == '\n') ++lnum;
    }

    if (l < 0) {
        hts_log_error("Error while reading %s: %s", fname, strerror(errno));
        goto fail;
    }
    free(buf);
    return fqi;

 fail:
    free(buf);
    fqi_destroy(fqi);
    return NULL;
}

void fqi_destroy(fqidx_t *fqi)
{
    int i;
    if (!fqi) return;
    for (i = 0; i < fqi->n; ++i) free(fqi->name[i]);
    free(fqi->name);
    kh_destroy(s, fqi->hash);
    if (fqi->bgzf) bgzf_close(fqi->bgzf);
    free(fqi);
}

int fqi_build3(const char *fn, const char *fnfqi, const char *fngzi)
{
    kstring_t fqi_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    BGZF *bgzf = NULL;
    hFILE *fp = NULL;
    fqidx_t *fqi = NULL;
    int save_errno, res;

    if (!fnfqi) {
        if (ksprintf(&fqi_kstr, "%s.fqi", fn) < 0) goto fail;
        fnfqi = fqi_kstr.s;
    }
    if (!fngzi) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto fail;
        fngzi = gzi_kstr.s;
    }

    bgzf = bgzf_open(fn, "r");
    if ( !bgzf ) {
        hts_log_error("Failed to open the FASTQ file %s", fn);
        goto fail;
    }
    if ( bgzf->is_compressed ) {
        if (bgzf_index_build_init(bgzf) != 0) {
            hts_log_error("Failed to allocate bgzf index");
            goto fail;
        }
    }
    fqi = fqi_build_core(bgzf);
    if ( !fqi ) {
        if (bgzf->is_compressed && bgzf->is_gzip) {
            hts_log_error("Cannot index files compressed with gzip, please use bgzip");
        }
        goto fail;
    }
    if ( bgzf->is_compressed ) {
        if (bgzf_index_dump(bgzf, fngzi, NULL) < 0) {
            hts_log_error("Failed to make bgzf index %s", fngzi);
            goto fail;
        }
    }
    res = bgzf_close(bgzf);
    bgzf = NULL;
    if (res < 0) {
        hts_log_error("Error on closing %s : %s", fn, strerror(errno));
        goto fail;
    }
    fp = hopen(fnfqi, "wb");
    if ( !fp ) {
        hts_log_error("Failed to open FASTQ index %s : %s", fnfqi, strerror(errno));
        goto fail;
    }
    if (fqi_save(fqi, fp) != 0) {
        hts_log_error("Failed to write FASTQ index %s : %s", fnfqi, strerror(errno));
        goto fail;
    }
    if (hclose(fp) != 0) {
        hts_log_error("Failed on closing FASTQ index %s : %s", fnfqi, strerror(errno));
        goto fail;
    }

    free(fqi_kstr.s);
    free(gzi_kstr.s);
    fqi_destroy(fqi);
    return 0;

 fail:
    save_errno = errno;
    free(fqi_kstr.s);
    free(gzi_kstr.s);
    bgzf_close(bgzf);
    fqi_destroy(fqi);
    errno = save_errno;
    return -1;
}

int fqi_build(const char *fn) {
    return fqi_build3(fn, NULL, NULL);
}

fqidx_t *fqi_load3(const char* fn, const char *fnfqi, const char *fngzi,
                   int flags)
{
    kstring_t fqi_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    hFILE *fp = NULL;
    fqidx_t *fqi = NULL;
    int res;

    if (fn == NULL)
        return NULL;

    if (fnfqi == NULL) {
        if (ksprintf(&fqi_kstr, "%s.fqi", fn) < 0) goto fail;
        fnfqi = fqi_kstr.s;
    }
    if (fngzi == NULL) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto fail;
        fngzi = gzi_kstr.s;
    }

    fp = hopen(fnfqi, "rb");

    if (fp == 0) {
        if (!(flags & FQI_CREATE) || errno != ENOENT) {
            hts_log_error("Failed to open FASTQ index %s: %s", fnfqi, strerror(errno));
            goto fail;
        }

        hts_log_info("Build FASTQ index");

        if (fqi_build3(fn, fnfqi, fngzi) < 0) {
            goto fail;
        }

        fp = hopen(fnfqi, "rb");
        if (fp == 0) {
            hts_log_error("Failed to open FASTQ index %s: %s", fnfqi, strerror(errno));
            goto fail;
        }
    }

    fqi = fqi_read(fp, fnfqi);
    if (fqi == NULL) {
        hts_log_error("Failed to read FASTQ index %s", fnfqi);
        goto fail;
    }

    res = hclose(fp);
    fp = NULL;
    if (res < 0) {
        hts_log_error("Failed on closing FASTQ index %s : %s", fnfqi, strerror(errno));
        goto fail;
    }

    fqi->bgzf = bgzf_open(fn, "rb");
    if (fqi->bgzf == 0) {
        hts_log_error("Failed to open FASTQ file %s", fn);
        goto fail;
    }
    if ( fqi->bgzf->is_compressed==1 ) {
        if ( bgzf_index_load(fqi->bgzf, fngzi, NULL) < 0 ) {
            hts_log_error("Failed to load .gzi index: %s", fngzi);
            goto fail;
        }
    }
    free(fqi_kstr.s);
    free(gzi_kstr.s);
    return fqi;

 fail:
    if (fqi) fqi_destroy(fqi);
    if (fp) hclose_abruptly(fp);
    free(fqi_kstr.s);
    free(gzi_kstr.s);
    return NULL;
}

fqidx_t *fqi_load(const char* fn)
{
    return fqi_load3(fn, NULL, NULL, FQI_CREATE);
}

static char *fqi_retrieve(const fqidx_t *fqi, const fqidx1_t *val,
                          const long beg, const long end, int *len) {
    char *s;
    size_t l;
    int c = 0;
    int ret = bgzf_useek(fqi->bgzf,
                         val->seq_offset
                         + beg / val->line_blen * val->line_len
                         + beg % val->line_blen, SEEK_SET);

    if (ret < 0) {
        *len = -1;
        hts_log_error("Failed to retrieve block. (Seeking in a compressed, .gzi unindexed, file?)");
        return NULL;
    }

    l = 0;
    s = (char*)malloc((size_t) end - beg + 2);
    if (!s) {
        *len = -1;
        return NULL;
    }

    while ( l < end - beg && (c=bgzf_getc(fqi->bgzf))>=0 )
        if (isgraph(c)) s[l++] = c;
    if (c < 0) {
        hts_log_error("Failed to retrieve block: %s",
            c == -1 ? "unexpected end of file" : "error reading file");
        free(s);
        *len = -1;
        return NULL;
    }

    s[l] = '\0';
    *len = l < INT_MAX ? l : INT_MAX;
    return s;
}

static char *fqi_retrievequal(const fqidx_t *fqi, const fqidx1_t *val,
                          const long beg, const long end, int *len) {
    char *s;
    size_t l;
    int c = 0;
    int ret = bgzf_useek(fqi->bgzf,
                         val->qual_offset
                         + beg / val->line_blen * val->line_len
                         + beg % val->line_blen, SEEK_SET);
    
    if (ret < 0) {
        *len = -1;
        hts_log_error("Failed to retrieve block. (Seeking in a compressed, .gzi unindexed, file?)");
        return NULL;
    }
    
    l = 0;
    s = (char*)malloc((size_t) end - beg + 2);
    if (!s) {
        *len = -1;
        return NULL;
    }
    
    while ( l < end - beg && (c=bgzf_getc(fqi->bgzf))>=0 )
        if (isgraph(c)) s[l++] = c;
    if (c < 0) {
        hts_log_error("Failed to retrieve block: %s",
                      c == -1 ? "unexpected end of file" : "error reading file");
        free(s);
        *len = -1;
        return NULL;
    }
    
    s[l] = '\0';
    *len = l < INT_MAX ? l : INT_MAX;
    return s;
}

char *fqi_fetch(const fqidx_t *fqi, const char *str, int *len)
{
    char *s, *ep;
    size_t i, l, k, name_end;
    khiter_t iter;
    fqidx1_t val;
    khash_t(s) *h;
    long beg, end;

    beg = end = -1;
    h = fqi->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    if (!s) {
        *len = -1;
        return NULL;
    }

    // remove space
    for (i = k = 0; i < l; ++i)
        if (!isspace_c(str[i])) s[k++] = str[i];
    s[k] = 0;
    name_end = l = k;
    // determine the sequence name
    for (i = l; i > 0; --i) if (s[i - 1] == ':') break; // look for colon from the end
    if (i > 0) name_end = i - 1;
    if (name_end < l) { // check if this is really the end
        int n_hyphen = 0;
        for (i = name_end + 1; i < l; ++i) {
            if (s[i] == '-') ++n_hyphen;
            else if (!isdigit_c(s[i]) && s[i] != ',') break;
        }
        if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
        s[name_end] = 0;
        iter = kh_get(s, h, s);
        if (iter == kh_end(h)) { // cannot find the sequence name
            iter = kh_get(s, h, str); // try str as the name
            if (iter != kh_end(h)) {
                s[name_end] = ':';
                name_end = l;
            }
        }
    } else iter = kh_get(s, h, str);
    if(iter == kh_end(h)) {
        hts_log_warning("Reference %s not found in FASTQ file, returning empty sequence", str);
        free(s);
        *len = -2;
        return 0;
    }
    val = kh_value(h, iter);
    // parse the interval
    if (name_end < l) {
        int save_errno = errno;
        errno = 0;
        for (i = k = name_end + 1; i < l; ++i)
            if (s[i] != ',') s[k++] = s[i];
        s[k] = 0;
        if (s[name_end + 1] == '-') {
            beg = 0;
            i = name_end + 2;
        } else {
            beg = strtol(s + name_end + 1, &ep, 10);
            for (i = ep - s; i < k;) if (s[i++] == '-') break;
        }
        end = i < k? strtol(s + i, &ep, 10) : val.len;
        if (beg > 0) --beg;
        // Check for out of range numbers.  Only going to be a problem on
        // 32-bit platforms with >2Gb sequence length.
        if (errno == ERANGE && (uint64_t) val.len > LONG_MAX) {
            hts_log_error("Positions in range %s are too large for this platform", s);
            free(s);
            *len = -2;
            return NULL;
        }
        errno = save_errno;
    } else beg = 0, end = val.len;
    if (beg >= val.len) beg = val.len;
    if (end >= val.len) end = val.len;
    if (beg > end) beg = end;
    free(s);

    // now retrieve the sequence
    return fqi_retrieve(fqi, &val, beg, end, len);
}

char *fqi_fetchqual(const fqidx_t *fqi, const char *str, int *len)
{
    char *s, *ep;
    size_t i, l, k, name_end;
    khiter_t iter;
    fqidx1_t val;
    khash_t(s) *h;
    long beg, end;
    
    beg = end = -1;
    h = fqi->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    if (!s) {
        *len = -1;
        return NULL;
    }
    
    // remove space
    for (i = k = 0; i < l; ++i)
        if (!isspace_c(str[i])) s[k++] = str[i];
    s[k] = 0;
    name_end = l = k;
    // determine the sequence name
    for (i = l; i > 0; --i) if (s[i - 1] == ':') break; // look for colon from the end
    if (i > 0) name_end = i - 1;
    if (name_end < l) { // check if this is really the end
        int n_hyphen = 0;
        for (i = name_end + 1; i < l; ++i) {
            if (s[i] == '-') ++n_hyphen;
            else if (!isdigit_c(s[i]) && s[i] != ',') break;
        }
        if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
        s[name_end] = 0;
        iter = kh_get(s, h, s);
        if (iter == kh_end(h)) { // cannot find the sequence name
            iter = kh_get(s, h, str); // try str as the name
            if (iter != kh_end(h)) {
                s[name_end] = ':';
                name_end = l;
            }
        }
    } else iter = kh_get(s, h, str);
    if(iter == kh_end(h)) {
        hts_log_warning("Reference %s not found in FASTQ file, returning empty quality", str);
        free(s);
        *len = -2;
        return 0;
    }
    val = kh_value(h, iter);
    // parse the interval
    if (name_end < l) {
        int save_errno = errno;
        errno = 0;
        for (i = k = name_end + 1; i < l; ++i)
            if (s[i] != ',') s[k++] = s[i];
        s[k] = 0;
        if (s[name_end + 1] == '-') {
            beg = 0;
            i = name_end + 2;
        } else {
            beg = strtol(s + name_end + 1, &ep, 10);
            for (i = ep - s; i < k;) if (s[i++] == '-') break;
        }
        end = i < k? strtol(s + i, &ep, 10) : val.len;
        if (beg > 0) --beg;
        // Check for out of range numbers.  Only going to be a problem on
        // 32-bit platforms with >2Gb sequence length.
        if (errno == ERANGE && (uint64_t) val.len > LONG_MAX) {
            hts_log_error("Positions in range %s are too large for this platform", s);
            free(s);
            *len = -2;
            return NULL;
        }
        errno = save_errno;
    } else beg = 0, end = val.len;
    if (beg >= val.len) beg = val.len;
    if (end >= val.len) end = val.len;
    if (beg > end) beg = end;
    free(s);
    
    // now retrieve the sequence
    return fqi_retrievequal(fqi, &val, beg, end, len);
}

int fqidx_fetch_nseq(const fqidx_t *fqi)
{
    return fqi->n;
}

int fqidx_nseq(const fqidx_t *fqi)
{
    return fqi->n;
}

const char *fqidx_iseq(const fqidx_t *fqi, int i)
{
    return fqi->name[i];
}

int fqidx_seq_len(const fqidx_t *fqi, const char *seq)
{
    khint_t k = kh_get(s, fqi->hash, seq);
    if ( k == kh_end(fqi->hash) ) return -1;
    return kh_val(fqi->hash, k).len;
}

char *fqidx_fetch_seq(const fqidx_t *fqi, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    khiter_t iter;
    fqidx1_t val;

    // Adjust position
    iter = kh_get(s, fqi->hash, c_name);
    if (iter == kh_end(fqi->hash))
    {
        *len = -2;
        hts_log_error("The sequence \"%s\" not found", c_name);
        return NULL;
    }
    val = kh_value(fqi->hash, iter);
    if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence
    return fqi_retrieve(fqi, &val, p_beg_i, (long) p_end_i + 1, len);
}

int fqidx_has_seq(const fqidx_t *fqi, const char *seq)
{
    khiter_t iter = kh_get(s, fqi->hash, seq);
    if (iter == kh_end(fqi->hash)) return 0;
    return 1;
}

