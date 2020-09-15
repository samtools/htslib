/*  vcf.c -- VCF/BCF API functions.

    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2020 Genome Research Ltd.
    Portions copyright (C) 2014 Intel Corporation.

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

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>

#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "htslib/hfile.h"
#include "hts_internal.h"
#include "htslib/hts_endian.h"
#include "htslib/khash_str2int.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "htslib/kseq.h"
HTSLIB_EXPORT
uint32_t bcf_float_missing    = 0x7F800001;

HTSLIB_EXPORT
uint32_t bcf_float_vector_end = 0x7F800002;

HTSLIB_EXPORT
uint8_t bcf_type_shift[] = { 0, 0, 1, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static bcf_idinfo_t bcf_idinfo_def = { .info = { 15, 15, 15 }, .hrec = { NULL, NULL, NULL}, .id = -1 };

/*
    Partial support for 64-bit POS and Number=1 INFO tags.
    Notes:
     - the support for 64-bit values is motivated by POS and INFO/END for large genomes
     - the use of 64-bit values does not conform to the specification
     - cannot output 64-bit BCF and if it does, it is not compatible with anything
     - experimental, use at your risk
*/
#ifdef VCF_ALLOW_INT64
    #define BCF_MAX_BT_INT64 (0x7fffffffffffffff)       /* INT64_MAX, for internal use only */
    #define BCF_MIN_BT_INT64 -9223372036854775800LL     /* INT64_MIN + 8, for internal use only */
#endif

#define BCF_IS_64BIT (1<<30)


static char *find_chrom_header_line(char *s)
{
    char *nl;
    if (strncmp(s, "#CHROM\t", 7) == 0) return s;
    else if ((nl = strstr(s, "\n#CHROM\t")) != NULL) return nl+1;
    else return NULL;
}

/*************************
 *** VCF header parser ***
 *************************/

static int bcf_hdr_add_sample_len(bcf_hdr_t *h, const char *s, size_t len)
{
    if ( !s ) return 0;
    if (len == 0) len = strlen(s);

    const char *ss = s;
    while ( *ss && isspace_c(*ss) && ss - s < len) ss++;
    if ( !*ss || ss - s == len)
    {
        hts_log_error("Empty sample name: trailing spaces/tabs in the header line?");
        return -1;
    }

    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    int ret;
    char *sdup = malloc(len + 1);
    if (!sdup) return -1;
    memcpy(sdup, s, len);
    sdup[len] = 0;

    // Ensure space is available in h->samples
    size_t n = kh_size(d);
    char **new_samples = realloc(h->samples, sizeof(char*) * (n + 1));
    if (!new_samples) {
        free(sdup);
        return -1;
    }
    h->samples = new_samples;

    int k = kh_put(vdict, d, sdup, &ret);
    if (ret < 0) {
        free(sdup);
        return -1;
    }
    if (ret) { // absent
        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = n;
    } else {
        hts_log_error("Duplicated sample name '%s'", s);
        free(sdup);
        return -1;
    }
    h->samples[n] = sdup;
    h->dirty = 1;
    return 0;
}

int bcf_hdr_add_sample(bcf_hdr_t *h, const char *s)
{
    return bcf_hdr_add_sample_len(h, s, 0);
}

int HTS_RESULT_USED bcf_hdr_parse_sample_line(bcf_hdr_t *h, const char *str)
{
    int ret = 0;
    int i = 0;
    const char *p, *q;
    // add samples
    for (p = q = str;; ++q) {
        if (*q > '\n') continue;
        if (++i > 9) {
            if ( bcf_hdr_add_sample_len(h, p, q - p) < 0 ) ret = -1;
        }
        if (*q == 0 || *q == '\n' || ret < 0) break;
        p = q + 1;
    }

    return ret;
}

int bcf_hdr_sync(bcf_hdr_t *h)
{
    int i;
    for (i = 0; i < 3; i++)
    {
        vdict_t *d = (vdict_t*)h->dict[i];
        khint_t k;
        if ( h->n[i] < kh_size(d) )
        {
            bcf_idpair_t *new_idpair;
            // this should be true only for i=2, BCF_DT_SAMPLE
            new_idpair = (bcf_idpair_t*) realloc(h->id[i], kh_size(d)*sizeof(bcf_idpair_t));
            if (!new_idpair) return -1;
            h->n[i] = kh_size(d);
            h->id[i] = new_idpair;
        }
        for (k=kh_begin(d); k<kh_end(d); k++)
        {
            if (!kh_exist(d,k)) continue;
            h->id[i][kh_val(d,k).id].key = kh_key(d,k);
            h->id[i][kh_val(d,k).id].val = &kh_val(d,k);
        }
    }
    h->dirty = 0;
    return 0;
}

void bcf_hrec_destroy(bcf_hrec_t *hrec)
{
    if (!hrec) return;
    free(hrec->key);
    if ( hrec->value ) free(hrec->value);
    int i;
    for (i=0; i<hrec->nkeys; i++)
    {
        free(hrec->keys[i]);
        free(hrec->vals[i]);
    }
    free(hrec->keys);
    free(hrec->vals);
    free(hrec);
}

// Copies all fields except IDX.
bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec)
{
    int save_errno;
    bcf_hrec_t *out = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
    if (!out) return NULL;

    out->type = hrec->type;
    if ( hrec->key ) {
        out->key = strdup(hrec->key);
        if (!out->key) goto fail;
    }
    if ( hrec->value ) {
        out->value = strdup(hrec->value);
        if (!out->value) goto fail;
    }
    out->nkeys = hrec->nkeys;
    out->keys = (char**) malloc(sizeof(char*)*hrec->nkeys);
    if (!out->keys) goto fail;
    out->vals = (char**) malloc(sizeof(char*)*hrec->nkeys);
    if (!out->vals) goto fail;
    int i, j = 0;
    for (i=0; i<hrec->nkeys; i++)
    {
        if ( hrec->keys[i] && !strcmp("IDX",hrec->keys[i]) ) continue;
        if ( hrec->keys[i] ) {
            out->keys[j] = strdup(hrec->keys[i]);
            if (!out->keys[j]) goto fail;
        }
        if ( hrec->vals[i] ) {
            out->vals[j] = strdup(hrec->vals[i]);
            if (!out->vals[j]) goto fail;
        }
        j++;
    }
    if ( i!=j ) out->nkeys -= i-j;   // IDX was omitted
    return out;

 fail:
    save_errno = errno;
    hts_log_error("%s", strerror(errno));
    bcf_hrec_destroy(out);
    errno = save_errno;
    return NULL;
}

void bcf_hrec_debug(FILE *fp, bcf_hrec_t *hrec)
{
    fprintf(fp, "key=[%s] value=[%s]", hrec->key, hrec->value?hrec->value:"");
    int i;
    for (i=0; i<hrec->nkeys; i++)
        fprintf(fp, "\t[%s]=[%s]", hrec->keys[i],hrec->vals[i]);
    fprintf(fp, "\n");
}

void bcf_header_debug(bcf_hdr_t *hdr)
{
    int i, j;
    for (i=0; i<hdr->nhrec; i++)
    {
        if ( !hdr->hrec[i]->value )
        {
            fprintf(stderr, "##%s=<", hdr->hrec[i]->key);
            fprintf(stderr,"%s=%s", hdr->hrec[i]->keys[0], hdr->hrec[i]->vals[0]);
            for (j=1; j<hdr->hrec[i]->nkeys; j++)
                fprintf(stderr,",%s=%s", hdr->hrec[i]->keys[j], hdr->hrec[i]->vals[j]);
            fprintf(stderr,">\n");
        }
        else
            fprintf(stderr,"##%s=%s\n", hdr->hrec[i]->key,hdr->hrec[i]->value);
    }
}

int bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, size_t len)
{
    char **tmp;
    size_t n = hrec->nkeys + 1;
    assert(len > 0 && len < SIZE_MAX);
    tmp = realloc(hrec->keys, sizeof(char*)*n);
    if (!tmp) return -1;
    hrec->keys = tmp;
    tmp = realloc(hrec->vals, sizeof(char*)*n);
    if (!tmp) return -1;
    hrec->vals = tmp;

    hrec->keys[hrec->nkeys] = (char*) malloc((len+1)*sizeof(char));
    if (!hrec->keys[hrec->nkeys]) return -1;
    memcpy(hrec->keys[hrec->nkeys],str,len);
    hrec->keys[hrec->nkeys][len] = 0;
    hrec->vals[hrec->nkeys] = NULL;
    hrec->nkeys = n;
    return 0;
}

int bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, size_t len, int is_quoted)
{
    if ( hrec->vals[i] ) {
        free(hrec->vals[i]);
        hrec->vals[i] = NULL;
    }
    if ( !str ) return 0;
    if ( is_quoted )
    {
        if (len >= SIZE_MAX - 3) {
            errno = ENOMEM;
            return -1;
        }
        hrec->vals[i] = (char*) malloc((len+3)*sizeof(char));
        if (!hrec->vals[i]) return -1;
        hrec->vals[i][0] = '"';
        memcpy(&hrec->vals[i][1],str,len);
        hrec->vals[i][len+1] = '"';
        hrec->vals[i][len+2] = 0;
    }
    else
    {
        if (len == SIZE_MAX) {
            errno = ENOMEM;
            return -1;
        }
        hrec->vals[i] = (char*) malloc((len+1)*sizeof(char));
        if (!hrec->vals[i]) return -1;
        memcpy(hrec->vals[i],str,len);
        hrec->vals[i][len] = 0;
    }
    return 0;
}

int hrec_add_idx(bcf_hrec_t *hrec, int idx)
{
    int n = hrec->nkeys + 1;
    char **tmp = (char**) realloc(hrec->keys, sizeof(char*)*n);
    if (!tmp) return -1;
    hrec->keys = tmp;

    tmp = (char**) realloc(hrec->vals, sizeof(char*)*n);
    if (!tmp) return -1;
    hrec->vals = tmp;

    hrec->keys[hrec->nkeys] = strdup("IDX");
    if (!hrec->keys[hrec->nkeys]) return -1;

    kstring_t str = {0,0,0};
    if (kputw(idx, &str) < 0) {
        free(hrec->keys[hrec->nkeys]);
        return -1;
    }
    hrec->vals[hrec->nkeys] = str.s;
    hrec->nkeys = n;
    return 0;
}

int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key)
{
    int i;
    for (i=0; i<hrec->nkeys; i++)
        if ( !strcasecmp(key,hrec->keys[i]) ) return i;
    return -1;
}

static inline int is_escaped(const char *min, const char *str)
{
    int n = 0;
    while ( --str>=min && *str=='\\' ) n++;
    return n%2;
}

bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len)
{
    bcf_hrec_t *hrec = NULL;
    const char *p = line;
    if (p[0] != '#' || p[1] != '#') { *len = 0; return NULL; }
    p += 2;

    const char *q = p;
    while ( *q && *q!='=' && *q != '\n' ) q++;
    ptrdiff_t n = q-p;
    if ( *q!='=' || !n ) // wrong format
        goto malformed_line;

    hrec = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
    if (!hrec) { *len = -1; return NULL; }
    hrec->key = (char*) malloc(sizeof(char)*(n+1));
    if (!hrec->key) goto fail;
    memcpy(hrec->key,p,n);
    hrec->key[n] = 0;

    p = ++q;
    if ( *p!='<' ) // generic field, e.g. ##samtoolsVersion=0.1.18-r579
    {
        while ( *q && *q!='\n' ) q++;
        hrec->value = (char*) malloc((q-p+1)*sizeof(char));
        if (!hrec->value) goto fail;
        memcpy(hrec->value, p, q-p);
        hrec->value[q-p] = 0;
        *len = q - line + (*q ? 1 : 0); // Skip \n but not \0
        return hrec;
    }

    // structured line, e.g.
    // ##INFO=<ID=PV1,Number=1,Type=Float,Description="P-value for baseQ bias">
    // ##PEDIGREE=<Name_0=G0-ID,Name_1=G1-ID,Name_3=GN-ID>
    int nopen = 1;
    while ( *q && *q!='\n' && nopen>0 )
    {
        p = ++q;
        while ( *q && *q==' ' ) { p++; q++; }
        // ^[A-Za-z_][0-9A-Za-z_.]*$
        if (p==q && *q && (isalpha_c(*q) || *q=='_'))
        {
            q++;
            while ( *q && (isalnum_c(*q) || *q=='_' || *q=='.') ) q++;
        }
        n = q-p;
        int m = 0;
        while ( *q && *q==' ' ) { q++; m++; }
        if ( *q!='=' || !n )
            goto malformed_line;

        if (bcf_hrec_add_key(hrec, p, q-p-m) < 0) goto fail;
        p = ++q;
        while ( *q && *q==' ' ) { p++; q++; }
        int quoted = *p=='"' ? 1 : 0;
        if ( quoted ) p++, q++;
        while ( *q && *q != '\n' )
        {
            if ( quoted ) { if ( *q=='"' && !is_escaped(p,q) ) break; }
            else
            {
                if ( *q=='<' ) nopen++;
                if ( *q=='>' ) nopen--;
                if ( !nopen ) break;
                if ( *q==',' && nopen==1 ) break;
            }
            q++;
        }
        const char *r = q;
        while ( r > p && r[-1] == ' ' ) r--;
        if (bcf_hrec_set_val(hrec, hrec->nkeys-1, p, r-p, quoted) < 0)
            goto fail;
        if ( quoted && *q=='"' ) q++;
        if ( *q=='>' ) { nopen--; q++; }
    }

    // Skip to end of line
    int nonspace = 0;
    p = q;
    while ( *q && *q!='\n' ) { nonspace |= !isspace_c(*q); q++; }
    if (nonspace) {
        char buffer[320];
        hts_log_warning("Dropped trailing junk from header line '%s'",
                        hts_strprint(buffer, sizeof(buffer),
                                     '"', line, q - line));
    }

    *len = q - line + (*q ? 1 : 0);
    return hrec;

 fail:
    *len = -1;
    bcf_hrec_destroy(hrec);
    return NULL;

 malformed_line:
    {
        char buffer[320];
        while ( *q && *q!='\n' ) q++;  // Ensure *len includes full line
        hts_log_error("Could not parse the header line: %s",
                      hts_strprint(buffer, sizeof(buffer),
                                   '"', line, q - line));
        *len = q - line + (*q ? 1 : 0);
        bcf_hrec_destroy(hrec);
        return NULL;
    }
}

static int bcf_hdr_set_idx(bcf_hdr_t *hdr, int dict_type, const char *tag, bcf_idinfo_t *idinfo)
{
    size_t new_n;

    // If available, preserve existing IDX
    if ( idinfo->id==-1 )
        idinfo->id = hdr->n[dict_type];
    else if ( idinfo->id < hdr->n[dict_type] && hdr->id[dict_type][idinfo->id].key )
    {
        hts_log_error("Conflicting IDX=%d lines in the header dictionary, the new tag is %s",
            idinfo->id, tag);
        errno = EINVAL;
        return -1;
    }

    new_n = idinfo->id >= hdr->n[dict_type] ? idinfo->id+1 : hdr->n[dict_type];
    if (hts_resize(bcf_idpair_t, new_n, &hdr->m[dict_type],
                   &hdr->id[dict_type], HTS_RESIZE_CLEAR)) {
        return -1;
    }
    hdr->n[dict_type] = new_n;

    // NB: the next kh_put call can invalidate the idinfo pointer, therefore
    // we leave it unassigned here. It must be set explicitly in bcf_hdr_sync.
    hdr->id[dict_type][idinfo->id].key = tag;

    return 0;
}

// returns: 1 when hdr needs to be synced, -1 on error, 0 otherwise
static int bcf_hdr_register_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)
{
    // contig
    int i, ret, replacing = 0;
    khint_t k;
    char *str = NULL;

    if ( !strcmp(hrec->key, "contig") )
    {
        hts_pos_t len = 0;
        hrec->type = BCF_HL_CTG;

        // Get the contig ID ($str) and length ($j)
        i = bcf_hrec_find_key(hrec,"length");
        if ( i<0 ) len = 0;
        else {
            char *end = hrec->vals[i];
            len = strtoll(hrec->vals[i], &end, 10);
            if (end == hrec->vals[i] || len < 0) return 0;
        }

        i = bcf_hrec_find_key(hrec,"ID");
        if ( i<0 ) return 0;
        str = strdup(hrec->vals[i]);
        if (!str) return -1;

        // Register in the dictionary
        vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
        khint_t k = kh_get(vdict, d, str);
        if ( k != kh_end(d) ) { // already present
            free(str); str=NULL;
            if (kh_val(d, k).hrec[0] != NULL) // and not removed
                return 0;
            replacing = 1;
        } else {
            k = kh_put(vdict, d, str, &ret);
            if (ret < 0) { free(str); return -1; }
        }

        int idx = bcf_hrec_find_key(hrec,"IDX");
        if ( idx!=-1 )
        {
            char *tmp = hrec->vals[idx];
            idx = strtol(hrec->vals[idx], &tmp, 10);
            if ( *tmp || idx < 0 || idx >= INT_MAX - 1)
            {
                if (!replacing) {
                    kh_del(vdict, d, k);
                    free(str);
                }
                hts_log_warning("Error parsing the IDX tag, skipping");
                return 0;
            }
        }

        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = idx;
        kh_val(d, k).info[0] = len;
        kh_val(d, k).hrec[0] = hrec;
        if (bcf_hdr_set_idx(hdr, BCF_DT_CTG, kh_key(d,k), &kh_val(d,k)) < 0) {
            if (!replacing) {
                kh_del(vdict, d, k);
                free(str);
            }
            return -1;
        }
        if ( idx==-1 ) {
            if (hrec_add_idx(hrec, kh_val(d,k).id) < 0) {
               return -1;
            }
        }

        return 1;
    }

    if ( !strcmp(hrec->key, "INFO") ) hrec->type = BCF_HL_INFO;
    else if ( !strcmp(hrec->key, "FILTER") ) hrec->type = BCF_HL_FLT;
    else if ( !strcmp(hrec->key, "FORMAT") ) hrec->type = BCF_HL_FMT;
    else if ( hrec->nkeys>0 ) { hrec->type = BCF_HL_STR; return 1; }
    else return 0;

    // INFO/FILTER/FORMAT
    char *id = NULL;
    uint32_t type = UINT32_MAX, var = UINT32_MAX;
    int num = -1, idx = -1;
    for (i=0; i<hrec->nkeys; i++)
    {
        if ( !strcmp(hrec->keys[i], "ID") ) id = hrec->vals[i];
        else if ( !strcmp(hrec->keys[i], "IDX") )
        {
            char *tmp = hrec->vals[i];
            idx = strtol(hrec->vals[i], &tmp, 10);
            if ( *tmp || idx < 0 || idx >= INT_MAX - 1)
            {
                hts_log_warning("Error parsing the IDX tag, skipping");
                return 0;
            }
        }
        else if ( !strcmp(hrec->keys[i], "Type") )
        {
            if ( !strcmp(hrec->vals[i], "Integer") ) type = BCF_HT_INT;
            else if ( !strcmp(hrec->vals[i], "Float") ) type = BCF_HT_REAL;
            else if ( !strcmp(hrec->vals[i], "String") ) type = BCF_HT_STR;
            else if ( !strcmp(hrec->vals[i], "Character") ) type = BCF_HT_STR;
            else if ( !strcmp(hrec->vals[i], "Flag") ) type = BCF_HT_FLAG;
            else
            {
                hts_log_warning("The type \"%s\" is not supported, assuming \"String\"", hrec->vals[i]);
                type = BCF_HT_STR;
            }
        }
        else if ( !strcmp(hrec->keys[i], "Number") )
        {
            if ( !strcmp(hrec->vals[i],"A") ) var = BCF_VL_A;
            else if ( !strcmp(hrec->vals[i],"R") ) var = BCF_VL_R;
            else if ( !strcmp(hrec->vals[i],"G") ) var = BCF_VL_G;
            else if ( !strcmp(hrec->vals[i],".") ) var = BCF_VL_VAR;
            else
            {
                sscanf(hrec->vals[i],"%d",&num);
                var = BCF_VL_FIXED;
            }
            if (var != BCF_VL_FIXED) num = 0xfffff;
        }
    }
    if (hrec->type == BCF_HL_INFO || hrec->type == BCF_HL_FMT) {
        if (type == -1) {
            hts_log_warning("%s %s field has no Type defined. Assuming String",
                *hrec->key == 'I' ? "An" : "A", hrec->key);
            type = BCF_HT_STR;
        }
        if (var == -1) {
            hts_log_warning("%s %s field has no Number defined. Assuming '.'",
                *hrec->key == 'I' ? "An" : "A", hrec->key);
            var = BCF_VL_VAR;
        }
    }
    uint32_t info = ((((uint32_t)num) & 0xfffff)<<12 |
                     (var & 0xf) << 8 |
                     (type & 0xf) << 4 |
                     (((uint32_t) hrec->type) & 0xf));

    if ( !id ) return 0;
    str = strdup(id);
    if (!str) return -1;

    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_ID];
    k = kh_get(vdict, d, str);
    if ( k != kh_end(d) )
    {
        // already present
        free(str);
        if ( kh_val(d, k).hrec[info&0xf] ) return 0;
        kh_val(d, k).info[info&0xf] = info;
        kh_val(d, k).hrec[info&0xf] = hrec;
        if ( idx==-1 ) {
            if (hrec_add_idx(hrec, kh_val(d, k).id) < 0) {
                return -1;
            }
        }
        return 1;
    }
    k = kh_put(vdict, d, str, &ret);
    if (ret < 0) {
        free(str);
        return -1;
    }
    kh_val(d, k) = bcf_idinfo_def;
    kh_val(d, k).info[info&0xf] = info;
    kh_val(d, k).hrec[info&0xf] = hrec;
    kh_val(d, k).id = idx;
    if (bcf_hdr_set_idx(hdr, BCF_DT_ID, kh_key(d,k), &kh_val(d,k)) < 0) {
        kh_del(vdict, d, k);
        free(str);
        return -1;
    }
    if ( idx==-1 ) {
        if (hrec_add_idx(hrec, kh_val(d,k).id) < 0) {
            return -1;
        }
    }

    return 1;
}

int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)
{
    int res;
    if ( !hrec ) return 0;

    hrec->type = BCF_HL_GEN;
    res = bcf_hdr_register_hrec(hdr,hrec);
    if (res < 0) return -1;
    if ( !res )
    {
        // If one of the hashed field, then it is already present
        if ( hrec->type != BCF_HL_GEN )
        {
            bcf_hrec_destroy(hrec);
            return 0;
        }

        // Is one of the generic fields and already present?
        int i;
        for (i=0; i<hdr->nhrec; i++)
        {
            if ( hdr->hrec[i]->type!=BCF_HL_GEN ) continue;
            if ( !strcmp(hdr->hrec[i]->key,hrec->key) && !strcmp(hrec->key,"fileformat") ) break;
            if ( !strcmp(hdr->hrec[i]->key,hrec->key) && !strcmp(hdr->hrec[i]->value,hrec->value) ) break;
        }
        if ( i<hdr->nhrec )
        {
            bcf_hrec_destroy(hrec);
            return 0;
        }
    }

    // New record, needs to be added
    int n = hdr->nhrec + 1;
    bcf_hrec_t **new_hrec = realloc(hdr->hrec, n*sizeof(bcf_hrec_t*));
    if (!new_hrec) return -1;
    hdr->hrec = new_hrec;
    hdr->hrec[hdr->nhrec] = hrec;
    hdr->dirty = 1;
    hdr->nhrec = n;

    return hrec->type==BCF_HL_GEN ? 0 : 1;
}

/*
 *  Note that while querying of FLT,INFO,FMT,CTG lines is fast (the keys are hashed),
 *  the STR,GEN lines are searched for linearly in a linked list of all header lines.
 *  This may become a problem for VCFs with huge headers, we might need to build a
 *  dictionary for these lines as well.
 */
bcf_hrec_t *bcf_hdr_get_hrec(const bcf_hdr_t *hdr, int type, const char *key, const char *value, const char *str_class)
{
    int i;
    if ( type==BCF_HL_GEN )
    {
        for (i=0; i<hdr->nhrec; i++)
        {
            if ( hdr->hrec[i]->type!=type ) continue;
            if ( strcmp(hdr->hrec[i]->key,key) ) continue;
            if ( !value || !strcmp(hdr->hrec[i]->value,value) ) return hdr->hrec[i];
        }
        return NULL;
    }
    else if ( type==BCF_HL_STR )
    {
        for (i=0; i<hdr->nhrec; i++)
        {
            if ( hdr->hrec[i]->type!=type ) continue;
            if ( strcmp(hdr->hrec[i]->key,str_class) ) continue;
            int j = bcf_hrec_find_key(hdr->hrec[i],key);
            if ( j>=0 && !strcmp(hdr->hrec[i]->vals[j],value) ) return hdr->hrec[i];
        }
        return NULL;
    }
    vdict_t *d = type==BCF_HL_CTG ? (vdict_t*)hdr->dict[BCF_DT_CTG] : (vdict_t*)hdr->dict[BCF_DT_ID];
    khint_t k = kh_get(vdict, d, value);
    if ( k == kh_end(d) ) return NULL;
    return kh_val(d, k).hrec[type==BCF_HL_CTG?0:type];
}

void bcf_hdr_check_sanity(bcf_hdr_t *hdr)
{
    static int PL_warned = 0, GL_warned = 0;

    if ( !PL_warned )
    {
        int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "PL");
        if ( bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) && bcf_hdr_id2length(hdr,BCF_HL_FMT,id)!=BCF_VL_G )
        {
            hts_log_warning("PL should be declared as Number=G");
            PL_warned = 1;
        }
    }
    if ( !GL_warned )
    {
        int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GL");
        if ( bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) && bcf_hdr_id2length(hdr,BCF_HL_FMT,id)!=BCF_VL_G )
        {
            hts_log_warning("GL should be declared as Number=G");
            GL_warned = 1;
        }
    }
}

int bcf_hdr_parse(bcf_hdr_t *hdr, char *htxt)
{
    int len, done = 0;
    char *p = htxt;

    // Check sanity: "fileformat" string must come as first
    bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr,p,&len);
    if ( !hrec || !hrec->key || strcasecmp(hrec->key,"fileformat") )
        hts_log_warning("The first line should be ##fileformat; is the VCF/BCF header broken?");
    if (bcf_hdr_add_hrec(hdr, hrec) < 0) {
        bcf_hrec_destroy(hrec);
        return -1;
    }

    // The filter PASS must appear first in the dictionary
    hrec = bcf_hdr_parse_line(hdr,"##FILTER=<ID=PASS,Description=\"All filters passed\">",&len);
    if (!hrec || bcf_hdr_add_hrec(hdr, hrec) < 0) {
        bcf_hrec_destroy(hrec);
        return -1;
    }

    // Parse the whole header
    do {
        while (NULL != (hrec = bcf_hdr_parse_line(hdr, p, &len))) {
            if (bcf_hdr_add_hrec(hdr, hrec) < 0) {
                bcf_hrec_destroy(hrec);
                return -1;
            }
            p += len;
        }
        assert(hrec == NULL);
        if (len < 0) {
            // len < 0 indicates out-of-memory, or similar error
            hts_log_error("Could not parse header line: %s", strerror(errno));
            return -1;
        } else if (len > 0) {
            // Bad header line.  bcf_hdr_parse_line() will have logged it.
            // Skip and try again on the next line (p + len will be the start
            // of the next one).
            p += len;
            continue;
        }

        // Next should be the sample line.  If not, it was a malformed
        // header, in which case print a warning and skip (many VCF
        // operations do not really care about a few malformed lines).
        // In the future we may want to add a strict mode that errors in
        // this case.
        if ( strncmp("#CHROM\tPOS",p,10) != 0 ) {
            char *eol = strchr(p, '\n');
            if (*p != '\0') {
                char buffer[320];
                hts_log_warning("Could not parse header line: %s",
                                hts_strprint(buffer, sizeof(buffer),
                                               '"', p,
                                               eol ? (eol - p) : SIZE_MAX));
            }
            if (eol) {
                p = eol + 1; // Try from the next line.
            } else {
                done = -1; // No more lines left, give up.
            }
        } else {
            done = 1; // Sample line found
        }
    } while (!done);

    if (done < 0) {
        // No sample line is fatal.
        hts_log_error("Could not parse the header, sample line not found");
        return -1;
    }

    if (bcf_hdr_parse_sample_line(hdr,p) < 0)
        return -1;
    if (bcf_hdr_sync(hdr) < 0)
        return -1;
    bcf_hdr_check_sanity(hdr);
    return 0;
}

int bcf_hdr_append(bcf_hdr_t *hdr, const char *line)
{
    int len;
    bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr, (char*) line, &len);
    if ( !hrec ) return -1;
    if (bcf_hdr_add_hrec(hdr, hrec) < 0)
        return -1;
    return 0;
}

void bcf_hdr_remove(bcf_hdr_t *hdr, int type, const char *key)
{
    int i = 0;
    bcf_hrec_t *hrec;
    if ( !key )
    {
        while ( i<hdr->nhrec )
        {
            if ( hdr->hrec[i]->type!=type ) { i++; continue; }
            hrec = hdr->hrec[i];

            if ( type==BCF_HL_FLT || type==BCF_HL_INFO || type==BCF_HL_FMT || type== BCF_HL_CTG )
            {
                int j = bcf_hrec_find_key(hdr->hrec[i], "ID");
                if ( j>=0 )
                {
                    vdict_t *d = type==BCF_HL_CTG ? (vdict_t*)hdr->dict[BCF_DT_CTG] : (vdict_t*)hdr->dict[BCF_DT_ID];
                    khint_t k = kh_get(vdict, d, hdr->hrec[i]->vals[j]);
                    kh_val(d, k).hrec[type==BCF_HL_CTG?0:type] = NULL;
                }
            }

            hdr->dirty = 1;
            hdr->nhrec--;
            if ( i < hdr->nhrec )
                memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
            bcf_hrec_destroy(hrec);
        }
        return;
    }
    while (1)
    {
        if ( type==BCF_HL_FLT || type==BCF_HL_INFO || type==BCF_HL_FMT || type== BCF_HL_CTG )
        {
            hrec = bcf_hdr_get_hrec(hdr, type, "ID", key, NULL);
            if ( !hrec ) return;

            for (i=0; i<hdr->nhrec; i++)
                if ( hdr->hrec[i]==hrec ) break;
            assert( i<hdr->nhrec );

            vdict_t *d = type==BCF_HL_CTG ? (vdict_t*)hdr->dict[BCF_DT_CTG] : (vdict_t*)hdr->dict[BCF_DT_ID];
            khint_t k = kh_get(vdict, d, key);
            kh_val(d, k).hrec[type==BCF_HL_CTG?0:type] = NULL;
        }
        else
        {
            for (i=0; i<hdr->nhrec; i++)
            {
                if ( hdr->hrec[i]->type!=type ) continue;
                if ( type==BCF_HL_GEN )
                {
                    if ( !strcmp(hdr->hrec[i]->key,key) ) break;
                }
                else
                {
                    // not all structured lines have ID, we could be more sophisticated as in bcf_hdr_get_hrec()
                    int j = bcf_hrec_find_key(hdr->hrec[i], "ID");
                    if ( j>=0 && !strcmp(hdr->hrec[i]->vals[j],key) ) break;
                }
            }
            if ( i==hdr->nhrec ) return;
            hrec = hdr->hrec[i];
        }

        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
        hdr->dirty = 1;
    }
}

int bcf_hdr_printf(bcf_hdr_t *hdr, const char *fmt, ...)
{
    char tmp[256], *line = tmp;
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(line, sizeof(tmp), fmt, ap);
    va_end(ap);

    if (n >= sizeof(tmp)) {
        n++; // For trailing NUL
        line = (char*)malloc(n);
        if (!line)
            return -1;

        va_start(ap, fmt);
        vsnprintf(line, n, fmt, ap);
        va_end(ap);
    }

    int ret = bcf_hdr_append(hdr, line);

    if (line != tmp) free(line);
    return ret;
}


/**********************
 *** BCF header I/O ***
 **********************/

const char *bcf_hdr_get_version(const bcf_hdr_t *hdr)
{
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_GEN, "fileformat", NULL, NULL);
    if ( !hrec )
    {
        hts_log_warning("No version string found, assuming VCFv4.2");
        return "VCFv4.2";
    }
    return hrec->value;
}

int bcf_hdr_set_version(bcf_hdr_t *hdr, const char *version)
{
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_GEN, "fileformat", NULL, NULL);
    if ( !hrec )
    {
        int len;
        kstring_t str = {0,0,0};
        ksprintf(&str,"##fileformat=%s", version);
        hrec = bcf_hdr_parse_line(hdr, str.s, &len);
        free(str.s);
    }
    else
    {
        free(hrec->value);
        hrec->value = strdup(version);
    }
    hdr->dirty = 1;
    return 0; // FIXME: check for errs in this function (return < 0 if so)
}

bcf_hdr_t *bcf_hdr_init(const char *mode)
{
    int i;
    bcf_hdr_t *h;
    h = (bcf_hdr_t*)calloc(1, sizeof(bcf_hdr_t));
    if (!h) return NULL;
    for (i = 0; i < 3; ++i)
        if ((h->dict[i] = kh_init(vdict)) == NULL) goto fail;
    if ( strchr(mode,'w') )
    {
        bcf_hdr_append(h, "##fileformat=VCFv4.2");
        // The filter PASS must appear first in the dictionary
        bcf_hdr_append(h, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    }
    return h;

 fail:
    for (i = 0; i < 3; ++i)
        kh_destroy(vdict, h->dict[i]);
    free(h);
    return NULL;
}

void bcf_hdr_destroy(bcf_hdr_t *h)
{
    int i;
    khint_t k;
    if (!h) return;
    for (i = 0; i < 3; ++i) {
        vdict_t *d = (vdict_t*)h->dict[i];
        if (d == 0) continue;
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
        kh_destroy(vdict, d);
        free(h->id[i]);
    }
    for (i=0; i<h->nhrec; i++)
        bcf_hrec_destroy(h->hrec[i]);
    if (h->nhrec) free(h->hrec);
    if (h->samples) free(h->samples);
    free(h->keep_samples);
    free(h->transl[0]); free(h->transl[1]);
    free(h->mem.s);
    free(h);
}

bcf_hdr_t *bcf_hdr_read(htsFile *hfp)
{
    if (hfp->format.format == vcf)
        return vcf_hdr_read(hfp);
    if (hfp->format.format != bcf) {
        hts_log_error("Input is not detected as bcf or vcf format");
        return NULL;
    }

    assert(hfp->is_bgzf);

    BGZF *fp = hfp->fp.bgzf;
    uint8_t magic[5];
    bcf_hdr_t *h;
    h = bcf_hdr_init("r");
    if (!h) {
        hts_log_error("Failed to allocate bcf header");
        return NULL;
    }
    if (bgzf_read(fp, magic, 5) != 5)
    {
        hts_log_error("Failed to read the header (reading BCF in text mode?)");
        bcf_hdr_destroy(h);
        return NULL;
    }
    if (strncmp((char*)magic, "BCF\2\2", 5) != 0)
    {
        if (!strncmp((char*)magic, "BCF", 3))
            hts_log_error("Invalid BCF2 magic string: only BCFv2.2 is supported");
        else
            hts_log_error("Invalid BCF2 magic string");
        bcf_hdr_destroy(h);
        return NULL;
    }
    uint8_t buf[4];
    size_t hlen;
    char *htxt = NULL;
    if (bgzf_read(fp, buf, 4) != 4) goto fail;
    hlen = buf[0] | (buf[1] << 8) | (buf[2] << 16) | ((size_t) buf[3] << 24);
    if (hlen >= SIZE_MAX) { errno = ENOMEM; goto fail; }
    htxt = (char*)malloc(hlen + 1);
    if (!htxt) goto fail;
    if (bgzf_read(fp, htxt, hlen) != hlen) goto fail;
    htxt[hlen] = '\0'; // Ensure htxt is terminated
    if ( bcf_hdr_parse(h, htxt) < 0 ) goto fail;
    free(htxt);
    return h;
 fail:
    hts_log_error("Failed to read BCF header");
    free(htxt);
    bcf_hdr_destroy(h);
    return NULL;
}

int bcf_hdr_write(htsFile *hfp, bcf_hdr_t *h)
{
    if (!h) {
        errno = EINVAL;
        return -1;
    }
    if ( h->dirty ) {
        if (bcf_hdr_sync(h) < 0) return -1;
    }
    hfp->format.category = variant_data;
    if (hfp->format.format == vcf || hfp->format.format == text_format) {
        hfp->format.format = vcf;
        return vcf_hdr_write(hfp, h);
    }

    if (hfp->format.format == binary_format)
        hfp->format.format = bcf;

    kstring_t htxt = {0,0,0};
    if (bcf_hdr_format(h, 1, &htxt) < 0) {
        free(htxt.s);
        return -1;
    }
    kputc('\0', &htxt); // include the \0 byte

    BGZF *fp = hfp->fp.bgzf;
    if ( bgzf_write(fp, "BCF\2\2", 5) !=5 ) return -1;
    uint8_t hlen[4];
    u32_to_le(htxt.l, hlen);
    if ( bgzf_write(fp, hlen, 4) !=4 ) return -1;
    if ( bgzf_write(fp, htxt.s, htxt.l) != htxt.l ) return -1;

    free(htxt.s);
    return 0;
}

/********************
 *** BCF site I/O ***
 ********************/

bcf1_t *bcf_init()
{
    bcf1_t *v;
    v = (bcf1_t*)calloc(1, sizeof(bcf1_t));
    return v;
}

void bcf_clear(bcf1_t *v)
{
    int i;
    for (i=0; i<v->d.m_info; i++)
    {
        if ( v->d.info[i].vptr_free )
        {
            free(v->d.info[i].vptr - v->d.info[i].vptr_off);
            v->d.info[i].vptr_free = 0;
        }
    }
    for (i=0; i<v->d.m_fmt; i++)
    {
        if ( v->d.fmt[i].p_free )
        {
            free(v->d.fmt[i].p - v->d.fmt[i].p_off);
            v->d.fmt[i].p_free = 0;
        }
    }
    v->rid = v->pos = v->rlen = v->unpacked = 0;
    bcf_float_set_missing(v->qual);
    v->n_info = v->n_allele = v->n_fmt = v->n_sample = 0;
    v->shared.l = v->indiv.l = 0;
    v->d.var_type = -1;
    v->d.shared_dirty = 0;
    v->d.indiv_dirty  = 0;
    v->d.n_flt = 0;
    v->errcode = 0;
    if (v->d.m_als) v->d.als[0] = 0;
    if (v->d.m_id) v->d.id[0] = 0;
}

void bcf_empty(bcf1_t *v)
{
    bcf_clear1(v);
    free(v->d.id);
    free(v->d.als);
    free(v->d.allele); free(v->d.flt); free(v->d.info); free(v->d.fmt);
    if (v->d.var ) free(v->d.var);
    free(v->shared.s); free(v->indiv.s);
    memset(&v->d,0,sizeof(v->d));
    memset(&v->shared,0,sizeof(v->shared));
    memset(&v->indiv,0,sizeof(v->indiv));
}

void bcf_destroy(bcf1_t *v)
{
    if (!v) return;
    bcf_empty1(v);
    free(v);
}

static inline int bcf_read1_core(BGZF *fp, bcf1_t *v)
{
    uint8_t x[32];
    ssize_t ret;
    uint32_t shared_len, indiv_len;
    if ((ret = bgzf_read(fp, x, 32)) != 32) {
        if (ret == 0) return -1;
        return -2;
    }
    bcf_clear1(v);
    shared_len = le_to_u32(x);
    if (shared_len < 24) return -2;
    shared_len -= 24; // to exclude six 32-bit integers
    if (ks_resize(&v->shared, shared_len) != 0) return -2;
    indiv_len = le_to_u32(x + 4);
    if (ks_resize(&v->indiv, indiv_len) != 0) return -2;
    v->rid  = le_to_i32(x + 8);
    v->pos  = le_to_u32(x + 12);
    v->rlen = le_to_i32(x + 16);
    v->qual = le_to_float(x + 20);
    v->n_info = le_to_u16(x + 24);
    v->n_allele = le_to_u16(x + 26);
    v->n_sample = le_to_u32(x + 28) & 0xffffff;
    v->n_fmt = x[31];
    v->shared.l = shared_len;
    v->indiv.l = indiv_len;
    // silent fix of broken BCFs produced by earlier versions of bcf_subset, prior to and including bd6ed8b4
    if ( (!v->indiv.l || !v->n_sample) && v->n_fmt ) v->n_fmt = 0;

    if (bgzf_read(fp, v->shared.s, v->shared.l) != v->shared.l) return -2;
    if (bgzf_read(fp, v->indiv.s, v->indiv.l) != v->indiv.l) return -2;
    return 0;
}

#define bit_array_size(n) ((n)/8+1)
#define bit_array_set(a,i)   ((a)[(i)/8] |=   1 << ((i)%8))
#define bit_array_clear(a,i) ((a)[(i)/8] &= ~(1 << ((i)%8)))
#define bit_array_test(a,i)  ((a)[(i)/8] &   (1 << ((i)%8)))

static int bcf_dec_typed_int1_safe(uint8_t *p, uint8_t *end, uint8_t **q,
                                   int32_t *val) {
    uint32_t t;
    if (end - p < 2) return -1;
    t = *p++ & 0xf;
    /* Use if .. else if ... else instead of switch to force order.  Assumption
       is that small integers are more frequent than big ones. */
    if (t == BCF_BT_INT8) {
        *val = *(int8_t *) p++;
    } else {
        if (end - p < (1<<bcf_type_shift[t])) return -1;
        if (t == BCF_BT_INT16) {
            *val = le_to_i16(p);
            p += 2;
        } else if (t == BCF_BT_INT32) {
            *val = le_to_i32(p);
            p += 4;
#ifdef VCF_ALLOW_INT64
        } else if (t == BCF_BT_INT64) {
            // This case should never happen because there should be no
            // 64-bit BCFs at all, definitely not coming from htslib
            *val = le_to_i64(p);
            p += 8;
#endif
        } else {
            return -1;
        }
    }
    *q = p;
    return 0;
}

static int bcf_dec_size_safe(uint8_t *p, uint8_t *end, uint8_t **q,
                             int *num, int *type) {
    int r;
    if (p >= end) return -1;
    *type = *p & 0xf;
    if (*p>>4 != 15) {
        *q = p + 1;
        *num = *p >> 4;
        return 0;
    }
    r = bcf_dec_typed_int1_safe(p + 1, end, q, num);
    if (r) return r;
    return *num >= 0 ? 0 : -1;
}

static const char *get_type_name(int type) {
    const char *types[9] = {
        "null", "int (8-bit)", "int (16 bit)", "int (32 bit)",
        "unknown", "float", "unknown", "char", "unknown"
    };
    int t = (type >= 0 && type < 8) ? type : 8;
    return types[t];
}

static void bcf_record_check_err(const bcf_hdr_t *hdr, bcf1_t *rec,
                                 char *type, uint32_t *reports, int i) {
    if (*reports == 0 || hts_verbose >= HTS_LOG_DEBUG)
        hts_log_warning("Bad BCF record at %s:%"PRIhts_pos
                        ": Invalid FORMAT %s %d",
                        bcf_seqname_safe(hdr,rec), rec->pos+1, type, i);
    (*reports)++;
}

static int bcf_record_check(const bcf_hdr_t *hdr, bcf1_t *rec) {
    uint8_t *ptr, *end;
    size_t bytes;
    uint32_t err = 0;
    int type = 0;
    int num  = 0;
    int reflen = 0;
    uint32_t i, reports;
    const uint32_t is_integer = ((1 << BCF_BT_INT8)  |
                                 (1 << BCF_BT_INT16) |
#ifdef VCF_ALLOW_INT64
                                 (1 << BCF_BT_INT64) |
#endif
                                 (1 << BCF_BT_INT32));
    const uint32_t is_valid_type = (is_integer          |
                                    (1 << BCF_BT_NULL)  |
                                    (1 << BCF_BT_FLOAT) |
                                    (1 << BCF_BT_CHAR));
    int32_t max_id = hdr ? hdr->n[BCF_DT_ID] : 0;

    // Check for valid contig ID
    if (rec->rid < 0
        || (hdr && (rec->rid >= hdr->n[BCF_DT_CTG]
                    || hdr->id[BCF_DT_CTG][rec->rid].key == NULL))) {
        hts_log_warning("Bad BCF record at %"PRIhts_pos": Invalid %s id %d", rec->pos+1, "CONTIG", rec->rid);
        err |= BCF_ERR_CTG_INVALID;
    }

    // Check ID
    ptr = (uint8_t *) rec->shared.s;
    end = ptr + rec->shared.l;
    if (bcf_dec_size_safe(ptr, end, &ptr, &num, &type) != 0) goto bad_shared;
    if (type != BCF_BT_CHAR) {
        hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s type %d (%s)", bcf_seqname_safe(hdr,rec), rec->pos+1, "ID", type, get_type_name(type));
        err |= BCF_ERR_TAG_INVALID;
    }
    bytes = (size_t) num << bcf_type_shift[type];
    if (end - ptr < bytes) goto bad_shared;
    ptr += bytes;

    // Check REF and ALT
    reports = 0;
    for (i = 0; i < rec->n_allele; i++) {
        if (bcf_dec_size_safe(ptr, end, &ptr, &num, &type) != 0) goto bad_shared;
        if (type != BCF_BT_CHAR) {
            if (!reports++ || hts_verbose >= HTS_LOG_DEBUG)
                hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s type %d (%s)", bcf_seqname_safe(hdr,rec), rec->pos+1, "REF/ALT", type, get_type_name(type));
            err |= BCF_ERR_CHAR;
        }
        if (i == 0) reflen = num;
        bytes = (size_t) num << bcf_type_shift[type];
        if (end - ptr < bytes) goto bad_shared;
        ptr += bytes;
    }

    // Check FILTER
    reports = 0;
    if (bcf_dec_size_safe(ptr, end, &ptr, &num, &type) != 0) goto bad_shared;
    if (num > 0) {
        bytes = (size_t) num << bcf_type_shift[type];
        if (((1 << type) & is_integer) == 0) {
            hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s type %d (%s)", bcf_seqname_safe(hdr,rec), rec->pos+1, "FILTER", type, get_type_name(type));
            err |= BCF_ERR_TAG_INVALID;
            if (end - ptr < bytes) goto bad_shared;
            ptr += bytes;
        } else {
            if (end - ptr < bytes) goto bad_shared;
            for (i = 0; i < num; i++) {
                int32_t key = bcf_dec_int1(ptr, type, &ptr);
                if (key < 0
                    || (hdr && (key >= max_id
                                || hdr->id[BCF_DT_ID][key].key == NULL))) {
                    if (!reports++ || hts_verbose >= HTS_LOG_DEBUG)
                        hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s id %d", bcf_seqname_safe(hdr,rec), rec->pos+1, "FILTER", key);
                    err |= BCF_ERR_TAG_UNDEF;
                }
            }
        }
    }

    // Check INFO
    reports = 0;
    bcf_idpair_t *id_tmp = hdr ? hdr->id[BCF_DT_ID] : NULL;
    for (i = 0; i < rec->n_info; i++) {
        int32_t key = -1;
        if (bcf_dec_typed_int1_safe(ptr, end, &ptr, &key) != 0) goto bad_shared;
        if (key < 0 || (hdr && (key >= max_id
                                || id_tmp[key].key == NULL))) {
            if (!reports++ || hts_verbose >= HTS_LOG_DEBUG)
                hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s id %d", bcf_seqname_safe(hdr,rec), rec->pos+1, "INFO", key);
            err |= BCF_ERR_TAG_UNDEF;
        }
        if (bcf_dec_size_safe(ptr, end, &ptr, &num, &type) != 0) goto bad_shared;
        if (((1 << type) & is_valid_type) == 0) {
            if (!reports++ || hts_verbose >= HTS_LOG_DEBUG)
                hts_log_warning("Bad BCF record at %s:%"PRIhts_pos": Invalid %s type %d (%s)", bcf_seqname_safe(hdr,rec), rec->pos+1, "INFO", type, get_type_name(type));
            err |= BCF_ERR_TAG_INVALID;
        }
        bytes = (size_t) num << bcf_type_shift[type];
        if (end - ptr < bytes) goto bad_shared;
        ptr += bytes;
    }

    // Check FORMAT and individual information
    ptr = (uint8_t *) rec->indiv.s;
    end = ptr + rec->indiv.l;
    reports = 0;
    for (i = 0; i < rec->n_fmt; i++) {
        int32_t key = -1;
        if (bcf_dec_typed_int1_safe(ptr, end, &ptr, &key) != 0) goto bad_indiv;
        if (key < 0
            || (hdr && (key >= max_id
                        || id_tmp[key].key == NULL))) {
            bcf_record_check_err(hdr, rec, "id", &reports, key);
            err |= BCF_ERR_TAG_UNDEF;
        }
        if (bcf_dec_size_safe(ptr, end, &ptr, &num, &type) != 0) goto bad_indiv;
        if (((1 << type) & is_valid_type) == 0) {
            bcf_record_check_err(hdr, rec, "type", &reports, type);
            err |= BCF_ERR_TAG_INVALID;
        }
        bytes = ((size_t) num << bcf_type_shift[type]) * rec->n_sample;
        if (end - ptr < bytes) goto bad_indiv;
        ptr += bytes;
    }

    if (!err && rec->rlen < 0) {
        // Treat bad rlen as a warning instead of an error, and try to
        // fix up by using the length of the stored REF allele.
        static int warned = 0;
        if (!warned) {
            hts_log_warning("BCF record at %s:%"PRIhts_pos" has invalid RLEN (%"PRIhts_pos"). "
                            "Only one invalid RLEN will be reported.",
                            bcf_seqname_safe(hdr,rec), rec->pos+1, rec->rlen);
            warned = 1;
        }
        rec->rlen = reflen >= 0 ? reflen : 0;
    }

    rec->errcode |= err;

    return err ? -2 : 0; // Return -2 so bcf_read() reports an error

 bad_shared:
    hts_log_error("Bad BCF record at %s:%"PRIhts_pos" - shared section malformed or too short", bcf_seqname_safe(hdr,rec), rec->pos+1);
    return -2;

 bad_indiv:
    hts_log_error("Bad BCF record at %s:%"PRIhts_pos" - individuals section malformed or too short", bcf_seqname_safe(hdr,rec), rec->pos+1);
    return -2;
}

static inline uint8_t *bcf_unpack_fmt_core1(uint8_t *ptr, int n_sample, bcf_fmt_t *fmt);
int bcf_subset_format(const bcf_hdr_t *hdr, bcf1_t *rec)
{
    if ( !hdr->keep_samples ) return 0;
    if ( !bcf_hdr_nsamples(hdr) )
    {
        rec->indiv.l = rec->n_sample = 0;
        return 0;
    }

    int i, j;
    uint8_t *ptr = (uint8_t*)rec->indiv.s, *dst = NULL, *src;
    bcf_dec_t *dec = &rec->d;
    hts_expand(bcf_fmt_t, rec->n_fmt, dec->m_fmt, dec->fmt);
    for (i=0; i<dec->m_fmt; ++i) dec->fmt[i].p_free = 0;

    for (i=0; i<rec->n_fmt; i++)
    {
        ptr = bcf_unpack_fmt_core1(ptr, rec->n_sample, &dec->fmt[i]);
        src = dec->fmt[i].p - dec->fmt[i].size;
        if ( dst )
        {
            memmove(dec->fmt[i-1].p + dec->fmt[i-1].p_len, dec->fmt[i].p - dec->fmt[i].p_off, dec->fmt[i].p_off);
            dec->fmt[i].p = dec->fmt[i-1].p + dec->fmt[i-1].p_len + dec->fmt[i].p_off;
        }
        dst = dec->fmt[i].p;
        for (j=0; j<hdr->nsamples_ori; j++)
        {
            src += dec->fmt[i].size;
            if ( !bit_array_test(hdr->keep_samples,j) ) continue;
            memmove(dst, src, dec->fmt[i].size);
            dst += dec->fmt[i].size;
        }
        rec->indiv.l -= dec->fmt[i].p_len - (dst - dec->fmt[i].p);
        dec->fmt[i].p_len = dst - dec->fmt[i].p;
    }
    rec->unpacked |= BCF_UN_FMT;

    rec->n_sample = bcf_hdr_nsamples(hdr);
    return 0;
}

int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
    if (fp->format.format == vcf) return vcf_read(fp,h,v);
    int ret = bcf_read1_core(fp->fp.bgzf, v);
    if (ret == 0) ret = bcf_record_check(h, v);
    if ( ret!=0 || !h->keep_samples ) return ret;
    return bcf_subset_format(h,v);
}

int bcf_readrec(BGZF *fp, void *null, void *vv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    bcf1_t *v = (bcf1_t *) vv;
    int ret = bcf_read1_core(fp, v);
    if (ret == 0) ret = bcf_record_check(NULL, v);
    if (ret  >= 0)
        *tid = v->rid, *beg = v->pos, *end = v->pos + v->rlen;
    return ret;
}

static inline int bcf1_sync_id(bcf1_t *line, kstring_t *str)
{
    // single typed string
    if ( line->d.id && strcmp(line->d.id, ".") ) {
        return bcf_enc_vchar(str, strlen(line->d.id), line->d.id);
    } else {
        return bcf_enc_size(str, 0, BCF_BT_CHAR);
    }
}
static inline int bcf1_sync_alleles(bcf1_t *line, kstring_t *str)
{
    // list of typed strings
    int i;
    for (i=0; i<line->n_allele; i++) {
        if (bcf_enc_vchar(str, strlen(line->d.allele[i]), line->d.allele[i]) < 0)
            return -1;
    }
    if ( !line->rlen && line->n_allele ) line->rlen = strlen(line->d.allele[0]);
    return 0;
}
static inline int bcf1_sync_filter(bcf1_t *line, kstring_t *str)
{
    // typed vector of integers
    if ( line->d.n_flt ) {
        return bcf_enc_vint(str, line->d.n_flt, line->d.flt, -1);
    } else {
        return bcf_enc_vint(str, 0, 0, -1);
    }
}

static inline int bcf1_sync_info(bcf1_t *line, kstring_t *str)
{
    // pairs of typed vectors
    int i, irm = -1, e = 0;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *info = &line->d.info[i];
        if ( !info->vptr )
        {
            // marked for removal
            if ( irm < 0 ) irm = i;
            continue;
        }
        e |= kputsn_(info->vptr - info->vptr_off, info->vptr_len + info->vptr_off, str) < 0;
        if ( irm >=0 )
        {
            bcf_info_t tmp = line->d.info[irm]; line->d.info[irm] = line->d.info[i]; line->d.info[i] = tmp;
            while ( irm<=i && line->d.info[irm].vptr ) irm++;
        }
    }
    if ( irm>=0 ) line->n_info = irm;
    return e == 0 ? 0 : -1;
}

static int bcf1_sync(bcf1_t *line)
{
    char *shared_ori = line->shared.s;
    size_t prev_len;

    kstring_t tmp = {0,0,0};
    if ( !line->shared.l )
    {
        // New line created via API, BCF data blocks do not exist. Get it ready for BCF output
        tmp = line->shared;
        bcf1_sync_id(line, &tmp);
        line->unpack_size[0] = tmp.l; prev_len = tmp.l;

        bcf1_sync_alleles(line, &tmp);
        line->unpack_size[1] = tmp.l - prev_len; prev_len = tmp.l;

        bcf1_sync_filter(line, &tmp);
        line->unpack_size[2] = tmp.l - prev_len;

        bcf1_sync_info(line, &tmp);
        line->shared = tmp;
    }
    else if ( line->d.shared_dirty )
    {
        // The line was edited, update the BCF data block.

        if ( !(line->unpacked & BCF_UN_STR) ) bcf_unpack(line,BCF_UN_STR);

        // ptr_ori points to the original unchanged BCF data.
        uint8_t *ptr_ori = (uint8_t *) line->shared.s;

        // ID: single typed string
        if ( line->d.shared_dirty & BCF1_DIRTY_ID )
            bcf1_sync_id(line, &tmp);
        else
            kputsn_(ptr_ori, line->unpack_size[0], &tmp);
        ptr_ori += line->unpack_size[0];
        line->unpack_size[0] = tmp.l; prev_len = tmp.l;

        // REF+ALT: list of typed strings
        if ( line->d.shared_dirty & BCF1_DIRTY_ALS )
            bcf1_sync_alleles(line, &tmp);
        else
        {
            kputsn_(ptr_ori, line->unpack_size[1], &tmp);
            if ( !line->rlen && line->n_allele ) line->rlen = strlen(line->d.allele[0]);
        }
        ptr_ori += line->unpack_size[1];
        line->unpack_size[1] = tmp.l - prev_len; prev_len = tmp.l;

        if ( line->unpacked & BCF_UN_FLT )
        {
            // FILTER: typed vector of integers
            if ( line->d.shared_dirty & BCF1_DIRTY_FLT )
                bcf1_sync_filter(line, &tmp);
            else if ( line->d.n_flt )
                kputsn_(ptr_ori, line->unpack_size[2], &tmp);
            else
                bcf_enc_vint(&tmp, 0, 0, -1);
            ptr_ori += line->unpack_size[2];
            line->unpack_size[2] = tmp.l - prev_len;

            if ( line->unpacked & BCF_UN_INFO )
            {
                // INFO: pairs of typed vectors
                if ( line->d.shared_dirty & BCF1_DIRTY_INF )
                {
                    bcf1_sync_info(line, &tmp);
                    ptr_ori = (uint8_t*)line->shared.s + line->shared.l;
                }
            }
        }

        int size = line->shared.l - (size_t)ptr_ori + (size_t)line->shared.s;
        if ( size ) kputsn_(ptr_ori, size, &tmp);

        free(line->shared.s);
        line->shared = tmp;
    }
    if ( line->shared.s != shared_ori && line->unpacked & BCF_UN_INFO )
    {
        // Reallocated line->shared.s block invalidated line->d.info[].vptr pointers
        size_t off_new = line->unpack_size[0] + line->unpack_size[1] + line->unpack_size[2];
        int i;
        for (i=0; i<line->n_info; i++)
        {
            uint8_t *vptr_free = line->d.info[i].vptr_free ? line->d.info[i].vptr - line->d.info[i].vptr_off : NULL;
            line->d.info[i].vptr = (uint8_t*) line->shared.s + off_new + line->d.info[i].vptr_off;
            off_new += line->d.info[i].vptr_len + line->d.info[i].vptr_off;
            if ( vptr_free )
            {
                free(vptr_free);
                line->d.info[i].vptr_free = 0;
            }
        }
    }

    if ( line->n_sample && line->n_fmt && (!line->indiv.l || line->d.indiv_dirty) )
    {
        // The genotype fields changed or are not present
        tmp.l = tmp.m = 0; tmp.s = NULL;
        int i, irm = -1;
        for (i=0; i<line->n_fmt; i++)
        {
            bcf_fmt_t *fmt = &line->d.fmt[i];
            if ( !fmt->p )
            {
                // marked for removal
                if ( irm < 0 ) irm = i;
                continue;
            }
            kputsn_(fmt->p - fmt->p_off, fmt->p_len + fmt->p_off, &tmp);
            if ( irm >=0 )
            {
                bcf_fmt_t tfmt = line->d.fmt[irm]; line->d.fmt[irm] = line->d.fmt[i]; line->d.fmt[i] = tfmt;
                while ( irm<=i && line->d.fmt[irm].p ) irm++;
            }

        }
        if ( irm>=0 ) line->n_fmt = irm;
        free(line->indiv.s);
        line->indiv = tmp;

        // Reallocated line->indiv.s block invalidated line->d.fmt[].p pointers
        size_t off_new = 0;
        for (i=0; i<line->n_fmt; i++)
        {
            uint8_t *p_free = line->d.fmt[i].p_free ? line->d.fmt[i].p - line->d.fmt[i].p_off : NULL;
            line->d.fmt[i].p = (uint8_t*) line->indiv.s + off_new + line->d.fmt[i].p_off;
            off_new += line->d.fmt[i].p_len + line->d.fmt[i].p_off;
            if ( p_free )
            {
                free(p_free);
                line->d.fmt[i].p_free = 0;
            }
        }
    }
    if ( !line->n_sample ) line->n_fmt = 0;
    line->d.shared_dirty = line->d.indiv_dirty = 0;
    return 0;
}

bcf1_t *bcf_copy(bcf1_t *dst, bcf1_t *src)
{
    bcf1_sync(src);

    bcf_clear(dst);
    dst->rid  = src->rid;
    dst->pos  = src->pos;
    dst->rlen = src->rlen;
    dst->qual = src->qual;
    dst->n_info = src->n_info; dst->n_allele = src->n_allele;
    dst->n_fmt = src->n_fmt; dst->n_sample = src->n_sample;

    if ( dst->shared.m < src->shared.l )
    {
        dst->shared.s = (char*) realloc(dst->shared.s, src->shared.l);
        dst->shared.m = src->shared.l;
    }
    dst->shared.l = src->shared.l;
    memcpy(dst->shared.s,src->shared.s,dst->shared.l);

    if ( dst->indiv.m < src->indiv.l )
    {
        dst->indiv.s = (char*) realloc(dst->indiv.s, src->indiv.l);
        dst->indiv.m = src->indiv.l;
    }
    dst->indiv.l = src->indiv.l;
    memcpy(dst->indiv.s,src->indiv.s,dst->indiv.l);

    return dst;
}
bcf1_t *bcf_dup(bcf1_t *src)
{
    bcf1_t *out = bcf_init1();
    return bcf_copy(out, src);
}

int bcf_write(htsFile *hfp, bcf_hdr_t *h, bcf1_t *v)
{
    if ( h->dirty ) {
        if (bcf_hdr_sync(h) < 0) return -1;
    }
    if ( bcf_hdr_nsamples(h)!=v->n_sample )
    {
        hts_log_error("Broken VCF record, the number of columns at %s:%"PRIhts_pos" does not match the number of samples (%d vs %d)",
            bcf_seqname_safe(h,v), v->pos+1, v->n_sample, bcf_hdr_nsamples(h));
        return -1;
    }

    if ( hfp->format.format == vcf || hfp->format.format == text_format )
        return vcf_write(hfp,h,v);

    if ( v->errcode )
    {
        // vcf_parse1() encountered a new contig or tag, undeclared in the
        // header.  At this point, the header must have been printed,
        // proceeding would lead to a broken BCF file. Errors must be checked
        // and cleared by the caller before we can proceed.
        hts_log_error("Unchecked error (%d) at %s:%"PRIhts_pos, v->errcode, bcf_seqname_safe(h,v), v->pos+1);
        return -1;
    }
    bcf1_sync(v);   // check if the BCF record was modified

    if ( v->unpacked & BCF_IS_64BIT )
    {
        hts_log_error("Data at %s:%"PRIhts_pos" contains 64-bit values not representable in BCF. Please use VCF instead", bcf_seqname_safe(h,v), v->pos+1);
        return -1;
    }

    BGZF *fp = hfp->fp.bgzf;
    uint8_t x[32];
    u32_to_le(v->shared.l + 24, x); // to include six 32-bit integers
    u32_to_le(v->indiv.l, x + 4);
    i32_to_le(v->rid, x + 8);
    u32_to_le(v->pos, x + 12);
    u32_to_le(v->rlen, x + 16);
    float_to_le(v->qual, x + 20);
    u16_to_le(v->n_info, x + 24);
    u16_to_le(v->n_allele, x + 26);
    u32_to_le((uint32_t)v->n_fmt<<24 | (v->n_sample & 0xffffff), x + 28);
    if ( bgzf_write(fp, x, 32) != 32 ) return -1;
    if ( bgzf_write(fp, v->shared.s, v->shared.l) != v->shared.l ) return -1;
    if ( bgzf_write(fp, v->indiv.s, v->indiv.l) != v->indiv.l ) return -1;

    if (hfp->idx) {
        if (hts_idx_push(hfp->idx, v->rid, v->pos, v->pos + v->rlen, bgzf_tell(fp), 1) < 0)
            return -1;
    }

    return 0;
}

/**********************
 *** VCF header I/O ***
 **********************/

static int add_missing_contig_hrec(bcf_hdr_t *h, const char *name) {
    bcf_hrec_t *hrec = calloc(1, sizeof(bcf_hrec_t));
    int save_errno;
    if (!hrec) goto fail;

    hrec->key = strdup("contig");
    if (!hrec->key) goto fail;

    if (bcf_hrec_add_key(hrec, "ID", strlen("ID")) < 0) goto fail;
    if (bcf_hrec_set_val(hrec, hrec->nkeys-1, name, strlen(name), 0) < 0)
        goto fail;
    if (bcf_hdr_add_hrec(h, hrec) < 0)
        goto fail;
    return 0;

 fail:
    save_errno = errno;
    hts_log_error("%s", strerror(errno));
    if (hrec) bcf_hrec_destroy(hrec);
    errno = save_errno;
    return -1;
}

bcf_hdr_t *vcf_hdr_read(htsFile *fp)
{
    kstring_t txt, *s = &fp->line;
    int ret;
    bcf_hdr_t *h;
    tbx_t *idx = NULL;
    const char **names = NULL;
    h = bcf_hdr_init("r");
    if (!h) {
        hts_log_error("Failed to allocate bcf header");
        return NULL;
    }
    txt.l = txt.m = 0; txt.s = 0;
    while ((ret = hts_getline(fp, KS_SEP_LINE, s)) >= 0) {
        int e = 0;
        if (s->l == 0) continue;
        if (s->s[0] != '#') {
            hts_log_error("No sample line");
            goto error;
        }
        if (s->s[1] != '#' && fp->fn_aux) { // insert contigs here
            kstring_t tmp = { 0, 0, NULL };
            hFILE *f = hopen(fp->fn_aux, "r");
            if (f == NULL) {
                hts_log_error("Couldn't open \"%s\"", fp->fn_aux);
                goto error;
            }
            while (tmp.l = 0, kgetline(&tmp, (kgets_func *) hgets, f) >= 0) {
                char *tab = strchr(tmp.s, '\t');
                if (tab == NULL) continue;
                e |= (kputs("##contig=<ID=", &txt) < 0);
                e |= (kputsn(tmp.s, tab - tmp.s, &txt) < 0);
                e |= (kputs(",length=", &txt) < 0);
                e |= (kputl(atol(tab), &txt) < 0);
                e |= (kputsn(">\n", 2, &txt) < 0);
            }
            free(tmp.s);
            if (hclose(f) != 0) {
                hts_log_error("Error on closing %s", fp->fn_aux);
                goto error;
            }
            if (e) goto error;
        }
        if (kputsn(s->s, s->l, &txt) < 0) goto error;
        if (kputc('\n', &txt) < 0) goto error;
        if (s->s[1] != '#') break;
    }
    if ( ret < -1 ) goto error;
    if ( !txt.s )
    {
        hts_log_error("Could not read the header");
        goto error;
    }
    if ( bcf_hdr_parse(h, txt.s) < 0 ) goto error;

    // check tabix index, are all contigs listed in the header? add the missing ones
    idx = tbx_index_load3(fp->fn, NULL, HTS_IDX_SAVE_REMOTE|HTS_IDX_SILENT_FAIL);
    if ( idx )
    {
        int i, n, need_sync = 0;
        names = tbx_seqnames(idx, &n);
        if (!names) goto error;
        for (i=0; i<n; i++)
        {
            bcf_hrec_t *hrec = bcf_hdr_get_hrec(h, BCF_HL_CTG, "ID", (char*) names[i], NULL);
            if ( hrec ) continue;
            if (add_missing_contig_hrec(h, names[i]) < 0) goto error;
            need_sync = 1;
        }
        if ( need_sync ) {
            if (bcf_hdr_sync(h) < 0) goto error;
        }
        free(names);
        tbx_destroy(idx);
    }
    free(txt.s);
    return h;

 error:
    if (idx) tbx_destroy(idx);
    free(names);
    free(txt.s);
    if (h) bcf_hdr_destroy(h);
    return NULL;
}

int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname)
{
    int i = 0, n = 0, save_errno;
    char **lines = hts_readlines(fname, &n);
    if ( !lines ) return 1;
    for (i=0; i<n-1; i++)
    {
        int k;
        bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr,lines[i],&k);
        if (!hrec) goto fail;
        if (bcf_hdr_add_hrec(hdr, hrec) < 0) {
            bcf_hrec_destroy(hrec);
            goto fail;
        }
        free(lines[i]);
        lines[i] = NULL;
    }
    if (bcf_hdr_parse_sample_line(hdr, lines[n-1]) < 0) goto fail;
    if (bcf_hdr_sync(hdr) < 0) goto fail;
    free(lines[n-1]);
    free(lines);
    return 0;

 fail:
    save_errno = errno;
    for (; i < n; i++)
        free(lines[i]);
    free(lines);
    errno = save_errno;
    return 1;
}

static int _bcf_hrec_format(const bcf_hrec_t *hrec, int is_bcf, kstring_t *str)
{
    uint32_t e = 0;
    if ( !hrec->value )
    {
        int j, nout = 0;
        e |= ksprintf(str, "##%s=<", hrec->key) < 0;
        for (j=0; j<hrec->nkeys; j++)
        {
            // do not output IDX if output is VCF
            if ( !is_bcf && !strcmp("IDX",hrec->keys[j]) ) continue;
            if ( nout ) e |= kputc(',',str) < 0;
            e |= ksprintf(str,"%s=%s", hrec->keys[j], hrec->vals[j]) < 0;
            nout++;
        }
        e |= ksprintf(str,">\n") < 0;
    }
    else
        e |= ksprintf(str,"##%s=%s\n", hrec->key,hrec->value) < 0;

    return e == 0 ? 0 : -1;
}

int bcf_hrec_format(const bcf_hrec_t *hrec, kstring_t *str)
{
    return _bcf_hrec_format(hrec,0,str);
}

int bcf_hdr_format(const bcf_hdr_t *hdr, int is_bcf, kstring_t *str)
{
    int i, r = 0;
    for (i=0; i<hdr->nhrec; i++)
        r |= _bcf_hrec_format(hdr->hrec[i], is_bcf, str) < 0;

    r |= ksprintf(str, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") < 0;
    if ( bcf_hdr_nsamples(hdr) )
    {
        r |= ksprintf(str, "\tFORMAT") < 0;
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
            r |= ksprintf(str, "\t%s", hdr->samples[i]) < 0;
    }
    r |= ksprintf(str, "\n") < 0;

    return r ? -1 : 0;
}

char *bcf_hdr_fmt_text(const bcf_hdr_t *hdr, int is_bcf, int *len)
{
    kstring_t txt = {0,0,0};
    if (bcf_hdr_format(hdr, is_bcf, &txt) < 0)
        return NULL;
    if ( len ) *len = txt.l;
    return txt.s;
}

const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *n)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    khint_t k;
    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        assert( tid<m );
        names[tid] = kh_key(d,k);
    }
    // sanity check: there should be no gaps
    for (tid=0; tid<m; tid++)
        assert(names[tid]);
    *n = m;
    return names;
}

int vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h)
{
    kstring_t htxt = {0,0,0};
    if (bcf_hdr_format(h, 0, &htxt) < 0) {
        free(htxt.s);
        return -1;
    }
    while (htxt.l && htxt.s[htxt.l-1] == '\0') --htxt.l; // kill trailing zeros
    int ret;
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, htxt.s, htxt.l);
    else
        ret = hwrite(fp->fp.hfile, htxt.s, htxt.l);
    free(htxt.s);
    return ret<0 ? -1 : 0;
}

/***********************
 *** Typed value I/O ***
 ***********************/

int bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
{
    int32_t max = INT32_MIN, min = INT32_MAX;
    int i;
    if (n <= 0) bcf_enc_size(s, 0, BCF_BT_NULL);
    else if (n == 1) bcf_enc_int1(s, a[0]);
    else {
        if (wsize <= 0) wsize = n;
        for (i = 0; i < n; ++i) {
            if (a[i] == bcf_int32_missing || a[i] == bcf_int32_vector_end ) continue;
            if (max < a[i]) max = a[i];
            if (min > a[i]) min = a[i];
        }
        if (max <= BCF_MAX_BT_INT8 && min >= BCF_MIN_BT_INT8) {
            bcf_enc_size(s, wsize, BCF_BT_INT8);
            for (i = 0; i < n; ++i)
                if ( a[i]==bcf_int32_vector_end ) kputc(bcf_int8_vector_end, s);
                else if ( a[i]==bcf_int32_missing ) kputc(bcf_int8_missing, s);
                else kputc(a[i], s);
        } else if (max <= BCF_MAX_BT_INT16 && min >= BCF_MIN_BT_INT16) {
            uint8_t *p;
            bcf_enc_size(s, wsize, BCF_BT_INT16);
            ks_resize(s, s->l + n * sizeof(int16_t));
            p = (uint8_t *) s->s + s->l;
            for (i = 0; i < n; ++i)
            {
                int16_t x;
                if ( a[i]==bcf_int32_vector_end ) x = bcf_int16_vector_end;
                else if ( a[i]==bcf_int32_missing ) x = bcf_int16_missing;
                else x = a[i];
                i16_to_le(x, p);
                p += sizeof(int16_t);
            }
            s->l += n * sizeof(int16_t);
        } else {
            uint8_t *p;
            bcf_enc_size(s, wsize, BCF_BT_INT32);
            ks_resize(s, s->l + n * sizeof(int32_t));
            p = (uint8_t *) s->s + s->l;
            for (i = 0; i < n; ++i) {
                i32_to_le(a[i], p);
                p += sizeof(int32_t);
            }
            s->l += n * sizeof(int32_t);
        }
    }

    return 0; // FIXME: check for errs in this function
}

#ifdef VCF_ALLOW_INT64
static int bcf_enc_long1(kstring_t *s, int64_t x) {
    uint32_t e = 0;
    if (x <= BCF_MAX_BT_INT32 && x >= BCF_MIN_BT_INT32)
        return bcf_enc_int1(s, x);
    if (x == bcf_int64_vector_end) {
        e |= bcf_enc_size(s, 1, BCF_BT_INT8);
        e |= kputc(bcf_int8_vector_end, s) < 0;
    } else if (x == bcf_int64_missing) {
        e |= bcf_enc_size(s, 1, BCF_BT_INT8);
        e |= kputc(bcf_int8_missing, s) < 0;
    } else {
        e |= bcf_enc_size(s, 1, BCF_BT_INT64);
        e |= ks_expand(s, 8);
        if (e == 0) { u64_to_le(x, (uint8_t *) s->s + s->l); s->l += 8; }
    }
    return e == 0 ? 0 : -1;
}
#endif

static inline int serialize_float_array(kstring_t *s, size_t n, const float *a) {
    uint8_t *p;
    size_t i;
    size_t bytes = n * sizeof(float);

    if (bytes / sizeof(float) != n) return -1;
    if (ks_resize(s, s->l + bytes) < 0) return -1;

    p = (uint8_t *) s->s + s->l;
    for (i = 0; i < n; i++) {
        float_to_le(a[i], p);
        p += sizeof(float);
    }
    s->l += bytes;

    return 0;
}

int bcf_enc_vfloat(kstring_t *s, int n, float *a)
{
    assert(n >= 0);
    bcf_enc_size(s, n, BCF_BT_FLOAT);
    serialize_float_array(s, n, a);
    return 0; // FIXME: check for errs in this function
}

int bcf_enc_vchar(kstring_t *s, int l, const char *a)
{
    bcf_enc_size(s, l, BCF_BT_CHAR);
    kputsn(a, l, s);
    return 0; // FIXME: check for errs in this function
}

int bcf_fmt_array(kstring_t *s, int n, int type, void *data)
{
    int j = 0;
    uint32_t e = 0;
    if (n == 0) {
        return kputc('.', s) >= 0 ? 0 : -1;
    }
    if (type == BCF_BT_CHAR)
    {
        char *p = (char*)data;
        for (j = 0; j < n && *p; ++j, ++p)
        {
            if ( *p==bcf_str_missing ) e |= kputc('.', s) < 0;
            else e |= kputc(*p, s) < 0;
        }
    }
    else
    {
        #define BRANCH(type_t, convert, is_missing, is_vector_end, kprint) { \
            uint8_t *p = (uint8_t *) data; \
            for (j=0; j<n; j++, p += sizeof(type_t))    \
            { \
                type_t v = convert(p); \
                if ( is_vector_end ) break; \
                if ( j ) kputc(',', s); \
                if ( is_missing ) kputc('.', s); \
                else e |= kprint < 0; \
            } \
        }
        switch (type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, v==bcf_int8_missing,  v==bcf_int8_vector_end,  kputw(v, s)); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, v==bcf_int16_missing, v==bcf_int16_vector_end, kputw(v, s)); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, v==bcf_int32_missing, v==bcf_int32_vector_end, kputw(v, s)); break;
            case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, v==bcf_float_missing, v==bcf_float_vector_end, kputd(le_to_float(p), s)); break;
            default: hts_log_error("Unexpected type %d", type); exit(1); break;
        }
        #undef BRANCH
    }
    return e == 0 ? 0 : -1;
}

uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr)
{
    int x, type;
    x = bcf_dec_size(ptr, &ptr, &type);
    bcf_fmt_array(s, x, type, ptr);
    return ptr + (x << bcf_type_shift[type]);
}

/********************
 *** VCF site I/O ***
 ********************/

typedef struct {
    int key, max_m, size, offset;
    uint32_t is_gt:1, max_g:31;
    uint32_t max_l;
    uint32_t y;
    uint8_t *buf;
} fmt_aux_t;

static inline int align_mem(kstring_t *s)
{
    int e = 0;
    if (s->l&7) {
        uint64_t zero = 0;
        e = kputsn((char*)&zero, 8 - (s->l&7), s) < 0;
    }
    return e == 0 ? 0 : -1;
}

// p,q is the start and the end of the FORMAT field
#define MAX_N_FMT 255   /* Limited by size of bcf1_t n_fmt field */
static int vcf_parse_format(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v, char *p, char *q)
{
    if ( !bcf_hdr_nsamples(h) ) return 0;

    static int extreme_val_warned = 0;
    char *r, *t;
    int j, l, m, g, overflow = 0;
    khint_t k;
    ks_tokaux_t aux1;
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
    kstring_t *mem = (kstring_t*)&h->mem;
    fmt_aux_t fmt[MAX_N_FMT];
    mem->l = 0;

    char *end = s->s + s->l;
    if ( q>=end )
    {
        hts_log_error("FORMAT column with no sample columns starting at %s:%"PRIhts_pos"", bcf_seqname_safe(h,v), v->pos+1);
        v->errcode |= BCF_ERR_NCOLS;
        return -1;
    }

    v->n_fmt = 0;
    if ( p[0]=='.' && p[1]==0 ) // FORMAT field is empty "."
    {
        v->n_sample = bcf_hdr_nsamples(h);
        return 0;
    }

    // get format information from the dictionary
    for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
        if (j >= MAX_N_FMT) {
            v->errcode |= BCF_ERR_LIMITS;
            hts_log_error("FORMAT column at %s:%"PRIhts_pos" lists more identifiers than htslib can handle",
                bcf_seqname_safe(h,v), v->pos+1);
            return -1;
        }

        *(char*)aux1.p = 0;
        k = kh_get(vdict, d, t);
        if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_FMT] == 15) {
            if ( t[0]=='.' && t[1]==0 )
            {
                hts_log_error("Invalid FORMAT tag name '.' at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_TAG_INVALID;
                return -1;
            }
            hts_log_warning("FORMAT '%s' at %s:%"PRIhts_pos" is not defined in the header, assuming Type=String", t, bcf_seqname_safe(h,v), v->pos+1);
            kstring_t tmp = {0,0,0};
            int l;
            ksprintf(&tmp, "##FORMAT=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", t);
            bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
            free(tmp.s);
            int res = hrec ? bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) : -1;
            if (res < 0) bcf_hrec_destroy(hrec);
            if (res > 0) res = bcf_hdr_sync((bcf_hdr_t*)h);

            k = kh_get(vdict, d, t);
            v->errcode = BCF_ERR_TAG_UNDEF;
            if (res || k == kh_end(d)) {
                hts_log_error("Could not add dummy header for FORMAT '%s' at %s:%"PRIhts_pos, t, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_TAG_INVALID;
                return -1;
            }
        }
        fmt[j].max_l = fmt[j].max_m = fmt[j].max_g = 0;
        fmt[j].key = kh_val(d, k).id;
        fmt[j].is_gt = !strcmp(t, "GT");
        fmt[j].y = h->id[0][fmt[j].key].val->info[BCF_HL_FMT];
        v->n_fmt++;
    }
    // compute max
    int n_sample_ori = -1;
    r = q + 1;  // r: position in the format string
    l = 0, m = g = 1, v->n_sample = 0;  // m: max vector size, l: max field len, g: max number of alleles
    while ( r<end )
    {
        // can we skip some samples?
        if ( h->keep_samples )
        {
            n_sample_ori++;
            if ( !bit_array_test(h->keep_samples,n_sample_ori) )
            {
                while ( *r!='\t' && r<end ) r++;
                if ( *r=='\t' ) { *r = 0; r++; }
                continue;
            }
        }

        // collect fmt stats: max vector size, length, number of alleles
        j = 0;  // j-th format field
        fmt_aux_t *f = fmt;
        for (;;) {
            switch (*r) {
            case ',':
                m++;
                break;

            case '|':
            case '/':
                if (f->is_gt) g++;
                break;

            case '\t':
                *r = 0; // fall through

            case '\0':
            case ':':
                if (f->max_m < m) f->max_m = m;
                if (f->max_l < l) f->max_l = l;
                if (f->is_gt && f->max_g < g) f->max_g = g;
                l = 0, m = g = 1;
                if ( *r==':' ) {
                    j++; f++;
                    if ( j>=v->n_fmt ) {
                        hts_log_error("Incorrect number of FORMAT fields at %s:%"PRIhts_pos"",
                                      h->id[BCF_DT_CTG][v->rid].key, v->pos+1);
                        v->errcode |= BCF_ERR_NCOLS;
                        return -1;
                    }
                } else goto end_for;
                break;
            }
            if ( r>=end ) break;
            r++; l++;
        }
    end_for:
        v->n_sample++;
        if ( v->n_sample == bcf_hdr_nsamples(h) ) break;
        r++;
    }

    // allocate memory for arrays
    for (j = 0; j < v->n_fmt; ++j) {
        fmt_aux_t *f = &fmt[j];
        if ( !f->max_m ) f->max_m = 1;  // omitted trailing format field
        if ((f->y>>4&0xf) == BCF_HT_STR) {
            f->size = f->is_gt? f->max_g << 2 : f->max_l;
        } else if ((f->y>>4&0xf) == BCF_HT_REAL || (f->y>>4&0xf) == BCF_HT_INT) {
            f->size = f->max_m << 2;
        } else
        {
            hts_log_error("The format type %d at %s:%"PRIhts_pos" is currently not supported", f->y>>4&0xf, bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_TAG_INVALID;
            return -1;
        }
        if (align_mem(mem) < 0) {
            hts_log_error("Memory allocation failure at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_LIMITS;
            return -1;
        }

        // Limit the total memory to ~2Gb per VCF row.  This should mean
        // malformed VCF data is less likely to take excessive memory and/or
        // time.
        if ((uint64_t) mem->l + v->n_sample * (uint64_t)f->size > INT_MAX) {
            hts_log_error("Excessive memory required by FORMAT fields at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_LIMITS;
            return -1;
        }

        f->offset = mem->l;
        if (ks_resize(mem, mem->l + v->n_sample * (size_t)f->size) < 0) {
            hts_log_error("Memory allocation failure at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_LIMITS;
            return -1;
        }
        mem->l += v->n_sample * f->size;
    }
    for (j = 0; j < v->n_fmt; ++j)
        fmt[j].buf = (uint8_t*)mem->s + fmt[j].offset;
    // fill the sample fields; at beginning of the loop, t points to the first char of a format
    n_sample_ori = -1;
    t = q + 1; m = 0;   // m: sample id
    while ( t<end )
    {
        // can we skip some samples?
        if ( h->keep_samples )
        {
            n_sample_ori++;
            if ( !bit_array_test(h->keep_samples,n_sample_ori) )
            {
                while ( *t && t<end ) t++;
                t++;
                continue;
            }
        }
        if ( m == bcf_hdr_nsamples(h) ) break;

        j = 0; // j-th format field, m-th sample
        while ( t < end )
        {
            fmt_aux_t *z = &fmt[j++];
            if (!z->buf) {
                hts_log_error("Memory allocation failure for FORMAT field type %d at %s:%"PRIhts_pos,
                              z->y>>4&0xf, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_LIMITS;
                return -1;
            }
            if ((z->y>>4&0xf) == BCF_HT_STR) {
                if (z->is_gt) { // genotypes
                    int32_t is_phased = 0;
                    uint32_t *x = (uint32_t*)(z->buf + z->size * (size_t)m);
                    uint32_t unreadable = 0;
                    uint32_t max = 0;
                    overflow = 0;
                    for (l = 0;; ++t) {
                        if (*t == '.') {
                            ++t, x[l++] = is_phased;
                        } else {
                            char *tt = t;
                            uint32_t val = hts_str2uint(t, &t, sizeof(val) * CHAR_MAX - 2, &overflow);
                            unreadable |= tt == t;
                            if (max < val) max = val;
                            x[l++] = (val + 1) << 1 | is_phased;
                        }
                        is_phased = (*t == '|');
                        if (*t != '|' && *t != '/') break;
                    }
                    // Possibly check max against v->n_allele instead?
                    if (overflow || max > (INT32_MAX >> 1) - 1) {
                        hts_log_error("Couldn't read GT data: value too large at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
                        return -1;
                    }
                    if (unreadable) {
                        hts_log_error("Couldn't read GT data: value not a number or '.' at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
                        return -1;
                    }
                    if ( !l ) x[l++] = 0;   // An empty field, insert missing value
                    for (; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
                } else {
                    char *x = (char*)z->buf + z->size * (size_t)m;
                    for (r = t, l = 0; *t != ':' && *t; ++t) x[l++] = *t;
                    for (; l < z->size; ++l) x[l] = 0;
                }
            } else if ((z->y>>4&0xf) == BCF_HT_INT) {
                int32_t *x = (int32_t*)(z->buf + z->size * (size_t)m);
                for (l = 0;; ++t) {
                    if (*t == '.') {
                        x[l++] = bcf_int32_missing, ++t; // ++t to skip "."
                    } else {
                        overflow = 0;
                        char *te;
                        long int tmp_val = hts_str2int(t, &te, sizeof(tmp_val)*CHAR_BIT, &overflow);
                        if ( te==t || overflow || tmp_val<BCF_MIN_BT_INT32 || tmp_val>BCF_MAX_BT_INT32 )
                        {
                            if ( !extreme_val_warned )
                            {
                                hts_log_warning("Extreme FORMAT/%s value encountered and set to missing at %s:%"PRIhts_pos, h->id[BCF_DT_ID][fmt[j-1].key].key, bcf_seqname_safe(h,v), v->pos+1);
                                extreme_val_warned = 1;
                            }
                            tmp_val = bcf_int32_missing;
                        }
                        x[l++] = tmp_val;
                        t = te;
                    }
                    if (*t != ',') break;
                }
                if ( !l ) x[l++] = bcf_int32_missing;
                for (; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
            } else if ((z->y>>4&0xf) == BCF_HT_REAL) {
                float *x = (float*)(z->buf + z->size * (size_t)m);
                for (l = 0;; ++t) {
                    if (*t == '.' && !isdigit_c(t[1])) {
                        bcf_float_set_missing(x[l++]), ++t; // ++t to skip "."
                    } else {
                        overflow = 0;
                        char *te;
                        float tmp_val = hts_str2dbl(t, &te, &overflow);
                        if ( (te==t || overflow) && !extreme_val_warned )
                        {
                            hts_log_warning("Extreme FORMAT/%s value encountered at %s:%"PRIhts_pos, h->id[BCF_DT_ID][fmt[j-1].key].key, bcf_seqname(h,v), v->pos+1);
                            extreme_val_warned = 1;
                        }
                        x[l++] = tmp_val;
                        t = te;
                    }
                    if (*t != ',') break;
                }
                if ( !l ) bcf_float_set_missing(x[l++]);    // An empty field, insert missing value
                for (; l < z->size>>2; ++l) bcf_float_set_vector_end(x[l]);
            } else {
                hts_log_error("Unknown FORMAT field type %d at %s:%"PRIhts_pos, z->y>>4&0xf, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_TAG_INVALID;
                return -1;
            }

            if (*t == '\0') {
                break;
            }
            else if (*t == ':') {
                t++;
            }
            else {
                char buffer[8];
                hts_log_error("Invalid character %s in '%s' FORMAT field at %s:%"PRIhts_pos"",
                    hts_strprint(buffer, sizeof buffer, '\'', t, 1),
                    h->id[BCF_DT_ID][z->key].key, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_CHAR;
                return -1;
            }
        }

        for (; j < v->n_fmt; ++j) { // fill end-of-vector values
            fmt_aux_t *z = &fmt[j];
            if ((z->y>>4&0xf) == BCF_HT_STR) {
                if (z->is_gt) {
                    int32_t *x = (int32_t*)(z->buf + z->size * (size_t)m);
                    if (z->size) x[0] = bcf_int32_missing;
                    for (l = 1; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
                } else {
                    char *x = (char*)z->buf + z->size * (size_t)m;
                    if ( z->size ) x[0] = '.';
                    for (l = 1; l < z->size; ++l) x[l] = 0;
                }
            } else if ((z->y>>4&0xf) == BCF_HT_INT) {
                int32_t *x = (int32_t*)(z->buf + z->size * (size_t)m);
                x[0] = bcf_int32_missing;
                for (l = 1; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
            } else if ((z->y>>4&0xf) == BCF_HT_REAL) {
                float *x = (float*)(z->buf + z->size * (size_t)m);
                bcf_float_set_missing(x[0]);
                for (l = 1; l < z->size>>2; ++l) bcf_float_set_vector_end(x[l]);
            }
        }

        m++; t++;
    }

    // write individual genotype information
    kstring_t *str = &v->indiv;
    int i;
    if (v->n_sample > 0) {
        for (i = 0; i < v->n_fmt; ++i) {
            fmt_aux_t *z = &fmt[i];
            bcf_enc_int1(str, z->key);
            if ((z->y>>4&0xf) == BCF_HT_STR && !z->is_gt) {
                bcf_enc_size(str, z->size, BCF_BT_CHAR);
                kputsn((char*)z->buf, z->size * (size_t)v->n_sample, str);
            } else if ((z->y>>4&0xf) == BCF_HT_INT || z->is_gt) {
                bcf_enc_vint(str, (z->size>>2) * v->n_sample, (int32_t*)z->buf, z->size>>2);
            } else {
                bcf_enc_size(str, z->size>>2, BCF_BT_FLOAT);
                if (serialize_float_array(str, (z->size>>2) * (size_t)v->n_sample,
                                          (float *) z->buf) != 0) {
                    v->errcode |= BCF_ERR_LIMITS;
                    hts_log_error("Out of memory at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
                    return -1;
                }
            }
        }
    }

    if ( v->n_sample!=bcf_hdr_nsamples(h) )
    {
        hts_log_error("Number of columns at %s:%"PRIhts_pos" does not match the number of samples (%d vs %d)",
            bcf_seqname_safe(h,v), v->pos+1, v->n_sample, bcf_hdr_nsamples(h));
        v->errcode |= BCF_ERR_NCOLS;
        return -1;
    }
    if ( v->indiv.l > 0xffffffff )
    {
        hts_log_error("The FORMAT at %s:%"PRIhts_pos" is too long", bcf_seqname_safe(h,v), v->pos+1);
        v->errcode |= BCF_ERR_LIMITS;

        // Error recovery: return -1 if this is a critical error or 0 if we want to ignore the FORMAT and proceed
        v->n_fmt = 0;
        return -1;
    }

    return 0;
}

static khint_t fix_chromosome(const bcf_hdr_t *h, vdict_t *d, const char *p) {
    // Simple error recovery for chromosomes not defined in the header. It will not help when VCF header has
    // been already printed, but will enable tools like vcfcheck to proceed.

    kstring_t tmp = {0,0,0};
    khint_t k;
    int l;
    if (ksprintf(&tmp, "##contig=<ID=%s>", p) < 0)
        return kh_end(d);
    bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
    free(tmp.s);
    int res = hrec ? bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) : -1;
    if (res < 0) bcf_hrec_destroy(hrec);
    if (res > 0) res = bcf_hdr_sync((bcf_hdr_t*)h);
    k = kh_get(vdict, d, p);

    return k;
}

static int vcf_parse_filter(kstring_t *str, const bcf_hdr_t *h, bcf1_t *v, char *p, char *q) {
    int i, n_flt = 1, max_n_flt = 0;
    char *r, *t;
    int32_t *a_flt = NULL;
    ks_tokaux_t aux1;
    khint_t k;
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
    // count the number of filters
    if (*(q-1) == ';') *(q-1) = 0;
    for (r = p; *r; ++r)
        if (*r == ';') ++n_flt;
    if (n_flt > max_n_flt) {
        a_flt = malloc(n_flt * sizeof(*a_flt));
        if (!a_flt) {
            hts_log_error("Could not allocate memory at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_LIMITS; // No appropriate code?
            return -1;
        }
        max_n_flt = n_flt;
    }
    // add filters
    for (t = kstrtok(p, ";", &aux1), i = 0; t; t = kstrtok(0, 0, &aux1)) {
        *(char*)aux1.p = 0;
        k = kh_get(vdict, d, t);
        if (k == kh_end(d))
        {
            // Simple error recovery for FILTERs not defined in the header. It will not help when VCF header has
            // been already printed, but will enable tools like vcfcheck to proceed.
            hts_log_warning("FILTER '%s' is not defined in the header", t);
            kstring_t tmp = {0,0,0};
            int l;
            ksprintf(&tmp, "##FILTER=<ID=%s,Description=\"Dummy\">", t);
            bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
            free(tmp.s);
            int res = hrec ? bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) : -1;
            if (res < 0) bcf_hrec_destroy(hrec);
            if (res > 0) res = bcf_hdr_sync((bcf_hdr_t*)h);
            k = kh_get(vdict, d, t);
            v->errcode |= BCF_ERR_TAG_UNDEF;
            if (res || k == kh_end(d)) {
                hts_log_error("Could not add dummy header for FILTER '%s' at %s:%"PRIhts_pos, t, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_TAG_INVALID;
                free(a_flt);
                return -1;
            }
        }
        a_flt[i++] = kh_val(d, k).id;
    }

    bcf_enc_vint(str, n_flt, a_flt, -1);
    free(a_flt);

    return 0;
}

static int vcf_parse_info(kstring_t *str, const bcf_hdr_t *h, bcf1_t *v, char *p, char *q) {
    static int extreme_int_warned = 0, negative_rlen_warned = 0;
    int max_n_val = 0, overflow = 0;
    char *r, *key;
    khint_t k;
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
    int32_t *a_val = NULL;

    v->n_info = 0;
    if (*(q-1) == ';') *(q-1) = 0;
    for (r = key = p;; ++r) {
        int c;
        char *val, *end;
        if (*r != ';' && *r != '=' && *r != 0) continue;
        if (v->n_info == UINT16_MAX) {
            hts_log_error("Too many INFO entries at %s:%"PRIhts_pos,
                          bcf_seqname_safe(h,v), v->pos+1);
            v->errcode |= BCF_ERR_LIMITS;
            return -1;
        }
        val = end = 0;
        c = *r; *r = 0;
        if (c == '=') {
            val = r + 1;
            for (end = val; *end != ';' && *end != 0; ++end);
            c = *end; *end = 0;
        } else end = r;
        if ( !*key ) { if (c==0) break; r = end; key = r + 1; continue; }  // faulty VCF, ";;" in the INFO
        k = kh_get(vdict, d, key);
        if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_INFO] == 15)
        {
            hts_log_warning("INFO '%s' is not defined in the header, assuming Type=String", key);
            kstring_t tmp = {0,0,0};
            int l;
            ksprintf(&tmp, "##INFO=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", key);
            bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
            free(tmp.s);
            int res = hrec ? bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) : -1;
            if (res < 0) bcf_hrec_destroy(hrec);
            if (res > 0) res = bcf_hdr_sync((bcf_hdr_t*)h);
            k = kh_get(vdict, d, key);
            v->errcode = BCF_ERR_TAG_UNDEF;
            if (res || k == kh_end(d)) {
                hts_log_error("Could not add dummy header for INFO '%s' at %s:%"PRIhts_pos, key, bcf_seqname_safe(h,v), v->pos+1);
                v->errcode |= BCF_ERR_TAG_INVALID;
                return -1;
            }
        }
        uint32_t y = kh_val(d, k).info[BCF_HL_INFO];
        ++v->n_info;
        bcf_enc_int1(str, kh_val(d, k).id);
        if (val == 0) {
            bcf_enc_size(str, 0, BCF_BT_NULL);
        } else if ((y>>4&0xf) == BCF_HT_FLAG || (y>>4&0xf) == BCF_HT_STR) { // if Flag has a value, treat it as a string
            bcf_enc_vchar(str, end - val, val);
        } else { // int/float value/array
            int i, n_val;
            char *t, *te;
            for (t = val, n_val = 1; *t; ++t) // count the number of values
                if (*t == ',') ++n_val;
            // Check both int and float size in one step for simplicity
            if (n_val > max_n_val) {
                int32_t *a_tmp = (int32_t *)realloc(a_val, n_val * sizeof(*a_val));
                if (!a_tmp) {
                    hts_log_error("Could not allocate memory at %s:%"PRIhts_pos, bcf_seqname_safe(h,v), v->pos+1);
                    v->errcode |= BCF_ERR_LIMITS; // No appropriate code?
                    return -1;
                }
                a_val = a_tmp;
                max_n_val = n_val;
            }
            if ((y>>4&0xf) == BCF_HT_INT) {
                i = 0, t = val;
                int64_t val1;
                int is_int64 = 0;
#ifdef VCF_ALLOW_INT64
                if ( n_val==1 )
                {
                    overflow = 0;
                    long long int tmp_val = hts_str2int(val, &te, sizeof(tmp_val)*CHAR_BIT, &overflow);
                    if ( te==val ) tmp_val = bcf_int32_missing;
                    else if ( overflow || tmp_val<BCF_MIN_BT_INT64 || tmp_val>BCF_MAX_BT_INT64 )
                    {
                        if ( !extreme_int_warned )
                        {
                            hts_log_warning("Extreme INFO/%s value encountered and set to missing at %s:%"PRIhts_pos,key,bcf_seqname_safe(h,v), v->pos+1);
                            extreme_int_warned = 1;
                        }
                        tmp_val = bcf_int32_missing;
                    }
                    else
                        is_int64 = 1;
                    val1 = tmp_val;
                    t = te;
                    i = 1;  // this is just to avoid adding another nested block...
                }
#endif
                for (; i < n_val; ++i, ++t)
                {
                    overflow = 0;
                    long int tmp_val = hts_str2int(t, &te, sizeof(tmp_val)*CHAR_BIT, &overflow);
                    if ( te==t ) tmp_val = bcf_int32_missing;
                    else if ( overflow || tmp_val<BCF_MIN_BT_INT32 || tmp_val>BCF_MAX_BT_INT32 )
                    {
                        if ( !extreme_int_warned )
                        {
                            hts_log_warning("Extreme INFO/%s value encountered and set to missing at %s:%"PRIhts_pos,key,bcf_seqname_safe(h,v), v->pos+1);
                            extreme_int_warned = 1;
                        }
                        tmp_val = bcf_int32_missing;
                    }
                    a_val[i] = tmp_val;
                    for (t = te; *t && *t != ','; t++);
                }
                if (n_val == 1) {
#ifdef VCF_ALLOW_INT64
                    if ( is_int64 )
                    {
                        v->unpacked |= BCF_IS_64BIT;
                        bcf_enc_long1(str, val1);
                    }
                    else
                        bcf_enc_int1(str, (int32_t)val1);
#else
                    val1 = a_val[0];
                    bcf_enc_int1(str, (int32_t)val1);
#endif
                } else {
                    bcf_enc_vint(str, n_val, a_val, -1);
                }
                if (n_val==1 && (val1!=bcf_int32_missing || is_int64) && strcmp(key, "END") == 0)
                {
                    if ( val1 <= v->pos )
                    {
                        if ( !negative_rlen_warned )
                        {
                            hts_log_warning("INFO/END=%"PRIhts_pos" is smaller than POS at %s:%"PRIhts_pos,val1,bcf_seqname_safe(h,v),v->pos+1);
                            negative_rlen_warned = 1;
                        }
                    }
                    else
                        v->rlen = val1 - v->pos;
                }
            } else if ((y>>4&0xf) == BCF_HT_REAL) {
                float *val_f = (float *)a_val;
                for (i = 0, t = val; i < n_val; ++i, ++t)
                {
                    overflow = 0;
                    val_f[i] = hts_str2dbl(t, &te, &overflow);
                    if ( te==t || overflow ) // conversion failed
                        bcf_float_set_missing(val_f[i]);
                    for (t = te; *t && *t != ','; t++);
                }
                bcf_enc_vfloat(str, n_val, val_f);
            }
        }
        if (c == 0) break;
        r = end;
        key = r + 1;
    }

    free(a_val);
    return 0;
}

int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
{
    int i = 0, ret = -2, overflow = 0;
    char *p, *q, *r, *t;
    kstring_t *str;
    khint_t k;
    ks_tokaux_t aux;

    if (!s || !h || !v || !(s->s))
        return ret;

    // Assumed in lots of places, but we may as well spot this early
    assert(sizeof(float) == sizeof(int32_t));

    bcf_clear1(v);
    str = &v->shared;
    memset(&aux, 0, sizeof(ks_tokaux_t));
    for (p = kstrtok(s->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
        q = (char*)aux.p;
        *q = 0;
        if (i == 0) { // CHROM
            vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
            k = kh_get(vdict, d, p);
            if (k == kh_end(d))
            {
                hts_log_warning("Contig '%s' is not defined in the header. (Quick workaround: index the file with tabix.)", p);
                v->errcode = BCF_ERR_CTG_UNDEF;
                if ((k = fix_chromosome(h, d, p)) == kh_end(d)) {
                    hts_log_error("Could not add dummy header for contig '%s'", p);
                    v->errcode |= BCF_ERR_CTG_INVALID;
                    goto err;
                }
            }
            v->rid = kh_val(d, k).id;
        } else if (i == 1) { // POS
            overflow = 0;
            v->pos = hts_str2uint(p, &p, 63, &overflow);
            if (overflow) {
                hts_log_error("Position value '%s' is too large", p);
                goto err;
            } else {
                v->pos -= 1;
            }
            if (v->pos >= INT32_MAX)
                v->unpacked |= BCF_IS_64BIT;
        } else if (i == 2) { // ID
            if (strcmp(p, ".")) bcf_enc_vchar(str, q - p, p);
            else bcf_enc_size(str, 0, BCF_BT_CHAR);
        } else if (i == 3) { // REF
            bcf_enc_vchar(str, q - p, p);
            v->n_allele = 1, v->rlen = q - p;
        } else if (i == 4) { // ALT
            if (strcmp(p, ".")) {
                for (r = t = p;; ++r) {
                    if (*r == ',' || *r == 0) {
                        if (v->n_allele == UINT16_MAX) {
                            hts_log_error("Too many ALT alleles at %s:%"PRIhts_pos,
                                          bcf_seqname_safe(h,v), v->pos+1);
                            v->errcode |= BCF_ERR_LIMITS;
                            goto err;
                        }
                        bcf_enc_vchar(str, r - t, t);
                        t = r + 1;
                        ++v->n_allele;
                    }
                    if (r == q) break;
                }
            }
        } else if (i == 5) { // QUAL
            if (strcmp(p, ".")) v->qual = atof(p);
            else bcf_float_set_missing(v->qual);
            if ( v->max_unpack && !(v->max_unpack>>1) ) goto end; // BCF_UN_STR
        } else if (i == 6) { // FILTER
            if (strcmp(p, ".")) {
                if (vcf_parse_filter(str, h, v, p, q)) goto err;
            } else bcf_enc_vint(str, 0, 0, -1);
            if ( v->max_unpack && !(v->max_unpack>>2) ) goto end; // BCF_UN_FLT
        } else if (i == 7) { // INFO
            if (strcmp(p, ".")) {
                if (vcf_parse_info(str, h, v, p, q)) goto err;
            }
            if ( v->max_unpack && !(v->max_unpack>>3) ) goto end;
        } else if (i == 8) {// FORMAT
            return vcf_parse_format(s, h, v, p, q) == 0 ? 0 : -2;
        }
    }

 end:
    ret = 0;

 err:
    return ret;
}

int vcf_open_mode(char *mode, const char *fn, const char *format)
{
    if (format == NULL) {
        // Try to pick a format based on the filename extension
        char extension[HTS_MAX_EXT_LEN];
        if (find_file_extension(fn, extension) < 0) return -1;
        return vcf_open_mode(mode, fn, extension);
    }
    else if (strcasecmp(format, "bcf") == 0) strcpy(mode, "b");
    else if (strcasecmp(format, "vcf") == 0) strcpy(mode, "");
    else if (strcasecmp(format, "vcf.gz") == 0 || strcasecmp(format, "vcf.bgz") == 0) strcpy(mode, "z");
    else return -1;

    return 0;
}

int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
    int ret;
    ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
    if (ret < 0) return ret;
    return vcf_parse1(&fp->line, h, v);
}

static inline uint8_t *bcf_unpack_fmt_core1(uint8_t *ptr, int n_sample, bcf_fmt_t *fmt)
{
    uint8_t *ptr_start = ptr;
    fmt->id = bcf_dec_typed_int1(ptr, &ptr);
    fmt->n = bcf_dec_size(ptr, &ptr, &fmt->type);
    fmt->size = fmt->n << bcf_type_shift[fmt->type];
    fmt->p = ptr;
    fmt->p_off  = ptr - ptr_start;
    fmt->p_free = 0;
    ptr += n_sample * fmt->size;
    fmt->p_len = ptr - fmt->p;
    return ptr;
}

static inline uint8_t *bcf_unpack_info_core1(uint8_t *ptr, bcf_info_t *info)
{
    uint8_t *ptr_start = ptr;
    info->key = bcf_dec_typed_int1(ptr, &ptr);
    info->len = bcf_dec_size(ptr, &ptr, &info->type);
    info->vptr = ptr;
    info->vptr_off  = ptr - ptr_start;
    info->vptr_free = 0;
    info->v1.i = 0;
    if (info->len == 1) {
        if (info->type == BCF_BT_INT8 || info->type == BCF_BT_CHAR) info->v1.i = *(int8_t*)ptr;
        else if (info->type == BCF_BT_INT32) info->v1.i = le_to_i32(ptr);
        else if (info->type == BCF_BT_FLOAT) info->v1.f = le_to_float(ptr);
        else if (info->type == BCF_BT_INT16) info->v1.i = le_to_i16(ptr);
        else if (info->type == BCF_BT_INT64) info->v1.i = le_to_i64(ptr);
    }
    ptr += info->len << bcf_type_shift[info->type];
    info->vptr_len = ptr - info->vptr;
    return ptr;
}

int bcf_unpack(bcf1_t *b, int which)
{
    if ( !b->shared.l ) return 0; // Building a new BCF record from scratch
    uint8_t *ptr = (uint8_t*)b->shared.s, *ptr_ori;
    int i;
    bcf_dec_t *d = &b->d;
    if (which & BCF_UN_FLT) which |= BCF_UN_STR;
    if (which & BCF_UN_INFO) which |= BCF_UN_SHR;
    if ((which&BCF_UN_STR) && !(b->unpacked&BCF_UN_STR))
    {
        kstring_t tmp;

        // ID
        tmp.l = 0; tmp.s = d->id; tmp.m = d->m_id;
        ptr_ori = ptr;
        ptr = bcf_fmt_sized_array(&tmp, ptr);
        b->unpack_size[0] = ptr - ptr_ori;
        kputc('\0', &tmp);
        d->id = tmp.s; d->m_id = tmp.m;

        // REF and ALT are in a single block (d->als) and d->alleles are pointers into this block
        hts_expand(char*, b->n_allele, d->m_allele, d->allele); // NM: hts_expand() is a macro
        tmp.l = 0; tmp.s = d->als; tmp.m = d->m_als;
        ptr_ori = ptr;
        for (i = 0; i < b->n_allele; ++i) {
            // Use offset within tmp.s as realloc may change pointer
            d->allele[i] = (char *)(intptr_t)tmp.l;
            ptr = bcf_fmt_sized_array(&tmp, ptr);
            kputc('\0', &tmp);
        }
        b->unpack_size[1] = ptr - ptr_ori;
        d->als = tmp.s; d->m_als = tmp.m;

        // Convert our offsets within tmp.s back to pointers again
        for (i = 0; i < b->n_allele; ++i)
            d->allele[i] = d->als + (ptrdiff_t)d->allele[i];
        b->unpacked |= BCF_UN_STR;
    }
    if ((which&BCF_UN_FLT) && !(b->unpacked&BCF_UN_FLT)) { // FILTER
        ptr = (uint8_t*)b->shared.s + b->unpack_size[0] + b->unpack_size[1];
        ptr_ori = ptr;
        if (*ptr>>4) {
            int type;
            d->n_flt = bcf_dec_size(ptr, &ptr, &type);
            hts_expand(int, d->n_flt, d->m_flt, d->flt);
            for (i = 0; i < d->n_flt; ++i)
                d->flt[i] = bcf_dec_int1(ptr, type, &ptr);
        } else ++ptr, d->n_flt = 0;
        b->unpack_size[2] = ptr - ptr_ori;
        b->unpacked |= BCF_UN_FLT;
    }
    if ((which&BCF_UN_INFO) && !(b->unpacked&BCF_UN_INFO)) { // INFO
        ptr = (uint8_t*)b->shared.s + b->unpack_size[0] + b->unpack_size[1] + b->unpack_size[2];
        hts_expand(bcf_info_t, b->n_info, d->m_info, d->info);
        for (i = 0; i < d->m_info; ++i) d->info[i].vptr_free = 0;
        for (i = 0; i < b->n_info; ++i)
            ptr = bcf_unpack_info_core1(ptr, &d->info[i]);
        b->unpacked |= BCF_UN_INFO;
    }
    if ((which&BCF_UN_FMT) && b->n_sample && !(b->unpacked&BCF_UN_FMT)) { // FORMAT
        ptr = (uint8_t*)b->indiv.s;
        hts_expand(bcf_fmt_t, b->n_fmt, d->m_fmt, d->fmt);
        for (i = 0; i < d->m_fmt; ++i) d->fmt[i].p_free = 0;
        for (i = 0; i < b->n_fmt; ++i)
            ptr = bcf_unpack_fmt_core1(ptr, b->n_sample, &d->fmt[i]);
        b->unpacked |= BCF_UN_FMT;
    }
    return 0;
}

int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)
{
    int i;
    int32_t max_dt_id = h->n[BCF_DT_ID];
    const char *chrom = bcf_seqname(h, v);
    if (!chrom) {
        hts_log_error("Invalid BCF, CONTIG id=%d not present in the header",
                      v->rid);
        errno = EINVAL;
        return -1;
    }
    bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
    kputs(chrom, s); // CHROM
    kputc('\t', s); kputll(v->pos + 1, s); // POS
    kputc('\t', s); kputs(v->d.id ? v->d.id : ".", s); // ID
    kputc('\t', s); // REF
    if (v->n_allele > 0) kputs(v->d.allele[0], s);
    else kputc('.', s);
    kputc('\t', s); // ALT
    if (v->n_allele > 1) {
        for (i = 1; i < v->n_allele; ++i) {
            if (i > 1) kputc(',', s);
            kputs(v->d.allele[i], s);
        }
    } else kputc('.', s);
    kputc('\t', s); // QUAL
    if ( bcf_float_is_missing(v->qual) ) kputc('.', s); // QUAL
    else kputd(v->qual, s);
    kputc('\t', s); // FILTER
    if (v->d.n_flt) {
        for (i = 0; i < v->d.n_flt; ++i) {
            int32_t idx = v->d.flt[i];
            if (idx < 0 || idx >= max_dt_id
                || h->id[BCF_DT_ID][idx].key == NULL) {
                hts_log_error("Invalid BCF, the FILTER tag id=%d at %s:%"PRIhts_pos" not present in the header",
                              idx, bcf_seqname_safe(h, v), v->pos + 1);
                errno = EINVAL;
                return -1;
            }
            if (i) kputc(';', s);
            kputs(h->id[BCF_DT_ID][idx].key, s);
        }
    } else kputc('.', s);
    kputc('\t', s); // INFO
    if (v->n_info) {
        int first = 1;
        for (i = 0; i < v->n_info; ++i) {
            bcf_info_t *z = &v->d.info[i];
            if ( !z->vptr ) continue;
            if ( !first ) kputc(';', s);
            first = 0;
            if (z->key < 0 || z->key >= max_dt_id
                || h->id[BCF_DT_ID][z->key].key == NULL) {
                hts_log_error("Invalid BCF, the INFO tag id=%d is %s at %s:%"PRIhts_pos,
                              z->key,
                              z->key < 0 ? "negative"
                              : (z->key >= max_dt_id ? "too large" : "not present in the header"),
                              bcf_seqname_safe(h, v), v->pos+1);
                errno = EINVAL;
                return -1;
            }
            kputs(h->id[BCF_DT_ID][z->key].key, s);
            if (z->len <= 0) continue;
            kputc('=', s);
            if (z->len == 1)
            {
                switch (z->type)
                {
                    case BCF_BT_INT8:  if ( z->v1.i==bcf_int8_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_INT16: if ( z->v1.i==bcf_int16_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_INT32: if ( z->v1.i==bcf_int32_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_INT64: if ( z->v1.i==bcf_int64_missing ) kputc('.', s); else kputll(z->v1.i, s); break;
                    case BCF_BT_FLOAT: if ( bcf_float_is_missing(z->v1.f) ) kputc('.', s); else kputd(z->v1.f, s); break;
                    case BCF_BT_CHAR:  kputc(z->v1.i, s); break;
                    default:
                        hts_log_error("Unexpected type %d at %s:%"PRIhts_pos, z->type, bcf_seqname_safe(h, v), v->pos+1);
                        errno = EINVAL;
                        return -1;
                }
            }
            else bcf_fmt_array(s, z->len, z->type, z->vptr);
        }
        if ( first ) kputc('.', s);
    } else kputc('.', s);
    // FORMAT and individual information
    if (v->n_sample)
    {
        int i,j;
        if ( v->n_fmt)
        {
            int gt_i = -1;
            bcf_fmt_t *fmt = v->d.fmt;
            int first = 1;
            for (i = 0; i < (int)v->n_fmt; ++i) {
                if ( !fmt[i].p ) continue;
                kputc(!first ? ':' : '\t', s); first = 0;
                if (fmt[i].id < 0 || fmt[i].id >= max_dt_id
                    || h->id[BCF_DT_ID][fmt[i].id].key == NULL) //!bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) )
                {
                    hts_log_error("Invalid BCF, the FORMAT tag id=%d at %s:%"PRIhts_pos" not present in the header", fmt[i].id, bcf_seqname_safe(h, v), v->pos+1);
                    errno = EINVAL;
                    return -1;
                }
                kputs(h->id[BCF_DT_ID][fmt[i].id].key, s);
                if (strcmp(h->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0) gt_i = i;
            }
            if ( first ) kputs("\t.", s);
            for (j = 0; j < v->n_sample; ++j) {
                kputc('\t', s);
                first = 1;
                for (i = 0; i < (int)v->n_fmt; ++i) {
                    bcf_fmt_t *f = &fmt[i];
                    if ( !f->p ) continue;
                    if (!first) kputc(':', s);
                    first = 0;
                    if (gt_i == i)
                        bcf_format_gt(f,j,s);
                    else
                        bcf_fmt_array(s, f->n, f->type, f->p + j * (size_t)f->size);
                }
                if ( first ) kputc('.', s);
            }
        }
        else
            for (j=0; j<=v->n_sample; j++)
                kputs("\t.", s);
    }
    kputc('\n', s);
    return 0;
}

int vcf_write_line(htsFile *fp, kstring_t *line)
{
    int ret;
    if ( line->s[line->l-1]!='\n' ) kputc('\n',line);
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, line->s, line->l);
    else
        ret = hwrite(fp->fp.hfile, line->s, line->l);
    return ret==line->l ? 0 : -1;
}

int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
    int ret;
    fp->line.l = 0;
    if (vcf_format1(h, v, &fp->line) != 0)
        return -1;
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l);
    else
        ret = hwrite(fp->fp.hfile, fp->line.s, fp->line.l);

    if (fp->idx) {
        int tid;
        if ((tid = hts_idx_tbi_name(fp->idx, v->rid, bcf_seqname_safe(h, v))) < 0)
            return -1;

        if (hts_idx_push(fp->idx, tid, v->pos, v->pos + v->rlen, bgzf_tell(fp->fp.bgzf), 1) < 0)
            return -1;
    }

    return ret==fp->line.l ? 0 : -1;
}

/************************
 * Data access routines *
 ************************/

int bcf_hdr_id2int(const bcf_hdr_t *h, int which, const char *id)
{
    khint_t k;
    vdict_t *d = (vdict_t*)h->dict[which];
    k = kh_get(vdict, d, id);
    return k == kh_end(d)? -1 : kh_val(d, k).id;
}


/********************
 *** BCF indexing ***
 ********************/

// Calculate number of index levels given min_shift and the header contig
// list.  Also returns number of contigs in *nids_out.
static int idx_calc_n_lvls_ids(const bcf_hdr_t *h, int min_shift,
                               int starting_n_lvls, int *nids_out)
{
    int n_lvls, i, nids = 0;
    int64_t max_len = 0, s;

    for (i = 0; i < h->n[BCF_DT_CTG]; ++i)
    {
        if ( !h->id[BCF_DT_CTG][i].val ) continue;
        if ( max_len < h->id[BCF_DT_CTG][i].val->info[0] )
            max_len = h->id[BCF_DT_CTG][i].val->info[0];
        nids++;
    }
    if ( !max_len ) max_len = (1LL<<31) - 1;  // In case contig line is broken.
    max_len += 256;
    s = 1LL << (min_shift + starting_n_lvls * 3);
    for (n_lvls = starting_n_lvls; max_len > s; ++n_lvls, s <<= 3);

    if (nids_out) *nids_out = nids;
    return n_lvls;
}

hts_idx_t *bcf_index(htsFile *fp, int min_shift)
{
    int n_lvls;
    bcf1_t *b = NULL;
    hts_idx_t *idx = NULL;
    bcf_hdr_t *h;
    int r;
    h = bcf_hdr_read(fp);
    if ( !h ) return NULL;
    int nids = 0;
    n_lvls = idx_calc_n_lvls_ids(h, min_shift, 0, &nids);
    idx = hts_idx_init(nids, HTS_FMT_CSI, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
    if (!idx) goto fail;
    b = bcf_init1();
    if (!b) goto fail;
    while ((r = bcf_read1(fp,h, b)) >= 0) {
        int ret;
        ret = hts_idx_push(idx, b->rid, b->pos, b->pos + b->rlen, bgzf_tell(fp->fp.bgzf), 1);
        if (ret < 0) goto fail;
    }
    if (r < -1) goto fail;
    hts_idx_finish(idx, bgzf_tell(fp->fp.bgzf));
    bcf_destroy1(b);
    bcf_hdr_destroy(h);
    return idx;

 fail:
    hts_idx_destroy(idx);
    bcf_destroy1(b);
    bcf_hdr_destroy(h);
    return NULL;
}

hts_idx_t *bcf_index_load2(const char *fn, const char *fnidx)
{
    return fnidx? hts_idx_load2(fn, fnidx) : bcf_index_load(fn);
}

hts_idx_t *bcf_index_load3(const char *fn, const char *fnidx, int flags)
{
    return hts_idx_load3(fn, fnidx, HTS_FMT_CSI, flags);
}

int bcf_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads)
{
    htsFile *fp;
    hts_idx_t *idx;
    tbx_t *tbx;
    int ret;
    if ((fp = hts_open(fn, "rb")) == 0) return -2;
    if (n_threads)
        hts_set_threads(fp, n_threads);
    if ( fp->format.compression!=bgzf ) { hts_close(fp); return -3; }
    switch (fp->format.format) {
        case bcf:
            if (!min_shift) {
                hts_log_error("TBI indices for BCF files are not supported");
                ret = -1;
            } else {
                idx = bcf_index(fp, min_shift);
                if (idx) {
                    ret = hts_idx_save_as(idx, fn, fnidx, HTS_FMT_CSI);
                    if (ret < 0) ret = -4;
                    hts_idx_destroy(idx);
                }
                else ret = -1;
            }
            break;

        case vcf:
            tbx = tbx_index(hts_get_bgzfp(fp), min_shift, &tbx_conf_vcf);
            if (tbx) {
                ret = hts_idx_save_as(tbx->idx, fn, fnidx, min_shift > 0 ? HTS_FMT_CSI : HTS_FMT_TBI);
                if (ret < 0) ret = -4;
                tbx_destroy(tbx);
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

int bcf_index_build2(const char *fn, const char *fnidx, int min_shift)
{
    return bcf_index_build3(fn, fnidx, min_shift, 0);
}

int bcf_index_build(const char *fn, int min_shift)
{
    return bcf_index_build3(fn, NULL, min_shift, 0);
}

// Initialise fp->idx for the current format type.
// This must be called after the header has been written but no other data.
static int vcf_idx_init(htsFile *fp, bcf_hdr_t *h, int min_shift, const char *fnidx) {
    int n_lvls, fmt;

    if (min_shift == 0) {
        min_shift = 14;
        n_lvls = 5;
        fmt = HTS_FMT_TBI;
    } else {
        // Set initial n_lvls to match tbx_index()
        int starting_n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3;
        // Increase if necessary
        n_lvls = idx_calc_n_lvls_ids(h, min_shift, starting_n_lvls, NULL);
        fmt = HTS_FMT_CSI;
    }

    fp->idx = hts_idx_init(0, fmt, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
    if (!fp->idx) return -1;

    // Tabix meta data, added even in CSI for VCF
    uint8_t conf[4*7];
    u32_to_le(TBX_VCF, conf+0);  // fmt
    u32_to_le(1,       conf+4);  // name col
    u32_to_le(2,       conf+8);  // beg col
    u32_to_le(0,       conf+12); // end col
    u32_to_le('#',     conf+16); // comment
    u32_to_le(0,       conf+20); // n.skip
    u32_to_le(0,       conf+24); // ref name len
    if (hts_idx_set_meta(fp->idx, sizeof(conf)*sizeof(*conf), (uint8_t *)conf, 1) < 0) {
        hts_idx_destroy(fp->idx);
        fp->idx = NULL;
        return -1;
    }
    fp->fnidx = fnidx;

    return 0;
}

// Initialise fp->idx for the current format type.
// This must be called after the header has been written but no other data.
int bcf_idx_init(htsFile *fp, bcf_hdr_t *h, int min_shift, const char *fnidx) {
    int n_lvls, nids = 0;

    if (fp->format.format == vcf)
        return vcf_idx_init(fp, h, min_shift, fnidx);

    if (!min_shift)
        min_shift = 14;

    n_lvls = idx_calc_n_lvls_ids(h, min_shift, 0, &nids);

    fp->idx = hts_idx_init(nids, HTS_FMT_CSI, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
    if (!fp->idx) return -1;
    fp->fnidx = fnidx;

    return 0;
}

// Finishes an index. Call after the last record has been written.
// Returns 0 on success, <0 on failure.
//
// NB: same format as SAM/BAM as it uses bgzf.
int bcf_idx_save(htsFile *fp) {
    return sam_idx_save(fp);
}

/*****************
 *** Utilities ***
 *****************/

int bcf_hdr_combine(bcf_hdr_t *dst, const bcf_hdr_t *src)
{
    int i, ndst_ori = dst->nhrec, need_sync = 0, ret = 0, res;
    for (i=0; i<src->nhrec; i++)
    {
        if ( src->hrec[i]->type==BCF_HL_GEN && src->hrec[i]->value )
        {
            int j;
            for (j=0; j<ndst_ori; j++)
            {
                if ( dst->hrec[j]->type!=BCF_HL_GEN ) continue;

                // Checking only the key part of generic lines, otherwise
                // the VCFs are too verbose. Should we perhaps add a flag
                // to bcf_hdr_combine() and make this optional?
                if ( !strcmp(src->hrec[i]->key,dst->hrec[j]->key) ) break;
            }
            if ( j>=ndst_ori ) {
                res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                if (res < 0) return -1;
                need_sync += res;
            }
        }
        else if ( src->hrec[i]->type==BCF_HL_STR )
        {
            // NB: we are ignoring fields without ID
            int j = bcf_hrec_find_key(src->hrec[i],"ID");
            if ( j>=0 )
            {
                bcf_hrec_t *rec = bcf_hdr_get_hrec(dst, src->hrec[i]->type, "ID", src->hrec[i]->vals[j], src->hrec[i]->key);
                if ( !rec ) {
                    res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                    if (res < 0) return -1;
                    need_sync += res;
                }
            }
        }
        else
        {
            int j = bcf_hrec_find_key(src->hrec[i],"ID");
            assert( j>=0 ); // this should always be true for valid VCFs

            bcf_hrec_t *rec = bcf_hdr_get_hrec(dst, src->hrec[i]->type, "ID", src->hrec[i]->vals[j], NULL);
            if ( !rec ) {
                res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                if (res < 0) return -1;
                need_sync += res;
            } else if ( src->hrec[i]->type==BCF_HL_INFO || src->hrec[i]->type==BCF_HL_FMT )
            {
                // Check that both records are of the same type. The bcf_hdr_id2length
                // macro cannot be used here because dst header is not synced yet.
                vdict_t *d_src = (vdict_t*)src->dict[BCF_DT_ID];
                vdict_t *d_dst = (vdict_t*)dst->dict[BCF_DT_ID];
                khint_t k_src  = kh_get(vdict, d_src, src->hrec[i]->vals[0]);
                khint_t k_dst  = kh_get(vdict, d_dst, src->hrec[i]->vals[0]);
                if ( (kh_val(d_src,k_src).info[rec->type]>>8 & 0xf) != (kh_val(d_dst,k_dst).info[rec->type]>>8 & 0xf) )
                {
                    hts_log_warning("Trying to combine \"%s\" tag definitions of different lengths",
                        src->hrec[i]->vals[0]);
                    ret |= 1;
                }
                if ( (kh_val(d_src,k_src).info[rec->type]>>4 & 0xf) != (kh_val(d_dst,k_dst).info[rec->type]>>4 & 0xf) )
                {
                    hts_log_warning("Trying to combine \"%s\" tag definitions of different types",
                        src->hrec[i]->vals[0]);
                    ret |= 1;
                }
            }
        }
    }
    if ( need_sync ) {
        if (bcf_hdr_sync(dst) < 0) return -1;
    }
    return ret;
}

bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const bcf_hdr_t *src)
{
    if ( !dst )
    {
        // this will effectively strip existing IDX attributes from src to become dst
        dst = bcf_hdr_init("r");
        kstring_t htxt = {0,0,0};
        if (bcf_hdr_format(src, 0, &htxt) < 0) {
            free(htxt.s);
            return NULL;
        }
        if ( bcf_hdr_parse(dst, htxt.s) < 0 ) {
            bcf_hdr_destroy(dst);
            dst = NULL;
        }
        free(htxt.s);
        return dst;
    }

    int i, ndst_ori = dst->nhrec, need_sync = 0, ret = 0, res;
    for (i=0; i<src->nhrec; i++)
    {
        if ( src->hrec[i]->type==BCF_HL_GEN && src->hrec[i]->value )
        {
            int j;
            for (j=0; j<ndst_ori; j++)
            {
                if ( dst->hrec[j]->type!=BCF_HL_GEN ) continue;

                // Checking only the key part of generic lines, otherwise
                // the VCFs are too verbose. Should we perhaps add a flag
                // to bcf_hdr_combine() and make this optional?
                if ( !strcmp(src->hrec[i]->key,dst->hrec[j]->key) ) break;
            }
            if ( j>=ndst_ori ) {
                res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                if (res < 0) return NULL;
                need_sync += res;
            }
        }
        else if ( src->hrec[i]->type==BCF_HL_STR )
        {
            // NB: we are ignoring fields without ID
            int j = bcf_hrec_find_key(src->hrec[i],"ID");
            if ( j>=0 )
            {
                bcf_hrec_t *rec = bcf_hdr_get_hrec(dst, src->hrec[i]->type, "ID", src->hrec[i]->vals[j], src->hrec[i]->key);
                if ( !rec ) {
                    res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                    if (res < 0) return NULL;
                    need_sync += res;
                }
            }
        }
        else
        {
            int j = bcf_hrec_find_key(src->hrec[i],"ID");
            assert( j>=0 ); // this should always be true for valid VCFs

            bcf_hrec_t *rec = bcf_hdr_get_hrec(dst, src->hrec[i]->type, "ID", src->hrec[i]->vals[j], NULL);
            if ( !rec ) {
                res = bcf_hdr_add_hrec(dst, bcf_hrec_dup(src->hrec[i]));
                if (res < 0) return NULL;
                need_sync += res;
            } else if ( src->hrec[i]->type==BCF_HL_INFO || src->hrec[i]->type==BCF_HL_FMT )
            {
                // Check that both records are of the same type. The bcf_hdr_id2length
                // macro cannot be used here because dst header is not synced yet.
                vdict_t *d_src = (vdict_t*)src->dict[BCF_DT_ID];
                vdict_t *d_dst = (vdict_t*)dst->dict[BCF_DT_ID];
                khint_t k_src  = kh_get(vdict, d_src, src->hrec[i]->vals[0]);
                khint_t k_dst  = kh_get(vdict, d_dst, src->hrec[i]->vals[0]);
                if ( (kh_val(d_src,k_src).info[rec->type]>>8 & 0xf) != (kh_val(d_dst,k_dst).info[rec->type]>>8 & 0xf) )
                {
                    hts_log_warning("Trying to combine \"%s\" tag definitions of different lengths",
                        src->hrec[i]->vals[0]);
                    ret |= 1;
                }
                if ( (kh_val(d_src,k_src).info[rec->type]>>4 & 0xf) != (kh_val(d_dst,k_dst).info[rec->type]>>4 & 0xf) )
                {
                    hts_log_warning("Trying to combine \"%s\" tag definitions of different types",
                        src->hrec[i]->vals[0]);
                    ret |= 1;
                }
            }
        }
    }
    if ( need_sync ) {
        if (bcf_hdr_sync(dst) < 0) return NULL;
    }
    return dst;
}

int bcf_translate(const bcf_hdr_t *dst_hdr, bcf_hdr_t *src_hdr, bcf1_t *line)
{
    int i;
    if ( line->errcode )
    {
        hts_log_error("Unchecked error (%d) at %s:%"PRIhts_pos", exiting", line->errcode, bcf_seqname_safe(src_hdr,line), line->pos+1);
        exit(1);
    }
    if ( src_hdr->ntransl==-1 ) return 0;    // no need to translate, all tags have the same id
    if ( !src_hdr->ntransl )  // called for the first time, see what needs translating
    {
        int dict;
        for (dict=0; dict<2; dict++)    // BCF_DT_ID and BCF_DT_CTG
        {
            src_hdr->transl[dict] = (int*) malloc(src_hdr->n[dict]*sizeof(int));
            for (i=0; i<src_hdr->n[dict]; i++)
            {
                if ( !src_hdr->id[dict][i].key ) // gap left after removed BCF header lines
                {
                    src_hdr->transl[dict][i] = -1;
                    continue;
                }
                src_hdr->transl[dict][i] = bcf_hdr_id2int(dst_hdr,dict,src_hdr->id[dict][i].key);
                if ( src_hdr->transl[dict][i]!=-1 && i!=src_hdr->transl[dict][i] ) src_hdr->ntransl++;
            }
        }
        if ( !src_hdr->ntransl )
        {
            free(src_hdr->transl[0]); src_hdr->transl[0] = NULL;
            free(src_hdr->transl[1]); src_hdr->transl[1] = NULL;
            src_hdr->ntransl = -1;
        }
        if ( src_hdr->ntransl==-1 ) return 0;
    }
    bcf_unpack(line,BCF_UN_ALL);

    // CHROM
    if ( src_hdr->transl[BCF_DT_CTG][line->rid] >=0 ) line->rid = src_hdr->transl[BCF_DT_CTG][line->rid];

    // FILTER
    for (i=0; i<line->d.n_flt; i++)
    {
        int src_id = line->d.flt[i];
        if ( src_hdr->transl[BCF_DT_ID][src_id] >=0 )
            line->d.flt[i] = src_hdr->transl[BCF_DT_ID][src_id];
        line->d.shared_dirty |= BCF1_DIRTY_FLT;
    }

    // INFO
    for (i=0; i<line->n_info; i++)
    {
        int src_id = line->d.info[i].key;
        int dst_id = src_hdr->transl[BCF_DT_ID][src_id];
        if ( dst_id<0 ) continue;
        line->d.info[i].key = dst_id;
        if ( !line->d.info[i].vptr ) continue;  // skip deleted
        int src_size = src_id>>7 ? ( src_id>>15 ? BCF_BT_INT32 : BCF_BT_INT16) : BCF_BT_INT8;
        int dst_size = dst_id>>7 ? ( dst_id>>15 ? BCF_BT_INT32 : BCF_BT_INT16) : BCF_BT_INT8;
        if ( src_size==dst_size )   // can overwrite
        {
            uint8_t *vptr = line->d.info[i].vptr - line->d.info[i].vptr_off;
            if ( dst_size==BCF_BT_INT8 ) { vptr[1] = (uint8_t)dst_id; }
            else if ( dst_size==BCF_BT_INT16 ) { *(uint16_t*)vptr = (uint16_t)dst_id; }
            else { *(uint32_t*)vptr = (uint32_t)dst_id; }
        }
        else    // must realloc
        {
            bcf_info_t *info = &line->d.info[i];
            kstring_t str = {0,0,0};
            bcf_enc_int1(&str, dst_id);
            bcf_enc_size(&str, info->len,info->type);
            uint32_t vptr_off = str.l;
            kputsn((char*)info->vptr, info->vptr_len, &str);
            if( info->vptr_free ) free(info->vptr - info->vptr_off);
            info->vptr_off = vptr_off;
            info->vptr = (uint8_t*)str.s + info->vptr_off;
            info->vptr_free = 1;
            line->d.shared_dirty |= BCF1_DIRTY_INF;
        }
    }

    // FORMAT
    for (i=0; i<line->n_fmt; i++)
    {
        int src_id = line->d.fmt[i].id;
        int dst_id = src_hdr->transl[BCF_DT_ID][src_id];
        if ( dst_id<0 ) continue;
        line->d.fmt[i].id = dst_id;
        if( !line->d.fmt[i].p ) continue;  // skip deleted
        int src_size = src_id>>7 ? ( src_id>>15 ? BCF_BT_INT32 : BCF_BT_INT16) : BCF_BT_INT8;
        int dst_size = dst_id>>7 ? ( dst_id>>15 ? BCF_BT_INT32 : BCF_BT_INT16) : BCF_BT_INT8;
        if ( src_size==dst_size )   // can overwrite
        {
            uint8_t *p = line->d.fmt[i].p - line->d.fmt[i].p_off;    // pointer to the vector size (4bits) and BT type (4bits)
            if ( dst_size==BCF_BT_INT8 ) { p[1] = dst_id; }
            else if ( dst_size==BCF_BT_INT16 ) { i16_to_le(dst_id, p + 1); }
            else { i32_to_le(dst_id, p + 1); }
        }
        else    // must realloc
        {
            bcf_fmt_t *fmt = &line->d.fmt[i];
            kstring_t str = {0,0,0};
            bcf_enc_int1(&str, dst_id);
            bcf_enc_size(&str, fmt->n, fmt->type);
            uint32_t p_off = str.l;
            kputsn((char*)fmt->p, fmt->p_len, &str);
            if( fmt->p_free ) free(fmt->p - fmt->p_off);
            fmt->p_off = p_off;
            fmt->p = (uint8_t*)str.s + fmt->p_off;
            fmt->p_free = 1;
            line->d.indiv_dirty = 1;
        }
    }
    return 0;
}

bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr)
{
    bcf_hdr_t *hout = bcf_hdr_init("r");
    if (!hout) {
        hts_log_error("Failed to allocate bcf header");
        return NULL;
    }
    kstring_t htxt = {0,0,0};
    if (bcf_hdr_format(hdr, 1, &htxt) < 0) {
        free(htxt.s);
        return NULL;
    }
    if ( bcf_hdr_parse(hout, htxt.s) < 0 ) {
        bcf_hdr_destroy(hout);
        hout = NULL;
    }
    free(htxt.s);
    return hout;
}

bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap)
{
    void *names_hash = khash_str2int_init();
    kstring_t htxt = {0,0,0};
    kstring_t str = {0,0,0};
    bcf_hdr_t *h = bcf_hdr_init("w");
    int r = 0;
    if (!h || !names_hash) {
        hts_log_error("Failed to allocate bcf header");
        goto err;
    }
    if (bcf_hdr_format(h0, 1, &htxt) < 0) {
        hts_log_error("Failed to get header text");
        goto err;
    }
    bcf_hdr_set_version(h,bcf_hdr_get_version(h0));
    int j;
    for (j=0; j<n; j++) imap[j] = -1;
    if ( bcf_hdr_nsamples(h0) > 0) {
        char *p = find_chrom_header_line(htxt.s);
        int i = 0, end = n? 8 : 7;
        while ((p = strchr(p, '\t')) != 0 && i < end) ++i, ++p;
        if (i != end) {
            hts_log_error("Wrong number of columns in header #CHROM line");
            goto err;
        }
        r |= kputsn(htxt.s, p - htxt.s, &str) < 0;
        for (i = 0; i < n; ++i) {
            if ( khash_str2int_has_key(names_hash,samples[i]) )
            {
                hts_log_error("Duplicate sample name \"%s\"", samples[i]);
                goto err;
            }
            imap[i] = bcf_hdr_id2int(h0, BCF_DT_SAMPLE, samples[i]);
            if (imap[i] < 0) continue;
            r |= kputc('\t', &str) < 0;
            r |= kputs(samples[i], &str) < 0;
            r |= khash_str2int_inc(names_hash,samples[i]) < 0;
        }
    } else r |= kputsn(htxt.s, htxt.l, &str) < 0;
    while (str.l && (!str.s[str.l-1] || str.s[str.l-1]=='\n') ) str.l--; // kill trailing zeros and newlines
    r |= kputc('\n',&str) < 0;
    if (r) {
        hts_log_error("%s", strerror(errno));
        goto err;
    }
    if ( bcf_hdr_parse(h, str.s) < 0 ) {
        bcf_hdr_destroy(h);
        h = NULL;
    }
    free(str.s);
    free(htxt.s);
    khash_str2int_destroy(names_hash);
    return h;

 err:
    ks_free(&str);
    ks_free(&htxt);
    khash_str2int_destroy(names_hash);
    bcf_hdr_destroy(h);
    return NULL;
}

int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file)
{
    if ( samples && !strcmp("-",samples) ) return 0;            // keep all samples

    int i, narr = bit_array_size(bcf_hdr_nsamples(hdr));
    hdr->keep_samples = (uint8_t*) calloc(narr,1);
    if (!hdr->keep_samples) return -1;

    hdr->nsamples_ori = bcf_hdr_nsamples(hdr);
    if ( !samples )
    {
        // exclude all samples
        khint_t k;
        vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_SAMPLE], *new_dict;
        new_dict = kh_init(vdict);
        if (!new_dict) return -1;

        bcf_hdr_nsamples(hdr) = 0;

        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
        kh_destroy(vdict, d);
        hdr->dict[BCF_DT_SAMPLE] = new_dict;
        if (bcf_hdr_sync(hdr) < 0) return -1;

        return 0;
    }

    if ( samples[0]=='^' )
        for (i=0; i<bcf_hdr_nsamples(hdr); i++) bit_array_set(hdr->keep_samples,i);

    int idx, n, ret = 0;
    char **smpls = hts_readlist(samples[0]=='^'?samples+1:samples, is_file, &n);
    if ( !smpls ) return -1;
    for (i=0; i<n; i++)
    {
        idx = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,smpls[i]);
        if ( idx<0 )
        {
            if ( !ret ) ret = i+1;
            continue;
        }
        assert( idx<bcf_hdr_nsamples(hdr) );
        if (  samples[0]=='^' )
            bit_array_clear(hdr->keep_samples, idx);
        else
            bit_array_set(hdr->keep_samples, idx);
    }
    for (i=0; i<n; i++) free(smpls[i]);
    free(smpls);

    bcf_hdr_nsamples(hdr) = 0;
    for (i=0; i<hdr->nsamples_ori; i++)
        if ( bit_array_test(hdr->keep_samples,i) ) bcf_hdr_nsamples(hdr)++;

    if ( !bcf_hdr_nsamples(hdr) ) { free(hdr->keep_samples); hdr->keep_samples=NULL; }
    else
    {
        // Make new list and dictionary with desired samples
        char **samples = (char**) malloc(sizeof(char*)*bcf_hdr_nsamples(hdr));
        vdict_t *new_dict, *d;
        int k, res;
        if (!samples) return -1;

        new_dict = kh_init(vdict);
        if (!new_dict) {
            free(samples);
            return -1;
        }
        idx = 0;
        for (i=0; i<hdr->nsamples_ori; i++) {
            if ( bit_array_test(hdr->keep_samples,i) ) {
                samples[idx] = hdr->samples[i];
                k = kh_put(vdict, new_dict, hdr->samples[i], &res);
                if (res < 0) {
                    free(samples);
                    kh_destroy(vdict, new_dict);
                    return -1;
                }
                kh_val(new_dict, k) = bcf_idinfo_def;
                kh_val(new_dict, k).id = idx;
                idx++;
            }
        }

        // Delete desired samples from old dictionary, so we don't free them
        d = (vdict_t*)hdr->dict[BCF_DT_SAMPLE];
        for (i=0; i < idx; i++) {
            int k = kh_get(vdict, d, samples[i]);
            if (k < kh_end(d)) kh_del(vdict, d, k);
        }

        // Free everything else
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
        kh_destroy(vdict, d);
        hdr->dict[BCF_DT_SAMPLE] = new_dict;

        free(hdr->samples);
        hdr->samples = samples;

        if (bcf_hdr_sync(hdr) < 0)
            return -1;
    }

    return ret;
}

int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap)
{
    kstring_t ind;
    ind.s = 0; ind.l = ind.m = 0;
    if (n) {
        bcf_fmt_t fmt[MAX_N_FMT];
        int i, j;
        uint8_t *ptr = (uint8_t*)v->indiv.s;
        for (i = 0; i < v->n_fmt; ++i)
            ptr = bcf_unpack_fmt_core1(ptr, v->n_sample, &fmt[i]);
        for (i = 0; i < (int)v->n_fmt; ++i) {
            bcf_fmt_t *f = &fmt[i];
            bcf_enc_int1(&ind, f->id);
            bcf_enc_size(&ind, f->n, f->type);
            for (j = 0; j < n; ++j)
                if (imap[j] >= 0) kputsn((char*)(f->p + imap[j] * f->size), f->size, &ind);
        }
        for (i = j = 0; j < n; ++j) if (imap[j] >= 0) ++i;
        v->n_sample = i;
    } else v->n_sample = 0;
    if ( !v->n_sample ) v->n_fmt = 0;
    free(v->indiv.s);
    v->indiv = ind;
    v->unpacked &= ~BCF_UN_FMT;    // only BCF is ready for output, VCF will need to unpack again
    return 0;
}

int bcf_is_snp(bcf1_t *v)
{
    int i;
    bcf_unpack(v, BCF_UN_STR);
    for (i = 0; i < v->n_allele; ++i)
    {
        if ( v->d.allele[i][1]==0 && v->d.allele[i][0]!='*' ) continue;

        // mpileup's <X> allele, see also below. This is not completely satisfactory,
        // a general library is here narrowly tailored to fit samtools.
        if ( v->d.allele[i][0]=='<' && v->d.allele[i][1]=='X' && v->d.allele[i][2]=='>' ) continue;
        if ( v->d.allele[i][0]=='<' && v->d.allele[i][1]=='*' && v->d.allele[i][2]=='>' ) continue;

        break;
    }
    return i == v->n_allele;
}

static void bcf_set_variant_type(const char *ref, const char *alt, bcf_variant_t *var)
{
    if ( *alt == '*' && !alt[1] ) { var->n = 0; var->type = VCF_OVERLAP; return; }  // overlapping variant

    // The most frequent case
    if ( !ref[1] && !alt[1] )
    {
        if ( *alt == '.' || *ref==*alt ) { var->n = 0; var->type = VCF_REF; return; }
        if ( *alt == 'X' ) { var->n = 0; var->type = VCF_REF; return; }  // mpileup's X allele shouldn't be treated as variant
        var->n = 1; var->type = VCF_SNP; return;
    }
    if ( alt[0]=='<' )
    {
        if ( alt[1]=='X' && alt[2]=='>' ) { var->n = 0; var->type = VCF_REF; return; }  // mpileup's X allele shouldn't be treated as variant
        if ( alt[1]=='*' && alt[2]=='>' ) { var->n = 0; var->type = VCF_REF; return; }
        if ( !strcmp("NON_REF>",alt+1) ) { var->n = 0; var->type = VCF_REF; return; }
        var->type = VCF_OTHER;
        return;
    }

    const char *r = ref, *a = alt;
    while (*r && *a && toupper_c(*r)==toupper_c(*a) ) { r++; a++; }     // unfortunately, matching REF,ALT case is not guaranteed

    if ( *a && !*r )
    {
        if ( *a==']' || *a=='[' ) { var->type = VCF_BND; return; }
        while ( *a ) a++;
        var->n = (a-alt)-(r-ref); var->type = VCF_INDEL; return;
    }
    else if ( *r && !*a )
    {
        while ( *r ) r++;
        var->n = (a-alt)-(r-ref); var->type = VCF_INDEL; return;
    }
    else if ( !*r && !*a )
    {
        var->n = 0; var->type = VCF_REF; return;
    }

    const char *re = r, *ae = a;
    while ( re[1] ) re++;
    while ( ae[1] ) ae++;
    while ( re>r && ae>a && toupper_c(*re)==toupper_c(*ae) ) { re--; ae--; }
    if ( ae==a )
    {
        if ( re==r ) { var->n = 1; var->type = VCF_SNP; return; }
        var->n = -(re-r);
        if ( toupper_c(*re)==toupper_c(*ae) ) { var->type = VCF_INDEL; return; }
        var->type = VCF_OTHER; return;
    }
    else if ( re==r )
    {
        var->n = ae-a;
        if ( toupper_c(*re)==toupper_c(*ae) ) { var->type = VCF_INDEL; return; }
        var->type = VCF_OTHER; return;
    }

    var->type = ( re-r == ae-a ) ? VCF_MNP : VCF_OTHER;
    var->n = ( re-r > ae-a ) ? -(re-r+1) : ae-a+1;

    // should do also complex events, SVs, etc...
}

static int bcf_set_variant_types(bcf1_t *b)
{
    if ( !(b->unpacked & BCF_UN_STR) ) bcf_unpack(b, BCF_UN_STR);
    bcf_dec_t *d = &b->d;
    if ( d->n_var < b->n_allele )
    {
        d->var = (bcf_variant_t *) realloc(d->var, sizeof(bcf_variant_t)*b->n_allele);
        d->n_var = b->n_allele;
    }
    int i;
    b->d.var_type = 0;
    d->var[0].type = VCF_REF;
    d->var[0].n    = 0;
    for (i=1; i<b->n_allele; i++)
    {
        bcf_set_variant_type(d->allele[0],d->allele[i], &d->var[i]);
        b->d.var_type |= d->var[i].type;
        //fprintf(stderr,"[set_variant_type] %d   %s %s -> %d %d .. %d\n", b->pos+1,d->allele[0],d->allele[i],d->var[i].type,d->var[i].n, b->d.var_type);
    }
    return 0;
}

int bcf_get_variant_types(bcf1_t *rec)
{
    if ( rec->d.var_type==-1 ) bcf_set_variant_types(rec);
    return rec->d.var_type;
}
int bcf_get_variant_type(bcf1_t *rec, int ith_allele)
{
    if ( rec->d.var_type==-1 ) bcf_set_variant_types(rec);
    return rec->d.var[ith_allele].type;
}

int bcf_update_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
{
    static int negative_rlen_warned = 0;
    int is_end_tag;

    // Is the field already present?
    int i, inf_id = bcf_hdr_id2int(hdr,BCF_DT_ID,key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,inf_id) ) return -1;    // No such INFO field in the header
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    is_end_tag = strcmp(key, "END") == 0;

    for (i=0; i<line->n_info; i++)
        if ( inf_id==line->d.info[i].key ) break;
    bcf_info_t *inf = i==line->n_info ? NULL : &line->d.info[i];

    if ( !n || (type==BCF_HT_STR && !values) )
    {
        if ( n==0 && is_end_tag )
            line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;
        if ( inf )
        {
            // Mark the tag for removal, free existing memory if necessary
            if ( inf->vptr_free )
            {
                free(inf->vptr - inf->vptr_off);
                inf->vptr_free = 0;
            }
            line->d.shared_dirty |= BCF1_DIRTY_INF;
            inf->vptr = NULL;
            inf->vptr_off = inf->vptr_len = 0;
        }
        return 0;
    }

    if (is_end_tag)
    {
        if (n != 1)
        {
            hts_log_error("END info tag should only have one value at %s:%"PRIhts_pos, bcf_seqname_safe(hdr,line), line->pos+1);
            line->errcode |= BCF_ERR_TAG_INVALID;
            return -1;
        }
        if (type != BCF_HT_INT && type != BCF_HT_LONG)
        {
            hts_log_error("Wrong type (%d) for END info tag at %s:%"PRIhts_pos, type, bcf_seqname_safe(hdr,line), line->pos+1);
            line->errcode |= BCF_ERR_TAG_INVALID;
            return -1;
        }
    }

    // Encode the values and determine the size required to accommodate the values
    kstring_t str = {0,0,0};
    bcf_enc_int1(&str, inf_id);
    if ( type==BCF_HT_INT )
        bcf_enc_vint(&str, n, (int32_t*)values, -1);
    else if ( type==BCF_HT_REAL )
        bcf_enc_vfloat(&str, n, (float*)values);
    else if ( type==BCF_HT_FLAG || type==BCF_HT_STR )
    {
        if ( values==NULL )
            bcf_enc_size(&str, 0, BCF_BT_NULL);
        else
            bcf_enc_vchar(&str, strlen((char*)values), (char*)values);
    }
#ifdef VCF_ALLOW_INT64
    else if ( type==BCF_HT_LONG )
    {
        if (n != 1) {
            hts_log_error("Only storing a single BCF_HT_LONG value is supported at %s:%"PRIhts_pos, bcf_seqname_safe(hdr,line), line->pos+1);
            abort();
        }
        bcf_enc_long1(&str, *(int64_t *) values);
    }
#endif
    else
    {
        hts_log_error("The type %d not implemented yet at %s:%"PRIhts_pos, type, bcf_seqname_safe(hdr,line), line->pos+1);
        abort();
    }

    // Is the INFO tag already present
    if ( inf )
    {
        // Is it big enough to accommodate new block?
        if ( str.l <= inf->vptr_len + inf->vptr_off )
        {
            if ( str.l != inf->vptr_len + inf->vptr_off ) line->d.shared_dirty |= BCF1_DIRTY_INF;
            uint8_t *ptr = inf->vptr - inf->vptr_off;
            memcpy(ptr, str.s, str.l);
            free(str.s);
            int vptr_free = inf->vptr_free;
            bcf_unpack_info_core1(ptr, inf);
            inf->vptr_free = vptr_free;
        }
        else
        {
            if ( inf->vptr_free )
                free(inf->vptr - inf->vptr_off);
            bcf_unpack_info_core1((uint8_t*)str.s, inf);
            inf->vptr_free = 1;
            line->d.shared_dirty |= BCF1_DIRTY_INF;
        }
    }
    else
    {
        // The tag is not present, create new one
        line->n_info++;
        hts_expand0(bcf_info_t, line->n_info, line->d.m_info , line->d.info);
        inf = &line->d.info[line->n_info-1];
        bcf_unpack_info_core1((uint8_t*)str.s, inf);
        inf->vptr_free = 1;
        line->d.shared_dirty |= BCF1_DIRTY_INF;
    }
    line->unpacked |= BCF_UN_INFO;

   if ( n==1 && is_end_tag) {
        hts_pos_t end = type == BCF_HT_INT ? *(int32_t *) values : *(int64_t *) values;
        if ( (type == BCF_HT_INT && end!=bcf_int32_missing) || (type == BCF_HT_LONG && end!=bcf_int64_missing) )
        {
            if ( end <= line->pos )
            {
                if ( !negative_rlen_warned )
                {
                    hts_log_warning("INFO/END=%"PRIhts_pos" is smaller than POS at %s:%"PRIhts_pos,end,bcf_seqname_safe(hdr,line),line->pos+1);
                    negative_rlen_warned = 1;
                }
                line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;
            }
            else
                line->rlen = end - line->pos;
        }
    }
    return 0;
}

int bcf_update_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char **values, int n)
{
    if ( !n )
        return bcf_update_format(hdr,line,key,NULL,0,BCF_HT_STR);

    int i, max_len = 0;
    for (i=0; i<n; i++)
    {
        int len = strlen(values[i]);
        if ( len > max_len ) max_len = len;
    }
    char *out = (char*) malloc(max_len*n);
    if ( !out ) return -2;
    for (i=0; i<n; i++)
    {
        char *dst = out+i*max_len;
        const char *src = values[i];
        int j = 0;
        while ( src[j] ) { dst[j] = src[j]; j++; }
        for (; j<max_len; j++) dst[j] = 0;
    }
    int ret = bcf_update_format(hdr,line,key,out,max_len*n,BCF_HT_STR);
    free(out);
    return ret;
}

int bcf_update_format(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
{
    // Is the field already present?
    int i, fmt_id = bcf_hdr_id2int(hdr,BCF_DT_ID,key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
    {
        if ( !n ) return 0;
        return -1;  // the key not present in the header
    }

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==fmt_id ) break;
    bcf_fmt_t *fmt = i==line->n_fmt ? NULL : &line->d.fmt[i];

    if ( !n )
    {
        if ( fmt )
        {
            // Mark the tag for removal, free existing memory if necessary
            if ( fmt->p_free )
            {
                free(fmt->p - fmt->p_off);
                fmt->p_free = 0;
            }
            line->d.indiv_dirty = 1;
            fmt->p = NULL;
        }
        return 0;
    }

    line->n_sample = bcf_hdr_nsamples(hdr);
    int nps = n / line->n_sample;  // number of values per sample
    assert( nps && nps*line->n_sample==n );     // must be divisible by n_sample

    // Encode the values and determine the size required to accommodate the values
    kstring_t str = {0,0,0};
    bcf_enc_int1(&str, fmt_id);
    assert(values != NULL);
    if ( type==BCF_HT_INT )
        bcf_enc_vint(&str, n, (int32_t*)values, nps);
    else if ( type==BCF_HT_REAL )
    {
        bcf_enc_size(&str, nps, BCF_BT_FLOAT);
        serialize_float_array(&str, nps*line->n_sample, (float *) values);
    }
    else if ( type==BCF_HT_STR )
    {
        bcf_enc_size(&str, nps, BCF_BT_CHAR);
        kputsn((char*)values, nps*line->n_sample, &str);
    }
    else
    {
        hts_log_error("The type %d not implemented yet at %s:%"PRIhts_pos, type, bcf_seqname_safe(hdr,line), line->pos+1);
        abort();
    }

    if ( !fmt )
    {
        // Not present, new format field
        line->n_fmt++;
        hts_expand0(bcf_fmt_t, line->n_fmt, line->d.m_fmt, line->d.fmt);

        // Special case: VCF specification requires that GT is always first
        if ( line->n_fmt > 1 && key[0]=='G' && key[1]=='T' && !key[2] )
        {
            for (i=line->n_fmt-1; i>0; i--)
                line->d.fmt[i] = line->d.fmt[i-1];
            fmt = &line->d.fmt[0];
        }
        else
            fmt = &line->d.fmt[line->n_fmt-1];
        bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
        line->d.indiv_dirty = 1;
        fmt->p_free = 1;
    }
    else
    {
        // The tag is already present, check if it is big enough to accommodate the new block
        if ( str.l <= fmt->p_len + fmt->p_off )
        {
            // good, the block is big enough
            if ( str.l != fmt->p_len + fmt->p_off ) line->d.indiv_dirty = 1;
            uint8_t *ptr = fmt->p - fmt->p_off;
            memcpy(ptr, str.s, str.l);
            free(str.s);
            int p_free = fmt->p_free;
            bcf_unpack_fmt_core1(ptr, line->n_sample, fmt);
            fmt->p_free = p_free;
        }
        else
        {
            if ( fmt->p_free )
                free(fmt->p - fmt->p_off);
            bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
            fmt->p_free = 1;
            line->d.indiv_dirty = 1;
        }
    }
    line->unpacked |= BCF_UN_FMT;
    return 0;
}


int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n)
{
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    line->d.shared_dirty |= BCF1_DIRTY_FLT;
    line->d.n_flt = n;
    if ( !n ) return 0;
    hts_expand(int, line->d.n_flt, line->d.m_flt, line->d.flt);
    int i;
    for (i=0; i<n; i++)
        line->d.flt[i] = flt_ids[i];
    return 0;
}

int bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id)
{
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    int i;
    for (i=0; i<line->d.n_flt; i++)
        if ( flt_id==line->d.flt[i] ) break;
    if ( i<line->d.n_flt ) return 0;    // this filter is already set
    line->d.shared_dirty |= BCF1_DIRTY_FLT;
    if ( flt_id==0 )    // set to PASS
        line->d.n_flt = 1;
    else if ( line->d.n_flt==1 && line->d.flt[0]==0 )
        line->d.n_flt = 1;
    else
        line->d.n_flt++;
    hts_expand(int, line->d.n_flt, line->d.m_flt, line->d.flt);
    line->d.flt[line->d.n_flt-1] = flt_id;
    return 1;
}
int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int pass)
{
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    int i;
    for (i=0; i<line->d.n_flt; i++)
        if ( flt_id==line->d.flt[i] ) break;
    if ( i==line->d.n_flt ) return 0;   // the filter is not present
    line->d.shared_dirty |= BCF1_DIRTY_FLT;
    if ( i!=line->d.n_flt-1 ) memmove(line->d.flt+i,line->d.flt+i+1,(line->d.n_flt-i-1)*sizeof(*line->d.flt));
    line->d.n_flt--;
    if ( !line->d.n_flt && pass ) bcf_add_filter(hdr,line,0);
    return 0;
}

int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter)
{
    if ( filter[0]=='.' && !filter[1] ) filter = "PASS";
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, filter);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FLT,id) ) return -1;  // not defined in the header

    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( id==0 && !line->d.n_flt) return 1; // PASS

    int i;
    for (i=0; i<line->d.n_flt; i++)
        if ( line->d.flt[i]==id ) return 1;
    return 0;
}

static inline int _bcf1_sync_alleles(const bcf_hdr_t *hdr, bcf1_t *line, int nals)
{
    line->d.shared_dirty |= BCF1_DIRTY_ALS;

    line->n_allele = nals;
    hts_expand(char*, line->n_allele, line->d.m_allele, line->d.allele);

    char *als = line->d.als;
    int n = 0;
    while (n<nals)
    {
        line->d.allele[n] = als;
        while ( *als ) als++;
        als++;
        n++;
    }

    // Update REF length. Note that END is 1-based while line->pos 0-based
    bcf_info_t *end_info = bcf_get_info(hdr,line,"END");
    if ( end_info )
    {
        if ( end_info->type==BCF_HT_INT && end_info->v1.i==bcf_int32_missing ) end_info = NULL;
        else if ( end_info->type==BCF_HT_LONG && end_info->v1.i==bcf_int64_missing ) end_info = NULL;
    }
    if ( end_info && end_info->v1.i > line->pos )
        line->rlen = end_info->v1.i - line->pos;
    else if ( nals > 0 )
        line->rlen = strlen(line->d.allele[0]);
    else
        line->rlen = 0;

    return 0;
}
int bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
{
    if ( !(line->unpacked & BCF_UN_STR) ) bcf_unpack(line, BCF_UN_STR);
    kstring_t tmp = {0,0,0};
    char *free_old = NULL;

    // If the supplied alleles are not pointers to line->d.als, the existing block can be reused.
    int i;
    for (i=0; i<nals; i++)
        if ( alleles[i]>=line->d.als && alleles[i]<line->d.als+line->d.m_als ) break;
    if ( i==nals )
    {
        // all alleles point elsewhere, reuse the existing block
        tmp.l = 0; tmp.s = line->d.als; tmp.m = line->d.m_als;
    }
    else
        free_old = line->d.als;

    for (i=0; i<nals; i++)
    {
        kputs(alleles[i], &tmp);
        kputc(0, &tmp);
    }
    line->d.als = tmp.s; line->d.m_als = tmp.m;
    free(free_old);
    return _bcf1_sync_alleles(hdr,line,nals);
}

int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string)
{
    if ( !(line->unpacked & BCF_UN_STR) ) bcf_unpack(line, BCF_UN_STR);
    kstring_t tmp;
    tmp.l = 0; tmp.s = line->d.als; tmp.m = line->d.m_als;
    kputs(alleles_string, &tmp);
    line->d.als = tmp.s; line->d.m_als = tmp.m;

    int nals = 1;
    char *t = line->d.als;
    while (*t)
    {
        if ( *t==',' ) { *t = 0; nals++; }
        t++;
    }
    return _bcf1_sync_alleles(hdr, line, nals);
}

int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)
{
    if ( !(line->unpacked & BCF_UN_STR) ) bcf_unpack(line, BCF_UN_STR);
    kstring_t tmp;
    tmp.l = 0; tmp.s = line->d.id; tmp.m = line->d.m_id;
    if ( id )
        kputs(id, &tmp);
    else
        kputs(".", &tmp);
    line->d.id = tmp.s; line->d.m_id = tmp.m;
    line->d.shared_dirty |= BCF1_DIRTY_ID;
    return 0;
}

int bcf_add_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)
{
    if ( !id ) return 0;
    if ( !(line->unpacked & BCF_UN_STR) ) bcf_unpack(line, BCF_UN_STR);

    kstring_t tmp;
    tmp.l = 0; tmp.s = line->d.id; tmp.m = line->d.m_id;

    int len = strlen(id);
    char *dst = line->d.id;
    while ( *dst && (dst=strstr(dst,id)) )
    {
        if ( dst[len]!=0 && dst[len]!=';' ) dst++;              // a prefix, not a match
        else if ( dst==line->d.id || dst[-1]==';' ) return 0;   // already present
        dst++;  // a suffix, not a match
    }
    if ( line->d.id && (line->d.id[0]!='.' || line->d.id[1]) )
    {
        tmp.l = strlen(line->d.id);
        kputc(';',&tmp);
    }
    kputs(id,&tmp);

    line->d.id = tmp.s; line->d.m_id = tmp.m;
    line->d.shared_dirty |= BCF1_DIRTY_ID;
    return 0;

}

bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) ) return NULL;   // no such FMT field in the header
    return bcf_get_fmt_id(line, id);
}

bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,id) ) return NULL;   // no such INFO field in the header
    return bcf_get_info_id(line, id);
}

bcf_fmt_t *bcf_get_fmt_id(bcf1_t *line, const int id)
{
    int i;
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);
    for (i=0; i<line->n_fmt; i++)
    {
        if ( line->d.fmt[i].id==id ) return &line->d.fmt[i];
    }
    return NULL;
}

bcf_info_t *bcf_get_info_id(bcf1_t *line, const int id)
{
    int i;
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);
    for (i=0; i<line->n_info; i++)
    {
        if ( line->d.info[i].key==id ) return &line->d.info[i];
    }
    return NULL;
}


int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)
{
    int i, ret = -4, tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id) ) return -1;    // no such INFO field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_INFO,tag_id)!=(type & 0xff) ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key==tag_id ) break;
    if ( i==line->n_info ) return ( type==BCF_HT_FLAG ) ? 0 : -3;       // the tag is not present in this record
    if ( type==BCF_HT_FLAG ) return 1;

    bcf_info_t *info = &line->d.info[i];
    if ( !info->vptr ) return -3;           // the tag was marked for removal
    if ( type==BCF_HT_STR )
    {
        if ( *ndst < info->len+1 )
        {
            *ndst = info->len + 1;
            *dst  = realloc(*dst, *ndst);
        }
        memcpy(*dst,info->vptr,info->len);
        ((uint8_t*)*dst)[info->len] = 0;
        return info->len;
    }

    // Make sure the buffer is big enough
    int size1;
    switch (type) {
        case BCF_HT_INT:  size1 = sizeof(int32_t); break;
        case BCF_HT_LONG: size1 = sizeof(int64_t); break;
        case BCF_HT_REAL: size1 = sizeof(float); break;
        default:
            hts_log_error("Unexpected output type %d at %s:%"PRIhts_pos, type, bcf_seqname_safe(hdr,line), line->pos+1);
            return -2;
    }
    if ( *ndst < info->len )
    {
        *ndst = info->len;
        *dst  = realloc(*dst, *ndst * size1);
    }

    #define BRANCH(type_t, convert, is_missing, is_vector_end, set_missing, set_regular, out_type_t) do { \
        out_type_t *tmp = (out_type_t *) *dst; \
        int j; \
        for (j=0; j<info->len; j++) \
        { \
            type_t p = convert(info->vptr + j * sizeof(type_t)); \
            if ( is_vector_end ) break; \
            if ( is_missing ) set_missing; \
            else set_regular; \
            tmp++; \
        } \
        ret = j; \
    } while (0)
    switch (info->type) {
        case BCF_BT_INT8:
            if (type == BCF_HT_LONG) {
                BRANCH(int8_t,  le_to_i8,  p==bcf_int8_missing,  p==bcf_int8_vector_end,  *tmp=bcf_int64_missing, *tmp=p, int64_t);
            } else {
                BRANCH(int8_t,  le_to_i8,  p==bcf_int8_missing,  p==bcf_int8_vector_end,  *tmp=bcf_int32_missing, *tmp=p, int32_t);
            }
            break;
        case BCF_BT_INT16:
            if (type == BCF_HT_LONG) {
                BRANCH(int16_t, le_to_i16, p==bcf_int16_missing, p==bcf_int16_vector_end, *tmp=bcf_int64_missing, *tmp=p, int64_t);
            } else {
                BRANCH(int16_t, le_to_i16, p==bcf_int16_missing, p==bcf_int16_vector_end, *tmp=bcf_int32_missing, *tmp=p, int32_t);
            }
            break;
        case BCF_BT_INT32:
            if (type == BCF_HT_LONG) {
                BRANCH(int32_t, le_to_i32, p==bcf_int32_missing, p==bcf_int32_vector_end, *tmp=bcf_int64_missing, *tmp=p, int64_t); break;
            } else {
                BRANCH(int32_t, le_to_i32, p==bcf_int32_missing, p==bcf_int32_vector_end, *tmp=bcf_int32_missing, *tmp=p, int32_t); break;
            }
        case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, p==bcf_float_missing, p==bcf_float_vector_end, bcf_float_set_missing(*tmp), bcf_float_set(tmp, p), float); break;
        default: hts_log_error("Unexpected type %d at %s:%"PRIhts_pos, info->type, bcf_seqname_safe(hdr,line), line->pos+1); return -2;
    }
    #undef BRANCH
    return ret;  // set by BRANCH
}

int bcf_get_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst)
{
    int i,tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) ) return -1;    // no such FORMAT field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=BCF_HT_STR ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==tag_id ) break;
    if ( i==line->n_fmt ) return -3;                               // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];
    if ( !fmt->p ) return -3;                                      // the tag was marked for removal

    int nsmpl = bcf_hdr_nsamples(hdr);
    if ( !*dst )
    {
        *dst = (char**) malloc(sizeof(char*)*nsmpl);
        if ( !*dst ) return -4;     // could not alloc
        (*dst)[0] = NULL;
    }
    int n = (fmt->n+1)*nsmpl;
    if ( *ndst < n )
    {
        (*dst)[0] = realloc((*dst)[0], n);
        if ( !(*dst)[0] ) return -4;    // could not alloc
        *ndst = n;
    }
    for (i=0; i<nsmpl; i++)
    {
        uint8_t *src = fmt->p + i*fmt->n;
        uint8_t *tmp = (uint8_t*)(*dst)[0] + i*(fmt->n+1);
        memcpy(tmp,src,fmt->n);
        tmp[fmt->n] = 0;
        (*dst)[i] = (char*) tmp;
    }
    return n;
}

int bcf_get_format_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)
{
    int i,j, tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) ) return -1;    // no such FORMAT field in the header
    if ( tag[0]=='G' && tag[1]=='T' && tag[2]==0 )
    {
        // Ugly: GT field is considered to be a string by the VCF header but BCF represents it as INT.
        if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=BCF_HT_STR ) return -2;
    }
    else if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=type ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==tag_id ) break;
    if ( i==line->n_fmt ) return -3;                               // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];
    if ( !fmt->p ) return -3;                                      // the tag was marked for removal

    if ( type==BCF_HT_STR )
    {
        int n = fmt->n*bcf_hdr_nsamples(hdr);
        if ( *ndst < n )
        {
            *dst  = realloc(*dst, n);
            if ( !*dst ) return -4;     // could not alloc
            *ndst = n;
        }
        memcpy(*dst,fmt->p,n);
        return n;
    }

    // Make sure the buffer is big enough
    int nsmpl = bcf_hdr_nsamples(hdr);
    int size1 = type==BCF_HT_INT ? sizeof(int32_t) : sizeof(float);
    if ( *ndst < fmt->n*nsmpl )
    {
        *ndst = fmt->n*nsmpl;
        *dst  = realloc(*dst, *ndst*size1);
        if ( !*dst ) return -4;     // could not alloc
    }

    #define BRANCH(type_t, convert, is_missing, is_vector_end, set_missing, set_vector_end, set_regular, out_type_t) { \
        out_type_t *tmp = (out_type_t *) *dst; \
        uint8_t *fmt_p = fmt->p; \
        for (i=0; i<nsmpl; i++) \
        { \
            for (j=0; j<fmt->n; j++) \
            { \
                type_t p = convert(fmt_p + j * sizeof(type_t)); \
                if ( is_missing ) set_missing; \
                else if ( is_vector_end ) { set_vector_end; break; } \
                else set_regular; \
                tmp++; \
            } \
            for (; j<fmt->n; j++) { set_vector_end; tmp++; } \
            fmt_p += fmt->size; \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, p==bcf_int8_missing,  p==bcf_int8_vector_end,  *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, p==bcf_int16_missing, p==bcf_int16_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, p==bcf_int32_missing, p==bcf_int32_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, p==bcf_float_missing, p==bcf_float_vector_end, bcf_float_set_missing(*tmp), bcf_float_set_vector_end(*tmp), bcf_float_set(tmp, p), float); break;
        default: hts_log_error("Unexpected type %d at %s:%"PRIhts_pos, fmt->type, bcf_seqname_safe(hdr,line), line->pos+1); exit(1);
    }
    #undef BRANCH
    return nsmpl*fmt->n;
}
