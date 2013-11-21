#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "hfile.h"

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "htslib/kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

uint32_t bcf_float_missing    = 0x7F800001;
uint32_t bcf_float_vector_end = 0x7F800002;
uint8_t bcf_type_shift[] = { 0, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static bcf_idinfo_t bcf_idinfo_def = { .info = { 15, 15, 15 }, .hrec = { NULL, NULL, NULL}, .id = -1 };

/*************************
 *** VCF header parser ***
 *************************/

int bcf_hdr_sync(bcf_hdr_t *h);

int bcf_hdr_add_sample(bcf_hdr_t *h, const char *s)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    int ret;
    char *sdup = strdup(s);
    int k = kh_put(vdict, d, sdup, &ret);
    if (ret) { // absent
        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = kh_size(d) - 1;
    } else {
        if (hts_verbose >= 2)
            fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
        free(sdup);
        return -1;
    }
    int n = kh_size(d);
    h->samples = (char**) realloc(h->samples,sizeof(char*)*n);
    h->samples[n-1] = sdup;
    bcf_hdr_sync(h);
    return 0;
}

void bcf_hdr_parse_sample_line(bcf_hdr_t *h, const char *str)
{
    int i = 0;
    const char *p, *q;
    // add samples
    for (p = q = str;; ++q) {
        if (*q != '\t' && *q != 0 && *q != '\n') continue;
        if (++i > 9) {
            char *s = (char*)malloc(q - p + 1);
            strncpy(s, p, q - p);
            s[q - p] = 0;
            bcf_hdr_add_sample(h,s);
            free(s);
        }
        if (*q == 0 || *q == '\n') break;
        p = q + 1;
    }
}

int bcf_hdr_sync(bcf_hdr_t *h)
{
	int i;
	for (i = 0; i < 3; ++i) {
		khint_t k;
		vdict_t *d = (vdict_t*)h->dict[i];
        if ( h->n[i] < kh_size(d) )
            h->id[i] = (bcf_idpair_t*)realloc(h->id[i], kh_size(d) * sizeof(bcf_idpair_t));
		h->n[i] = kh_size(d);
		for (k = kh_begin(d); k != kh_end(d); ++k) {
			if (!kh_exist(d, k)) continue;
			h->id[i][kh_val(d, k).id].key = kh_key(d, k);
			h->id[i][kh_val(d, k).id].val = &kh_val(d, k);
		}
	}
	return 0;
}

void bcf_hrec_destroy(bcf_hrec_t *hrec)
{
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

bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec)
{
    bcf_hrec_t *out = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
    out->type = hrec->type;
    if ( hrec->key ) out->key = strdup(hrec->key);
    if ( hrec->value ) out->value = strdup(hrec->value);
    out->nkeys = hrec->nkeys;
    out->keys = (char**) malloc(sizeof(char*)*hrec->nkeys);
    out->vals = (char**) malloc(sizeof(char*)*hrec->nkeys);
    int i;
    for (i=0; i<hrec->nkeys; i++)
    {
        if ( hrec->keys[i] ) out->keys[i] = strdup(hrec->keys[i]);
        if ( hrec->vals[i] ) out->vals[i] = strdup(hrec->vals[i]);
    }
    return out;
}

void bcf_hrec_debug(bcf_hrec_t *hrec)
{
    printf("key=[%s] value=[%s]", hrec->key, hrec->value?hrec->value:"");
    int i;
    for (i=0; i<hrec->nkeys; i++)
        printf("\t[%s]=[%s]", hrec->keys[i],hrec->vals[i]);
    printf("\n");
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

void bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, int len)
{
    int n = ++hrec->nkeys;
    hrec->keys = (char**) realloc(hrec->keys, sizeof(char*)*n);
    hrec->vals = (char**) realloc(hrec->vals, sizeof(char*)*n);
    assert( len );
    hrec->keys[n-1] = (char*) malloc((len+1)*sizeof(char));
    memcpy(hrec->keys[n-1],str,len);
    hrec->keys[n-1][len] = 0;
    hrec->vals[n-1] = NULL;
}

void bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, int len, int is_quoted)
{
    if ( !str ) { hrec->vals[i] = NULL; return; }
    if ( hrec->vals[i] ) free(hrec->vals[i]);
    if ( is_quoted )
    {
        hrec->vals[i] = (char*) malloc((len+3)*sizeof(char));
        hrec->vals[i][0] = '"';
        memcpy(&hrec->vals[i][1],str,len);
        hrec->vals[i][len+1] = '"';
        hrec->vals[i][len+2] = 0;
    }
    else
    {
        hrec->vals[i] = (char*) malloc((len+1)*sizeof(char));
        memcpy(hrec->vals[i],str,len);
        hrec->vals[i][len] = 0;
    }
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
    const char *p = line;
    if (p[0] != '#' || p[1] != '#') { *len = 0; return NULL; }
    p += 2;

    const char *q = p;
    while ( *q && *q!='=' ) q++;
    int n = q-p;
    if ( *q!='=' || !n ) { *len = q-line+1; return NULL; } // wrong format

    bcf_hrec_t *hrec = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
    hrec->key = (char*) malloc(sizeof(char)*(n+1));
    memcpy(hrec->key,p,n);
    hrec->key[n] = 0;

    p = ++q;
    if ( *p!='<' ) // generic field, e.g. ##samtoolsVersion=0.1.18-r579
    {
        while ( *q && *q!='\n' ) q++;
        hrec->value = (char*) malloc((q-p+1)*sizeof(char));
        memcpy(hrec->value, p, q-p);
        hrec->value[q-p] = 0;
        *len = q-line+1;
        return hrec;
    }

    // structured line, e.g. ##INFO=<ID=PV1,Number=1,Type=Float,Description="P-value for baseQ bias">
    int nopen = 1;
    while ( *q && *q!='\n' && nopen )
    {
        p = ++q;
        while ( *q && *q!='=' ) q++;
        n = q-p;
        if ( *q!='=' || !n ) { *len = q-line+1; bcf_hrec_destroy(hrec); return NULL; } // wrong format
        bcf_hrec_add_key(hrec, p, q-p);
        p = ++q;
        int quoted = *p=='"' ? 1 : 0;
        if ( quoted ) p++, q++;
        while (1)
        {
            if ( !*q ) break;
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
        bcf_hrec_set_val(hrec, hrec->nkeys-1, p, q-p, quoted);
        if ( quoted ) q++;
        if ( *q=='>' ) { nopen--; q++; }
    }
    *len = q-line+1;
    return hrec;
}

// returns: 1 when hdr needs to be synced, 0 otherwise
int bcf_hdr_register_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)
{
    // contig
    int i,j,k, ret;
    char *str;
    if ( !strcmp(hrec->key, "contig") ) 
    {
        hrec->type = BCF_HL_CTG;

        // Get the contig ID ($str) and length ($j)
        i = bcf_hrec_find_key(hrec,"length"); 
        if ( i<0 ) return 0;
        if ( sscanf(hrec->vals[i],"%d",&j)!=1 ) return 0;

        i = bcf_hrec_find_key(hrec,"ID"); 
        if ( i<0 ) return 0; 
        str = strdup(hrec->vals[i]);

        // Register in the dictionary
        vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
        k = kh_put(vdict, d, str, &ret);
        if ( !ret ) { free(str); return 0; }    // already present

        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = kh_size(d) - 1;
        kh_val(d, k).info[0] = i;
        kh_val(d, k).hrec[0] = hrec;
        return 1;
    }

    if ( !strcmp(hrec->key, "INFO") ) hrec->type = BCF_HL_INFO;
    else if ( !strcmp(hrec->key, "FILTER") ) hrec->type = BCF_HL_FLT;
    else if ( !strcmp(hrec->key, "FORMAT") ) hrec->type = BCF_HL_FMT;
    else if ( hrec->nkeys>0 ) { hrec->type = BCF_HL_STR; return 1; }
    else return 0;
    
    // INFO/FILTER/FORMAT
    char *id = NULL;
    int type = -1, num = -1, var = -1; 
    for (i=0; i<hrec->nkeys; i++)
    {
        if ( !strcmp(hrec->keys[i], "ID") ) id = hrec->vals[i];
        else if ( !strcmp(hrec->keys[i], "Type") )
        {
            if ( !strcmp(hrec->vals[i], "Integer") ) type = BCF_HT_INT;
            else if ( !strcmp(hrec->vals[i], "Float") ) type = BCF_HT_REAL;
            else if ( !strcmp(hrec->vals[i], "String") ) type = BCF_HT_STR;
            else if ( !strcmp(hrec->vals[i], "Flag") ) type = BCF_HT_FLAG;
            else
            {
                fprintf(stderr, "[E::%s] The type \"%s\" not supported, assuming \"String\"\n", __func__, hrec->vals[i]);
                type = BCF_HT_STR;
            }
        }
        else if ( !strcmp(hrec->keys[i], "Number") )
        {
            if ( !strcmp(hrec->vals[i],"A") ) var = BCF_VL_A;
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
    uint32_t info = (uint32_t)num<<12 | var<<8 | type<<4 | hrec->type;

    if ( !id ) return 0;
    str = strdup(id);

    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_ID];
    k = kh_put(vdict, d, str, &ret);
    if ( !ret ) 
    { 
        // already present
        free(str);
        if ( kh_val(d, k).hrec[info&0xf] ) return 0;
        kh_val(d, k).info[info&0xf] = info;
        kh_val(d, k).hrec[info&0xf] = hrec;
        return 1;
    }
    kh_val(d, k) = bcf_idinfo_def;
    kh_val(d, k).info[info&0xf] = info;
    kh_val(d, k).hrec[info&0xf] = hrec;
    kh_val(d, k).id = kh_size(d) - 1;
    return 1;
}

int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)
{
    hrec->type = BCF_HL_GEN;
    if ( !bcf_hdr_register_hrec(hdr,hrec) )
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
            if ( !strcmp(hdr->hrec[i]->key,hrec->key) && !strcmp(hdr->hrec[i]->value,hrec->value) ) break;
        }
        if ( i<hdr->nhrec ) 
        {
            bcf_hrec_destroy(hrec);
            return 0;
        }
    }

    // New record, needs to be added
    int n = ++hdr->nhrec;
    hdr->hrec = (bcf_hrec_t**) realloc(hdr->hrec, n*sizeof(bcf_hrec_t*));
    hdr->hrec[n-1] = hrec;

    return hrec->type==BCF_HL_GEN ? 0 : 1;
}

bcf_hrec_t *bcf_hdr_get_hrec(bcf_hdr_t *hdr, int type, const char *id)
{
    vdict_t *d = type==BCF_HL_CTG ? (vdict_t*)hdr->dict[BCF_DT_CTG] : (vdict_t*)hdr->dict[BCF_DT_ID];
    khint_t k = kh_get(vdict, d, id);
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
            fprintf(stderr,"[W::%s] PL should be declared as Number=G\n", __func__);
            PL_warned = 1;
        }
    }
    if ( !GL_warned )
    {
        int id = bcf_hdr_id2int(hdr, BCF_HL_FMT, "GL");
        if ( bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) && bcf_hdr_id2length(hdr,BCF_HL_FMT,id)!=BCF_VL_G )
        {
            fprintf(stderr,"[W::%s] GL should be declared as Number=G\n", __func__);
            PL_warned = 1;
        }
    }
}

int bcf_hdr_parse(bcf_hdr_t *hdr)
{
    int len, needs_sync = 0;
    char *p = hdr->text;

    // Check sanity: "fileformat" string must come as first
    bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr,p,&len);
    if ( !hrec->key || strcasecmp(hrec->key,"fileformat") )
        fprintf(stderr, "[W::%s] The first line should be ##fileformat; is the VCF/BCF header broken?\n", __func__);
    needs_sync += bcf_hdr_add_hrec(hdr, hrec);

    // The filter PASS must appear first in the dictionary
    hrec = bcf_hdr_parse_line(hdr,"##FILTER=<ID=PASS,Description=\"All filters passed\">",&len);
    needs_sync += bcf_hdr_add_hrec(hdr, hrec);

    // Parse the whole header
    while ( (hrec=bcf_hdr_parse_line(hdr,p,&len)) )
    {
        // bcf_hrec_debug(hrec);
        needs_sync += bcf_hdr_add_hrec(hdr, hrec);
        p += len;
    }
    bcf_hdr_parse_sample_line(hdr,p);
    if ( needs_sync ) bcf_hdr_sync(hdr);
    bcf_hdr_check_sanity(hdr);
	return 0;
}

int bcf_hdr_append(bcf_hdr_t *hdr, const char *line)
{
    int len;
    bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr, (char*) line, &len);
    if ( bcf_hdr_add_hrec(hdr, hrec) )
        bcf_hdr_sync(hdr);
    return 0;
}

int bcf_hdr_printf(bcf_hdr_t *hdr, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *line = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(line, n, fmt, ap);
    va_end(ap);

    bcf_hdr_append(hdr, line);

    free(line);
    return 0;
}


/**********************
 *** BCF header I/O ***
 **********************/

bcf_hdr_t *bcf_hdr_init(const char *mode)
{
	int i;
	bcf_hdr_t *h;
	h = (bcf_hdr_t*)calloc(1, sizeof(bcf_hdr_t));
	for (i = 0; i < 3; ++i)
		h->dict[i] = kh_init(vdict);
    if ( strchr(mode,'w') )
    {
        bcf_hdr_append(h, "##fileformat=VCFv4.1");
        // The filter PASS must appear first in the dictionary
        bcf_hdr_append(h, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    }
	return h;
}

void bcf_hdr_destroy(bcf_hdr_t *h)
{
	int i;
	khint_t k;
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
	free(h->mem.s); free(h->text);
	free(h);
}

bcf_hdr_t *bcf_hdr_read(htsFile *hfp)
{
	if (!hfp->is_bin) 
        return vcf_hdr_read(hfp);

    BGZF *fp = hfp->fp.bgzf;
	uint8_t magic[5];
	bcf_hdr_t *h;
	h = bcf_hdr_init("r");
	if ( bgzf_read(fp, magic, 5)<0 ) 
    {
        fprintf(stderr,"[%s:%d %s] Failed to read the header (reading BCF in text mode?)\n", __FILE__,__LINE__,__FUNCTION__);
        return NULL;
    }
	if (strncmp((char*)magic, "BCF\2\2", 5) != 0) 
    {
        if (!strncmp((char*)magic, "BCF", 3)) 
            fprintf(stderr,"[%s:%d %s] invalid BCF2 magic string: only BCFv2.2 is supported.\n", __FILE__,__LINE__,__FUNCTION__);
		else if (hts_verbose >= 2)
			fprintf(stderr, "[E::%s] invalid BCF2 magic string\n", __func__);
		bcf_hdr_destroy(h);
		return 0;
	}
	bgzf_read(fp, &h->l_text, 4);
	h->text = (char*)malloc(h->l_text);
	bgzf_read(fp, h->text, h->l_text);
	bcf_hdr_parse(h);
	return h;
}

int bcf_hdr_write(htsFile *hfp, const bcf_hdr_t *h)
{
    if (!hfp->is_bin) return vcf_hdr_write(hfp, h);
    BGZF *fp = hfp->fp.bgzf;
	if ( bgzf_write(fp, "BCF\2\2", 5) !=5 ) return -1;
	if ( bgzf_write(fp, &h->l_text, 4) !=4 ) return -1;
	if ( bgzf_write(fp, h->text, h->l_text) != h->l_text ) return -1;
    return 0;
}

/********************
 *** BCF site I/O ***
 ********************/

bcf1_t *bcf_init1()
{
	bcf1_t *v;
	v = (bcf1_t*)calloc(1, sizeof(bcf1_t));
	return v;
}

void bcf_clear1(bcf1_t *v)
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
    v->unpack_ptr = NULL;
	v->qual = 0;
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

void bcf_empty1(bcf1_t *v)
{
    bcf_clear1(v);
    free(v->d.id);
    free(v->d.als);
	free(v->d.allele); free(v->d.flt); free(v->d.info); free(v->d.fmt);
	if (v->d.var ) free(v->d.var);
	free(v->shared.s); free(v->indiv.s);
}

void bcf_destroy1(bcf1_t *v)
{
    bcf_empty1(v);
	free(v);
}

static inline int bcf_read1_core(BGZF *fp, bcf1_t *v)
{
	uint32_t x[8];
	int ret;
	if ((ret = bgzf_read(fp, x, 32)) != 32) {
		if (ret == 0) return -1;
		return -2;
	}
    bcf_clear1(v);
	x[0] -= 24; // to exclude six 32-bit integers
	ks_resize(&v->shared, x[0]);
	ks_resize(&v->indiv, x[1]);
	memcpy(v, x + 2, 16);
	v->n_allele = x[6]>>16; v->n_info = x[6]&0xffff;
	v->n_fmt = x[7]>>24; v->n_sample = x[7]&0xffffff;
	v->shared.l = x[0], v->indiv.l = x[1];
	bgzf_read(fp, v->shared.s, v->shared.l);
	bgzf_read(fp, v->indiv.s, v->indiv.l);
	return 0;
}

int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v) 
{ 
    if (!fp->is_bin) return vcf_read(fp,h,v);
    return bcf_read1_core(fp->fp.bgzf, v); 
}

int bcf_readrec(BGZF *fp, void *null, bcf1_t *v, int *tid, int *beg, int *end)
{
	int ret;
	if ((ret = bcf_read1_core(fp, v)) >= 0)
		*tid = v->rid, *beg = v->pos, *end = v->pos + v->rlen;
	return ret;
}

static inline void bcf1_sync_id(bcf1_t *line, kstring_t *str)
{
    // single typed string
    if ( line->d.id && strcmp(line->d.id, ".") ) bcf_enc_vchar(str, strlen(line->d.id), line->d.id);
    else bcf_enc_size(str, 0, BCF_BT_CHAR);
}
static inline void bcf1_sync_alleles(bcf1_t *line, kstring_t *str)
{
    // list of typed strings
    int i;
    for (i=0; i<line->n_allele; i++)
        bcf_enc_vchar(str, strlen(line->d.allele[i]), line->d.allele[i]);
    line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;     // beware: this neglects SV's END tag
}
static inline void bcf1_sync_filter(bcf1_t *line, kstring_t *str)
{
    // typed vector of integers
    if ( line->d.n_flt ) bcf_enc_vint(str, line->d.n_flt, line->d.flt, -1);
    else bcf_enc_vint(str, 0, 0, -1);
}
static inline void bcf1_sync_info(bcf1_t *line, kstring_t *str)
{
    // pairs of typed vectors
    int i, irm = -1;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *info = &line->d.info[i];
        if ( !info->vptr )
        { 
            // marked for removal
            if ( irm < 0 ) irm = i;
            continue;
        }
        kputsn_(info->vptr - info->vptr_off, info->vptr_len + info->vptr_off, str);
        if ( irm >=0 )
        {
            bcf_info_t tmp = line->d.info[irm]; line->d.info[irm] = line->d.info[i]; line->d.info[i] = tmp;
            while ( irm<=i && line->d.info[irm].vptr ) irm++;
        }
    }
    if ( irm>=0 ) line->n_info = irm;
}

static int bcf1_sync(bcf1_t *line)
{
    int i;
    kstring_t tmp = {0,0,0};
    if ( !line->shared.l )
    {
        // New line, get ready for BCF output
        tmp = line->shared;
        bcf1_sync_id(line, &tmp);
        bcf1_sync_alleles(line, &tmp);
        bcf1_sync_filter(line, &tmp);
        bcf1_sync_info(line, &tmp);
        line->shared = tmp;
    }
    else if ( line->d.shared_dirty )
    {
        // The line was edited, update the BCF data block
        uint8_t *ptr, *ptr_ori = (uint8_t *) line->shared.s;
        int type, size;

        // ID: single typed string
        size = bcf_dec_size(ptr_ori, &ptr, &type);
        size = (size << bcf_type_shift[type]) + (ptr - ptr_ori);
        if ( line->d.shared_dirty & BCF1_DIRTY_ID ) 
            bcf1_sync_id(line, &tmp);
        else 
            kputsn_(ptr_ori, size, &tmp);
        ptr_ori += size;

        // REF+ALT: list of typed strings
        uint8_t *_ptr, *_ptr_ori = ptr_ori;
        size = 0;
        for (i=0; i<line->n_allele; i++) 
        {
            int _size = bcf_dec_size(_ptr_ori, &_ptr, &type);
            size += (_size << bcf_type_shift[type]) + (_ptr - _ptr_ori);
            _ptr_ori = _ptr + (_size << bcf_type_shift[type]);
        }
        if ( line->d.shared_dirty & BCF1_DIRTY_ALS ) 
            bcf1_sync_alleles(line, &tmp);
        else
        {
            line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;     // beware: this neglects SV's END tag
            kputsn_(ptr_ori, size, &tmp);
        }
        ptr_ori += size;

        // FILTER: typed vector of integers
        ptr = ptr_ori;
        size = bcf_dec_size(ptr_ori, &ptr, &type);
        for (i=0; i<size; i++) 
            bcf_dec_int1(ptr, type, &ptr);
        size = ptr - ptr_ori;
        if ( line->d.shared_dirty & BCF1_DIRTY_FLT )
            bcf1_sync_filter(line, &tmp);
        else if ( line->d.n_flt ) 
            kputsn_(ptr_ori, size, &tmp);
        else
            bcf_enc_vint(&tmp, 0, 0, -1);
        ptr_ori += size;

        // INFO: pairs of typed vectors
        if ( line->d.shared_dirty & BCF1_DIRTY_INF ) 
            bcf1_sync_info(line, &tmp);
        else
        {
            size = line->shared.l - (size_t)ptr_ori + (size_t)line->shared.s;
            kputsn_(ptr_ori, size, &tmp);
        }
        free(line->shared.s);
        line->shared = tmp;
    }
    if ( line->n_sample && line->n_fmt && (!line->indiv.l || line->d.indiv_dirty) )
    {
        // The genotype fields has changed or are not present
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
    }
    return 0;
}

int bcf_write(htsFile *hfp, const bcf_hdr_t *h, bcf1_t *v)
{
    if ( !hfp->is_bin ) return vcf_write(hfp,h,v);

    if ( v->errcode ) 
    {
        // vcf_parse1() encountered a new contig or tag, undeclared in the
        // header.  At this point, the header must have been printed,
        // proceeding would lead to a broken BCF file. Errors must be checked
        // and cleared by the caller before we can proceed.
        fprintf(stderr,"[%s:%d %s] Unchecked error (%d), exiting.\n", __FILE__,__LINE__,__FUNCTION__,v->errcode);
        exit(1);
    }
    bcf1_sync(v);   // check if the BCF record was modified

    BGZF *fp = hfp->fp.bgzf;
	uint32_t x[8];
	x[0] = v->shared.l + 24; // to include six 32-bit integers
	x[1] = v->indiv.l;
	memcpy(x + 2, v, 16);
	x[6] = (uint32_t)v->n_allele<<16 | v->n_info;
	x[7] = (uint32_t)v->n_fmt<<24 | v->n_sample;
	if ( bgzf_write(fp, x, 32) != 32 ) return -1;
    if ( bgzf_write(fp, v->shared.s, v->shared.l) != v->shared.l ) return -1;
	if ( bgzf_write(fp, v->indiv.s, v->indiv.l) != v->indiv.l ) return -1;
	return 0;
}

/**********************
 *** VCF header I/O ***
 **********************/

bcf_hdr_t *vcf_hdr_read(htsFile *fp)
{
    kstring_t txt, *s = &fp->line;
    bcf_hdr_t *h;
    h = bcf_hdr_init("r");
    txt.l = txt.m = 0; txt.s = 0;
    while (hts_getline(fp, KS_SEP_LINE, s) >= 0) {
        if (s->l == 0) continue;
        if (s->s[0] != '#') {
            if (hts_verbose >= 2)
                fprintf(stderr, "[E::%s] no sample line\n", __func__);
            free(txt.s);
            bcf_hdr_destroy(h);
            return 0;
        }
        if (s->s[1] != '#' && fp->fn_aux) { // insert contigs here
            int dret;
            gzFile f;
            kstream_t *ks;
            kstring_t tmp;
            tmp.l = tmp.m = 0; tmp.s = 0;
            f = gzopen(fp->fn_aux, "r");
            ks = ks_init(f);
            while (ks_getuntil(ks, 0, &tmp, &dret) >= 0) {
                int c;
                kputs("##contig=<ID=", &txt); kputs(tmp.s, &txt);
                ks_getuntil(ks, 0, &tmp, &dret);
                kputs(",length=", &txt); kputw(atol(tmp.s), &txt);
                kputsn(">\n", 2, &txt);
                if (dret != '\n')
                    while ((c = ks_getc(ks)) != '\n' && c != -1); // skip the rest of the line
            }
            free(tmp.s);
            ks_destroy(ks);
            gzclose(f);
        }
        kputsn(s->s, s->l, &txt);
        kputc('\n', &txt);
        if (s->s[1] != '#') break;
    }
    h->l_text = txt.l + 1; // including NULL
    h->text = txt.s;
    if ( !h->text ) 
    {
        fprintf(stderr,"[%s:%d %s] Could not read the header\n", __FILE__,__LINE__,__FUNCTION__);
        return NULL;
    }
    bcf_hdr_parse(h);
    // check tabix index, are all contigs listed in the header? add the missing ones
    tbx_t *idx = tbx_index_load(fp->fn);
    if ( idx )
    {
        int i, n, need_sync = 0;
        const char **names = tbx_seqnames(idx, &n);
        for (i=0; i<n; i++)
        {
            bcf_hrec_t *hrec = bcf_hdr_get_hrec(h, BCF_DT_CTG, (char*) names[i]);
            if ( hrec ) continue;
            hrec = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
            hrec->key = strdup("contig");
            bcf_hrec_add_key(hrec, "ID", strlen("ID"));
            bcf_hrec_set_val(hrec, hrec->nkeys-1, (char*) names[i], strlen(names[i]), 0);
            bcf_hrec_add_key(hrec, "length", strlen("length"));
            bcf_hrec_set_val(hrec, hrec->nkeys-1, "2147483647", strlen("2147483647"), 0);
            bcf_hdr_add_hrec(h, hrec);
            need_sync = 1;
        }
        free(names);
        tbx_destroy(idx);
        if ( need_sync )
        {
            bcf_hdr_sync(h);
            bcf_hdr_fmt_text(h);
        }
    }
    return h;
}

int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname)
{
    int i, n;
    char **lines = hts_readlines(fname, &n);
    if ( !lines ) return 1;
    for (i=0; i<n-1; i++)
    {
        int k;
        bcf_hrec_t *hrec = bcf_hdr_parse_line(hdr,lines[i],&k);
        bcf_hdr_add_hrec(hdr, hrec);
        free(lines[i]);
    }
    bcf_hdr_parse_sample_line(hdr,lines[n-1]);
    free(lines[n-1]);
    free(lines);
    bcf_hdr_sync(hdr);
    bcf_hdr_fmt_text(hdr);
    return 0;
}

void bcf_hdr_fmt_text(bcf_hdr_t *hdr)
{
    int i,j;
    kstring_t txt = {0,0,0};
    for (i=0; i<hdr->nhrec; i++)
    {
        if ( !hdr->hrec[i]->value )
        {
            ksprintf(&txt, "##%s=<", hdr->hrec[i]->key);
            ksprintf(&txt,"%s=%s", hdr->hrec[i]->keys[0], hdr->hrec[i]->vals[0]);
            for (j=1; j<hdr->hrec[i]->nkeys; j++)
                ksprintf(&txt,",%s=%s", hdr->hrec[i]->keys[j], hdr->hrec[i]->vals[j]);
            ksprintf(&txt,">\n");
        }
        else
            ksprintf(&txt,"##%s=%s\n", hdr->hrec[i]->key,hdr->hrec[i]->value);
    }
    ksprintf(&txt,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if ( hdr->n[BCF_DT_SAMPLE] )
    {
        ksprintf(&txt,"\tFORMAT");
        for (i=0; i<hdr->n[BCF_DT_SAMPLE]; i++)
            ksprintf(&txt,"\t%s", hdr->samples[i]);
    }
    ksprintf(&txt,"\n");

    if ( hdr->text ) free(hdr->text);
    hdr->text = txt.s;
    hdr->l_text = txt.l + 1;    // the terminating \0 must be included
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
    int l = h->l_text;
    while (l && h->text[l-1] == 0) --l; // kill the trailing zeros
    int ret;
    if ( fp->is_compressed==1 )
        ret = bgzf_write(fp->fp.bgzf, h->text, l);
    else
        ret = hwrite(fp->fp.hfile, h->text, l);
    return ret<0 ? -1 : 0;
}

/***********************
 *** Typed value I/O ***
 ***********************/

void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
{
	int32_t max = INT32_MIN + 1, min = INT32_MAX;
	int i;
	if (n == 0) bcf_enc_size(s, 0, BCF_BT_NULL);
	else if (n == 1) bcf_enc_int1(s, a[0]);
	else {
		if (wsize <= 0) wsize = n;
		for (i = 0; i < n; ++i) {
			if (a[i] == bcf_int32_missing || a[i] == bcf_int32_vector_end ) continue;
			if (max < a[i]) max = a[i];
			if (min > a[i]) min = a[i];
		}
		if (max <= INT8_MAX && min > bcf_int8_vector_end) {
			bcf_enc_size(s, wsize, BCF_BT_INT8);
			for (i = 0; i < n; ++i)
                if ( a[i]==bcf_int32_vector_end ) kputc(bcf_int8_vector_end, s);
                else if ( a[i]==bcf_int32_missing ) kputc(bcf_int8_missing, s);
                else kputc(a[i], s);
		} else if (max <= INT16_MAX && min > bcf_int16_vector_end) {
			bcf_enc_size(s, wsize, BCF_BT_INT16);
			for (i = 0; i < n; ++i) 
            {
                int16_t x;
                if ( a[i]==bcf_int32_vector_end ) x = bcf_int16_vector_end;
                else if ( a[i]==bcf_int32_missing ) x = bcf_int16_missing;
                else x = a[i];
				kputsn((char*)&x, 2, s);
			}
		} else {
			bcf_enc_size(s, wsize, BCF_BT_INT32);
			for (i = 0; i < n; ++i) {
				int32_t x = a[i];
				kputsn((char*)&x, 4, s);
			}
		}
	}
}

void bcf_enc_vfloat(kstring_t *s, int n, float *a)
{
	bcf_enc_size(s, n, BCF_BT_FLOAT);
	kputsn((char*)a, n << 2, s);
}

void bcf_enc_vchar(kstring_t *s, int l, const char *a)
{
	bcf_enc_size(s, l, BCF_BT_CHAR);
	kputsn(a, l, s);
}

void bcf_fmt_array(kstring_t *s, int n, int type, void *data)
{
	int j = 0;
	if (n == 0) {
		kputc('.', s);
		return;
	}
    if (type == BCF_BT_CHAR) 
    {
        char *p = (char*)data;
        for (j = 0; j < n && *p; ++j, ++p) 
        {
            if ( *p==bcf_str_missing ) kputc('.', s);
            else kputc(*p, s);
        }
    }
    else
    {
        #define BRANCH(type_t, is_missing, is_vector_end, kprint) { \
            type_t *p = (type_t *) data; \
            for (j=0; j<n; j++) \
            { \
                if ( is_vector_end ) break; \
                if ( j ) kputc(',', s); \
                if ( is_missing ) kputc('.', s); \
                else kprint; \
            } \
        }
        switch (type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  kputw(p[j], s)); break;
            case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, kputw(p[j], s)); break;
            case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, kputw(p[j], s)); break;
            case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), ksprintf(s, "%g", p[j])); break;
            default: fprintf(stderr,"todo: type %d\n", type); exit(1); break;
        }
        #undef BRANCH
    }
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
	uint32_t is_gt:1, max_g:15, max_l:16;
	uint32_t y;
	uint8_t *buf;
} fmt_aux_t;

static inline void align_mem(kstring_t *s)
{
	if (s->l&7) {
		uint64_t zero = 0;
		int l = ((s->l + 7)>>3<<3) - s->l;
		kputsn((char*)&zero, l, s);
	}
}

int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
{
	int i = 0;
	char *p, *q, *r, *t;
	fmt_aux_t *fmt = 0;
	kstring_t *str, *mem = (kstring_t*)&h->mem;
	khint_t k;
	ks_tokaux_t aux;

    bcf_clear1(v);
	mem->l = 0;
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
                // Simple error recovery for chromosomes not defined in the header. It will not help when VCF header has
                // been already printed, but will enable tools like vcfcheck to proceed.
                fprintf(stderr, "[W::%s] contig '%s' is not defined in the header\n", __func__, p);
                kstring_t tmp = {0,0,0};
                int l;
                ksprintf(&tmp, "##contig=<ID=%s,length=2147483647>", p);
                bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
                free(tmp.s);
                if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
                k = kh_get(vdict, d, p);
                v->errcode = BCF_ERR_CTG_UNDEF;
			}
            v->rid = kh_val(d, k).id;
		} else if (i == 1) { // POS
			v->pos = atoi(p) - 1;
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
						bcf_enc_vchar(str, r - t, t);
						t = r + 1;
						++v->n_allele;
					}
					if (r == q) break;
				}
			}
		} else if (i == 5) { // QUAL
			if (strcmp(p, ".")) v->qual = atof(p);
			else memcpy(&v->qual, &bcf_float_missing, 4);
            if ( v->max_unpack && !(v->max_unpack>>1) ) return 0; // BCF_UN_STR
		} else if (i == 6) { // FILTER
			if (strcmp(p, ".")) {
				int32_t *a;
				int n_flt = 1, i;
				ks_tokaux_t aux1;
				vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
				// count the number of filters
				if (*(q-1) == ';') *(q-1) = 0;
				for (r = p; *r; ++r)
					if (*r == ';') ++n_flt;
				a = (int32_t*)alloca(n_flt * 4);
				// add filters
				for (t = kstrtok(p, ";", &aux1), i = 0; t; t = kstrtok(0, 0, &aux1)) {
					*(char*)aux1.p = 0;
					k = kh_get(vdict, d, t);
					if (k == kh_end(d)) 
                    {
                        // Simple error recovery for FILTERs not defined in the header. It will not help when VCF header has
                        // been already printed, but will enable tools like vcfcheck to proceed.
                        fprintf(stderr, "[W::%s] FILTER '%s' is not defined in the header\n", __func__, t);
                        kstring_t tmp = {0,0,0};
                        int l;
                        ksprintf(&tmp, "##FILTER=<ID=%s,Description=\"Dummy\">", t);
                        bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
                        free(tmp.s);
                        if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
                        k = kh_get(vdict, d, t);
                        v->errcode = BCF_ERR_TAG_UNDEF;
                    }
					a[i++] = kh_val(d, k).id;
				}
				n_flt = i;
				bcf_enc_vint(str, n_flt, a, -1);
			} else bcf_enc_vint(str, 0, 0, -1);
            if ( v->max_unpack && !(v->max_unpack>>2) ) return 0;    // BCF_UN_FLT
		} else if (i == 7) { // INFO
			char *key;
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
			v->n_info = 0;
			if (strcmp(p, ".")) {
				if (*(q-1) == ';') *(q-1) = 0;
				for (r = key = p;; ++r) {
					int c;
					char *val, *end;
					if (*r != ';' && *r != '=' && *r != 0) continue;
					val = end = 0;
					c = *r; *r = 0;
					if (c == '=') {
						val = r + 1;
						for (end = val; *end != ';' && *end != 0; ++end);
						c = *end; *end = 0;
					} else end = r;
					k = kh_get(vdict, d, key);
					if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_INFO] == 15) 
                    {
                        fprintf(stderr, "[W::%s] INFO '%s' is not defined in the header, assuming Type=String\n", __func__, key);
                        kstring_t tmp = {0,0,0};
                        int l;
                        ksprintf(&tmp, "##INFO=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", key);
                        bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
                        free(tmp.s);
                        if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
                        k = kh_get(vdict, d, key);
                        v->errcode = BCF_ERR_TAG_UNDEF;
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
                        char *t;
                        for (t = val, n_val = 1; *t; ++t) // count the number of values
                            if (*t == ',') ++n_val;
                        if ((y>>4&0xf) == BCF_HT_INT) {
                            int32_t *z;
                            z = (int32_t*)alloca(n_val<<2);
                            for (i = 0, t = val; i < n_val; ++i, ++t)
                                z[i] = strtol(t, &t, 10);
                            bcf_enc_vint(str, n_val, z, -1);
                            if (strcmp(key, "END") == 0) v->rlen = z[0] - v->pos;
                        } else if ((y>>4&0xf) == BCF_HT_REAL) {
                            float *z;
                            z = (float*)alloca(n_val<<2);
                            for (i = 0, t = val; i < n_val; ++i, ++t)
                                z[i] = strtod(t, &t);
                            bcf_enc_vfloat(str, n_val, z);
                        }
                    }
                    if (c == 0) break;
                    r = end;
                    key = r + 1;
                }
            }
            if ( v->max_unpack && !(v->max_unpack>>3) ) return 0; 
		} else if (i == 8) { // FORMAT
			int j, l, m, g;
			ks_tokaux_t aux1;
			vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
			char *end = s->s + s->l;
			// count the number of format fields
			for (r = p, v->n_fmt = 1; *r; ++r)
				if (*r == ':') ++v->n_fmt;
			fmt = (fmt_aux_t*)alloca(v->n_fmt * sizeof(fmt_aux_t));
			// get format information from the dictionary
			for (j = 0, t = kstrtok(p, ":", &aux1); t; t = kstrtok(0, 0, &aux1), ++j) {
				*(char*)aux1.p = 0;
				k = kh_get(vdict, d, t);
				if (k == kh_end(d) || kh_val(d, k).info[BCF_HL_FMT] == 15) {
                    fprintf(stderr, "[W::%s] FORMAT '%s' is not defined in the header, assuming Type=String\n", __func__, t);
                    kstring_t tmp = {0,0,0};
                    int l;
                    ksprintf(&tmp, "##FORMAT=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", t);
                    bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
                    free(tmp.s);
                    if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
                    k = kh_get(vdict, d, t);
                    v->errcode = BCF_ERR_TAG_UNDEF;
                }
                fmt[j].max_l = fmt[j].max_m = fmt[j].max_g = 0;
                fmt[j].key = kh_val(d, k).id;
                fmt[j].is_gt = !strcmp(t, "GT");
                fmt[j].y = h->id[0][fmt[j].key].val->info[BCF_HL_FMT];
			}
            if ( q>=end ) 
            {
                fprintf(stderr,"[%s:%d %s] Error: FORMAT column with no sample columns starting at %s:%d\n", __FILE__,__LINE__,__FUNCTION__,s->s,v->pos+1);
                return -1;
            }
			// compute max
			for (r = q + 1, j = 0, m = l = g = 1, v->n_sample = 0;; ++r, ++l) {
				if (*r == '\t') *r = 0;
				if (*r == ':' || *r == '\0') { // end of a sample
					if (fmt[j].max_m < m) fmt[j].max_m = m;
					if (fmt[j].max_l < l - 1) fmt[j].max_l = l - 1;
					if (fmt[j].is_gt && fmt[j].max_g < g) fmt[j].max_g = g;
					l = 0, m = g = 1;
					if (*r) ++j;
					else j = 0, ++v->n_sample;
				} else if (*r == ',') ++m;
				else if (*r == '|' || *r == '/') ++g;
				if (r == end) break;
			}
			// allocate memory for arrays
			for (j = 0; j < v->n_fmt; ++j) {
				fmt_aux_t *f = &fmt[j];
				if ((f->y>>4&0xf) == BCF_HT_STR) {
					f->size = f->is_gt? f->max_g << 2 : f->max_l;
				} else if ((f->y>>4&0xf) == BCF_HT_REAL || (f->y>>4&0xf) == BCF_HT_INT) {
					f->size = f->max_m << 2;
				} else 
                {
                    fprintf(stderr, "[E::%s] the format type %d currently not supported\n", __func__, f->y>>4&0xf);
                    abort(); // I do not know how to do with Flag in the genotype fields
                }
				align_mem(mem);
				f->offset = mem->l;
				ks_resize(mem, mem->l + v->n_sample * f->size);
				mem->l += v->n_sample * f->size;
			}
			for (j = 0; j < v->n_fmt; ++j)
				fmt[j].buf = (uint8_t*)mem->s + fmt[j].offset;
			// fill the sample fields; at beginning of the loop, t points to the first char of a format
			for (t = q + 1, j = m = 0;; ++t) { // j: fmt id, m: sample id
				fmt_aux_t *z = &fmt[j];
				if ((z->y>>4&0xf) == BCF_HT_STR) {
					if (z->is_gt) { // genotypes
						int32_t is_phased = 0, *x = (int32_t*)(z->buf + z->size * m);
						for (l = 0;; ++t) {
							if (*t == '.') ++t, x[l++] = is_phased;
							else x[l++] = (strtol(t, &t, 10) + 1) << 1 | is_phased;
							is_phased = (*t == '|');
							if (*t == ':' || *t == 0) break;
						}
                        if ( !l ) x[l++] = 0;   // An empty field, insert missing value
                        for (; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
					} else {
						char *x = (char*)z->buf + z->size * m;
						for (r = t, l = 0; *t != ':' && *t; ++t) x[l++] = *t;
						for (; l < z->size; ++l) x[l] = 0;
					}
				} else if ((z->y>>4&0xf) == BCF_HT_INT) {
					int32_t *x = (int32_t*)(z->buf + z->size * m);
					for (l = 0;; ++t) {
						if (*t == '.') x[l++] = bcf_int32_missing, ++t; // ++t to skip "."
						else x[l++] = strtol(t, &t, 10);
						if (*t == ':' || *t == 0) break;
					}
                    if ( !l ) x[l++] = bcf_int32_missing;
					for (; l < z->size>>2; ++l) x[l] = bcf_int32_vector_end;
				} else if ((z->y>>4&0xf) == BCF_HT_REAL) {
					float *x = (float*)(z->buf + z->size * m);
					for (l = 0;; ++t) {
						if (*t == '.' && !isdigit(t[1])) bcf_float_set_missing(x[l++]), ++t; // ++t to skip "."
						else x[l++] = strtod(t, &t);
						if (*t == ':' || *t == 0) break;
					}
                    if ( !l ) bcf_float_set_missing(x[l++]);    // An empty field, insert missing value 
					for (; l < z->size>>2; ++l) bcf_float_set_vector_end(x[l]);
				} else abort();
				if (*t == 0) {
					for (++j; j < v->n_fmt; ++j) { // fill end-of-vector values
						z = &fmt[j];
						if ((z->y>>4&0xf) == BCF_HT_STR) {
							if (z->is_gt) {
								int32_t *x = (int32_t*)(z->buf + z->size * m);
                                x[0] = bcf_int32_missing;
								for (l = 1; l != z->size>>2; ++l) x[l] = bcf_int32_vector_end;
							} else {
								char *x = (char*)z->buf + z->size * m;
								for (l = 0; l != z->size; ++l) x[l] = 0;
							}
						} else if ((z->y>>4&0xf) == BCF_HT_INT) {
							int32_t *x = (int32_t*)(z->buf + z->size * m);
                            x[0] = bcf_int32_missing;
							for (l = 1; l != z->size>>2; ++l) x[l] = bcf_int32_vector_end;
						} else if ((z->y>>4&0xf) == BCF_HT_REAL) {
							float *x = (float*)(z->buf + z->size * m);
                            bcf_float_set_missing(x[0]);
							for (l = 1; l != z->size>>2; ++l) bcf_float_set_vector_end(x[l]);
						}
					}
					if (t == end) break;
					++m, j = 0;
				} else if (*t == ':') ++j;
			}
			break;
		}
	}
	// write individual genotype information
	str = &v->indiv;
	if (v->n_sample > 0) {
		for (i = 0; i < v->n_fmt; ++i) {
			fmt_aux_t *z = &fmt[i];
			bcf_enc_int1(str, z->key);
			if ((z->y>>4&0xf) == BCF_HT_STR && !z->is_gt) {
				bcf_enc_size(str, z->size, BCF_BT_CHAR);
				kputsn((char*)z->buf, z->size * v->n_sample, str);
			} else if ((z->y>>4&0xf) == BCF_HT_INT || z->is_gt) {
				bcf_enc_vint(str, (z->size>>2) * v->n_sample, (int32_t*)z->buf, z->size>>2);
			} else {
				bcf_enc_size(str, z->size>>2, BCF_BT_FLOAT);
				kputsn((char*)z->buf, z->size * v->n_sample, str);
			}
		}
	}
	return 0;
}

int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
    int ret;
    ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
    if (ret < 0) return -1;
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
        else if (info->type == BCF_BT_INT32) info->v1.i = *(int32_t*)ptr;
        else if (info->type == BCF_BT_FLOAT) info->v1.f = *(float*)ptr;
        else if (info->type == BCF_BT_INT16) info->v1.i = *(int16_t*)ptr;
    }
    ptr += info->len << bcf_type_shift[info->type];
    info->vptr_len = ptr - info->vptr;
    return ptr;
}

int bcf_unpack(bcf1_t *b, int which)
{
    if ( !b->shared.l ) return 0; // Building a new BCF record from scratch
	uint8_t *ptr = (uint8_t*)b->shared.s;
	int *offset, i;
	bcf_dec_t *d = &b->d;
	if (which & BCF_UN_FLT) which |= BCF_UN_STR;
	if (which & BCF_UN_INFO) which |= BCF_UN_SHR;
	if ((which&BCF_UN_STR) && !(b->unpacked&BCF_UN_STR)) 
    { 
        kstring_t tmp; 

        // ID
        tmp.l = 0; tmp.s = d->id; tmp.m = d->m_id;
        ptr = bcf_fmt_sized_array(&tmp, ptr); 
        kputc('\0', &tmp);
        d->id = tmp.s; d->m_id = tmp.m;

        // REF and ALT are in a single block (d->als) and d->alleles are pointers into this block
        tmp.l = 0; tmp.s = d->als; tmp.m = d->m_als;
        offset = (int*)alloca(b->n_allele * sizeof(int));
        offset = (int*)alloca(b->n_allele * sizeof(int));
        for (i = 0; i < b->n_allele; ++i) {
            offset[i] = tmp.l;
            ptr = bcf_fmt_sized_array(&tmp, ptr);
            kputc('\0', &tmp);
        }
        d->als = tmp.s; d->m_als = tmp.m;

        hts_expand(char*, b->n_allele, d->m_allele, d->allele); // NM: hts_expand() is a macro
        for (i = 0; i < b->n_allele; ++i)
            d->allele[i] = d->als + offset[i];
        b->unpack_ptr = ptr;
        b->unpacked |= BCF_UN_STR;
	}
	if ((which&BCF_UN_FLT) && !(b->unpacked&BCF_UN_FLT)) { // FILTER
		ptr = b->unpack_ptr;
		if (*ptr>>4) {
			int type;
			d->n_flt = bcf_dec_size(ptr, &ptr, &type);
			hts_expand(int, d->n_flt, d->m_flt, d->flt);
			for (i = 0; i < d->n_flt; ++i)
				d->flt[i] = bcf_dec_int1(ptr, type, &ptr);
		} else ++ptr, d->n_flt = 0;
		b->unpack_ptr = ptr;
		b->unpacked |= BCF_UN_FLT;
	}
	if ((which&BCF_UN_INFO) && !(b->unpacked&BCF_UN_INFO)) { // INFO
		ptr = b->unpack_ptr;
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
	bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
	kputs(h->id[BCF_DT_CTG][v->rid].key, s); // CHROM
	kputc('\t', s); kputw(v->pos + 1, s); // POS
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
	if (memcmp(&v->qual, &bcf_float_missing, 4) == 0) kputc('.', s); // QUAL
	else ksprintf(s, "%g", v->qual);
	kputc('\t', s); // FILTER
	if (v->d.n_flt) {
		for (i = 0; i < v->d.n_flt; ++i) {
			if (i) kputc(';', s);
			kputs(h->id[BCF_DT_ID][v->d.flt[i]].key, s);
		}
	} else kputc('.', s);
	kputc('\t', s); // INFO
	if (v->n_info) {
        int first = 1;
		for (i = 0; i < v->n_info; ++i) {
			bcf_info_t *z = &v->d.info[i];
            if ( !z->vptr ) continue;
			if ( !first ) kputc(';', s); first = 0;
			kputs(h->id[BCF_DT_ID][z->key].key, s);
			if (z->len <= 0) continue;
			kputc('=', s);
			if (z->len == 1) {
				if (z->type == BCF_BT_FLOAT) ksprintf(s, "%g", z->v1.f);
				else if (z->type != BCF_BT_CHAR) kputw(z->v1.i, s);
				else kputc(z->v1.i, s);
			} else bcf_fmt_array(s, z->len, z->type, z->vptr);
		}
	} else kputc('.', s);
	// FORMAT and individual information
	if (v->n_sample && v->n_fmt) { // FORMAT
		int i, j, gt_i = -1;
		bcf_fmt_t *fmt = v->d.fmt;
        int first = 1;
		for (i = 0; i < (int)v->n_fmt; ++i) {
            if ( !fmt[i].p ) continue;
			kputc(!first ? ':' : '\t', s); first = 0;
            if ( fmt[i].id<0 ) //!bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) ) 
            {
                fprintf(stderr, "[E::%s] invalid BCF, the FORMAT tag id=%d not present in the header.\n", __func__, fmt[i].id);
                abort();
            }
			kputs(h->id[BCF_DT_ID][fmt[i].id].key, s);
			if (strcmp(h->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0) gt_i = i;
		}
		for (j = 0; j < v->n_sample; ++j) {
			kputc('\t', s);
            first = 1;
			for (i = 0; i < (int)v->n_fmt; ++i) {
				bcf_fmt_t *f = &fmt[i];
                if ( !f->p ) continue;
				if (!first) kputc(':', s); first = 0;
                if (gt_i == i)
                    bcf_format_gt(f,j,s);
                else
                    bcf_fmt_array(s, f->n, f->type, f->p + j * f->size);
			}
		}
	}
    kputc('\n', s);
	return 0;
}

int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
    int ret;
    fp->line.l = 0;
    vcf_format1(h, v, &fp->line);
    if ( fp->is_compressed==1 )
        ret = bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l);
    else
        ret = hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
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

hts_idx_t *bcf_index(htsFile *fp, int min_shift)
{
	int n_lvls, i;
	bcf1_t *b;
	hts_idx_t *idx;
	bcf_hdr_t *h;
	int64_t max_len = 0, s;
	h = bcf_hdr_read(fp);
    if ( !h ) return NULL;
	for (i = 0; i < h->n[BCF_DT_CTG]; ++i)
		if (max_len < h->id[BCF_DT_CTG][i].val->info[0])
			max_len = h->id[BCF_DT_CTG][i].val->info[0];
    if ( !max_len ) max_len = ((int64_t)1<<31) - 1;  // In case contig line is broken.
	max_len += 256;
	for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);
	idx = hts_idx_init(h->n[BCF_DT_CTG], HTS_FMT_CSI, bgzf_tell(fp->fp.bgzf), min_shift, n_lvls);
	b = bcf_init1();
	while (bcf_read1(fp,h, b) >= 0) {
		int ret;
		ret = hts_idx_push(idx, b->rid, b->pos, b->pos + b->rlen, bgzf_tell(fp->fp.bgzf), 1);
		if (ret < 0) break;
	}
	hts_idx_finish(idx, bgzf_tell(fp->fp.bgzf));
	bcf_destroy1(b);
	bcf_hdr_destroy(h);
	return idx;
}

int bcf_index_build(const char *fn, int min_shift)
{
    htsFile *fp;
	hts_idx_t *idx;
	if ((fp = hts_open(fn, "rb")) == 0) return -1;
    if ( !fp->fp.bgzf->is_compressed ) { hts_close(fp); return -1; }
	idx = bcf_index(fp, min_shift);
	hts_close(fp);
    if ( !idx ) return -1;
	hts_idx_save(idx, fn, HTS_FMT_CSI);
	hts_idx_destroy(idx);
	return 0;
}

/*****************
 *** Utilities ***
 *****************/

bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr)
{
    bcf_hdr_t *hout = bcf_hdr_init("r");
    hout->text = strdup(hdr->text);
    hout->l_text = hdr->l_text;
    bcf_hdr_parse(hout);
    return hout;
}

bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap)
{
	kstring_t str;
	bcf_hdr_t *h;
	str.l = str.m = 0; str.s = 0;
	h = bcf_hdr_init("w");
	if (h0->n[BCF_DT_SAMPLE] > 0) {
		char *p;
		int i = 0, end = n? 8 : 7;
		while ((p = strstr(h0->text, "#CHROM\t")) != 0)
			if (p > h0->text && *(p-1) == '\n') break;
		while ((p = strchr(p, '\t')) != 0 && i < end) ++i, ++p;
		if (i != end) {
			free(h); free(str.s);
			return 0; // malformated header
		}
		kputsn(h0->text, p - h0->text, &str);
		for (i = 0; i < n; ++i) {
			imap[i] = bcf_hdr_id2int(h0, BCF_DT_SAMPLE, samples[i]);
			if (imap[i] < 0) continue;
			kputc('\t', &str);
			kputs(samples[i], &str);
		}
	} else kputsn(h0->text, h0->l_text, &str);
    while (str.l && (!str.s[str.l-1] || str.s[str.l-1]=='\n') ) str.l--; // kill trailing zeros and newlines
    kputc('\n',&str);
	h->text = str.s;
	h->l_text = str.l;
	bcf_hdr_parse(h);
	return h;
}

int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap)
{
	kstring_t ind;
	ind.s = 0; ind.l = ind.m = 0;
	if (n) {
		bcf_fmt_t *fmt;
		int i, j;
		fmt = (bcf_fmt_t*)alloca(v->n_fmt * sizeof(bcf_fmt_t));
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
		if (strlen(v->d.allele[i]) != 1) break;
	return i == v->n_allele;
}

static void bcf_set_variant_type(const char *ref, const char *alt, variant_t *var)
{
	// The most frequent case
	if ( !ref[1] && !alt[1] )
	{
		if ( *alt == '.' || *ref==*alt ) { var->n = 0; var->type = VCF_REF; return; }
        if ( *alt == 'X' ) { var->n = 0; var->type = VCF_REF; return; }  // mpileup's X allele shouldn't be treated as variant
		var->n = 1; var->type = VCF_SNP; return;
	}

	const char *r = ref, *a = alt;
	while (*r && *a && *r==*a ) { r++; a++; }

	if ( *a && !*r )
	{
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
    while ( *re==*ae && re>r && ae>a ) { re--; ae--; }
    if ( ae==a ) 
    { 
        if ( re==r ) { var->n = 1; var->type = VCF_SNP; return; }
        var->n = -(re-r);
        if ( *re==*ae ) { var->type = VCF_INDEL; return; }
        var->type = VCF_OTHER; return;
    }
    else if ( re==r ) 
    { 
        var->n = ae-a;
        if ( *re==*ae ) { var->type = VCF_INDEL; return; }
        var->type = VCF_OTHER; return;
    }

	var->type = ( re-r == ae-a ) ? VCF_MNP : VCF_OTHER;
    var->n = ( re-r > ae-a ) ? -(re-r+1) : ae-a+1;

	// should do also complex events, SVs, etc...
}

static void bcf_set_variant_types(bcf1_t *b)
{
    if ( !(b->unpacked & BCF_UN_STR) ) bcf_unpack(b, BCF_UN_STR);
	bcf_dec_t *d = &b->d;
	if ( d->n_var < b->n_allele ) 
	{
		d->var = (variant_t *) realloc(d->var, sizeof(variant_t)*b->n_allele);
		d->n_var = b->n_allele;
	}
	int i;
	b->d.var_type = 0;
	for (i=1; i<b->n_allele; i++)
	{
		bcf_set_variant_type(d->allele[0],d->allele[i], &d->var[i]);
		b->d.var_type |= d->var[i].type;
		//fprintf(stderr,"[set_variant_type] %d   %s %s -> %d %d .. %d\n", b->pos+1,d->allele[0],d->allele[i],d->var[i].type,d->var[i].n, b->d.var_type);
	}
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

int bcf_update_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
{
    // Is the field already present?
    int i, inf_id = bcf_hdr_id2int(hdr,BCF_DT_ID,key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,inf_id) )
    {
        fprintf(stderr,"[%s:%d %s] The tag was not defined in the header: %s\n", __FILE__,__LINE__,__FUNCTION__,key);
        exit(-1);
    }

    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    for (i=0; i<line->n_info; i++)
        if ( inf_id==line->d.info[i].key ) break;
    bcf_info_t *inf = i==line->n_info ? NULL : &line->d.info[i];

    if ( !n || (type==BCF_HT_STR && !values) )
    {
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
        }
        return 0;
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
    else
    {
        fprintf(stderr, "[E::%s] the type %d not implemented yet\n", __func__, type);
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
            bcf_unpack_info_core1(ptr, inf);
        }
        else
        {
            assert( !inf->vptr_free );  // fix the caller or improve here: this has been modified before 
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
    return 0;
}


int bcf_update_format(bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
{
    // Is the field already present?
    int i, fmt_id = bcf_hdr_id2int(hdr,BCF_DT_ID,key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,fmt_id) )
    {
        if ( !n ) return 0;
        fprintf(stderr,"[%s:%d] Wrong usage of bcf_update_format: The key \"%s\" not present in the header.\n",  __FILE__, __LINE__, key);
        exit(-1);
    }

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==fmt_id ) break;
    bcf_fmt_t *fmt = i==line->n_fmt ? NULL : &line->d.fmt[i];

    if ( !n )
    {
        // Mark the tag for removal, free existing memory if necessary
        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
        return 0;
    }

    int nps = n / line->n_sample;  // number of values per sample
    assert( nps && nps*line->n_sample==n );     // must be divisible by n_sample

    // Encode the values and determine the size required to accommodate the values
    kstring_t str = {0,0,0};
    bcf_enc_int1(&str, fmt_id);
    if ( type==BCF_HT_INT )
        bcf_enc_vint(&str, n, (int32_t*)values, nps);
    else if ( type==BCF_HT_REAL )
    {
        bcf_enc_size(&str, nps, BCF_BT_FLOAT);
        kputsn((char*)values, nps*line->n_sample*sizeof(float), &str);
    }
    else if ( type==BCF_HT_STR )
    {
        bcf_enc_size(&str, nps, BCF_BT_CHAR);
        kputsn((char*)values, nps*line->n_sample, &str);
    }
    else
    {
        fprintf(stderr, "[E::%s] the type %d not implemented yet\n", __func__, type);
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
        // The tag is already present, check if it is big enough to accomodate the new block
        if ( str.l <= fmt->p_len + fmt->p_off )
        {
            // good, the block is big enough
            if ( str.l != fmt->p_len + fmt->p_off ) line->d.indiv_dirty = 1;
            uint8_t *ptr = fmt->p - fmt->p_off;
            memcpy(ptr, str.s, str.l);
            free(str.s);
            bcf_unpack_fmt_core1(ptr, line->n_sample, fmt);
        }
        else
        {
            assert( !fmt->p_free );  // fix the caller or improve here: this has been modified before 
            bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
            fmt->p_free = 1;
            line->d.indiv_dirty = 1;
        }
    }
    line->unpacked |= BCF_UN_FMT;
    return 0;
}


int bcf_update_filter(bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n)
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

int bcf_add_filter(bcf_hdr_t *hdr, bcf1_t *line, int flt_id)
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

static inline int _bcf1_sync_alleles(bcf_hdr_t *hdr, bcf1_t *line, int nals)
{
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
    line->n_allele = nals;
    line->d.shared_dirty |= BCF1_DIRTY_ALS;
    return 0;
}
int bcf_update_alleles(bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
{
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

int bcf_update_alleles_str(bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string)
{
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

int bcf_update_id(bcf_hdr_t *hdr, bcf1_t *line, const char *id)
{
    kstring_t tmp;
    tmp.l = 0; tmp.s = line->d.id; tmp.m = line->d.m_id;
    kputs(id, &tmp);
    line->d.id = tmp.s; line->d.m_id = tmp.m;
    line->d.shared_dirty |= BCF1_DIRTY_ID;
    return 0;
}

bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
{
    int i, id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) ) return NULL;   // no such FMT field in the header
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);
    for (i=0; i<line->n_fmt; i++)  
    {
        if ( line->d.fmt[i].id==id ) return &line->d.fmt[i];
    }
    return NULL;
}

bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
{
    int i, id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,id) ) return NULL;   // no such INFO field in the header
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);
    for (i=0; i<line->n_info; i++)  
    {
        if ( line->d.info[i].key==id ) return &line->d.info[i];
    }
    return NULL;
}

int bcf_get_info_values(bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)
{
    int i,j, tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id) ) return -1;    // no such INFO field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_INFO,tag_id)!=type ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key==tag_id ) break;
    if ( i==line->n_info ) return -3;                               // the tag is not present in this record
    bcf_info_t *info = &line->d.info[i];

    // Make sure the buffer is big enough
    int size1 = type==BCF_HT_INT ? sizeof(int) : sizeof(float);
    if ( *ndst < info->len )
    {
        *ndst = info->len;
        *dst  = realloc(*dst, *ndst * size1);
    }

    if ( info->len == 1 )
    {
        if ( info->type==BCF_BT_FLOAT ) *((float*)*dst) = info->v1.f;
        else *((int*)*dst) = info->v1.i;
        return 1;
    }

    #define BRANCH(type_t, is_missing, is_vector_end, set_missing, out_type_t) { \
        out_type_t *tmp = (out_type_t *) *dst; \
        type_t *p = (type_t *) info->vptr; \
        for (j=0; j<info->len; j++) \
        { \
            if ( is_vector_end ) return j; \
            if ( is_missing ) set_missing; \
            else *tmp = p[j]; \
            tmp++; \
        } \
        return j; \
    }
    switch (info->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  *tmp=bcf_int32_missing, int); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, *tmp=bcf_int32_missing, int); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, *tmp=bcf_int32_missing, int); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), bcf_float_set_missing(*tmp), float); break;
        default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
    }
    #undef BRANCH
    return -4;  // this can never happen
}

int bcf_get_format_values(bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)
{
    int i,j, tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) ) return -1;    // no such FORMAT field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=type ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==tag_id ) break;
    if ( i==line->n_fmt ) return -3;                               // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];

    // Make sure the buffer is big enough
    int nsmpl = hdr->n[BCF_DT_SAMPLE];
    int size1 = type==BCF_HT_INT ? sizeof(int) : sizeof(float);
    if ( *ndst < fmt->n*nsmpl )
    {
        *ndst = fmt->n*nsmpl;
        *dst  = realloc(*dst, *ndst*size1);
    }

    #define BRANCH(type_t, is_missing, is_vector_end, set_missing, set_vector_end, out_type_t) { \
        out_type_t *tmp = (out_type_t *) *dst; \
        type_t *p = (type_t*) fmt->p; \
        for (i=0; i<nsmpl; i++) \
        { \
            for (j=0; j<fmt->n; j++) \
            { \
                if ( is_missing ) set_missing; \
                else if ( is_vector_end ) { set_vector_end; break; } \
                else *tmp = p[j]; \
                tmp++; \
            } \
            for (; j<fmt->n; j++) { set_vector_end; tmp++; } \
            p = (type_t *)((char *)p + fmt->size); \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, int); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, int); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, int); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), bcf_float_set_missing(*tmp), bcf_float_set_vector_end(*tmp), float); break;
        default: fprintf(stderr,"TODO: %s:%d .. fmt->type=%d\n", __FILE__,__LINE__, fmt->type); exit(1);
    }
    #undef BRANCH
    return nsmpl*fmt->n;
}

