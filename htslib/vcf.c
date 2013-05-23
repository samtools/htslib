#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "kstring.h"
#include "bgzf.h"
#include "vcf.h"
#include "tbx.h"

#include "khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

uint32_t bcf_float_missing    = 0x7F800001;
uint32_t bcf_float_vector_end = 0x7F800002;
uint8_t bcf_type_shift[] = { 0, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static bcf_idinfo_t bcf_idinfo_def = { .info = { 15, 15, 15 }, .hrec = { NULL, NULL, NULL}, .id = -1 };

/*************************
 *** VCF header parser ***
 *************************/

void bcf_hdr_add_sample(bcf_hdr_t *h, char *s)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    int ret;
    int k = kh_put(vdict, d, s, &ret);
    if (ret) { // absent
        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = kh_size(d) - 1;
    } else {
        if (hts_verbose >= 2)
            fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
    }
    int n = kh_size(d);
    h->samples = (char**) realloc(h->samples,sizeof(char*)*n);
    h->samples[n-1] = s;
}

void bcf_hdr_parse_sample_line(bcf_hdr_t *h, const char *str)
{
    int i = 0;
    const char *p, *q;
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_ID];
    // add samples
    d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    for (p = q = str;; ++q) {
        if (*q != '\t' && *q != 0 && *q != '\n') continue;
        if (++i > 9) {
            char *s;
            s = (char*)malloc(q - p + 1);
            strncpy(s, p, q - p);
            s[q - p] = 0;
            bcf_hdr_add_sample(h,s);
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

inline int is_escaped(const char *min, const char *str)
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
    while ( *q && *q!='\n' )
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
            else { if ( *q==',' || *q=='>' ) break; }
            q++;
        }
        bcf_hrec_set_val(hrec, hrec->nkeys-1, p, q-p, quoted);
        if ( quoted ) q++;
        if ( *q=='>' ) { q++; break; }
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
        int id = bcf_id2int(hdr, BCF_DT_ID, "PL");
        if ( id>=0 && bcf_idinfo_exists(hdr,BCF_HL_FMT,id) && bcf_id2length(hdr,BCF_HL_FMT,id)!=BCF_VL_G )
        {
            fprintf(stderr,"[W::%s] PL should be declared as Number=G\n", __func__);
            PL_warned = 1;
        }
    }
    if ( !GL_warned )
    {
        int id = bcf_id2int(hdr, BCF_HL_FMT, "GL");
        if ( id>=0 && bcf_idinfo_exists(hdr,BCF_HL_FMT,id) && bcf_id2length(hdr,BCF_HL_FMT,id)!=BCF_VL_G )
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


/**********************
 *** BCF header I/O ***
 **********************/

bcf_hdr_t *bcf_hdr_init(void)
{
	int i;
	bcf_hdr_t *h;
	h = (bcf_hdr_t*)calloc(1, sizeof(bcf_hdr_t));
	for (i = 0; i < 3; ++i)
		h->dict[i] = kh_init(vdict);
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

bcf_hdr_t *bcf_hdr_read(BGZF *fp)
{
	uint8_t magic[5];
	bcf_hdr_t *h;
	h = bcf_hdr_init();
	bgzf_read(fp, magic, 5);
	if (strncmp((char*)magic, "BCF\2\2", 5) != 0) {
		if (hts_verbose >= 2)
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

void bcf_hdr_write(BGZF *fp, const bcf_hdr_t *h)
{
	bgzf_write(fp, "BCF\2\2", 5);
	bgzf_write(fp, &h->l_text, 4);
	bgzf_write(fp, h->text, h->l_text);
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
}

void bcf_destroy1(bcf1_t *v)
{
    bcf_clear1(v);
	free(v->d.id); free(v->d.als); free(v->d.allele); free(v->d.flt); free(v->d.info); free(v->d.fmt);
	if (v->d.var ) free(v->d.var);
	free(v->shared.s); free(v->indiv.s);
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

int bcf_read1(BGZF *fp, bcf1_t *v) { return bcf_read1_core(fp, v); }

int bcf_readrec(BGZF *fp, void *null, bcf1_t *v, int *tid, int *beg, int *end)
{
	int ret;
	if ((ret = bcf_read1_core(fp, v)) >= 0)
		*tid = v->rid, *beg = v->pos, *end = v->pos + v->rlen;
	return ret;
}

int bcf_write1(BGZF *fp, const bcf1_t *v)
{
	uint32_t x[8];
	x[0] = v->shared.l + 24; // to include six 32-bit integers
	x[1] = v->indiv.l;
	memcpy(x + 2, v, 16);
	x[6] = (uint32_t)v->n_allele<<16 | v->n_info;
	x[7] = (uint32_t)v->n_fmt<<24 | v->n_sample;
	bgzf_write(fp, x, 32);
	bgzf_write(fp, v->shared.s, v->shared.l);
	bgzf_write(fp, v->indiv.s, v->indiv.l);
	return 0;
}

/**********************
 *** VCF header I/O ***
 **********************/

bcf_hdr_t *vcf_hdr_read(htsFile *fp)
{
	if (fp->is_bin) 
        return bcf_hdr_read((BGZF*)fp->fp);

    kstring_t txt, *s = &fp->line;
    bcf_hdr_t *h;
    h = bcf_hdr_init();
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
        if (s->s[1] != '#') break;
        kputc('\n', &txt);
    }
    h->l_text = txt.l + 1; // including NULL
    h->text = txt.s;
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
    hdr->l_text = txt.l;
}

const char **bcf_seqnames(const bcf_hdr_t *h, int *n)
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

void vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h)
{
	if (!fp->is_bin) {
		int l = h->l_text;
		while (l && h->text[l-1] == 0) --l; // kill the trailing zeros
		if (l && h->text[l-1] == '\n') --l;
		fwrite(h->text, 1, l, (FILE*)fp->fp);
		fputc('\n', (FILE*)fp->fp);
	} else bcf_hdr_write((BGZF*)fp->fp, h);
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

int *bcf_set_iarray(bcf_fmt_t *fmt, int nsmpl, int *arr, int *narr)
{
    if ( nsmpl*fmt->n > *narr )
    {
        *narr = nsmpl*fmt->n;
        arr = (int*) realloc(arr, sizeof(int)*(*narr));
    }

    #define BRANCH(type_t, missing, vector_end) { \
        type_t *ptr = (type_t*) fmt->p; \
        int i, j; \
        if ( fmt->n==1 ) \
        { \
            for (i=0; i<nsmpl; i++) \
            { \
                if ( ptr[0]==missing ) arr[i] = bcf_int32_missing; \
                else if ( ptr[0]==vector_end ) arr[i] = bcf_int32_vector_end; \
                else arr[i] = ptr[0]; \
                ptr = (type_t *)((void*)ptr + fmt->size); \
            } \
        } \
        else \
        { \
            int *p_arr = arr; \
            for (i=0; i<nsmpl; i++) \
            { \
                for (j=0; j<fmt->n; j++) \
                { \
                    if ( ptr[j]==missing ) p_arr[j] = bcf_int32_missing; \
                    else if ( ptr[j]==vector_end ) { p_arr[j] = bcf_int32_vector_end; break; } \
                    else p_arr[j] = ptr[j]; \
                } \
                for (; j<fmt->n; j++) p_arr[j] = bcf_int32_vector_end; \
                ptr = (type_t *)((void*)ptr + fmt->size); \
                p_arr += fmt->n; \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr,"fixme: type %d in bcf_set_iarray?\n", fmt->type); abort(); break;
    }
    #undef BRANCH
    return arr;
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
        for (j = 0; j < n && *p; ++j, ++p) kputc(*p, s);
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

int vcf_parse1(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
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
                }
                fmt[j].max_l = fmt[j].max_m = fmt[j].max_g = 0;
                fmt[j].key = kh_val(d, k).id;
                fmt[j].is_gt = !strcmp(t, "GT");
                fmt[j].y = h->id[0][fmt[j].key].val->info[BCF_HL_FMT];
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

int vcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
{
	if (!fp->is_bin) {
		int ret;
		ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
		if (ret < 0) return -1;
		ret = vcf_parse1(&fp->line, h, v);
		return 0;
	} else return bcf_read1((BGZF*)fp->fp, v);
}

inline uint8_t *bcf_unpack_fmt_core1(uint8_t *ptr, int n_sample, bcf_fmt_t *fmt)
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

inline uint8_t *bcf_unpack_info_core1(uint8_t *ptr, bcf_info_t *info)
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
	uint8_t *ptr = (uint8_t*)b->shared.s;
	int *offset, i;
	bcf_dec_t *d = &b->d;
	if (which & BCF_UN_FLT) which |= BCF_UN_STR;
	if (which & BCF_UN_INFO) which |= BCF_UN_SHR;
	if ((which&BCF_UN_STR) && !(b->unpacked&BCF_UN_STR)) { // ID
        kstring_t tmp;
        tmp.l = 0; tmp.m = d->m_str; tmp.s = d->id;
        ptr = bcf_fmt_sized_array(&tmp, ptr); kputc('\0', &tmp);
        d->m_str = tmp.m; d->id = tmp.s; // write tmp back
        // REF and ALT
        tmp.l = 0; tmp.m = d->m_als; tmp.s = d->als;
        offset = (int*)alloca(b->n_allele * sizeof(int));
        for (i = 0; i < b->n_allele; ++i) {
            offset[i] = tmp.l;
            ptr = bcf_fmt_sized_array(&tmp, ptr);
            kputc('\0', &tmp);
        }
        hts_expand(char*, b->n_allele, d->m_allele, d->allele); // NM: hts_expand() is a macro
        for (i = 0; i < b->n_allele; ++i)
            d->allele[i] = tmp.s + offset[i];
        d->m_als = tmp.m; d->als = tmp.s; // write tmp back
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

int vcf_format1(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)
{
	uint8_t *ptr = (uint8_t*)v->shared.s;
	int i;
	bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
	kputs(h->id[BCF_DT_CTG][v->rid].key, s); // CHROM
	kputc('\t', s); kputw(v->pos + 1, s); // POS
	kputc('\t', s); kputs(v->d.id, s); // ID
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
	ptr = (uint8_t*)v->indiv.s;
	if (v->n_sample && v->n_fmt) { // FORMAT
		int i, j, gt_i = -1;
		bcf_fmt_t *fmt = v->d.fmt;
        // This way it is not possible to modify the output line. Also,
        // the bcf_unpack_fmt_core has been already called above.
		//      fmt = (bcf_fmt_t*)alloca(v->n_fmt * sizeof(bcf_fmt_t));
		//      ptr = bcf_unpack_fmt_core(ptr, v->n_sample, v->n_fmt, fmt);
        int first = 1;
		for (i = 0; i < (int)v->n_fmt; ++i) {
            if ( !fmt[i].p ) continue;
			kputc(!first ? ':' : '\t', s); first = 0;
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
	return 0;
}

int vcf_write1(htsFile *fp, const bcf_hdr_t *h, const bcf1_t *v)
{
	if (!fp->is_bin) {
	    fp->line.l = 0;
		vcf_format1(h, v, &fp->line);
		fwrite(fp->line.s, 1, fp->line.l, (FILE*)fp->fp);
		fputc('\n', (FILE*)fp->fp);
	} else return bcf_write1((BGZF*)fp->fp, v);
	return 0;
}

/************************
 * Data access routines *
 ************************/

int bcf_id2int(const bcf_hdr_t *h, int which, const char *id)
{
	khint_t k;
	vdict_t *d = (vdict_t*)h->dict[which];
	k = kh_get(vdict, d, id);
	return k == kh_end(d)? -1 : kh_val(d, k).id;
}

int bcf_name2id(const bcf_hdr_t *h, const char *id)
{
	return bcf_id2int(h, BCF_DT_CTG, id);
}

/********************
 *** BCF indexing ***
 ********************/

hts_idx_t *bcf_index(BGZF *fp, int min_shift)
{
	int n_lvls, i;
	bcf1_t *b;
	hts_idx_t *idx;
	bcf_hdr_t *h;
	int64_t max_len = 0, s;
	h = bcf_hdr_read(fp);
	for (i = 0; i < h->n[BCF_DT_CTG]; ++i)
		if (max_len < h->id[BCF_DT_CTG][i].val->info[0])
			max_len = h->id[BCF_DT_CTG][i].val->info[0];
    if ( !max_len ) max_len = ((int64_t)1<<31) - 1;  // In case contig line is broken.
	max_len += 256;
	for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);
	idx = hts_idx_init(h->n[BCF_DT_CTG], HTS_FMT_CSI, bgzf_tell(fp), min_shift, n_lvls);
	bcf_hdr_destroy(h);
	b = bcf_init1();
	while (bcf_read1(fp, b) >= 0) {
		int ret;
		ret = hts_idx_push(idx, b->rid, b->pos, b->pos + b->rlen, bgzf_tell(fp), 1);
		if (ret < 0) break;
	}
	hts_idx_finish(idx, bgzf_tell(fp));
	bcf_destroy1(b);
	return idx;
}

int bcf_index_build(const char *fn, int min_shift)
{
	BGZF *fp;
	hts_idx_t *idx;
	if ((fp = bgzf_open(fn, "r")) == 0) return -1;
	idx = bcf_index(fp, min_shift);
	bgzf_close(fp);
	hts_idx_save(idx, fn, HTS_FMT_CSI);
	hts_idx_destroy(idx);
	return 0;
}

/*****************
 *** Utilities ***
 *****************/

bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap)
{
	kstring_t str;
	bcf_hdr_t *h;
	str.l = str.m = 0; str.s = 0;
	h = bcf_hdr_init();
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
			imap[i] = bcf_id2int(h0, BCF_DT_SAMPLE, samples[i]);
			if (imap[i] < 0) continue;
			kputc('\t', &str);
			kputs(samples[i], &str);
		}
	} else kputsn(h0->text, h0->l_text, &str);
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
    v->unpacked = 0;    // only BCF is ready for output, VCF will need to unpack again
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

void bcf_set_variant_types(bcf1_t *b)
{
	if ( b->d.var_type!=-1 ) return;	// already set

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

int bcf1_update_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key, void *values, int n, int type)
{
    // Get the number and type from the bcf_idinfo_t header line
    vdict_t *d = (vdict_t*) hdr->dict[BCF_DT_ID];
    khint_t k = kh_get(vdict, d, key);
    int inf_id = kh_val(d,k).id;

    // Is the field already present?
    int i;
    for (i=0; i<line->n_info; i++)
        if ( inf_id==line->d.info[i].key ) break;
    bcf_info_t *inf = i==line->n_info ? NULL : &line->d.info[i];

    if ( !n )
    {
        if ( inf )
        {
            // Mark the tag for removal
            line->d.shared_dirty = 1;
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
            if ( str.l != inf->vptr_len + inf->vptr_off ) line->d.shared_dirty = 1;
            uint8_t *ptr = inf->vptr - inf->vptr_off;
            memcpy(ptr, str.s, str.l);
            free(str.s);
            bcf_unpack_info_core1(ptr, inf);
        }
        else
        {
            bcf_unpack_info_core1((uint8_t*)str.s, inf);
            inf->vptr_free = 1;
            line->d.shared_dirty = 1;
        }
    }
    else
    {
        // The tag is not present, create new one
        line->n_info++;
        hts_expand(bcf_info_t, line->n_info, line->d.m_info , line->d.info);
        inf = &line->d.info[line->n_info-1];
        bcf_unpack_info_core1((uint8_t*)str.s, inf);
        inf->vptr_free = 1;
        line->d.shared_dirty = 1;
    }
    return 0;
}


int bcf1_update_format(bcf_hdr_t *hdr, bcf1_t *line, const char *key, void *values, int n, int type)
{
    assert( n );        // todo: n=0

    // Is the field already present?
    int i, fmt_id = bcf_id2int(hdr,BCF_DT_ID,key);
    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==fmt_id ) break;
    bcf_fmt_t *fmt = i==line->n_fmt ? NULL : &line->d.fmt[i];

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
    else
    {
        fprintf(stderr, "[E::%s] the type %d not implemented yet\n", __func__, type);
        abort();
    }

    if ( !fmt )
    {
        line->n_fmt++;
        hts_expand0(bcf_fmt_t, line->n_fmt, line->d.m_fmt, line->d.fmt);
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
            bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
            fmt->p_free = 1;
            line->d.indiv_dirty = 1;
        }
    }
    return 0;
}



