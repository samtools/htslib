#include "vcfutils.h"

int calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)
{
	int i;
	for (i=0; i<line->n_allele; i++) ac[i]=0;

	// Use INFO/AC,AN field only when asked
	if ( which&BCF_UN_INFO )
	{
		bcf_unpack(line, BCF_UN_INFO);
		int an_id = bcf_id2int(header, BCF_DT_ID, "AN");
		int ac_id = bcf_id2int(header, BCF_DT_ID, "AC");
		if ( an_id>=0 && ac_id>=0 )
		{
			int i, an=0, ac_len=0, ac_type=0;
			uint8_t *ac_ptr=NULL;
			for (i=0; i<line->n_info; i++)
			{
				bcf_info_t *z = &line->d.info[i];
				if ( z->key == an_id ) an = z->v1.i;
				else if ( z->key == ac_id ) { ac_ptr = z->vptr; ac_len = z->len; ac_type = z->type; }
			}
			int nac = 0;
            #define BRANCH_INT(type_t) {        \
                type_t *p = (type_t *) ac_ptr;  \
                for (i=0; i<ac_len; i++)        \
                {                               \
                    ac[i+1] = p[i];             \
                    nac += p[i];                \
                }                               \
            }
            switch (ac_type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t); break;
                default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, ac_type); exit(1); break;
            }
            #undef BRANCH_INT
			ac[0] = an - nac;
			return 1;
		}
	}

	// Split genotype fields only when asked
	if ( which&BCF_UN_FMT )
	{
		int i, gt_id = bcf_id2int(header,BCF_DT_ID,"GT");
		if ( gt_id<0 ) return 0;
		bcf_unpack(line, BCF_UN_FMT);
		bcf_fmt_t *gt = NULL;
		for (i=0; i<(int)line->n_fmt; i++) 
			if ( line->d.fmt[i].id==gt_id ) { gt = &line->d.fmt[i]; break; }
		if ( !gt ) return 0;
		uint8_t *p = gt->p;
		for (i=0; i<line->n_sample; i++)
		{
			int ial;
			for (ial=0; ial<gt->size; ial++)
			{
                if ( !*p || !(*p)>>1 || *p==(uint8_t)INT8_MIN ) { p += gt->size - ial; break; }
				ac[((*p)>>1)-1]++;
				p++;
			}
		}
		return 1;
	}

	return 0;
}

inline int gt_type(bcf_fmt_t *fmt_ptr, int isample, int *ial)
{
	uint8_t *p = &fmt_ptr->p[isample*fmt_ptr->size];
	int i, a = p[0]>>1, b = a, min = a, nref = a>1 ? a : 255;
	for (i=1; i<fmt_ptr->size; i++)
	{
        if ( p[i] == (uint8_t)INT8_MIN ) break;   // smaller ploidy
		int tmp = p[i]>>1;
		if ( tmp < min ) min = tmp;
		if ( tmp > 1 && nref > tmp ) nref = tmp;
		a |= tmp;
		b &= tmp;
	}
	if ( min==0 ) return GT_UNKN;       // missing GT
	if ( ial ) *ial = nref-1;
	if ( a==b ) return min==1 ? GT_HOM_RR : GT_HOM_AA;
	return min==1 ? GT_HET_RA : GT_HET_AA;
}

bcf_fmt_t *get_fmt_ptr(const bcf_hdr_t *header, bcf1_t *line, char *tag)
{
    bcf_unpack(line, BCF_UN_FMT);

    int i, id = bcf_id2int(header, BCF_DT_ID, tag);
    if ( id<0 ) return NULL;

    bcf_unpack(line, BCF_UN_FMT);
    for (i=0; i<(int)line->n_fmt; i++)
        if ( line->d.fmt[i].id==id ) return &line->d.fmt[i];

    return NULL;
}

int trim_alleles(const bcf_hdr_t *header, bcf1_t *line)
{
    int i;
    bcf_fmt_t *gt = get_fmt_ptr(header, line, "GT");
    if ( !gt ) return 0;

    int *ac = (int*) calloc(line->n_allele,sizeof(int));

    // check if all alleles are populated
    uint8_t *p = gt->p;
    for (i=0; i<line->n_sample; i++)
    {
        int ial;
        for (ial=0; ial<gt->size; ial++)
        {
            if ( !*p || !(*p)>>1 || *p==(uint8_t)INT8_MIN ) { p += gt->size - ial; break; }
            ac[((*p)>>1)-1]++;
            p++;
        }
    }

    int rm_als = 0, nrm = 0;
    for (i=1; i<line->n_allele; i++) 
    {
        if ( !ac[i] ) { rm_als |= 1<<i; nrm++; }
    }
    free(ac);

    if ( nrm ) remove_alleles(header, line, rm_als);
    return nrm;
}

extern uint32_t bcf_missing_float;

void remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int rm_mask)
{
    int *map = (int*) calloc(line->n_allele, sizeof(int));

    // create map of indexes from old to new ALT numbering and modify ALT
    int nrm = 0, i, j;
    for (i=1, j=1; i<line->n_allele; i++) 
    {
        if ( rm_mask & 1<<i )
        {
            // remove this allele
            line->d.allele[i] = NULL;
            nrm++;
            continue;
        }
        if ( line->d.allele[j] ) 
        {
            map[i] = j++;
            continue;
        }
        line->d.allele[j] = line->d.allele[i];
        line->d.allele[i] = NULL;
        map[i] = j;
        j++;
    }
    if ( !nrm ) { free(map); return; }

    // remove from GT fields 
    bcf_fmt_t *gt = get_fmt_ptr(header, line, "GT");
    if ( gt )
    {
        for (i=1; i<line->n_allele; i++) if ( map[i]!=i ) break;
        if ( i<line->n_allele )
        {
            uint8_t *p = gt->p;
            for (i=0; i<line->n_sample; i++)
            {
                int ial;
                for (ial=0; ial<gt->size; ial++)
                {
                    if ( !*p || !(*p)>>1 || *p==(uint8_t)INT8_MIN ) { p += gt->size - ial; break; }
                    *p = (map[((*p)>>1)-1] + 1) <<1 | ((*p)&1);
                    p++;
                }
            }
        }
    }

    // remove from Number=G fields. Assuming haploid or diploid GTs
    int nG_ori = line->n_allele*(line->n_allele + 1)/2;
    int nG_new = (line->n_allele - nrm)*(line->n_allele - nrm + 1)/2;
    for (i=0; i<(int)line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];

        if ( ((header->id[BCF_DT_ID][fmt->id].val->info[BCF_HL_FMT]>>8)&0xf) == BCF_VL_G )
        {
            assert( fmt->n==nG_ori );

            #define BRANCH_INT(type_t,missing) { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *p = (type_t *) (fmt->p + j*fmt->size); \
                    int k, nset = 0; \
                    for (k=0; k<nG_ori; k++) if ( p[k] != missing ) nset++; \
                    if ( nset==nG_ori ) \
                    { \
                        /* diploid */ \
                        int ia, ib, k_ori = 0, k_new = 0; \
                        for (ia=0; ia<line->n_allele; ia++) \
                        { \
                            for (ib=0; ib<=ia; ib++) \
                            { \
                                if ( rm_mask & 1<<ia || rm_mask & 1<<ib ) { k_ori++; continue; } \
                                p[k_new] = p[k_ori]; \
                                k_ori++; \
                                k_new++; \
                            } \
                        } \
                    } \
                    else if ( nset==line->n_allele ) \
                    { \
                        /* haploid */ \
                        int k_ori, k_new = 0; \
                        for (k_ori=0; k_ori<line->n_allele; k_ori++) \
                            if ( !(rm_mask & 1<<k_ori) ) p[k_new++] = p[k_ori]; \
                        for (; k_new<line->n_allele; k_new++) \
                            p[k_new] = missing; \
                    } \
                    else { fprintf(stderr, "[E::%s] todo, missing values: %d %d\n", __func__, nset,nG_ori); exit(1); } \
                } \
            }
            switch (fmt->type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t,INT8_MIN); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t,INT16_MIN); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t,INT32_MIN); break;
                case BCF_BT_FLOAT: BRANCH_INT(float,bcf_missing_float); break;  // fixme: float will not work
                default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt->type); exit(1); break;
            }
            #undef BRANCH_INT

            fmt->n = nG_new;
        }
        else if ( ((header->id[BCF_DT_ID][fmt->id].val->info[BCF_HL_FMT]>>8)&0xf) == BCF_VL_A )
        {
            fprintf(stderr, "[E::%s] todo A\n", __func__); exit(1);
        }
    }
    line->n_allele -= nrm;

    free(map);
    return;
}



