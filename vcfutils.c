#include "htslib/vcfutils.h"

int bcf_calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)
{
	int i;
	for (i=0; i<line->n_allele; i++) ac[i]=0;

	// Use INFO/AC,AN field only when asked
	if ( which&BCF_UN_INFO )
	{
		bcf_unpack(line, BCF_UN_INFO);
		int an_id = bcf_hdr_id2int(header, BCF_DT_ID, "AN");
		int ac_id = bcf_hdr_id2int(header, BCF_DT_ID, "AC");
        int i, an=0, ac_len=0, ac_type=0;
        uint8_t *ac_ptr=NULL;
		if ( an_id>=0 && ac_id>=0 )
		{
			for (i=0; i<line->n_info; i++)
			{
				bcf_info_t *z = &line->d.info[i];
				if ( z->key == an_id ) an = z->v1.i;
				else if ( z->key == ac_id ) { ac_ptr = z->vptr; ac_len = z->len; ac_type = z->type; }
			}
        }
        if ( ac_ptr )
        {
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
                default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, ac_type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
            }
            #undef BRANCH_INT
            assert( an>=nac );  // sanity check for missing values
			ac[0] = an - nac;
			return 1;
        }
	}

	// Split genotype fields only when asked
	if ( which&BCF_UN_FMT )
	{
		int i, gt_id = bcf_hdr_id2int(header,BCF_DT_ID,"GT");
		if ( gt_id<0 ) return 0;
		bcf_unpack(line, BCF_UN_FMT);
		bcf_fmt_t *fmt_gt = NULL;
		for (i=0; i<(int)line->n_fmt; i++) 
			if ( line->d.fmt[i].id==gt_id ) { fmt_gt = &line->d.fmt[i]; break; }
		if ( !fmt_gt ) return 0;
        #define BRANCH_INT(type_t,missing,vector_end) { \
		    for (i=0; i<line->n_sample; i++) \
		    { \
                type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
		    	int ial; \
		    	for (ial=0; ial<fmt_gt->n; ial++) \
		    	{ \
                    if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                    if ( !(p[ial]>>1) || p[ial]==missing ) continue; /* missing allele */ \
		    		ac[(p[ial]>>1)-1]++; \
		    	} \
		    } \
        }
        switch (fmt_gt->type) {
            case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, fmt_gt->type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
        }
        #undef BRANCH_INT
		return 1;
	}
	return 0;
}

int bcf_gt_type(bcf_fmt_t *fmt_ptr, int isample, int *_ial, int *_jal)
{
    int i, nals = 0, has_ref = 0, has_alt = 0, ial = 0, jal = 0;
    #define BRANCH_INT(type_t,missing,vector_end) { \
        type_t *p = (type_t*) (fmt_ptr->p + isample*fmt_ptr->size); \
        for (i=0; i<fmt_ptr->n; i++) \
        { \
            if ( p[i] == vector_end ) break; /* smaller ploidy */ \
            if ( !p[i] || p[i] == missing ) continue; /* missing allele */ \
            int tmp = p[i]>>1; \
            if ( tmp>1 ) \
            { \
                if ( !ial ) { ial = tmp; has_alt = 1; } \
                else if ( tmp!=ial ) \
                { \
                    if ( tmp<ial ) \
                    { \
                        jal = ial; \
                        ial = tmp; \
                    } \
                    else \
                    { \
                        jal = tmp; \
                    } \
                    has_alt = 2; \
                } \
            } \
            else has_ref = 1; \
            nals++; \
        } \
    }
    switch (fmt_ptr->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr, "[E::%s] todo: fmt_type %d\n", __func__, fmt_ptr->type); exit(1); break;
    }
    #undef BRANCH_INT

    if ( _ial ) *_ial = ial>0 ? ial-1 : ial;
    if ( _jal ) *_jal = jal>0 ? jal-1 : jal;
    if ( !nals ) return GT_UNKN;
    if ( nals==1 )
        return has_ref ? GT_HAPL_R : GT_HAPL_A;
    if ( !has_ref ) 
        return has_alt==1 ? GT_HOM_AA : GT_HET_AA;
    if ( !has_alt )
        return GT_HOM_RR;
    return GT_HET_RA;
}

int bcf_trim_alleles(const bcf_hdr_t *header, bcf1_t *line)
{
    int i;
    bcf_fmt_t *gt = bcf_get_fmt(header, line, "GT");
    if ( !gt ) return 0;

    int *ac = (int*) calloc(line->n_allele,sizeof(int));

    // check if all alleles are populated
    #define BRANCH(type_t,missing,vector_end) { \
        for (i=0; i<line->n_sample; i++) \
        { \
            type_t *p = (type_t*) (gt->p + i*gt->size); \
            int ial; \
            for (ial=0; ial<gt->size; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( !(p[ial]>>1) || p[ial]==missing ) continue; /* missing allele */ \
                ac[(p[ial]>>1)-1]++; \
            } \
        } \
    }
    switch (gt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, gt->type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
    }
    #undef BRANCH

    int rm_als = 0, nrm = 0;
    for (i=1; i<line->n_allele; i++) 
    {
        if ( !ac[i] ) { rm_als |= 1<<i; nrm++; }
    }
    free(ac);

    if ( nrm ) bcf_remove_alleles(header, line, rm_als);
    return nrm;
}

void bcf_remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int rm_mask)
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
    bcf_fmt_t *gt = bcf_get_fmt(header, line, "GT");
    if ( gt )
    {
        for (i=1; i<line->n_allele; i++) if ( map[i]!=i ) break;
        if ( i<line->n_allele )
        {
            #define BRANCH(type_t,missing,vector_end) { \
                for (i=0; i<line->n_sample; i++) \
                { \
                    type_t *p = (type_t*) (gt->p + i*gt->size); \
                    int ial; \
                    for (ial=0; ial<gt->n; ial++) \
                    { \
                        if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                        if ( !(p[ial]>>1) || p[ial]==missing ) continue; /* missing allele */ \
                        p[ial] = (map[(p[ial]>>1)-1] + 1) <<1 | (p[ial]&1); \
                    } \
                } \
            }
            switch (gt->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
                default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, gt->type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
            }
            #undef BRANCH
        }
    }

    int nG_ori = line->n_allele*(line->n_allele + 1)/2;
    int nG_new = (line->n_allele - nrm)*(line->n_allele - nrm + 1)/2;
    int nA_ori = line->n_allele - 1;
    int nA_new = line->n_allele-nrm - 1;

    // remove from Number=G and Number=A FORMAT fields. Assuming haploid or diploid GTs
    for (i=0; i<(int)line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];

        if ( ((header->id[BCF_DT_ID][fmt->id].val->info[BCF_HL_FMT]>>8)&0xf) == BCF_VL_G ) // FORMAT/Number==G
        {
            if (fmt->type==BCF_BT_CHAR) { continue; } // todo: support for strings
            assert( fmt->n==nG_ori || fmt->n==line->n_allele ); // diploid or haploid

            #define BRANCH(type_t,is_missing,is_vector_end,set_vector_end) { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *p = (type_t *) (fmt->p + j*fmt->size); \
                    int k, nset = 0; \
                    for (k=0; k<nG_ori; k++) \
                        if ( is_vector_end ) break; \
                        else if ( !(is_missing) ) nset++; \
                    if ( nset==0 ) set_vector_end; \
                    else if ( nset==nG_ori ) \
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
                        int k_ori; k = 0; \
                        for (k_ori=0; k_ori<line->n_allele; k_ori++) \
                            if ( !(rm_mask & 1<<k_ori) ) p[k++] = p[k_ori]; \
                        for (; k<line->n_allele; k++) set_vector_end;  \
                    } \
                    else { fprintf(stderr, "[E::%s] todo, missing values: %d %d\n", __func__, nset,nG_ori); exit(1); } \
                } \
            }
            switch (fmt->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  p[k]==bcf_int8_missing,  p[k]==bcf_int8_vector_end,  p[k]=bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH(int16_t, p[k]==bcf_int16_missing, p[k]==bcf_int16_vector_end, p[k]=bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH(int32_t, p[k]==bcf_int32_missing, p[k]==bcf_int32_vector_end, p[k]=bcf_int32_vector_end); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[k]), bcf_float_is_vector_end(p[k]), bcf_float_set_vector_end(p[k]) ); break;
                default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt->type); exit(1); break;
            }
            #undef BRANCH

            fmt->n = nG_new;
        }
        else if ( ((header->id[BCF_DT_ID][fmt->id].val->info[BCF_HL_FMT]>>8)&0xf) == BCF_VL_A )  // FORMAT/Number==A
        {
            if (fmt->type==BCF_BT_CHAR) { continue; } // todo: support for strings
            assert( fmt->n==nA_ori );
            
            if (nA_new == 0)
            {
                // delete FORMAT field since it now has no length
                int k;
                for (k=i; k<(int)line->n_fmt; k++)
                {
                    line->d.fmt[k]=line->d.fmt[k+1];
                }
                line->n_fmt--; i--;
                continue;
            }
            else
            {
                #define BRANCH_INT(type_t,missing) {         \
                    type_t *p = (type_t *) fmt->p;  \
                    int ia, j, k; \
                    for (j=0; j<line->n_sample; j++) \
                    { \
                        for (ia=0, k=0; ia<line->n_allele; ia++) \
                        {                                \
                            if ( rm_mask & 1<<(ia+1) ) \
                            { \
                                p[ia] = missing; \
                                continue; \
                            } \
                            if ( p[k] != missing ) \
                            { \
                                k++; \
                                continue; \
                            } \
                            p[k] = p[ia]; \
                            p[ia] = missing; \
                            k++; \
                        } \
                        p += fmt->n; \
                    } \
                }
                switch (fmt->type) {
                    case BCF_BT_INT8:  BRANCH_INT(int8_t,INT8_MIN); break;
                    case BCF_BT_INT16: BRANCH_INT(int16_t,INT16_MIN); break;
                    case BCF_BT_INT32: BRANCH_INT(int32_t,INT32_MIN); break;
                    case BCF_BT_FLOAT: BRANCH_INT(float,bcf_float_missing); break;
                    default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt->type); exit(1); break;
                }
                #undef BRANCH_INT
            }
            fmt->n = nA_new;
        }
    }
    // remove from Number=G and Number=A INFO fields.
    for (i=0; i<(int)line->n_info; i++)
    {
        bcf_info_t *info = &line->d.info[i];
        
        if ( ((header->id[BCF_DT_ID][info->key].val->info[BCF_HL_INFO]>>8)&0xf) == BCF_VL_G ) // INFO/Number==G
        {
            if (info->type==BCF_BT_CHAR) { continue; } // todo: support for strings
            assert( info->len==nG_ori );

            #define BRANCH(type_t,is_missing,is_vector_end,set_vector_end) { \
                type_t *p = (type_t *) (type_t *) info->vptr; \
                int k, nset = 0; \
                for (k=0; k<nG_ori; k++) \
                    if ( is_vector_end ) break; \
                    else if ( !(is_missing) ) nset++; \
                if ( nset==0 ) return; \
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
                    for (; k_new<line->n_allele; k_new++) set_vector_end; \
                } \
                else { fprintf(stderr, "[E::%s] todo, missing values: %d %d\n", __func__, nset, nG_ori); exit(1); } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  p[k]==bcf_int8_missing,  p[k]==bcf_int8_vector_end,  p[k]=bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH(int16_t, p[k]==bcf_int16_missing, p[k]==bcf_int16_vector_end, p[k]=bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH(int32_t, p[k]==bcf_int32_missing, p[k]==bcf_int32_vector_end, p[k]=bcf_int32_vector_end); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[k]), bcf_float_is_vector_end(p[k]), bcf_float_set_vector_end(p[k]) ); break;
                default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, info->type); exit(1); break;
            }
            #undef BRANCH
            
            info->len = nG_new;
        }
        else if ( ((header->id[BCF_DT_ID][info->key].val->info[BCF_HL_INFO]>>8)&0xf) == BCF_VL_A) // INFO/Number==A
        {
            if (info->type==BCF_BT_CHAR) { continue; } // todo: support for strings
            assert( info->len==nA_ori );
            
            if (nA_new == 0)
            {
                // delete INFO field since it now has no length
                int k;
                for (k=i; k<(int)line->n_info; k++)
                {
                    line->d.info[k]=line->d.info[k+1];
                }
                line->n_info--; i--;
                continue;
            }
            else
            {
                #define BRANCH(type_t,missing,vector_end) {         \
                    type_t *p = (type_t *) info->vptr;  \
                    int ia, j; \
                    for (ia=0, j=0; ia<line->n_allele; ia++) \
                    {                                \
                        if ( rm_mask & 1<<(ia+1) ) \
                        { \
                            p[ia] = missing; \
                            continue; \
                        } \
                        if ( p[j] != missing ) \
                        { \
                            j++; \
                            continue; \
                        } \
                        p[j] = p[ia]; \
                        p[ia] = missing; \
                        j++; \
                    } \
                    if (nA_new == 1) { \
                        info->v1.i = p[0]; \
                    } \
                }
                switch (info->type) {
                    case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
                    case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
                    case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
                    case BCF_BT_FLOAT: BRANCH(float, bcf_float_missing, bcf_float_vector_end); break;
                    default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, info->type); exit(1); break;
                }
                #undef BRANCH
            }
            info->len = nA_new;
        }
    }
    line->n_allele -= nrm;

    free(map);
    return;
}
