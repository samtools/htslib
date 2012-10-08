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
            if ( ac_type==BCF_BT_INT8 ) { BRANCH_INT(uint8_t) }
            else if ( ac_type==BCF_BT_INT16 ) { BRANCH_INT(uint16_t) }
            else if ( ac_type==BCF_BT_INT32 ) { BRANCH_INT(uint32_t) }
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

int trim_alleles(const bcf_hdr_t *header, bcf1_t *line)
{
    int i, gt_id = bcf_id2int(header,BCF_DT_ID,"GT");
    if ( gt_id<0 ) return 0;
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *gt = NULL;
    for (i=0; i<(int)line->n_fmt; i++) 
        if ( line->d.fmt[i].id==gt_id ) { gt = &line->d.fmt[i]; break; }
    if ( !gt ) return 0;

    int *ac = (int*) calloc(line->n_sample,sizeof(int));

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

    // fprintf(stderr,"before: ");
    // for (i=0; i<line->n_allele; i++) fprintf(stderr,"\t(%d) %s", i,line->d.allele[i]);
    // fprintf(stderr,"\n");

    int nrm = 0, j;
    for (i=1, j=1; i<line->n_allele; i++) 
    {
        if ( !ac[i] )
        {
            // remove this allele
            line->d.allele[i] = NULL;
            nrm++;
        }
        else
        {
            if ( line->d.allele[j] ) 
            {
                ac[i] = j++;  // ac now serves as a map from old to new ALT indexes
                continue;
            }
            line->d.allele[j] = line->d.allele[i];
            line->d.allele[i] = NULL;
            while ( j<i && line->d.allele[j] ) j++;
        }
    }

    line->n_allele -= nrm;

    // fprintf(stderr,"after:  ");
    // for (i=0; i<line->n_allele; i++) fprintf(stderr,"\t(%d) %s", i,line->d.allele[i]);
    // fprintf(stderr,"\tnrm=%d\n", nrm);

    if ( nrm )
    {
        for (i=1; i<line->n_allele; i++) 
        {
            if ( ac[i]==i ) continue;
            
            fprintf(stderr,"here is one: %d\n", line->pos+1);
            exit(1);

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

            break;
        }
    }

    free(ac);
    return nrm;
}

