#include "vcfutils.h"

int calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)
{
	// Use INFO/AC,AN field only when asked
	if ( which&BCF_UN_INFO )
	{
		bcf_unpack(line, BCF_UN_INFO);
		int an_id = bcf_id2int(header, BCF_DT_ID, "AN");
		int ac_id = bcf_id2int(header, BCF_DT_ID, "AC");
		if ( an_id>=0 && ac_id>=0 )
		{
			int i, an=0, ac_len=0;
			uint8_t *ac_ptr=NULL;
			for (i=0; i<line->n_info; i++)
			{
				bcf_info_t *z = &line->d.info[i];
				if ( z->key == an_id ) an = z->v1.i;
				else if ( z->key == ac_id ) { ac_ptr = z->vptr; ac_len = z->len; }
			}
			int nac = 0;
			for (i=0; i<ac_len; i++)
			{
				ac[i+1] = ac_ptr[i]; 
				nac += ac_ptr[i];
			}
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
				if ( *p ) ac[((*p)>>1)-1]++;
				p++;
			}
		}
		return 1;
	}

	return 0;
}

