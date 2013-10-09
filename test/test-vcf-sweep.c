#include <stdio.h>
#include <htslib/vcf_sweep.h>

int main(int argc, char **argv)
{
    if ( argc!=2 ) 
    {
        fprintf(stderr,"Usage: test-vcf-sweep <file.bcf|file.vcf>\n");
        return 1;
    }

    // Init variables. The checksum is just for this test program to output
    // something and verify that all sites are read in both passes - fwd and
    // bwd.
    bcf_sweep_t *sw = bcf_sweep_init(argv[1]);
    bcf_hdr_t *hdr  = bcf_sweep_hdr(sw);
    int chksum = 0;

    // First we must sweep forward and read the whole file to build an index.
    // If this is undesirable, we can require the presence of a .gzi index
    // which can be created with `bgzip -r` from the samtools/htslib package
    bcf1_t *rec;
    while ( (rec = bcf_sweep_fwd(sw)) ) chksum += rec->pos+1;
    printf("fwd position chksum: %d\n", chksum);

    // Now sweep backward. 
    chksum = 0;
    while ( (rec = bcf_sweep_bwd(sw)) ) chksum += rec->pos+1;
    printf("bwd position chksum: %d\n", chksum);

    // And forward and backward again, this time summing the PL vectors
    int i,j, *PLs = NULL, mPLs = 0, nPLs;
    chksum = 0;
    while ( (rec = bcf_sweep_fwd(sw)) ) 
    {
        // get copy of the PL vectors
        nPLs = bcf_get_format_int(hdr, rec, "PL", &PLs, &mPLs);
        if ( !nPLs ) continue;  // PL not present

        // how many values are there per sample
        int nvals = nPLs / bcf_hdr_nsamples(hdr);

        int *ptr = PLs;
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
        {
            for (j=0; j<nvals; j++)
            {
                // check for shorter vectors (haploid genotypes amongst diploids)
                if ( ptr[j]==bcf_int32_vector_end ) break;

                // skip missing values
                if ( ptr[j]==bcf_int32_missing ) continue;

                chksum += ptr[j];
            }
            ptr += nvals;
        }
    }
    printf("fwd PL chksum: %d\n", chksum);

    // And the same backwards..
    chksum = 0;
    while ( (rec = bcf_sweep_bwd(sw)) )
    {
        nPLs = bcf_get_format_int(hdr, rec, "PL", &PLs, &mPLs);
        if ( !nPLs ) continue;
        int nvals = nPLs / bcf_hdr_nsamples(hdr);
        int *ptr = PLs;
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
        {
            for (j=0; j<nvals; j++)
            {
                if ( ptr[j]==bcf_int32_vector_end ) break;
                if ( ptr[j]==bcf_int32_missing ) continue;
                chksum += ptr[j];
            }
            ptr += nvals;
        }
    }
    printf("bwd PL chksum: %d\n", chksum);

    // Clean up
    bcf_sweep_destroy(sw);
    return 0;
}


