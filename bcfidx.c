#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "vcf.h"

int main_bcfidx(int argc, char *argv[])
{
	int c, min_shift = 14;
	while ((c = getopt(argc, argv, "s:")) >= 0)
		if (c == 's') min_shift = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: bamidx [-s minShift] <in.bam>\n");
		return 1;
	}
	bcf_index_build(argv[optind], min_shift);
	return 0;
}
