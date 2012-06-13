#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "sam.h"

int main_bamidx(int argc, char *argv[])
{
	int c, min_shift = -1;
	while ((c = getopt(argc, argv, "s:")) >= 0)
		if (c == 's') min_shift = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: bamidx [-s minShift] <in.bam>\n");
		return 1;
	}
	bam_index_build(argv[optind], 0, min_shift);
	return 0;
}
