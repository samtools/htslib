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
		fprintf(stderr, "\nUsage: bamidx [-s minBits] <in.bam>\n\n");
		fprintf(stderr, "Note:  The minimal interval size equals 1<<minBits. If '-s' is unset, bamidx\n\
       creates an old BAM index with file extension '.bai'; otherwise it\n\
       writes a new index with extension '.csi'.\n\n");
		return 1;
	}
	bam_index_build(argv[optind], min_shift);
	return 0;
}
