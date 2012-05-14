#include <stdio.h>
#include "sam.h"

int main_bamidx(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: bamidx <in.bam>\n");
		return 1;
	}
	sam_index_build(argv[1], 0);
	return 0;
}
