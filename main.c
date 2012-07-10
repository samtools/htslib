#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "vcf.h"

int main_samview(int argc, char *argv[]);
int main_vcfview(int argc, char *argv[]);
int main_bamidx(int argc, char *argv[]);
int main_bcfidx(int argc, char *argv[]);
int main_bamshuf(int argc, char *argv[]);
int main_bam2fq(int argc, char *argv[]);
int main_tabix(int argc, char *argv[]);
int main_sam2break(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "\nUsage:   htscmd <command> <argument>\n\n");
	fprintf(stderr, "Command: samview      SAM<->BAM conversion\n");
	fprintf(stderr, "         vcfview      VCF<->BCF conversion (unfinished)\n");
	fprintf(stderr, "         tabix        tabix for BGZF'd BED, GFF, SAM, VCF and more\n");
	fprintf(stderr, "         bamidx       index BAM\n");
	fprintf(stderr, "         bcfidx       index BCF\n\n");
	fprintf(stderr, "         bamshuf      shuffle BAM and group alignments by query name\n");
	fprintf(stderr, "         bam2fq       convert name grouped BAM to interleaved fastq\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "samview") == 0) return main_samview(argc-1, argv+1);
	else if (strcmp(argv[1], "vcfview") == 0) return main_vcfview(argc-1, argv+1);
	else if (strcmp(argv[1], "bamidx") == 0) return main_bamidx(argc-1, argv+1);
	else if (strcmp(argv[1], "bcfidx") == 0) return main_bcfidx(argc-1, argv+1);
	else if (strcmp(argv[1], "bamshuf") == 0) return main_bamshuf(argc-1, argv+1);
	else if (strcmp(argv[1], "bam2fq") == 0) return main_bam2fq(argc-1, argv+1);
	else if (strcmp(argv[1], "tabix") == 0) return main_tabix(argc-1, argv+1);
	else if (strcmp(argv[1], "sam2break") == 0) return main_sam2break(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
}
