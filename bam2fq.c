#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "sam.h"
#include "kstring.h"

static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int main_bam2fq(int argc, char *argv[])
{
	BGZF *fp, *fpse = 0;
	bam1_t *b;
	uint8_t *buf;
	int max_buf, c, has12 = 0;
	kstring_t str;
	int64_t n_singletons = 0, n_reads = 0;
	char last[512], *fnse = 0;

	while ((c = getopt(argc, argv, "as:")) > 0)
		if (c == 'a') has12 = 1;
		else if (c == 's') fnse = optarg;
	if (argc == optind) {
		fprintf(stderr, "\nUsage:   bam2fq [-a] [-s outSE] <in.bam>\n\n");
		fprintf(stderr, "Options: -a        append /1 and /2 to the read name\n");
		fprintf(stderr, "         -s FILE   write singleton reads to FILE [assume single-end]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bgzf_open(argv[optind], "r") : bgzf_dopen(fileno(stdin), "r");
	assert(fp);
	bam_hdr_destroy(bam_hdr_read(fp));
	buf = 0;
	max_buf = 0;
	str.l = str.m = 0; str.s = 0;
	last[0] = 0;
	if (fnse) fpse = bgzf_open(fnse, "w1");

	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		int i, qlen = b->core.l_qseq, is_print = 0;
		uint8_t *qual, *seq;
		++n_reads;
		if (fpse) {
			if (str.l && strcmp(last, bam_get_qname(b))) {
				bgzf_write(fpse, str.s, str.l);
				str.l = 0;
				++n_singletons;
			}
			if (str.l) is_print = 1;
			strcpy(last, bam_get_qname(b));
		} else is_print = 1;
		qual = bam_get_qual(b);
		kputc(qual[0] == 0xff? '>' : '@', &str);
		kputsn(bam_get_qname(b), b->core.l_qname - 1, &str);
		if (has12) {
			kputc('/', &str);
			kputw(b->core.flag>>6&3, &str);
		}
		kputc('\n', &str);
		if (max_buf < qlen + 1) {
			max_buf = qlen + 1;
			kroundup32(max_buf);
			buf = (uint8_t*)realloc(buf, max_buf);
		}
		buf[qlen] = 0;
		seq = bam_get_seq(b);
		for (i = 0; i < qlen; ++i) buf[i] = bam_seqi(seq, i); // copy the sequence
		if (bam_is_rev(b)) { // reverse complement
			for (i = 0; i < qlen>>1; ++i) {
				int8_t t = seq_comp_table[buf[qlen - 1 - i]];
				buf[qlen - 1 - i] = seq_comp_table[buf[i]];
				buf[i] = t;
			}
			if (qlen&1) buf[i] = seq_comp_table[buf[i]];
		}
		for (i = 0; i < qlen; ++i) buf[i] = seq_nt16_str[buf[i]];
		kputsn((char*)buf, qlen, &str); kputc('\n', &str);
		if (qual[0] != 0xff) {
			kputsn("+\n", 2, &str);
			for (i = 0; i < qlen; ++i) buf[i] = 33 + qual[i];
			if (bam_is_rev(b)) { // reverse
				for (i = 0; i < qlen>>1; ++i) {
					uint8_t t = buf[qlen - 1 - i];
					buf[qlen - 1 - i] = buf[i];
					buf[i] = t;
				}
			}
		}
		kputsn((char*)buf, qlen, &str); kputc('\n', &str);
		if (is_print) {
			fwrite(str.s, 1, str.l, stdout);
			str.l = 0;
		}
	}
	if (fpse) {
		if (str.l) {
			bgzf_write(fpse, str.s, str.l);
			++n_singletons;
		}
		fprintf(stderr, "[M::%s] discarded %lld singletons\n", __func__, (long long)n_singletons);
		bgzf_close(fpse);
	}
	fprintf(stderr, "[M::%s] processed %lld reads\n", __func__, (long long)n_reads);
	free(buf); free(str.s);
	bam_destroy1(b);
	bgzf_close(fp);
	return 0;
}
