/*
    Copyright (C) 2018-2019 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Test region description parser.
    Usage: test-parse-reg [-c] file.bam region
           test-parse-reg [-c] -m file.bam region,region...
           test-parse-reg -t

    -c is chr:pos is a single base coordinate, ie chr:pos-pos,
       otherwise it is chr:pos-<end>
    -m is multi-region list.
    -t runs built-in tests

    ./test/test-parse-reg -c -m test/colons.bam "{chr1:100-200},{chr1}:100-200,{chr1:100-200}:100,{chr1,chr3},chr1:"
*/

#include <config.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

void reg_expected(sam_hdr_t *hdr, const char *reg, int flags,
                 char *reg_exp, int tid_exp, hts_pos_t beg_exp, hts_pos_t end_exp) {
    const char *reg_out;
    int tid_out = -1;
    hts_pos_t beg_out = -1, end_out = -1;

    reg_out = sam_parse_region(hdr, reg, &tid_out, &beg_out, &end_out, flags);

    if ((reg_out != NULL) != (reg_exp != NULL) ||
        (reg_out && reg_exp && strcmp(reg_out, reg_exp) != 0) ||
        (reg_exp && tid_out != tid_exp) ||
        (reg_exp && beg_out != beg_exp) ||
        (reg_exp && end_out != end_exp)) {
        fprintf(stderr, "Parsing \"%s\" expected return \"%s\", %d:%"PRIhts_pos"-%"PRIhts_pos", "
                "but got \"%s\", %d:%"PRIhts_pos"-%"PRIhts_pos"\n",
                reg,
                reg_exp?reg_exp:"(null)", tid_exp, beg_exp, end_exp,
                reg_out?reg_out:"(null)", tid_out, beg_out, end_out);
        exit(1);
    }
}

int reg_test(char *fn) {
    samFile *fp;
    sam_hdr_t *hdr;

    if (!(fp = sam_open(fn, "r")))
        return 1;

    if (!(hdr = sam_hdr_read(fp)))
        return 1;

    // 0 chr1
    // 1 chr1:100
    // 2 chr1:100-200
    // 3 chr2:100-200
    // 4 chr3
    // 5 chr1,chr3

    // Check range extensions.
    reg_expected(hdr, "chr1", 0, "",  0, 0, HTS_POS_MAX);
    reg_expected(hdr, "chr1:50", 0, "",  0, 49, HTS_POS_MAX);
    reg_expected(hdr, "chr1:50", HTS_PARSE_ONE_COORD, "",  0, 49, 50);
    reg_expected(hdr, "chr1:50-100", 0, "",  0, 49, 100);
    reg_expected(hdr, "chr1:50-", 0, "",  0, 49, HTS_POS_MAX);
    reg_expected(hdr, "chr1:-50", 0, "",  0, 0, 50);

    // Check quoting
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "chr1:100-200", 0, NULL,  0, 0, 0); // ambiguous
    reg_expected(hdr, "{chr1}:100-200", 0, "",  0, 99, 200);
    reg_expected(hdr, "{chr1:100-200}", 0, "",  2, 0, HTS_POS_MAX);
    reg_expected(hdr, "{chr1:100-200}:100-200", 0, "",  2, 99, 200);
    reg_expected(hdr, "{chr2:100-200}:100-200", 0, "",  3, 99, 200);
    reg_expected(hdr, "chr2:100-200:100-200", 0, "",  3, 99, 200);
    reg_expected(hdr, "chr2:100-200", 0, "",  3, 0, HTS_POS_MAX);

    // Check numerics
    reg_expected(hdr, "chr3", 0, "",  4, 0, HTS_POS_MAX);
    reg_expected(hdr, "chr3:", 0, "",  4, 0, HTS_POS_MAX);
    reg_expected(hdr, "chr3:1000-1500", 0, "",  4, 999, 1500);
    reg_expected(hdr, "chr3:1,000-1,500", 0, "",  4, 999, 1500);
    reg_expected(hdr, "chr3:1k-1.5K", 0, "",  4, 999, 1500);
    reg_expected(hdr, "chr3:1e3-1.5e3", 0, "",  4, 999, 1500);
    reg_expected(hdr, "chr3:1e3-15e2", 0, "",  4, 999, 1500);

    // Check list mode
    reg_expected(hdr, "chr1,chr3", HTS_PARSE_LIST, "chr3", 0, 0, HTS_POS_MAX);
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "chr1:100-200,chr3", HTS_PARSE_LIST, NULL,  0, 0, 0); // ambiguous
    reg_expected(hdr, "{chr1,chr3}", HTS_PARSE_LIST, "", 5, 0, HTS_POS_MAX);
    reg_expected(hdr, "{chr1,chr3},chr1", HTS_PARSE_LIST, "chr1", 5, 0, HTS_POS_MAX);
    // incorrect usage; first reg is valid (but not what user expects).
    reg_expected(hdr, "chr3:1,000-1,500", HTS_PARSE_LIST | HTS_PARSE_ONE_COORD, "000-1,500",  4, 0, 1);

    // More expected failures
    reg_expected(hdr, "chr2", 0, NULL, 0, 0, 0);
    reg_expected(hdr, "chr1,", 0, NULL, 0, 0, 0);
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "{chr1", 0, NULL, 0, 0, 0);
    reg_expected(hdr, "chr1:10-10", 0, "", 0, 9, 10); // OK
    reg_expected(hdr, "chr1:10-9", 0, NULL, 0, 0, 0); // Issue#353
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "chr1:x", 0, NULL, 0, 0, 0);
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "chr1:1-y", 0, NULL, 0, 0, 0);
    fprintf(stderr, "Expected error: ");
    reg_expected(hdr, "chr1:1,chr3", 0, NULL, 0, 0, 0);

    sam_hdr_destroy(hdr);
    sam_close(fp);

    exit(0);
}

int main(int argc, char **argv) {
    sam_hdr_t *hdr;
    samFile *fp;
    int flags = 0;

    while (argc > 1) {
        if (strcmp(argv[1], "-m") == 0) {
            flags |= HTS_PARSE_LIST;
            argc--; argv++;
            continue;
        }

        if (strcmp(argv[1], "-c") == 0) {
            flags |= HTS_PARSE_ONE_COORD;
            argc--; argv++;
            continue;
        }

        // Automatic mode for test harness
        if (strcmp(argv[1], "-t") == 0)
            reg_test(argv[2]);

        break;
    }

    // Interactive mode for debugging
    if (argc != 3) {
        fprintf(stderr, "Usage: test-parse-reg [-m] [-c] region[,region]...\n");
        exit(1);
    }

    if (!(fp = sam_open(argv[1], "r"))) {
        perror(argv[1]);
        exit(1);
    }

    if (!(hdr = sam_hdr_read(fp))) {
        fprintf(stderr, "Couldn't read header\n");
        exit(1);
    }

    const char *reg = argv[2];
    while (*reg) {
        int tid;
        hts_pos_t beg, end;
        reg = sam_parse_region(hdr, reg, &tid, &beg, &end, flags);
        if (!reg) {
            fprintf(stderr, "Failed to parse region\n");
            exit(1);
        }
        printf("%-20s %12"PRIhts_pos" %12"PRIhts_pos"\n",
               tid == -1 ? "*" : hdr->target_name[tid],
               beg, end);
    }

    sam_hdr_destroy(hdr);
    sam_close(fp);

    return 0;
}
