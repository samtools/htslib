/*  test/test_index_load.c -- Test loading of malformed index files.

    Copyright (C) 2024 Genome Research Ltd.

    Author: Mohammed Adib <adib@bugqore.com>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "../htslib/hts.h"
#include "../htslib/bgzf.h"

static void u32(uint8_t *p, uint32_t v) {
    p[0] = v; p[1] = v >> 8; p[2] = v >> 16; p[3] = v >> 24;
}

// Write a minimal CSI index with the given min_shift and depth (n_lvls) and
// no references, used to exercise the header validity checks in idx_read().
static int write_csi(const char *fn, uint32_t min_shift, uint32_t depth) {
    uint8_t buf[20];
    memcpy(buf, "CSI\1", 4);
    u32(buf + 4,  min_shift);
    u32(buf + 8,  depth);
    u32(buf + 12, 0);           // l_aux
    u32(buf + 16, 0);           // n_ref
    BGZF *fp = bgzf_open(fn, "w");
    if (!fp) return -1;
    if (bgzf_write(fp, buf, sizeof(buf)) != (ssize_t) sizeof(buf)) {
        bgzf_close(fp);
        return -1;
    }
    return bgzf_close(fp);
}

static int check_rejected(const char *fn, uint32_t min_shift, uint32_t depth) {
    if (write_csi(fn, min_shift, depth) != 0) {
        fprintf(stderr, "Failed to write %s\n", fn);
        return 1;
    }
    hts_idx_t *idx = hts_idx_load3("test_index_load.dummy", fn,
                                   HTS_FMT_CSI, 0);
    if (idx) {
        fprintf(stderr, "Loaded malformed CSI (min_shift=%u, depth=%u) "
                "that should have been rejected\n", min_shift, depth);
        hts_idx_destroy(idx);
        return 1;
    }
    return 0;
}

static int check_accepted(const char *fn, uint32_t min_shift, uint32_t depth) {
    if (write_csi(fn, min_shift, depth) != 0) {
        fprintf(stderr, "Failed to write %s\n", fn);
        return 1;
    }
    hts_idx_t *idx = hts_idx_load3("test_index_load.dummy", fn,
                                   HTS_FMT_CSI, 0);
    if (!idx) {
        fprintf(stderr, "Rejected valid CSI (min_shift=%u, depth=%u)\n",
                min_shift, depth);
        return 1;
    }
    hts_idx_destroy(idx);
    return 0;
}

int main(int argc, char **argv) {
    char fn[1024];
    const char *dir = argc > 1 ? argv[1] : ".";
    int fail = 0;
    enum htsLogLevel logging = hts_get_log_level();

    snprintf(fn, sizeof(fn), "%s/test_index_load.csi", dir);

    // Silence the expected "Index file corrupted" errors.
    hts_set_log_level(HTS_LOG_OFF);

    // depth drives 1 << (3*depth + 3) in hts_idx_init(); > 9 overflows.
    fail |= check_rejected(fn, 14, 30);
    fail |= check_rejected(fn, 14, 0xffffffff);
    // min_shift + 3*depth feeds hts_bin_maxpos() and the query shifts.
    fail |= check_rejected(fn, 63, 5);
    fail |= check_rejected(fn, 0xffffffff, 5);

    // A normal CSI header must still load.
    fail |= check_accepted(fn, 14, 5);

    hts_set_log_level(logging);
    return fail;
}
