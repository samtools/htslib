/*  test/usepublic.cpp -- Test compiling public headers with a C++ compiler.

    Copyright (C) 2023 Centre for Population Genomics.

    Author: John Marshall <jmarshall@hey.com>

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

// Include *all* the public HTSlib headers.

#include "../htslib/bgzf.h"
#include "../htslib/cram.h"
#include "../htslib/faidx.h"
#include "../htslib/hfile.h"
#include "../htslib/hts.h"
#include "../htslib/hts_defs.h"
#include "../htslib/hts_endian.h"
#include "../htslib/hts_expr.h"
#include "../htslib/hts_log.h"
#include "../htslib/hts_os.h"
#include "../htslib/kbitset.h"
#include "../htslib/kfunc.h"
#include "../htslib/khash.h"
#include "../htslib/khash_str2int.h"
#include "../htslib/klist.h"
#include "../htslib/knetfile.h"
#include "../htslib/kroundup.h"
#include "../htslib/kseq.h"
#include "../htslib/ksort.h"
#include "../htslib/kstring.h"
#include "../htslib/regidx.h"
#include "../htslib/sam.h"
#include "../htslib/synced_bcf_reader.h"
#include "../htslib/tbx.h"
#include "../htslib/thread_pool.h"
#include "../htslib/vcf.h"
#include "../htslib/vcf_sweep.h"
#include "../htslib/vcfutils.h"

// Instantiate macro-based klib facilities so the resulting function
// definitions are seen by the C++ compiler.

KHASH_SET_INIT_STR(strhash)

#define noop_free(ptr)
KLIST_INIT(intlist, int, noop_free)

KSORT_INIT_STR

struct myFILE;
extern int myread(struct myFILE *, void *, int);
KSEQ_INIT2(, struct myFILE *, myread)

int main()
{
    return 0;
}
