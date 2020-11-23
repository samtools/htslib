/*  expr.c -- filter expression parsing and processing.

    Copyright (C) 2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTS_EXPR_H
#define HTS_EXPR_H

#include <htslib/hts.h>
#include <htslib/kstring.h>

// fexpr_t is our return type and the type for elements within the expr.
// Note we cope with zero-but-true in order to implement a basic
// "exists(something)" check where "something" may even be zero.
//
// Eg in the aux tag searching syntax, "[NM]" should return true if
// NM tag exists even if zero.
// Take care when negating this. "[NM] != 0" will be true when
// [NM] is absent, thus consider "[NM] && [NM] != 0".
typedef struct {
    char is_str;  // Use .s vs .d
    char is_true; // Force true if even zero
    kstring_t s;  // is_str and empty s permitted (eval as false)
    double d;     // otherwise this
} fexpr_t;

#define FEXPR_INIT {0, 0, KS_INITIALIZE, 0}

// Create a SAM filter for expression "str".
//
// Returns a pointer on success,
//         NULL on failure
sam_filter_t *sam_filter_init(const char *str);

// Frees a sam_filter_t created via sam_filter_init
void sam_filter_free(sam_filter_t *filt);

typedef int (sym_func)(void *data, char *str, char **end, fexpr_t *res);
int sam_filter_eval(sam_filter_t *filt, void *data, sym_func *f, fexpr_t *res);

static inline void fexpr_free(fexpr_t *f) {
    ks_free(&f->s);
}

#endif /* HTS_EXPR_H */
