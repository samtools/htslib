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

#include "kstring.h"
#include "hts_defs.h"

/// Holds a filter variable.  This is also used to return the results.
/**
 * Note we cope with zero-but-true in order to implement a basic
 * "exists(something)" check where "something" may even be zero.
 *
 * Eg in the aux tag searching syntax, "[NM]" should return true if
 * NM tag exists even if zero.
 * Take care when negating this. "[NM] != 0" will be true when
 * [NM] is absent, thus consider "[NM] && [NM] != 0".
 */
typedef struct hts_expr_val_t {
    char is_str;  // Use .s vs .d
    char is_true; // Force true if even zero
    kstring_t s;  // is_str and empty s permitted (eval as false)
    double d;     // otherwise this
} hts_expr_val_t;

/// Frees a hts_expr_val_t type.
static inline void hts_expr_val_free(hts_expr_val_t *f) {
    ks_free(&f->s);
}

/// Opaque hts_filter_t type.  Definition in hts_expr.c
typedef struct hts_filter_t hts_filter_t;

/// For static initialisation of hts_expr_val_t values
#define HTS_EXPR_VAL_INIT {0, 0, KS_INITIALIZE, 0}

/// Creates a filter for expression "str".
/** @param str    The filter expression
 *  @return       A pointer on success, NULL on failure
 */
HTSLIB_EXPORT
hts_filter_t *hts_filter_init(const char *str);

/// Frees an hts_filter_t created via hts_filter_init
/** @param filt    The filter pointer.
 */
HTSLIB_EXPORT
void hts_filter_free(hts_filter_t *filt);

/// Type for expression symbol lookups; name -> value.
typedef int (hts_expr_sym_func)(void *data, char *str, char **end,
                                hts_expr_val_t *res);

/// Evaluates a filter expression and returns the value
/** @param filt      The filter, produced by hts_filter_init
 *  @param data      Arbitrary caller data, passed into sym_func
 *  @param sym_func  Callback function to lookup variables.
 *  @param res       Filled out with the result of the filter evaluation
 *  @return          Returns 0 on success, -1 on failure
 *
 *  sym_func and data may be NULL if the caller does not need its own data
 *  pointer or if it has no variables to lookup.
 *
 *  The type of the returned result may be numeric of string, as defined by
 *  the is_str member.  It can also be explicitly defined to be true even
 *  for a null value.  This may be used to check for the existence of
 *  something, irrespective of whether that something evaluates to zero.
 */
HTSLIB_EXPORT
int hts_filter_eval(hts_filter_t *filt,
                    void *data, hts_expr_sym_func *sym_func,
                    hts_expr_val_t *res);


#endif /* HTS_EXPR_H */
