/*  expr.c -- filter expression parsing and processing.

    Copyright (C) 2020, 2022 Genome Research Ltd.

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

#include <math.h>
#include "kstring.h"
#include "hts_defs.h"

/// Holds a filter variable.  This is also used to return the results.
/**
 * The expression language has 3-states of string, numeric, and unknown.
 * The unknown state is either a NaN numeric or a null string, with both
 * internally considered to have the same "unknown" meaning.
 *
 * These largely match the IEE 754 semantics for NaN comparisons: <, >, ==,
 * != all fail, (even NaN == NaN).  Similarly arithmetic (+,-,/,*,%) with
 * unknown values are still unknown (and false).
 *
 * The departure from NaN semantics though is that our unknown/null state is
 * considered to be false while NaN in C is true.  Similarly the false nature
 * of our unknown state meants !val becomes true, !!val is once again false,
 * val && 1 is false, val || 0 is false, and val || 1 is true along with
 * !val || 0 and !val && 1.
 *
 * Note it is possible for empty strings and zero numbers to also be true.
 * An example of this is the aux string '[NM]' which returns true if the
 * NM tag is found, regardless of whether it is also zero.  However the
 * better approach added in 1.16 is 'exists([NM])'.
 */
typedef struct hts_expr_val_t {
    char is_str;  // Use .s vs .d
    char is_true; // Force true if even zero
    kstring_t s;  // is_str and empty s permitted (eval as false)
    double d;     // otherwise this
} hts_expr_val_t;

/// Returns true if an hts_expr_val_t is defined.
/* An example usage of this is in the SAM expression filter where an
 * [X0] aux tag will be the value of X0 (string or numeric) if set, or
 * a false nul-string (not the same as an empty one) when not set.
 */
static inline int hts_expr_val_exists(hts_expr_val_t *v) {
    return v && !(v->is_str == 1 && v->s.s == NULL)
             && !(v->is_str == 0 && isnan(v->d));
}

/// Returns true if an hts_expr_val_t is defined or is undef-but-true
static inline int hts_expr_val_existsT(hts_expr_val_t *v) {
    return (v && v->is_true) || hts_expr_val_exists(v);
}

/// Set a value to be undefined (nan).
static inline void hts_expr_val_undef(hts_expr_val_t *v) {
    ks_clear(&v->s);
    v->is_true = 0;
    v->is_str = 0;
    v->d = NAN;
}

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
 *
 *  @p res must be initialized using HTS_EXPR_VAL_INIT before passing it
 *  to this function for the first time.
 */
HTSLIB_EXPORT
int hts_filter_eval2(hts_filter_t *filt,
                     void *data, hts_expr_sym_func *sym_func,
                     hts_expr_val_t *res);

/// Evaluate a filter expression (derecated API)
/**
 *  @copydetails hts_filter_eval2()
 *
 *  If calling this function more than once with the same @p res
 *  parameter, hts_expr_val_free(res) must be used between invocations
 *  to clear any allocated memory prior to reuse.
 *
 *  @deprecated This function has been replaced by hts_filter_eval2(),
 *              which clears @p res properly itself.
 */
HTSLIB_EXPORT
int hts_filter_eval(hts_filter_t *filt,
                    void *data, hts_expr_sym_func *sym_func,
                    hts_expr_val_t *res)
    HTS_DEPRECATED("Please use hts_filter_eval2 instead");


#endif /* HTS_EXPR_H */
