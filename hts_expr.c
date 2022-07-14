/*  hts_expr.c -- filter expression parsing and processing.

    Copyright (C) 2020-2022 Genome Research Ltd.

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

// TODO:
// - ?: operator for conditionals?

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <regex.h>
#include <math.h>

#include "htslib/hts_expr.h"
#include "htslib/hts_log.h"
#include "textutils_internal.h"

// Could also cache hts_expr_val_t stack here for kstring reuse?
#define MAX_REGEX 10
struct hts_filter_t {
    char *str;
    int parsed;
    int curr_regex, max_regex;
    regex_t preg[MAX_REGEX];
};

/*
 * This is designed to be mostly C like with mostly same the precedence rules,
 * with the exception of bit operators (widely considered as a mistake in C).
 * It's not full C (eg no bit-shifting), but good enough for our purposes.
 *
 * Supported syntax, in order of precedence:
 *
 * Grouping:      (, ),   eg "(1+2)*3"
 * Values:        integers, floats, strings or variables
 * Unary ops:     +, -, !, ~  eg -10 +10, !10 (0), ~5 (bitwise not)
 * Math ops:      *, /, %  [TODO: add // for floor division?]
 * Math ops:      +, -
 * Bit-wise:      &, ^, |  [NB as 3 precedence levels, in that order]
 * Conditionals:  >, >=, <, <=,
 * Equality:      ==, !=, =~, !~
 * Boolean:       &&, ||
 */

// Skip to start of term
static char *ws(char *str) {
    while (*str && (*str == ' ' || *str == '\t'))
        str++;
    return str;
}

static int expression(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                      char *str, char **end, hts_expr_val_t *res);

/*
 * Simple functions operating on strings only.
 * length, min, max, avg.
 *
 * All return 0 on success,
 *           -1 on failure
 */
static int expr_func_length(hts_expr_val_t *res) {
    if (!res->is_str)
        return -1;

    res->is_str = 0;
    res->d = res->s.l;
    return 0;
}

static int expr_func_min(hts_expr_val_t *res) {
    if (!res->is_str)
        return -1;

    size_t l = res->s.l;
    int v = INT_MAX;
    const uint8_t *x = (uint8_t *)res->s.s;
    for (l = 0; l < res->s.l; l++)
        if (v > x[l])
            v = x[l];

    res->is_str = 0;
    res->d = v == INT_MAX ? NAN : v;

    return 0;
}

static int expr_func_max(hts_expr_val_t *res) {
    if (!res->is_str)
        return -1;

    size_t l = res->s.l;
    int v = INT_MIN;
    const uint8_t *x = (uint8_t *)res->s.s;
    for (l = 0; l < res->s.l; l++)
        if (v < x[l])
            v = x[l];

    res->is_str = 0;
    res->d = v == INT_MIN ? NAN : v;

    return 0;
}

static int expr_func_avg(hts_expr_val_t *res) {
    if (!res->is_str)
        return -1;

    size_t l = res->s.l;
    double v = 0;
    const uint8_t *x = (uint8_t *)res->s.s;
    for (l = 0; l < res->s.l; l++)
        v += x[l];
    if (l)
        v /= l;

    res->is_str = 0;
    res->d = v;

    return 0;
}

/*
 * functions:  FUNC(expr).
 * Note for simplicity of parsing, the "(" must immediately follow FUNC,
 * so "FUNC (x)" is invalid.
 */
static int func_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                     char *str, char **end, hts_expr_val_t *res) {
    int func_ok = -1;
    switch (*str) {
    case 'a':
        if (strncmp(str, "avg(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = expr_func_avg(res);
        }
        break;

    case 'd':
        if (strncmp(str, "default(", 8) == 0) {
            if (expression(filt, data, fn, str+8, end, res)) return -1;
            if (**end != ',')
                return -1;
            (*end)++;
            hts_expr_val_t val = HTS_EXPR_VAL_INIT;
            if (expression(filt, data, fn, ws(*end), end, &val)) return -1;
            func_ok = 1;
            if (!hts_expr_val_existsT(res)) {
                kstring_t swap = res->s;
                *res = val;
                val.s = swap;
                hts_expr_val_free(&val);
            }
        }
        break;

    case 'e':
        if (strncmp(str, "exists(", 7) == 0) {
            if (expression(filt, data, fn, str+7, end, res)) return -1;
            func_ok = 1;
            res->is_true = res->d = hts_expr_val_existsT(res);
            res->is_str = 0;
        } else if (strncmp(str, "exp(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = 1;
            res->d = exp(res->d);
            res->is_str = 0;
            if (isnan(res->d))
                hts_expr_val_undef(res);
        }

        break;

    case 'l':
        if (strncmp(str, "length(", 7) == 0) {
            if (expression(filt, data, fn, str+7, end, res)) return -1;
            func_ok = expr_func_length(res);
        } else if (strncmp(str, "log(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = 1;
            res->d = log(res->d);
            res->is_str = 0;
            if (isnan(res->d))
                hts_expr_val_undef(res);
        }
        break;

    case 'm':
        if (strncmp(str, "min(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = expr_func_min(res);
        } else if (strncmp(str, "max(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = expr_func_max(res);
        }
        break;

    case 'p':
        if (strncmp(str, "pow(", 4) == 0) {
            if (expression(filt, data, fn, str+4, end, res)) return -1;
            func_ok = 1;

            if (**end != ',')
                return -1;
            (*end)++;
            hts_expr_val_t val = HTS_EXPR_VAL_INIT;
            if (expression(filt, data, fn, ws(*end), end, &val)) return -1;
            if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
                hts_expr_val_undef(res);
            } else if (res->is_str || val.is_str) {
                hts_expr_val_free(&val); // arith on strings
                return -1;
            } else {
                func_ok = 1;
                res->d = pow(res->d, val.d);
                hts_expr_val_free(&val);
                res->is_str = 0;
            }

            if (isnan(res->d))
                hts_expr_val_undef(res);
        }
        break;

    case 's':
        if (strncmp(str, "sqrt(", 5) == 0) {
            if (expression(filt, data, fn, str+5, end, res)) return -1;
            func_ok = 1;
            res->d = sqrt(res->d);
            res->is_str = 0;
            if (isnan(res->d))
                hts_expr_val_undef(res);
        }
        break;
    }

    if (func_ok < 0)
        return -1;

    str = ws(*end);
    if (*str != ')') {
        fprintf(stderr, "Missing ')'\n");
        return -1;
    }
    *end = str+1;

    return 0;
}

/*
 * simple_expr
 *     : identifier
 *     | constant
 *     | string
 *     | func_expr
 *     | '(' expression ')'
*/
static int simple_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                       char *str, char **end, hts_expr_val_t *res) {
    // Main recursion step
    str = ws(str);
    if (*str == '(') {
        if (expression(filt, data, fn, str+1, end, res)) return -1;
        str = ws(*end);
        if (*str != ')') {
            fprintf(stderr, "Missing ')'\n");
            return -1;
        }
        *end = str+1;

        return 0;
    }

    // Otherwise a basic element.
    int fail = 0;
    double d = hts_str2dbl(str, end, &fail);
    if (str != *end) {
        res->is_str = 0;
        res->d = d;
    } else {
        // Not valid floating point syntax.
        // TODO: add function call names in here; len(), sqrt(), pow(), etc
        if (*str == '"') {
            res->is_str = 1;
            char *e = str+1;
            int backslash = 0;
            while (*e && *e != '"') {
                if (*e == '\\')
                    backslash=1, e+=1+(e[1]!='\0');
                else
                    e++;
            }

            kputsn(str+1, e-(str+1), ks_clear(&res->s));
            if (backslash) {
                size_t i, j;
                for (i = j = 0; i < res->s.l; i++) {
                    res->s.s[j++] = res->s.s[i];
                    if (res->s.s[i] == '\\') {
                        switch (res->s.s[++i]) {
                        case '"': res->s.s[j-1] = '"'; break;
                        case '\\':res->s.s[j-1] = '\\'; break;
                        case 't': res->s.s[j-1] = '\t'; break;
                        case 'n': res->s.s[j-1] = '\n'; break;
                        case 'r': res->s.s[j-1] = '\r'; break;
                        default:  res->s.s[j++] = res->s.s[i];
                        }
                    }
                }
                res->s.s[j] = 0;
                res->s.l = j;
            }
            if (*e != '"')
                return -1;
            *end = e+1;
        } else if (fn) {
            // Try lookup as variable, if not as function
            if (fn(data, str, end, res) == 0)
                return 0;
            else
                return func_expr(filt, data, fn, str, end, res);
        } else {
            return -1;
        }
    }

    return 0;
}

/*
 * unary_expr
 *     : simple_expr
 *     | '+' simple_expr
 *     | '-' simple_expr
 *     | '!' unary_expr // higher precedence
 *     | '~' unary_expr // higher precedence
 */
static int unary_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                      char *str, char **end, hts_expr_val_t *res) {
    int err;
    str = ws(str);
    if (*str == '+' || *str == '-') {
        err = simple_expr(filt, data, fn, str+1, end, res);
        if (!hts_expr_val_exists(res)) {
            hts_expr_val_undef(res);
        } else {
            err |= res->is_str;
            if (*str == '-')
                res->d = -res->d;
            res->is_true = res->d != 0;
        }
    } else if (*str == '!') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        if (res->is_true) {
            // Any explicitly true value becomes false
            res->d = res->is_true = 0;
        } else if (!hts_expr_val_exists(res)) {
            // We can also still negate undef values by toggling the
            // is_true override value.
            res->d = res->is_true = !res->is_true;
        } else if (res->is_str) {
            // !null = true, !"foo" = false, NOTE: !"" = false also
            res->d = res->is_true = (res->s.s == NULL);
        } else {
            res->d = !(int64_t)res->d;
            res->is_true = res->d != 0;
        }
        res->is_str = 0;
    } else if (*str == '~') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        if (!hts_expr_val_exists(res)) {
            hts_expr_val_undef(res);
        } else {
            err |= res->is_str;
            if (!hts_expr_val_exists(res)) {
                hts_expr_val_undef(res);
            } else {
                res->d = ~(int64_t)res->d;
                res->is_true = res->d != 0;
            }
        }
    } else {
        err = simple_expr(filt, data, fn, str, end, res);
    }
    return err ? -1 : 0;
}


/*
 * mul_expr
 *     : unary_expr (
 *           '*' unary_expr
 *         | '/' unary_expr
 *         | '%' unary_expr
 *       )*
 */
static int mul_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                    char *str, char **end, hts_expr_val_t *res) {
    if (unary_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    while (*str) {
        str = ws(str);
        if (*str == '*' || *str == '/' || *str == '%') {
            if (unary_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)) {
                hts_expr_val_undef(res);
            } else if (val.is_str || res->is_str) {
                hts_expr_val_free(&val);
                return -1; // arith on strings
            }
        }

        if (*str == '*')
            res->d *= val.d;
        else if (*str == '/')
            res->d /= val.d;
        else if (*str == '%') {
            if (val.d)
                res->d = (int64_t)res->d % (int64_t)val.d;
            else
                hts_expr_val_undef(res);
        } else
            break;

        res->is_true = hts_expr_val_exists(res) && (res->d != 0);
        str = *end;
    }

    hts_expr_val_free(&val);

    return 0;
}

/*
 * add_expr
 *     : mul_expr (
 *           '+' mul_expr
 *         | '-' mul_expr
 *       )*
 */
static int add_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                    char *str, char **end, hts_expr_val_t *res) {
    if (mul_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    while (*str) {
        str = ws(str);
        int undef = 0;
        if (*str == '+' || *str == '-') {
            if (mul_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)) {
                undef = 1;
            } else if (val.is_str || res->is_str) {
                hts_expr_val_free(&val);
                return -1; // arith on strings
            }
        }

        if (*str == '+')
            res->d += val.d;
        else if (*str == '-')
            res->d -= val.d;
        else
            break;

        if (undef)
            hts_expr_val_undef(res);
        else
            res->is_true = res->d != 0;

        str = *end;
    }

    hts_expr_val_free(&val);

    return 0;
}

/*
 * bitand_expr
 *     : add_expr
 *     | bitand_expr '&' add_expr
 */
static int bitand_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                       char *str, char **end, hts_expr_val_t *res) {
    if (add_expr(filt, data, fn, str, end, res)) return -1;

    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    int undef = 0;
    for (;;) {
        str = ws(*end);
        if (*str == '&' && str[1] != '&') {
            if (add_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)) {
                undef = 1;
            } else if (res->is_str || val.is_str) {
                hts_expr_val_free(&val);
                return -1;
            }
            res->is_true = (res->d = ((int64_t)res->d & (int64_t)val.d)) != 0;
        } else {
            break;
        }
    }
    hts_expr_val_free(&val);
    if (undef)
        hts_expr_val_undef(res);

    return 0;
}

/*
 * bitxor_expr
 *     : bitand_expr
 *     | bitxor_expr '^' bitand_expr
 */
static int bitxor_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                       char *str, char **end, hts_expr_val_t *res) {
    if (bitand_expr(filt, data, fn, str, end, res)) return -1;

    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    int undef = 0;
    for (;;) {
        str = ws(*end);
        if (*str == '^') {
            if (bitand_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)) {
                undef = 1;
            } else if (res->is_str || val.is_str) {
                hts_expr_val_free(&val);
                return -1;
            }
            res->is_true = (res->d = ((int64_t)res->d ^ (int64_t)val.d)) != 0;
        } else {
            break;
        }
    }
    hts_expr_val_free(&val);
    if (undef)
        hts_expr_val_undef(res);

    return 0;
}

/*
 * bitor_expr
 *     : bitxor_expr
 *     | bitor_expr '|' bitxor_expr
 */
static int bitor_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                      char *str, char **end, hts_expr_val_t *res) {
    if (bitxor_expr(filt, data, fn, str, end, res)) return -1;

    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    int undef = 0;
    for (;;) {
        str = ws(*end);
        if (*str == '|' && str[1] != '|') {
            if (bitxor_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)) {
                undef = 1;
            } else if (res->is_str || val.is_str) {
                hts_expr_val_free(&val);
                return -1;
            }
            res->is_true = (res->d = ((int64_t)res->d | (int64_t)val.d)) != 0;
        } else {
            break;
        }
    }
    hts_expr_val_free(&val);
    if (undef)
        hts_expr_val_undef(res);

    return 0;
}

/*
 * cmp_expr
 *     : bitor_expr
 *     | cmp_expr '<=' bitor_expr
 *     | cmp_expr '<'  bitor_expr
 *     | cmp_expr '>=' bitor_expr
 *     | cmp_expr '>'  bitor_expr
 */
static int cmp_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                    char *str, char **end, hts_expr_val_t *res) {
    if (bitor_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);
    hts_expr_val_t val = HTS_EXPR_VAL_INIT;
    int err = 0, cmp_done = 0;

    if (*str == '>' && str[1] == '=') {
        cmp_done = 1;
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
            hts_expr_val_undef(res);
        } else {
            res->is_true=res->d
                = res->is_str && res->s.s && val.is_str && val.s.s
                ? strcmp(res->s.s, val.s.s) >= 0
                : !res->is_str && !val.is_str && res->d >= val.d;
            res->is_str = 0;
        }
    } else if (*str == '>') {
        cmp_done = 1;
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
            hts_expr_val_undef(res);
        } else {
            res->is_true=res->d
                = res->is_str && res->s.s && val.is_str && val.s.s
                ? strcmp(res->s.s, val.s.s) > 0
                : !res->is_str && !val.is_str && res->d > val.d;
            res->is_str = 0;
        }
    } else if (*str == '<' && str[1] == '=') {
        cmp_done = 1;
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
            hts_expr_val_undef(res);
        } else {
            res->is_true=res->d
                = res->is_str && res->s.s && val.is_str && val.s.s
                ? strcmp(res->s.s, val.s.s) <= 0
                : !res->is_str && !val.is_str && res->d <= val.d;
            res->is_str = 0;
        }
    } else if (*str == '<') {
        cmp_done = 1;
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
            hts_expr_val_undef(res);
        } else {
            res->is_true=res->d
                = res->is_str && res->s.s && val.is_str && val.s.s
                ? strcmp(res->s.s, val.s.s) < 0
                : !res->is_str && !val.is_str && res->d < val.d;
            res->is_str = 0;
        }
    }

    if (cmp_done && (!hts_expr_val_exists(&val) || !hts_expr_val_exists(res)))
        hts_expr_val_undef(res);
    hts_expr_val_free(&val);

    return err ? -1 : 0;
}

/*
 * eq_expr
 *     : cmp_expr
 *     | eq_expr '==' cmp_expr
 *     | eq_expr '!=' cmp_expr
 *     | eq_expr '=~' cmp_expr
 *     | eq_expr '!~' cmp_expr
 */
static int eq_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                   char *str, char **end, hts_expr_val_t *res) {
    if (cmp_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);

    int err = 0, eq_done = 0;
    hts_expr_val_t val = HTS_EXPR_VAL_INIT;

    // numeric vs numeric comparison is as expected
    // string vs string comparison is as expected
    // numeric vs string is false
    if (str[0] == '=' && str[1] == '=') {
        eq_done = 1;
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
                hts_expr_val_undef(res);
            } else {
                res->is_true = res->d = res->is_str
                    ? (res->s.s && val.s.s ?strcmp(res->s.s, val.s.s)==0 :0)
                    : !res->is_str && !val.is_str && res->d == val.d;
            }
        }
        res->is_str = 0;

    } else if (str[0] == '!' && str[1] == '=') {
        eq_done = 1;
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            if (!hts_expr_val_exists(res) || !hts_expr_val_exists(&val)) {
                hts_expr_val_undef(res);
            } else {
                res->is_true = res->d = res->is_str
                    ? (res->s.s && val.s.s ?strcmp(res->s.s, val.s.s) != 0 :1)
                    : res->is_str != val.is_str || res->d != val.d;
            }
        }
        res->is_str = 0;

    } else if ((str[0] == '=' && str[1] == '~') ||
               (str[0] == '!' && str[1] == '~')) {
        eq_done = 1;
        err = eq_expr(filt, data, fn, str+2, end, &val);
        if (!val.is_str || !res->is_str) {
            hts_expr_val_free(&val);
            return -1;
        }
        if (val.s.s && res->s.s && val.is_true >= 0 && res->is_true >= 0) {
            regex_t preg_, *preg;
            if (filt->curr_regex >= filt->max_regex) {
                // Compile regex if not seen before
                if (filt->curr_regex >= MAX_REGEX) {
                    preg = &preg_;
                } else {
                    preg = &filt->preg[filt->curr_regex];
                    filt->max_regex++;
                }

                int ec = regcomp(preg, val.s.s, REG_EXTENDED | REG_NOSUB);
                if (ec != 0) {
                    char errbuf[1024];
                    regerror(ec, preg, errbuf, 1024);
                    fprintf(stderr, "Failed regex: %.1024s\n", errbuf);
                    hts_expr_val_free(&val);
                    return -1;
                }
            } else {
                preg = &filt->preg[filt->curr_regex];
            }
            res->is_true = res->d = regexec(preg, res->s.s, 0, NULL, 0) == 0
                ? *str == '='  // matcn
                : *str == '!'; // no-match
            if (preg == &preg_)
                regfree(preg);

            filt->curr_regex++;
        } else {
            // nul regexp or input is considered false
            res->is_true = 0;
        }
        res->is_str = 0;
    }

    if (eq_done && ((!hts_expr_val_exists(&val)) || !hts_expr_val_exists(res)))
        hts_expr_val_undef(res);
    hts_expr_val_free(&val);

    return err ? -1 : 0;
}

/*
 * and_expr
 *     : eq_expr
 *     | and_expr 'and' eq_expr
 *     | and_expr 'or'  eq_expr
 */
static int and_expr(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                    char *str, char **end, hts_expr_val_t *res) {
    if (eq_expr(filt, data, fn, str, end, res)) return -1;

    for (;;) {
        hts_expr_val_t val = HTS_EXPR_VAL_INIT;
        str = ws(*end);
        if (str[0] == '&' && str[1] == '&') {
            if (eq_expr(filt, data, fn, str+2, end, &val)) return -1;
            if (!hts_expr_val_existsT(res) || !hts_expr_val_existsT(&val)) {
                hts_expr_val_undef(res);
                res->d = 0;
            } else {
                res->is_true = res->d =
                    (res->is_true || (res->is_str && res->s.s) || res->d) &&
                    (val.is_true  || (val.is_str && val.s.s) || val.d);
                res->is_str = 0;
            }
        } else if (str[0] == '|' && str[1] == '|') {
            if (eq_expr(filt, data, fn, str+2, end, &val)) return -1;
            if (!hts_expr_val_existsT(res) && !hts_expr_val_existsT(&val)) {
                // neither defined
                hts_expr_val_undef(res);
                res->d = 0;
            } else if (!hts_expr_val_existsT(res) &&
                       !(val.is_true  || (val.is_str  && val.s.s ) || val.d)) {
                // LHS undef and RHS false
                hts_expr_val_undef(res);
                res->d = 0;
            } else if (!hts_expr_val_existsT(&val) &&
                       !(res->is_true || (res->is_str && res->s.s) || res->d)){
                // RHS undef and LHS false
                hts_expr_val_undef(res);
                res->d = 0;
            } else {
                res->is_true = res->d =
                    res->is_true || (res->is_str && res->s.s) || res->d ||
                    val.is_true  || (val.is_str  && val.s.s ) || val.d;
                res->is_str = 0;
            }
        } else {
            break;
        }
        hts_expr_val_free(&val);
    }

    return 0;
}

static int expression(hts_filter_t *filt, void *data, hts_expr_sym_func *fn,
                      char *str, char **end, hts_expr_val_t *res) {
    return and_expr(filt, data, fn, str, end, res);
}

hts_filter_t *hts_filter_init(const char *str) {
    hts_filter_t *f = calloc(1, sizeof(*f));
    if (!f) return NULL;

    // Oversize to permit faster comparisons with memcmp over strcmp
    size_t len = strlen(str)+100;
    if (!(f->str = malloc(len))) {
        free(f);
        return NULL;
    }
    strcpy(f->str, str);
    return f;
}

void hts_filter_free(hts_filter_t *filt) {
    if (!filt)
        return;

    int i;
    for (i = 0; i < filt->max_regex; i++)
        regfree(&filt->preg[i]);

    free(filt->str);
    free(filt);
}

static int hts_filter_eval_(hts_filter_t *filt,
                            void *data, hts_expr_sym_func *fn,
                            hts_expr_val_t *res) {
    char *end = NULL;

    filt->curr_regex = 0;
    if (expression(filt, data, fn, filt->str, &end, res))
        return -1;

    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", filt->str);
        return -1;
    }

    // Strings evaluate to true.  An empty string is also true, but an
    // absent (null) string is false, unless overriden by is_true.  An
    // empty string has kstring length of zero, but a pointer as it's
    // nul-terminated.
    if (res->is_str) {
        res->is_true |= res->s.s != NULL;
        res->d = res->is_true;
    } else if (hts_expr_val_exists(res)) {
        res->is_true |= res->d != 0;
    }

    return 0;
}

int hts_filter_eval(hts_filter_t *filt,
                    void *data, hts_expr_sym_func *fn,
                    hts_expr_val_t *res) {
    if (res->s.l != 0 || res->s.m != 0 || res->s.s != NULL) {
        // As *res is cleared below, it's not safe to call this function
        // with res->s.s set, as memory would be leaked.  It's also not
        // possible to know is res was initialised correctly, so in
        // either case we fail.
        hts_log_error("Results structure must be cleared before calling this function");
        return -1;
    }

    memset(res, 0, sizeof(*res));

    return hts_filter_eval_(filt, data, fn, res);
}

int hts_filter_eval2(hts_filter_t *filt,
                     void *data, hts_expr_sym_func *fn,
                     hts_expr_val_t *res) {
    ks_free(&res->s);
    memset(res, 0, sizeof(*res));

    return hts_filter_eval_(filt, data, fn, res);
}
