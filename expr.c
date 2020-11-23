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

// TODO:
// - add maths functions.  pow, sqrt, log, min, max, ?
// - ?: operator for conditionals?

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <regex.h> // may need configure rule for this

#include "expr.h"
#include "textutils_internal.h"

// Could also cache fexpr_t stack here for kstring reuse?
#define MAX_REGEX 10
struct sam_filter_t {
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
 * Bit-wise:      &, |, ^  [NB as 3 precedence levels, in that order]
 * Conditionals:  >, >=, <, <=,
 * Equality:      ==, !=, =~ !~
 * Boolean:       &&, ||
 */

// Skip to start of term
static char *ws(char *str) {
    while (*str && (*str == ' ' || *str == '\t'))
        str++;
    return str;
}

static int expression(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res);

/*
 * simple_expr
 *     : identifier
 *     | constant
 * //  | string ?
 *     | '(' expression ')'
*/
static int simple_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
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
        } else if (fn)
            // Look up variable.
            return fn(data, str, end, res);
        else
            return -1;
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
static int unary_expr(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    int err;
    str = ws(str);
    if (*str == '+') {
        err = simple_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->is_true = res->d != 0;
    } else if (*str == '-') {
        err = simple_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->d = -res->d;
        res->is_true = res->d != 0;
    } else if (*str == '!') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        if (res->is_str) {
            res->is_str = 0;
            res->d = 0;
            res->is_true = !res->is_true;
        } else {
            res->d = !(int64_t)res->d;
            res->is_true = res->d != 0;
        }
    } else if (*str == '~') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->d = ~(int64_t)res->d;
        res->is_true = res->d != 0;
    } else {
        err = simple_expr(filt, data, fn, str, end, res);
    }
    return err ? -1 : 0;
}


/*
 * mul_expr
 *     : unary_expr (
 *           unary_expr '*' unary_expr
 *         | unary_expr '/' unary_expr
 *         | unary_expr '%' unary_expr
 *       )*
 */
static int mul_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (unary_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    fexpr_t val = FEXPR_INIT;
    while (*str) {
        str = ws(str);
        if (*str == '*' || *str == '/' || *str == '%') {
            if (unary_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (val.is_str || res->is_str) {
                fexpr_free(&val);
                return -1; // arith on strings
            }
        }

        if (*str == '*')
            res->d *= val.d;
        else if (*str == '/')
            res->d /= val.d;
        else if (*str == '%')
            res->d = (int64_t)res->d % (int64_t)val.d;
        else
            break;

        str = *end;
    }
    fexpr_free(&val);

    return 0;
}

/*
 * add_expr
 *     : mul_expr (
 *           mul_expr '+' mul_expr
 *         | mul_expr '-' mul_expr
 *       )*
 */
static int add_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (mul_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    fexpr_t val = FEXPR_INIT;
    while (*str) {
        str = ws(str);
        if (*str == '+' || *str == '-') {
            if (mul_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (val.is_str || res->is_str) {
                fexpr_free(&val);
                return -1; // arith on strings
            }
        }

        if (*str == '+')
            res->d += val.d;
        else if (*str == '-')
            res->d -= val.d;
        else
            break;

        str = *end;
    }
    fexpr_free(&val);

    return 0;
}

/*
 * bitand_expr
 *     : add_expr
 *     | bitand_expr '&' add_expr
 */
static int bitand_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
    if (add_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '&' && str[1] != '&') {
            if (add_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d & (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

/*
 * bitxor_expr
 *     : bitand_expr
 *     | bitxor_expr '^' bitand_expr
 */
static int bitxor_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
    if (bitand_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '^') {
            if (bitand_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d ^ (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

/*
 * bitor_expr
 *     : xor_expr
 *     | bitor_expr '|' xor_expr
 */
static int bitor_expr(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    if (bitxor_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '|' && str[1] != '|') {
            if (bitxor_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d | (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

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
static int cmp_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (bitor_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);
    fexpr_t val = FEXPR_INIT;
    int err = 0;

    if (*str == '>' && str[1] == '=') {
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) >= 0
            : !res->is_str && !val.is_str && res->d >= val.d;
        res->is_str = 0;
    } else if (*str == '>') {
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) > 0
            : !res->is_str && !val.is_str && res->d > val.d;
        res->is_str = 0;
    } else if (*str == '<' && str[1] == '=') {
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) <= 0
            : !res->is_str && !val.is_str && res->d <= val.d;
        res->is_str = 0;
    } else if (*str == '<') {
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) < 0
            : !res->is_str && !val.is_str && res->d < val.d;
        res->is_str = 0;
    }
    fexpr_free(&val);

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
static int eq_expr(sam_filter_t *filt, void *data, sym_func *fn,
                   char *str, char **end, fexpr_t *res) {
    if (cmp_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);

    int err = 0;
    fexpr_t val = FEXPR_INIT;

    // numeric vs numeric comparison is as expected
    // string vs string comparison is as expected
    // numeric vs string is false
    if (str[0] == '=' && str[1] == '=') {
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            res->is_true = res->d = res->is_str
                ? (res->s.s && val.s.s ? strcmp(res->s.s, val.s.s)==0 : 0)
                : !res->is_str && !val.is_str && res->d == val.d;
        }
        res->is_str = 0;

    } else if (str[0] == '!' && str[1] == '=') {
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            res->is_true = res->d = res->is_str
                ? (res->s.s && val.s.s ? strcmp(res->s.s, val.s.s) != 0 : 1)
                : res->is_str != val.is_str || res->d != val.d;
        }
        res->is_str = 0;

    } else if ((str[0] == '=' && str[1] == '~') ||
               (str[0] == '!' && str[1] == '~')) {
        err = eq_expr(filt, data, fn, str+2, end, &val);
        if (!val.is_str || !res->is_str) {
            fexpr_free(&val);
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
                    fexpr_free(&val);
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
    fexpr_free(&val);

    return err ? -1 : 0;
}

/*
 * and_expr
 *     : eq_expr
 *     | and_expr 'and' eq_expr
 *     | and_expr 'or'  eq_expr
 */
static int and_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (eq_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (str[0] == '&' && str[1] == '&') {
            if (eq_expr(filt, data, fn, str+2, end, &val)) return -1;
            res->is_true = res->d =
                (res->is_true || (res->is_str && res->s.s) || res->d) &&
                (val.is_true  || (val.is_str && val.s.s) || val.d);
            res->is_str = 0;
        } else if (str[0] == '|' && str[1] == '|') {
            if (eq_expr(filt, data, fn, str+2, end, &val)) return -1;
            res->is_true = res->d =
                res->is_true || (res->is_str && res->s.s) || res->d ||
                val.is_true  || (val.is_str  && val.s.s ) || val.d;
            res->is_str = 0;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

static int expression(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    return and_expr(filt, data, fn, str, end, res);
}

sam_filter_t *sam_filter_init(const char *str) {
    sam_filter_t *f = calloc(1, sizeof(*f));
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

void sam_filter_free(sam_filter_t *filt) {
    if (!filt)
        return;

    int i;
    for (i = 0; i < filt->max_regex; i++)
        regfree(&filt->preg[i]);

    free(filt->str);
    free(filt);
}

int sam_filter_eval(sam_filter_t *filt, void *data, sym_func *fn,
                    fexpr_t *res) {
    char *end = NULL;

    memset(res, 0, sizeof(*res));

    filt->curr_regex = 0;
    if (expression(filt, data, fn, filt->str, &end, res))
        return -1;

    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", filt->str);
        return -1;
    }

    // Strings evaluate to true.  An empty string is also true, but an
    // absent (null) string is false.  An empty string has kstring length
    // of zero, but a pointer as it's nul-terminated.
    if (res->is_str)
        res->is_true = res->d = res->s.s != NULL;
    else
        res->is_true |= res->d != 0;

    return 0;
}
