/*  test-expr.c -- Testing: filter expression parsing and processing.

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

#include <config.h>

#include <stdio.h>
#include <string.h>
#include "../htslib/hts_expr.h"

int lookup(void *data, char *str, char **end, hts_expr_val_t *res) {
    int foo = 15551; // my favourite palindromic prime
    int a = 1;
    int b = 2;
    int c = 3;
    res->is_str = 0;
    if (strncmp(str, "foo", 3) == 0) {
        *end = str+3;
        res->d = foo;
    } else if (*str == 'a') {
        *end = str+1;
        res->d = a;
    } else if (*str == 'b') {
        *end = str+1;
        res->d = b;
    } else if (*str == 'c') {
        *end = str+1;
        res->d = c;
    } else if (strncmp(str, "magic", 5) == 0) {
        // non-empty string
        *end = str+5;
        res->is_str = 1;
        kputs("plugh", ks_clear(&res->s));
    } else if (strncmp(str, "empty", 5) == 0) {
        // empty string
        *end = str+5;
        res->is_str = 1;
        kputs("", ks_clear(&res->s));
    } else if (strncmp(str, "null", 4) == 0) {
        // null string (eg aux:Z tag is absent)
        *end = str+4;
        res->is_str = 1;
        ks_clear(&res->s);

    } else {
        return -1;
    }

    return 0;
}

typedef struct {
    double dval;
    char *sval;
    char *str;
} test_ev;

int test(void) {
    // These are all valid expressions that should work
    test_ev tests[] = {
        {  1, NULL, "1"},
        {  1, NULL, "+1"},
        { -1, NULL, "-1"},
        {  0, NULL, "!7"},
        {  1, NULL, "!0"},
        {  1, NULL, "!(!7)"},
        {  1, NULL, "!!7"},

        {  5, NULL, "2+3"},
        { -1, NULL, "2+-3"},
        {  6, NULL, "1+2+3"},
        {  1, NULL, "-2+3"},

        {  6, NULL, "2*3"},
        {  6, NULL, "1*2*3"},
        {  0, NULL, "2*0"},

        {  7, NULL, "(7)"},
        {  7, NULL, "((7))"},
        { 21, NULL, "(1+2)*(3+4)"},
        { 14, NULL, "(4*5)-(-2*-3)"},

        {  1, NULL, "(1+2)*3==9"},
        {  1, NULL, "(1+2)*3!=8"},
        {  0, NULL, "(1+2)*3!=9"},
        {  0, NULL, "(1+2)*3==8"},

        {  0, NULL, "1>2"},
        {  1, NULL, "1<2"},
        {  0, NULL, "3<3"},
        {  0, NULL, "3>3"},
        {  1, NULL, "9<=9"},
        {  1, NULL, "9>=9"},
        {  1, NULL, "2*4==8"},
        {  1, NULL, "16==0x10"},
        {  1, NULL, "15<0x10"},
        {  1, NULL, "17>0x10"},
        {  0, NULL, "2*4!=8"},
        {  1, NULL, "4+2<3+4"},
        {  0, NULL, "4*2<3+4"},
        {  8, NULL, "4*(2<3)+4"},  // boolean; 4*(1)+4

        {  1, NULL, "(1<2) == (3>2)"},
        {  1, NULL, "1<2 == 3>2"},

        {  1, NULL, "2 && 1"},
        {  0, NULL, "2 && 0"},
        {  0, NULL, "0 && 2"},
        {  1, NULL, "2 || 1"},
        {  1, NULL, "2 || 0"},
        {  1, NULL, "0 || 2"},
        {  1, NULL, "1 || 2 && 3"},
        {  1, NULL, "2 && 3 || 1"},
        {  1, NULL, "0 && 3 || 2"},
        {  0, NULL, "0 && 3 || 0"},

        {  1, NULL, "3 & 1"},
        {  2, NULL, "3 & 2"},
        {  3, NULL, "1 | 2"},
        {  3, NULL, "1 | 3"},
        {  7, NULL, "1 | 6"},
        {  2, NULL, "1 ^ 3"},

        {  1, NULL, "(1^0)&(4^3)"},
        {  2, NULL, "1 ^(0&4)^ 3"},
        {  2, NULL, "1 ^ 0&4 ^ 3"},  // precedence, & before ^

        {  6, NULL, "(1|0)^(4|3)"},
        {  7, NULL, "1 |(0^4)| 3"},
        {  7, NULL, "1 | 0^4 | 3"},  // precedence, ^ before |

        {  1, NULL, "4 & 2 || 1"},
        {  1, NULL, "(4 & 2) || 1"},
        {  0, NULL, "4 & (2 || 1)"},
        {  1, NULL, "1 || 4 & 2"},
        {  1, NULL, "1 || (4 & 2)"},
        {  0, NULL, "(1 || 4) & 2"},

        {  1, NULL, " (2*3)&7  > 4"},
        {  0, NULL, " (2*3)&(7 > 4)"}, // C precedence equiv
        {  1, NULL, "((2*3)&7) > 4"},  // Python precedence equiv
        {  1, NULL, "((2*3)&7) > 4 && 2*2 <= 4"},

        {  1, "plugh", "magic"},
        {  1, "",   "empty"},
        {  1, NULL, "magic == \"plugh\""},
        {  1, NULL, "magic != \"xyzzy\""},

        {  1, NULL, "\"abc\" < \"def\""},
        {  1, NULL, "\"abc\" <= \"abc\""},
        {  0, NULL, "\"abc\" < \"ab\""},
        {  0, NULL, "\"abc\" <= \"ab\""},

        {  0, NULL, "\"abc\" > \"def\""},
        {  1, NULL, "\"abc\" >= \"abc\""},
        {  1, NULL, "\"abc\" > \"ab\""},
        {  1, NULL, "\"abc\" >= \"ab\""},

        {  1, NULL, "\"abbc\" =~ \"^a+b+c+$\""},
        {  0, NULL, "\"aBBc\" =~ \"^a+b+c+$\""},
        {  1, NULL, "\"aBBc\" !~ \"^a+b+c+$\""},
        {  1, NULL, "\"xyzzy plugh abracadabra\" =~ magic"},
    };

    int i;
    hts_expr_val_t r;
    for (i = 0; i < sizeof(tests) / sizeof(*tests); i++) {
        hts_filter_t *filt = hts_filter_init(tests[i].str);
        if (!filt)
            return 1;
        if (hts_filter_eval(filt, NULL, lookup, &r)) {
            fprintf(stderr, "Failed to parse filter string %s\n",
                    tests[i].str);
            return 1;
        }

        if (r.is_str && (strcmp(r.s.s, tests[i].sval) != 0
                         || r.d != tests[i].dval)) {
            fprintf(stderr, "Failed test: %s == %s, got %s, %f\n",
                    tests[i].str, tests[i].sval, r.s.s, r.d);
            return 1;
        } else if (!r.is_str && r.d != tests[i].dval) {
            fprintf(stderr, "Failed test: %s == %f, got %f\n",
                    tests[i].str, tests[i].dval, r.d);
            return 1;
        }

        hts_expr_val_free(&r);
        hts_filter_free(filt);
    }

    return 0;
}

int main(int argc, char **argv) {
    if (argc > 1) {
        hts_expr_val_t v;
        hts_filter_t *filt = hts_filter_init(argv[1]);
        if (hts_filter_eval(filt, NULL, lookup, &v))
            return 1;

        if (v.is_str)
            puts(v.s.s);
        else
            printf("%g\n", v.d);

        hts_expr_val_free(&v);
        hts_filter_free(filt);
        return 0;
    }

    return test();
}
