/*  test_kfunc.c -- kt_fisher_exact() unit tests

    Copyright (C) 2020 University of Glasgow.

    Author: John Marshall <John.W.Marshall@glasgow.ac.uk>

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../htslib/kfunc.h"

int differ(double obs, double expected)
{
    return fabs(obs - expected) > 1e-8;
}

int nfailed = 0;

void fail(const char *test, double obs, double expected,
          int n11, int n12, int n21, int n22)
{
    fprintf(stderr, "[%d %d | %d %d] %s: %g (expected %g)\n",
            n11, n12, n21, n22, test, obs, expected);
    nfailed++;
}

void test_fisher(int n11, int n12, int n21, int n22,
                 double eleft, double eright, double etwo, double eprob)
{
    double prob, left, right, two;
    prob = kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two);
    if (differ(left, eleft)) fail("LEFT", left, eleft, n11, n12, n21, n22);
    if (differ(right, eright)) fail("RIGHT", right, eright, n11, n12, n21, n22);
    if (differ(two, etwo)) fail("TWO-TAIL", two, etwo, n11, n12, n21, n22);
    if (differ(prob, eprob)) fail("RESULT", prob, eprob, n11, n12, n21, n22);
}

int main(int argc, char **argv)
{
    test_fisher(2, 1, 0, 31,   1.0, 0.005347593583, 0.005347593583, 0.005347593583);
    test_fisher(2, 1, 0, 1,    1.0, 0.5, 1.0, 0.5);
    test_fisher(3, 1, 0, 0,    1.0, 1.0, 1.0, 1.0);
    test_fisher(3, 15, 37, 45, 0.021479750169, 0.995659202564, 0.033161943699, 0.017138952733);
    test_fisher(12, 5, 29, 2,  0.044554737835, 0.994525206022, 0.080268552074, 0.039079943857);

    test_fisher(781, 23171, 4963, 2455001,   1.0, 0.0, 0.0, 0.0);
    test_fisher(333, 381, 801722, 7664285,   1.0, 0.0, 0.0, 0.0);
    test_fisher(4155, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0);
    test_fisher(4455, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0);
    test_fisher(5455, 4903, 805463, 8507517, 1.0, 0.0, 0.0, 0.0);

    test_fisher(1, 1, 100000, 1000000, 0.991735477166, 0.173555146661, 0.173555146661, 0.165290623827);
    test_fisher(1000, 1000, 100000, 1000000, 1.0, 0.0, 0.0, 0.0);
    test_fisher(1000, 1000, 1000000, 100000, 0.0, 1.0, 0.0, 0.0);

    test_fisher(49999, 10001,  90001, 49999, 1.0, 0.0, 0.0, 0.0);
    test_fisher(50000, 10000,  90000, 50000, 1.0, 0.0, 0.0, 0.0);
    test_fisher(50001,  9999,  89999, 50001, 1.0, 0.0, 0.0, 0.0);
    test_fisher(10000, 50000, 130000, 10000, 0.0, 1.0, 0.0, 0.0);

    if (nfailed > 0) {
        const char *plural = (nfailed == 1)? "" : "s";
        fprintf(stderr, "Failed %d test case%s\n", nfailed, plural);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
