/*  textutils.c -- non-bioinformatics utility routines for text etc.

    Copyright (C) 2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#include "hts_internal.h"

static int dehex(char c)
{
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    else if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    else if (c >= '0' && c <= '9') return c - '0';
    else return -1;  // Hence dehex('\0') = -1
}

int hts_decode_percent(char *dest, size_t *destlen, const char *s)
{
    char *d = dest;
    int hi, lo;

    while (*s) {
        if (*s == '%' && (hi = dehex(s[1])) >= 0 && (lo = dehex(s[2])) >= 0) {
            *d++ = (hi << 4) | lo;
            s += 3;
        }
        else *d++ = *s++;
    }

    *d = '\0';
    *destlen = d - dest;
    return 0;
}

static int debase64(char c)
{
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    else if (c >= 'A' && c <= 'Z') return c - 'A';
    else if (c >= '0' && c <= '9') return c - '0' + 52;
    else if (c == '/') return 63;
    else if (c == '+') return 62;
    else return -1;  // Hence debase64('\0') = -1
}

size_t hts_base64_decoded_length(size_t len)
{
    size_t nquartets = (len + 2) / 4;
    return 3 * nquartets;
}

int hts_decode_base64(char *dest, size_t *destlen, const char *s)
{
    char *d = dest;
    int x0, x1, x2, x3;

    while (1) {
        x0 = debase64(*s++);
        x1 = (x0 >= 0)? debase64(*s++) : -1;
        x2 = (x1 >= 0)? debase64(*s++) : -1;
        x3 = (x2 >= 0)? debase64(*s++) : -1;
        if (x3 < 0) break;

        *d++ = (x0 << 2) | (x1 >> 4);
        *d++ = (x1 << 4) | (x2 >> 2);
        *d++ = (x2 << 6) | x3;
    }

    if (x1 >= 0) *d++ = (x0 << 2) | (x1 >> 4);
    if (x2 >= 0) *d++ = (x1 << 4) | (x2 >> 2);

    *destlen = d - dest;
    return 0;
}
