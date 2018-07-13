/// @file ref.h
/// Reference genome fetching
/*
    Copyright (C) 2017-2018 Genome Research Ltd.

    Author: Thomas Hickman <th10@sanger.ac.uk>

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


#ifndef HTSLIB_REF_H
#define HTSLIB_REF_H

#include "htslib/hfile.h"

#ifdef __cplusplus
extern "C" {
#endif

/* m5_to_ref() - returns the reference genome that has a given MD5 string
 * @param m5_str: The m5 string to query
 * @returns: A hFILE containing a file pointer to a reference genome.
 *           NULL on failure
 *
 * Note: This function is not currently thread safe, so locks
 * need to be acquired before calling this, in a multi-threaded
 * enviroment
 */
hFILE* m5_to_ref(const char *m5_str);

/* m5_to_path() - returns a path to a reference genome that has a
 * given MD5 string
 * @param m5_str: The m5 string to query
 * @returns: A path to the correct reference genome.
 *           NULL on failure
 *
 * Note: This function is not currently thread safe, so locks
 * need to be acquired before calling this, in a multi-threaded
 * enviroment
 */
char* m5_to_path(const char *m5_str);

#ifdef __cplusplus
}
#endif

#endif
