/*  fuzz_settings.h -- fuzz-tester specific definitions

    Copyright (C) 2023 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

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

#ifndef HTSLIB_FUZZ_SETTINGS_H
#define HTSLIB_FUZZ_SETTINGS_H
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION

#ifndef FUZZ_ALLOC_LIMIT
// By default libfuzzer reports out-of-memory on allocations > 2 Gbytes
#define FUZZ_ALLOC_LIMIT 2000000000ULL
#endif

#endif
#endif
