/*  log_files.h -- ref-cache log files interface

    Copyright (C) 2025 Genome Research Ltd.

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


#ifndef LOG_FILES_H_INCLUDED
#define LOG_FILES_H_INCLUDED

#include <stddef.h>
#include "../htslib/hts_defs.h"
#include "types.h"

int redirect_stdio(const Logfiles *logfiles, const Options *opts)
HTS_ACCESS(read_only, 1)
HTS_ACCESS(read_only, 2)
HTS_RESULT_USED;

void close_logs(Logfiles *logfiles);

Logfiles * open_logs(const Options *opts)
HTS_ACCESS(read_only, 1);

int write_to_log(Logfiles *logfiles, const Options *opts,
                 const char *msg, size_t bytes)
HTS_ACCESS(read_only, 2)
HTS_ACCESS(read_only, 3, 4)
HTS_RESULT_USED;

#endif
