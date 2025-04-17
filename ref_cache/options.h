/*  options.h -- ref-cache command-line options

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


#ifndef OPTIONS_H_INCLUDED
#define OPTIONS_H_INCLUDED

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "types.h"

typedef struct MatchAddr {
    sa_family_t family;
    uint8_t mask_bytes;
    uint8_t mask;
    unsigned char addr[16];
} MatchAddr;

struct Options {
    const char *cache_dir;
    const char *log_dir;
    const char *error_log_file;
    FILE       *log;
    const char *upstream_url;
    size_t      upstream_url_len;
    MatchAddr  *match_addrs;
    size_t      num_match_addrs;
    size_t      match_addrs_size;
    size_t      first_ip6;
    off_t       max_log_sz;
    int         cache_fd;
    int         log_dir_fd;
    uint16_t    port;
    uint16_t    nlogs;
    uint16_t    max_kids;
    uint8_t     verbosity;
    uint8_t     daemon;
    uint8_t     no_log;
};

#endif /* OPTIONS_H_INCLUDED */
