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

#define FIRST_SD_LISTEN_FD 3

typedef enum DaemonType {
    not_a_daemon = 0,
    sysv_daemon = 1,
    systemd_socket_service = 2
} DaemonType;

typedef struct MatchAddr {
    sa_family_t family;
    uint8_t mask_bytes;
    uint8_t mask;
    unsigned char addr[16];
} MatchAddr;

struct Options {
    const char *cache_dir;          // Directory for cached reference files
    const char *log_dir;            // Directory for log files
    const char *error_log_file;     // Error log file name
    FILE       *log;
    const char *upstream_url;       // URL for upstream server
    size_t      upstream_url_len;   // Cached length of upstream_url
    MatchAddr  *match_addrs;        // Client CIDR network allow list
    size_t      num_match_addrs;    // Number of items in match_addrs
    size_t      match_addrs_size;   // match_addrs allocation
    size_t      first_ip6;          // First ip6 range in sorted match_addrs
    off_t       max_log_sz;         // Size limit for log files in bytes
    int         cache_fd;           // File descriptor for cache_dir
    int         listen_fds;         // Number of fds passed in by systemd
    DaemonType  daemon;             // Run as a daemon or socket service
    uint16_t    port;               // Port to listen on
    uint16_t    nlogs;              // Number of log files to keep
    uint16_t    max_kids;           // Number of server processes to run
    uint8_t     verbosity;          // Print debugging messages
    uint8_t     no_log;             // Turn off logging
};

#endif /* OPTIONS_H_INCLUDED */
