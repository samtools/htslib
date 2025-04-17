/*  upstream.h -- download ref-cache files from upstream host

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

#ifndef UPSTREAM_H_INCLUDED
#define UPSTREAM_H_INCLUDED

#include <stdint.h>

#include "types.h"

typedef enum {
  US_START,
  US_CONTENT_LENGTH,
  US_PARTIAL_LENGTH,
  US_RESULT
} Upstream_msg_code;

typedef struct {
  unsigned int        id;
  Upstream_msg_code code;
  int64_t           val;
} Upstream_msg;

int upstream_send_cmd(int cmd_fd, const char hexmd5[32], unsigned int id);
int upstream_recv_msg(int cmd_fd, Upstream_msg *msg, int *fd);
int run_upstream_handler(Options *opts, int *sockets, int liveness_fd);
#endif /* UPSTREAM_H_INCLUDED */
