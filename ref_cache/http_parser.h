/*  http_parser.h -- ref-cache HTTP protocol handler interface

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

#ifndef HTTP_PARSER_H_INCLUDED
#define HTTP_PARSER_H_INCLUDED

#include <stdint.h>
#include "../htslib/hts_defs.h"
#include "types.h"

/* Maximum lengths for log file */
#define MAX_UA_LEN       128
#define MAX_REFERRER_LEN 128
#define MAX_REQUEST_LEN  64

typedef enum {
  READING_REQUEST_LINE = 0,
  READING_HEADERS,
  READING_CHUNK_HEADER,
  READING_CHUNK,
  READING_CHUNK_TRAILER,
  READING_BODY,
  GOT_REQUEST,
  SHUTTING_DOWN,
  ERR_BAD_REQUEST   = 400,
  ERR_NOT_FOUND     = 404,
  ERR_TOO_LARGE     = 413,
  ERR_LONG_URI      = 414,
  ERR_INTERNAL      = 500,
  ERR_UNIMPLEMENTED = 501,
  ERR_HTTP_VERS     = 505,
} Parse_state;

typedef enum {
  REQ_OPTIONS = 0,
  REQ_GET,
  REQ_HEAD,
  REQ_POST,
  REQ_PUT,
  REQ_DELETE,
  REQ_TRACE,
  REQ_CONNECT,
  REQ_OTHER
} Request_type;

typedef enum {
  TE_IDENT = 0,
  TE_CHUNKED = 1,
  TE_OTHER = 0x80
} Transfer_encoding;

static const unsigned int REQ_EOF          =  1U;
static const unsigned int REQ_KEEP_ALIVE   =  2U;
static const unsigned int REQ_DEBUG        =  4U;
static const unsigned int REQ_RANGE_FROM   =  8U;
static const unsigned int REQ_RANGE_TO     = 16U;
static const unsigned int REQ_RANGE_SUFFIX = 32U;

typedef enum {
  HTTP_0_9 = 0,
  HTTP_1_0,
  HTTP_1_1,
  HTTP_OTHER
} Http_version;

struct Http_Parser {
  Parse_state       state;
  Request_type      req_type;
  Http_version      http_vers;
  Transfer_encoding trans_enc;
  unsigned long     content_length;
  unsigned long     bytes;
  char             *uri;
  char             *key;
  char             *val;
  char             *buffer;
  char             *user_agent;
  char             *referrer;
  off_t             range_from;
  off_t             range_to;
  size_t            key_sz;
  size_t            key_used;
  size_t            val_sz;
  size_t            val_used;
  int               upstream;
  unsigned int      flags;
  unsigned int      in;
  unsigned int      out;
  unsigned int      pos;
  unsigned int      used;
};

typedef enum {
  READ_BLOCKED,
  READ_MORE,
  READ_EOF,
  READ_ERROR
} Read_result;

int init_http_parser(Http_Parser *parser, int upstream);

void cleanup_http_parser(Http_Parser *parser);

Read_result parser_read_data(const Options *opts, Client *client,
                             Http_Parser *parser, int fd)
HTS_ACCESS(read_only, 1);

char * steal_user_agent_from_parser(Http_Parser *parser);

char * steal_referrer_from_parser(Http_Parser *parser);

#endif
