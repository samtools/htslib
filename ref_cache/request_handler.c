/*  request_handler.c -- ref-cache http requests

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

#include <config.h>

#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "request_handler.h"
#include "http_parser.h"
#include "misc.h"
#include "options.h"
#include "ref_files.h"
#include "transaction.h"
#include "upstream.h"

static const char *text_plain = "text/plain";
/* static const char *text_html  = "text/html"; */

static inline int is_hexmd5(char *str) {
    int i;
    for (i = 0; i < 32; i++)
        if (hexval(str[i]) == -1) return 0;
    return hexval(str[32]) == -1 ? 1 : 0;
}

static char * decode_uri(Http_Parser *parser) {
    char *uri = parser->uri, *querypart, *out, last;
    const char *in;
    int d1, d2;
    if (uri == NULL) return NULL;

    /* Deal with absolute URLs */
    if (strncasecmp(uri, "http://", 7) == 0) {
        char *localpart = strchr(uri + 7, '/');
        if (localpart == NULL) return NULL;
        uri = localpart;
    }

    /* Should always start with / now */
    if (*uri != '/') return NULL;

    /* Hack off query part */
    querypart = strchr(uri, '?');
    if (querypart != NULL) *querypart = '\0';

    /* Deal with multiple // and % decoding. URI will always shrink. */
    last = '\0';
    for (in = out = uri; *in; in++) {
        if (*in == '/' && last == '/') continue;
        if (*in == '%'
            && ((d1 = hexval(*(in + 1))) >= 0)
            && ((d2 = hexval(*(in + 2))) >= 0)) {
            char v = (char) (d1 << 4 | d2);
            in += 2;
            if (v == '/' && last == '/') continue;
            *out++ = last = v; continue;
        }
        *out++ = last = *in;
    }
    *out = '\0';
    return uri;
}

static void handle_hello(Transaction *transact) {
    const char *resp = "Hello\r\n";
    set_message_response(transact, text_plain, resp, strlen(resp));
}

static void handle_md5(const Options *opts, Http_Parser *parser,
                       Transaction *transact, char *md5) {
    // off_t range_start = -1, range_end = 0;
    // int have_range = 0;
    RefFile *ref_file = get_ref_file(opts, md5, parser->upstream);
    RefFileStatus status;
    off_t size;

    if (ref_file == NULL) {
        set_error_response(transact, ERR_INTERNAL);
        return;
    }

    status = get_ref_status(ref_file);
    if (status == REF_NOT_FOUND) {
        set_error_response(transact, ERR_NOT_FOUND);
        release_ref_file(ref_file);
        return;
    }

    size = get_ref_size(ref_file);

    transaction_set_ref(transact, ref_file);

    set_transaction_file_range(transact, size, status >= REF_DOWNLOAD_STARTED);
}

static void handle_get(const Options *opts, Http_Parser *parser,
                       Transaction *transact) {
    char *requested = decode_uri(parser);

    if (opts->verbosity > 1) {
        fprintf(stderr, "Request: GET %s\n", requested ? requested : "");
    }

    transaction_set_req_str(transact, requested ? requested : "");

    if (requested == NULL || *requested != '/') {
        set_error_response(transact, ERR_BAD_REQUEST); return;
    }
    requested++;
    if (is_hexmd5(requested) && requested[32] == '\0') {
        handle_md5(opts, parser, transact, requested); return;
    }
    if (strcmp(requested, "hello") == 0) { handle_hello(transact); return; }
    set_error_response(transact, ERR_NOT_FOUND);
    return;
}

void handle_request(const Options *opts, Client *client,
                    Http_Parser *parser, Transaction **transact_out) {
    Transaction *transact = new_transaction(client, parser);
    if (transact == NULL) { parser->state = ERR_INTERNAL; return; }

    switch (parser->req_type) {
    case REQ_GET:  handle_get(opts, parser, transact);  break;
        /* case REQ_HEAD: handle_head(parser, transact); break; */
    default:       set_error_response(transact, ERR_UNIMPLEMENTED);
    }
    *transact_out = transact;

    if ((parser->flags & REQ_KEEP_ALIVE) != 0 && parser->http_vers == HTTP_1_1) {
        parser->state = READING_REQUEST_LINE;
    } else {
        parser->state = SHUTTING_DOWN;
    }
}

void handle_error(Client *client, Http_Parser *parser, int code, Transaction **transact_out) {
    Transaction *transact = new_transaction(client, parser);
    if (transact == NULL) { parser->state = ERR_INTERNAL; return; }
    set_error_response(transact,
                       code >= ERR_BAD_REQUEST ? (unsigned) code : 500U);
    *transact_out = transact;
}
