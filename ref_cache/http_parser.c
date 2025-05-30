/*  http_parser.c -- ref-cache HTTP protocol handler

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <unistd.h>
#include <errno.h>
#include <sys/uio.h>

#include "http_parser.h"
#include "misc.h"
#include "options.h"
#include "request_handler.h"
#include "server.h"
#include "../cram/pooled_alloc.h"

#if defined(REF_CACHE_NO_POOLED_ALLOC)
#  define pool_alloc(x) malloc(x->dsize)
#  define pool_free(p, x) free(x)
#endif

#define BUF_SZ   0x400U
#define BUF_MASK (BUF_SZ - 1)

#define OFF_MAX ((off_t) ((sizeof(off_t) < 8) ? INT32_MAX : INT64_MAX))

static const unsigned char lws_chars[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
#define is_lws(X) (lws_chars[(unsigned char) (X)])

static const unsigned char token_chars[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
#define is_token(X) (token_chars[(unsigned char) (X)])

static inline size_t lws_spn(char *s) {
    size_t c = 0;
    while (is_lws(s[c])) c++;
    return c;
}

static inline size_t lws_cspn(char *s) {
    size_t c = 0;
    while (s[c] && !is_lws(s[c])) c++;
    return c;
}

static inline size_t tok_spn(char *s) {
    size_t c = 0;
    while (is_token(s[c])) c++;
    return c;
}

int init_http_parser(Http_Parser *parser, int upstream) {

    memset(parser, 0, sizeof(*parser));

    parser->uri = NULL;
    parser->key = NULL;
    parser->val = NULL;
    parser->user_agent = NULL;
    parser->referrer = NULL;

    parser->buffer = malloc(BUF_SZ);
    if (parser->buffer == NULL) {
        return -1;
    }

    parser->upstream = upstream;

    return 0;
}

void cleanup_http_parser(Http_Parser *parser) {
    if (parser->uri) {
        free(parser->uri);
        parser->uri = NULL;
    }
    if (parser->key) {
        free(parser->key);
        parser->key = NULL;
    }
    if (parser->val) {
        free(parser->val);
        parser->val = NULL;
    }
    if (parser->user_agent) {
        free(parser->user_agent);
        parser->user_agent = NULL;
    }
    if (parser->referrer) {
        free(parser->referrer);
        parser->referrer = NULL;
    }
    if (parser->buffer) {
        free(parser->buffer);
        parser->buffer = NULL;
    }
}

static char * parser_get_line(Http_Parser *parser, size_t *len) {
    static char line[BUF_SZ + 1];
    size_t lpos;

    if (parser->used == 0) {
        assert(parser->in == parser->out);
        assert(parser->pos == parser->out);
        return NULL;
    }

    /* Search for end of line */
    do {
        if (parser->buffer[parser->pos] == '\x0a') break;
        parser->pos = (parser->pos + 1) & BUF_MASK;
    } while (parser->pos != parser->in);
    if (parser->pos == parser->in && parser->buffer[parser->pos] != '\x0a') {
        if (parser->used == BUF_SZ) {
            parser->state = ERR_TOO_LARGE;
            return NULL;
        }
        return NULL;
    }

    lpos = 0;
    if (parser->pos <= parser->out) {
        memcpy(&line[lpos], parser->buffer + parser->out, BUF_SZ - parser->out);
        lpos += BUF_SZ - parser->out;
        parser->out = 0;
    }
    memcpy(&line[lpos],
           parser->buffer + parser->out, parser->pos - parser->out);
    lpos += parser->pos - parser->out;
    line[lpos] = '\0';
    if (lpos > 0 && line[lpos - 1] == '\x0d') line[--lpos] = '\0';
    *len = lpos;
    parser->pos = (parser->pos + 1) & BUF_MASK; /* Get past \n */
    parser->out = parser->pos; /* Eat the data */
    parser->used = (parser->out > parser->in
                    ? BUF_SZ - parser->out + parser->in
                    : parser->in - parser->out);
    return line;
}

static Read_result parser_read_input(Http_Parser *parser, int fd) {
    struct iovec iov[2] = {
        { parser->buffer + parser->in, 0U },
        { parser->buffer, 0U }
    };
    ssize_t res;
    int nio = 1;

    /* Read into circular buffer */

    assert(parser->used <= BUF_SZ);
    if (parser->used == BUF_SZ) // Full
        return READ_MORE;

    if (parser->in > parser->out || parser->used == 0) {
        // Unused space from in to end
        iov[0].iov_len = BUF_SZ - parser->in;
        if (parser->out > 0) {
            // and start to out
            iov[1].iov_len = parser->out;
            nio = 2;
        }
    } else {
        // Unused space from in to out
        assert(parser->in < parser->out);
        iov[0].iov_len = parser->out - parser->in;
    }

    res = readv(fd, iov, nio);
    if (res < 0) {
        if (errno != EAGAIN && errno != EWOULDBLOCK && errno != EINTR) {
            fprintf(stderr, "Error from fd #%d : %s\n", fd, strerror(errno));
            return READ_ERROR;
        }
        return READ_BLOCKED;
    } else {
        parser->in = ((parser->in + (unsigned int) res) & BUF_MASK);
        parser->used += (unsigned int) res;
    }

    if (res == 0) return READ_EOF;
    return ((size_t) res < iov[0].iov_len + iov[1].iov_len
            ? READ_BLOCKED : READ_MORE);
}

static void read_request_line(const Options *opts, Http_Parser *parser) {
    char *line;
    size_t len, reqlen, uripos, urilen, verpos;

    do {
        line = parser_get_line(parser, &len);
    } while (line != NULL && len == 0);
    if (line == NULL) return;

    if (opts->verbosity > 2) fprintf(stderr, "RECV'D: %s\n", line);

    reqlen = lws_cspn(line);
    uripos = lws_spn(line + reqlen) + reqlen;
    urilen = lws_cspn(line + uripos);
    verpos = lws_spn(line + uripos + urilen) + uripos + urilen;

    if (reqlen == 0 || urilen == 0) { parser->state = ERR_BAD_REQUEST; return; }
    if (urilen > 128) { parser->state = ERR_LONG_URI; return; }

    /* Decode request type */

    switch (reqlen) {
    case 3:
        if (strncmp(line, "GET", 3) == 0)    { parser->req_type = REQ_GET; break; }
        if (strncmp(line, "PUT", 3) == 0)    { parser->req_type = REQ_PUT; break; }
        parser->req_type = REQ_OTHER; break;
    case 4:
        if (strncmp(line, "HEAD", 4) == 0)   { parser->req_type = REQ_HEAD; break; }
        if (strncmp(line, "POST", 4) == 0)   { parser->req_type = REQ_POST; break; }
        parser->req_type = REQ_OTHER; break;
    case 5:
        if (strncmp(line, "TRACE", 5) == 0)  { parser->req_type = REQ_TRACE; break; }
        parser->req_type = REQ_OTHER; break;
    case 6:
        if (strncmp(line, "DELETE", 6) == 0) { parser->req_type = REQ_DELETE; break; }
        parser->req_type = REQ_OTHER; break;
    case 7:
        if (strncmp(line, "OPTIONS", 7) == 0){ parser->req_type = REQ_OPTIONS; break; }
        if (strncmp(line, "CONNECT", 7) == 0){ parser->req_type = REQ_CONNECT; break; }
        // fall through
    default:
        parser->req_type = REQ_OTHER; break;
    }

    /* Copy URI */

    parser->uri = malloc(urilen + 1);
    if (parser->uri == NULL) { parser->state = ERR_INTERNAL; return; }
    memcpy(parser->uri, line + uripos, urilen);
    parser->uri[urilen] = '\0';

    /* Get HTTP version */

    if (line[verpos] == '\0') {
        parser->http_vers = HTTP_0_9;
    } else {
        long major, minor;
        char *p;
        if (strncmp(line + verpos, "HTTP/", 5) != 0) {
            parser->state = ERR_BAD_REQUEST;
            return;
        }
        major = strtol(line + verpos + 5, &p, 10);
        if (*p != '.') { parser->state = ERR_BAD_REQUEST; return; }
        minor = strtol(p + 1, NULL, 10);
        if (major == 0) {
            parser->http_vers = HTTP_0_9;
        } else if (major == 1) {
            parser->http_vers = minor == 0 ? HTTP_1_0 : HTTP_1_1;
        } else {
            parser->state = ERR_HTTP_VERS; return;
        }
    }

    if (parser->http_vers == HTTP_1_1) parser->flags |= REQ_KEEP_ALIVE;

    parser->state = READING_HEADERS;
}

/* parse_range output states:
   flags set                         desired range
   REQ_RANGE_FROM and REQ_RANGE_TO   range_from to range_to, inclusive
   REQ_RANGE_FROM only               range_from to end
   REQ_RANGE_SUFFIX only             last range_to bytes
   none                              everything
*/

static void parse_range(Http_Parser *parser) {
    char *v = parser->val, *end;
    long long int ll;

    /* bytes= */
    if (strncasecmp(v, "bytes", 5) != 0) goto fail;
    v += 5 + lws_spn(v + 5);
    if (*v != '=') goto fail;
    v += 1 + lws_spn(v + 1);

    /* From rfc2068:
       byte-range-set  = 1#( byte-range-spec | suffix-byte-range-spec )
       byte-range-spec = first-byte-pos "-" [last-byte-pos]
       first-byte-pos  = 1*DIGIT
       last-byte-pos   = 1*DIGIT
       suffix-byte-range-spec = "-" suffix-length
       suffix-length = 1*DIGIT */

    parser->range_from = OFF_MAX;
    parser->range_to   = 0;

    for (;;) {
        if (*v == '-') { /* suffix-byte-range-spec */
            ll = strtoll(v + 1, &end, 10);  /* suffix-length */
            if (ll < 0)
                goto fail;
            parser->flags |= REQ_RANGE_SUFFIX;
            parser->range_to = parser->range_to > ll ? parser->range_to : ll;
            v = end + lws_spn(end);
        } else { /* byte-range-spec */
            ll = strtoll(v, &end, 10);  /* first-byte-pos */
            if (ll < 0)
                goto fail;
            parser->flags |= REQ_RANGE_FROM;
            parser->range_from = parser->range_from < ll ? parser->range_from : ll;

            v = end + lws_spn(end);
            if (*v != '-') goto fail;
            v += 1 + lws_spn(v + 1);
            if (*v == '\0') break;

            ll = strtoll(v, &end, 10); /* last-byte-pos */
            parser->flags |= REQ_RANGE_TO;
            parser->range_to = parser->range_to > ll ? parser->range_to : ll;
            v = end + lws_spn(end);
        }
        if (*v == '\0') break;
        if (*v != ',') goto fail;
    }
    /* Convert from + suffix into just from */
    if ((parser->flags & (REQ_RANGE_FROM|REQ_RANGE_SUFFIX))
        == (REQ_RANGE_FROM|REQ_RANGE_SUFFIX)) {
        parser->flags &= ~(REQ_RANGE_TO|REQ_RANGE_SUFFIX);
    }

    return;

 fail:
    parser->flags &= ~(REQ_RANGE_FROM|REQ_RANGE_TO|REQ_RANGE_SUFFIX);
}

static int parser_parse_header(Http_Parser *parser) {
    int res = 0;

    if (strcasecmp(parser->key, "Content-Length") == 0) {
        size_t len;
        char *end;

        if (!parser->val_used) {
            parser->state = ERR_BAD_REQUEST; res = 1;
        } else {
            len = strtoul(parser->val, &end, 10);
            if (*parser->val == '\0' || *end != '\0') {
                parser->state = ERR_BAD_REQUEST; res = 1;
            } else {
                parser->content_length = len;
            }
        }
    } else if (strcasecmp(parser->key, "Transfer-Encoding") == 0) {
        if (!parser->val_used) {
            parser->state = ERR_BAD_REQUEST; res = 1;
        } else if (strcasecmp(parser->val, "identity") == 0) {
            parser->trans_enc = TE_IDENT;
        } else if (strcasecmp(parser->val, "chunked") == 0) {
            parser->trans_enc = TE_CHUNKED;
        } else {
            parser->state = ERR_UNIMPLEMENTED; res = 1;
        }
    } else if (strcasecmp(parser->key, "Connection") == 0) {
        if ((parser->flags & REQ_KEEP_ALIVE) != 0) {
            char *v = parser->val;

            while (*v != '\0') {
                size_t l = lws_cspn(v);
                if (l == 5 && strncasecmp(v, "close", l) == 0) {
                    parser->flags &= ~REQ_KEEP_ALIVE;
                    break;
                }
                v += l + lws_spn(v + l);
            }
        }
    } else if (strcasecmp(parser->key, "User-Agent") == 0) {
        if (parser->val_used) {
            free(parser->user_agent);
            parser->user_agent = lim_strdup(parser->val, parser->val_used, MAX_UA_LEN);
        }
    } else if (strcasecmp(parser->key, "Referer") == 0) {
        if (parser->val_used) {
            free(parser->referrer);
            parser->referrer = lim_strdup(parser->val, parser->val_used, MAX_REFERRER_LEN);
        }
    } else if (strcasecmp(parser->key, "Range") == 0) {
        if (parser->val_used) parse_range(parser);
    }
    parser->key_used = 0; *parser->key = '\0';
    parser->val_used = 0; if (parser->val) *parser->val = '\0';
    return res;
}

static void read_headers(Http_Parser *parser) {
    char *line;
    size_t len, keylen, spaces, valpos;

    line = parser_get_line(parser, &len);

    while (line != NULL) {

        if (len == 0) { /* End of headers */
            /* Store last key/value if present */
            if (parser->key_used) {
                if (parser_parse_header(parser) != 0) return;
            }

            /* Go on to next stage */
            switch (parser->trans_enc) {
            case TE_IDENT:
                parser->state = parser->content_length ? READING_BODY : GOT_REQUEST;
                parser->bytes = parser->content_length;
                return;
            case TE_CHUNKED:
                parser->state = READING_CHUNK_HEADER; return;
            default:
                parser->state = ERR_UNIMPLEMENTED; return;
            }
        }

        keylen = tok_spn(line);
        spaces = lws_spn(line + keylen);
        if (keylen > 0) {
            /* Check for colon */
            if (line[keylen + spaces] != ':') {parser->state = ERR_BAD_REQUEST; return; }
            spaces += 1 + lws_spn(line + keylen + spaces + 1);

            /* Store previous key/value if present */
            if (parser->key_used) {
                if (parser_parse_header(parser) != 0) return;
            }

            /* Copy key */
            if (parser->key_sz < keylen + 1) {
                if (parser->key) { free(parser->key); parser->key_sz = 0; }
                parser->key = malloc(keylen + 1);
                if (parser->key == NULL) { parser->state = ERR_INTERNAL; return; }
                parser->key_sz = keylen + 1;
            }
            memcpy(parser->key, line, keylen); parser->key[keylen] = '\0';
            parser->key_used = keylen;
            parser->val_used = 0;
        } else if (spaces == 0) {
            /* Bad char in key */
            parser->state = ERR_BAD_REQUEST; return;
        }

        /* Copy value, if present */
        valpos = keylen + spaces;
        if (line[valpos] != '\0') {
            if (parser->val_used + len - valpos + 2 > parser->val_sz) {
                char *new_val;
                size_t new_sz = parser->val_sz ? parser->val_sz * 2 : 64;
                while (new_sz < parser->val_used + len - valpos + 2) new_sz *= 2;
                new_val = realloc(parser->val, new_sz);
                if (new_val == NULL) { parser->state = ERR_INTERNAL; return; }
                parser->val = new_val;
                parser->val_sz = new_sz;
            }
            if (parser->val_used > 0) parser->val[parser->val_used++] = ' '; /* continuation */
            assert(parser->val != NULL);
            memcpy(parser->val + parser->val_used, line + valpos, len - valpos);
            parser->val_used += len - valpos;
            parser->val[parser->val_used] = '\0';
        }
        line = parser_get_line(parser, &len);
    }
    return; /* Wait for more data */
}

static void read_chunk_header(Http_Parser *parser) {
    char *line, *end;
    size_t len, bytes;

    line = parser_get_line(parser, &len);
    if (line == NULL) return;

    bytes = strtoul(line, &end, 16);
    if (*end != '\0' && !is_lws(*end) && *end != ';') {
        /* Bad hex string */
        parser->state = ERR_BAD_REQUEST; return;
    }

    parser->bytes = bytes;
    parser->state = bytes ? READING_CHUNK : READING_CHUNK_TRAILER;
}

static inline void eat_data(Http_Parser *parser) {
    unsigned int l;

    if (parser->used == 0) {
        assert(parser->in == parser->out);
        assert(parser->pos == parser->out);
        return;
    }

    if (parser->in <= parser->out) {
        // Data available from out to end
        l = BUF_SZ - parser->out;
        if (l > parser->bytes)
            l = (unsigned int) parser->bytes;
        assert(l <= parser->used);
        parser->out = (parser->out + l) & BUF_MASK;
        parser->bytes -= l;
        parser->used -= l;
    }
    if (parser->out < parser->in) {
        // Data available from out to in
        l = parser->in - parser->out;
        if (l > parser->bytes)
            l = (unsigned int) parser->bytes;
        assert(l <= parser->used);
        parser->out   += l;
        parser->bytes -= l;
        parser->used -= l;
    }
    parser->pos = parser->out;
}

static void read_chunk(Http_Parser *parser) {
    /* Throw away for now */
    eat_data(parser);
    if (parser->bytes == 0) parser->state = READING_CHUNK_TRAILER;
}

static void read_body(Http_Parser *parser) {
    /* Throw away for now */
    eat_data(parser);
    if (parser->bytes == 0) parser->state = GOT_REQUEST;
}

static void read_chunk_trailer(Http_Parser *parser) {
    size_t len;
    /* Throw away for now */
    while (parser_get_line(parser, &len) != NULL) {
        if (len == 0) { /* Finished */
            parser->state = GOT_REQUEST; return;
        }
    }
    return;  /* Wait for more lines */
}

Read_result parser_read_data(const Options *opts, Client *client,
                             Http_Parser *parser, int fd) {
    Read_result res;
    Parse_state last_state;

    res = parser_read_input(parser, fd);
    if (res == READ_EOF || res == READ_ERROR) return res;

    do {
        Transaction *transact = NULL;
        last_state = parser->state;
        switch (parser->state) {
        case READING_REQUEST_LINE:  read_request_line(opts, parser);  break;
        case READING_HEADERS:       read_headers(parser);             break;
        case READING_CHUNK_HEADER:  read_chunk_header(parser);        break;
        case READING_CHUNK:         read_chunk(parser);               break;
        case READING_BODY:          read_body(parser);                break;
        case READING_CHUNK_TRAILER: read_chunk_trailer(parser);       break;
        case GOT_REQUEST:
            handle_request(opts, client, parser, &transact);
            break;
        case SHUTTING_DOWN:         return READ_EOF;
        default:
            handle_error(client, parser, (int) parser->state, &transact);
            break;
        }

        if (transact != NULL) {
            client_add_transaction(client, transact);
        }

    } while (last_state != parser->state);

    return res;
}

char * steal_user_agent_from_parser(Http_Parser *parser) {
    char *ua = parser->user_agent;
    parser->user_agent = NULL;
    return ua;
}

char * steal_referrer_from_parser(Http_Parser *parser) {
    char *referrer = parser->referrer;
    parser->referrer = NULL;
    return referrer;
}
