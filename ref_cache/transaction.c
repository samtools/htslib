/*  transaction.c -- ref-cache http transactions

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/uio.h>
#include <time.h>
#include <assert.h>

#include "transaction.h"
#include "http_parser.h"
#include "options.h"
#include "poll_wrap.h"
#include "ref_files.h"
#include "sendfile_wrap.h"
#include "server.h"
#include "../cram/pooled_alloc.h"

#if defined(REF_CACHE_NO_POOLED_ALLOC)
#  define pool_alloc(x) malloc(x->dsize)
#  define pool_free(p, x) free(x)
#endif

/* Transaction state */
typedef enum {
    TRANSACT_WAITING_TEXT = 0,
    TRANSACT_GOT_TEXT,
    TRANSACT_SENDING_TEXT,
    TRANSACT_SENDING_FILE,
    TRANSACT_FINISHED
} Transaction_state;

/* Transaction flags */

const unsigned int TRANSACT_TEXT_CONST       =  1;
const unsigned int TRANSACT_KEEP_ALIVE       =  2;
const unsigned int TRANSACT_READING_PIPE     =  4;
const unsigned int TRANSACT_RANGE_FROM       =  8;
const unsigned int TRANSACT_RANGE_TO         = 16;
const unsigned int TRANSACT_RANGE_SUFFIX     = 32;
const unsigned int TRANSACT_SEND_BUFFER_FULL = 64;

struct Transaction {
    Transaction_state state;
    unsigned int  flags;
    struct Transaction *next;
    struct Transaction *next_id;
    Client       *client;
    char         *user_agent;
    char         *referrer;
    char         *req_str;
    char         *text;
    off_t         range_from;
    off_t         range_to;
    size_t        sz;
    size_t        out;
    RefFile      *ref;
    off_t         fd_sz;
    off_t         fd_sent;
    unsigned int  rc;
    Http_version  http_vers;
#ifdef HAVE_SENDFILE
#define HAVE_SEND_BUFFER 0
#else
#define SEND_BUFFER_SIZE 0x40000
#define SEND_BUFFER_MASK (SEND_BUFFER_SIZE - 1)
#define HAVE_SEND_BUFFER 1
    unsigned char *buffer;
    unsigned int   buffer_in;
    unsigned int   buffer_out;
#endif
};

#define ERROR_1_0_TXT(S, L) "HTTP/1.0 " S "\r\n"             \
                            "Connection: close\r\n"          \
                            "Content-Type: text/plain\r\n"   \
                            "Content-Length: " #L "\r\n\r\n" \
                            S "\r\n"

#define ERROR_1_1_TXT(S, L) "HTTP/1.1 " S "\r\n"             \
                            "Connection: close\r\n"          \
                            "Content-Type: text/plain\r\n"   \
                            "Content-Length: " #L "\r\n\r\n" \
                            S "\r\n"

#define ERROR_TXT(S, L) { \
    ERROR_1_0_TXT(S, L),  \
    ERROR_1_1_TXT(S, L)   \
}

/* Static error texts so we don't have to malloc them */
static char *err_400_txt[2] = ERROR_TXT("400 Bad Request", 17);
static char *err_404_txt[2] = ERROR_TXT("404 Not found", 15);
static char *err_413_txt[2] = ERROR_TXT("413 Request Entity Too Large", 30);
static char *err_414_txt[2] = ERROR_TXT("414 Request-URI Too Large", 27);
static char *err_500_txt[2] = ERROR_TXT("500 Internal Server Error", 27);
static char *err_501_txt[2] = ERROR_TXT("501 Not Implemented", 21);
static char *err_502_txt[2] = ERROR_TXT("502 Bad Gateway", 17);
static char *err_505_txt[2] = ERROR_TXT("505 HTTP Version not supported", 32);

/* Lookup table of transactions by id, for communicating with upstream */
#define TRANSACTIONS 0x400
#define TRANSACT_MASK (TRANSACTIONS - 1)
Transaction *transactions[TRANSACTIONS] = { 0 };

/* Pool for transaction structs */
static pool_alloc_t *transaction_pool = NULL;

Transaction *new_transaction(Client *client, Http_Parser *parser) {
    Transaction *transact;

    if (transaction_pool == NULL) {
        transaction_pool = pool_create(sizeof(Transaction));
        if (transaction_pool == NULL) return NULL;
    }

    transact = pool_alloc(transaction_pool);
    if (transact == NULL) return NULL;
    memset(transact, 0, sizeof(*transact));

    transact->client     = client;
    transact->ref        = NULL;
    transact->user_agent = steal_user_agent_from_parser(parser);
    transact->referrer   = steal_referrer_from_parser(parser);
    transact->req_str    = NULL;
    transact->range_from = parser->range_from;
    transact->range_to   = parser->range_to;
    transact->flags      = parser->flags & (REQ_KEEP_ALIVE|REQ_RANGE_FROM|REQ_RANGE_TO|REQ_RANGE_SUFFIX);
    transact->rc         = 0;
    transact->next_id    = NULL;
    transact->http_vers  = parser->http_vers;
#if HAVE_SEND_BUFFER
    transact->buffer     = NULL;
#endif

    return transact;
}

static void transaction_clear_ref(Transaction *transact) {
    if (transact->ref == NULL) return;
    unsigned int id = get_ref_id(transact->ref);
    if (transactions[id & TRANSACT_MASK] == transact) {
        transactions[id & TRANSACT_MASK] = transact->next_id;
    } else {
        Transaction *t = transactions[id & TRANSACT_MASK];
        while (t != NULL && t->next_id != transact) t = t->next_id;
        if (t) t->next_id = transact->next_id;
    }
    release_ref_file(transact->ref);
    transact->ref = NULL;
    transact->next_id = NULL;
}

void free_transaction(Transaction *transact) {
    transaction_clear_ref(transact);
    if ((transact->flags & TRANSACT_TEXT_CONST) == 0 && transact->text != NULL)
        free(transact->text);
    free(transact->user_agent);
    free(transact->referrer);
    free(transact->req_str);
#if HAVE_SEND_BUFFER
    free(transact->buffer);
#endif
    memset(transact, 0, sizeof(*transact));

    pool_free(transaction_pool, transact);
}

void free_transaction_list(Transaction *head) {
    while (head != NULL) {
        Transaction *next = head->next;
        free_transaction(head);
        head = next;
    }
}

Transaction * switch_to_next_transaction(Transaction *transact) {
    Transaction *next = transact->next;
    free_transaction(transact);
    return next;
}

Client * transaction_get_client(Transaction *transact) {
    return transact->client;
}

int transaction_get_keep_alive(Transaction *transact) {
    return (transact->flags & TRANSACT_KEEP_ALIVE) != 0;
}

void transaction_set_ref(Transaction *transact, RefFile *ref) {
    transact->ref = ref;
    unsigned int id = get_ref_id(ref);
    transact->next_id = transactions[id & TRANSACT_MASK];
    transactions[id & TRANSACT_MASK] = transact;
}

static void calculate_range_available(Transaction *transact, off_t size,
                                      off_t *range_start_out,
                                      off_t *range_end_out) {
    /*
      Calculate the file range to be returned, based on the request in
      the http headers, and the given size.
      If the requested range is not valid (e.g. it starts beyond size),
      we revert to sending the whole file, as permitted by the http
      specification.
    */
    off_t range_start = -1, range_end = -1;
    int have_range = (0 != (transact->flags & (TRANSACT_RANGE_FROM
                                               | TRANSACT_RANGE_SUFFIX)));

    if (have_range) {
        if (0 != (transact->flags & REQ_RANGE_SUFFIX)) {
            // Request for the last suffix-length bytes in the file
            range_end = size;
            range_start = (transact->range_to < range_end
                           ? range_end - transact->range_to : 0);
            if (range_start == 0 || transact->range_to == 0) have_range = 0;
        } else {
            // Request for bytes range-from ... range-to or range-from ... end
            if (transact->range_from > transact->range_to
                || transact->range_from >= size) {
                have_range = 0;
            } else {
                range_start = transact->range_from;
                range_end   = (0 != (transact->flags & TRANSACT_RANGE_TO)
                               ? (transact->range_to + 1 < size
                                  ? transact->range_to + 1 : size)
                               : size);
            }
        }
    }
    *range_start_out = have_range ? range_start : -1;
    *range_end_out   = have_range ? range_end   : -1;
}

void transaction_set_req_str(Transaction *transact, const char *requested) {
    size_t len = strlen(requested);
    transact->req_str = strndup(requested,
                                len < MAX_REQUEST_LEN ? len : MAX_REQUEST_LEN);
}

static Transaction * transaction_by_id(unsigned int id, Transaction *start) {
    Transaction *r = start == NULL ? transactions[id & TRANSACT_MASK] : start->next_id;
    while (r != NULL && r->ref != NULL && get_ref_id(r->ref) != id)
        r = r->next_id;
    return r;
}

void set_error_response(Transaction *transact, unsigned int code) {
    unsigned int vers;

    if (transact->state >= TRANSACT_SENDING_TEXT) {
        // To late to send an error.  Best we can do is abandon the transaction
        // and hang up.
        transact->state =  TRANSACT_FINISHED;
        transact->flags &= ~TRANSACT_KEEP_ALIVE;
        return;
    }

    if ((transact->flags & TRANSACT_TEXT_CONST) == 0 && transact->text != NULL)
        free(transact->text);

    vers = transact->http_vers == HTTP_1_1 ? 1 : 0;

    switch (code) {
      case 400: transact->text = err_400_txt[vers]; break;
      case 404: transact->text = err_404_txt[vers]; break;
      case 413: transact->text = err_413_txt[vers]; break;
      case 414: transact->text = err_414_txt[vers]; break;
      case 500: transact->text = err_500_txt[vers]; break;
      case 501: transact->text = err_501_txt[vers]; break;
      case 502: transact->text = err_502_txt[vers]; break;
      case 505: transact->text = err_505_txt[vers]; break;
      default:  transact->text = err_500_txt[vers]; break;
    }

    transact->rc = code;
    transact->sz = strlen(transact->text);
    transaction_clear_ref(transact);
    transact->flags &= ~TRANSACT_KEEP_ALIVE;
    transact->flags |= TRANSACT_TEXT_CONST;
    transact->state = TRANSACT_GOT_TEXT;
}

static inline const char * http_vers_string(Transaction *transact) {
    switch (transact->http_vers) {
      case HTTP_1_0: return "HTTP/1.0";
      case HTTP_1_1: return "HTTP/1.1";
      default: break;
    }
    return NULL;
}

static int set_ref_file_response(Transaction *transact, int64_t len,
                                 off_t range_start, off_t range_end) {
    size_t text_len = 128;
    char *text;
    const char *vers = http_vers_string(transact);
    const char *content_type = "text/plain";
    int keep_alive = ((transact->flags & TRANSACT_KEEP_ALIVE) != 0
                      && transact->http_vers == HTTP_1_1);
    unsigned int rc = 200;
    int l;

    assert(transact->state < TRANSACT_SENDING_TEXT);

    if (!vers) {
        set_error_response(transact, 505);
        return -1;
    }

    if ((transact->flags & TRANSACT_TEXT_CONST) == 0
        && transact->text != NULL) {
        free(transact->text);
        transact->text = NULL;
    }

    if (range_start >= 0)
        text_len += 100;

    text = malloc(text_len);
    if (text == NULL) {
        set_error_response(transact, 500);
        return -1;
    }

    if (len) {
        if (range_start >= 0) {
            assert(range_end > range_start && range_end > 0);
            rc = 206;
            l = snprintf(text, text_len,
                         "%s 206 Partial content\r\n"
                         "Content-Type: %s\r\n"
                         "Content-Range: bytes %lld-%lld/%lld\r\n"
                         "Content-Length: %lld\r\n"
                         "%s\r\n",
                         vers, content_type,
                         (long long) range_start,
                         (long long) range_end - 1,
                         (long long) len,
                         (long long) (range_end - range_start),
                         keep_alive ? "" : "Connection: close\r\n");
        } else {
            l = snprintf(text, text_len,
                         "%s 200 OK\r\n"
                         "Content-Type: %s\r\n"
                         "Content-Length: %lld\r\n"
                         "%s\r\n",
                         vers, content_type, (long long) len,
                         keep_alive ? "" : "Connection: close\r\n");
        }
    } else {
        l = snprintf(text, text_len,
                     "%s 200 OK\r\n"
                     "Content-Type: %s\r\n"
                     "Connection: close\r\n\r\n",
                     vers, content_type);
        keep_alive = 0;
    }

    if ((size_t) l > text_len) { // Shouldn't happen
        free(text);
        set_error_response(transact, 500);
        return -1;
    }

    transact->rc = rc;
    transact->text = text;
    transact->sz   = (size_t) l;
    transact->out  = 0;
    transact->flags &= ~TRANSACT_TEXT_CONST;
    transact->state = TRANSACT_GOT_TEXT;

    if (!keep_alive) {
        transact->flags &= ~TRANSACT_KEEP_ALIVE;
    }
    return 0;
}

void set_message_response(Transaction *transact, const char *content_type,
                          const char *message, size_t len) {
    size_t text_len = 128 + len, l;
    char *text;
    const char *vers = http_vers_string(transact);
    int keep_alive = ((transact->flags & TRANSACT_KEEP_ALIVE) != 0
                      && transact->http_vers == HTTP_1_1);

    assert(transact->state < TRANSACT_SENDING_TEXT);

    if (!vers) {
        set_error_response(transact, 505);
    }

    if ((transact->flags & TRANSACT_TEXT_CONST) == 0
        && transact->text != NULL) {
        free(transact->text);
        transact->text = NULL;
    }

    text = malloc(text_len);
    if (text == NULL) {
        set_error_response(transact, 500);
    }

    l = (size_t) snprintf(text, text_len,
                          "%s 200 OK\r\n"
                          "Content-Type: %s\r\n"
                          "Content-Length: %zu\r\n"
                          "%s\r\n",
                          vers, content_type, len,
                          keep_alive ? "" : "Connection: close\r\n");

    if (l > text_len || text_len - l < len) { // Shouldn't happen
        free(text);
        set_error_response(transact, 500);
    }

    memcpy(text + l, message, len);

    transact->rc = 200;
    transact->text = text;
    transact->sz   = l + len;
    transact->out  = 0;
    transact->flags &= ~TRANSACT_TEXT_CONST;
    transact->state = TRANSACT_GOT_TEXT;

    if (!keep_alive) {
        transact->flags &= ~TRANSACT_KEEP_ALIVE;
    }
}

static inline ssize_t send_file(Transaction *transact, int out_fd, int in_fd,
                                off_t end) {
#ifdef HAVE_SENDFILE
    assert(end >= transact->fd_sent);
    return sendfile_wrap(out_fd, in_fd,
                         &transact->fd_sent,
                         (size_t) (end - transact->fd_sent));
#else
    struct iovec iov[2];
    int iovcnt = 1;
    off_t count = end - transact->fd_sent;
    unsigned int buffered;
    ssize_t bytes = 0;

    if (!transact->buffer) {
        transact->buffer = malloc(SEND_BUFFER_SIZE);
        if (!transact->buffer)
            return -1;
        transact->buffer_in = transact->buffer_out = 0;
        transact->flags &= ~TRANSACT_SEND_BUFFER_FULL;
    }

    // Adjust count for bytes already in buffer
    if ((transact->flags & TRANSACT_SEND_BUFFER_FULL) != 0) {
        buffered = SEND_BUFFER_SIZE;
    } else {
        if (transact->buffer_out <= transact->buffer_in) {
            buffered = transact->buffer_in - transact->buffer_out;
        } else {
            buffered = SEND_BUFFER_SIZE - transact->buffer_out + transact->buffer_in;
        }
    }
    count -= buffered;
    assert(count >= 0);

    // Only read if there's data available and the buffer has a reasonable
    // amount of space free (to avoid lots of short reads).
    while (count > 0 && buffered < SEND_BUFFER_SIZE / 2) {
        size_t len = (transact->buffer_in < transact->buffer_out
                      ? transact->buffer_out - transact->buffer_in
                      : SEND_BUFFER_SIZE - transact->buffer_in);
        assert(len <= SEND_BUFFER_SIZE);
        if ((off_t) len > count)
            len = (size_t) count;

        do {
            bytes = pread(in_fd, transact->buffer + transact->buffer_in,
                          len, transact->fd_sent + buffered);
        } while (bytes < 0
                 && (errno == EINTR
                     || errno == EAGAIN || errno == EWOULDBLOCK));
        if (bytes < 0)
            return -1;
        transact->buffer_in = ((transact->buffer_in + (unsigned int) bytes)
                               & SEND_BUFFER_MASK);
        buffered += (unsigned int) bytes;
        count -= bytes;
        assert(count >= 0);

        if (transact->buffer_in == transact->buffer_out)
            transact->flags |= TRANSACT_SEND_BUFFER_FULL;

        if (bytes == 0 || transact->buffer_in != 0)
            break;
        // Otherwise try one more read, to the front of the buffer
    }

    if ((transact->flags & TRANSACT_SEND_BUFFER_FULL) != 0
        || transact->buffer_out != transact->buffer_in) {
        iovcnt = 1;
        iov[0].iov_base = transact->buffer + transact->buffer_out;
        if (transact->buffer_in > transact->buffer_out) {
            iov[0].iov_len  = transact->buffer_in - transact->buffer_out;
        } else {
            iov[0].iov_len  = SEND_BUFFER_SIZE - transact->buffer_out;
            if (transact->buffer_in > 0) {
                iov[1].iov_base = transact->buffer;
                iov[1].iov_len  = transact->buffer_in;
                iovcnt = 2;
            }
        }
        bytes = writev(out_fd, iov, iovcnt);
        if (bytes < 0)
            return -1;
        transact->fd_sent += bytes;
        transact->buffer_out = ((transact->buffer_out + (unsigned int) bytes)
                                & SEND_BUFFER_MASK);
        if (bytes > 0)
            transact->flags &= ~TRANSACT_SEND_BUFFER_FULL;

        return bytes;
    }
    return 0; // Nothing to send
#endif
}

Write_result transaction_send_data(Transaction *transact, int fd) {
    switch (transact->state) {
    case TRANSACT_GOT_TEXT:
        transact->state = TRANSACT_SENDING_TEXT;
        // fall through
    case TRANSACT_SENDING_TEXT: {
        ssize_t bytes;
        assert(transact->text != NULL);
        bytes = write(fd, transact->text, transact->sz - transact->out);
        if (bytes < 0) {
            if (errno != EAGAIN || errno != EWOULDBLOCK || errno != EINTR) {
                fprintf(stderr, "Error from fd #%d : %s\n", fd, strerror(errno));
                return WRITE_ERROR;
            }
            return WRITE_BLOCKED;
        }
        transact->out += (size_t) bytes;
        if (transact->out < transact->sz) return WRITE_BLOCKED;
        if (transact->ref == NULL) {
            transact->state = TRANSACT_FINISHED;
            return WRITE_COMPLETE;
        }
        transact->state = TRANSACT_SENDING_FILE;
    }
        // fall through
    case TRANSACT_SENDING_FILE: {
        off_t end, available;
        ssize_t sent;
        int ref_fd;
        assert(transact->ref != NULL);
        ref_fd = get_ref_fd(transact->ref);
        available = get_ref_available(transact->ref);
        assert(ref_fd >= 0);

        end = get_ref_size(transact->ref) > 0 ? transact->fd_sz : available;

        if (available <= transact->fd_sent) return WRITE_BLOCKED_UPSTREAM;
        if (available < end) end = available;
        sent = send_file(transact, fd, ref_fd, end);
        if (sent < 0) {
            if (errno != EAGAIN || errno != EWOULDBLOCK) {
                fprintf(stderr, "sendfile fd #%d : %s\n", fd, strerror(errno));
                return WRITE_ERROR;
            }
            return WRITE_BLOCKED;
        }
        return (transact->fd_sent < end || !get_ref_complete(transact->ref)
                ? (transact->fd_sent == end ? WRITE_BLOCKED_UPSTREAM : WRITE_BLOCKED)
                : WRITE_COMPLETE);
    }

    case TRANSACT_FINISHED:
        return WRITE_COMPLETE;

    default:
        fprintf(stderr,
                "resp_send_data entered when transaction in wrong state (%d)\n",
                transact->state);
        abort();
    }
    return WRITE_MORE;
}

int transaction_have_content(Transaction *transact) {
    return transact->state > TRANSACT_WAITING_TEXT;
}

void transaction_set_next(Transaction *transact, Transaction *next) {
    assert(transact->next == NULL);
    transact->next = next;
}

int transaction_has_data_to_send(Transaction *transact) {
    int have_data = 0;

    switch (transact->state) {
    case TRANSACT_WAITING_TEXT:
        break;
    case TRANSACT_GOT_TEXT:
    case TRANSACT_SENDING_TEXT:
        if (transact->out < transact->sz) have_data = 1;
        break;
    case TRANSACT_SENDING_FILE: {
        off_t available = get_ref_available(transact->ref);
        off_t end = get_ref_size(transact->ref) > 0 ? transact->fd_sz : available;
        if (available < end) end = available;
        if (transact->fd_sent < end) have_data = 1;
        break;
    }
    case TRANSACT_FINISHED:
        have_data = 1;
        break;
    default:
        abort(); /* Shouldn't happen */
    }
    return have_data;
}

void set_transaction_file_range(Transaction *transact, int64_t size,
                                int ref_data_available) {
    int res = 0;
    off_t range_start = -1, range_end = -1;

    if (size > 0) {
        calculate_range_available(transact, size, &range_start, &range_end);
    }
    if (ref_data_available)
        res = set_ref_file_response(transact, size, range_start, range_end);
    if (res == 0) {
        transact->fd_sz = range_end >= 0
            ? range_end : get_ref_size(transact->ref);
        transact->fd_sent = range_start >= 0 ? range_start : 0;
    }
}

static void update_with_initial_size(Transaction *transact, int64_t size,
                                     unsigned int ref_id,
                                     Client **write_stack) {
    for (; transact != NULL; transact = transaction_by_id(ref_id, transact)) {
        set_transaction_file_range(transact, size, 1);
        queue_transaction_write(transact, write_stack);
    }
}

void got_download_started(unsigned int id, int64_t val, int fd,
                          Client **write_stack) {
    /* Handle a message from upstream to say that a refernce download has
       started. */
    Transaction *transact = transaction_by_id(id, NULL);
    if (transact == NULL) return;

    assert(transact->ref != NULL);
    update_ref_download_started(transact->ref, fd, val);

    if (val >= 0) {
        /* We got a size.  This means the file must already have been downloaded
           and we can treat it as an ordinary file.  Go through the list of
           transactions for this ref and update them. */
        update_with_initial_size(transact, val, id, write_stack);
    }
}

void got_download_part(unsigned int id, int64_t val, Client **write_stack) {
    /* A bit more downloading has happened. Record the new size. */
    assert(val >= 0);

    Transaction *transact = transaction_by_id(id, NULL);
    if (transact == NULL) return;

    assert(transact->ref != NULL);
    update_ref_available(transact->ref, val);

    /* Add any clients that can now make progress to the write list in the
       server. */
    for (; transact != NULL; transact = transaction_by_id(id, transact)) {
        queue_transaction_write(transact, write_stack);
    }
}

void got_download_clen(unsigned int id, int64_t val, Client **write_stack) {
    /* Got content length (or 0 if unknown).  This means a download has
       actually started.  Update the transaction structs for this download. */

    assert(val >= 0);

    Transaction *transact = transaction_by_id(id, NULL);
    if (transact == NULL) return;

    assert(transact->ref != NULL);
    update_ref_with_content_len(transact->ref, val);

    update_with_initial_size(transact, val, id, write_stack);
}

void got_download_result(unsigned int id, int64_t val, Client **write_stack) {

    assert(val >= 0 && val < 1000);

    Transaction *transact = transaction_by_id(id, NULL);
    if (!transact)
        return;

    if (val != 200) {
        // Download errored, pass on the result code
        for (; transact != NULL; transact = transaction_by_id(id, transact)) {
            set_error_response(transact, (unsigned int) val);
            queue_transaction_write(transact, write_stack);
        }
        return;
    }


    assert(transact->ref != NULL);

    int no_initial_content_length = set_ref_complete(transact->ref);
    off_t available = get_ref_available(transact->ref);

    for (; transact != NULL; transact = transaction_by_id(id, transact)) {
        if (no_initial_content_length) {
            transact->fd_sz = available;
        } else if (transact->fd_sz > available) {
            /* Incorrect content-length again ... */
            transact->fd_sz = available;
        }
        if (transact->fd_sent >= transact->fd_sz) {
            transact->state = TRANSACT_FINISHED;
        } else if (transact->state < TRANSACT_SENDING_TEXT) {
            set_transaction_file_range(transact, available, 1);
        }
        queue_transaction_write(transact, write_stack);
    }
}

size_t make_log_message(Transaction *transact, char *buffer, size_t size) {
    int bytes;
    char timestamp[32];
    time_t t = time(NULL);
    strftime(timestamp, sizeof(timestamp), "%d/%b/%Y:%H:%M:%S +0000",
             gmtime(&t));
    bytes = snprintf(buffer, size,
                     "REQ %s - - %s \"%s\" %d %zu \"%s\" \"%s\"\n",
                     client_host(transact->client), timestamp,
                     transact->req_str != NULL ? transact->req_str : "",
                     transact->rc,
                     transact->out + (size_t) transact->fd_sent,
                     transact->user_agent != NULL ? transact->user_agent : "",
                     transact->referrer != NULL ? transact->referrer : "");
    return (size_t) bytes;
}
