/*  server.c -- ref-cache http server

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
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#define MAX_EVENTS 128

#include "server.h"
#include "http_parser.h"
#include "listener.h"
#include "misc.h"
#include "options.h"
#include "poll_wrap.h"
#include "ref_files.h"
#include "request_handler.h"
#include "transaction.h"
#include "upstream.h"
#include "../cram/pooled_alloc.h"

#if defined(REF_CACHE_NO_POOLED_ALLOC)
#  define pool_alloc(x) malloc(x->dsize)
#  define pool_free(p, x) free(x)
#endif

// May not have been defined, depending on compiler options
#ifndef NI_MAXHOST
#  define NI_MAXHOST 1025
#endif

/* Client flags */
static const unsigned int ON_READ_LIST  = 0x01U;
static const unsigned int ON_WRITE_LIST = 0x02U;
static const unsigned int READ_CLOSED   = 0x04U;
static const unsigned int READ_DRAIN    = 0x08U;
static const unsigned int CAN_WRITE     = 0x10U;


struct Client {
    Pw_item *polled;
    struct Client *prev;
    struct Client *next;
    struct Client *prev_write;
    struct Client *next_write;
    char *host;
    Transaction *transact;
    Transaction *last_transact;
    Http_Parser  parser;
    uint32_t flags;
};

typedef struct {
    pool_alloc_t *client_pool;
    Client *can_read;
    Client *can_write;
    size_t  count;
} Clients;

#define LOG_BUF_SZ   0x10000U
#define LOG_BUF_MASK (LOG_BUF_SZ - 1)

typedef struct {
    char buffer[LOG_BUF_SZ];
    size_t in;
    size_t out;
} Log_buffer;

union SocketAddress {
    struct sockaddr     sa;
    struct sockaddr_in  in4;
    struct sockaddr_in6 in6;
    struct sockaddr_storage storage;
};

int init_clients(Clients *clients) {
    clients->client_pool = pool_create(sizeof(Client));
    if (clients->client_pool == NULL) return -1;
    clients->can_read = clients->can_write = NULL;
    clients->count = 0;
    return 0;
}

static inline void read_stack_push(Client *client, Client **stack) {
    if (client->flags & ON_READ_LIST) // Client may already be on the list
        return;
    client->flags |= ON_READ_LIST;
    if (*stack != NULL) (*stack)->prev = client;
    client->next = *stack;
    client->prev = NULL;
    *stack = client;
}

static inline Client * read_stack_pop(Clients *clients) {
    Client *can_read = clients->can_read;
    if (can_read == NULL)
        return NULL;
    clients->can_read = can_read->next;
    if (clients->can_read != NULL) clients->can_read->prev = NULL;
    can_read->prev = can_read->next = NULL;
    assert(can_read->flags & ON_READ_LIST);
    can_read->flags &= ~ON_READ_LIST;
    return can_read;
}

static inline void read_stack_remove(Clients *clients, Client *client) {
    client->flags &= ~ON_READ_LIST;
    if (clients->can_read == client) {
        clients->can_read = client->next;
        if (clients->can_read != NULL) clients->can_read->prev = NULL;
    } else {
        Client *prev = client->prev;
        Client *next = client->next;
        prev->next = next;
        if (next != NULL) next->prev = prev;
    }
    client->prev = client->next = NULL;
}

static inline void write_stack_push(Client *client, Client **stack) {
    if (client->flags & ON_WRITE_LIST)
        return;
    client->flags |= ON_WRITE_LIST;
    if (*stack != NULL) (*stack)->prev_write = client;
    client->next_write = *stack;
    client->prev_write = NULL;
    *stack = client;
}

static inline Client * write_stack_pop(Clients *clients) {
    Client *can_write = clients->can_write;
    if (can_write == NULL)
        return NULL;
    clients->can_write = can_write->next_write;
    if (clients->can_write != NULL) clients->can_write->prev_write = NULL;
    can_write->prev_write = can_write->next_write = NULL;
    assert(can_write->flags & ON_WRITE_LIST);
    can_write->flags &= ~ON_WRITE_LIST;
    return can_write;
}

static inline void write_stack_remove(Clients *clients, Client *client) {
    client->flags &= ~ON_WRITE_LIST;
    if (clients->can_write == client) {
        clients->can_write = client->next_write;
        if (clients->can_write != NULL) clients->can_write->prev_write = NULL;
    } else {
        Client *prev = client->prev_write;
        Client *next = client->next_write;
        prev->next_write = next;
        if (next != NULL) next->prev_write = prev;
    }
    client->prev_write = client->next_write = NULL;
}

static inline Client * get_free_client(Clients *clients) {
    Client *client = pool_alloc(clients->client_pool);
    memset(client, 0, sizeof(*client));
    ++clients->count;
    return client;
}

static inline void close_client(Clients *clients, Client *client,
                                Poll_wrap *pw) {
    assert(clients->count > 0);
    if (client->polled != NULL) {
        pw_remove(pw, client->polled, 1);
        client->polled = NULL;
    }

    if ((client->flags & ON_READ_LIST) != 0)
        read_stack_remove(clients, client);

    if ((client->flags & ON_WRITE_LIST) != 0)
        write_stack_remove(clients, client);

    cleanup_http_parser(&client->parser);

    free_transaction_list(client->transact);
    client->transact = NULL;
    client->last_transact = NULL;
    if (client->host != NULL) {
        free(client->host);
        client->host = NULL;
    }

    pool_free(clients->client_pool, client);
    --clients->count;
}

static int check_addr_allowed(const Options *opts, sa_family_t family,
                              const union SocketAddress *addr,
                              socklen_t addrlen) {
    const uint8_t *addrp;
    size_t i, start, end;

    if (!opts->match_addrs)
        return 0;

    if (family == AF_INET && addrlen == sizeof(struct sockaddr_in)) {
        addrp = (const uint8_t *) &addr->in4.sin_addr;
        start = 0;
        end = opts->first_ip6;
    } else if (family == AF_INET6 && addrlen == sizeof(struct sockaddr_in6)) {
        const uint8_t v4mapped_prefix[12] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0xff, 0xff
        };
        addrp = (const uint8_t *) &addr->in6.sin6_addr;
        if (memcmp(addrp, v4mapped_prefix, sizeof(v4mapped_prefix)) == 0) {
            // V4MAPPED addresses are really IPv4, with the mapped address
            // in the last four bytes
            family = AF_INET;
            addrp += 12;
            start = 0;
            end = opts->first_ip6;
        } else {
            // Ordinary IP6
            start = opts->first_ip6;
            end = opts->num_match_addrs;
        }
    } else {
        return -1;
    }

    for (i = start; i < end; i++) {
        MatchAddr *ma = &opts->match_addrs[i];
        if (family != ma->family)
            continue;
        if (memcmp(addrp, ma->addr, ma->mask_bytes) == 0
            && (addrp[ma->mask_bytes] & ma->mask) == (ma->addr[ma->mask_bytes] & ma->mask)) {
            return 0;
        }
    }
    return -1;
}

static inline void handle_incoming(const Options *opts, int listener,
                                   int upstream, Poll_wrap *pw,
                                   Clients *clients) {
    union SocketAddress addr;
    socklen_t addrlen = sizeof(addr);
    Client *client;
    int incoming;
    int res;
    char host[NI_MAXHOST];

    while ((incoming = accept(listener, &addr.sa, &addrlen))>=0){

        res = getnameinfo(&addr.sa, addrlen, host, sizeof(host),
                          NULL, 0, NI_NUMERICHOST);
        if (res != 0) {
            fprintf(stderr, "Couldn't resolve host for incoming client: %s\n",
                    gai_strerror(res));
            strcpy(host, "<unknown>");
        }

        if (check_addr_allowed(opts, addr.sa.sa_family, &addr, addrlen) != 0) {
            if (opts->verbosity > 1) {
                fprintf(stderr, "Rejected incoming client from %s\n", host);
            }
            close(incoming);
            continue;
        }

        if (opts->verbosity > 1) {
            fprintf(stderr, "Accepted incoming client from %s on fd #%d\n",
                    host, incoming);
        }

        if (setnonblock(incoming) != 0) {
            close(incoming);
            return;
        }

        client = get_free_client(clients);
        if (client == NULL) {
            perror("Getting free client struct");
            close(incoming);
            return;
        }
#if defined(PW_HAVE_EDGE) && PW_HAVE_EDGE
        client->polled = pw_register(pw, incoming, SV_CLIENT,
                                     PW_IN|PW_OUT|PW_ERR|PW_HUP|PW_EDGE, client);
#else
        client->polled = pw_register(pw, incoming, SV_CLIENT,
                                     PW_IN|PW_ERR|PW_HUP, client);
#endif
        if (client->polled == NULL) {
            perror("Registering client with poll");
            close_client(clients, client, pw);
            return;
        }
        client->flags = 0;
        client->host = strdup(host);

        client->transact = client->last_transact = NULL;
        if (init_http_parser(&client->parser, upstream) < 0) {
            perror("Setting up Http_Parser struct");
            close_client(clients, client, pw);
            return;
        }
    }

    switch (errno) {
    case EAGAIN: case EINTR:
    case ENETDOWN: case EPROTO: case ENOPROTOOPT:
    case EHOSTUNREACH: case EOPNOTSUPP: case ENETUNREACH:
#if defined(HAVE_DECL_EHOSTDOWN) && HAVE_DECL_EHOSTDOWN
    case EHOSTDOWN:
#endif
#if defined(HAVE_DECL_ENONET) && HAVE_DECL_ENONET
    case ENONET:
#endif
      break;
    default:
        perror("Error while accepting incoming client");
    }
    return;
}

void client_add_transaction(Client *client, Transaction *transact) {
    if (client->transact == NULL) {
        assert(client->last_transact == NULL);
        client->transact = client->last_transact = transact;
    } else {
        transaction_set_next(client->last_transact, transact);
        client->last_transact = transact;
    }
}

static void handle_upstream(int upstream, Clients *clients, int verbosity) {
    int res;
    Upstream_msg msg;
    int       fd = -1;

    res = upstream_recv_msg(upstream, &msg, &fd);
    if (res <= 0) return;

    if (verbosity > 2) {
        fprintf(stderr, "Message from upstream: %u, %d, %lld\n",
                msg.id, msg.code, (long long) msg.val);
    }

    switch (msg.code) {
    case US_START:
        got_download_started(msg.id, msg.val, fd, &clients->can_write);
        break;
    case US_PARTIAL_LENGTH:
        got_download_part(msg.id, msg.val, &clients->can_write);
        break;
    case US_CONTENT_LENGTH:
        got_download_clen(msg.id, msg.val, &clients->can_write);
        break;
    case US_RESULT:
        got_download_result(msg.id, msg.val, &clients->can_write);
        break;
    default:
        abort(); /* Shouldn't happen */
    }
}

/* Callback so got_download_* functions called above can push
   transactions that can make more progress on to the write_stack */
void queue_transaction_write(Transaction *transact, Client **write_stack) {
    Client *client = transaction_get_client(transact);

    /* Check if already on queue */
    if ((client->flags & ON_WRITE_LIST) != 0)
        return;

    /* Check if at the head of the transaction list */
    if (transact != client->transact)
        return;

    /* See if there is any data waiting to be sent */
    if (!transaction_has_data_to_send(transact))
        return;

#if defined(PW_HAVE_EDGE) && PW_HAVE_EDGE
    /* Check if blocked.  Only if poll is edge triggered. If it isn't
       we may need to re-enable polling on output, which will happen
       if the attempted write blocks. */
    if ((client->flags & CAN_WRITE) == 0)
        return;
#endif

    /* Put onto write_stack */
    write_stack_push(client, write_stack);
}

static Read_result eat_input(int fd) {
    char buffer[65536];
    ssize_t bytes;

    do {
        bytes = read(fd, buffer, sizeof(buffer));
    } while (bytes < 0 && errno == EINTR);
    if (bytes < 0 && (errno == EAGAIN || errno == EWOULDBLOCK))
        return READ_BLOCKED;
    if (bytes < 0) return READ_ERROR;
    if (bytes == 0) return READ_EOF;
    return ((size_t) bytes < sizeof(buffer)) ? READ_BLOCKED : READ_MORE;
}

static void queue_log_msg(Log_buffer *log_buf, Transaction *transact) {
    char buf[MAX_UA_LEN + MAX_REFERRER_LEN + MAX_REQUEST_LEN + NI_MAXHOST + 128];
    char *b = buf;
    size_t bytes = make_log_message(transact, buf, sizeof(buf));
    size_t space;
#if (!defined(NDEBUG))
    assert(bytes < sizeof(buf));
#else
    if (bytes >= sizeof(buf)) {
        bytes = sizeof(buf) - 1;
    }
#endif
    space = (log_buf->in >= log_buf->out
             ? LOG_BUF_SZ - log_buf->in + log_buf->out - 1U
             : log_buf->out - log_buf->in - 1U);
    if (bytes > space) return;

    if (log_buf->in >= log_buf->out) {
        size_t l = LOG_BUF_SZ - log_buf->in - (log_buf->out == 0 ? 1 : 0);
        if (l > bytes) l = bytes;
        memcpy(log_buf->buffer + log_buf->in, b, l);
        b += l;
        bytes -= l;
        log_buf->in = (log_buf->in + l) & LOG_BUF_MASK;
    }
    if (bytes > 0 && log_buf->in < log_buf->out) {
        assert(bytes < log_buf->out - log_buf->in - 1);
        memcpy(log_buf->buffer + log_buf->in, b, bytes);
        log_buf->in += bytes;
    }
    assert(log_buf->in < LOG_BUF_SZ);
    assert(log_buf->in != log_buf->out);
}

static Write_result send_to_log(Log_buffer *log_buf, int log_fd) {
    size_t len;
    ssize_t bytes;
    if (log_buf->out == log_buf->in) return WRITE_MORE;
    if (log_buf->out > log_buf->in) {
        len = LOG_BUF_SZ - log_buf->out;
        do {
            bytes = write(log_fd, log_buf->buffer + log_buf->out, len);
        } while (bytes < 0 && errno == EINTR);

        if (bytes < 0) {
            if (errno == EAGAIN || errno == EWOULDBLOCK) return WRITE_BLOCKED;
            fprintf(stderr, "Couldn't write to log: %s\n", strerror(errno));
            return WRITE_ERROR;
        }

        log_buf->out = (log_buf->out + (size_t) bytes) & LOG_BUF_MASK;
        if ((size_t) bytes < len) return WRITE_BLOCKED;
    }

    if (log_buf->out == log_buf->in) return WRITE_MORE;
    assert(log_buf->out < log_buf->in);

    len = log_buf->in - log_buf->out;
    do {
        bytes = write(log_fd, log_buf->buffer + log_buf->out, len);
    } while (bytes < 0 && errno == EINTR);

    if (bytes < 0) {
        if (errno == EAGAIN || errno == EWOULDBLOCK) return WRITE_BLOCKED;
        fprintf(stderr, "Couldn't write to log: %s\n", strerror(errno));
        return WRITE_ERROR;
    }

    log_buf->out += (size_t) bytes;
    assert(log_buf->out < LOG_BUF_SZ);
    return ((size_t) bytes < len) ? WRITE_BLOCKED : WRITE_MORE;
}

static void handle_client_events(Clients *clients, const Options *opts,
                                 Poll_wrap *pw, Pw_item *item,
                                 unsigned int evts) {
    Client *client = (Client *) item->userp;
    assert(client != NULL);

    if (opts->verbosity > 2) {
        fprintf(stderr, "Events %04x for fd #%d\n", evts, client->polled->fd);
    }
    if ((evts & PW_HUP) != 0) {
        if (opts->verbosity > 1) {
            fprintf(stderr, "Got hangup on fd#%d\n", client->polled->fd);
        }
        close_client(clients, client, pw);
        return;
    } else if ((evts & (PW_IN|PW_ERR)) != 0) {
#if !defined(PW_HAVE_EDGE) || !PW_HAVE_EDGE
        uint32_t out = ((client->flags & ON_WRITE_LIST) == 0
                        && client->transact != NULL
                        && transaction_have_content(client->transact))
            ? PW_OUT : 0;
        if (pw_mod(pw, item, out|PW_ERR|PW_HUP) != 0) {
            fprintf(stderr, "Altering poll settings for fd#%d : %s\n",
                    item->fd, strerror(errno));
            close_client(clients, client, pw);
            return;
        }
#endif /* PW_HAVE_EDGE */
        /* Add to can_read list (if not already there) */
        if ((client->flags & READ_CLOSED) == 0) {
            read_stack_push(client, &clients->can_read);
        } else {
            assert((client->flags & ON_READ_LIST) == 0);
        }
    }

    if ((evts & PW_OUT) != 0) {
#if !defined(PW_HAVE_EDGE) || !PW_HAVE_EDGE
        uint32_t in = (client->flags & ON_READ_LIST) == 0 ? PW_IN : 0;
        if (pw_mod(pw, item, in|PW_ERR|PW_HUP) != 0) {
            fprintf(stderr, "Altering poll settings for fd#%d : %s\n",
                    item->fd, strerror(errno));
            close_client(clients, client, pw);
            return;
        }
#endif /* PW_HAVE_EDGE */
        client->flags |= CAN_WRITE;
        if ((client->flags & ON_WRITE_LIST) == 0
            && client->transact != NULL
            && transaction_have_content(client->transact)) {
            /* Add to the can_write list */
            write_stack_push(client, &clients->can_write);
        }
    }
}

static void do_reads(Clients *clients, const Options *opts, Poll_wrap *pw) {
    Client *can_read = NULL;
    Client *client = NULL;

    while ((client = read_stack_pop(clients)) != NULL) {
        Read_result res;

        if (opts->verbosity > 2) {
            fprintf(stderr, "Reading fd #%d\n", client->polled->fd);
        }
        if ((client->flags & READ_DRAIN) == 0) {
            res = parser_read_data(opts, client, &client->parser,
                                   client->polled->fd);
            if (opts->verbosity > 2) {
                fprintf(stderr, "parser_read_data returned %d\n", res);
            }
            // parser_read_data() may have added new transactions,
            // so check if this client can be added to the write list
            if ((client->flags & (CAN_WRITE|ON_WRITE_LIST)) == CAN_WRITE
                && client->transact
                && transaction_have_content(client->transact)) {
                write_stack_push(client, &clients->can_write);
            }
        } else {
            res = eat_input(client->polled->fd);
            if (opts->verbosity > 2) {
                fprintf(stderr, "eat_input returned %d\n", res);
            }
        }

        switch (res) {
        case READ_MORE:
            if (opts->verbosity > 2) {
                fprintf(stderr, "fd #%d has more data...\n", client->polled->fd);
            }
            read_stack_push(client, &can_read);
            break;
        case READ_EOF:
            client->flags |= READ_CLOSED;
            /* Falls through */
        case READ_BLOCKED:
            if ((client->flags & (READ_CLOSED|READ_DRAIN)) != 0
                && (client->transact == NULL
                    || !transaction_have_content(client->transact))) {
                if (opts->verbosity > 1) {
                    fprintf(stderr, "Closing fd #%d (EOF)\n", client->polled->fd);
                }
                assert((client->flags & ON_READ_LIST) == 0);
                close_client(clients, client, pw);
            } else {
#if !defined(PW_HAVE_EDGE) || !PW_HAVE_EDGE
                uint32_t out = ((client->flags & ON_WRITE_LIST) == 0
                                && client->transact != NULL
                                && transaction_have_content(client->transact))
                    ? PW_OUT : 0;
                uint32_t in = (client->flags & READ_CLOSED) == 0 ? PW_IN : 0;
                if (pw_mod(pw, client->polled, out|in|PW_ERR|PW_HUP) != 0) {
                    fprintf(stderr, "Altering poll settings for fd#%d : %s\n",
                            client->polled->fd, strerror(errno));
                    assert((client->flags & ON_READ_LIST) == 0);
                    close_client(clients, client, pw);
                }
#endif /* PW_HAVE_EDGE */
            }
            break;
        case READ_ERROR:
            if (opts->verbosity > 1) {
                fprintf(stderr, "Closing fd #%d (read error)\n", client->polled->fd);
            }
            assert((client->flags & ON_READ_LIST) == 0);
            close_client(clients, client, pw);
            break;
        }
    }
    clients->can_read = can_read;
}

static void do_writes(Clients *clients, const Options *opts, Poll_wrap *pw,
                      Log_buffer *log_buf) {
    Client *can_write = NULL;
    Client *client = NULL;

    while ((client = write_stack_pop(clients)) != NULL) {
        Write_result res;

        if (opts->verbosity > 1) {
            fprintf(stderr, "Writing to fd #%d\n", client->polled->fd);
        }
        res = transaction_send_data(client->transact, client->polled->fd);
        if (opts->verbosity > 2) {
            fprintf(stderr, "transact_send_data returned %d\n", res);
        }

        switch (res) {
        case WRITE_BLOCKED:
            client->flags &= ~CAN_WRITE;
            // fall through
        case WRITE_BLOCKED_UPSTREAM: {
#if !defined(PW_HAVE_EDGE) || !PW_HAVE_EDGE
            uint32_t in = (client->flags & ON_READ_LIST) == 0 ? PW_IN : 0;
            if (pw_mod(pw, client->polled, in|PW_OUT|PW_ERR|PW_HUP) != 0) {
                fprintf(stderr, "Altering poll settings for fd#%d : %s\n",
                        client->polled->fd, strerror(errno));
                close_client(clients, client, pw);
                break;
            }
#endif /* PW_HAVE_EDGE */
            break;
        }
        case WRITE_COMPLETE: {
            queue_log_msg(log_buf, client->transact);
            if (transaction_get_keep_alive(client->transact)) {
                client->transact = switch_to_next_transaction(client->transact);
            } else {
                // Drop any remaining transactions in queue
                free_transaction_list(client->transact);
                client->transact = NULL;
                client->flags |= READ_DRAIN;
            }

            if (client->transact == NULL) {
                client->last_transact = NULL;
                if ((client->flags & (READ_CLOSED|READ_DRAIN)) != 0) {
                    if (opts->verbosity > 2) {
                        fprintf(stderr, "Shutting down fd #%d\n", client->polled->fd);
                    }
                    if (shutdown(client->polled->fd, SHUT_WR) == -1) {
                        fprintf(stderr, "Error from shutdown(%d, SHUT_WR) : %s\n",
                                client->polled->fd, strerror(errno));
                        close_client(clients, client, pw);
                    }
                }
                break;
            }

            if (!transaction_have_content(client->transact)) break;
        }
            /* else fall through */
        case WRITE_MORE:
            if (opts->verbosity > 2) {
                fprintf(stderr, "fd #%d can send more data...\n", client->polled->fd);
            }
            write_stack_push(client, &can_write);
            break;
        case WRITE_EOF:
        case WRITE_ERROR:
            if (opts->verbosity > 1) {
                fprintf(stderr, "Closing fd #%d\n", client->polled->fd);
            }
            close_client(clients, client, pw);
            break;
        }
    }
    clients->can_write = can_write;
}

int run_poll_loop(const Options *opts, Listeners *lsocks,
                  int upstream, int log_fd) {
    Clients clients;
    Poll_wrap *pw;
    Pw_item  **polled_listeners;
    Pw_item   *polled_upstream = NULL;
    Pw_item   *polled_log = NULL;
    Pw_events  events[MAX_EVENTS];
    Log_buffer log_buf;
    int log_can_write = 0;
    int running = 1;

    log_buf.in = log_buf.out = 0;

    if (setnonblock(log_fd) != 0) {
        perror("Setting nonblocking on log pipe");
        return -1;
    }

    if (init_clients(&clients) != 0) {
        perror("Initializing clients array");
        return -1;
    }

    pw = pw_init(opts->verbosity > 2);
    if (pw == NULL) {
        perror("Initializing poller");
        return -1;
    }

    polled_listeners = register_listener_pollers(lsocks, pw);
    if (polled_listeners == NULL)
        return -1;

    if (upstream != -1) {
        polled_upstream = pw_register(pw, upstream, SV_UPSTREAM,
                                      PW_IN|PW_ERR|PW_HUP, NULL);
        if (polled_upstream == NULL) {
            perror("Adding upstream socket to poller");
            return -1;
        }
    }

    polled_log = pw_register(pw, log_fd, SV_LOG, PW_OUT|PW_ERR|PW_HUP, NULL);
    if (polled_log == NULL) {
        perror("Adding log pipe to poller");
        return -1;
    }

    while (running || clients.count > 0) {
        int e;
        int nevts;

        nevts = pw_wait(pw, events, MAX_EVENTS,
                        clients.can_read == 0 && clients.can_write == 0 ? -1 : 0);
        if (nevts == -1) {
            if (errno == EINTR) /* Signal */
                continue;
            perror("Error from poller");
            return -1;
        }

        // Handle poll events
        for (e = 0; e < nevts; e++) {
            unsigned int  evts = PWE(events[e]);
            Pw_item      *item = PWI(events[e]);
            assert(item != NULL);

            switch (item->fd_type) {
            case SV_LISTENER:
                if (running)
                    handle_incoming(opts, item->fd, upstream, pw, &clients);
                break;
            case SV_UPSTREAM:
                handle_upstream(upstream, &clients, opts->verbosity);
                break;
            case SV_LOG: {
                if (pw_mod(pw, item, PW_ERR|PW_HUP) != 0) {
                    perror("Disarming log poller");
                    return -1;
                }
                if (evts & (PW_HUP | PW_ERR)) {
                    running = 0;
                } else {
                    log_can_write = 1;
                }
                break;
            }
            case SV_CLIENT:
                handle_client_events(&clients, opts, pw, item, evts);
                break;
            default:
                fprintf(stderr, "Unexpected polled item type %d\n", item->fd_type);
                return -1;
            }
        }

        /* Do some reading */
        do_reads(&clients, opts, pw);

        /* Do some writing */
        do_writes(&clients, opts, pw, &log_buf);

        /* Write log output */
        if (log_can_write && log_buf.in != log_buf.out) {
            Write_result res = send_to_log(&log_buf, log_fd);
            if (res == WRITE_ERROR) return -1;
            if (res == WRITE_BLOCKED) {
                if (pw_mod(pw, polled_log, PW_OUT|PW_ERR|PW_HUP) != 0) {
                    perror("Re-arming log poller");
                    return -1;
                }
                log_can_write = 0;
            }
        }
    }
    return 0;
}

const char * client_host(Client *client) {
    return client->host ? client->host : "<NULL>";
}
