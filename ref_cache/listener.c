/*  listener.c -- ref-cache socket listener interface

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

#include "../config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <fcntl.h>
#include <assert.h>

#include "listener.h"
#include "misc.h"
#include "options.h"
#include "poll_wrap.h"

struct Listeners {
  size_t nsocks;
  int *sockets;
};

static inline int open_socket(struct addrinfo *addr, char **cause) {
    int s, serrno, val = 1;

    s = socket(addr->ai_family, addr->ai_socktype, addr->ai_protocol);
    if (s == -1) {
        *cause = "socket";
        return -1;
    }

    if (setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &val, sizeof(val)) < 0) {
        *cause = "setsockopt";
        goto fail;
    }

    if (addr->ai_family == AF_INET6) {
        if (setsockopt(s, IPPROTO_IPV6, IPV6_V6ONLY, &val, sizeof(val)) < 0) {
            *cause = "setsockopt";
            goto fail;
        }
    }

    if (bind(s, addr->ai_addr, addr->ai_addrlen) == -1) {
        *cause = "bind";
        goto fail;
    }
    /* Set non-blocking so we don't get stuck in accept() */
    if (setnonblock(s) != 0) {
        *cause = "setnonblock";
        goto fail;
    }

    if (listen(s, SOMAXCONN) != 0) {
        *cause = "listen";
        goto fail;
    }

    return s;
 fail:
    serrno = errno;
    close(s);
    errno = serrno;
    return -1;
}

Listeners * get_listen_sockets(int port) {
    Listeners *lsocks = calloc(1, sizeof(Listeners));
    struct addrinfo hints, *addr, *addr_list;
    char pnum[20];
    char *cause = NULL;
    size_t count;
    int res;

    if (lsocks == NULL)
        return NULL;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE|AI_NUMERICSERV|AI_ADDRCONFIG;

    snprintf(pnum, sizeof(pnum), "%d", port);
    res = getaddrinfo(NULL, pnum, &hints, &addr_list);
    if (res != 0) {
        fprintf(stderr, "getaddrinfo failed: %s\n", gai_strerror(res));
        free(lsocks);
        return NULL;
    }

    count = 0;
    for (addr = addr_list; addr != NULL; addr = addr->ai_next)
        count++;

    if (count == 0) {
        fprintf(stderr, "getaddrinfo returned nothing.\n");
        free(lsocks);
        return NULL;
    }

    lsocks->sockets = malloc(count * sizeof(*lsocks->sockets));
    if (lsocks->sockets == NULL) {
        perror("Allocating socket list");
        freeaddrinfo(addr_list);
        free(lsocks);
        return NULL;
    }

    for (addr = addr_list; addr != NULL; addr = addr->ai_next) {
        assert(lsocks->nsocks < count);
        lsocks->sockets[lsocks->nsocks] = open_socket(addr, &cause);
        if (lsocks->sockets[lsocks->nsocks] != -1)
            lsocks->nsocks++;
    }

    freeaddrinfo(addr_list);

    if (lsocks->nsocks == 0) {
        fprintf(stderr, "Failure in %s while getting socket: %s\n",
                cause, strerror(errno));
        free(lsocks->sockets);
        free(lsocks);
        return NULL;
    }

    return lsocks;
}

void close_listen_sockets(Listeners *lsocks) {
    size_t i;
    assert(lsocks != NULL);
    for (i = 0; i < lsocks->nsocks; i++)
        close(lsocks->sockets[i]);
}

Pw_item ** register_listener_pollers(Listeners *lsocks, Poll_wrap *pw) {
    Pw_item **polled_listeners = calloc(lsocks->nsocks,
                                        sizeof(*polled_listeners));
    size_t i;

    if (polled_listeners == NULL) {
        perror("Allocating listener poll structs");
        return NULL;
    }

    for (i = 0; i < lsocks->nsocks; i++) {
        polled_listeners[i] = pw_register(pw, lsocks->sockets[i], SV_LISTENER,
                                          PW_IN|PW_ERR|PW_HUP, NULL);
        if (polled_listeners[i] == NULL) {
            perror("Adding listener socket to poller");
            goto fail;
        }
    }

    return polled_listeners;

 fail:
    assert(polled_listeners[i] == NULL);
    if (i > 0) {
        do {
            --i;
            if (polled_listeners[i])
                pw_remove(pw, polled_listeners[i], 0);
        } while (i > 0);
    }
    free(polled_listeners);
    return NULL;
}
