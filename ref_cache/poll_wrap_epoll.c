/*  poll_wrap_epoll.c -- ref-cache wrapper around epoll

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

#if defined(HAVE_EPOLL) && HAVE_EPOLL

#include <sys/epoll.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "poll_wrap.h"
#include "../cram/pooled_alloc.h"

#if defined(REF_CACHE_NO_POOLED_ALLOC)
#  define pool_alloc(x) malloc(x->dsize)
#  define pool_free(p, x) free(x)
#endif

#define INIT_EPOLL_SIZE 128

struct Poll_wrap {
    pool_alloc_t *pool;
    int epfd;
    int debug;
};

Poll_wrap *pw_init(int debug) {
    Poll_wrap *pw = calloc(1, sizeof(Poll_wrap));
    if (pw == NULL) return NULL;

    pw->pool = pool_create(sizeof(Pw_item));
    if (pw->pool == NULL) {
        free(pw);
        return NULL;
    }

    pw->epfd = epoll_create(INIT_EPOLL_SIZE);
    if (pw->epfd < 0) {
        perror("epoll_create");
        pool_destroy(pw->pool);
        free(pw);
        return NULL;
    }

    pw->debug = debug;

    return pw;
}

void pw_close(Poll_wrap *pw) {
    close(pw->epfd);
    pool_destroy(pw->pool);
    free(pw);
}

Pw_item * pw_register(Poll_wrap *pw, int fd, Pw_fd_type fd_type,
                      uint32_t init_events, void *userp) {
    struct epoll_event event;
    Pw_item *item = pool_alloc(pw->pool);
    if (item == NULL)
        return NULL;

    if (pw->debug) {
        fprintf(stderr, "pw_register(%p, %d, %d, 0x%04x, %p)\n",
                (void *) pw, fd, (int) fd_type, init_events, userp);
    }

    item->fd      = fd;
    item->fd_type = fd_type;
    item->userp   = userp;

    event.events = init_events;
    event.data.ptr = item;

    if (epoll_ctl(pw->epfd, EPOLL_CTL_ADD, fd, &event) != 0) {
        perror("epoll_ctl");
        pool_free(pw->pool, item);
        return NULL;
    }

    return item;
}

int pw_mod(Poll_wrap *pw, Pw_item *item, uint32_t events) {
    struct epoll_event event;

    if (pw->debug) {
        fprintf(stderr, "pw_mod(%p, %d, 0x%04x)\n",
                (void *) pw, item->fd, events);
    }

    event.events   = events;
    event.data.ptr = item;

    return epoll_ctl(pw->epfd, EPOLL_CTL_MOD, item->fd, &event);
}

int pw_wait(Poll_wrap *pw, Pw_events *events,
            int max_events, int timeout) {
    return epoll_wait(pw->epfd, (struct epoll_event *) events,
                      max_events, timeout);
}

int pw_remove(Poll_wrap *pw, Pw_item *item, int do_close) {
    struct epoll_event dummy;
    int res;

    if (pw->debug) {
        fprintf(stderr, "pw_remove(%p, %d%s)\n",
                (void *) pw, item->fd, do_close ? ", close" : "");
    }

    if (do_close) {
        res = close(item->fd);
    } else {
        res = epoll_ctl(pw->epfd, EPOLL_CTL_DEL, item->fd, &dummy);
    }
    pool_free(pw->pool, item);
    return res;
}
#else
// Prevent "empty translation unit" errors
const char poll_wrap_epoll_disabled[] = "poll_wrap_epoll not built";
#endif /* HAVE_EPOLL */
