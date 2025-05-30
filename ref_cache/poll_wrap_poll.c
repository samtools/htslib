/*  poll_wrap_poll.c -- ref-cache wrapper around poll

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

#if !(defined(HAVE_EPOLL) && HAVE_EPOLL)

#include <poll.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include "poll_wrap.h"
#include "../cram/pooled_alloc.h"

#if defined(REF_CACHE_NO_POOLED_ALLOC)
#  define pool_alloc(x) malloc(x->dsize)
#  define pool_free(p, x) free(x)
#endif

#define INIT_POLLED_SZ 16
#define INIT_IDX_SZ    16

struct Poll_wrap {
    pool_alloc_t  *pool;
    struct pollfd *polled;
    unsigned int   npolled;
    unsigned int   polled_sz;
    unsigned int  *fd_index;
    Pw_item      **item_index;
    unsigned int   idx_sz;
    unsigned int   last_out;
    int            need_compact;
    int            debug;
};

void pw_close(Poll_wrap *pw) {
    if (pw == NULL) return;
    if (pw->pool != NULL)       pool_destroy(pw->pool);
    if (pw->polled != NULL)     free(pw->polled);
    if (pw->fd_index != NULL)   free(pw->fd_index);
    if (pw->item_index != NULL) free(pw->item_index);
    free(pw);
}

Poll_wrap *pw_init(int debug) {
    Poll_wrap *pw = calloc(1, sizeof(Poll_wrap));
    if (pw == NULL) return NULL;

    pw->pool = pool_create(sizeof(Pw_item));
    if (pw->pool == NULL) goto fail;

    pw->polled_sz = INIT_POLLED_SZ;
    pw->polled = malloc(pw->polled_sz * sizeof(struct pollfd));
    if (pw->polled == NULL) goto fail;

    pw->idx_sz = INIT_IDX_SZ;
    pw->fd_index = malloc(pw->idx_sz * sizeof(unsigned int));
    if (pw->fd_index == NULL) goto fail;

    pw->item_index = calloc(pw->idx_sz, sizeof(Pw_item *));
    if (pw->item_index == NULL) goto fail;

    pw->debug = debug;

    return pw;
 fail:
    pw_close(pw);
    return NULL;
}

Pw_item * pw_register(Poll_wrap *pw, int fd, Pw_fd_type fd_type,
                      uint32_t init_events, void *userp) {
    Pw_item *item;

    if (pw->debug) {
        fprintf(stderr, "pw_register(%p, %d, %d, 0x%04x, %p)\n",
                (void *) pw, fd, (int) fd_type, init_events, userp);
    }

    if (fd < 0) {
        errno = EBADF;
        return NULL;
    }

    if ((unsigned) fd < pw->idx_sz && pw->item_index[fd] != NULL) {
        errno = EEXIST;
        return NULL;
    }

    if (pw->idx_sz <= (unsigned int) fd) {
        unsigned int new_sz = (unsigned int) fd + 1;
        unsigned int *new_index = realloc(pw->fd_index,
                                          new_sz * sizeof(unsigned int));
        Pw_item **new_items;
        if (new_index == NULL)
            return NULL;
        pw->fd_index = new_index;
        new_items = realloc(pw->item_index, new_sz * sizeof(Pw_item *));
        if (new_items == NULL)
            return NULL;
        memset(new_items + pw->idx_sz,0,(new_sz - pw->idx_sz) * sizeof(Pw_item *));
        pw->item_index = new_items;
        pw->idx_sz = new_sz;
    }
    if (pw->npolled == pw->polled_sz) {
        unsigned int new_sz = pw->polled_sz * 2;
        struct pollfd *new_polled = realloc(pw->polled,
                                            new_sz * sizeof(*new_polled));
        if (new_polled == NULL)
            return NULL;
        pw->polled = new_polled;
        pw->polled_sz = new_sz;
    }

    item = pool_alloc(pw->pool);
    if (item == NULL)
        return NULL;

    item->fd      = fd;
    item->fd_type = fd_type;
    item->userp   = userp;

    pw->fd_index[fd]   = pw->npolled;
    pw->item_index[fd] = item;

    pw->polled[pw->npolled].fd     = fd;
    pw->polled[pw->npolled].events = (short) init_events;
    pw->npolled++;

    return item;
}

int pw_mod(Poll_wrap *pw, Pw_item *item, uint32_t events) {
    if (pw->debug) {
        fprintf(stderr, "pw_mod(%p, %d, 0x%04x)\n",
                (void *) pw, item->fd, events);
    }

    if (item->fd < 0 || (unsigned int) item->fd >= pw->idx_sz
        || pw->item_index[item->fd] == NULL) {
        errno = ENOENT;
        return -1;
    }

    pw->polled[pw->fd_index[item->fd]].events = (short) events;
    return 0;
}

int pw_wait(Poll_wrap *pw, Pw_events *events,
            int max_events, int timeout) {
    unsigned int i;
    unsigned int j;
    unsigned int end;
    int          out = 0;
    int          res;

    /* Remove deleted items */
    if (pw->need_compact) {
        for (i = 0, j = 0; i < pw->npolled; i++) {
            if (i == pw->last_out) pw->last_out = j;
            if (pw->item_index[pw->polled[i].fd] == NULL)
                continue;
            if (i != j) {
                pw->polled[j] = pw->polled[i];
                pw->fd_index[pw->polled[j].fd] = j;
            }
            j++;
        }
        pw->need_compact = 0;
        pw->npolled = j;
    }

    /* Do the poll */
    res = poll(pw->polled, pw->npolled, out == 0 ? timeout : 0);
    if (res < 0) return res;

    /* Start copying events from where the previous pw_wait left off, to
       prevent later events from never making it to the list. */
    end = pw->last_out;
    for (; pw->last_out < pw->npolled && out < max_events; pw->last_out++) {
        if (pw->polled[pw->last_out].revents == 0)
            continue;
        events[out].events = (uint32_t) pw->polled[pw->last_out].revents;
        events[out++].item = pw->item_index[pw->polled[pw->last_out].fd];
    }
    for (pw->last_out = 0;
         pw->last_out < end && out < max_events; pw->last_out++) {
        if (pw->polled[pw->last_out].revents == 0) continue;
        events[out].events = (uint32_t) pw->polled[pw->last_out].revents;
        events[out++].item = pw->item_index[pw->polled[pw->last_out].fd];
    }

    return out;
}

int pw_remove(Poll_wrap *pw, Pw_item *item, int do_close) {
    int fd = item->fd;

    if (pw->debug) {
        fprintf(stderr, "pw_remove(%p, %d%s)\n",
                (void *) pw, item->fd, do_close ? ", close" : "");
    }

    if (item->fd < 0 || (unsigned int) item->fd >= pw->idx_sz
        || pw->item_index[item->fd] == NULL) {
        errno = ENOENT;
        return -1;
    }
    pw->item_index[item->fd] = NULL;
    pw->polled[pw->fd_index[item->fd]].events = 0;
    pw->need_compact = 1;
    pool_free(pw->pool, item);
    if (!do_close) return 0;
    return close(fd);
}
#else
// Prevent "empty translation unit" errors
const char poll_wrap_poll_disabled[] = "poll_wrap_poll not built";
#endif /* HAVE_EPOLL */
