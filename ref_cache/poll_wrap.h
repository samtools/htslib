/*  poll_wrap.h -- ref-cache wrapper around poll()-like interfaces

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

#ifndef POLL_WRAP_H
#define POLL_WRAP_H

#include <stddef.h>
#include <stdint.h>

#if defined(HAVE_EPOLL) && HAVE_EPOLL
#  include <sys/epoll.h>
#else
#  include <poll.h>
#endif

#include "types.h"

typedef enum {
    SV_LISTENER = 1,
    SV_UPSTREAM,
    SV_LOG,
    SV_CLIENT,
    SV_US_PIPE,
    US_COMMAND = 0x10000,
    US_CURL,
    US_PIPE,
    US_LIVE,
    MAIN_LOG_RD = 0x20000,
    MAIN_SIG
} Pw_fd_type;

#if defined(HAVE_EPOLL) && HAVE_EPOLL

#  define PW_IN  EPOLLIN
#  define PW_OUT EPOLLOUT
#  define PW_ERR EPOLLERR
#  define PW_HUP EPOLLHUP
#  if defined(PW_HAVE_EDGE) && PW_HAVE_EDGE // Using configure
#    define PW_EDGE EPOLLET
#  elif !defined(PW_HAVE_EDGE) && defined(EPOLLET)
#    define PW_EDGE EPOLLET
#    define PW_HAVE_EDGE 1
#  else
#    define PW_EDGE 0
#  endif /* EPOLLET */

typedef struct epoll_event Pw_events;
#  define PWE(E) (E.events)
#  define PWI(E) ((Pw_item *) (E.data.ptr))

#else /* HAVE_EPOLL */

#  define PW_IN  POLLIN
#  define PW_OUT POLLOUT
#  define PW_ERR POLLERR
#  define PW_HUP POLLHUP
#  define PW_EDGE 0

typedef struct {
    uint32_t events;
    Pw_item *item;
} Pw_events;
#  define PWE(E) (E.events)
#  define PWI(E) (E.item)

#endif /* HAVE_EPOLL */

struct Pw_item {
    int fd;
    Pw_fd_type fd_type;
    void *userp;
};

Poll_wrap *pw_init(int debug);
void pw_close(Poll_wrap *pw);
Pw_item * pw_register(Poll_wrap *pw, int fd, Pw_fd_type fd_type,
                      uint32_t init_events, void *userp);
int pw_mod(Poll_wrap *pw, Pw_item *item, uint32_t events);
int pw_wait(Poll_wrap *pw, Pw_events *events, int max_events, int timeout);
int pw_remove(Poll_wrap *pw, Pw_item *item, int do_close);

#endif /* POLL_WRAP_H */
