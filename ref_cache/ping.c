/*  ping.c -- check if ref-cache is already running

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
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

#include "ping.h"
#include "misc.h"
#include "options.h"

int check_running(int port_num) {
    char port[16];
    char reply[256];
    const char *request = "GET /hello HTTP/1.0\r\n\r\n";
    const char *expected = "HTTP/1.0 200 OK\r\n";
    const char *fails[3] = {
        "No addresses returned by getaddrinfo.",
        "Error getting socket: ",
        "Couldn't connect socket: "
    };
    size_t len = strlen(expected);
    struct addrinfo hints;
    struct addrinfo *addrs, *addr;
    int res, s, where = 0, got_connection = 0;
    ssize_t bytes;

    snprintf(port, sizeof(port), "%d", port_num);
    memset(&hints, 0, sizeof(hints));
    hints.ai_family   = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags    = AI_NUMERICSERV;
#if defined(HAVE_DECL_AI_V4MAPPED) && HAVE_DECL_AI_V4MAPPED
    hints.ai_flags    |= AI_V4MAPPED;
#endif
#if defined(HAVE_DECL_AI_ADDRCONFIG) && HAVE_DECL_AI_ADDRCONFIG
    hints.ai_flags    |= AI_ADDRCONFIG;
#endif
    res = getaddrinfo(NULL, port, &hints, &addrs);
    if (res != 0) {
        fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(res));
        return -1;
    }

    for (addr = addrs; addr != NULL; addr = addr->ai_next) {
        where = 0;
        do {
            s = socket(addr->ai_family, addr->ai_socktype, addr->ai_protocol);
        } while (s == -1 && errno == EINTR);
        if (s == -1) {
            where = 1;
            continue;
        }
        do {
            res = connect(s, addr->ai_addr, addr->ai_addrlen);
        } while (res == -1 && errno == EINTR);
        if (res == 0) {
            got_connection = 1;
            break;
        }

        where = 2;
        close(s);
    }

    freeaddrinfo(addrs);

    if (!got_connection) {
        // This is expected if the server isn't running yet.
        if (where == 2 && errno == ECONNREFUSED)
            return 0;

        fprintf(stderr, "%s%s\n", fails[where],
                where > 0 ? strerror(errno) : "");
        return -1;
    }

    // Try sending a simple request and check what comes back.
    // If it's the expected answer, the server is likely already running.
    if (do_write_all(s, request, strlen(request)) != 0) {
        fprintf(stderr, "Writing to socket: %s\n", strerror(errno));
        close(s);
        return -1;
    }

    bytes = do_read_all(s, reply, sizeof(reply));
    close(s);
    if (bytes < 0) {
        fprintf(stderr, "Reading from socket: %s\n", strerror(errno));
        return -1;
    }

    if ((size_t) bytes < len || memcmp(reply, expected, len) != 0) {
        fprintf(stderr, "Unexpected reply from server:\n%.*s\n",
                (int) bytes, reply);
        return -1;
    }
    return 1;
}
