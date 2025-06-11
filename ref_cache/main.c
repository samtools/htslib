/*  main.c -- ref-cache main entry point

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

#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
// Work around older MacOS SDKs which hid O_DIRECTORY if _POSIX_C_SOURCE
// or _XOPEN_SOURCE are defined.
#define _DARWIN_C_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <assert.h>
#ifdef HAVE_PRCTL
#include <sys/prctl.h>
#endif
#include <time.h>
#include <limits.h>

#include "listener.h"
#include "log_files.h"
#include "misc.h"
#include "options.h"
#include "ping.h"
#include "poll_wrap.h"
#include "server.h"
#include "upstream.h"

// May not have been defined, depending on compiler options
#ifndef NI_MAXHOST
#  define NI_MAXHOST 1025
#endif


typedef enum {
    CHLD_SERVER,
    CHLD_UPSTREAM,
} Child_type;

typedef struct {
    Child_type type;
    pid_t      pid;
    Pw_item   *polled_rd;
    int        log_rd;
    int        log_wr;
    int        upstream;
} Child_proc;

static int *upstream = NULL;
static Child_proc *kids = NULL;
static size_t nkids = 0;
static int sig_fds[2];
Pw_item *polled_sig = NULL;


static volatile sig_atomic_t got_chld = 0;

static char *argv0 = NULL;
static size_t argv0_len = 0;



static void change_name(char *name) {
    if (argv0_len < strlen(name)) return;
    strncpy(argv0, name, argv0_len);
#if HAVE_PRCTL
    prctl(PR_SET_NAME, (unsigned long) name, 0, 0, 0);
#endif
}

static int init_children(Options *opts) {
    int k;

    upstream = malloc((opts->max_kids + 1U) * sizeof(int));
    if (upstream == NULL) return -1;

    kids = malloc((opts->max_kids + 1U) * sizeof(Child_proc));
    if (kids == NULL) return -1;

    for (k = 0; k <= opts->max_kids; k++) {
        kids[k].type = CHLD_SERVER;
        kids[k].pid  = 0;
        kids[k].polled_rd = NULL;
        kids[k].log_rd = kids[k].log_wr = kids[k].upstream = upstream[k] = -1;
    }

    kids[opts->max_kids].type = CHLD_UPSTREAM;

    for (k = 0; k < opts->max_kids; k++) {
        int sv[2];
        int pipefd[2];

        /* Make a pipe for the log file */
        if (pipe(pipefd) != 0) {
            perror("pipe");
            goto fail;
        }

        kids[k].log_rd = pipefd[0];
        kids[k].log_wr = pipefd[1];

        /* Sockets for upstream */
        if (opts->upstream_url != NULL) {
            if (socketpair(AF_UNIX, SOCK_DGRAM, 0, sv) != 0) {
                perror("socketpair");
                goto fail;
            }

            kids[k].upstream = sv[0];
            upstream[k]      = sv[1];
        }
    }

    /* Pipe for signals */
    if (pipe(sig_fds) != 0) {
        perror("pipe");
        goto fail;
    }

    return 0;

 fail:
    for (k = 0; k < opts->max_kids; k++) {
        if (kids[k].log_rd   >= 0) close(kids[k].log_rd);
        if (kids[k].log_wr   >= 0) close(kids[k].log_wr);
        if (kids[k].upstream >= 0) close(kids[k].upstream);
        if (upstream[k] >= 0) close(upstream[k]);
    }
    return -1;
}

static int set_up_child(Options *opts, int k, int is_upstream, Poll_wrap *pw) {
    struct sigaction sigact;
    int i;

    if (is_upstream) k = opts->max_kids;

    /* Restore default signal handler */
    sigact.sa_handler = SIG_DFL;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = 0;
    if (sigaction(SIGCHLD, &sigact, NULL) != 0) {
        perror("Resetting SIGCHLD handler");
        return -1;
    }

    assert(kids != NULL);

    /* Close all the file descriptors we don't need */
    for (i = 0; i < opts->max_kids; i++) {
        close(kids[i].log_rd);
        if (!is_upstream && upstream[i] >= 0) close(upstream[i]);
        if (i != k) {
            close(kids[i].log_wr);
            if (kids[i].upstream >= 0) close(kids[i].upstream);
        }
    }
    close(sig_fds[0]);
    close(sig_fds[1]);

    if (pw != NULL) pw_close(pw);

    if (opts->log != stdout) {
        close(fileno(opts->log));
    }

    free(kids);
    if (!is_upstream) free(upstream);

    change_name(is_upstream ? "refc[dl]" : "refc[svr]");

    return 0;
}

static int make_new_child(Options *opts, Listeners *lsocks, Poll_wrap *pw) {
    int k;
    pid_t pid;

    assert(kids != NULL);

    /* Find a free slot */
    for (k = 0; k < opts->max_kids; k++) {
        if (kids[k].pid == 0) break;
    }
    assert(k < opts->max_kids);

    /* Start the child process */
    pid = fork();
    if (pid < 0) {
        perror("fork");
        return -1;
    }

    if (pid == 0) {
        /* Copy file descriptors as set_up_child frees kids[] */
        int upstr  = kids[k].upstream;
        int log_wr = kids[k].log_wr;
        int res;
        if (set_up_child(opts, k, 0, pw) != 0) _exit(EXIT_FAILURE);
        res = run_poll_loop(opts, lsocks, upstr, log_wr);
        _exit(res == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
    }

    kids[k].pid = pid;
    nkids++;
    return 0;
}

int start_upstream(Options *opts, Poll_wrap *pw) {
    pid_t upstream_pid;
    int liveness_pipe[2] = { -1, -1 };

    // Make pipe so child can detect parent going away
    if (pipe(liveness_pipe) < 0) {
        perror("Opening pipe");
        return -1;
    }

    upstream_pid = fork();
    if (upstream_pid == -1) {
        perror("start_upstream couldn't fork");
        close(liveness_pipe[0]);
        close(liveness_pipe[1]);
        return -1;
    }

    if (upstream_pid == 0) {
        int res;
        close(liveness_pipe[0]);
        if (set_up_child(opts, 0, 1, pw) != 0) _exit(EXIT_FAILURE);

        res = run_upstream_handler(opts, upstream, liveness_pipe[1]);

        _exit(res == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
    }

    assert(kids != NULL);

    close(liveness_pipe[1]);
    kids[opts->max_kids].pid = upstream_pid;
    kids[opts->max_kids].type = CHLD_UPSTREAM;
    kids[opts->max_kids].log_rd = liveness_pipe[0];

    return 0;
}

static void sig_handler(int signal) {
    ssize_t bytes;
    char c;

    switch (signal) {
    case SIGCHLD:
        if (got_chld) return; /* save writing repeatedly to the pipe */
        got_chld = 1;
        break;
    default:
        return;
    }

    if (sig_fds[1] < 0) return;
    c = '*';
    do {
        bytes = write(sig_fds[1], &c, 1);
    } while (bytes < 0 && errno == EINTR);
    if (bytes < 0 && (errno != EAGAIN && errno != EWOULDBLOCK)) {
        close(sig_fds[1]); /* Should get the attention of the other end... */
    }
}

static int handle_sigchld(Options *opts, Poll_wrap *pw) {
    char buffer[16];
    ssize_t bytes;
    pid_t pid;

    /* Drain the pipe */
    do {
        bytes = read(sig_fds[0], buffer, sizeof(buffer));
    } while (bytes < 0 && errno == EINTR);
    if (bytes <= 0) {
        fprintf(stderr, "%s on signal fd #%d\n",
                bytes == 0 ? "EOF" : strerror(errno), sig_fds[0]);
        return -1;
    }

    got_chld = 0; /* Allow more writes to the pipe */

    /* Reap the children */
    do {
        int status;
        pid = waitpid(-1, &status, WNOHANG);
        if (pid < 0) {
            if (errno == ECHILD || errno == EINTR) continue;
            perror("waitpid");
            return -1;
        }

        if (pid > 0) {
            int i;

            if (WIFEXITED(status)) {
                fprintf(stderr, "Child PID %d exited with status %d.\n",
                        pid, WEXITSTATUS(status));
            } else if (WIFSIGNALED(status)) {
                fprintf(stderr, "Child PID %d terminated by signal %d.\n",
                        pid, WTERMSIG(status));
            } else {
                fprintf(stderr, "Child PID %d terminated\n", pid);
            }
            for (i = 0; i < opts->max_kids + 1; i++) {
                if (kids[i].pid == pid) {
                    kids[i].pid = 0;
                    if (kids[i].type == CHLD_UPSTREAM) {
                        if (start_upstream(opts, pw) != 0) return -1;
                    } else {
                        --nkids;
                    }
                    break;
                }
            }
        }
    } while (pid > 0);

    return 0;
}

static Poll_wrap * init_poller(Options *opts) {
    Poll_wrap *pw = pw_init(opts->verbosity > 2);
    size_t k;

    if (pw == NULL) return NULL;

    for (k = 0; k < opts->max_kids; k++) {
        kids[k].polled_rd = pw_register(pw, kids[k].log_rd, MAIN_LOG_RD,
                                        PW_IN|PW_ERR|PW_HUP, &kids[k]);
        if (kids[k].polled_rd == NULL) {
            perror("Registering log pipe with poller");
            goto fail;
        }
    }

    polled_sig = pw_register(pw, sig_fds[0], MAIN_SIG, PW_IN|PW_HUP|PW_ERR, NULL);
    if (polled_sig == NULL) {
        perror("Registering signal pipe with poller");
        goto fail;
    }

    return pw;
 fail:

    for (k = 0; k < nkids; k++) {
        if (kids[k].polled_rd != NULL) pw_remove(pw, kids[k].polled_rd, 0);
    }
    return NULL;
}

#define MAX_EVENTS 16
static int run_server_population(Options *opts, Listeners *lsocks,
                                 Logfiles *logfiles) {
    struct sigaction sigact;
    Poll_wrap *pw;
    int logged = 0;

    /* Set up signal handler */
    sigact.sa_handler = sig_handler;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = SA_NOCLDSTOP;
    if (sigaction(SIGCHLD, &sigact, NULL) != 0) {
        perror("Setting up SIGCHLD handler");
        return -1;
    }

    pw = init_poller(opts);
    if (pw == NULL) return -1;

    /* Run poll loop */
    for (;;) {
        int e;
        int ret;
        Pw_events events[MAX_EVENTS];

        /* Make new children if necessary */
        while (nkids < opts->max_kids) {
            if (make_new_child(opts, lsocks, pw) != 0) {
                if (nkids == 0) {
                    fprintf(stderr, "Unable to make server processes, giving up.\n");
                    return -1;
                } else {
                    break;
                }
            }
        }

        /* Wait for events */
        ret = pw_wait(pw, events, MAX_EVENTS, logged ? 100 : -1);
        if (ret < 0) {
            if (errno != EINTR) {
                perror("Waiting for poller");
                return -1;
            }
            continue;
        }

        if (ret == 0 && logged) {
            /* Flush the log */
            fflush(opts->log);
            logged = 0;
            continue;
        }

        for (e = 0; e < ret; e++) {
            Pw_item *item = PWI(events[e]);
            assert(item != NULL);

            switch (item->fd_type) {

            case MAIN_SIG:
                /* Got a signal */
                if (handle_sigchld(opts, pw) != 0) return -1;
                break;

            case MAIN_LOG_RD: {
                /* Deal with log messages */
                char buffer[65536];
                ssize_t bytes;
                Child_proc *kid = (Child_proc *) item->userp;
                assert(kid != NULL);
                assert(kid->pid != 0);

                do {
                    bytes = read(kid->log_rd, buffer, sizeof(buffer));
                } while (bytes < 0 && errno == EINTR);

                if (bytes <= 0) {
                    fprintf(stderr, "%s reading log fd #%d\n",
                            bytes == 0 ? "EOF" : strerror(errno), kid->log_rd);
                    close(kid->log_rd);
                    kid->log_rd = -1;
                } else {
                    if (!opts->no_log) {
                        if (write_to_log(logfiles, opts,
                                         buffer, (size_t) bytes) < 0)
                            return -1;
                        if (opts->log_dir != NULL)
                            logged = 1;
                    }
                }
                break;
            }
            default:
                fprintf(stderr, "Unexpected item type %x from poll\n", item->fd_type);
                abort();
            }
        }
    }
}

static int daemonise(int daemon_fds[2], Options *opts) {
    const char *error_log_file = opts->error_log_file;
    int i, save_errno, pid1, pid2, devnull;
    struct rlimit fd_limit;
    struct sigaction sigact;
    sigset_t all_sigs;

    if (getrlimit(RLIMIT_NOFILE, &fd_limit) != 0) {
        perror("Getting max file descriptor count");
        return -1;
    }

    // Close any open file descriptors above 3 (apart from the pipe passed in)
    save_errno = errno;
    for (i = 3; i < (int) fd_limit.rlim_cur; i++) {
        if (i != daemon_fds[0] && i != daemon_fds[1])
            close(i);
    }

    // Reset signal handlers to default
    // Unfortunately there's no portable way to count them, so take
    // a punt that it's going to be less than sizeof(sigset_t) * 8
    sigfillset(&all_sigs);
    sigdelset(&all_sigs, SIGKILL);  // Can't be changed
    sigdelset(&all_sigs, SIGSTOP);  // Can't be changed
    for (i = 1; i < (int) sizeof(all_sigs) * 8; i++) {
        if (!sigismember(&all_sigs, i))
            continue;
        memset(&sigact, 0, sizeof(sigact));
        sigact.sa_handler = SIG_DFL;
        sigemptyset(&sigact.sa_mask);
        if (sigaction(i, &sigact, NULL) != 0) {
            // Ideally, sigfillset() would only fill valid signals,
            // but sadly on Linux, at least, that doesn't appear to be
            // the case so we should expect to get EINVAL
            if (errno != EINVAL) {
                fprintf(stderr, "Resetting signal handler %d: %s\n",
                        i, strerror(errno));
                return -1;
            } else {
                break;
            }
        }
    }
    errno = save_errno;

    // Unblock all signals
    sigemptyset(&all_sigs);
    if (sigprocmask(SIG_SETMASK, &all_sigs, NULL) != 0) {
        perror("Setting signal mask");
        return -1;
    }

    pid1 = fork();
    if (pid1 < 0) {
        perror("Couldn't fork");
        return -1;
    }
    if (pid1 != 0) {
        // Check that the daemon started up successfully, by waiting for
        // a message from the pipe in daemon_fds[].
        char msg[2] = { '\0', '\0' };

        // Ensure a clean exit from initial process
        free(opts->match_addrs);

        close(daemon_fds[1]); // Writing end is for the daemon process
        ssize_t res = do_read_all(daemon_fds[0], &msg, sizeof(msg));
        close(daemon_fds[0]);
        if (res < 0 || !(msg[0] == 'o' && msg[1] == 'k')) {
            fprintf(stderr, "Daemon failed to start, sorry.\n");
            _exit(EXIT_FAILURE);
        }
        if (!opts->error_log_file && !opts->no_log) {
            // Last chance to tell the user...
            fprintf(stderr,
                    "Error messages will be unavailable after this one.\n");
        }
        _exit(EXIT_SUCCESS);
    }

    // Close the reading end of the pipe.  Writing end will be closed on
    // successful start.
    close(daemon_fds[0]);
    daemon_fds[0] = -1;

    // Become session leader
    if (setsid() < 0) {
        perror("Couldn't become session leader");
        _exit(EXIT_FAILURE);
    }

    // Fork again to make the daemon
    pid2 = fork();
    if (pid2 < 0) {
        perror("Couldn't fork");
        _exit(EXIT_FAILURE);
    }
    if (pid2 != 0) {
        // Exit so the daemon is inherited by PID 1
        exit(EXIT_SUCCESS);
    }

    // Redirect stdin and stdout to /dev/null
    devnull = open("/dev/null", O_RDWR);
    if (devnull < 0) {
        perror("Opening /dev/null");
        exit(EXIT_FAILURE);
    }
    if (dup2(devnull, 0) < 0
        || dup2(devnull, 1) < 0
        || dup2(devnull, 2) < 0) {
        perror("Redirecting stdin and stdout to /dev/null");
        exit(EXIT_FAILURE);
    }

    // Redirect stderr to error_log_file, or /dev/null if not set
    if (error_log_file != NULL) {
        if (freopen(error_log_file, "a", stderr) == NULL) {
            /* Does it still exist? */
            fprintf(stderr, "Couldn't redirect stderr to %s: %s\n",
                    error_log_file, strerror(errno));
            exit(EXIT_FAILURE);
        }
    } else {
        if (dup2(devnull, 2) < 0) {
            perror("Redirecting stderr to /dev/null");
            exit(EXIT_FAILURE);
        }
    }

    if (devnull > 2)
        close(devnull);

    // Reset umask
    umask(0);

    return 0;
}

static int get_systemd_listen_fds(Options *opts) {
    struct rlimit fd_limit;
    pid_t our_pid = getpid();
    char *env_listen_pid = getenv("LISTEN_PID");
    char *env_listen_fds = getenv("LISTEN_FDS");
    char *end = env_listen_pid;
    long long listen_pid;
    long listen_fds;
    if (env_listen_pid == NULL) {
        fprintf(stderr, "LISTEN_PID environment variable is not set\n");
        return -1;
    }
    if (env_listen_fds == NULL) {
        fprintf(stderr, "LISTEN_FDS environment variable is not set\n");
        return -1;
    }
    if (getrlimit(RLIMIT_NOFILE, &fd_limit) != 0) {
        perror("Getting max file descriptor count");
        return -1;
    }

    listen_pid = strtoll(env_listen_pid, &end, 10);
    if (*env_listen_pid == '\0' || *end != '\0'
        || listen_pid < 0 || listen_pid != (long long) our_pid) {
        fprintf(stderr, "LISTEN_PID is incorrect\n");
        return -1;
    }
    end = env_listen_fds;
    listen_fds = strtol(env_listen_fds, &end, 10);
    if (*env_listen_fds == '\0' || *end != '\0'
        || listen_fds <= 0 || listen_fds > INT_MAX
        || listen_fds > (long) fd_limit.rlim_cur - FIRST_SD_LISTEN_FD) {
        fprintf(stderr, "LISTEN_FDS is not valid\n");
        return -1;
    }
    opts->listen_fds = (int) listen_fds;
    unsetenv("LISTEN_PID");
    unsetenv("LISTEN_FDS");
    unsetenv("LISTEN_FDNAMES");
    return 0;
}

// Append supplied IP address ranges to the list we will accept
// connections from.
static int add_match_addr(Options *opts, const char *addr_list) {
    char host[NI_MAXHOST + 1];
    unsigned long netmask_bits;
    size_t host_start, host_len, netmask_end, p;
    size_t addr_list_len = strlen(addr_list);
    struct addrinfo hints, *addrs, *addr;
    int res;

    for (host_start = p = 0; host_start < addr_list_len;
         host_start = netmask_end) {
        // Look for the IP address part
        size_t p = host_start;
        host_len = strcspn(addr_list + p, "/,");
        if (host_len >= NI_MAXHOST) {
            fprintf(stderr, "IP address \"%.20s...\" too long\n",
                    addr_list + host_start);
            return -1;
        }
        memcpy(host, addr_list + p, host_len);
        host[host_len] = '\0';
        p += host_len;

        // Check for CIDR-notation netmask
        if (addr_list[p] == '/') {
            size_t netmask_len = strcspn(addr_list + p, ",");
            const char *nm = addr_list + p + 1;
            char *endp = NULL;
            netmask_bits = strtoul(nm, &endp, 10);
            if (endp == nm) {
                fprintf(stderr, "Empty netmask for host \"%s\"\n", host);
                return -1;
            }
            if (netmask_bits > 128) {
                fprintf(stderr, "Netmask \"%s/%lu\" too large\n",
                        host, netmask_bits);
                return -1;
            }
            netmask_end = p + netmask_len;
        } else {
            netmask_bits = 128;
            netmask_end = p;
        }
        if (netmask_end < addr_list_len && addr_list[netmask_end] == ',')
            netmask_end++;

        // Use getaddrinfo() to convert the address to a struct sockaddr
        addrs = NULL;
        memset(&hints, 0, sizeof(hints));
        hints.ai_family = AF_UNSPEC;
        hints.ai_socktype = SOCK_STREAM;
        hints.ai_protocol = 0;
        hints.ai_canonname = NULL;
        hints.ai_addr = NULL;
        hints.ai_next = NULL;
        hints.ai_flags = AI_NUMERICHOST;
        res = getaddrinfo(host, NULL, &hints, &addrs);
        if (res != 0) {
            fprintf(stderr, "Couldn't resolve IP address \"%s\" : %s\n",
                    host, gai_strerror(res));
            return -1;
        }
        // In theory there could be more than one address, so iterate just
        // in case.
        for (addr = addrs; addr != NULL; addr = addr->ai_next) {
            MatchAddr *ma;
            size_t alen;
            uint8_t *addrp;
            unsigned long max_mask_bits;

            // Ignore unexpected types
            if (addr->ai_family != AF_INET && addr->ai_family != AF_INET6)
                continue;

            // Grow output array if necessary
            if (opts->match_addrs_size == opts->num_match_addrs) {
                size_t new_size = (opts->match_addrs_size > 0
                                   ? opts->match_addrs_size * 2 : 16);
                MatchAddr *new_addrs = realloc(opts->match_addrs,
                                               new_size * sizeof(*new_addrs));
                if (!new_addrs) {
                    perror("Allocating address match list");
                    freeaddrinfo(addrs);
                    return -1;
                }
                opts->match_addrs_size = new_size;
                opts->match_addrs = new_addrs;
            }

            ma = &opts->match_addrs[opts->num_match_addrs];
            memset(ma, 0, sizeof(*ma));
            ma->family = addr->ai_addr->sa_family;

            if (addr->ai_family == AF_INET
                && addr->ai_addrlen == sizeof(struct sockaddr_in)) {
                alen = sizeof(struct in_addr);
                addrp = ((uint8_t *) addr->ai_addr) + offsetof(struct sockaddr_in, sin_addr);
            } else if (addr->ai_family == AF_INET6
                       && addr->ai_addrlen == sizeof(struct sockaddr_in6)) {
                alen = sizeof(struct in6_addr);
                addrp = ((uint8_t *) addr->ai_addr) + offsetof(struct sockaddr_in6, sin6_addr);
            } else {
                fprintf(stderr,
                        "Unexpected address type/length! "
                        "Got %u/%zd expected %u/%zd or %u/%zd\n",
                        (unsigned) addr->ai_addr->sa_family,
                        (size_t) addr->ai_addrlen,
                        (unsigned) AF_INET,  sizeof(struct sockaddr_in),
                        (unsigned) AF_INET6, sizeof(struct sockaddr_in6));
                freeaddrinfo(addrs);
                return -1;
            }
            max_mask_bits = alen * 8;

            // Store netmask
            if (netmask_bits > max_mask_bits)
                netmask_bits = max_mask_bits;

            ma->mask_bytes = (uint8_t) (netmask_bits / 8);
            if (ma->mask_bytes == alen) {
                // This ensures ma->addr[ma->mask_bytes] is always valid
                // which makes checking easier
                --ma->mask_bytes;
                ma->mask = 0xffU;
            } else {
                ma->mask = (uint8_t) ((0xff00U >> (netmask_bits & 7U)) & 0xffU);
            }

            // Copy masked IP address
            memcpy(ma->addr, addrp, ma->mask_bytes);
            ma->addr[ma->mask_bytes]
                = (uint8_t) (addrp[ma->mask_bytes] & ma->mask);

            opts->num_match_addrs++;
        }
        freeaddrinfo(addrs);
    }
    return 0;
}

static int cmp_match_addrs(const void *av, const void *bv) {
    const MatchAddr *a = (const MatchAddr *) av;
    const MatchAddr *b = (const MatchAddr *) bv;
    int ip_order = AF_INET < AF_INET6 ? 1 : -1;

    if (a->family != b->family)
        return ip_order * (a->family < b->family ? -1 : 1);
    return memcmp(a->addr, b->addr, sizeof(a->addr));
}

static void sort_match_addrs(Options *opts) {
    size_t i, j;
    int ip6_seen = 0;

    if (!opts->match_addrs)
        return;

    qsort(opts->match_addrs, opts->num_match_addrs,
          sizeof(opts->match_addrs[0]), cmp_match_addrs);

    for (i = j = 0; j < opts->num_match_addrs;) {
        if (!ip6_seen && opts->match_addrs[i].family == AF_INET6) {
            opts->first_ip6 = i;
            ip6_seen = 1;
        }
        do {
            ++j;
        } while (j < opts->num_match_addrs
                 && memcmp(&opts->match_addrs[i], &opts->match_addrs[j],
                           sizeof(opts->match_addrs[0])) == 0);
        ++i;
        if (j < opts->num_match_addrs && i < j) {
            memcpy(&opts->match_addrs[i], &opts->match_addrs[j],
                   sizeof(opts->match_addrs[0]));
        }
    }
    opts->num_match_addrs = i;
    if (!ip6_seen)
        opts->first_ip6 = i;
}

static void usage(const char *prog, int help, const Options *opts) {
    fprintf(stderr, "Usage: %s [options] -d <dir>\n", prog);
    if (help) {
        fprintf(stderr,
                "Options:\n"
                "  -b         Run in background as a daemon\n"
                "  -d <dir>   Directory for cached reference files\n"
                "  -h         Show help\n"
                "  -l <dir>   Directory for log files.  Log to stdout if not set and running in\n"
                "             foreground\n"
                "  -L         Don't log\n"
                "  -m <list>  Only respond to connections from these networks\n"
                "  -n <1-4>   Number of server processes to run [%u]\n"
                "  -p <num>   Port number to listen on [%u]\n"
                "  -r <num>   Number of request log files to keep [%u]\n"
                "  -R <num>   Maximum size of a single request log file (MiB) [%lld]\n"
                "  -s         Run as a systemd socket service\n"
                "  -u <url>   URL for upstream server%s%s%s\n"
                "  -U         Only serve local files, turn off upstream\n"
                "  -v         Turn on debugging output\n",
                opts->max_kids,
                opts->port,
                opts->nlogs,
                (long long) opts->max_log_sz>>20,
                (opts->upstream_url
                 ? (opts->upstream_url_len > 40 ? "\n             [" : " [")
                 : ""),
                opts->upstream_url ? opts->upstream_url : "",
                opts->upstream_url ? "]" : "");
    }
}

static inline uint16_t get_opt_val(const char *arg, const char *prog,
                                   const char *opt,
                                   uint16_t min, uint16_t max,
                                   int *badarg) {
    char *endp = NULL;
    long lmin = (long) min;
    long lmax = (long) max;
    long val = strtol(arg, &endp, 0);
    if (*arg == '\0' || *endp != '\0') {
        fprintf(stderr, "%s : %s option value \"%s\" is not a number\n",
                prog, opt, arg);
        *badarg = 1;
        return 0;
    }
    if (val < lmin || val > lmax) {
        fprintf(stderr,
                "%s : %s option value \"%s\" should be between %ld and %ld\n",
                prog, opt, arg, lmin, lmax);
        *badarg = 1;
        return 0;
    }
    return (uint16_t) val;
}

int main(int argc, char **argv) {
    Options opts;
    Listeners *lsocks = NULL;
    Logfiles *logfiles = NULL;
    int c, res, badarg = 0, retval = EXIT_FAILURE, show_help = 0;
    int daemon_pipe[2] = { -1, -1 };
    const char *ip_ranges_all = "0.0.0.0/0,::/0";
    const char *ip_ranges_localhost = "127.0.0.0/8,::1/128";
    const char *ip_ranges_default = "10.0.0.0/8,172.16.0.0/12,192.168.0.0/16,fc00::/7,fe80::/10";

    /* Copy argv[0] for change_name() */
    argv0 = argv[0];
    argv0_len = strlen(argv0);

    /* Options */

    memset(&opts, 0, sizeof(opts));

    opts.port = 8000;
    opts.cache_dir = NULL;
    opts.cache_fd  = -1;
    opts.max_kids = 1;
    opts.verbosity = 0;
    opts.log_dir = NULL;
    opts.log = stdout;
    opts.match_addrs = NULL;
    opts.num_match_addrs = 0;
    opts.match_addrs_size = 0;
    opts.error_log_file = NULL;
    opts.nlogs = 5;
    opts.max_log_sz = 10<<20;
    opts.daemon = not_a_daemon;
    opts.no_log = 0;
    opts.upstream_url = "https://www.ebi.ac.uk/ena/cram/md5/";

    while ((c = getopt(argc, argv, "be:d:hl:Lm:n:p:r:R:su:Uv")) >= 0) {
        switch (c) {
        case 'b':
            if (opts.daemon != not_a_daemon) {
                fprintf(stderr,
                        "The -b and -s options cannot be used together\n");
                badarg = 1;
            } else {
                opts.daemon = sysv_daemon;
            }
            break;
        case 'e': opts.error_log_file = optarg; break;
        case 'd': opts.cache_dir = optarg; break;
        case 'h': show_help = 1; break;
        case 'l':
            opts.log_dir = optarg;
            opts.no_log = 0;
            break;
        case 'L':
            if (!opts.log_dir)
                opts.no_log = 1;
            break;
        case 'm': {
            const char *to_add = optarg;
            if (strcmp(optarg, "all") == 0) {
                to_add = ip_ranges_all;
            } else if (strcmp(optarg, "default") == 0) {
                to_add = ip_ranges_default;
            } else if (strcmp(optarg, "localhost") == 0) {
                to_add = ip_ranges_localhost;
            }
            if (add_match_addr(&opts, to_add) != 0) {
                goto cleanup;
            }
            break;
        }
        case 'n':
            opts.max_kids = get_opt_val(optarg, argv[0], "-n", 1, 4, &badarg);
            break;
        case 'p':
            opts.port = get_opt_val(optarg, argv[0], "-p", 1, 65535, &badarg);
            break;
        case 'r':
            opts.nlogs = get_opt_val(optarg, argv[0], "-r", 1, 100, &badarg);
            break;
        case 'R':
            opts.max_log_sz = (off_t) get_opt_val(optarg, argv[0], "-R",
                                                  1, 1000, &badarg) * (1<<20);
            break;
        case 's':
            if (opts.daemon != not_a_daemon) {
                fprintf(stderr,
                        "The -b and -s options cannot be used together\n");
                badarg = 1;
            } else {
                opts.daemon = systemd_socket_service;
            }
            break;
        case 'u': opts.upstream_url = optarg; break;
        case 'U': opts.upstream_url = NULL; break;
        case 'v': opts.verbosity++; break;
        default:
            usage(argv[0], 0, &opts);
            goto cleanup;
        }
    }

    if (show_help) {
        usage(argv[0], 1, &opts);
        retval = EXIT_SUCCESS;
        goto cleanup;
    }

    if (badarg)
        goto cleanup;

    if (opts.cache_dir == NULL) {
        usage(argv[0], 0, &opts);
        goto cleanup;
    }

    if (opts.match_addrs == NULL) {
        if (add_match_addr(&opts, ip_ranges_default) != 0)
            goto cleanup;
    }
    if (add_match_addr(&opts, ip_ranges_localhost) != 0)
        goto cleanup;

    sort_match_addrs(&opts);

    if (opts.daemon != systemd_socket_service) {
        /* See if we're already running */
        res = check_running(opts.port);
        if (res != 0) {
            retval = res < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
            goto cleanup;
        }
    }

    /* Turn into a daemon, if requested */
    switch (opts.daemon) {
    case sysv_daemon:
        if (!opts.log_dir && !opts.no_log) {
            fprintf(stderr,
                    "Warning: Running as a daemon without setting a request "
                    "log directory!\nRequest logs will be unavailable.\n");
        }
        if (!opts.error_log_file && !opts.no_log) {
            fprintf(stderr,
                    "Warning: Running as a daemon without setting an error "
                    "log file!\n");
        }

        if (pipe(daemon_pipe) < 0) {
            perror("Opening pipe");
            goto cleanup;
        }
        if (daemonise(daemon_pipe, &opts) != 0) {
            goto cleanup;
        }
        break;

    case systemd_socket_service:
        if (get_systemd_listen_fds(&opts) != 0)
            goto cleanup;
        break;

    default: // not_a_daemon
        break;
    }

    /* Log files */
    logfiles = open_logs(&opts);
    if (logfiles == NULL)
        goto cleanup;

    /* Open cache directory */
    opts.cache_fd = open(opts.cache_dir, O_RDONLY|O_DIRECTORY);
    if (opts.cache_fd < 0) {
        fprintf(stderr, "Couldn't open directory %s: %s\n",
                opts.cache_dir, strerror(errno));
        goto cleanup;
    }

    /* Allocate data structures and plumbing for child processes */

    if (init_children(&opts) != 0)
        goto cleanup;

    /* Start the upstream handler if needed */

    if (opts.upstream_url != NULL) {
        opts.upstream_url_len = strlen(opts.upstream_url);
        if (start_upstream(&opts, NULL) != 0)
            goto cleanup;
    }

    /* Open the socket to listen on */

    switch (opts.daemon) {
    case systemd_socket_service:
        lsocks = adopt_listen_sockets(FIRST_SD_LISTEN_FD, opts.listen_fds);
        break;
    default:
        lsocks = get_listen_sockets(opts.port);
        break;
    }
    if (lsocks == NULL) {
        fprintf(stderr, "Couldn't start up.  Sorry.\n");
        goto cleanup;
    }

    if (opts.daemon == sysv_daemon) {
        if (do_write_all(daemon_pipe[1], "ok", 2) < 0) {
            perror("Couldn't report successful start back to parent");
            goto cleanup;
        }
        close(daemon_pipe[1]);
        daemon_pipe[1] = -1;
    } else if (opts.error_log_file) {
        if (freopen(opts.error_log_file, "a", stderr) == NULL) {
            /* Does it still exist? */
            fprintf(stderr, "Couldn't redirect stderr to %s: %s\n",
                    opts.error_log_file, strerror(errno));
            goto cleanup;
        }
    }

    /* Run the servers */

    res = run_server_population(&opts, lsocks, logfiles);
    if (res == -1) {
        fprintf(stderr, "Server died.\n");
    } else {
        retval = EXIT_SUCCESS;
    }

    /* Clean up and finish */

 cleanup:
    if (lsocks)
        close_listen_sockets(lsocks);
    if (opts.cache_fd >= 0)
        close(opts.cache_fd);
    if (logfiles)
        close_logs(logfiles);

    free(opts.match_addrs);

    return retval;
}
