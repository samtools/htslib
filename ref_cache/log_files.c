/*  log_files.c -- ref-cache log file handling

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

#if defined(__NetBSD__) && !defined(_NETBSD_SOURCE)
// Needed for dirfd()
#define _NETBSD_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "log_files.h"
#include "options.h"

#define LOG_NAME_LEN 80
typedef struct {
    char name[LOG_NAME_LEN];
    off_t size;
} Logfile;

struct Logfiles {
    DIR        *dir_handle;
    FILE       *curr_log;
    Logfile    *logs;
    size_t      nlogs;
    size_t      sz;
    int         log_dir_fd;
};

static int log_compare(const void *av, const void *bv) {
    const Logfile *a = (const Logfile *) av;
    const Logfile *b = (const Logfile *) bv;

    return strcmp(b->name, a->name);
}

static int rotate_logs(Logfiles *logfiles, const Options *opts) {
    size_t i;
    char name[LOG_NAME_LEN];
    time_t now = time(NULL);
    struct tm *gmt = gmtime(&now);
    int log_fd = -1;
    FILE *file = NULL;

    /* Open a new log file */
    for (i = 0; i < 99; i++) {
        snprintf(name, LOG_NAME_LEN,
                 "ref_cache_%04d%02d%02d%02d%02d%02d_%02u.log",
                 gmt->tm_year + 1900, gmt->tm_mon + 1, gmt->tm_mday,
                 gmt->tm_hour, gmt->tm_min, gmt->tm_sec, (unsigned int) i);
        do {
            log_fd = openat(logfiles->log_dir_fd, name, O_WRONLY | O_CREAT | O_EXCL, 0644);
        } while (log_fd < 0 && errno == EINTR);
        if (log_fd >= 0) break;
        if (errno != EEXIST) break;
    }

    if (log_fd < 0) {
        fprintf(stderr, "Couldn't open %s/%s for writing: %s\n",
                opts->log_dir, name, strerror(errno));
        return -1;
    }

    file = fdopen(log_fd, "w");
    if (file == NULL) {
        fprintf(stderr, "Couldn't fdopen %s/%s : %s\n",
                opts->log_dir, name, strerror(errno));
        close(log_fd);
        unlinkat(logfiles->log_dir_fd, name, 0);
        return -1;
    }

    /* Remove stale logs */
    assert(opts->nlogs > 0);
    assert(logfiles->sz > opts->nlogs);
    for (i = opts->nlogs - 1U; i < logfiles->nlogs; i++) {
        if (unlinkat(logfiles->log_dir_fd, logfiles->logs[i].name, 0) != 0) {
            fprintf(stderr, "Warning: Couldn't remove old log file %s/%s: %s\n",
                    opts->log_dir, logfiles->logs[i].name, strerror(errno));
        }
    }

    if (logfiles->nlogs > (size_t) opts->nlogs - 1U)
        logfiles->nlogs = opts->nlogs - 1U;

    /* Put the new file at the front of the list */
    if (logfiles->nlogs > 0) {
        memmove(&logfiles->logs[1], &logfiles->logs[0],
                logfiles->nlogs * sizeof(Logfile));
    }
    memcpy(logfiles->logs[0].name, name, LOG_NAME_LEN);
    logfiles->logs[0].size = 0;
    logfiles->nlogs++;

    if (logfiles->curr_log != NULL && logfiles->curr_log != stdout)
        fclose(logfiles->curr_log);
    logfiles->curr_log = file;

    return 0;
}

void close_logs(Logfiles *logfiles) {
    if (logfiles == NULL)
        return;

    if (logfiles->dir_handle != NULL) {
        closedir(logfiles->dir_handle);
    }
    if (logfiles->curr_log != NULL && logfiles->curr_log != stdout) {
        fclose(logfiles->curr_log);
    }
    free(logfiles->logs);
    free(logfiles);
}

Logfiles * open_logs(const Options *opts) {
    struct dirent *ent;
    Logfiles *logfiles = calloc(1, sizeof(*logfiles));

    if (logfiles == NULL) {
        perror("Allocating logfiles");
        return NULL;
    }

    logfiles->logs = NULL;

    if (opts->log_dir == NULL) {
        /* Logging to stdout */
        logfiles->dir_handle = NULL;
        logfiles->curr_log = stdout;
        return logfiles;
    }

    assert(opts->nlogs > 0);

    logfiles->curr_log = NULL;
    logfiles->dir_handle = opendir(opts->log_dir);

    if (logfiles->dir_handle == NULL) {
        fprintf(stderr, "Couldn't open directory %s: %s\n",
                opts->log_dir, strerror(errno));
        goto fail;
    }

    logfiles->sz = (size_t) opts->nlogs + 1;
    logfiles->logs = calloc(logfiles->sz, sizeof(Logfile));
    if (logfiles->logs == NULL) {
        perror(NULL);
        goto fail;
    }

    logfiles->log_dir_fd = dirfd(logfiles->dir_handle);
    if (logfiles->log_dir_fd < 0) {
        fprintf(stderr, "Couldn't get descriptor for %s : %s",
                opts->log_dir, strerror(errno));
        goto fail;
    }

    /* Find existing log files */
    errno = 0;
    while ((ent = readdir(logfiles->dir_handle)) != NULL) {
        char dt[15];
        unsigned int idx;
        char suff[8];
        struct stat st;

        if (sscanf(ent->d_name, "ref_cache_%14[0-9]_%u.%7s", dt, &idx, suff)
            && (strcmp(suff, "log") == 0 || strcmp(suff, "log.gz") == 0)) {
            if (logfiles->nlogs == logfiles->sz) {
                size_t new_sz = logfiles->sz * 2;
                Logfile *new_logs = realloc(logfiles->logs, new_sz * sizeof(Logfile));
                if (new_logs == NULL) {
                    perror(NULL);
                    goto fail;
                }
                memset(&new_logs[logfiles->nlogs], 0,
                       (new_sz - logfiles->nlogs) * sizeof(Logfile));
                logfiles->logs = new_logs;
                logfiles->sz = new_sz;
            }

            int l = snprintf(logfiles->logs[logfiles->nlogs].name, LOG_NAME_LEN,
                             "%s", ent->d_name);
            if (l >= LOG_NAME_LEN)
                abort();  // should never happen
            if (fstatat(logfiles->log_dir_fd,
                        ent->d_name, &st, AT_SYMLINK_NOFOLLOW) != 0) {
                fprintf(stderr, "Warning: Couldn't stat %s/%s : %s\n",
                        opts->log_dir, ent->d_name, strerror(errno));
                continue;
            }
            if (S_ISREG(st.st_mode)) {
                logfiles->logs[logfiles->nlogs].size = st.st_size;
                logfiles->nlogs++;
            }
        }
    }

    // Sort into reverse date order
    if (logfiles->nlogs > 0) {
        qsort(logfiles->logs, logfiles->nlogs, sizeof(Logfile), log_compare);
    }

    // Call rotate_logs() to trim old logs and start a new one
    if (rotate_logs(logfiles, opts) < 0)
        goto fail;

    return logfiles;

 fail:
    close_logs(logfiles);
    return NULL;
}

static const uint8_t needs_escape[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

int write_to_log(Logfiles *logfiles, const Options *opts,
                 const char *msg, size_t len) {
    size_t written = 0;
    // Escape any non-ascii characters
    while (written < len) {
        size_t ok = written;
        while (ok < len && !needs_escape[(uint8_t) msg[ok]])
            ok++;

        if (ok > written) {
            size_t wrote = fwrite(msg, 1, ok - written, logfiles->curr_log);
            if (wrote != ok - written)
                break;
        }
        if (ok < len && needs_escape[(uint8_t) msg[ok]]) {
            if (msg[ok] == '\\') {
                if (fprintf(logfiles->curr_log, "\\\\") < 0) break;
            } else {
                if (fprintf(logfiles->curr_log, "\\x%02x", (uint8_t) msg[ok]) < 0) break;
            }
            ok++;
        }
        written = ok;
    }

    if (written < len) {
        if (opts->log_dir != NULL) {
            fprintf(stderr, "Error writing to %s/%s : %s\n",
                    opts->log_dir, logfiles->logs[0].name, strerror(errno));
        } else {
            fprintf(stderr, "Error writing to stdout: %s\n", strerror(errno));
        }
        return -1;
    }
    if (opts->log_dir != NULL) {
        logfiles->logs[0].size += (off_t) written;
        if (logfiles->logs[0].size > opts->max_log_sz) {
            if (rotate_logs(logfiles, opts) != 0)
                return -1;
        }
    }
    return 0;
}
