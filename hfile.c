/*  hfile.c -- buffered low-level input/output streams.

    Copyright (C) 2013-2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <pthread.h>

#include "htslib/hfile.h"
#include "hfile_internal.h"

#ifndef ENOTSUP
#define ENOTSUP EINVAL
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif
#ifndef EPROTONOSUPPORT
#define EPROTONOSUPPORT ENOSYS
#endif

/* hFILE fields are used as follows:

   char *buffer;     // Pointer to the start of the I/O buffer
   char *begin;      // First not-yet-read character / unused position
   char *end;        // First unfilled/unfillable position
   char *limit;      // Pointer to the first position past the buffer

   const hFILE_backend *backend;  // Methods to refill/flush I/O buffer

   off_t offset;     // Offset within the stream of buffer position 0
   unsigned at_eof:1;// For reading, whether EOF has been seen
   int has_errno;    // Error number from the last failure on this stream

For reading, begin is the first unread character in the buffer and end is the
first unfilled position:

   -----------ABCDEFGHIJKLMNO---------------
   ^buffer    ^begin         ^end           ^limit

For writing, begin is the first unused position and end is unused so remains
equal to buffer:

   ABCDEFGHIJKLMNOPQRSTUVWXYZ---------------
   ^buffer                   ^begin         ^limit
   ^end

Thus if begin > end then there is a non-empty write buffer, if begin < end
then there is a non-empty read buffer, and if begin == end then both buffers
are empty.  In all cases, the stream's file position indicator corresponds
to the position pointed to by begin.  */

hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity)
{
    hFILE *fp = (hFILE *) malloc(struct_size);
    if (fp == NULL) goto error;

    if (capacity == 0) capacity = 32768;
    // FIXME For now, clamp input buffer sizes so mpileup doesn't eat memory
    if (strchr(mode, 'r') && capacity > 32768) capacity = 32768;

    fp->buffer = (char *) malloc(capacity);
    if (fp->buffer == NULL) goto error;

    fp->begin = fp->end = fp->buffer;
    fp->limit = &fp->buffer[capacity];

    fp->offset = 0;
    fp->at_eof = 0;
    fp->has_errno = 0;
    return fp;

error:
    hfile_destroy(fp);
    return NULL;
}

void hfile_destroy(hFILE *fp)
{
    int save = errno;
    if (fp) free(fp->buffer);
    free(fp);
    errno = save;
}

static inline int writebuffer_is_nonempty(hFILE *fp)
{
    return fp->begin > fp->end;
}

/* Refills the read buffer from the backend (once, so may only partially
   fill the buffer), returning the number of additional characters read
   (which might be 0), or negative when an error occurred.  */
static ssize_t refill_buffer(hFILE *fp)
{
    ssize_t n;

    // Move any unread characters to the start of the buffer
    if (fp->begin > fp->buffer) {
        fp->offset += fp->begin - fp->buffer;
        memmove(fp->buffer, fp->begin, fp->end - fp->begin);
        fp->end = &fp->buffer[fp->end - fp->begin];
        fp->begin = fp->buffer;
    }

    // Read into the available buffer space at fp->[end,limit)
    if (fp->at_eof || fp->end == fp->limit) n = 0;
    else {
        n = fp->backend->read(fp, fp->end, fp->limit - fp->end);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
    }

    fp->end += n;
    return n;
}

/* Called only from hgetc(), when our buffer is empty.  */
int hgetc2(hFILE *fp)
{
    return (refill_buffer(fp) > 0)? (unsigned char) *(fp->begin++) : EOF;
}

ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)
{
    size_t n = fp->end - fp->begin;
    while (n < nbytes) {
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;
        else if (ret == 0) break;
        else n += ret;
    }

    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    return n;
}

/* Called only from hread(); when called, our buffer is empty and nread bytes
   have already been placed in the destination buffer.  */
ssize_t hread2(hFILE *fp, void *destv, size_t nbytes, size_t nread)
{
    const size_t capacity = fp->limit - fp->buffer;
    char *dest = (char *) destv;
    dest += nread, nbytes -= nread;

    // Read large requests directly into the destination buffer
    while (nbytes * 2 >= capacity && !fp->at_eof) {
        ssize_t n = fp->backend->read(fp, dest, nbytes);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
        fp->offset += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    while (nbytes > 0 && !fp->at_eof) {
        size_t n;
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;

        n = fp->end - fp->begin;
        if (n > nbytes) n = nbytes;
        memcpy(dest, fp->begin, n);
        fp->begin += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    return nread;
}

/* Flushes the write buffer, fp->[buffer,begin), out through the backend
   returning 0 on success or negative if an error occurred.  */
static ssize_t flush_buffer(hFILE *fp)
{
    const char *buffer = fp->buffer;
    while (buffer < fp->begin) {
        ssize_t n = fp->backend->write(fp, buffer, fp->begin - buffer);
        if (n < 0) { fp->has_errno = errno; return n; }
        buffer += n;
        fp->offset += n;
    }

    fp->begin = fp->buffer;  // Leave the buffer empty
    return 0;
}

int hflush(hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    if (fp->backend->flush) {
        if (fp->backend->flush(fp) < 0) { fp->has_errno = errno; return EOF; }
    }
    return 0;
}

/* Called only from hputc(), when our buffer is already full.  */
int hputc2(int c, hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    *(fp->begin++) = c;
    return c;
}

/* Called only from hwrite() and hputs2(); when called, our buffer is full and
   ncopied bytes from the source have already been copied to our buffer.  */
ssize_t hwrite2(hFILE *fp, const void *srcv, size_t totalbytes, size_t ncopied)
{
    const char *src = (const char *) srcv;
    ssize_t ret;
    const size_t capacity = fp->limit - fp->buffer;
    size_t remaining = totalbytes - ncopied;
    src += ncopied;

    ret = flush_buffer(fp);
    if (ret < 0) return ret;

    // Write large blocks out directly from the source buffer
    while (remaining * 2 >= capacity) {
        ssize_t n = fp->backend->write(fp, src, remaining);
        if (n < 0) { fp->has_errno = errno; return n; }
        fp->offset += n;
        src += n, remaining -= n;
    }

    // Just buffer any remaining characters
    memcpy(fp->begin, src, remaining);
    fp->begin += remaining;

    return totalbytes;
}

/* Called only from hputs(), when our buffer is already full.  */
int hputs2(const char *text, size_t totalbytes, size_t ncopied, hFILE *fp)
{
    return (hwrite2(fp, text, totalbytes, ncopied) >= 0)? 0 : EOF;
}

off_t hseek(hFILE *fp, off_t offset, int whence)
{
    off_t curpos, pos;

    if (writebuffer_is_nonempty(fp)) {
        int ret = flush_buffer(fp);
        if (ret < 0) return ret;
    }

    curpos = htell(fp);

    // Relative offsets are given relative to the hFILE's stream position,
    // which may differ from the backend's physical position due to buffering
    // read-ahead.  Correct for this by converting to an absolute position.
    if (whence == SEEK_CUR) {
        if (curpos + offset < 0) {
            // Either a negative offset resulted in a position before the
            // start of the file, or we overflowed when given a positive offset
            fp->has_errno = errno = (offset < 0)? EINVAL : EOVERFLOW;
            return -1;
        }

        whence = SEEK_SET;
        offset = curpos + offset;
    }

    // TODO Avoid seeking if the desired position is within our read buffer

    pos = fp->backend->seek(fp, offset, whence);
    if (pos < 0) { fp->has_errno = errno; return pos; }

    // Seeking succeeded, so discard any non-empty read buffer
    fp->begin = fp->end = fp->buffer;
    fp->at_eof = 0;

    fp->offset = pos;
    return pos;
}

int hclose(hFILE *fp)
{
    int err = fp->has_errno;

    if (writebuffer_is_nonempty(fp) && hflush(fp) < 0) err = fp->has_errno;
    if (fp->backend->close(fp) < 0) err = errno;
    hfile_destroy(fp);

    if (err) {
        errno = err;
        return EOF;
    }
    else return 0;
}

void hclose_abruptly(hFILE *fp)
{
    int save = errno;
    if (fp->backend->close(fp) < 0) { /* Ignore subsequent errors */ }
    hfile_destroy(fp);
    errno = save;
}


/***************************
 * File descriptor backend *
 ***************************/

#ifndef _WIN32
#include <sys/socket.h>
#include <sys/stat.h>
#define HAVE_STRUCT_STAT_ST_BLKSIZE
#else
#include <winsock2.h>
#define HAVE_CLOSESOCKET
#define HAVE_SETMODE
#endif
#include <fcntl.h>
#include <unistd.h>

/* For Unix, it doesn't matter whether a file descriptor is a socket.
   However Windows insists on send()/recv() and its own closesocket()
   being used when fd happens to be a socket.  */

typedef struct {
    hFILE base;
    int fd;
    unsigned is_socket:1;
} hFILE_fd;

static ssize_t fd_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket? recv(fp->fd, buffer, nbytes, 0)
                         : read(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static ssize_t fd_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket?  send(fp->fd, buffer, nbytes, 0)
                         : write(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static off_t fd_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    return lseek(fp->fd, offset, whence);
}

static int fd_flush(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;
    do {
#ifdef HAVE_FDATASYNC
        ret = fdatasync(fp->fd);
#else
        ret = fsync(fp->fd);
#endif
        // Ignore invalid-for-fsync(2) errors due to being, e.g., a pipe,
        // and operation-not-supported errors (Mac OS X)
        if (ret < 0 && (errno == EINVAL || errno == ENOTSUP)) ret = 0;
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static int fd_close(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;
    do {
#ifdef HAVE_CLOSESOCKET
        ret = fp->is_socket? closesocket(fp->fd) : close(fp->fd);
#else
        ret = close(fp->fd);
#endif
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static const struct hFILE_backend fd_backend =
{
    fd_read, fd_write, fd_seek, fd_flush, fd_close
};

static size_t blksize(int fd)
{
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
    struct stat sbuf;
    if (fstat(fd, &sbuf) != 0) return 0;
    return sbuf.st_blksize;
#else
    return 0;
#endif
}

static hFILE *hopen_fd(const char *filename, const char *mode)
{
    hFILE_fd *fp = NULL;
    int fd = open(filename, hfile_oflags(mode), 0666);
    if (fd < 0) goto error;

    fp = (hFILE_fd *) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) goto error;

    fp->fd = fd;
    fp->is_socket = 0;
    fp->base.backend = &fd_backend;
    return &fp->base;

error:
    if (fd >= 0) { int save = errno; (void) close(fd); errno = save; }
    hfile_destroy((hFILE *) fp);
    return NULL;
}

hFILE *hdopen(int fd, const char *mode)
{
    hFILE_fd *fp = (hFILE_fd*) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) return NULL;

    fp->fd = fd;
    fp->is_socket = (strchr(mode, 's') != NULL);
    fp->base.backend = &fd_backend;
    return &fp->base;
}

static hFILE *hopen_fd_fileuri(const char *url, const char *mode)
{
    if (strncmp(url, "file://localhost/", 17) == 0) url += 16;
    else if (strncmp(url, "file:///", 8) == 0) url += 7;
    else { errno = EPROTONOSUPPORT; return NULL; }

    return hopen_fd(url, mode);
}

static hFILE *hopen_fd_stdinout(const char *mode)
{
    int fd = (strchr(mode, 'r') != NULL)? STDIN_FILENO : STDOUT_FILENO;
#if defined HAVE_SETMODE && defined O_BINARY
    if (setmode(fd, O_BINARY) < 0) return NULL;
#endif
    return hdopen(fd, mode);
}

int hfile_oflags(const char *mode)
{
    int rdwr = 0, flags = 0;
    const char *s;
    for (s = mode; *s; s++)
        switch (*s) {
        case 'r': rdwr = O_RDONLY;  break;
        case 'w': rdwr = O_WRONLY; flags |= O_CREAT | O_TRUNC;  break;
        case 'a': rdwr = O_WRONLY; flags |= O_CREAT | O_APPEND;  break;
        case '+': rdwr = O_RDWR;  break;
#ifdef O_CLOEXEC
        case 'e': flags |= O_CLOEXEC;  break;
#endif
#ifdef O_EXCL
        case 'x': flags |= O_EXCL;  break;
#endif
        default:  break;
        }

#ifdef O_BINARY
    flags |= O_BINARY;
#endif

    return rdwr | flags;
}


/*********************
 * In-memory backend *
 *********************/

typedef struct {
    hFILE base;
    const char *buffer;
    size_t length, pos;
} hFILE_mem;

static ssize_t mem_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    size_t avail = fp->length - fp->pos;
    if (nbytes > avail) nbytes = avail;
    memcpy(buffer, fp->buffer + fp->pos, nbytes);
    fp->pos += nbytes;
    return nbytes;
}

static off_t mem_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    size_t absoffset = (offset >= 0)? offset : -offset;
    size_t origin;

    switch (whence) {
    case SEEK_SET: origin = 0; break;
    case SEEK_CUR: origin = fp->pos; break;
    case SEEK_END: origin = fp->length; break;
    default: errno = EINVAL; return -1;
    }

    if ((offset  < 0 && absoffset > origin) ||
        (offset >= 0 && absoffset > fp->length - origin)) {
        errno = EINVAL;
        return -1;
    }

    fp->pos = origin + offset;
    return fp->pos;
}

static int mem_close(hFILE *fpv)
{
    return 0;
}

static const struct hFILE_backend mem_backend =
{
    mem_read, NULL, mem_seek, NULL, mem_close
};

static hFILE *hopen_mem(const char *data, const char *mode)
{
    if (strncmp(data, "data:", 5) == 0) data += 5;

    // TODO Implement write modes, which will require memory allocation
    if (strchr(mode, 'r') == NULL) { errno = EINVAL; return NULL; }

    hFILE_mem *fp = (hFILE_mem *) hfile_init(sizeof (hFILE_mem), mode, 0);
    if (fp == NULL) return NULL;

    fp->buffer = data;
    fp->length = strlen(data);
    fp->pos = 0;
    fp->base.backend = &mem_backend;
    return &fp->base;
}


/*****************************************
 * Plugin and hopen() backend dispatcher *
 *****************************************/

#include "hts_internal.h"
#include "htslib/khash.h"

KHASH_MAP_INIT_STR(scheme_string, const struct hFILE_scheme_handler *);
static khash_t(scheme_string) *schemes = NULL;

struct hFILE_plugin_list {
    struct hFILE_plugin plugin;
    struct hFILE_plugin_list *next;
};

static struct hFILE_plugin_list *plugins = NULL;
static pthread_mutex_t plugins_lock = PTHREAD_MUTEX_INITIALIZER;

static void hfile_exit()
{
    pthread_mutex_lock(&plugins_lock);

    kh_destroy(scheme_string, schemes);

    while (plugins != NULL) {
        struct hFILE_plugin_list *p = plugins;
        if (p->plugin.destroy) p->plugin.destroy();
#ifdef ENABLE_PLUGINS
        if (p->plugin.obj) close_plugin(p->plugin.obj);
#endif
        plugins = p->next;
        free(p);
    }

    pthread_mutex_unlock(&plugins_lock);
    pthread_mutex_destroy(&plugins_lock);
}

void hfile_add_scheme_handler(const char *scheme,
                              const struct hFILE_scheme_handler *handler)
{
    int absent;
    khint_t k = kh_put(scheme_string, schemes, scheme, &absent);
    if (absent || handler->priority > kh_value(schemes, k)->priority) {
        kh_value(schemes, k) = handler;
    }
}

static int init_add_plugin(void *obj, int (*init)(struct hFILE_plugin *),
                           const char *pluginname)
{
    struct hFILE_plugin_list *p = malloc (sizeof (struct hFILE_plugin_list));
    if (p == NULL) abort();

    p->plugin.api_version = 1;
    p->plugin.obj = obj;
    p->plugin.name = NULL;
    p->plugin.destroy = NULL;

    int ret = (*init)(&p->plugin);

    if (ret != 0) {
        if (hts_verbose >= 4)
            fprintf(stderr, "[W::load_hfile_plugins] "
                    "initialisation failed for plugin \"%s\": %d\n",
                    pluginname, ret);
        free(p);
        return ret;
    }

    if (hts_verbose >= 5)
        fprintf(stderr, "[M::load_hfile_plugins] loaded \"%s\"\n", pluginname);

    p->next = plugins, plugins = p;
    return 0;
}

static void load_hfile_plugins()
{
    static const struct hFILE_scheme_handler
        data = { hopen_mem, hfile_always_local, "built-in", 80 },
        file = { hopen_fd_fileuri, hfile_always_local, "built-in", 80 };

    schemes = kh_init(scheme_string);
    if (schemes == NULL) abort();

    hfile_add_scheme_handler("data", &data);
    hfile_add_scheme_handler("file", &file);
    init_add_plugin(NULL, hfile_plugin_init_net, "knetfile");

#ifdef ENABLE_PLUGINS
    struct hts_path_itr path;
    const char *pluginname;
    hts_path_itr_setup(&path, NULL, NULL, "hfile_", 6, NULL, 0);
    while ((pluginname = hts_path_itr_next(&path)) != NULL) {
        void *obj;
        int (*init)(struct hFILE_plugin *) = (int (*)(struct hFILE_plugin *))
            load_plugin(&obj, pluginname, "hfile_plugin_init");

        if (init) {
            if (init_add_plugin(obj, init, pluginname) != 0)
                close_plugin(obj);
        }
    }
#else

#ifdef HAVE_IRODS
    init_add_plugin(NULL, hfile_plugin_init_irods, "iRODS");
#endif
#ifdef HAVE_LIBCURL
    init_add_plugin(NULL, hfile_plugin_init_libcurl, "libcurl");
#endif

#endif

    // In the unlikely event atexit() fails, it's better to succeed here and
    // carry on; then eventually when the program exits, we'll merely close
    // down the plugins uncleanly, as if we had aborted.
    (void) atexit(hfile_exit);
}

/* A filename like "foo:bar" in which we don't recognise the scheme is
   either an ordinary file or an indication of a missing or broken plugin.
   Try to open it as an ordinary file; but if there's no such file, set
   errno distinctively to make the plugin issue apparent.  */
static hFILE *hopen_unknown_scheme(const char *fname, const char *mode)
{
    hFILE *fp = hopen_fd(fname, mode);
    if (fp == NULL && errno == ENOENT) errno = EPROTONOSUPPORT;
    return fp;
}

/* Returns the appropriate handler, or NULL if the string isn't an URL.  */
static const struct hFILE_scheme_handler *find_scheme_handler(const char *s)
{
    static const struct hFILE_scheme_handler unknown_scheme =
        { hopen_unknown_scheme, hfile_always_local, "built-in", 0 };

    char scheme[12];
    int i;

    for (i = 0; i < sizeof scheme; i++)
        if (isalnum_c(s[i]) || s[i] == '+' || s[i] == '-' || s[i] == '.')
            scheme[i] = tolower_c(s[i]);
        else if (s[i] == ':') break;
        else return NULL;

    if (i == 0 || i >= sizeof scheme) return NULL;
    scheme[i] = '\0';

    pthread_mutex_lock(&plugins_lock);
    if (! schemes) load_hfile_plugins();
    pthread_mutex_unlock(&plugins_lock);

    khint_t k = kh_get(scheme_string, schemes, scheme);
    return (k != kh_end(schemes))? kh_value(schemes, k) : &unknown_scheme;
}

hFILE *hopen(const char *fname, const char *mode)
{
    const struct hFILE_scheme_handler *handler = find_scheme_handler(fname);
    if (handler) return handler->open(fname, mode);
    else if (strcmp(fname, "-") == 0) return hopen_fd_stdinout(mode);
    else return hopen_fd(fname, mode);
}

int hfile_always_local (const char *fname) { return 0; }
int hfile_always_remote(const char *fname) { return 1; }

int hisremote(const char *fname)
{
    const struct hFILE_scheme_handler *handler = find_scheme_handler(fname);
    return handler? handler->isremote(fname) : 0;
}
