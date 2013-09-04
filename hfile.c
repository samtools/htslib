/*  hfile.c -- buffered low-level input/output streams.

    Copyright (C) 2013 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "hfile.h"
#include "hfile_internal.h"

/* hFILE fields are used as follows:

   char *buffer;     // Pointer to the start of the I/O buffer
   char *begin;      // First not-yet-read character / unused position
   char *end;        // First unfilled/unfillable position
   size_t capacity;  // Size of the I/O buffer

   const hFILE_backend *backend;  // Methods to refill/flush I/O buffer

   int writing:1;    // Whether the buffer is being used for writing
   int at_eof:1;     // For reading, whether EOF has been seen
   int has_errno;    // Error number from the last failure on this stream

For reading, begin is the first unread character in the buffer and end is the
first unfilled position:

   -----------ABCDEFGHIJKLMNO---------------
   ^buffer    ^begin         ^end

For writing, begin is the first unused position and end is the end of the
buffer (so end == &buffer[capacity]):

   ABCDEFGHIJKLMNOPQRSTUVWXYZ---------------
   ^buffer                   ^begin         ^end

In both cases, available characters/positions are in [begin,end).  */

int hinit_buffer(hFILE *fp, const char *mode, size_t capacity)
{
    fp->capacity = capacity? capacity : 32768;
    fp->buffer = (char *) malloc(fp->capacity);
    if (fp->buffer == NULL) return -1;

    fp->writing = (strchr(mode, 'w') != NULL || strchr(mode, 'a') != NULL);

    fp->begin = fp->buffer;
    fp->end = fp->writing? &fp->buffer[fp->capacity] : fp->buffer;

    fp->at_eof = 0;
    fp->has_errno = 0;
    return 0;
}

/* Refills the read buffer from the backend (once, so may only partially
   fill the buffer), returning the number of additional characters read
   (which might be 0), or negative when an error occurred.  */
static ssize_t refill_buffer(hFILE *fp)
{
    size_t avail;
    ssize_t n;

    // Move any unread characters to the start of the buffer
    if (fp->begin > fp->buffer) {
        memmove(fp->buffer, fp->begin, fp->end - fp->begin);
        fp->end = &fp->buffer[fp->end - fp->begin];
        fp->begin = fp->buffer;
    }

    // Read into the available buffer space at fp->[end,capacity)
    avail = &fp->buffer[fp->capacity] - fp->end;
    if (fp->at_eof || avail == 0) n = 0;
    else {
        n = fp->backend->read(fp, fp->end, avail);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
    }

    fp->end += n;
    return n;
}

/* Called only from hgetc(), when our buffer is empty.  */
int hgetc2(hFILE *fp)
{
    return (refill_buffer(fp) > 0)? *(fp->begin++) : EOF;
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
    char *dest = (char *) destv;
    dest += nread, nbytes -= nread;

    // Read large requests directly into the destination buffer
    while (nbytes * 2 >= fp->capacity && !fp->at_eof) {
        ssize_t n = fp->backend->read(fp, dest, nbytes);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
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
    }

    fp->begin = fp->buffer;  // Leave the buffer empty
    return 0;
}

int hflush(hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    if (fp->backend->flush(fp) < 0) { fp->has_errno = errno; return EOF; }
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
    size_t remaining = totalbytes - ncopied;
    src += ncopied;

    ret = flush_buffer(fp);
    if (ret < 0) return ret;

    // Write large blocks out directly from the source buffer
    while (remaining * 2 >= fp->capacity) {
        ssize_t n = fp->backend->write(fp, src, remaining);
        if (n < 0) { fp->has_errno = errno; return n; }
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
    off_t pos = fp->backend->seek(fp, offset, whence);
    if (pos < 0) fp->has_errno = errno;
    return pos;
}

off_t htell(hFILE *fp)
{
    off_t pos = fp->backend->tell(fp);
    if (pos < 0) fp->has_errno = errno;
    return pos;
}

int hclose(hFILE *fp)
{
    int err = fp->has_errno;

    if (fp->writing && hflush(fp) < 0) err = fp->has_errno;
    free(fp->buffer);
    if (fp->backend->close(fp) < 0) err = errno;

    if (err) {
        errno = err;
        return EOF;
    }
    else return 0;
}


/***************************
 * File descriptor backend *
 ***************************/

#include <sys/socket.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef _WIN32
#define HAVE_CLOSESOCKET
#endif

/* For Unix, it doesn't matter whether a file descriptor is a socket.
   However Windows insists on send()/recv() and its own closesocket()
   being used when fd happens to be a socket.  */

typedef struct {
    hFILE base;
    int fd;
    int is_socket:1;
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

static off_t fd_tell(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    return lseek(fp->fd, 0, SEEK_CUR);
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
        // Ignore invalid-for-fsync(2) errors due to being, e.g., a pipe
        if (ret < 0 && errno == EINVAL) ret = 0;
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
    fd_read, fd_write, fd_seek, fd_tell, fd_flush, fd_close
};

static size_t blksize(int fd)
{
    struct stat sbuf;
    if (fstat(fd, &sbuf) != 0) return 0;
    return sbuf.st_blksize;
}

static hFILE *hopen_fd(const char *filename, const char *mode)
{
    hFILE_fd *fp = (hFILE_fd *) malloc(sizeof (hFILE_fd));
    if (fp == NULL) goto error;

    fp->fd = open(filename, hinit_oflags(mode), 0666);
    if (fp->fd < 0) goto error;

    if (hinit_buffer(&fp->base, mode, blksize(fp->fd)) < 0) goto error;

    fp->is_socket = 0;
    fp->base.backend = &fd_backend;
    return &fp->base;

error: {
    int save = errno;
    if (fp && fp->fd >= 0) close(fp->fd);
    free(fp);
    errno = save;
    }
    return NULL;
}

hFILE *hdopen(int fd, const char *mode)
{
    hFILE_fd *fp = (hFILE_fd *) malloc(sizeof (hFILE_fd));
    if (fp == NULL) return NULL;

    if (hinit_buffer(&fp->base, mode, blksize(fd)) < 0) {
        int save = errno; free(fp); errno = save;
        return NULL;
    }

    fp->fd = fd;
    fp->is_socket = (strchr(mode, 's') != NULL);
    fp->base.backend = &fd_backend;
    return &fp->base;
}

static hFILE *hopen_fd_stdinout(const char *mode)
{
    int fd = (strchr(mode, 'r') != NULL)? STDIN_FILENO : STDOUT_FILENO;
    // TODO Set binary mode (for Windows)
    return hdopen(fd, mode);
}

int hinit_oflags(const char *mode)
{
    int rdwr = 0, flags = 0;
    const char *s;
    for (s = mode; *s; s++)
        switch (*s) {
        case 'r': rdwr = O_RDONLY;  break;
        case 'w': rdwr = O_WRONLY; flags |= O_CREAT | O_TRUNC;  break;
        case 'a': rdwr = O_WRONLY; flags |= O_CREAT | O_APPEND;  break;
        case '+': rdwr = O_RDWR;  break;
        default:  break;
        }

#ifdef O_BINARY
    flags |= O_BINARY;
#endif

    return rdwr | flags;
}


/******************************
 * hopen() backend dispatcher *
 ******************************/

hFILE *hopen(const char *fname, const char *mode)
{
    if (strncmp(fname, "http://", 7) == 0 ||
        strncmp(fname, "ftp://", 6) == 0) return hopen_net(fname, mode);
    else if (strcmp(fname, "-") == 0) return hopen_fd_stdinout(mode);
    else return hopen_fd(fname, mode);
}
