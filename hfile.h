/*  hfile.h -- buffered low-level input/output streams.

    Copyright (C) 2013 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>
*/

#ifndef HFILE_H
#define HFILE_H

#include <string.h>

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* These fields are declared here solely for the benefit of the inline functions
   below.  They may change in future releases.  User code should not use them
   directly; you should imagine that hFILE is an opaque incomplete type.  */
struct hFILE_backend;
typedef struct hFILE {
    char *buffer, *begin, *end;
    size_t capacity;
    const struct hFILE_backend *backend;
    int writing:1, at_eof:1;
    int has_errno;
} hFILE;

/*!
  @abstract  Open the named file or URL as a stream
  @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
*/
hFILE *hopen(const char *filename, const char *mode);

/*!
  @abstract  Associate a stream with an existing open file descriptor
  @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
  @notes     For socket descriptors (on Windows), mode should contain 's'.
*/
hFILE *hdopen(int fd, const char *mode);

/*!
  @abstract  Flush (for output streams) and close the stream
  @return    0 if successful, or EOF (with errno set) if an error occurred.
*/
int hclose(hFILE *fp);

/*!
  @abstract  Return the stream's error indicator
  @return    Non-zero (in fact, an errno value) if an error has occurred.
  @notes     This would be called herror() and return true/false to parallel
    ferror(3), but a networking-related herror(3) function already exists.  */
static inline int herrno(hFILE *fp)
{
    return fp->has_errno;
}

/*!
  @abstract  Reposition the read/write stream offset
  @return    The resulting offset within the stream (as per lseek(2)),
    or negative if an error occurred.
*/
off_t hseek(hFILE *fp, off_t offset, int whence);

/*!
  @abstract  Report the current stream offset
  @return    The offset within the stream, or negative if an error occurred.
*/
off_t htell(hFILE *fp);

/*!
  @abstract  Read one character from the stream
  @return    The character read, or EOF on end-of-file or error
*/
static inline int hgetc(hFILE *fp)
{
    extern int hgetc2(hFILE *);
    return (fp->end > fp->begin)? *(fp->begin++) : hgetc2(fp);
}

/*!
  @abstract  Peek at characters to be read without removing them from buffers
  @param fp      The file stream
  @param buffer  The buffer to which the peeked bytes will be written
  @param nbytes  The number of bytes to peek at; limited by the size of the
    internal buffer, which could be as small as 4K.
  @return    The number of bytes peeked, which may be less than nbytes if EOF
    is encountered; or negative, if there was an I/O error.
  @notes  The characters peeked at remain in the stream's internal buffer,
    and will be returned by later hread() etc calls.
*/
ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes);

/*!
  @abstract  Read a block of characters from the file
  @return    The number of bytes read, or negative if an error occurred.
  @notes     The full nbytes requested will be returned, except as limited
    by EOF or I/O errors.
*/
static inline ssize_t hread(hFILE *fp, void *buffer, size_t nbytes)
{
    extern ssize_t hread2(hFILE *, void *, size_t, size_t);

    size_t n = fp->end - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    fp->begin += n;
    return (n == nbytes)? (ssize_t) n : hread2(fp, buffer, nbytes, n);
}

/*!
  @abstract  Write a character to the stream
  @return    The character written, or EOF if an error occurred.
*/
static inline int hputc(int c, hFILE *fp)
{
    extern int hputc2(int, hFILE *);
    if (fp->begin < fp->end) *(fp->begin++) = c;
    else c = hputc2(c, fp);
    return c;
}

/*!
  @abstract  Write a string to the stream
  @return    0 if successful, or EOF if an error occurred.
*/
static inline int hputs(const char *text, hFILE *fp)
{
    extern int hputs2(const char *, size_t, size_t, hFILE *);

    size_t nbytes = strlen(text), n = fp->end - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, text, n);
    fp->begin += n;
    return (n == nbytes)? 0 : hputs2(text, nbytes, n, fp);
}

/*!
  @abstract  Write a block of characters to the file
  @return    Either nbytes, or negative if an error occurred.
  @notes     In the absence of I/O errors, the full nbytes will be written.
*/
static inline ssize_t hwrite(hFILE *fp, const void *buffer, size_t nbytes)
{
    extern ssize_t hwrite2(hFILE *, const void *, size_t, size_t);

    size_t n = fp->end - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, buffer, n);
    fp->begin += n;
    return (n==nbytes)? (ssize_t) n : hwrite2(fp, buffer, nbytes, n);
}

/*!
  @abstract  For writing streams, flush buffered output to the underlying stream
  @return    0 if successful, or EOF if an error occurred.
*/
int hflush(hFILE *fp);

#ifdef __cplusplus
}
#endif

#endif
