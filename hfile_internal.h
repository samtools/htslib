/*  hfile_internal.h -- internal parts of low-level input/output streams.

    Copyright (C) 2013 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>
*/

#ifndef HFILE_INTERNAL_H
#define HFILE_INTERNAL_H

#include "hfile.h"

struct hFILE_backend {
    /* As per read(2), returning the number of bytes read (possibly 0) or
       negative (and setting errno) on errors.  Front-end code will call this
       repeatedly if necessary to attempt to get the desired byte count.  */
    ssize_t (*read)(hFILE *fp, void *buffer, size_t nbytes) HTS_RESULT_USED;

    /* As per write(2), returning the number of bytes written or negative (and
       setting errno) on errors.  Front-end code will call this repeatedly if
       necessary until the desired block is written or an error occurs.  */
    ssize_t (*write)(hFILE *fp, const void *buffer, size_t nbytes)
        HTS_RESULT_USED;

    /* As per lseek(2), returning the resulting offset within the stream or
       negative (and setting errno) on errors.  */
    off_t (*seek)(hFILE *fp, off_t offset, int whence) HTS_RESULT_USED;

    /* Performs low-level flushing, if any, e.g., fsync(2); for writing streams
       only.  Returns 0 for success or negative (and sets errno) on errors. */
    int (*flush)(hFILE *fp) HTS_RESULT_USED;

    /* Closes the underlying stream (for output streams, the buffer will
       already have been flushed), returning 0 for success or negative (and
       setting errno) on errors, as per close(2).  */
    int (*close)(hFILE *fp) HTS_RESULT_USED;
};

/* These are called from the hopen() dispatcher, and should call hfile_init()
   to malloc a struct "derived" from hFILE and initialise it appropriately,
   including setting base.backend to their own backend vector.  */
hFILE *hopen_net(const char *filename, const char *mode);

/* May be called by hopen_*() functions to decode a fopen()-style mode into
   open(2)-style flags.  */
int hfile_oflags(const char *mode);

/* Must be called by hopen_*() functions to allocate the hFILE struct and set
   up its base.  Capacity is a suggested buffer size (e.g., via fstat(2))
   or 0 for a default-sized buffer.  */
hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity);

/* May be called by hopen_*() functions to undo the effects of hfile_init()
   in the event opening the stream subsequently fails.  (This is safe to use
   even if fp is NULL.  This takes care to preserve errno.)  */
void hfile_destroy(hFILE *fp);

#endif
