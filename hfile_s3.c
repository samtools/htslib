/*  hfile_s3.c -- network backend for low-level input/output streams.

    Copyright (C) 2013-2014 Genome Research Ltd.

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

#include <stdlib.h>
#include <errno.h>

#include "hfile_internal.h"

#include "htslib/ks3file.h"

typedef struct {
    hFILE base;
    kurl_t *netfp;
} hFILE_s3;

static int s3_inited = 0;

#ifdef _WIN32
static void s3_exit(void)
{
    kurl_destroy();
}
#endif

static int s3_init(void)
{
#ifdef _WIN32
    if (kurl_init() != 0) return -1;

    // In the unlikely event atexit() fails, it's better to succeed here and
    // carry on and do the I/O; then eventually when the program exits, we'll
    // merely have failed to clean up properly, as if we had aborted.
    (void) atexit(s3_exit);
#endif

    s3_inited = 1;
    return 0;
}

static ssize_t s3_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_s3 *fp = (hFILE_s3 *) fpv;
    return kurl_read(fp->netfp, buffer, nbytes);
}

static off_t s3_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_s3 *fp = (hFILE_s3 *) fpv;
    return kurl_seek(fp->netfp, offset, whence);
}

static int s3_close(hFILE *fpv)
{
    hFILE_s3 *fp = (hFILE_s3 *) fpv;
    return kurl_close(fp->netfp);
}

static const struct hFILE_backend s3_backend =
{
    s3_read, NULL, s3_seek, NULL, s3_close
};

hFILE *hopen_s3(const char *filename, const char *mode)
{
    hFILE_s3 *fp;

    // Do any networking initialisation if this is the first use.
    if (! s3_inited) { if (s3_init() < 0) return NULL; }

    fp = (hFILE_s3 *) hfile_init(sizeof (hFILE_s3), mode, 0);
    if (fp == NULL) return NULL;

    fp->netfp = kurl_open(filename, 0);
    if (fp->netfp == NULL) { hfile_destroy((hFILE *) fp); return NULL; }

    fp->base.backend = &s3_backend;
    return &fp->base;
}
