/*  hfile_net.c -- network backend for low-level input/output streams.

    Copyright (C) 2013 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>
*/

#include <stdlib.h>
#include <errno.h>

#include "hfile_internal.h"

#include "htslib/knetfile.h"

typedef struct {
    hFILE base;
    knetFile *netfp;
} hFILE_net;

static ssize_t net_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_read(fp->netfp, buffer, nbytes);
}

static off_t net_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_seek(fp->netfp, offset, whence);
}

static int net_close(hFILE *fpv)
{
    hFILE_net *fp = (hFILE_net *) fpv;
    return knet_close(fp->netfp);
}

static const struct hFILE_backend net_backend =
{
    net_read, NULL, net_seek, NULL, net_close
};

hFILE *hopen_net(const char *filename, const char *mode)
{
    hFILE_net *fp = (hFILE_net *) malloc(sizeof (hFILE_net));
    if (fp == NULL) return NULL;

    if (hinit_buffer(&fp->base, mode, 0) < 0) goto error;

    fp->netfp = knet_open(filename, mode);
    if (fp->netfp == NULL) goto error;

    fp->base.backend = &net_backend;
    return &fp->base;

error: {
    int save = errno;
    free(fp);
    errno = save;
    }
    return NULL;
}
