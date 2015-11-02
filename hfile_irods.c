/*  hfile_irods.c -- iRODS backend for low-level file streams.

    Copyright (C) 2013, 2015 Genome Research Ltd.

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

#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "hfile_internal.h"

#include <rcConnect.h>
#include <dataObjOpen.h>
#include <dataObjRead.h>
#include <dataObjWrite.h>
#include <dataObjLseek.h>
#include <dataObjClose.h>

typedef struct {
    hFILE base;
    int descriptor;
} hFILE_irods;

static int status_errno(int status)
{
    switch (status) {
    case SYS_NO_API_PRIV: return EACCES;
    case SYS_MALLOC_ERR: return ENOMEM;
    case SYS_OUT_OF_FILE_DESC: return ENFILE;
    case SYS_BAD_FILE_DESCRIPTOR: return EBADF;
    case CAT_NO_ROWS_FOUND: return ENOENT;
    case CATALOG_ALREADY_HAS_ITEM_BY_THAT_NAME: return EEXIST;
    default: return EIO;
    }
}

static void set_errno(int status)
{
    int err = abs(status) % 1000;
    errno = err? err : status_errno(status);
}

static struct {
    rcComm_t *conn;
    rodsEnv env;
} irods = { NULL };

static void irods_exit()
{
    if (irods.conn) { (void) rcDisconnect(irods.conn); }
    irods.conn = NULL;
}

static int irods_init()
{
    struct sigaction pipehandler;
    rErrMsg_t err;
    int ret, pipehandler_ret;

    ret = getRodsEnv(&irods.env);
    if (ret < 0) goto error;

    // Prior to iRODS 4.1, rcConnect() (even if it fails) installs its own
    // SIGPIPE handler, which just prints a message and otherwise ignores the
    // signal.  Most actual SIGPIPEs encountered will pertain to e.g. stdout
    // rather than iRODS's connection, so we save and restore the existing
    // state (by default, termination; or as already set by our caller).
    pipehandler_ret = sigaction(SIGPIPE, NULL, &pipehandler);

    irods.conn = rcConnect(irods.env.rodsHost, irods.env.rodsPort,
                           irods.env.rodsUserName, irods.env.rodsZone,
                           NO_RECONN, &err);
    if (pipehandler_ret == 0) sigaction(SIGPIPE, &pipehandler, NULL);
    if (irods.conn == NULL) { ret = err.status; goto error; }

    if (strcmp(irods.env.rodsUserName, PUBLIC_USER_NAME) != 0) {
        ret = clientLogin(irods.conn);
        if (ret != 0) goto error;
    }

    return 0;

error:
    if (irods.conn) { (void) rcDisconnect(irods.conn); }
    irods.conn = NULL;
    set_errno(ret);
    return -1;
}

static ssize_t irods_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_irods *fp = (hFILE_irods *) fpv;
    openedDataObjInp_t args;
    bytesBuf_t buf;
    int ret;

    memset(&args, 0, sizeof args);
    args.l1descInx = fp->descriptor;
    args.len = nbytes;

    buf.buf = buffer;
    buf.len = nbytes;

    ret = rcDataObjRead(irods.conn, &args, &buf);
    if (ret < 0) set_errno(ret);
    return ret;
}

static ssize_t irods_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_irods *fp = (hFILE_irods *) fpv;
    openedDataObjInp_t args;
    bytesBuf_t buf;
    int ret;

    memset(&args, 0, sizeof args);
    args.l1descInx = fp->descriptor;
    args.len = nbytes;

    buf.buf = (void *) buffer; // ...the iRODS API is not const-correct here
    buf.len = nbytes;

    ret = rcDataObjWrite(irods.conn, &args, &buf);
    if (ret < 0) set_errno(ret);
    return ret;
}

static off_t irods_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_irods *fp = (hFILE_irods *) fpv;
    openedDataObjInp_t args;
    fileLseekOut_t *out = NULL;
    int ret;

    memset(&args, 0, sizeof args);
    args.l1descInx = fp->descriptor;
    args.offset = offset;
    args.whence = whence;

    ret = rcDataObjLseek(irods.conn, &args, &out);

    if (out) { offset = out->offset; free(out); }
    else offset = -1;
    if (ret < 0) { set_errno(ret); return -1; }
    return offset;
}

static int irods_close(hFILE *fpv)
{
    hFILE_irods *fp = (hFILE_irods *) fpv;
    openedDataObjInp_t args;
    int ret;

    memset(&args, 0, sizeof args);
    args.l1descInx = fp->descriptor;

    ret = rcDataObjClose(irods.conn, &args);
    if (ret < 0) set_errno(ret);
    return ret;
}

static const struct hFILE_backend irods_backend =
{
    irods_read, irods_write, irods_seek, NULL, irods_close
};

hFILE *hopen_irods(const char *filename, const char *mode)
{
    hFILE_irods *fp;
    rodsPath_t path;
    dataObjInp_t args;
    int ret;

    // Initialise the iRODS connection if this is the first use.
    if (irods.conn == NULL) { if (irods_init() < 0) return NULL; }

    if (strncmp(filename, "irods:", 6) == 0) filename += 6;
    else { errno = EINVAL; return NULL; }

    fp = (hFILE_irods *) hfile_init(sizeof (hFILE_irods), mode, 0);
    if (fp == NULL) return NULL;

    strncpy(path.inPath, filename, MAX_NAME_LEN-1);
    path.inPath[MAX_NAME_LEN-1] = '\0';

    ret = parseRodsPath(&path, &irods.env);
    if (ret < 0) goto error;

    memset(&args, 0, sizeof args);
    strcpy(args.objPath, path.outPath);
    args.openFlags = hfile_oflags(mode);
    if (args.openFlags & O_CREAT) {
        args.createMode = 0666;
        addKeyVal(&args.condInput, DEST_RESC_NAME_KW,irods.env.rodsDefResource);
    }

    ret = rcDataObjOpen(irods.conn, &args);
    if (ret < 0) goto error;
    fp->descriptor = ret;

    fp->base.backend = &irods_backend;
    return &fp->base;

error:
    hfile_destroy((hFILE *) fp);
    set_errno(ret);
    return NULL;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_irods)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { hopen_irods, hfile_always_remote, "iRODS", 50 };

    self->name = "iRODS";
    hfile_add_scheme_handler("irods", &handler);
    self->destroy = irods_exit;
    return 0;
}
