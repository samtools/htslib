/*  hfile_dx.cpp -- DNA Nexus low-level file streams.

    Copyright (C) 2015-2017, 2019-2020 Genome Research Ltd.

    Author: Pierre Lindenbaum

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

#include "dxcpp/api.h"
#include "config.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "hfile_internal.h"

using namespace std;
using namespace dx;

typedef struct {
    hFILE base;
    DXFile *dxfile;
    off_t file_size;
    off_t offset;
} hFILE_dx;


static ssize_t dx_read(hFILE *fpv, void *bufferv, size_t nbytes)
{
    hFILE_dx *fp = (hFILE_dx *) fpv;
    int64_t n = 0;
    ssize_t n_reads = 0;
    size_t remain = nbytes;
    while(!fp->dxfile->eof() && remain>0) {
      fp->dxfile-> read(&((char*)bufferv)[n_reads],nbytes);
      n = fp->dxfile->gcount();
      if(n==0) break;
      fp->offset+=n;
      remain-=n;
      n_reads+=n;
      }
    return n_reads;

}


static ssize_t dx_write(hFILE *fpv, const void *bufferv, size_t nbytes)
{
    return -1;
}


static off_t dx_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_dx* fp = (hFILE_dx*)fpv;
    off_t origin, pos;


    switch (whence) {
    case SEEK_SET:
        origin = 0;
        break;
    case SEEK_CUR:
        errno = ENOSYS;
        return -1;
    case SEEK_END:
        if (fp->file_size < 0) { errno = ESPIPE; return -1; }
        origin = fp->file_size;
        break;
    default:
        errno = EINVAL;
        return -1;
    }

    // Check 0 <= origin+offset < fp->file_size carefully, avoiding overflow
    if ((offset < 0)? origin + offset < 0
                : (fp->file_size >= 0 && offset > fp->file_size - origin)) {
        errno = EINVAL;
        return -1;
    }

  
    try {
      fp->dxfile->seek(origin + offset);
      }
    catch(...) {
      fprintf(stderr,"cannot seek\n");
      return EINVAL;
      }
   fp->offset = origin + offset;
   return fp->offset;
}


static int dx_close(hFILE *fpv)
{
    hFILE_dx *fp = (hFILE_dx *) fpv;
    fp->dxfile->close();
    delete fp->dxfile;
    fp->dxfile=NULL;
    else return 0;
}

static const struct hFILE_backend dx_backend =
{
    dx_read, dx_write, dx_seek, NULL, dx_close
};


static hFILE *dx_open(const char *fname, const char *mode)
{
if(strncmp(fname,"dx:")!=0) {
  fprintf(stderr,"illegal, %s doesn't start with 'dx:'\n",fname);
  goto early_error;
  }
hFILE_dx* fp = (hFILE_dx *) hfile_init(sizeof (hFILE_dx), mode, 0);
if (fp == NULL) goto early_error;
fp->offset = 0;
fp->dxfile = new DXFile(fname);
if(!fp->dxfile.is_open()) {
  fprintf(stderr,"Cannot open dx file \"%s\".\n",fname);
  goto error;
  }
fp->file_size =  fp->dxfile->describe()["size"].get<int64_t>();
fp->base.backend = &libdx_backend;

return fp;

error:
    if(fp->dxfile!=NULL) delete fp->dxfile;
    hfile_destroy((hFILE *) fp);
    return NULL;

early_error:
    return NULL;
}



int PLUGIN_GLOBAL(hfile_plugin_init,_dx)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        {
       /* Opens a stream when dispatched by hopen(); should call hfile_init()
       to malloc a struct "derived" from hFILE and initialise it appropriately,
       including setting base.backend to its own backend vector.  */
        dx_open,
       /* Returns whether the URL denotes remote storage when dispatched by
       hisremote().  For simple cases, use one of hfile_always_*() below.  */
        hfile_always_remote,
        /* The name of the plugin or other code providing this handler.  */
        "DX DNAnexus",
        /* priority */
        5000 + 50, 
        /* Same as the open() method, used when extra arguments have been given to hopen().  */
        NULL
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_dx plugin (htslib)\t" HTS_VERSION_TEXT;
#endif

    self->name = "DX DNA Nexus";
    hfile_add_scheme_handler("dx", &handler);
    return 0;
}
