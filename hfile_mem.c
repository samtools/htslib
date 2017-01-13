/* The MIT License

   Copyright (c) 2016 Illumina Cambridge Ltd.

   Author: Peter Krusche <pkrusche@illumina.com>

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "htslib/hfile.h"
#include "htslib/hfile_mem.h"
#include "hfile_internal.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdint.h>
#include <errno.h>


static buffer_lookup_fn hfile_mem_lookup_buffer = NULL;
void hfile_mem_set_lookup_function(buffer_lookup_fn fn)
{
    hfile_mem_lookup_buffer = fn;
}


typedef struct
{
    hFILE base;
    char *filename;
    char *mode;
    size_t buffer_size;
    size_t used_size;
    off_t offset;
    uint8_t *buffer;
    int buffer_is_mine;
    int write_flag;
} hFILE_mem;


/*
 * Implementation
 */


static ssize_t mem_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    const size_t max_read = fp->used_size - fp->offset;
    const size_t to_read = max_read < nbytes ? max_read : nbytes;

    if(fp->offset >= fp->buffer_size)
    {
        return 0;
    }
    memcpy(buffer, fp->buffer + fp->offset, to_read);
    fp->offset += to_read;
    return to_read;
}

static ssize_t mem_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    const ssize_t available = fp->buffer_size - fp->offset;
    const size_t round_mask = ((ssize_t) -1) << 10;
    void *tmp = NULL;
    size_t new_buffer_size;

    if(!fp->buffer_is_mine)
    {
        // Cannot write: we don't own the buffer and can only read
        errno = EROFS;
        return -1;
    }

    if(available < nbytes)
    {
        new_buffer_size = (fp->offset + nbytes + 1023) & round_mask;
        tmp = realloc(fp->buffer, new_buffer_size) ;
        if(tmp == NULL)
        {
            return -1;
        }
        fp->buffer_size = new_buffer_size;
        fp->buffer = tmp;
    }
    fp->write_flag = 1;
    memcpy(fp->buffer + fp->offset, buffer, nbytes);
    fp->offset += nbytes;
    if(fp->offset > fp->used_size)
    {
        fp->used_size = (size_t) fp->offset;
    }
    return nbytes;
}

static off_t mem_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    if(whence == SEEK_END)
    {
        fp->offset = (off_t) fp->buffer_size + offset;
        return fp->offset;
    }
    else if(whence == SEEK_CUR)
    {
        fp->offset += offset;
        return fp->offset;
    }
    else if(whence == SEEK_SET)
    {
        fp->offset = offset;
        return fp->offset;
    }
    else return -1;
}

static int mem_close(hFILE *fpv)
{
    hFILE_mem *fp = (hFILE_mem *) fpv;
    if(fp->filename)
    {
        free(fp->filename);
    }
    if(fp->mode)
    {
        free(fp->mode);
    }
    if(fp->buffer_is_mine && fp->buffer)
    {
        free(fp->buffer);
    }
    return 0;
}

static const struct hFILE_backend mem_backend = {
        mem_read, mem_write, mem_seek, NULL, mem_close
};

hFILE *hopen_mem(const char *filename, const char *mode)
{
    hFILE_mem *fp;
    FILE *fpr;
    size_t len;

    const char *realfilename = strchr(filename, ':');
    if(!realfilename)
    {
        realfilename = filename;
    }
    else
    {
        ++realfilename;
    }
    fp = (hFILE_mem *) hfile_init(sizeof(hFILE_mem), mode, 0);
    if(!fp)
    {
        return NULL;
    }

    fp->base.backend = &mem_backend;
    fp->buffer = NULL;
    fp->buffer_size = 0;
    fp->used_size = 0;
    fp->write_flag = 0;
    fp->offset = 0;
    fp->mode = strdup(mode);
    fp->buffer_is_mine = 0;

    if(realfilename[0] == '@')
    {
        if(hfile_mem_lookup_buffer == NULL)
        {
            free(fp);
            errno = EINVAL;
            return NULL;
        }
        ++realfilename;
        fp->filename = NULL;
        if(hfile_mem_lookup_buffer(realfilename, (void**)&fp->buffer, &fp->buffer_size))
        {
            free(fp);
            errno = EINVAL;
            return NULL;
        }

        fp->used_size = fp->buffer_size;
    }
    else
    {
        fp->filename = strdup(realfilename);

        if(strchr(mode, 'r'))
        {
            fpr = fopen(realfilename, mode);
            if(!fpr)
            {
                //            fprintf(stderr, "[E::mem_file] Cannot open %s for reading.\n", filename);
                // don't write an error, this happens all the time when htslib tries to open a
                // csi file that doesn't exist
                free(fp);
                return NULL;
            }
            fseek(fpr, 0, SEEK_END);
            len = ftell(fpr);
            fseek(fpr, 0, SEEK_SET);
            fp->buffer_is_mine = 1;
            fp->buffer = malloc(len);
            if(fp->buffer == NULL)
            {
                free(fp);
                fclose(fpr);
                errno = ENOMEM;
                return NULL;
            }
            if(fread(fp->buffer, 1, len, fpr) != len)
            {
                free(fp);
                fclose(fpr);
                errno = EIO;
                return NULL;
            }
            fp->buffer_size = len;
            fp->used_size = len;
            fclose(fpr);
        }
        else
        {
            fp->buffer = malloc(1024);
            fp->buffer_size = 1024;
            fp->buffer_is_mine = 1;
        }
    }
    return &fp->base;
}

int hfile_mem_get_buffer(hFILE * file, void ** buffer, size_t * length)
{
    if(file->backend != &mem_backend)
    {
        errno = EINVAL;
        return -1;
    }
    hFILE_mem *fp = (hFILE_mem *) file;

    if(fp->buffer)
    {
        *buffer = fp->buffer;
        *length = fp->used_size;
    }
    else
    {
        return -1;
    }
    return 0;
}

int hfile_plugin_init_mem(struct hFILE_plugin *self)
{
    // mem files are declared remote so they work with a tabix index
    static const struct hFILE_scheme_handler handler =
            {hopen_mem, hfile_always_remote, "mem", 0};
    self->name = "mem";
    hfile_add_scheme_handler("mem", &handler);
    return 0;
}

