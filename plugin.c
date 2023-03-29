/*  plugin.c -- low-level path parsing and plugin functions.

    Copyright (C) 2015-2016, 2020 Genome Research Ltd.

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

#include <dirent.h>
#include <dlfcn.h>

#include "hts_internal.h"
#include "htslib/kstring.h"

#ifndef PLUGINPATH
#define PLUGINPATH ""
#endif

static DIR *open_nextdir(struct hts_path_itr *itr)
{
    DIR *dir;

    while (1) {
        const char *colon = strchr(itr->pathdir, HTS_PATH_SEPARATOR_CHAR);
        if (colon == NULL) return NULL;

        itr->entry.l = 0;
        kputsn(itr->pathdir, colon - itr->pathdir, &itr->entry);
        itr->pathdir = &colon[1];
        if (itr->entry.l == 0) continue;

        dir = opendir(itr->entry.s);
        if (dir) break;

        if (hts_verbose >= 4)
            fprintf(stderr,
                    "[W::hts_path_itr] can't scan directory \"%s\": %s\n",
                    itr->entry.s, strerror(errno));
    }

    if (itr->entry.s[itr->entry.l-1] != '/') kputc('/', &itr->entry);
    itr->entry_dir_l = itr->entry.l;
    return dir;
}

void hts_path_itr_setup(struct hts_path_itr *itr, const char *path,
        const char *builtin_path, const char *prefix, size_t prefix_len,
        const char *suffix, size_t suffix_len)
{
    itr->prefix = prefix;
    itr->prefix_len = prefix_len;

    if (suffix) itr->suffix = suffix, itr->suffix_len = suffix_len;
    else itr->suffix = PLUGIN_EXT, itr->suffix_len = strlen(PLUGIN_EXT);

    itr->path.l = itr->path.m = 0; itr->path.s = NULL;
    itr->entry.l = itr->entry.m = 0; itr->entry.s = NULL;

    if (! builtin_path) builtin_path = PLUGINPATH;
    if (! path) {
        path = getenv("HTS_PATH");
        if (! path) path = "";
    }

    while (1) {
        size_t len = strcspn(path, HTS_PATH_SEPARATOR_STR);
        if (len == 0) kputs(builtin_path, &itr->path);
        else kputsn(path, len, &itr->path);
        kputc(HTS_PATH_SEPARATOR_CHAR, &itr->path);

        path += len;
        if (*path == HTS_PATH_SEPARATOR_CHAR) path++;
        else break;
    }

    // Note that ':' now terminates entries rather than separates them
    itr->pathdir = itr->path.s;
    itr->dirv = open_nextdir(itr);
}

const char *hts_path_itr_next(struct hts_path_itr *itr)
{
    while (itr->dirv) {
        struct dirent *e;
        while ((e = readdir((DIR *) itr->dirv)) != NULL) {
            size_t d_name_len = strlen(e->d_name);
            if (strncmp(e->d_name, itr->prefix, itr->prefix_len) == 0 &&
                d_name_len >= itr->suffix_len &&
                strncmp(e->d_name + d_name_len - itr->suffix_len, itr->suffix,
                        itr->suffix_len) == 0) {
                itr->entry.l = itr->entry_dir_l;
                kputs(e->d_name, &itr->entry);
                return itr->entry.s;
            }
        }

        closedir((DIR *) itr->dirv);
        itr->dirv = open_nextdir(itr);
    }

    itr->pathdir = NULL;
    free(itr->path.s); itr->path.s = NULL;
    free(itr->entry.s); itr->entry.s = NULL;
    return NULL;
}


#ifndef RTLD_NOLOAD
#define RTLD_NOLOAD 0
#endif

plugin_void_func *load_plugin(void **pluginp, const char *filename, const char *symbol)
{
    void *lib = dlopen(filename, RTLD_NOW | RTLD_LOCAL);
    if (lib == NULL) goto error;

    plugin_void_func *sym;
    *(void **) &sym = dlsym(lib, symbol);
    if (sym == NULL) {
        // Reopen the plugin with RTLD_GLOBAL and check for uniquified symbol
        void *libg = dlopen(filename, RTLD_NOLOAD | RTLD_NOW | RTLD_GLOBAL);
        if (libg == NULL) goto error;
        dlclose(lib);
        lib = libg;

        kstring_t symbolg = { 0, 0, NULL };
        kputs(symbol, &symbolg);
        kputc('_', &symbolg);
        const char *slash = strrchr(filename, '/');
        const char *basename = slash? slash+1 : filename;
        kputsn(basename, strcspn(basename, ".-+"), &symbolg);

        *(void **) &sym = dlsym(lib, symbolg.s);
        free(symbolg.s);
        if (sym == NULL) goto error;
    }

    *pluginp = lib;
    return sym;

error:
    if (hts_verbose >= 4)
        fprintf(stderr, "[W::%s] can't load plugin \"%s\": %s\n",
                __func__, filename, dlerror());
    if (lib) dlclose(lib);
    return NULL;
}

void *plugin_sym(void *plugin, const char *name, const char **errmsg)
{
    void *sym = dlsym(plugin, name);
    if (sym == NULL) *errmsg = dlerror();
    return sym;
}

plugin_void_func *plugin_func(void *plugin, const char *name, const char **errmsg)
{
    plugin_void_func *sym;
    *(void **) &sym = plugin_sym(plugin, name, errmsg);
    return sym;
}

void close_plugin(void *plugin)
{
    if (dlclose(plugin) != 0) {
        if (hts_verbose >= 4)
            fprintf(stderr, "[W::%s] dlclose() failed: %s\n",
                    __func__, dlerror());
    }
}

const char *hts_plugin_path(void) {
#ifdef ENABLE_PLUGINS
    char *path = getenv("HTS_PATH");
    if (!path) path = "";

    kstring_t ks = {0};
    while(1) {
        size_t len = strcspn(path, HTS_PATH_SEPARATOR_STR);
        if (len == 0) kputs(PLUGINPATH, &ks);
        else kputsn(path, len, &ks);
        kputc(HTS_PATH_SEPARATOR_CHAR, &ks);

        path += len;
        if (*path == HTS_PATH_SEPARATOR_CHAR) path++;
        else break;
    }

    static char s_path[1024];
    snprintf(s_path, sizeof(s_path), "%s", ks.s ? ks.s : "");
    free(ks.s);

    return s_path;
#else
    return NULL;
#endif
}
