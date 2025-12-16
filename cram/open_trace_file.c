/*
Author: James Bonfield

Copyright (c) 2000-2001 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
Copyright (c) 2008, 2009, 2013, 2014-2015, 2018-2020, 2025 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "os.h"
#ifndef PATH_MAX
#  define PATH_MAX 1024
#endif

#include "open_trace_file.h"
#include "misc.h"
#include "../htslib/hfile.h"
#include "../htslib/hts_log.h"
#include "../htslib/hts.h"

/*
 * Returns whether the path refers to a regular file.
 */
static int is_file(char *fn) {
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISREG(buf.st_mode);
}

/*
 * Tokenises the search path splitting on colons (unix) or semicolons
 * (windows).
 * We also  explicitly add a "./" to the end of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
char *tokenise_search_path(const char *searchpath) {
    char *newsearch;
    unsigned int i, j;
    size_t len;
    char path_sep = HTS_PATH_SEPARATOR_CHAR;

    if (!searchpath)
        searchpath="";

    newsearch = (char *)malloc((len = strlen(searchpath))+5);
    if (!newsearch)
        return NULL;

    for (i = 0, j = 0; i < len; i++) {
        /* "::" => ":". Used for escaping colons in http://foo */
        if (i < len-1 && searchpath[i] == ':' && searchpath[i+1] == ':') {
            newsearch[j++] = ':';
            i++;
            continue;
        }

        /* Handle http:// and ftp:// too without :: */
        if (path_sep == ':') {
            if ((i == 0 || (i > 0 && searchpath[i-1] == ':')) &&
                (!strncmp(&searchpath[i], "http:",     5) ||
                 !strncmp(&searchpath[i], "https:",    6) ||
                 !strncmp(&searchpath[i], "ftp:",      4) ||
                 !strncmp(&searchpath[i], "|http:",    6) ||
                 !strncmp(&searchpath[i], "|https:",   7) ||
                 !strncmp(&searchpath[i], "|ftp:",     5) ||
                 !strncmp(&searchpath[i], "URL=http:", 9) ||
                 !strncmp(&searchpath[i], "URL=https:",10)||
                 !strncmp(&searchpath[i], "URL=ftp:",  8))) {
                do {
                    newsearch[j++] = searchpath[i];
                } while (i<len && searchpath[i++] != ':');
                if (searchpath[i] == ':')
                    i++;
                if (searchpath[i]=='/')
                    newsearch[j++] = searchpath[i++];
                if (searchpath[i]=='/')
                    newsearch[j++] = searchpath[i++];
                // Look for host:port
                do {
                    newsearch[j++] = searchpath[i++];
                } while (i<len && searchpath[i] != ':' && searchpath[i] != '/');
                newsearch[j++] = searchpath[i++];
                if (searchpath[i] == ':')
                    i++;
            }
        }

        if (searchpath[i] == path_sep) {
            /* Skip blank path components */
            if (j && newsearch[j-1] != 0)
                newsearch[j++] = 0;
        } else {
            newsearch[j++] = searchpath[i];
        }
    }

    if (j)
        newsearch[j++] = 0;
    newsearch[j++] = '.';
    newsearch[j++] = '/';
    newsearch[j++] = 0;
    newsearch[j++] = 0;

    return newsearch;
}

static char *expand_path(const char *file, const char *dirname,
                         int max_s_digits);

mFILE *find_file_url(const char *file, char *url) {
    char *path = NULL, buf[8192];
    mFILE *mf = NULL;
    ssize_t len;
    hFILE *hf = NULL;

    /* Expand %s for the trace name.  Only one digit is allowed between
       The % and s to avoid ambiguity with percent-encoded URLs */

    path = expand_path(file, url, 1);
    if (!path)
        return NULL;

    if (!(hf = hopen(path, "r"))) {
        if (errno != ENOENT)
            hts_log_warning("Failed to open reference \"%s\": %s", path, strerror(errno));
        goto fail;
    }

    if (NULL == (mf = mfcreate(NULL, 0)))
        goto fail;
    while ((len = hread(hf, buf, sizeof(buf))) > 0) {
        if (mfwrite(buf, len, 1, mf) <= 0) {
            hclose_abruptly(hf);
            goto fail;
        }
    }
    if (hclose(hf) < 0 || len < 0) {
        hts_log_warning("Failed to read reference \"%s\": %s", path, strerror(errno));
        goto fail;
    }

    free(path);
    mrewind(mf);
    return mf;

 fail:
    mfdestroy(mf);
    free(path);
    return NULL;
}

/*
 * Takes a dirname possibly including % rules and appends the filename
 * to it.
 *
 * Returns expanded pathname or NULL for malloc failure.
 */
static char *expand_path(const char *file, const char *dirname,
                         int max_s_digits) {
    size_t len = strlen(dirname);
    size_t lenf = strlen(file);
    const char *end_dirname = dirname + len, *cp;
    char *path;

    path = malloc(len+lenf+2); // worst expansion DIR/FILE
    if (!path) {
        hts_log_error("Out of memory");
        return NULL;
    }

    // Remove trailing '/'s, unless the path matches "/"
    while (len > 1 && dirname[len-1] == '/') {
        len--;
        end_dirname--;
    }

    /* Special case for "./" or absolute filenames */
    if (*file == '/' || (len==1 && *dirname == '.')) {
        memcpy(path, file, lenf + 1);
    } else {
        /* Handle %[0-9]*s expansions, if required */
        char *path_end = path;

        while ((cp = strchr(dirname, '%'))) {
            char *endp;
            // Get optional length
            long l = strtol(cp+1, &endp, 10);
            if (*endp != 's' || l < 0 || endp - cp - 1 > max_s_digits) {
                // Not %[0-9]s. Copy over directly, taking care of edge cases
                // like the string ending with '%' or '%[0-9]*'.
                const char *e = MIN(endp+1, end_dirname);
                memcpy(path_end, dirname, e - dirname);
                path_end += e - dirname;
                dirname = e;
                continue;
            }

            // Copy part up to '%'
            memcpy(path_end, dirname, cp-dirname);
            path_end += cp-dirname;

            // Insert segment from file
            size_t to_copy = l > 0 ? MIN(lenf, l) : lenf;
            memcpy(path_end, file, to_copy);
            path_end += to_copy;
            file     += to_copy;
            lenf     -= to_copy;

            // Skip to part of dirname after the 's'
            dirname = endp+1;
        }

        // Add anything left in dirname
        if (dirname < end_dirname) {
            memcpy(path_end, dirname, end_dirname - dirname);
            path_end += end_dirname - dirname;
        }

        if (*file) {
            // Add remainder of file
            if (path_end > path && *(path_end - 1) != '/')
                *path_end++ = '/';
            memcpy(path_end, file, lenf);
            path_end += lenf;
        }
        // Terminate string
        *path_end = '\0';
    }

    // fprintf(stderr, "*PATH=\"%s\"\n", path);
    return path;
}

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns mFILE pointer if found
 *         NULL if not
 */
static mFILE *find_file_dir(const char *file, char *dirname) {
    char *path;
    mFILE *mf = NULL;

    path = expand_path(file, dirname, INT_MAX);
    if (!path)
        return NULL;

    if (is_file(path))
        mf = mfopen(path, "rbm");

    free(path);
    return mf;
}

/*
 * ------------------------------------------------------------------------
 * Public functions below.
 */

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in 'path' (which is a colon separated list).
 * If 'path' is NULL it uses the RAWDATA environment variable instead.
 *
 * If non-NULL *local is filled out to 1 for a local file and 0 for a remote
 * URL.
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_path_mfile(const char *file, char *path, char *relative_to,
                       int *local) {
    char *newsearch;
    char *ele;
    mFILE *fp;

    if (local)
        *local = 1;

    /* Use path first */
    if (!path)
        path = getenv("RAWDATA");
    if (NULL == (newsearch = tokenise_search_path(path)))
        return NULL;

    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
        char *ele2;

        /*
         * '|' prefixing a path component indicates that we do not
         * wish to perform the compression extension searching in that
         * location.
         *
         * NB: this has been removed from the htslib implementation.
         */
        if (*ele == '|') {
            ele2 = ele+1;
        } else {
            ele2 = ele;
        }

        if (0 == strncmp(ele2, "URL=", 4)) {
            if ((fp = find_file_url(file, ele2+4))) {
                if (local)
                    *local = strncmp(ele2+4, "file:", 5) == 0 ? 1 : 0;
                free(newsearch);
                return fp;
            }
        } else if (hisremote(ele2)) {
            if ((fp = find_file_url(file, ele2))) {
                free(newsearch);
                if (local)
                    *local = 0;
                return fp;
            }
        } else if ((fp = find_file_dir(file, ele2))) {
            free(newsearch);
            return fp;
        }
    }

    free(newsearch);

    /* Look in the same location as the incoming 'relative_to' filename */
    if (relative_to) {
        char *cp;
        char relative_path[PATH_MAX+1];
        strcpy(relative_path, relative_to);
        if ((cp = strrchr(relative_path, '/')))
            *cp = 0;
        if ((fp = find_file_dir(file, relative_path)))
            return fp;
    }

    return NULL;
}


/*
 * As per open_path_mfile, but searching only for local filenames.
 * This is useful as we may avoid doing a full mfopen and loading
 * the entire file into memory.
 *
 * Returns the expanded pathname if found.
 *         NULL if not
 */
char *find_path(const char *file, const char *path) {
    char *newsearch;
    char *ele;
    char *outpath = NULL;

    /* Use path first */
    if (!path)
        path = getenv("RAWDATA");
    if (NULL == (newsearch = tokenise_search_path(path)))
        return NULL;

    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
        char *ele2 = (*ele == '|') ? ele+1 : ele;

        if (!strncmp(ele2, "URL=", 4) ||
            !strncmp(ele2, "http:", 5) ||
            !strncmp(ele2, "https:", 6) ||
            !strncmp(ele2, "ftp:", 4)) {
            continue;
        } else {
            outpath = expand_path(file, ele2, INT_MAX);
            if (is_file(outpath)) {
                free(newsearch);
                return outpath;
            } else {
                free(outpath);
            }
        }
    }

    free(newsearch);

    return NULL;
}
