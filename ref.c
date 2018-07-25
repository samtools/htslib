/*
Copyright (c) 2008-2009, 2012-2018 Genome Research Ltd.
Authors: James Bonfield <jkb@sanger.ac.uk>, Thomas Hickman <th10@sanger.ac.uk>

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

#include "htslib/kstring.h"

#include "cram/cram.h"
#include "cram/cram_io.h"
#include "cram/os.h"
#include "cram/misc.h"
#include "htslib/ref.h"

#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>

#include <config.h>

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX // was 1024 in open_trace_file
#endif


/*
 * Return the cache directory to use, based on the first of these
 * environment variables to be set to a non-empty value.
 */
static const char *get_cache_basedir(const char **extra)
{
    char *base;

    *extra = "";

    base = getenv("XDG_CACHE_HOME");
    if (base && *base)
        return base;

    base = getenv("HOME");
    if (base && *base)
    {
        *extra = "/.cache";
        return base;
    }

    base = getenv("TMPDIR");
    if (base && *base)
        return base;

    base = getenv("TEMP");
    if (base && *base)
        return base;

    return "/tmp";
}

/*
 * Converts a directory and a filename into an expanded path, replacing %s
 * in directory with the filename and %[0-9]+s with portions of the filename
 * Any remaining parts of filename are added to the end with /%s.
 */
static char* expand_path(char *dir, char *fn)
{
    char *cp;
    char *out = malloc(strlen(dir) + strlen(fn) + 1);

    char *out_writer = out;

    while ((cp = strchr(dir, '%')))
    {
        strncpy(out_writer, dir, cp - dir);
        out_writer += cp - dir;

        if (*++cp == 's')
        {
            strcpy(out_writer, fn);
            out_writer += strlen(fn);
            fn += strlen(fn);
            cp++;
        }
        else if (*cp >= '0' && *cp <= '9')
        {
            char *endp;
            long l;

            l = strtol(cp, &endp, 10);
            l = MIN(l, strlen(fn));
            if (*endp == 's')
            {
                strncpy(out_writer, fn, l);
                out_writer += l;
                fn += l;
                *out_writer = 0;
                cp = endp + 1;
            }
            else
            {
                *out_writer++ = '%';
                *out_writer++ = *cp++;
            }
        }
        else
        {
            *out_writer++ = '%';
            *out_writer++ = *cp++;
        }
        dir = cp;
    }
    strcpy(out_writer, dir);
    out_writer += strlen(dir);
    if (*fn && out_writer[-1] != '/')
        *out_writer++ = '/';
    strcpy(out_writer, fn);

    return out;
}

/*
 * Like expand_path, but doesn't use %[0-9]+s rules.
 */
static char* expand_path_basic(char *dir, char *fn)
{
    char *cp;
    char *out = malloc(strlen(dir) + strlen(fn) + 1);

    char *out_writer = out;

    while ((cp = strchr(dir, '%')))
    {
        strncpy(out_writer, dir, cp - dir);
        out_writer += cp - dir;

        if (*++cp == 's')
        {
            strcpy(out_writer, fn);
            out_writer += strlen(fn);
            fn += strlen(fn);
            cp++;
        }
        dir = cp;
    }
    strcpy(out_writer, dir);
    out_writer += strlen(dir);
    if (*fn && out_writer[-1] != '/')
        *out_writer++ = '/';
    strcpy(out_writer, fn);

    return out;
}

/*
 * Return an integer representation of pthread_self().
 */
static unsigned get_int_threadid()
{
    pthread_t pt = pthread_self();
    unsigned char *s = (unsigned char *)&pt;
    size_t i;
    unsigned h = 0;
    for (i = 0; i < sizeof(pthread_t); i++)
        h = (h << 5) - h + s[i];
    return h;
}

/*
 * Make the directory containing path and any prefix directories.
 */
static void mkdir_prefix(char *path, int mode)
{
    char *cp = strrchr(path, '/');
    if (!cp)
        return;

    *cp = 0;
    if (is_directory(path))
    {
        *cp = '/';
        return;
    }

    if (mkdir(path, mode) == 0)
    {
        chmod(path, mode);
        *cp = '/';
        return;
    }

    mkdir_prefix(path, mode);
    mkdir(path, mode);
    chmod(path, mode);
    *cp = '/';
}

#define READ_CHUNK_SIZE 8192

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
static char *tokenise_search_path(char *searchpath)
{
    char *newsearch;
    unsigned int i, j;
    size_t len;
#ifdef _WIN32
    char path_sep = ';';
#else
    char path_sep = ':';
#endif

    if (!searchpath)
        searchpath = "";

    newsearch = (char *)malloc((len = strlen(searchpath)) + 5);
    if (!newsearch)
        return NULL;

    for (i = 0, j = 0; i < len; i++)
    {
        /* "::" => ":". Used for escaping colons in http://foo */
        if (i < len - 1 && searchpath[i] == ':' && searchpath[i + 1] == ':')
        {
            newsearch[j++] = ':';
            i++;
            continue;
        }

        /* Handle http:// and ftp:// too without :: */
        if (path_sep == ':')
        {
            if ((i == 0 || (i > 0 && searchpath[i - 1] == ':')) &&
                (!strncmp(&searchpath[i], "http:",      5) ||
                 !strncmp(&searchpath[i], "https:",     6) ||
                 !strncmp(&searchpath[i], "ftp:",       4) ||
                 !strncmp(&searchpath[i], "|http:",     6) ||
                 !strncmp(&searchpath[i], "|https:",    7) ||
                 !strncmp(&searchpath[i], "|ftp:",      5) ||
                 !strncmp(&searchpath[i], "URL=http:",  9) ||
                 !strncmp(&searchpath[i], "URL=https:", 10)||
                 !strncmp(&searchpath[i], "URL=ftp:",   8)))
            {
                do
                {
                    newsearch[j++] = searchpath[i];
                } while (i < len && searchpath[i++] != ':');
                if (searchpath[i] == ':')
                    i++;
                if (searchpath[i] == '/')
                    newsearch[j++] = searchpath[i++];
                if (searchpath[i] == '/')
                    newsearch[j++] = searchpath[i++];
                // Look for host:port
                do
                {
                    newsearch[j++] = searchpath[i++];
                } while (i < len && searchpath[i] != ':' && searchpath[i] != '/');
                newsearch[j++] = searchpath[i++];
                if (searchpath[i] == ':')
                    i++;
            }
        }

        if (searchpath[i] == path_sep)
        {
            /* Skip blank path components */
            if (j && newsearch[j - 1] != 0)
                newsearch[j++] = 0;
        }
        else
        {
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

/*
 * Looks in a colon (in non-windows enviroments) or semi-colon (windows) separated
 * list (path) to find a file.
 *
 * Any path can contain %s subtitutions and filesystem subsitutions can contain %Ns
 * like subsitutions
 *
 * Returns a hFILE pointer when found.
 *           NULL otherwise.
 */
static char* resolve_file_in_path(char *file, char *path)
{
    char *newsearch;
    char *proposed_path;
    char* resolved_path;
    hFILE *fp;

    if (NULL == (newsearch = tokenise_search_path(path)))
        return NULL;

    for (proposed_path = newsearch; *proposed_path; proposed_path += strlen(proposed_path) + 1)
    {
        if (strncmp(proposed_path, "URL=", 4) == 0)
            resolved_path = expand_path_basic(proposed_path+4, file);
        else if (!strncmp(proposed_path, "http:", 5) ||
                 !strncmp(proposed_path, "https:", 5) ||
                 !strncmp(proposed_path, "ftp:", 4))
            resolved_path = expand_path_basic(proposed_path, file);
        else
            resolved_path = expand_path(proposed_path, file);

        // Does the file exist? Use the hFILE logic to find by opening the file, then closing it.
        if((fp = hopen(resolved_path, "r")))
        {
            hclose_abruptly(fp);
            free(newsearch);

            return resolved_path;
        }

        free(resolved_path);
    }

    free(newsearch);
    return NULL;
}

// Public functions

char* m5_to_path(const char *m5_str){
    char *ref_path = getenv("REF_PATH");
    char *ref_cache = getenv("REF_CACHE");
    char path_tmp[PATH_MAX];
    char cache[PATH_MAX], cache_root[PATH_MAX];

    cache_root[0] = '\0';

    if (!ref_path || *ref_path == '\0') {
        /*
        * If we have no ref path, we use the EBI server.
        * However to avoid spamming it we require a local ref cache too.
        */

        ref_path = "https://www.ebi.ac.uk/ena/cram/md5/%s";
        if (!ref_cache || *ref_cache == '\0') {
            const char *extra;
            const char *base = get_cache_basedir(&extra);
            snprintf(cache_root, PATH_MAX, "%s%s/hts-ref", base, extra);
            snprintf(cache, PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
            ref_cache = cache;
            hts_log_info("Populating local cache: %s", ref_cache);
        }
    }

    /* Try in REF_CACHE */
    if (ref_cache && *ref_cache) {
        struct stat sb;
        char* found_path = expand_path(ref_cache, (char *)m5_str);

        if(0 == stat(found_path, &sb)){
            return found_path;
        }
    }

    char* found_path;
    /* Try in REF_PATH */
    if (!(found_path = resolve_file_in_path((char*)m5_str, ref_path))) {
        hts_log_info("Failed to fetch file. REF_PATH: '%s', M5: '%s'", ref_path, m5_str);
        return NULL;
    }

    /* If the REF_CACHE enviromental variable is set, populate the cache. */
    if (ref_cache && *ref_cache) {
        hFILE* ref;
        if(!(ref = hopen(found_path, "r"))){
            return NULL;
        }

        int pid = (int)getpid();
        unsigned thrid = get_int_threadid();
        hFILE *fp;

        if (*cache_root && !is_directory(cache_root)) {
            hts_log_warning("Creating reference cache directory %s\n"
                            "This may become large; see the samtools(1) manual page REF_CACHE discussion",
                            cache_root);
        }

        char* cache_path = expand_path(ref_cache, (char *)m5_str);
        hts_log_info("Writing cache file '%s'", cache_path);
        mkdir_prefix(cache_path, 01777);

        do {
            // Attempt to further uniquify the temporary filename
            unsigned t = ((unsigned)time(NULL)) ^ ((unsigned)clock());
            thrid++; // Ensure filename changes even if time/clock haven't

            sprintf(path_tmp, "%s.tmp_%d_%u_%u", cache_path, pid, thrid, t);
            fp = hopen(path_tmp, "wx");
        } while (fp == NULL && errno == EEXIST);
        if (!fp) {
            perror(path_tmp);
            free(cache_path);

            // Doesn't matter if we can't write to the temp file, just return the non
            // cached path. This argument is used many times below.
            return found_path;
        }

        // Stream the file into the cache and check the md5
        hts_md5_context *md5;
        char unsigned md5_buf1[16];
        char md5_buf2[33];

        if (!(md5 = hts_md5_init())) {
            hclose_abruptly(fp);
            unlink(path_tmp);
            free(cache_path);
            hts_log_error("Function hts_md5_init failed");

            return found_path;
        }

        int read_length;
        char buf[READ_CHUNK_SIZE];

        while ((read_length = hread(ref, buf, READ_CHUNK_SIZE)) > 0) {
            hts_md5_update(md5, buf, read_length);
            if(hwrite(fp, buf, read_length) != read_length){
                perror(cache_path);
                hclose_abruptly(ref);
                free(cache_path);

                return found_path;
            }
        }

        hts_md5_final(md5_buf1, md5);
        hts_md5_destroy(md5);
        hts_md5_hex(md5_buf2, md5_buf1);

        if (strncmp(m5_str, md5_buf2, 32) != 0) {
            hclose_abruptly(fp);
            free(found_path);
            free(cache_path);
            unlink(path_tmp);

            hts_log_error("Mismatching md5sum for downloaded reference");
            return NULL;
        }

        if (hclose(fp) < 0)
        {
            perror(cache_path);
            unlink(path_tmp);
            free(cache_path);

            return found_path;
        }
        else
        {
            if (0 == chmod(path_tmp, 0444))
                rename(path_tmp, cache_path);
            else{
                perror(cache_path);
                unlink(path_tmp);
                free(cache_path);

                return found_path;
            }
        }

        free(found_path);
        return cache_path;
    }

    return found_path;
}

hFILE* m5_to_ref(const char *m5_str){
    char* m5_path;

    if (!(m5_path = m5_to_path(m5_str)))
        return NULL;

    hFILE* hf = hopen(m5_path, "r");
    free(m5_path);

    return hf;
}
