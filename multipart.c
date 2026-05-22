/*  multipart.c -- GA4GH redirection and multipart backend for file streams.

    Copyright (C) 2016-2017 Genome Research Ltd.

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

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "htslib/kstring.h"

#include "hts_internal.h"
#include "hfile_internal.h"

#ifndef EPROTO
#define EPROTO ENOEXEC
#endif
#ifndef ENOTSUP
#define ENOTSUP EINVAL
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif

typedef struct hfile_part {
    char *url;
    char **headers;
} hfile_part;

typedef struct {
    hFILE base;
    hfile_part *parts;
    size_t nparts, maxparts, current;
    hFILE *currentfp;
} hFILE_multipart;

static void free_part(hfile_part *p)
{
    free(p->url);
    if (p->headers) {
        char **hdr;
        for (hdr = p->headers; *hdr; hdr++) free(*hdr);
        free(p->headers);
    }

    p->url = NULL;
    p->headers = NULL;
}

static void free_all_parts(hFILE_multipart *fp)
{
    size_t i;
    for (i = 0; i < fp->nparts; i++) free_part(&fp->parts[i]);
    free(fp->parts);
}

static ssize_t multipart_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_multipart *fp = (hFILE_multipart *) fpv;
    size_t n;

open_next:
    if (fp->currentfp == NULL) {
        if (fp->current < fp->nparts) {
            const hfile_part *p = &fp->parts[fp->current];
            hts_log_debug("Opening part #%zu of %zu: \"%.120s%s\"",
                fp->current+1, fp->nparts, p->url,
                (strlen(p->url) > 120)? "..." : "");

            fp->currentfp = p->headers?
                  hopen(p->url, "r:",
                        "httphdr:v", p->headers,
                        "auth_token_enabled", "false", NULL)
                : hopen(p->url, "r:", "auth_token_enabled", "false", NULL);

            if (fp->currentfp == NULL) return -1;
        }
        else return 0;  // No more parts, so we're truly at EOF
    }

    n = fp->currentfp->mobile?
          fp->currentfp->backend->read(fp->currentfp, buffer, nbytes)
        : hread(fp->currentfp, buffer, nbytes);

    if (n == 0) {
        // We're at EOF on this part, so set up the next part
        hFILE *prevfp = fp->currentfp;
        free_part(&fp->parts[fp->current]);
        fp->current++;
        fp->currentfp = NULL;
        if (hclose(prevfp) < 0) return -1;
        goto open_next;
    }

    return n;  // Number of bytes read by (or an error from) fp->currentfp
}

static ssize_t multipart_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    errno = EROFS;
    return -1;
}

static off_t multipart_seek(hFILE *fpv, off_t offset, int whence)
{
    errno = ESPIPE;
    return -1;
}

static int multipart_close(hFILE *fpv)
{
    hFILE_multipart *fp = (hFILE_multipart *) fpv;

    free_all_parts(fp);
    if (fp->currentfp) {
        if (hclose(fp->currentfp) < 0) return -1;
    }

    return 0;
}

static const struct hFILE_backend multipart_backend =
{
    multipart_read, multipart_write, multipart_seek, NULL, multipart_close
};

typedef struct hfile_chunked_part {
    char *url;
    off_t offset;
    off_t length;
} hfile_chunked_part;

typedef struct {
    hFILE base;
    hfile_chunked_part *parts;
    size_t nparts, maxparts, current;
    hFILE *currentfp;
    off_t length, position;
} hFILE_chunked;

static int chunked_has_uri_scheme(const char *s)
{
    int n = 0;

    while (isalnum_c(s[n]) || s[n] == '+' || s[n] == '-' || s[n] == '.')
        n++;

    return n > 1 && s[n] == ':';
}

static int chunked_is_absolute_path(const char *s)
{
    if (s[0] == '/')
        return 1;
#if defined(_WIN32) || defined(__MSYS__)
    if (isalpha_c(s[0]) && s[1] == ':' && (s[2] == '/' || s[2] == '\\'))
        return 1;
#endif
    return 0;
}

static char *chunked_manifest_dir(const char *manifest_name)
{
    const char *slash;
    char *dir;
    size_t len;

    if (chunked_has_uri_scheme(manifest_name))
        return NULL;

    slash = strrchr(manifest_name, '/');
    if (!slash)
        return NULL;

    len = slash - manifest_name + 1;
    dir = malloc(len + 1);
    if (!dir) {
        errno = ENOMEM;
        return NULL;
    }

    memcpy(dir, manifest_name, len);
    dir[len] = '\0';
    return dir;
}

static char *chunked_resolve_name(const char *base, const char *name)
{
    char *copy;

    if (!base || chunked_is_absolute_path(name) || chunked_has_uri_scheme(name))
        return strdup(name);

    size_t base_len = strlen(base), name_len = strlen(name);
    copy = malloc(base_len + name_len + 1);
    if (!copy)
        return NULL;

    memcpy(copy, base, base_len);
    memcpy(copy + base_len, name, name_len + 1);
    return copy;
}

static void chunked_free_parts(hFILE_chunked *fp)
{
    size_t i;

    for (i = 0; i < fp->nparts; i++)
        free(fp->parts[i].url);
    free(fp->parts);
}

static int chunked_close_current(hFILE_chunked *fp)
{
    if (fp->currentfp) {
        hFILE *part = fp->currentfp;
        fp->currentfp = NULL;
        if (hclose(part) < 0)
            return -1;
    }

    return 0;
}

static size_t chunked_find_part(hFILE_chunked *fp, off_t offset)
{
    size_t lo = 0, hi = fp->nparts;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        off_t start = fp->parts[mid].offset;
        off_t end = start + fp->parts[mid].length;

        if (offset < start)
            hi = mid;
        else if (offset >= end)
            lo = mid + 1;
        else
            return mid;
    }

    return lo;
}

static int chunked_open_current(hFILE_chunked *fp)
{
    hfile_chunked_part *part;
    off_t offset;

    if (fp->current >= fp->nparts)
        return 0;
    if (fp->currentfp)
        return 0;

    part = &fp->parts[fp->current];
    offset = fp->position - part->offset;

    fp->currentfp = hopen(part->url, "r");
    if (!fp->currentfp)
        return -1;

    if (offset != 0 && hseek(fp->currentfp, offset, SEEK_SET) < 0) {
        int save = errno;
        hclose_abruptly(fp->currentfp);
        fp->currentfp = NULL;
        errno = save;
        return -1;
    }

    return 0;
}

static ssize_t chunked_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_chunked *fp = (hFILE_chunked *) fpv;

    if (nbytes == 0)
        return 0;

    while (fp->current < fp->nparts) {
        hfile_chunked_part *part = &fp->parts[fp->current];
        off_t end = part->offset + part->length;
        off_t remaining;
        size_t request = nbytes;
        ssize_t n;

        if (part->length == 0 || fp->position >= end) {
            if (chunked_close_current(fp) < 0)
                return -1;
            fp->current++;
            continue;
        }

        if (fp->position < part->offset) {
            errno = EIO;
            return -1;
        }

        if (chunked_open_current(fp) < 0)
            return -1;

        remaining = end - fp->position;
        if (remaining < (off_t) request)
            request = remaining;

        n = hread(fp->currentfp, buffer, request);
        if (n < 0)
            return -1;
        if (n == 0) {
            errno = EIO;
            return -1;
        }

        fp->position += n;
        return n;
    }

    return 0;
}

static ssize_t chunked_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    errno = EROFS;
    return -1;
}

static off_t chunked_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_chunked *fp = (hFILE_chunked *) fpv;
    off_t base, target;
    size_t part;

    switch (whence) {
    case SEEK_SET:
        base = 0;
        break;
    case SEEK_CUR:
        base = htell(&fp->base);
        break;
    case SEEK_END:
        base = fp->length;
        break;
    default:
        errno = EINVAL;
        return -1;
    }

    if ((offset < 0 && offset < -base)
        || (offset > 0 && offset > fp->length - base)) {
        errno = EINVAL;
        return -1;
    }
    target = base + offset;

    part = target == fp->length ? fp->nparts : chunked_find_part(fp, target);
    if (part == fp->nparts && target != fp->length) {
        errno = EINVAL;
        return -1;
    }

    if (part != fp->current && chunked_close_current(fp) < 0)
        return -1;

    fp->current = part;
    fp->position = target;

    if (fp->currentfp) {
        off_t part_offset = target - fp->parts[fp->current].offset;
        if (hseek(fp->currentfp, part_offset, SEEK_SET) < 0)
            return -1;
    }

    return target;
}

static int chunked_close(hFILE *fpv)
{
    hFILE_chunked *fp = (hFILE_chunked *) fpv;
    int ret = chunked_close_current(fp);

    chunked_free_parts(fp);
    return ret;
}

static const struct hFILE_backend chunked_backend =
{
    chunked_read, chunked_write, chunked_seek, NULL, chunked_close
};

static int chunked_part_length(const char *url, off_t *length)
{
    hFILE *part = hopen(url, "r");
    off_t len;

    if (!part)
        return -1;

    len = hseek(part, 0, SEEK_END);
    if (len < 0) {
        int save = errno;
        hclose_abruptly(part);
        errno = save;
        return -1;
    }

    if (hclose(part) < 0)
        return -1;

    *length = len;
    return 0;
}

hFILE *hopen_chunked_manifest(const char *url, const char *mode)
{
    hFILE_chunked *fp = NULL;
    hFILE *manifest = NULL;
    kstring_t line = KS_INITIALIZE;
    const char *manifest_name = url + 8; // len("chunked:") = 8
    char *manifest_dir = NULL;

    if (!strchr(mode, 'r') || strchr(mode, '+') || strchr(mode, 'w')
        || strchr(mode, 'a')) {
        errno = ENOTSUP;
        return NULL;
    }

    if (*manifest_name == '\0') {
        errno = EINVAL;
        return NULL;
    }

    fp = (hFILE_chunked *) hfile_init(sizeof(*fp), mode, 0);
    if (!fp) return NULL;

    fp->parts = NULL;
    fp->nparts = fp->maxparts = 0;
    fp->current = 0;
    fp->currentfp = NULL;
    fp->length = fp->position = 0;

    manifest = hopen(manifest_name, "r");
    if (!manifest)
        goto fail;

    errno = 0;
    manifest_dir = chunked_manifest_dir(manifest_name);
    if (!manifest_dir && errno)
        goto fail;

    while (ks_clear(&line), khgetline(&line, manifest) == 0) {
        char *chunk_name;
        off_t chunk_length;
        hfile_chunked_part *part;

        if (line.l == 0 || line.s[0] == '#')
            continue;

        chunk_name = chunked_resolve_name(manifest_dir, line.s);
        if (!chunk_name)
            goto fail;

        if (chunked_part_length(chunk_name, &chunk_length) < 0) {
            free(chunk_name);
            goto fail;
        }
        if (chunk_length < 0 || chunk_length > INT64_MAX - fp->length) {
            free(chunk_name);
            errno = EOVERFLOW;
            goto fail;
        }

        hts_expand(hfile_chunked_part, fp->nparts + 1, fp->maxparts, fp->parts);
        part = &fp->parts[fp->nparts++];
        part->url = chunk_name;
        part->offset = fp->length;
        part->length = chunk_length;
        fp->length += chunk_length;
    }

    if (herrno(manifest))
        goto fail;
    if (hclose(manifest) < 0) {
        manifest = NULL;
        goto fail;
    }
    manifest = NULL;

    if (fp->nparts == 0) {
        errno = EINVAL;
        goto fail;
    }

    fp->base.backend = &chunked_backend;

    free(manifest_dir);
    ks_free(&line);
    return &fp->base;

 fail:
    {
        int save = errno;
        if (manifest)
            hclose_abruptly(manifest);
        free(manifest_dir);
        ks_free(&line);
        chunked_free_parts(fp);
        hfile_destroy((hFILE *) fp);
        errno = save;
    }
    return NULL;
}

// Returns 'v' (valid value), 'i' (invalid; required GA4GH field missing),
// or upon encountering an unexpected token, that token's type.
// Explicit `return '?'` means a JSON parsing error, typically a member key
// that is not a string.  An unexpected token may be a valid token that was
// not the type expected for a particular GA4GH field, or it may be '?' or
// '\0' which should be propagated.
static char
parse_ga4gh_body_json(hFILE_multipart *fp, hFILE *json,
                      kstring_t *b, kstring_t *header)
{
    hts_json_token t;

    if (hts_json_fnext(json, &t, b) != '{') return t.type;
    while (hts_json_fnext(json, &t, b) != '}') {
        if (t.type != 's') return '?';

        if (strcmp(t.str, "urls") == 0) {
            if (hts_json_fnext(json, &t, b) != '[') return t.type;

            while (hts_json_fnext(json, &t, b) != ']') {
                hfile_part *part;
                size_t n = 0, max = 0;

                hts_expand(hfile_part, fp->nparts+1, fp->maxparts, fp->parts);
                part = &fp->parts[fp->nparts++];
                part->url = NULL;
                part->headers = NULL;

                if (t.type != '{') return t.type;
                while (hts_json_fnext(json, &t, b) != '}') {
                    if (t.type != 's') return '?';

                    if (strcmp(t.str, "url") == 0) {
                        if (hts_json_fnext(json, &t, b) != 's') return t.type;
                        part->url = ks_release(b);
                    }
                    else if (strcmp(t.str, "headers") == 0) {
                        if (hts_json_fnext(json, &t, b) != '{') return t.type;

                        while (hts_json_fnext(json, &t, header) != '}') {
                            if (t.type != 's') return '?';

                            if (hts_json_fnext(json, &t, b) != 's')
                                return t.type;

                            kputs(": ", header);
                            kputs(t.str, header);
                            n++;
                            hts_expand(char *, n+1, max, part->headers);
                            part->headers[n-1] = ks_release(header);
                            part->headers[n] = NULL;
                        }
                    }
                    else if (hts_json_fskip_value(json, '\0') != 'v')
                        return '?';
                }

                if (! part->url) return 'i';
            }
        }
        else if (strcmp(t.str, "format") == 0) {
            if (hts_json_fnext(json, &t, b) != 's') return t.type;

            hts_log_debug("GA4GH JSON redirection to multipart %s data", t.str);
        }
        else if (hts_json_fskip_value(json, '\0') != 'v') return '?';
    }

    return 'v';
}

// Returns 'v' (valid value), 'i' (invalid; required GA4GH field missing),
// or upon encountering an unexpected token, that token's type.
// Explicit `return '?'` means a JSON parsing error, typically a member key
// that is not a string.  An unexpected token may be a valid token that was
// not the type expected for a particular GA4GH field, or it may be '?' or
// '\0' which should be propagated.
static char
parse_ga4gh_redirect_json(hFILE_multipart *fp, hFILE *json,
                          kstring_t *b, kstring_t *header) {
    hts_json_token t;

    if (hts_json_fnext(json, &t, b) != '{') return t.type;
    while (hts_json_fnext(json, &t, b) != '}') {
        if (t.type != 's') return '?';

        if (strcmp(t.str, "htsget") == 0) {
            char ret = parse_ga4gh_body_json(fp, json, b, header);
            if (ret != 'v') return ret;
        }
        else return '?';
    }

    if (hts_json_fnext(json, &t, b) != '\0') return '?';

    return 'v';
}

hFILE *hopen_htsget_redirect(hFILE *hfile, const char *mode)
{
    hFILE_multipart *fp;
    kstring_t s1 = { 0, 0, NULL }, s2 = { 0, 0, NULL };
    char ret;

    fp = (hFILE_multipart *) hfile_init(sizeof (hFILE_multipart), mode, 0);
    if (fp == NULL) return NULL;

    fp->parts = NULL;
    fp->nparts = fp->maxparts = 0;

    ret = parse_ga4gh_redirect_json(fp, hfile, &s1, &s2);
    free(s1.s);
    free(s2.s);
    if (ret != 'v') {
        free_all_parts(fp);
        hfile_destroy((hFILE *) fp);
        errno = (ret == '?' || ret == '\0')? EPROTO : EINVAL;
        return NULL;
    }

    fp->current = 0;
    fp->currentfp = NULL;
    fp->base.backend = &multipart_backend;
    return &fp->base;
}
