/*  hfile_tcga.c -- TCGA Genomic Data Commons backend for low-level file streams.

    Copyright (C) 2017 Lexent Bio, Inc.

    Author: Ram Yalamanchili <ramy@lexentbio.com>

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/knetfile.h"
#include "hfile_internal.h"
#include "hts_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif

#define GDC_DOWNLOAD_URL "https://gdc-api.nci.nih.gov/data/"
#define GDC_METADATA_URL "https://gdc-api.nci.nih.gov/v0/legacy/files/%s?expand=index_files"
#define INDEX_FMT ".bai"

static char *parse_tcga_index_id(char *json)
{
    hts_json_token token;
    size_t state = 0;
    char v;

    v = hts_json_snext(json, &state, &token);
    // Look for json['data']['index_files'][0]['file_id']
    while (v) {
      if (v == 's' && strcmp(token.str, "file_id") == 0) {
        v = hts_json_snext(json, &state, &token);
        return token.str;
      }
      v = hts_json_snext(json, &state, &token);
    }

    return NULL;
}

static hFILE *
tcga_rewrite(const char *tcga_url, const char *mode, int mode_has_colon,
            va_list *argsp)
{
    int url_len;
    const char *uuid, *path, *access_token, *bai_idx;
    kstring_t mode_colon = { 0, 0, NULL };
    kstring_t url = { 0, 0, NULL };
    kstring_t auth_hdr = { 0, 0, NULL };
    hFILE *fp = NULL;

    if (!strchr(mode, 'r')) {
        fprintf(stderr, "[M::tcga_open] unsupported mode: %s\b", mode);
        return NULL;
    }

    uuid = &tcga_url[5];
    while (*uuid == '/') uuid++;

    if (strlen(uuid) != strcspn(uuid, "/?#")) {
      fprintf(stderr, "[M::tcga_open] Invalid URL: %s, must be of format tcga://UUID\n", tcga_url);
      return NULL;
    }

    // TCGA URL format is tcga://UUID or tcga://UUID.bai
    url_len = strlen(tcga_url);
    bai_idx = strstr(tcga_url, INDEX_FMT);

    // Requested TCGA index file; Look up corresponding UUID and return TCGA URL
    if (bai_idx && (bai_idx - tcga_url) == (url_len - strlen(INDEX_FMT))) {
      kstring_t uuid_s = { 0, 0, NULL };
      kstring_t index_url = { 0, 0, NULL };
      hFILE *index_url_h;
      int json_buffer_len = 0x8000;

      kputsn(uuid, (bai_idx - uuid), &uuid_s);
      ksprintf(&index_url, GDC_METADATA_URL, uuid_s.s);
      index_url_h = hopen(index_url.s, "r", NULL);
      free(uuid_s.s);
      free(index_url.s);

      if (index_url_h == NULL) {
        fprintf(stderr, "[M::tcga_open] Unable to fetch index file for %s\n", tcga_url);
      } else {
        char *buffer = malloc(json_buffer_len);

        if(hread(index_url_h, buffer, json_buffer_len) > 0) {
          char *index_id = parse_tcga_index_id(buffer);
          if (index_id != NULL) {
            kputs(GDC_DOWNLOAD_URL, &url);
            kputs(index_id, &url);
          }
        }
        free(buffer);
      }
    } else {
      // Return download URL for corresponding TCGA UUID
      kputs(GDC_DOWNLOAD_URL, &url);
      path = uuid + strcspn(uuid, "/?#");
      kputsn(uuid, path - uuid, &url);
    }

    if (url.s == NULL) {
      goto error;
    }

    if (hts_verbose >= 8)
        fprintf(stderr, "[M::tcga_open] rewrote URL as %s\n", url.s);

    // TODO Find the access token in a more standard way
    access_token = getenv("TCGA_AUTH_TOKEN");

    if (access_token == NULL) {
      fprintf(stderr, "[M::tcga_open] FATAL: TCGA_AUTH_TOKEN environment variable not found.\n");
      goto error;
    }

    kputs("X-Auth-Token: ", &auth_hdr);
    kputs(access_token, &auth_hdr);

    if (argsp || auth_hdr.l > 0 || mode_has_colon) {
        if (! mode_has_colon) {
            kputs(mode, &mode_colon);
            kputc(':', &mode_colon);
            mode = mode_colon.s;
        }

        fp = hopen(url.s, mode, "va_list", argsp,
                   "httphdr", (auth_hdr.l > 0)? auth_hdr.s : NULL, NULL);
    }
    else
        fp = hopen(url.s, mode);

    free(url.s);
error:
    free(mode_colon.s);
    free(auth_hdr.s);
    return fp;
}

static hFILE *tcga_open(const char *url, const char *mode)
{
    return tcga_rewrite(url, mode, 0, NULL);
}

static hFILE *tcga_vopen(const char *url, const char *mode_colon, va_list args0)
{
    // Need to use va_copy() as we can only take the address of an actual
    // va_list object, not that of a parameter as its type may have decayed.
    va_list args;
    va_copy(args, args0);
    hFILE *fp = tcga_rewrite(url, mode_colon, 1, &args);
    va_end(args);
    return fp;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_tcga)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { tcga_open, hfile_always_remote, "TCGA Genomic Data Commons",
          2000 + 50, tcga_vopen
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_tcga plugin (htslib)\t" HTS_VERSION;
    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_tcga.init] version %s\n", strchr(id, '\t')+1);
#endif

    self->name = "TCGA Genomic Data Commons";
    hfile_add_scheme_handler("tcga", &handler);
    return 0;
}
