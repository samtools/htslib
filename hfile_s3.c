/*  hfile_s3.c -- Amazon S3 backend for low-level file streams.

    Copyright (C) 2015-2017 Genome Research Ltd.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "hts_internal.h"
#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"

#if defined HAVE_COMMONCRYPTO

#include <CommonCrypto/CommonHMAC.h>

#define DIGEST_BUFSIZ CC_SHA1_DIGEST_LENGTH

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    CCHmac(kCCHmacAlgSHA1, key->s, key->l, message->s, message->l, digest);
    return CC_SHA1_DIGEST_LENGTH;
}

#elif defined HAVE_HMAC

#include <openssl/hmac.h>

#define DIGEST_BUFSIZ EVP_MAX_MD_SIZE

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    unsigned int len;
    HMAC(EVP_sha1(), key->s, key->l,
         (unsigned char *) message->s, message->l, digest, &len);
    return len;
}

#else
#error No HMAC() routine found by configure
#endif

static void
urldecode_kput(const char *s, int len, kstring_t *str)
{
    char buf[3];
    int i = 0;

    while (i < len)
        if (s[i] == '%' && i+2 < len) {
            buf[0] = s[i+1], buf[1] = s[i+2], buf[2] = '\0';
            kputc(strtol(buf, NULL, 16), str);
            i += 3;
        }
        else kputc(s[i++], str);
}

static void base64_kput(const unsigned char *data, size_t len, kstring_t *str)
{
    static const char base64[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    size_t i = 0;
    unsigned x = 0;
    int bits = 0, pad = 0;

    while (bits || i < len) {
        if (bits < 6) {
            x <<= 8, bits += 8;
            if (i < len) x |= data[i++];
            else pad++;
        }

        bits -= 6;
        kputc(base64[(x >> bits) & 63], str);
    }

    str->l -= pad;
    kputsn("==", pad, str);
}

static int is_dns_compliant(const char *s0, const char *slim)
{
    int has_nondigit = 0, len = 0;
    const char *s;

    for (s = s0; s < slim; len++, s++)
        if (islower_c(*s))
            has_nondigit = 1;
        else if (*s == '-') {
            has_nondigit = 1;
            if (s == s0 || s+1 == slim) return 0;
        }
        else if (isdigit_c(*s))
            ;
        else if (*s == '.') {
            if (s == s0 || ! isalnum_c(s[-1])) return 0;
            if (s+1 == slim || ! isalnum_c(s[1])) return 0;
        }
        else return 0;

    return has_nondigit && len >= 3 && len <= 63;
}

static FILE *expand_tilde_open(const char *fname, const char *mode)
{
    FILE *fp;

    if (strncmp(fname, "~/", 2) == 0) {
        kstring_t full_fname = { 0, 0, NULL };
        const char *home = getenv("HOME");
        if (! home) return NULL;

        kputs(home, &full_fname);
        kputs(&fname[1], &full_fname);

        fp = fopen(full_fname.s, mode);
        free(full_fname.s);
    }
    else
        fp = fopen(fname, mode);

    return fp;
}

static void parse_ini(const char *fname, const char *section, ...)
{
    kstring_t line = { 0, 0, NULL };
    int active = 1;  // Start active, so global properties are accepted
    char *s;

    FILE *fp = expand_tilde_open(fname, "r");
    if (fp == NULL) return;

    while (line.l = 0, kgetline(&line, (kgets_func *) fgets, fp) >= 0)
        if (line.s[0] == '[' && (s = strchr(line.s, ']')) != NULL) {
            *s = '\0';
            active = (strcmp(&line.s[1], section) == 0);
        }
        else if (active && (s = strpbrk(line.s, ":=")) != NULL) {
            const char *key = line.s, *value = &s[1], *akey;
            va_list args;

            while (isspace_c(*key)) key++;
            while (s > key && isspace_c(s[-1])) s--;
            *s = '\0';

            while (isspace_c(*value)) value++;
            while (line.l > 0 && isspace_c(line.s[line.l-1]))
                line.s[--line.l] = '\0';

            va_start(args, section);
            while ((akey = va_arg(args, const char *)) != NULL) {
                kstring_t *avar = va_arg(args, kstring_t *);
                if (strcmp(key, akey) == 0) { kputs(value, avar); break; }
            }
            va_end(args);
        }

    fclose(fp);
    free(line.s);
}

static void parse_simple(const char *fname, kstring_t *id, kstring_t *secret)
{
    kstring_t text = { 0, 0, NULL };
    char *s;
    size_t len;

    FILE *fp = expand_tilde_open(fname, "r");
    if (fp == NULL) return;

    while (kgetline(&text, (kgets_func *) fgets, fp) >= 0)
        kputc(' ', &text);
    fclose(fp);

    s = text.s;
    while (isspace_c(*s)) s++;
    kputsn(s, len = strcspn(s, " \t"), id);

    s += len;
    while (isspace_c(*s)) s++;
    kputsn(s, strcspn(s, " \t"), secret);

    free(text.s);
}

static hFILE * s3_rewrite(const char *s3url, const char *mode, va_list *argsp)
{
    const char *bucket, *path;
    char date_hdr[40];
    char *header_list[4], **header = header_list;

    kstring_t message = { 0, 0, NULL };
    kstring_t url = { 0, 0, NULL };
    kstring_t profile = { 0, 0, NULL };
    kstring_t id = { 0, 0, NULL };
    kstring_t secret = { 0, 0, NULL };
    kstring_t host_base = { 0, 0, NULL };
    kstring_t token = { 0, 0, NULL };
    kstring_t token_hdr = { 0, 0, NULL };
    kstring_t auth_hdr = { 0, 0, NULL };

    time_t now = time(NULL);
#ifdef HAVE_GMTIME_R
    struct tm tm_buffer;
    struct tm *tm = gmtime_r(&now, &tm_buffer);
#else
    struct tm *tm = gmtime(&now);
#endif

    kputs(strchr(mode, 'r')? "GET\n" : "PUT\n", &message);
    kputc('\n', &message);
    kputc('\n', &message);
    strftime(date_hdr, sizeof date_hdr, "Date: %a, %d %b %Y %H:%M:%S GMT", tm);
    *header++ = date_hdr;
    kputs(&date_hdr[6], &message);
    kputc('\n', &message);

    // Our S3 URL format is s3[+SCHEME]://[ID[:SECRET[:TOKEN]]@]BUCKET/PATH

    if (s3url[2] == '+') {
        bucket = strchr(s3url, ':') + 1;
        kputsn(&s3url[3], bucket - &s3url[3], &url);
    }
    else {
        kputs("https:", &url);
        bucket = &s3url[3];
    }
    while (*bucket == '/') kputc(*bucket++, &url);

    path = bucket + strcspn(bucket, "/?#@");
    if (*path == '@') {
        const char *colon = strpbrk(bucket, ":@");
        if (*colon != ':') {
            urldecode_kput(bucket, colon - bucket, &profile);
        }
        else {
            const char *colon2 = strpbrk(&colon[1], ":@");
            urldecode_kput(bucket, colon - bucket, &id);
            urldecode_kput(&colon[1], colon2 - &colon[1], &secret);
            if (*colon2 == ':')
                urldecode_kput(&colon2[1], path - &colon2[1], &token);
        }

        bucket = &path[1];
        path = bucket + strcspn(bucket, "/?#");
    }
    else {
        // If the URL has no ID[:SECRET]@, consider environment variables.
        const char *v;
        if ((v = getenv("AWS_ACCESS_KEY_ID")) != NULL) kputs(v, &id);
        if ((v = getenv("AWS_SECRET_ACCESS_KEY")) != NULL) kputs(v, &secret);
        if ((v = getenv("AWS_SESSION_TOKEN")) != NULL) kputs(v, &token);

        if ((v = getenv("AWS_DEFAULT_PROFILE")) != NULL) kputs(v, &profile);
        else if ((v = getenv("AWS_PROFILE")) != NULL) kputs(v, &profile);
        else kputs("default", &profile);
    }

    if (id.l == 0) {
        const char *v = getenv("AWS_SHARED_CREDENTIALS_FILE");
        parse_ini(v? v : "~/.aws/credentials", profile.s,
                  "aws_access_key_id", &id, "aws_secret_access_key", &secret,
                  "aws_session_token", &token, NULL);
    }
    if (id.l == 0)
        parse_ini("~/.s3cfg", profile.s, "access_key", &id,
                  "secret_key", &secret, "access_token", &token,
                  "host_base", &host_base, NULL);
    if (id.l == 0)
        parse_simple("~/.awssecret", &id, &secret);

    if (host_base.l == 0)
        kputs("s3.amazonaws.com", &host_base);
    // Use virtual hosted-style access if possible, otherwise path-style.
    if (is_dns_compliant(bucket, path)) {
        kputsn(bucket, path - bucket, &url);
        kputc('.', &url);
        kputs(host_base.s, &url);
    }
    else {
        kputs(host_base.s, &url);
        kputc('/', &url);
        kputsn(bucket, path - bucket, &url);
    }
    kputs(path, &url);

    if (token.l > 0) {
        kputs("x-amz-security-token:", &message);
        kputs(token.s, &message);
        kputc('\n', &message);

        kputs("X-Amz-Security-Token: ", &token_hdr);
        kputs(token.s, &token_hdr);
        *header++ = token_hdr.s;
    }

    kputc('/', &message);
    kputs(bucket, &message); // CanonicalizedResource is '/' + bucket + path

    // If we have no id/secret, we can't sign the request but will
    // still be able to access public data sets.
    if (id.l > 0 && secret.l > 0) {
        unsigned char digest[DIGEST_BUFSIZ];
        size_t digest_len = s3_sign(digest, &secret, &message);

        kputs("Authorization: AWS ", &auth_hdr);
        kputs(id.s, &auth_hdr);
        kputc(':', &auth_hdr);
        base64_kput(digest, digest_len, &auth_hdr);

        *header++ = auth_hdr.s;
    }

    *header = NULL;
    hFILE *fp = hopen(url.s, mode, "va_list", argsp, "httphdr:v", header_list,
                      NULL);
    free(message.s);
    free(url.s);
    free(profile.s);
    free(id.s);
    free(secret.s);
    free(host_base.s);
    free(token.s);
    free(token_hdr.s);
    free(auth_hdr.s);
    return fp;
}

static hFILE *s3_open(const char *url, const char *mode)
{
    kstring_t mode_colon = { 0, 0, NULL };
    kputs(mode, &mode_colon);
    kputc(':', &mode_colon);
    hFILE *fp = s3_rewrite(url, mode_colon.s, NULL);
    free(mode_colon.s);
    return fp;
}

static hFILE *s3_vopen(const char *url, const char *mode_colon, va_list args0)
{
    // Need to use va_copy() as we can only take the address of an actual
    // va_list object, not that of a parameter whose type may have decayed.
    va_list args;
    va_copy(args, args0);
    hFILE *fp = s3_rewrite(url, mode_colon, &args);
    va_end(args);
    return fp;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_s3)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { s3_open, hfile_always_remote, "Amazon S3", 2000 + 50, s3_vopen
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_s3 plugin (htslib)\t" HTS_VERSION;
    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_s3.init] version %s\n", strchr(id, '\t')+1);
#endif

    self->name = "Amazon S3";
    hfile_add_scheme_handler("s3", &handler);
    hfile_add_scheme_handler("s3+http", &handler);
    hfile_add_scheme_handler("s3+https", &handler);
    return 0;
}
