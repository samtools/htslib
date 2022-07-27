/*  hfile_s3.c -- Amazon S3 backend for low-level file streams.

    Copyright (C) 2015-2017, 2019-2022 Genome Research Ltd.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#include <errno.h>

#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"
#include "hts_time_funcs.h"

typedef struct s3_auth_data {
    kstring_t id;
    kstring_t token;
    kstring_t secret;
    kstring_t region;
    kstring_t canonical_query_string;
    kstring_t user_query_string;
    kstring_t host;
    kstring_t profile;
    time_t creds_expiry_time;
    char *bucket;
    kstring_t auth_hdr;
    time_t auth_time;
    char date[40];
    char date_long[17];
    char date_short[9];
    kstring_t date_html;
    char mode;
    char *headers[5];
    int refcount;
} s3_auth_data;

#define AUTH_LIFETIME 60  // Regenerate auth headers if older than this
#define CREDENTIAL_LIFETIME 60 // Seconds before expiry to reread credentials

#if defined HAVE_COMMONCRYPTO

#include <CommonCrypto/CommonHMAC.h>

#define DIGEST_BUFSIZ CC_SHA1_DIGEST_LENGTH
#define SHA256_DIGEST_BUFSIZE CC_SHA256_DIGEST_LENGTH
#define HASH_LENGTH_SHA256 (SHA256_DIGEST_BUFSIZE * 2) + 1

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    CCHmac(kCCHmacAlgSHA1, key->s, key->l, message->s, message->l, digest);
    return CC_SHA1_DIGEST_LENGTH;
}


static void s3_sha256(const unsigned char *in, size_t length, unsigned char *out) {
    CC_SHA256(in, length, out);
}


static void s3_sign_sha256(const void *key, int key_len, const unsigned char *d, int n, unsigned char *md, unsigned int *md_len) {
    CCHmac(kCCHmacAlgSHA256, key, key_len, d, n, md);
    *md_len = CC_SHA256_DIGEST_LENGTH;
}


#elif defined HAVE_HMAC

#include <openssl/hmac.h>
#include <openssl/sha.h>

#define DIGEST_BUFSIZ EVP_MAX_MD_SIZE
#define SHA256_DIGEST_BUFSIZE SHA256_DIGEST_LENGTH
#define HASH_LENGTH_SHA256 (SHA256_DIGEST_BUFSIZE * 2) + 1

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    unsigned int len;
    HMAC(EVP_sha1(), key->s, key->l,
         (unsigned char *) message->s, message->l, digest, &len);
    return len;
}


static void s3_sha256(const unsigned char *in, size_t length, unsigned char *out) {
    SHA256(in, length, out);
}


static void s3_sign_sha256(const void *key, int key_len, const unsigned char *d, int n, unsigned char *md, unsigned int *md_len) {
    HMAC(EVP_sha256(), key, key_len, d, n, md, md_len);
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

static int is_dns_compliant(const char *s0, const char *slim, int is_https)
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
            if (is_https) return 0;
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
                if (strcmp(key, akey) == 0) {
                    avar->l = 0;
                    kputs(value, avar);
                    break; }
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

static int copy_auth_headers(s3_auth_data *ad, char ***hdrs) {
    char **hdr = &ad->headers[0];
    int idx = 0;
    *hdrs = hdr;

    hdr[idx] = strdup(ad->date);
    if (!hdr[idx]) return -1;
    idx++;

    if (ad->token.l) {
        kstring_t token_hdr = KS_INITIALIZE;
        kputs("X-Amz-Security-Token: ", &token_hdr);
        kputs(ad->token.s, &token_hdr);
        if (token_hdr.s) {
            hdr[idx++] = token_hdr.s;
        } else {
            goto fail;
        }
    }

    if (ad->auth_hdr.l) {
        hdr[idx] = strdup(ad->auth_hdr.s);
        if (!hdr[idx]) goto fail;
        idx++;
    }

    hdr[idx] = NULL;
    return 0;

 fail:
    for (--idx; idx >= 0; --idx)
        free(hdr[idx]);
    return -1;
}

static void free_auth_data(s3_auth_data *ad) {
    if (ad->refcount > 0) {
        --ad->refcount;
        return;
    }
    free(ad->profile.s);
    free(ad->id.s);
    free(ad->token.s);
    free(ad->secret.s);
    free(ad->region.s);
    free(ad->canonical_query_string.s);
    free(ad->user_query_string.s);
    free(ad->host.s);
    free(ad->bucket);
    free(ad->auth_hdr.s);
    free(ad->date_html.s);
    free(ad);
}

static time_t parse_rfc3339_date(kstring_t *datetime)
{
    int offset = 0;
    time_t when;
    int num;
    char should_be_t = '\0', timezone[10] = { '\0' };
    unsigned int year, mon, day, hour, min, sec;

    if (!datetime->s)
        return 0;

    // It should be possible to do this with strptime(), but it seems
    // to not get on with our feature definitions.
    num = sscanf(datetime->s, "%4u-%2u-%2u%c%2u:%2u:%2u%9s",
                 &year, &mon, &day, &should_be_t, &hour, &min, &sec, timezone);
    if (num < 8)
        return 0;
    if (should_be_t != 'T' && should_be_t != 't' && should_be_t != ' ')
        return 0;
    struct tm parsed = { sec, min, hour, day, mon - 1, year - 1900, 0, 0, 0 };

    switch (timezone[0]) {
      case 'Z':
      case 'z':
      case '\0':
          break;
      case '+':
      case '-': {
          unsigned hr_off, min_off;
          if (sscanf(timezone + 1, "%2u:%2u", &hr_off, &min_off)) {
              if (hr_off < 24 && min_off <= 60) {
                  offset = ((hr_off * 60 + min_off)
                            * (timezone[0] == '+' ? -60 : 60));
              }
          }
          break;
      }
      default:
          return 0;
    }

    when = hts_time_gm(&parsed);
    return when >= 0 ? when + offset : 0;
}

static void refresh_auth_data(s3_auth_data *ad) {
    // Basically a copy of the AWS_SHARED_CREDENTIALS_FILE part of
    // setup_auth_data(), but this only reads the authorisation parts.
    const char *v = getenv("AWS_SHARED_CREDENTIALS_FILE");
    kstring_t expiry_time = KS_INITIALIZE;
    parse_ini(v? v : "~/.aws/credentials", ad->profile.s,
              "aws_access_key_id", &ad->id,
              "aws_secret_access_key", &ad->secret,
              "aws_session_token", &ad->token,
              "expiry_time", &expiry_time);
    if (expiry_time.l) {
        ad->creds_expiry_time = parse_rfc3339_date(&expiry_time);
    }
    ks_free(&expiry_time);
}

static int auth_header_callback(void *ctx, char ***hdrs) {
    s3_auth_data *ad = (s3_auth_data *) ctx;

    time_t now = time(NULL);
#ifdef HAVE_GMTIME_R
    struct tm tm_buffer;
    struct tm *tm = gmtime_r(&now, &tm_buffer);
#else
    struct tm *tm = gmtime(&now);
#endif
    kstring_t message = { 0, 0, NULL };
    unsigned char digest[DIGEST_BUFSIZ];
    size_t digest_len;

    if (!hdrs) { // Closing connection
        free_auth_data(ad);
        return 0;
    }

    if (ad->creds_expiry_time > 0
        && ad->creds_expiry_time - now < CREDENTIAL_LIFETIME) {
        refresh_auth_data(ad);
    } else if (now - ad->auth_time < AUTH_LIFETIME) {
        // Last auth string should still be valid
        *hdrs = NULL;
        return 0;
    }

    strftime(ad->date, sizeof(ad->date), "Date: %a, %d %b %Y %H:%M:%S GMT", tm);
    if (!ad->id.l || !ad->secret.l) {
        ad->auth_time = now;
        return copy_auth_headers(ad, hdrs);
    }

    if (ksprintf(&message, "%s\n\n\n%s\n%s%s%s%s",
                 ad->mode == 'r' ? "GET" : "PUT", ad->date + 6,
                 ad->token.l ? "x-amz-security-token:" : "",
                 ad->token.l ? ad->token.s : "",
                 ad->token.l ? "\n" : "",
                 ad->bucket) < 0) {
        return -1;
    }

    digest_len = s3_sign(digest, &ad->secret, &message);
    ad->auth_hdr.l = 0;
    if (ksprintf(&ad->auth_hdr, "Authorization: AWS %s:", ad->id.s) < 0)
        goto fail;
    base64_kput(digest, digest_len, &ad->auth_hdr);

    free(message.s);
    ad->auth_time = now;
    return copy_auth_headers(ad, hdrs);

 fail:
    free(message.s);
    return -1;
}


/* like a escape path but for query strings '=' and '&' are untouched */
static char *escape_query(const char *qs) {
    size_t i, j = 0, length;
    char *escaped;

    length = strlen(qs);

    if ((escaped = malloc(length * 3 + 1)) == NULL) {
        return NULL;
    }

    for (i = 0; i < length; i++) {
        int c = qs[i];

        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
             c == '_' || c == '-' || c == '~' || c == '.' || c == '/' || c == '=' || c == '&') {
            escaped[j++] = c;
        } else {
            sprintf(escaped + j, "%%%02X", c);
            j += 3;
        }
    }

    if (i != length) {
        // in the case of a '?' copy the rest of the qs across unchanged
        strcpy(escaped + j, qs + i);
    } else {
        escaped[j] = '\0';
    }

    return escaped;
}


static char *escape_path(const char *path) {
    size_t i, j = 0, length;
    char *escaped;

    length = strlen(path);

    if ((escaped = malloc(length * 3 + 1)) == NULL) {
        return NULL;
    }

    for (i = 0; i < length; i++) {
        int c = path[i];

        if (c == '?') break; // don't escape ? or beyond

        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
             c == '_' || c == '-' || c == '~' || c == '.' || c == '/') {
            escaped[j++] = c;
        } else {
            sprintf(escaped + j, "%%%02X", c);
            j += 3;
        }
    }

    if (i != length) {
        // in the case of a '?' copy the rest of the path across unchanged
        strcpy(escaped + j, path + i);
    } else {
        escaped[j] = '\0';
    }

    return escaped;
}


static int is_escaped(const char *str) {
    const char *c = str;
    int escaped = 0;
    int needs_escape = 0;

    while (*c != '\0') {
        if (*c == '%' && c[1] != '\0' && c[2] != '\0') {
            if (isxdigit_c(c[1]) && isxdigit_c(c[2])) {
                escaped = 1;
                c += 3;
                continue;
            } else {
                // only escaped if all % signs are escaped
                escaped = 0;
            }
        }
        if (!((*c >= '0' && *c <= '9') || (*c >= 'A' && *c <= 'Z')
              || (*c >= 'a' && *c <= 'z') ||
              *c == '_' || *c == '-' || *c == '~' || *c == '.' || *c == '/')) {
            needs_escape = 1;
        }
        c++;
    }

    return escaped || !needs_escape;
}

static int redirect_endpoint_callback(void *auth, long response,
                                      kstring_t *header, kstring_t *url) {
    s3_auth_data *ad = (s3_auth_data *)auth;
    char *new_region;
    char *end;
    int ret = -1;

    // get the new region from the reply header
    if ((new_region = strstr(header->s, "x-amz-bucket-region: "))) {

        new_region += strlen("x-amz-bucket-region: ");
        end = new_region;

        while (isalnum_c(*end) || ispunct_c(*end)) end++;

        *end = 0;

        if (strstr(ad->host.s, "amazonaws.com")) {
            ad->region.l = 0;
            kputs(new_region, &ad->region);

            ad->host.l = 0;
            ksprintf(&ad->host, "s3.%s.amazonaws.com", new_region);

            if (ad->region.l && ad->host.l) {
               url->l = 0;
               kputs(ad->host.s, url);
               kputsn(ad->bucket, strlen(ad->bucket), url);
               if (ad->user_query_string.l) {
                   kputc('?', url);
                   kputsn(ad->user_query_string.s, ad->user_query_string.l, url);
               }
               ret = 0;
            }
        }
    }

    return ret;
}

static s3_auth_data * setup_auth_data(const char *s3url, const char *mode,
                                      int sigver, kstring_t *url)
{
    s3_auth_data *ad = calloc(1, sizeof(*ad));
    const char *bucket, *path;
    char *escaped = NULL;
    size_t url_path_pos;
    ptrdiff_t bucket_len;
    int is_https = 1, dns_compliant;
    char *query_start;
    enum {s3_auto, s3_virtual, s3_path} address_style = s3_auto;

    if (!ad)
        return NULL;
    ad->mode = strchr(mode, 'r') ? 'r' : 'w';

    // Our S3 URL format is s3[+SCHEME]://[ID[:SECRET[:TOKEN]]@]BUCKET/PATH

    if (s3url[2] == '+') {
        bucket = strchr(s3url, ':') + 1;
        if (bucket == NULL) {
            free(ad);
            return NULL;
        }
        kputsn(&s3url[3], bucket - &s3url[3], url);
        is_https = strncmp(url->s, "https:", 6) == 0;
    }
    else {
        kputs("https:", url);
        bucket = &s3url[3];
    }
    while (*bucket == '/') kputc(*bucket++, url);

    path = bucket + strcspn(bucket, "/?#@");

    if (*path == '@') {
        const char *colon = strpbrk(bucket, ":@");
        if (*colon != ':') {
            urldecode_kput(bucket, colon - bucket, &ad->profile);
        }
        else {
            const char *colon2 = strpbrk(&colon[1], ":@");
            urldecode_kput(bucket, colon - bucket, &ad->id);
            urldecode_kput(&colon[1], colon2 - &colon[1], &ad->secret);
            if (*colon2 == ':')
                urldecode_kput(&colon2[1], path - &colon2[1], &ad->token);
        }

        bucket = &path[1];
        path = bucket + strcspn(bucket, "/?#");
    }
    else {
        // If the URL has no ID[:SECRET]@, consider environment variables.
        const char *v;
        if ((v = getenv("AWS_ACCESS_KEY_ID")) != NULL) kputs(v, &ad->id);
        if ((v = getenv("AWS_SECRET_ACCESS_KEY")) != NULL) kputs(v, &ad->secret);
        if ((v = getenv("AWS_SESSION_TOKEN")) != NULL) kputs(v, &ad->token);
        if ((v = getenv("AWS_DEFAULT_REGION")) != NULL) kputs(v, &ad->region);
        if ((v = getenv("HTS_S3_HOST")) != NULL) kputs(v, &ad->host);

        if ((v = getenv("AWS_DEFAULT_PROFILE")) != NULL) kputs(v, &ad->profile);
        else if ((v = getenv("AWS_PROFILE")) != NULL) kputs(v, &ad->profile);
        else kputs("default", &ad->profile);

        if ((v = getenv("HTS_S3_ADDRESS_STYLE")) != NULL) {
            if (strcasecmp(v, "virtual") == 0) {
                address_style = s3_virtual;
            } else if (strcasecmp(v, "path") == 0) {
                address_style = s3_path;
            }
        }
    }

    if (ad->id.l == 0) {
        kstring_t url_style = KS_INITIALIZE;
        kstring_t expiry_time = KS_INITIALIZE;
        const char *v = getenv("AWS_SHARED_CREDENTIALS_FILE");
        parse_ini(v? v : "~/.aws/credentials", ad->profile.s,
                  "aws_access_key_id", &ad->id,
                  "aws_secret_access_key", &ad->secret,
                  "aws_session_token", &ad->token,
                  "region", &ad->region,
                  "addressing_style", &url_style,
                  "expiry_time", &expiry_time,
                  NULL);

        if (url_style.l) {
            if (strcmp(url_style.s, "virtual") == 0) {
                address_style = s3_virtual;
            } else if (strcmp(url_style.s, "path") == 0) {
                address_style = s3_path;
            } else {
                address_style = s3_auto;
            }
        }
        if (expiry_time.l) {
            // Not a real part of the AWS configuration file, but it allows
            // support for short-term credentials like those for the IAM
            // service.  The botocore library uses the key "expiry_time"
            // internally for this purpose.
            // See https://github.com/boto/botocore/blob/develop/botocore/credentials.py
            ad->creds_expiry_time = parse_rfc3339_date(&expiry_time);
        }

        ks_free(&url_style);
        ks_free(&expiry_time);
    }

    if (ad->id.l == 0) {
        kstring_t url_style = KS_INITIALIZE;
        const char *v = getenv("HTS_S3_S3CFG");
        parse_ini(v? v : "~/.s3cfg", ad->profile.s, "access_key", &ad->id,
                  "secret_key", &ad->secret, "access_token", &ad->token,
                  "host_base", &ad->host,
                  "bucket_location", &ad->region,
                  "host_bucket", &url_style,
                  NULL);

        if (url_style.l) {
            // Conforming to s3cmd's GitHub PR#416, host_bucket without the "%(bucket)s" string
            // indicates use of path style adressing.
            if (strstr(url_style.s, "%(bucket)s") == NULL) {
                address_style = s3_path;
            } else {
                address_style = s3_auto;
            }
        }

        ks_free(&url_style);
    }

    if (ad->id.l == 0)
        parse_simple("~/.awssecret", &ad->id, &ad->secret);


    // if address_style is set, force the dns_compliant setting
    if (address_style == s3_virtual) {
        dns_compliant = 1;
    } else if (address_style == s3_path) {
        dns_compliant = 0;
    } else {
        dns_compliant = is_dns_compliant(bucket, path, is_https);
    }

    if (ad->host.l == 0)
        kputs("s3.amazonaws.com", &ad->host);

    if (!dns_compliant && ad->region.l > 0
        && strcmp(ad->host.s, "s3.amazonaws.com") == 0) {
        // Can avoid a redirection by including the region in the host name
        // (assuming the right one has been specified)
        ad->host.l = 0;
        ksprintf(&ad->host, "s3.%s.amazonaws.com", ad->region.s);
    }

    if (ad->region.l == 0)
        kputs("us-east-1", &ad->region);

    if (!is_escaped(path)) {
        escaped = escape_path(path);
        if (escaped == NULL) {
            goto error;
        }
    }

    bucket_len = path - bucket;

    // Use virtual hosted-style access if possible, otherwise path-style.
    if (dns_compliant) {
        size_t url_host_pos = url->l;
        // Append "bucket.host" to url
        kputsn_(bucket, bucket_len, url);
        kputc('.', url);
        kputsn(ad->host.s, ad->host.l, url);
        url_path_pos = url->l;

        if (sigver == 4) {
            // Copy back to ad->host to use when making the signature
            ad->host.l = 0;
            kputsn(url->s + url_host_pos, url->l - url_host_pos, &ad->host);
        }
    }
    else {
        // Append "host/bucket" to url
        kputsn(ad->host.s, ad->host.l, url);
        url_path_pos = url->l;
        kputc('/', url);
        kputsn(bucket, bucket_len, url);
    }

    kputs(escaped == NULL ? path : escaped, url);

    if (sigver == 4 || !dns_compliant) {
        ad->bucket = malloc(url->l - url_path_pos + 1);
        if (ad->bucket == NULL) {
            goto error;
        }
        memcpy(ad->bucket, url->s + url_path_pos, url->l - url_path_pos + 1);
    }
    else {
        ad->bucket = malloc(url->l - url_path_pos + bucket_len + 2);
        if (ad->bucket == NULL) {
            goto error;
        }
        ad->bucket[0] = '/';
        memcpy(ad->bucket + 1, bucket, bucket_len);
        memcpy(ad->bucket + bucket_len + 1,
               url->s + url_path_pos, url->l - url_path_pos + 1);
    }

    // write any query strings to its own place to use later
    if ((query_start = strchr(ad->bucket, '?'))) {
        kputs(query_start + 1, &ad->user_query_string);
        *query_start = 0;
    }

    free(escaped);

    return ad;

 error:
    free(escaped);
    free_auth_data(ad);
    return NULL;
}

static hFILE * s3_rewrite(const char *s3url, const char *mode, va_list *argsp)
{
    kstring_t url = { 0, 0, NULL };
    s3_auth_data *ad = setup_auth_data(s3url, mode, 2, &url);

    if (!ad)
        return NULL;

    hFILE *fp = hopen(url.s, mode, "va_list", argsp,
                      "httphdr_callback", auth_header_callback,
                      "httphdr_callback_data", ad,
                      "redirect_callback", redirect_endpoint_callback,
                      "redirect_callback_data", ad,
                      NULL);
    if (!fp) goto fail;

    free(url.s);
    return fp;

 fail:
    free(url.s);
    free_auth_data(ad);
    return NULL;
}

/***************************************************************

AWS S3 sig version 4 writing code

****************************************************************/

static void hash_string(char *in, size_t length, char *out) {
    unsigned char hashed[SHA256_DIGEST_BUFSIZE];
    int i, j;

    s3_sha256((const unsigned char *)in, length, hashed);

    for (i = 0, j = 0; i < SHA256_DIGEST_BUFSIZE; i++, j+= 2) {
        sprintf(out + j, "%02x", hashed[i]);
    }
}

static void ksinit(kstring_t *s) {
    s->l = 0;
    s->m = 0;
    s->s = NULL;
}


static void ksfree(kstring_t *s) {
    free(s->s);
    ksinit(s);
}


static int make_signature(s3_auth_data *ad, kstring_t *string_to_sign, char *signature_string) {
    unsigned char date_key[SHA256_DIGEST_BUFSIZE];
    unsigned char date_region_key[SHA256_DIGEST_BUFSIZE];
    unsigned char date_region_service_key[SHA256_DIGEST_BUFSIZE];
    unsigned char signing_key[SHA256_DIGEST_BUFSIZE];
    unsigned char signature[SHA256_DIGEST_BUFSIZE];

    const unsigned char service[] = "s3";
    const unsigned char request[] = "aws4_request";

    kstring_t secret_access_key = {0, 0, NULL};
    unsigned int len;
    unsigned int i, j;

    ksprintf(&secret_access_key, "AWS4%s", ad->secret.s);

    if (secret_access_key.l == 0) {
        return -1;
    }

    s3_sign_sha256(secret_access_key.s, secret_access_key.l, (const unsigned char *)ad->date_short, strlen(ad->date_short), date_key, &len);
    s3_sign_sha256(date_key, len, (const unsigned char *)ad->region.s, ad->region.l, date_region_key, &len);
    s3_sign_sha256(date_region_key, len, service, 2, date_region_service_key, &len);
    s3_sign_sha256(date_region_service_key, len, request, 12, signing_key, &len);
    s3_sign_sha256(signing_key, len, (const unsigned char *)string_to_sign->s, string_to_sign->l, signature, &len);

    for (i = 0, j = 0; i < len; i++, j+= 2) {
        sprintf(signature_string + j, "%02x", signature[i]);
    }

    ksfree(&secret_access_key);

    return 0;
}


static int make_authorisation(s3_auth_data *ad, char *http_request, char *content, kstring_t *auth) {
    kstring_t signed_headers = {0, 0, NULL};
    kstring_t canonical_headers = {0, 0, NULL};
    kstring_t canonical_request = {0, 0, NULL};
    kstring_t scope = {0, 0, NULL};
    kstring_t string_to_sign = {0, 0, NULL};
    char cr_hash[HASH_LENGTH_SHA256];
    char signature_string[HASH_LENGTH_SHA256];
    int ret = -1;


    if (!ad->token.l) {
        kputs("host;x-amz-content-sha256;x-amz-date", &signed_headers);
    } else {
        kputs("host;x-amz-content-sha256;x-amz-date;x-amz-security-token", &signed_headers);
    }

    if (signed_headers.l == 0) {
        return -1;
    }


    if (!ad->token.l) {
        ksprintf(&canonical_headers, "host:%s\nx-amz-content-sha256:%s\nx-amz-date:%s\n",
        ad->host.s, content, ad->date_long);
    } else {
        ksprintf(&canonical_headers, "host:%s\nx-amz-content-sha256:%s\nx-amz-date:%s\nx-amz-security-token:%s\n",
        ad->host.s, content, ad->date_long, ad->token.s);
    }

    if (canonical_headers.l == 0) {
        goto cleanup;
    }

    // bucket == canonical_uri
    ksprintf(&canonical_request, "%s\n%s\n%s\n%s\n%s\n%s",
        http_request, ad->bucket, ad->canonical_query_string.s,
        canonical_headers.s, signed_headers.s, content);

    if (canonical_request.l == 0) {
        goto cleanup;
    }

    hash_string(canonical_request.s, canonical_request.l, cr_hash);

    ksprintf(&scope, "%s/%s/s3/aws4_request", ad->date_short, ad->region.s);

    if (scope.l == 0) {
        goto cleanup;
    }

    ksprintf(&string_to_sign, "AWS4-HMAC-SHA256\n%s\n%s\n%s", ad->date_long, scope.s, cr_hash);

    if (string_to_sign.l == 0) {
        goto cleanup;
    }

    if (make_signature(ad, &string_to_sign, signature_string)) {
        goto cleanup;
    }

    ksprintf(auth, "Authorization: AWS4-HMAC-SHA256 Credential=%s/%s/%s/s3/aws4_request,SignedHeaders=%s,Signature=%s",
                ad->id.s, ad->date_short, ad->region.s, signed_headers.s, signature_string);

    if (auth->l == 0) {
        goto cleanup;
    }

    ret = 0;

 cleanup:
    ksfree(&signed_headers);
    ksfree(&canonical_headers);
    ksfree(&canonical_request);
    ksfree(&scope);
    ksfree(&string_to_sign);

    return ret;
}


static int update_time(s3_auth_data *ad, time_t now) {
    int ret = -1;
#ifdef HAVE_GMTIME_R
    struct tm tm_buffer;
    struct tm *tm = gmtime_r(&now, &tm_buffer);
#else
    struct tm *tm = gmtime(&now);
#endif

    if (now - ad->auth_time > AUTH_LIFETIME) {
        // update timestamp
        ad->auth_time = now;

        if (strftime(ad->date_long, 17, "%Y%m%dT%H%M%SZ", tm) != 16) {
            return -1;
        }

        if (strftime(ad->date_short, 9, "%Y%m%d", tm) != 8) {
            return -1;;
        }

        ad->date_html.l = 0;
        ksprintf(&ad->date_html, "x-amz-date: %s", ad->date_long);
    }

    if (ad->date_html.l) ret = 0;

    return ret;
}


static int query_cmp(const void *p1, const void *p2) {
    char **q1 = (char **)p1;
    char **q2 = (char **)p2;

    return strcmp(*q1, *q2);
}


/* Query strings must be in alphabetical order for authorisation */

static int order_query_string(kstring_t *qs) {
    int *query_offset = NULL;
    int num_queries, i;
    char **queries = NULL;
    kstring_t ordered = {0, 0, NULL};
    char *escaped = NULL;
    int ret = -1;

    if ((query_offset = ksplit(qs, '&', &num_queries)) == NULL) {
        return -1;
    }

    if ((queries = malloc(num_queries * sizeof(char*))) == NULL)
        goto err;

    for (i = 0; i < num_queries; i++) {
        queries[i] = qs->s + query_offset[i];
    }

    qsort(queries, num_queries, sizeof(char *), query_cmp);

    for (i = 0; i < num_queries; i++) {
        if (i) {
            kputs("&", &ordered);
        }

        kputs(queries[i], &ordered);
    }

    if ((escaped = escape_query(ordered.s)) == NULL)
        goto err;

    qs->l = 0;
    kputs(escaped, qs);

    ret = 0;
 err:
    free(ordered.s);
    free(queries);
    free(query_offset);
    free(escaped);

    return ret;
}


static int write_authorisation_callback(void *auth, char *request, kstring_t *content, char *cqs,
                                        kstring_t *hash, kstring_t *auth_str, kstring_t *date,
                                        kstring_t *token, int uqs) {
    s3_auth_data *ad = (s3_auth_data *)auth;
    char content_hash[HASH_LENGTH_SHA256];
    time_t now;

    if (request == NULL) {
        // signal to free auth data
        free_auth_data(ad);
        return 0;
    }

    now = time(NULL);

    if (update_time(ad, now)) {
        return -1;
    }
    if (ad->creds_expiry_time > 0
        && ad->creds_expiry_time - now < CREDENTIAL_LIFETIME) {
        refresh_auth_data(ad);
    }

    if (content) {
        hash_string(content->s, content->l, content_hash);
    } else {
        // empty hash
        hash_string("", 0, content_hash);
    }

    ad->canonical_query_string.l = 0;
    kputs(cqs, &ad->canonical_query_string);

    if (ad->canonical_query_string.l == 0) {
        return -1;
    }

    /* add a user provided query string, normally only useful on upload initiation */
    if (uqs) {
        kputs("&", &ad->canonical_query_string);
        kputs(ad->user_query_string.s, &ad->canonical_query_string);

        if (order_query_string(&ad->canonical_query_string)) {
            return -1;
        }
    }

    if (make_authorisation(ad, request, content_hash, auth_str)) {
        return -1;
    }

    kputs(ad->date_html.s, date);
    kputsn(content_hash, HASH_LENGTH_SHA256, hash);

    if (date->l == 0 || hash->l == 0) {
        return -1;
    }

    if (ad->token.l) {
        ksprintf(token, "x-amz-security-token: %s", ad->token.s);
    }

    return 0;
}


static int v4_auth_header_callback(void *ctx, char ***hdrs) {
    s3_auth_data *ad = (s3_auth_data *) ctx;
    char content_hash[HASH_LENGTH_SHA256];
    kstring_t content = KS_INITIALIZE;
    kstring_t authorisation = KS_INITIALIZE;
    kstring_t token_hdr = KS_INITIALIZE;
    char *date_html = NULL;
    time_t now;
    int idx;

    if (!hdrs) { // Closing connection
        free_auth_data(ad);
        return 0;
    }

    now = time(NULL);

    if (update_time(ad, now)) {
        return -1;
    }

    if (ad->creds_expiry_time > 0
        && ad->creds_expiry_time - now < CREDENTIAL_LIFETIME) {
        refresh_auth_data(ad);
    }

    if (!ad->id.l || !ad->secret.l) {
        return copy_auth_headers(ad, hdrs);
    }

    hash_string("", 0, content_hash); // empty hash

    ad->canonical_query_string.l = 0;

    if (ad->user_query_string.l > 0) {
        kputs(ad->user_query_string.s, &ad->canonical_query_string);

        if (order_query_string(&ad->canonical_query_string)) {
            return -1;
        }
    } else {
        kputs("", &ad->canonical_query_string);
    }

    if (make_authorisation(ad, "GET", content_hash, &authorisation)) {
        return -1;
    }

    ksprintf(&content, "x-amz-content-sha256: %s", content_hash);
    date_html = strdup(ad->date_html.s);

    if (ad->token.l > 0) {
        kputs("X-Amz-Security-Token: ", &token_hdr);
        kputs(ad->token.s, &token_hdr);
    }

    if (content.l == 0 || date_html == NULL) {
        ksfree(&authorisation);
        ksfree(&content);
        ksfree(&token_hdr);
        free(date_html);
        return -1;
    }

    *hdrs = &ad->headers[0];
    idx = 0;
    ad->headers[idx++] = ks_release(&authorisation);
    ad->headers[idx++] = date_html;
    ad->headers[idx++] = ks_release(&content);
    if (token_hdr.s)
        ad->headers[idx++] = ks_release(&token_hdr);
    ad->headers[idx++] = NULL;

    return 0;
}

static int handle_400_response(hFILE *fp, s3_auth_data *ad) {
    // v4 signatures in virtual hosted mode return 400 Bad Request if the
    // wrong region is used to make the signature.  The response is an xml
    // document which includes the name of the correct region.  This can
    // be extracted and used to generate a corrected signature.
    // As the xml is fairly simple, go with something "good enough" instead
    // of trying to parse it properly.

    char buffer[1024], *region, *reg_end;
    ssize_t bytes;

    bytes = hread(fp, buffer, sizeof(buffer) - 1);
    if (bytes < 0) {
        return -1;
    }
    buffer[bytes] = '\0';
    region = strstr(buffer, "<Region>");
    if (region == NULL) {
        return -1;
    }
    region += 8;
    while (isspace((unsigned char) *region)) ++region;
    reg_end = strchr(region, '<');
    if (reg_end == NULL || strncmp(reg_end + 1, "/Region>", 8) != 0) {
        return -1;
    }
    while (reg_end > region && isspace((unsigned char) reg_end[-1])) --reg_end;
    ad->region.l = 0;
    kputsn(region, reg_end - region, &ad->region);
    if (ad->region.l == 0) {
        return -1;
    }

    return 0;
}

static int set_region(void *adv, kstring_t *region) {
    s3_auth_data *ad = (s3_auth_data *) adv;

    ad->region.l = 0;
    return kputsn(region->s, region->l, &ad->region) < 0;
}

static int http_status_errno(int status)
{
    if (status >= 500)
        switch (status) {
        case 501: return ENOSYS;
        case 503: return EBUSY;
        case 504: return ETIMEDOUT;
        default:  return EIO;
        }
    else if (status >= 400)
        switch (status) {
        case 401: return EPERM;
        case 403: return EACCES;
        case 404: return ENOENT;
        case 405: return EROFS;
        case 407: return EPERM;
        case 408: return ETIMEDOUT;
        case 410: return ENOENT;
        default:  return EINVAL;
        }
    else return 0;
}

static hFILE *s3_open_v4(const char *s3url, const char *mode, va_list *argsp) {
    kstring_t url = { 0, 0, NULL };

    s3_auth_data *ad = setup_auth_data(s3url, mode, 4, &url);
    hFILE *fp = NULL;

    if (ad == NULL) {
        return NULL;
    }

    if (ad->mode == 'r') {
        long http_response = 0;

        fp = hopen(url.s, mode, "va_list", argsp,
                   "httphdr_callback", v4_auth_header_callback,
                   "httphdr_callback_data", ad,
                   "redirect_callback", redirect_endpoint_callback,
                   "redirect_callback_data", ad,
                   "http_response_ptr", &http_response,
                   "fail_on_error", 0,
                   NULL);

        if (fp == NULL) goto error;

        if (http_response == 400) {
            ad->refcount = 1;
            if (handle_400_response(fp, ad) != 0) {
                goto error;
            }
            hclose_abruptly(fp);
            fp = hopen(url.s, mode, "va_list", argsp,
                       "httphdr_callback", v4_auth_header_callback,
                       "httphdr_callback_data", ad,
                       "redirect_callback", redirect_endpoint_callback,
                       "redirect_callback_data", ad,
                       NULL);
        } else if (http_response > 400) {
            ad->refcount = 1;
            errno = http_status_errno(http_response);
            goto error;
        }

        if (fp == NULL) goto error;
    } else {
        kstring_t final_url = {0, 0, NULL};

         // add the scheme marker
        ksprintf(&final_url, "s3w+%s", url.s);

        if(final_url.l == 0) goto error;

        fp = hopen(final_url.s, mode, "va_list", argsp,
                   "s3_auth_callback",  write_authorisation_callback,
                   "s3_auth_callback_data", ad,
                   "redirect_callback", redirect_endpoint_callback,
                   "set_region_callback", set_region,
                   NULL);
        free(final_url.s);

        if (fp == NULL) goto error;
    }

    free(url.s);

    return fp;

  error:

    if (fp) hclose_abruptly(fp);
    free(url.s);
    free_auth_data(ad);

    return NULL;
}


static hFILE *s3_open(const char *url, const char *mode)
{
    hFILE *fp;

    kstring_t mode_colon = { 0, 0, NULL };
    kputs(mode, &mode_colon);
    kputc(':', &mode_colon);

    if (getenv("HTS_S3_V2") == NULL) { // Force the v2 signature code
        fp = s3_open_v4(url, mode_colon.s, NULL);
    } else {
        fp = s3_rewrite(url, mode_colon.s, NULL);
    }

    free(mode_colon.s);

    return fp;
}

static hFILE *s3_vopen(const char *url, const char *mode_colon, va_list args0)
{
    hFILE *fp;
    // Need to use va_copy() as we can only take the address of an actual
    // va_list object, not that of a parameter whose type may have decayed.
    va_list args;
    va_copy(args, args0);

    if (getenv("HTS_S3_V2") == NULL) { // Force the v2 signature code
        fp = s3_open_v4(url, mode_colon, &args);
    } else {
        fp = s3_rewrite(url, mode_colon, &args);
    }

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
    static const char id[] = "@(#)hfile_s3 plugin (htslib)\t" HTS_VERSION_TEXT;
    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_s3.init] version %s\n", strchr(id, '\t')+1);
#endif

    self->name = "Amazon S3";
    hfile_add_scheme_handler("s3", &handler);
    hfile_add_scheme_handler("s3+http", &handler);
    hfile_add_scheme_handler("s3+https", &handler);
    return 0;
}
