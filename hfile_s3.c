/*  hfile_s3.c -- Amazon S3 backend for low-level file streams.

    Copyright (C) 2015-2017, 2019-2025 Genome Research Ltd.

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
#include <pthread.h>

#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"
#include "hts_time_funcs.h"

#include <curl/curl.h>

typedef struct s3_auth_data {
    kstring_t id;
    kstring_t token;
    kstring_t secret;
    kstring_t region;
    kstring_t canonical_query_string;
    kstring_t user_query_string;
    kstring_t host;
    kstring_t profile;
    enum {s3_auto, s3_virtual, s3_path} url_style;
    time_t creds_expiry_time;
    char *bucket;
    time_t auth_time;
    char date[40];
    char date_long[17];
    char date_short[9];
    kstring_t date_html;
    char mode;
    int is_v4;
} s3_auth_data;

typedef struct {
    hFILE base;
    CURL *curl;
    CURLcode ret;
    s3_auth_data *au;
    kstring_t buffer;
    kstring_t url;
    long verbose;
    int write;
    int part_size; // size for reading or writing

    kstring_t content_hash;
    kstring_t authorisation;
    kstring_t content;
    kstring_t date;
    kstring_t token;
    kstring_t range;

    // write variables
    kstring_t upload_id;
    kstring_t completion_message;
    int part_no;
    int aborted;
    size_t index;
    int expand;

    // read variables
    size_t last_read;               // last read position (remote)
    size_t last_read_buffer;        // last read (local buffer)
    int64_t file_size;              // size of the file being read
    int keep_going;

} hFILE_s3;

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


static void free_auth_data(s3_auth_data *ad) {
    free(ad->profile.s);
    free(ad->id.s);
    free(ad->token.s);
    free(ad->secret.s);
    free(ad->region.s);
    free(ad->canonical_query_string.s);
    free(ad->user_query_string.s);
    free(ad->host.s);
    free(ad->bucket);
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


/* like a escape path but for query strings '=' and '&' are untouched */
static char *escape_query(const char *qs) {
    size_t i, j = 0, length, alloced;
    char *escaped;

    length = strlen(qs);
    alloced = length * 3 + 1;
    if ((escaped = malloc(alloced)) == NULL) {
        return NULL;
    }

    for (i = 0; i < length; i++) {
        int c = qs[i];

        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
             c == '_' || c == '-' || c == '~' || c == '.' || c == '/' || c == '=' || c == '&') {
            escaped[j++] = c;
        } else {
            snprintf(escaped + j, alloced - j, "%%%02X", c);
            j += 3;
        }
    }

    escaped[j] = '\0';

    return escaped;
}


static char *escape_path(const char *path) {
    size_t i, j = 0, length, alloced;
    char *escaped;

    length = strlen(path);
    alloced = length * 3 + 1;

    if ((escaped = malloc(alloced)) == NULL) {
        return NULL;
    }

    for (i = 0; i < length; i++) {
        int c = path[i];

        if (c == '?') break; // don't escape ? or beyond

        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
             c == '_' || c == '-' || c == '~' || c == '.' || c == '/') {
            escaped[j++] = c;
        } else {
            snprintf(escaped + j, alloced - j, "%%%02X", c);
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


static int redirect_endpoint(hFILE_s3 *fp, kstring_t *header) {
    s3_auth_data *ad = fp->au;
    kstring_t *url = &fp->url;
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

            if (ad->url_style == s3_path) {
                // Path style https://s3.{region-code}.amazonaws.com/{bucket-name}/{key-name}
                ksprintf(&ad->host, "s3.%s.amazonaws.com", new_region);
            } else {
                // Virtual https://{bucket-name}.s3.{region-code}.amazonaws.com/{key-name}
                // Extract the {bucket-name} from {ad->host} to include in subdomain
                kstring_t url_prefix = KS_INITIALIZE;
                kputsn(ad->host.s, strcspn(ad->host.s, "."), &url_prefix);

                ksprintf(&ad->host, "%s.s3.%s.amazonaws.com", url_prefix.s, new_region);
                free(url_prefix.s);
            }
            if (ad->region.l && ad->host.l) {
               int e = 0;
               url->l = 0;
               e |= kputs("https://", url) < 0;
               e |= kputs(ad->host.s, url) < 0;
               e |= kputsn(ad->bucket, strlen(ad->bucket), url) < 0;

               if (!e)
                   ret = 0;
            }
            if (ad->user_query_string.l) {
                kputc('?', url);
                kputsn(ad->user_query_string.s, ad->user_query_string.l, url);
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

    if (!ad)
        return NULL;
    ad->mode = strchr(mode, 'r') ? 'r' : 'w';
    ad->url_style = s3_auto;

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
                ad->url_style = s3_virtual;
            } else if (strcasecmp(v, "path") == 0) {
                ad->url_style = s3_path;
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
                ad->url_style = s3_virtual;
            } else if (strcmp(url_style.s, "path") == 0) {
                ad->url_style = s3_path;
            } else {
                ad->url_style = s3_auto;
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
                ad->url_style = s3_path;
            } else {
                ad->url_style = s3_auto;
            }
        }

        ks_free(&url_style);
    }

    if (ad->id.l == 0)
        parse_simple("~/.awssecret", &ad->id, &ad->secret);


    // if address_style is set, force the dns_compliant setting
    if (ad->url_style == s3_virtual) {
        dns_compliant = 1;
    } else if (ad->url_style == s3_path) {
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
        ad->is_v4 = 1;
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
        ad->is_v4 = 0;
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


static int v2_authorisation(hFILE_s3 *fp, char *request) {
    s3_auth_data *ad = fp->au;
    time_t now = time(NULL);

#ifdef HAVE_GMTIME_R
    struct tm tm_buffer;
    struct tm *tm = gmtime_r(&now, &tm_buffer);
#else
    struct tm *tm = gmtime(&now);
#endif

    kstring_t message = KS_INITIALIZE;
    unsigned char digest[DIGEST_BUFSIZ];
    size_t digest_len;

    if (ad->creds_expiry_time > 0
        && ad->creds_expiry_time - now < CREDENTIAL_LIFETIME) {
        refresh_auth_data(ad);
    }

    // date format between v2 and v4 is different.

    strftime(ad->date, sizeof(ad->date), "Date: %a, %d %b %Y %H:%M:%S GMT", tm);

    kputs(ad->date, &fp->date);

    if (!ad->id.l || !ad->secret.l) {
        ad->auth_time = now;
        return 0;
    }

    if (ksprintf(&message, "%s\n\n\n%s\n%s%s%s%s",
                 request, ad->date + 6,
                 ad->token.l ? "x-amz-security-token:" : "",
                 ad->token.l ? ad->token.s : "",
                 ad->token.l ? "\n" : "",
                 ad->bucket) < 0) {
        return -1;
    }

    digest_len = s3_sign(digest, &ad->secret, &message);

    if (ksprintf(&fp->authorisation, "Authorization: AWS %s:", ad->id.s) < 0)
        goto fail;

    base64_kput(digest, digest_len, &fp->authorisation);

    free(message.s);
    ad->auth_time = now;
    return 0;

 fail:
    free(message.s);
    return -1;
}

/***************************************************************

AWS S3 sig version 4 writing code

****************************************************************/

static void hash_string(char *in, size_t length, char *out, size_t out_len) {
    unsigned char hashed[SHA256_DIGEST_BUFSIZE];
    int i, j;

    s3_sha256((const unsigned char *)in, length, hashed);

    for (i = 0, j = 0; i < SHA256_DIGEST_BUFSIZE; i++, j+= 2) {
        snprintf(out + j, out_len - j, "%02x", hashed[i]);
    }
}


static int make_signature(s3_auth_data *ad, kstring_t *string_to_sign, char *signature_string, size_t sig_string_len) {
    unsigned char date_key[SHA256_DIGEST_BUFSIZE];
    unsigned char date_region_key[SHA256_DIGEST_BUFSIZE];
    unsigned char date_region_service_key[SHA256_DIGEST_BUFSIZE];
    unsigned char signing_key[SHA256_DIGEST_BUFSIZE];
    unsigned char signature[SHA256_DIGEST_BUFSIZE];

    const unsigned char service[] = "s3";
    const unsigned char request[] = "aws4_request";

    kstring_t secret_access_key = KS_INITIALIZE;
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
        snprintf(signature_string + j, sig_string_len - j, "%02x", signature[i]);
    }

    ks_free(&secret_access_key);

    return 0;
}


static int make_authorisation(s3_auth_data *ad, char *http_request, char *content, kstring_t *auth) {
    kstring_t signed_headers = KS_INITIALIZE;
    kstring_t canonical_headers = KS_INITIALIZE;
    kstring_t canonical_request = KS_INITIALIZE;
    kstring_t scope = KS_INITIALIZE;
    kstring_t string_to_sign = KS_INITIALIZE;
    char cr_hash[HASH_LENGTH_SHA256];
    char signature_string[HASH_LENGTH_SHA256];
    int ret = -1;

    if (!ad->id.l || !ad->secret.l) {
        return 0;
    }

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

    hash_string(canonical_request.s, canonical_request.l, cr_hash, sizeof(cr_hash));

    ksprintf(&scope, "%s/%s/s3/aws4_request", ad->date_short, ad->region.s);

    if (scope.l == 0) {
        goto cleanup;
    }

    ksprintf(&string_to_sign, "AWS4-HMAC-SHA256\n%s\n%s\n%s", ad->date_long, scope.s, cr_hash);

    if (string_to_sign.l == 0) {
        goto cleanup;
    }

    if (make_signature(ad, &string_to_sign, signature_string, sizeof(signature_string))) {
        goto cleanup;
    }

    ksprintf(auth, "Authorization: AWS4-HMAC-SHA256 Credential=%s/%s/%s/s3/aws4_request,SignedHeaders=%s,Signature=%s",
                ad->id.s, ad->date_short, ad->region.s, signed_headers.s, signature_string);

    if (auth->l == 0) {
        goto cleanup;
    }

    ret = 0;

 cleanup:
    ks_free(&signed_headers);
    ks_free(&canonical_headers);
    ks_free(&canonical_request);
    ks_free(&scope);
    ks_free(&string_to_sign);

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
    kstring_t ordered = KS_INITIALIZE;
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


static int v4_authorisation(hFILE_s3 *fp, char *request, kstring_t *content, char *cqs, int uqs) {
    s3_auth_data *ad = fp->au;
    char content_hash[HASH_LENGTH_SHA256];
    time_t now;

    now = time(NULL);

    if (update_time(ad, now)) {
        return -1;
    }

    if (ad->creds_expiry_time > 0
        && ad->creds_expiry_time - now < CREDENTIAL_LIFETIME) {
        refresh_auth_data(ad);
    }

    if (content) {
        hash_string(content->s, content->l, content_hash, sizeof(content_hash));
    } else {
        // empty hash
        hash_string("", 0, content_hash, sizeof(content_hash));
    }

    ad->canonical_query_string.l = 0;

    if (cqs) {
        kputs(cqs, &ad->canonical_query_string);

        /* add a user provided query string, normally only useful on upload initiation */
        if (uqs) {
            kputs("&", &ad->canonical_query_string);
            kputs(ad->user_query_string.s, &ad->canonical_query_string);

            if (order_query_string(&ad->canonical_query_string)) {
                return -1;
            }
        }
    }

    if (make_authorisation(ad, request, content_hash, &fp->authorisation)) {
        return -1;
    }

    kputs(ad->date_html.s, &fp->date);
    kputsn(content_hash, HASH_LENGTH_SHA256, &fp->content_hash);

    if (fp->date.l == 0 || fp->content_hash.l == 0) {
        return -1;
    }

    if (ad->token.l) {
        ksprintf(&fp->token, "x-amz-security-token: %s", ad->token.s);
    }

    return 0;
}

static int set_region(s3_auth_data *ad, kstring_t *region) {
    ad->region.l = 0;
    return kputsn(region->s, region->l, &ad->region) < 0;
}

//
// Writing and reading handling
//

// Some common code

#define S3_MOVED_PERMANENTLY 301
#define S3_BAD_REQUEST 400


static struct {
    kstring_t useragent;
    CURLSH *share;
    pthread_mutex_t share_lock;
} curl = { { 0, 0, NULL }, NULL, PTHREAD_MUTEX_INITIALIZER };

static void share_lock(CURL *handle, curl_lock_data data,
                       curl_lock_access access, void *userptr) {
    pthread_mutex_lock(&curl.share_lock);
}

static void share_unlock(CURL *handle, curl_lock_data data, void *userptr) {
    pthread_mutex_unlock(&curl.share_lock);
}


static void initialise_authorisation_values(hFILE_s3 *fp) {
    ks_initialize(&fp->content_hash);
    ks_initialize(&fp->authorisation);
    ks_initialize(&fp->content);
    ks_initialize(&fp->date);
    ks_initialize(&fp->token);
    ks_initialize(&fp->range);
}


static void clear_authorisation_values(hFILE_s3 *fp) {
    ks_clear(&fp->content_hash);
    ks_clear(&fp->authorisation);
    ks_clear(&fp->content);
    ks_clear(&fp->date);
    ks_clear(&fp->token);
    ks_clear(&fp->range);
}


static void free_authorisation_values(hFILE_s3 *fp) {
    ks_free(&fp->content_hash);
    ks_free(&fp->authorisation);
    ks_free(&fp->content);
    ks_free(&fp->date);
    ks_free(&fp->token);
    ks_free(&fp->range);
}

/* As the response text is case insensitive we need a version of strstr that
   is also case insensitive.  The response is small so no need to get too
   complicated on the string search.
*/
static char *stristr(char *haystack, char *needle) {

    while (*haystack) {
        char *h = haystack;
        char *n = needle;

        while (toupper(*h) == toupper(*n)) {
            h++, n++;
            if (!*h || !*n) break;
        }

        if (!*n) break;

        haystack++;
    }

    if (!*haystack) return NULL;

    return haystack;
}


static int get_entry(char *in, char *start_tag, char *end_tag, kstring_t *out) {
    char *start;
    char *end;

    if (!in) {
        return EOF;
    }

    start = stristr(in, start_tag);
    if (!start) return EOF;

    start += strlen(start_tag);
    end = stristr(start, end_tag);

    if (!end) return EOF;

    return kputsn(start, end - start, out);
}


static int report_s3_error(kstring_t *body, long resp_code) {
    kstring_t entry = KS_INITIALIZE;

    if (get_entry(body->s, "<Code>", "</Code>", &entry) == EOF) {
        return -1;
    }

    fprintf(stderr, "hfile_s3: S3 error %ld: %s\n", resp_code, entry.s);

    ks_clear(&entry);

    if (get_entry(body->s, "<Message>", "</Message>", &entry) == EOF) {
        return -1;
    }

    if (entry.l)
        fprintf(stderr, "%s\n", entry.s);

    ks_free(&entry);

    return 0;
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
    else if (status >= 300)
        return EIO;
    else return 0;
}


static void initialise_local(hFILE_s3 *fp) {
    ks_initialize(&fp->buffer);
    ks_initialize(&fp->url);
    ks_initialize(&fp->upload_id);           // write only
    ks_initialize(&fp->completion_message);  // write only
}


static void cleanup_local(hFILE_s3 *fp) {
    ks_free(&fp->buffer);
    ks_free(&fp->url);
    ks_free(&fp->upload_id);
    ks_free(&fp->completion_message);
    curl_easy_cleanup(fp->curl);
    free_authorisation_values(fp);
}


static void cleanup(hFILE_s3 *fp) {
    // free up authorisation data
    free_auth_data(fp->au);
    cleanup_local(fp);
}

static size_t response_callback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    kstring_t *resp = (kstring_t *)userp;

    if (kputsn((const char *)contents, realsize, resp) == EOF) {
        return 0;
    }

    return realsize;
}


static int add_header(struct curl_slist **head, char *value) {
    int err = 0;
    struct curl_slist *tmp;

    if ((tmp = curl_slist_append(*head, value)) == NULL) {
        err = 1;
    } else {
        *head = tmp;
    }

    return err;
}


static struct curl_slist *set_html_headers(hFILE_s3 *fp, kstring_t *auth, kstring_t *date,
                 kstring_t *content, kstring_t *token, kstring_t *range) {
    struct curl_slist *headers = NULL;
    int err = 0;

    if (auth->l)
        if ((err = add_header(&headers, auth->s)))
            goto error;

    if ((err = add_header(&headers, date->s)))
        goto error;

    if (content->l)
        if ((err = add_header(&headers, content->s)))
            goto error;

    if (range)
        if ((err = add_header(&headers, range->s)))
            goto error;

    if (token->l)
        if ((err = add_header(&headers, token->s)))
            goto error;

    curl_easy_setopt(fp->curl, CURLOPT_HTTPHEADER, headers);

error:

    if (err) {
        curl_slist_free_all(headers);
        headers = NULL;
    }

    return headers;
}


/*

S3 Multipart Upload
-------------------

There are several steps in the Mulitipart upload.


1) Initiate Upload
------------------

Initiate the upload and get an upload ID.  This ID is used in all other steps.


2) Upload Part
--------------

Upload a part of the data.  5Mb minimum part size (except for the last part).
Each part is numbered and a successful upload returns an Etag header value that
needs to used for the completion step.

Step repeated till all data is uploaded.


3) Completion
-------------

Complete the upload by sending all the part numbers along with their associated
Etag values.


Optional - Abort
----------------

If something goes wrong this instructs the server to delete all the partial
uploads and abandon the upload process.
*/

/*
   This is the writing code.
*/

#define MINIMUM_S3_WRITE_SIZE 5242880

// Lets the part memory size grow to about 1Gb giving a 2.5Tb max file size.
// Max. parts allowed by AWS is 10000, so use ceil(10000.0/9.0)
#define EXPAND_ON 1112



/*
    The partially uploaded file will hang around unless the delete command is sent.
*/
static int abort_upload(hFILE_s3 *fp) {
    kstring_t url = KS_INITIALIZE;
    kstring_t canonical_query_string = KS_INITIALIZE;
    int ret = -1;
    struct curl_slist *headers = NULL;
    char http_request[] = "DELETE";
    CURLcode err;

    clear_authorisation_values(fp);

    if (ksprintf(&canonical_query_string, "uploadId=%s", fp->upload_id.s) < 0) {
        goto out;
    }

    if (v4_authorisation(fp,  http_request, NULL, canonical_query_string.s, 0) != 0) {
        goto out;
    }

    if (ksprintf(&url, "%s?%s", fp->url.s, canonical_query_string.s) < 0) {
        goto out;
    }

    if (ksprintf(&fp->content, "x-amz-content-sha256: %s", fp->content_hash.s) < 0) {
        goto out;
    }

    curl_easy_reset(fp->curl);

    err = curl_easy_setopt(fp->curl, CURLOPT_CUSTOMREQUEST, http_request);
    err |= curl_easy_setopt(fp->curl, CURLOPT_USERAGENT, curl.useragent.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_URL, url.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_VERBOSE, fp->verbose);

    if (err != CURLE_OK)
        goto out;

    headers = set_html_headers(fp, &fp->authorisation, &fp->date, &fp->content, &fp->token, NULL);

    if (!headers)
        goto out;

    fp->ret = curl_easy_perform(fp->curl);

    if (fp->ret == CURLE_OK) {
        ret = 0;
    }

 out:
    ks_free(&url);
    ks_free(&canonical_query_string);
    curl_slist_free_all(headers);

    fp->aborted = 1;
    cleanup(fp);

    return ret;
}


static int complete_upload(hFILE_s3 *fp, kstring_t *resp) {
    kstring_t url = KS_INITIALIZE;
    kstring_t canonical_query_string = KS_INITIALIZE;
    int ret = -1;
    struct curl_slist *headers = NULL;
    char http_request[] = "POST";
    CURLcode err;

    clear_authorisation_values(fp);

    if (ksprintf(&canonical_query_string, "uploadId=%s", fp->upload_id.s) < 0) {
        return -1;
    }

    // finish off the completion reply
    if (kputs("</CompleteMultipartUpload>\n", &fp->completion_message) < 0) {
        goto out;
    }

    if (v4_authorisation(fp,  http_request, &fp->completion_message, canonical_query_string.s, 0) != 0) {
        goto out;
    }

    if (ksprintf(&url, "%s?%s", fp->url.s, canonical_query_string.s) < 0) {
        goto out;
    }

    if (ksprintf(&fp->content, "x-amz-content-sha256: %s", fp->content_hash.s) < 0) {
        goto out;
    }

    curl_easy_reset(fp->curl);

    err = curl_easy_setopt(fp->curl, CURLOPT_POST, 1L);
    err |= curl_easy_setopt(fp->curl, CURLOPT_POSTFIELDS, fp->completion_message.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_POSTFIELDSIZE, (long) fp->completion_message.l);
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEFUNCTION, response_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEDATA, (void *)resp);
    err |= curl_easy_setopt(fp->curl, CURLOPT_URL, url.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_USERAGENT, curl.useragent.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_VERBOSE, fp->verbose);

    if (err != CURLE_OK)
        goto out;

    headers = set_html_headers(fp, &fp->authorisation, &fp->date, &fp->content, &fp->token, NULL);

    if (!headers)
        goto out;

    fp->ret = curl_easy_perform(fp->curl);

    if (fp->ret == CURLE_OK) {
        ret = 0;
    }

 out:
    ks_free(&url);
    ks_free(&canonical_query_string);
    curl_slist_free_all(headers);

    return ret;
}


static size_t upload_callback(void *ptr, size_t size, size_t nmemb, void *stream) {
    size_t realsize = size * nmemb;
    hFILE_s3 *fp = (hFILE_s3 *)stream;
    size_t read_length;

    if (realsize > (fp->buffer.l - fp->index)) {
        read_length = fp->buffer.l - fp->index;
    } else {
        read_length = realsize;
    }

    memcpy(ptr, fp->buffer.s + fp->index, read_length);
    fp->index += read_length;

    return read_length;
}


static int upload_part(hFILE_s3 *fp, kstring_t *resp) {
    kstring_t url = KS_INITIALIZE;
    kstring_t canonical_query_string = KS_INITIALIZE;
    int ret = -1;
    struct curl_slist *headers = NULL;
    char http_request[] = "PUT";
    CURLcode err;

    clear_authorisation_values(fp);

    if (ksprintf(&canonical_query_string, "partNumber=%d&uploadId=%s", fp->part_no, fp->upload_id.s) < 0) {
        return -1;
    }

    if (v4_authorisation(fp, http_request, &fp->buffer, canonical_query_string.s, 0) != 0) {
        goto out;
    }

    if (ksprintf(&url, "%s?%s", fp->url.s, canonical_query_string.s) < 0) {
        goto out;
    }

    fp->index = 0;
    if (ksprintf(&fp->content, "x-amz-content-sha256: %s", fp->content_hash.s) < 0) {
        goto out;
    }

    curl_easy_reset(fp->curl);

    err = curl_easy_setopt(fp->curl, CURLOPT_UPLOAD, 1L);
    err |= curl_easy_setopt(fp->curl, CURLOPT_READFUNCTION, upload_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_READDATA, fp);
    err |= curl_easy_setopt(fp->curl, CURLOPT_INFILESIZE_LARGE, (curl_off_t)fp->buffer.l);
    err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERFUNCTION, response_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERDATA, (void *)resp);
    err |= curl_easy_setopt(fp->curl, CURLOPT_URL, url.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_USERAGENT, curl.useragent.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_VERBOSE, fp->verbose);

    if (err != CURLE_OK)
        goto out;

    headers = set_html_headers(fp, &fp->authorisation, &fp->date, &fp->content, &fp->token, NULL);

    if (!headers)
        goto out;

    fp->ret = curl_easy_perform(fp->curl);

    if (fp->ret == CURLE_OK) {
        ret = 0;
    }

 out:
    ks_free(&url);
    ks_free(&canonical_query_string);
    curl_slist_free_all(headers);

    return ret;
}


static ssize_t s3_write(hFILE *fpv, const void *bufferv, size_t nbytes) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;
    const char *buffer  = (const char *)bufferv;
    CURLcode cret;

    if (kputsn(buffer, nbytes, &fp->buffer) == EOF) {
        return -1;
    }

    if (fp->buffer.l > fp->part_size) {
        // time to write out our data
        kstring_t response = {0, 0, NULL};
        int ret;

        ret = upload_part(fp, &response);

        if (!ret) {
            long response_code;
            kstring_t etag = {0, 0, NULL};

            cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);

            if (cret != CURLE_OK || response_code > 200) {
                errno = http_status_errno(response_code);
                ret = -1;
            } else {
                if (get_entry(response.s, "Etag: \"", "\"", &etag) == EOF) {
                    fprintf(stderr, "hfile_s3: Failed to read Etag\n");
                    ret = -1;
                } else {
                    ksprintf(&fp->completion_message, "\t<Part>\n\t\t<PartNumber>%d</PartNumber>\n\t\t<ETag>%s</ETag>\n\t</Part>\n",
                        fp->part_no, etag.s);

                    ks_free(&etag);
                }
            }
        }

        ks_free(&response);

        if (ret) {
            abort_upload(fp);
            return -1;
        }

        fp->part_no++;
        fp->buffer.l = 0;

        if (fp->expand && (fp->part_no % EXPAND_ON == 0)) {
            fp->part_size *= 2;
        }
    }

    return nbytes;
}


static int s3_write_close(hFILE *fpv) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;
    kstring_t response = {0, 0, NULL};
    int ret = 0;
    CURLcode cret;

    if (!fp->aborted) {

        if (fp->buffer.l) {
            // write the last part

            ret = upload_part(fp, &response);

            if (!ret) {
                long response_code;
                kstring_t etag = {0, 0, NULL};

                cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);

                if (cret != CURLE_OK || response_code > 200) {
                    errno = http_status_errno(response_code);
                    ret = -1;
                } else {
                    if (get_entry(response.s, "ETag: \"", "\"", &etag) == EOF) {
                        ret = -1;
                    } else {
                        ksprintf(&fp->completion_message, "\t<Part>\n\t\t<PartNumber>%d</PartNumber>\n\t\t<ETag>%s</ETag>\n\t</Part>\n",
                            fp->part_no, etag.s);

                        ks_free(&etag);
                    }
                }
            }

            ks_free(&response);

            if (ret) {
                abort_upload(fp);
                return -1;
            }

            fp->part_no++;
        }

        if (fp->part_no > 1) {
            ret = complete_upload(fp, &response);

            if (!ret) {
                if (strstr(response.s, "CompleteMultipartUploadResult") == NULL) {
                    ret = -1;
                }
            }
        } else {
            ret = -1;
        }

        if (ret) {
            abort_upload(fp);
        } else {
            cleanup(fp);
        }
    }

    ks_free(&response);

    return ret;
}


static int handle_bad_request(hFILE_s3 *fp, kstring_t *resp) {
    kstring_t region = {0, 0, NULL};
    int ret = -1;

    if (get_entry(resp->s, "<Region>", "</Region>", &region) == EOF) {
        return -1;
    }

    ret = set_region(fp->au, &region);

    ks_free(&region);

    return ret;
}

static int initialise_upload(hFILE_s3 *fp, kstring_t *head, kstring_t *resp, int user_query) {
    kstring_t url = KS_INITIALIZE;
    int ret = -1;
    struct curl_slist *headers = NULL;
    char http_request[] = "POST";
    char delimiter = '?';
    CURLcode err;

    clear_authorisation_values(fp);

    if (user_query) {
        delimiter = '&';
    }

    if (v4_authorisation(fp, http_request, NULL, "uploads=", user_query) != 0) {
        goto out;
    }

    if (ksprintf(&url, "%s%cuploads", fp->url.s, delimiter) < 0) {
        goto out;
    }

    if (ksprintf(&fp->content, "x-amz-content-sha256: %s", fp->content_hash.s) < 0) {
        goto out;
    }

    err = curl_easy_setopt(fp->curl, CURLOPT_URL, url.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_POST, 1L);
    err |= curl_easy_setopt(fp->curl, CURLOPT_POSTFIELDS, "");  // send no data
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEFUNCTION, response_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEDATA, (void *)resp);
    err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERFUNCTION, response_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERDATA, (void *)head);
    err |= curl_easy_setopt(fp->curl, CURLOPT_USERAGENT, curl.useragent.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_VERBOSE, fp->verbose);

    if (err != CURLE_OK)
        goto out;

    headers = set_html_headers(fp, &fp->authorisation, &fp->date, &fp->content, &fp->token, NULL);

    if (!headers)
        goto out;

    fp->ret = curl_easy_perform(fp->curl);

    if (fp->ret == CURLE_OK) {
        ret = 0;
    }

 out:
    curl_slist_free_all(headers);
    ks_free(&url);

    return ret;
}


static int get_upload_id(hFILE_s3 *fp, kstring_t *resp) {
    int ret = 0;

    if (get_entry(resp->s, "<UploadId>", "</UploadId>", &fp->upload_id) == EOF) {
        ret = -1;
    }

    return ret;
}


/*
    Now for the reading code
*/

#define READ_PART_SIZE 1048576

static size_t recv_callback(char *ptr, size_t size, size_t nmemb, void *fpv) {
    hFILE_s3 *fp = (hFILE_s3 *) fpv;
    size_t n = size * nmemb;

    if (n) {
        if (kputsn(ptr, n, &fp->buffer) == EOF) {
            fprintf(stderr, "hfile_s3: error: unable to allocate memory to read data.\n");
            return 0;
        }
    }

    return n;
}


static int s3_read_close(hFILE *fpv) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;

    cleanup(fp);

    return 0;
}


static int get_part(hFILE_s3 *fp, kstring_t *resp) {
    struct curl_slist *headers = NULL;
    int ret = -1;
    char http_request[] = "GET";
    char canonical_query_string = 0;
    CURLcode err;

    ks_clear(&fp->buffer); // reset storage buffer
    clear_authorisation_values(fp);

    if (fp->au->is_v4) {
        if (v4_authorisation(fp, http_request, NULL, &canonical_query_string, 0) != 0) {
            goto out;
        }

        if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: get_part: v4 auth done\n");

        if (ksprintf(&fp->content, "x-amz-content-sha256: %s", fp->content_hash.s) < 0) {
            goto out;
        }
    } else {
        if (v2_authorisation(fp, http_request) != 0) {
            goto out;
        }

        if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: get_part v2 auth done\n");
    }

    if (ksprintf(&fp->range, "Range: bytes=%zu-%zu", fp->last_read, fp->last_read + fp->part_size - 1) < 0) {
        goto out;
    }

    if (hts_verbose >= HTS_LOG_INFO) {
        fprintf(stderr, "hfile_s3: get_part: range set %s\n", fp->range.s);
        fprintf(stderr, "hfile_s3: url %s\n", fp->url.s);
    }

    curl_easy_reset(fp->curl);

    err = curl_easy_setopt(fp->curl, CURLOPT_URL, fp->url.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEFUNCTION, recv_callback);
    err |= curl_easy_setopt(fp->curl, CURLOPT_WRITEDATA, (void *)fp);
    err |= curl_easy_setopt(fp->curl, CURLOPT_USERAGENT, curl.useragent.s);
    err |= curl_easy_setopt(fp->curl, CURLOPT_VERBOSE, fp->verbose);

    if (resp) {
        err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERFUNCTION, response_callback);
        err |= curl_easy_setopt(fp->curl, CURLOPT_HEADERDATA, (void *)resp);
    }

    if (err != CURLE_OK)
        goto out;

    headers = set_html_headers(fp, &fp->authorisation, &fp->date, &fp->content, &fp->token, &fp->range);

    if (!headers)
        goto out;

    fp->ret = curl_easy_perform(fp->curl);

    if (fp->ret == CURLE_OK) {
        ret = 0;
    }

out:
    if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: get_part: ret %d\n", ret);
    curl_slist_free_all(headers);

    return ret;
}


static ssize_t s3_read(hFILE *fpv, void *bufferv, size_t nbytes) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;
    char *buffer = (char *)bufferv;
    size_t got = 0;

    /* Transfer data from the fp->buffer to the calling buffer.
       If there is no data left in the fp->buffer, grab another chunk of
       data from s3.
    */
    while (fp->keep_going && got < nbytes) {

        if (fp->buffer.l && fp->last_read_buffer < fp->buffer.l) {
            // copy data across
            size_t to_copy;
            size_t remaining = fp->buffer.l - fp->last_read_buffer;
            size_t bytes_left = nbytes - got;

            if (hts_verbose >  HTS_LOG_INFO) fprintf(stderr, "hfile_s3: read - remaining %zu read %zu bytes_left %zu, nbytes %zu\n", remaining, got, bytes_left, nbytes);

            if (bytes_left < remaining) {
                to_copy = bytes_left;
            } else {
                to_copy = remaining;
            }

            memcpy(buffer + got, fp->buffer.s + fp->last_read_buffer, to_copy);
            got += to_copy;
            fp->last_read_buffer += to_copy;

            if ((fp->buffer.l < fp->part_size) && (fp->last_read_buffer == fp->buffer.l)) {
                fp->keep_going = 0;
            }
        } else {
            int ret;

            ret = get_part(fp, NULL);

            if (!ret) {
                long response_code;
                CURLcode cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);

                if (cret != CURLE_OK || response_code > 300) {
                    errno = http_status_errno(response_code);
                    ret = -1;
                }
            }

            if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: read - read error %d\n", ret);

            if (ret < 0)
                return ret;

            if (fp->buffer.l == 0) {
                fp->keep_going = 0;
                break;
            }

            fp->last_read_buffer = 0;
            fp->last_read = fp->last_read + fp->buffer.l;
        }
    }

    return got;
}


static off_t s3_seek(hFILE *fpv, off_t offset, int whence) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;
    off_t origin;

    if (fp->write) {
        // lets not try and seek while writing
        errno = ESPIPE;
        return -1;
    }

    // I am not sure we handle any seek other than one from the beginning
    switch (whence) {
        case SEEK_SET:
            origin = 0;
            break;
        case SEEK_CUR:
            // hseek() should convert this to SEEK_SET
            errno = ENOSYS;
            return -1;
        case SEEK_END:
            if (fp->file_size < 0) {
                errno = ESPIPE;
                return -1;
            }

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

    fp->keep_going = 1;

    size_t pos = origin + offset; // origin is really only useful if we can make the other modes work

    if (pos <= fp->last_read && pos > (fp->last_read - fp->buffer.l)) {
        // within the current local buffer
        fp->last_read_buffer = pos - (fp->last_read - fp->buffer.l);
    } else {
        fp->last_read = pos;
        ks_clear(&fp->buffer); // resetting fp->buffer triggers a new remote read
    }

    return fp->last_read;
}


/*
    Unlike upload, download does not really need an initialisation.  Here we use it to
    get the size of the wanted files and as a test for redirects.
*/
static int initialise_download(hFILE_s3 *fp, kstring_t *resp) {

    fp->last_read = 0;
    ks_clear(resp);

    return get_part(fp, resp);
}


static int s3_close(hFILE *fpv) {
    hFILE_s3 *fp = (hFILE_s3 *)fpv;
    int ret;

    if (!fp->write) {
        ret = s3_read_close(fpv);
    } else {
        ret = s3_write_close(fpv);
    }

    return ret;
}


static const struct hFILE_backend s3_backend = {
    s3_read, s3_write, s3_seek, NULL, s3_close
};

/* Read and write open here, need to be after the s3_backend declaration. */
static hFILE *s3_write_open(const char *url, s3_auth_data *auth) {
    hFILE_s3 *fp;
    kstring_t response = {0, 0, NULL};
    kstring_t header   = {0, 0, NULL};
    int ret, has_user_query = 0;
    char *query_start;
    const char *env;
    CURLcode cret;
    long response_code;


    fp = (hFILE_s3 *)hfile_init(sizeof(hFILE_s3), "w", 0);

    if (fp == NULL) {
        return NULL;
    }

    if ((fp->curl = curl_easy_init()) == NULL) {
        errno = ENOMEM;
        goto error;
    }

    fp->au = auth;

    initialise_local(fp);
    initialise_authorisation_values(fp);
    fp->aborted = 0;
    fp->part_size = MINIMUM_S3_WRITE_SIZE;
    fp->expand = 1;
    fp->write = 1;

    if ((env = getenv("HTS_S3_PART_SIZE")) != NULL) {
        int part_size = atoi(env) * 1024 * 1024;

        if (part_size > fp->part_size)
            fp->part_size = part_size;

        fp->expand = 0;
    }

    if (hts_verbose >= 8) {
        fp->verbose = 1L;
    } else {
        fp->verbose = 0L;
    }

    kputs(url, &fp->url);

    if ((query_start = strchr(fp->url.s, '?'))) {
        has_user_query = 1;;
    }

    ret = initialise_upload(fp, &header, &response, has_user_query);
    cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);

    if (ret == 0) {
        if (cret == CURLE_OK) {
            if (response_code == S3_MOVED_PERMANENTLY) {
                if (redirect_endpoint(fp, &header) == 0) {
                    ks_clear(&response);
                    ks_clear(&header);

                    ret = initialise_upload(fp, &header, &response, has_user_query);
                }
            } else if (response_code == S3_BAD_REQUEST) {
                if (handle_bad_request(fp, &response) == 0) {
                    ks_clear(&response);
                    ks_clear(&header);

                    ret = initialise_upload(fp, &header, &response, has_user_query);
                }
            }
        } else {
            // unable to get a response code from curl
            ret = -1;
        }
    }

    if (response_code >= 300) {
        // something went wrong with the initialisation

        if (cret == CURLE_OK) {
            if (hts_verbose >= HTS_LOG_INFO) {
                if (report_s3_error(&response, response_code)) {
                    fprintf(stderr, "hfile_s3: warning, unable to report full S3 error status.\n");
                }
            }

            errno = http_status_errno(response_code);
        }

        ret = -1;
    }

    if (ret) goto error;

    if (get_upload_id(fp, &response)) goto error;

    // start the completion message (a formatted list of parts)
    if (kputs("<CompleteMultipartUpload>\n", &fp->completion_message) == EOF) {
        goto error;
    }

    fp->part_no = 1;

    // user query string no longer a useful part of the URL
    if (query_start)
         *query_start = '\0';

    fp->base.backend = &s3_backend;
    ks_free(&response);
    ks_free(&header);

    return &fp->base;

error:
    ks_free(&response);
    ks_free(&header);
    cleanup_local(fp);
    free_authorisation_values(fp);
    hfile_destroy((hFILE *)fp);
    return NULL;
}


static hFILE *s3_read_open(const char *url, s3_auth_data *auth) {
    hFILE_s3 *fp;
    const char *env;
    kstring_t response   = {0, 0, NULL};
    kstring_t file_range = {0, 0, NULL};
    int ret;
    CURLcode cret;
    long response_code = 0;

    fp = (hFILE_s3 *)hfile_init(sizeof(hFILE_s3), "r", 0);

    if (fp == NULL) {
        return NULL;
    }

    if ((fp->curl = curl_easy_init()) == NULL) {
        errno = ENOMEM;
        goto error;
    }

    fp->au = auth;

    initialise_local(fp);
    initialise_authorisation_values(fp);

    fp->last_read = 0; // ranges start at 0
    fp->write = 0;

    if ((env = getenv("HTS_S3_READ_PART_SIZE")) != NULL) {
        fp->part_size = atoi(env) * 1024 * 1024;
    } else {
        fp->part_size = READ_PART_SIZE;
    }

    if (hts_verbose >= 8) {
        fp->verbose = 1L;
    } else {
        fp->verbose = 0L;
    }

    kputs(url, &fp->url);

    ret = initialise_download(fp, &response);
    cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);

    if (ret == 0) {
        if (cret == CURLE_OK) {
            if (response_code == S3_MOVED_PERMANENTLY) {
                ks_clear(&response);

                if (redirect_endpoint(fp, &response) == 0) {
                    ret = initialise_download(fp, &response);
                }
            } else if (response_code == S3_BAD_REQUEST) {
                ks_clear(&response);

                if (handle_bad_request(fp, &fp->buffer) == 0) {
                    ret = initialise_download(fp, &response);
                }
            }

            // reget the response code (may not have changed)
            cret = curl_easy_getinfo(fp->curl, CURLINFO_RESPONSE_CODE, &response_code);
        } else {
            // unable to get a response code from curl
            ret = -1;
        }
    }

    if (response_code >= 300) {
        // something went wrong with the initialisation

        if (cret == CURLE_OK) {
            if (hts_verbose >= HTS_LOG_INFO) {
                if (report_s3_error(&fp->buffer, response_code)) {
                    fprintf(stderr, "hfile_s3: warning, unable to report full S3 error status.\n");
                }
            }

            errno = http_status_errno(response_code);
        }

        ret = -1;
    }

    if (ret) goto error;

    if (get_entry(response.s, "content-range: bytes ", "\n", &file_range) == EOF) {
        fprintf(stderr, "hfile_s3: warning: failed to read file size.\n");
        fp->file_size = -1;
    } else {
        char *s;
        if ((s = strchr(file_range.s, '/'))) {
            fp->file_size = strtoll(s + 1, NULL, 10);
        } else {
            fp->file_size = -1;
        }
    }

    fp->last_read_buffer = 0;
    fp->last_read = fp->last_read + fp->buffer.l;
    fp->base.backend = &s3_backend;
    fp->keep_going = 1;

    ks_free(&response);
    ks_free(&file_range);
    return &fp->base;


 error:
    ks_free(&response);
    ks_free(&file_range);
    cleanup_local(fp);
    free_authorisation_values(fp);
    hfile_destroy((hFILE *)fp);
    return NULL;
}


static hFILE *s3_open_v4(const char *s3url, const char *mode, va_list *argsp) {
    kstring_t url = { 0, 0, NULL };

    s3_auth_data *ad = setup_auth_data(s3url, mode, 4, &url);
    hFILE *fp = NULL;

    if (ad == NULL) {
        return NULL;
    }

    if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: s3_open_v4 url %s\n", url.s);

    if (*mode == 'r') {
        fp  = s3_read_open(url.s, ad);
    } else {
        fp =  s3_write_open(url.s, ad);
    }

    ks_free(&url);
    if (!fp)
        free_auth_data(ad);

    return fp;
}


static hFILE *s3_open_v2(const char *s3url, const char *mode, va_list *argsp) {
    kstring_t url = { 0, 0, NULL };

    s3_auth_data *ad = setup_auth_data(s3url, mode, 2, &url);
    hFILE *fp = NULL;

    if (ad == NULL) {
        return NULL;
    }

    if (hts_verbose >= HTS_LOG_INFO) fprintf(stderr, "hfile_s3: s3_open_v2 url %s\n", url.s);

    if (*mode == 'r') {
        fp  = s3_read_open(url.s, ad);
    } else {
        fprintf(stderr, "hfile_s3: error - signature v2 not handled for writing.\n.");
    }

    ks_free(&url);
    if (!fp)
        free_auth_data(ad);

    return fp;
}


static hFILE *hopen_s3(const char *url, const char *mode)
{
    hFILE *fp;

    if (getenv("HTS_S3_V2") == NULL) { // Force the v2 signature code
        fp = s3_open_v4(url, mode, NULL);
    } else {
        fp = s3_open_v2(url, mode, NULL);
    }

    return fp;
}


static hFILE *vhopen_s3(const char *url, const char *mode, va_list args0)
{
    hFILE *fp;

    // This should handle to vargs case.  Not sure what vargs we want
    // to handle
    fp = hopen_s3(url, mode);

    return fp;
}


static void s3_exit(void) {
    if (curl_share_cleanup(curl.share) == CURLSHE_OK)
        curl.share = NULL;

    free(curl.useragent.s);
    curl.useragent.l = curl.useragent.m = 0; curl.useragent.s = NULL;
    curl_global_cleanup();
}


int PLUGIN_GLOBAL(hfile_plugin_init,_s3)(struct hFILE_plugin *self) {

    static const struct hFILE_scheme_handler handler =
        { hopen_s3, hfile_always_remote, "Amazon S3",
          2000 + 50, vhopen_s3
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] =
        "@(#)hfile_s3 plugin (htslib)\t" HTS_VERSION_TEXT;
    const char *version = strchr(id, '\t') + 1;

    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_s3.init] version %s\n",
                version);
#else
    const char *version = hts_version();
#endif

    const curl_version_info_data *info;
    CURLcode err;
    CURLSHcode errsh;

    err = curl_global_init(CURL_GLOBAL_ALL);

    if (err != CURLE_OK) {
        // look at putting in an errno here
        return -1;
    }

    curl.share = curl_share_init();

    if (curl.share == NULL) {
        curl_global_cleanup();
        errno = EIO;
        return -1;
    }

    errsh  = curl_share_setopt(curl.share, CURLSHOPT_LOCKFUNC, share_lock);
    errsh |= curl_share_setopt(curl.share, CURLSHOPT_UNLOCKFUNC, share_unlock);
    errsh |= curl_share_setopt(curl.share, CURLSHOPT_SHARE, CURL_LOCK_DATA_DNS);

    if (errsh != 0) {
        curl_share_cleanup(curl.share);
        curl_global_cleanup();
        errno = EIO;
        return -1;
    }

    info = curl_version_info(CURLVERSION_NOW);
    ksprintf(&curl.useragent, "htslib/%s libcurl/%s", version, info->version);

    self->name = "Amazon S3";
    self->destroy = s3_exit;

    hfile_add_scheme_handler("s3",       &handler);
    hfile_add_scheme_handler("s3+http",  &handler);
    hfile_add_scheme_handler("s3+https", &handler);

    return 0;
}

