/*  hfile_libcurl.c -- libcurl backend for low-level file streams.

    Copyright (C) 2015 Genome Research Ltd.

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

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <sys/select.h>

#include "hfile_internal.h"
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"

#include <curl/curl.h>

typedef struct {
    hFILE base;
    CURL *easy;
    struct curl_slist *headers;
    off_t file_size;
    struct {
        union { char *rd; const char *wr; } ptr;
        size_t len;
    } buffer;
    CURLcode final_result;  // easy result code for finished transfers
    // Flags for communicating with libcurl callbacks:
    unsigned paused : 1;    // callback tells us that it has paused transfer
    unsigned closing : 1;   // informs callback that hclose() has been invoked
    unsigned finished : 1;  // wait_perform() tells us transfer is complete
} hFILE_libcurl;

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

static int easy_errno(CURL *easy, CURLcode err)
{
    long lval;

    switch (err) {
    case CURLE_OK:
        return 0;

    case CURLE_UNSUPPORTED_PROTOCOL:
    case CURLE_URL_MALFORMAT:
        return EINVAL;

    case CURLE_NOT_BUILT_IN:
        return ENOSYS;

    case CURLE_COULDNT_RESOLVE_PROXY:
    case CURLE_COULDNT_RESOLVE_HOST:
    case CURLE_FTP_CANT_GET_HOST:
        return EDESTADDRREQ; // Lookup failure

    case CURLE_COULDNT_CONNECT:
    case CURLE_SEND_ERROR:
    case CURLE_RECV_ERROR:
        if (curl_easy_getinfo(easy, CURLINFO_OS_ERRNO, &lval) == CURLE_OK)
            return lval;
        else
            return ECONNABORTED;

    case CURLE_REMOTE_ACCESS_DENIED:
    case CURLE_LOGIN_DENIED:
    case CURLE_TFTP_PERM:
        return EACCES;

    case CURLE_PARTIAL_FILE:
        return EPIPE;

    case CURLE_HTTP_RETURNED_ERROR:
        if (curl_easy_getinfo(easy, CURLINFO_RESPONSE_CODE, &lval) == CURLE_OK)
            return http_status_errno(lval);
        else
            return EIO;

    case CURLE_WRITE_ERROR:
    case CURLE_READ_ERROR:
        return ENOTRECOVERABLE; // Indicates bugs in our callback routines

    case CURLE_OUT_OF_MEMORY:
        return ENOMEM;

    case CURLE_OPERATION_TIMEDOUT:
        return ETIMEDOUT;

    case CURLE_RANGE_ERROR:
        return ESPIPE;

    case CURLE_SSL_CONNECT_ERROR:
        // TODO return SSL error buffer messages
        return ECONNABORTED;

    case CURLE_FILE_COULDNT_READ_FILE:
    case CURLE_TFTP_NOTFOUND:
        return ENOENT;

    case CURLE_TOO_MANY_REDIRECTS:
        return ELOOP;

    case CURLE_FILESIZE_EXCEEDED:
        return EFBIG;

    case CURLE_REMOTE_DISK_FULL:
        return ENOSPC;

    case CURLE_REMOTE_FILE_EXISTS:
        return EEXIST;

    default:
        return EIO;
    }
}

static int multi_errno(CURLMcode errm)
{
    switch (errm) {
    case CURLM_CALL_MULTI_PERFORM:
    case CURLM_OK:
        return 0;

    case CURLM_BAD_HANDLE:
    case CURLM_BAD_EASY_HANDLE:
    case CURLM_BAD_SOCKET:
        return EBADF;

    case CURLM_OUT_OF_MEMORY:
        return ENOMEM;

    default:
        return EIO;
    }
}


static struct {
    CURLM *multi;
    kstring_t useragent;
    int nrunning;
    unsigned perform_again : 1;
} curl = { NULL, { 0, 0, NULL }, 0, 0 };

static void libcurl_exit()
{
    (void) curl_multi_cleanup(curl.multi);
    curl.multi = NULL;

    free(curl.useragent.s);
    curl.useragent.l = curl.useragent.m = 0; curl.useragent.s = NULL;

    curl_global_cleanup();
}


static void process_messages()
{
    CURLMsg *msg;
    int remaining;

    while ((msg = curl_multi_info_read(curl.multi, &remaining)) != NULL) {
        hFILE_libcurl *fp = NULL;
        curl_easy_getinfo(msg->easy_handle, CURLINFO_PRIVATE, (char **) &fp);
        switch (msg->msg) {
        case CURLMSG_DONE:
            fp->finished = 1;
            fp->final_result = msg->data.result;
            break;

        default:
            break;
        }
    }
}

static int wait_perform()
{
    fd_set rd, wr, ex;
    int maxfd, nrunning;
    long timeout;
    CURLMcode errm;

    FD_ZERO(&rd);
    FD_ZERO(&wr);
    FD_ZERO(&ex);
    if (curl_multi_fdset(curl.multi, &rd, &wr, &ex, &maxfd) != CURLM_OK)
        maxfd = -1, timeout = 1000;
    else if (maxfd < 0)
        timeout = 100;  // as recommended by curl_multi_fdset(3)
    else {
        if (curl_multi_timeout(curl.multi, &timeout) != CURLM_OK)
            timeout = 1000;
        else if (timeout < 0)
            timeout = 10000;  // as recommended by curl_multi_timeout(3)
    }

    if (timeout > 0 && ! curl.perform_again) {
        struct timeval tval;
        tval.tv_sec  = (timeout / 1000);
        tval.tv_usec = (timeout % 1000) * 1000;

        if (select(maxfd + 1, &rd, &wr, &ex, &tval) < 0) return -1;
    }

    errm = curl_multi_perform(curl.multi, &nrunning);
    curl.perform_again = 0;
    if (errm == CURLM_CALL_MULTI_PERFORM) curl.perform_again = 1;
    else if (errm != CURLM_OK) { errno = multi_errno(errm); return -1; }

    if (nrunning < curl.nrunning) process_messages();
    return 0;
}


static size_t recv_callback(char *ptr, size_t size, size_t nmemb, void *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    size_t n = size * nmemb;

    if (n > fp->buffer.len) { fp->paused = 1; return CURL_WRITEFUNC_PAUSE; }
    else if (n == 0) return 0;

    memcpy(fp->buffer.ptr.rd, ptr, n);
    fp->buffer.ptr.rd += n;
    fp->buffer.len -= n;
    return n;
}

static ssize_t libcurl_read(hFILE *fpv, void *bufferv, size_t nbytes)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    char *buffer = (char *) bufferv;
    CURLcode err;

    fp->buffer.ptr.rd = buffer;
    fp->buffer.len = nbytes;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    while (! fp->paused && ! fp->finished)
        if (wait_perform() < 0) return -1;

    nbytes = fp->buffer.ptr.rd - buffer;
    fp->buffer.ptr.rd = NULL;
    fp->buffer.len = 0;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return nbytes;
}

static size_t send_callback(char *ptr, size_t size, size_t nmemb, void *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    size_t n = size * nmemb;

    if (fp->buffer.len == 0) {
        // Send buffer is empty; normally pause, or signal EOF if we're closing
        if (fp->closing) return 0;
        else { fp->paused = 1; return CURL_READFUNC_PAUSE; }
    }

    if (n > fp->buffer.len) n = fp->buffer.len;
    memcpy(ptr, fp->buffer.ptr.wr, n);
    fp->buffer.ptr.wr += n;
    fp->buffer.len -= n;
    return n;
}

static ssize_t libcurl_write(hFILE *fpv, const void *bufferv, size_t nbytes)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    const char *buffer = (const char *) bufferv;
    CURLcode err;

    fp->buffer.ptr.wr = buffer;
    fp->buffer.len = nbytes;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    while (! fp->paused && ! fp->finished)
        if (wait_perform() < 0) return -1;

    nbytes = fp->buffer.ptr.wr - buffer;
    fp->buffer.ptr.wr = NULL;
    fp->buffer.len = 0;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return nbytes;
}

static off_t libcurl_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;

    CURLcode err;
    CURLMcode errm;
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

    pos = origin + offset;

    errm = curl_multi_remove_handle(curl.multi, fp->easy);
    if (errm != CURLM_OK) { errno = multi_errno(errm); return -1; }
    curl.nrunning--;

    // TODO If we seem to be doing random access, use CURLOPT_RANGE to do
    // limited reads (e.g. about a BAM block!) so seeking can reuse the
    // existing connection more often.

    if (pos <= 2147483647) err = curl_easy_setopt(fp->easy, CURLOPT_RESUME_FROM, (long) pos);
    else err = curl_easy_setopt(fp->easy, CURLOPT_RESUME_FROM_LARGE, (curl_off_t) pos);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    fp->buffer.len = 0;
    fp->paused = fp->finished = 0;

    errm = curl_multi_add_handle(curl.multi, fp->easy);
    if (errm != CURLM_OK) { errno = multi_errno(errm); return -1; }
    curl.nrunning++;

    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    while (! fp->paused && ! fp->finished)
        if (wait_perform() < 0) return -1;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return pos;
}

static int libcurl_close(hFILE *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    CURLcode err;
    CURLMcode errm;
    int save_errno = 0;

    // Before closing the file, unpause it and perform on it so that uploads
    // have the opportunity to signal EOF to the server -- see send_callback().

    fp->buffer.len = 0;
    fp->closing = 1;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) save_errno = easy_errno(fp->easy, err);

    while (save_errno == 0 && ! fp->paused && ! fp->finished)
        if (wait_perform() < 0) save_errno = errno;

    if (fp->finished && fp->final_result != CURLE_OK)
        save_errno = easy_errno(fp->easy, fp->final_result);

    errm = curl_multi_remove_handle(curl.multi, fp->easy);
    if (errm != CURLM_OK && save_errno == 0) save_errno = multi_errno(errm);
    curl.nrunning--;

    curl_easy_cleanup(fp->easy);

    if (save_errno) { errno = save_errno; return -1; }
    else return 0;
}

static const struct hFILE_backend libcurl_backend =
{
    libcurl_read, libcurl_write, libcurl_seek, NULL, libcurl_close
};

static int add_header(hFILE_libcurl *fp, const char *header)
{
    struct curl_slist *list = curl_slist_append(fp->headers, header);
    if (list == NULL) { errno = ENOMEM; return -1; }
    fp->headers = list;
    return 0;
}

static int
add_s3_settings(hFILE_libcurl *fp, const char *url, kstring_t *message);

hFILE *hopen_libcurl(const char *url, const char *modes)
{
    hFILE_libcurl *fp;
    char mode;
    const char *s;
    CURLcode err;
    CURLMcode errm;
    int save;

    if ((s = strpbrk(modes, "rwa+")) != NULL) {
        mode = *s;
        if (strpbrk(&s[1], "rwa+")) mode = 'e';
    }
    else mode = '\0';

    if (mode != 'r' && mode != 'w') { errno = EINVAL; return NULL; }

    fp = (hFILE_libcurl *) hfile_init(sizeof (hFILE_libcurl), modes, 0);
    if (fp == NULL) return NULL;

    fp->easy = curl_easy_init();
    if (fp->easy == NULL) { errno = ENOMEM; goto error; }

    fp->headers = NULL;
    fp->file_size = -1;
    fp->buffer.ptr.rd = NULL;
    fp->buffer.len = 0;
    fp->final_result = (CURLcode) -1;
    fp->paused = fp->closing = fp->finished = 0;

    // Make a route to the hFILE_libcurl* given just a CURL* easy handle
    err = curl_easy_setopt(fp->easy, CURLOPT_PRIVATE, fp);

    if (mode == 'r') {
        err |= curl_easy_setopt(fp->easy, CURLOPT_WRITEFUNCTION, recv_callback);
        err |= curl_easy_setopt(fp->easy, CURLOPT_WRITEDATA, fp);
    }
    else {
        err |= curl_easy_setopt(fp->easy, CURLOPT_READFUNCTION, send_callback);
        err |= curl_easy_setopt(fp->easy, CURLOPT_READDATA, fp);
        err |= curl_easy_setopt(fp->easy, CURLOPT_UPLOAD, 1L);
        if (add_header(fp, "Transfer-Encoding: chunked") < 0) goto error;
    }

    if (tolower(url[0]) == 's' && url[1] == '3') {
        // Construct the HTTP-Method/Content-MD5/Content-Type part of the
        // message to be signed.  This will be destroyed by add_s3_settings().
        kstring_t message = { 0, 0, NULL };
        kputs((mode == 'r')? "GET\n" : "PUT\n", &message);
        kputc('\n', &message);
        kputc('\n', &message);
        if (add_s3_settings(fp, url, &message) < 0) goto error;
    }
    else
        err |= curl_easy_setopt(fp->easy, CURLOPT_URL, url);

    err |= curl_easy_setopt(fp->easy, CURLOPT_USERAGENT, curl.useragent.s);
    if (fp->headers)
        err |= curl_easy_setopt(fp->easy, CURLOPT_HTTPHEADER, fp->headers);
    err |= curl_easy_setopt(fp->easy, CURLOPT_FOLLOWLOCATION, 1L);
    err |= curl_easy_setopt(fp->easy, CURLOPT_FAILONERROR, 1L);
    if (hts_verbose >= 8)
        err |= curl_easy_setopt(fp->easy, CURLOPT_VERBOSE, 1L);

    if (err != 0) { errno = ENOSYS; goto error; }

    errm = curl_multi_add_handle(curl.multi, fp->easy);
    if (errm != CURLM_OK) { errno = multi_errno(errm); goto error; }
    curl.nrunning++;

    while (! fp->paused && ! fp->finished)
        if (wait_perform() < 0) goto error_remove;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        goto error_remove;
    }

    if (mode == 'r') {
        double dval;
        if (curl_easy_getinfo(fp->easy, CURLINFO_CONTENT_LENGTH_DOWNLOAD,
                              &dval) == CURLE_OK && dval >= 0.0)
            fp->file_size = (off_t) (dval + 0.1);
    }

    fp->base.backend = &libcurl_backend;
    return &fp->base;

error_remove:
    save = errno;
    (void) curl_multi_remove_handle(curl.multi, fp->easy);
    curl.nrunning--;
    errno = save;

error:
    save = errno;
    curl_easy_cleanup(fp->easy);
    if (fp->headers) curl_slist_free_all(fp->headers);
    hfile_destroy((hFILE *) fp);
    errno = save;
    return NULL;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_libcurl)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { hopen_libcurl, hfile_always_remote, "libcurl", 50 };

    const curl_version_info_data *info;
    const char * const *protocol;
    CURLcode err;

    err = curl_global_init(CURL_GLOBAL_ALL);
    if (err != CURLE_OK) { errno = easy_errno(NULL, err); return -1; }

    curl.multi = curl_multi_init();
    if (curl.multi == NULL) { curl_global_cleanup(); errno = EIO; return -1; }

    info = curl_version_info(CURLVERSION_NOW);
    ksprintf(&curl.useragent, "htslib/%s libcurl/%s",
             hts_version(), info->version);

    curl.nrunning = 0;
    curl.perform_again = 0;
    self->name = "libcurl";
    self->destroy = libcurl_exit;

    for (protocol = info->protocols; *protocol; protocol++)
        hfile_add_scheme_handler(*protocol, &handler);

    hfile_add_scheme_handler("s3", &handler);
    hfile_add_scheme_handler("s3+http", &handler);
    if (info->features & CURL_VERSION_SSL)
        hfile_add_scheme_handler("s3+https", &handler);

    return 0;
}


/*******************
 * Rewrite S3 URLs *
 *******************/

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
urldecode_kput(const char *s, int len, hFILE_libcurl *fp, kstring_t *str)
{
    if (memchr(s, '%', len) != NULL) {
        int len2;
        char *s2 = curl_easy_unescape(fp->easy, s, len, &len2);
        if (s2 == NULL) abort();
        kputsn(s2, len2, str);
        curl_free(s2);
    }
    else kputsn(s, len, str);
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
        if (islower(*s))
            has_nondigit = 1;
        else if (*s == '-') {
            has_nondigit = 1;
            if (s == s0 || s+1 == slim) return 0;
        }
        else if (isdigit(*s))
            ;
        else if (*s == '.') {
            if (s == s0 || ! isalnum(s[-1])) return 0;
            if (s+1 == slim || ! isalnum(s[1])) return 0;
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

            while (isspace(*key)) key++;
            while (s > key && isspace(s[-1])) s--;
            *s = '\0';

            while (isspace(*value)) value++;
            while (line.l > 0 && isspace(line.s[line.l-1]))
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
    while (isspace(*s)) s++;
    kputsn(s, len = strcspn(s, " \t"), id);

    s += len;
    while (isspace(*s)) s++;
    kputsn(s, strcspn(s, " \t"), secret);

    free(text.s);
}

static int
add_s3_settings(hFILE_libcurl *fp, const char *s3url, kstring_t *message)
{
    int ret, save;
    const char *bucket, *path;
    char date_hdr[40];
    CURLcode err;

    kstring_t url = { 0, 0, NULL };
    kstring_t profile = { 0, 0, NULL };
    kstring_t id = { 0, 0, NULL };
    kstring_t secret = { 0, 0, NULL };
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

    strftime(date_hdr, sizeof date_hdr, "Date: %a, %d %b %Y %H:%M:%S GMT", tm);
    if (add_header(fp, date_hdr) < 0) goto error;
    kputs(&date_hdr[6], message);
    kputc('\n', message);

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
            urldecode_kput(bucket, colon - bucket, fp, &profile);
        }
        else {
            const char *colon2 = strpbrk(&colon[1], ":@");
            urldecode_kput(bucket, colon - bucket, fp, &id);
            urldecode_kput(&colon[1], colon2 - &colon[1], fp, &secret);
            if (*colon2 == ':')
                urldecode_kput(&colon2[1], path - &colon2[1], fp, &token);
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

    // Use virtual hosted-style access if possible, otherwise path-style.
    if (is_dns_compliant(bucket, path)) {
        kputsn(bucket, path - bucket, &url);
        kputs(".s3.amazonaws.com", &url);
    }
    else {
        kputs("s3.amazonaws.com/", &url);
        kputsn(bucket, path - bucket, &url);
    }
    kputs(path, &url);

    if (id.l == 0) {
        const char *v = getenv("AWS_SHARED_CREDENTIALS_FILE");
        parse_ini(v? v : "~/.aws/credentials", profile.s,
                  "aws_access_key_id", &id, "aws_secret_access_key", &secret,
                  "aws_session_token", &token, NULL);
    }
    if (id.l == 0)
        parse_ini("~/.s3cfg", profile.s, "access_key", &id,
                  "secret_key", &secret, "access_token", &token, NULL);
    if (id.l == 0)
        parse_simple("~/.awssecret", &id, &secret);

    if (token.l > 0) {
        kputs("x-amz-security-token:", message);
        kputs(token.s, message);
        kputc('\n', message);

        kputs("X-Amz-Security-Token: ", &token_hdr);
        kputs(token.s, &token_hdr);
        if (add_header(fp, token_hdr.s) < 0) goto error;
    }

    kputc('/', message);
    kputs(bucket, message); // CanonicalizedResource is '/' + bucket + path

    err = curl_easy_setopt(fp->easy, CURLOPT_URL, url.s);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); goto error; }

    // If we have no id/secret, we can't sign the request but will
    // still be able to access public data sets.
    if (id.l > 0 && secret.l > 0) {
        unsigned char digest[DIGEST_BUFSIZ];
        size_t digest_len = s3_sign(digest, &secret, message);

        kputs("Authorization: AWS ", &auth_hdr);
        kputs(id.s, &auth_hdr);
        kputc(':', &auth_hdr);
        base64_kput(digest, digest_len, &auth_hdr);

        if (add_header(fp, auth_hdr.s) < 0) goto error;
    }

    ret = 0;
    goto free_and_return;

error:
    ret = -1;

free_and_return:
    save = errno;
    free(url.s);
    free(profile.s);
    free(id.s);
    free(secret.s);
    free(token.s);
    free(token_hdr.s);
    free(auth_hdr.s);
    free(message->s);
    errno = save;
    return ret;
}
