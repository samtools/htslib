/*  upstream.c -- download ref-cache files from upstream host

    Copyright (C) 2025 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <curl/curl.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <signal.h>
#include <limits.h>
#include <assert.h>

#include "upstream.h"
#include "cmsg_wrap.h"
#include "misc.h"
#include "options.h"
#include "poll_wrap.h"
#include "../htslib/hts_defs.h"

// This avoids having to include all of htslib/hts.h
typedef struct hts_md5_context hts_md5_context;
extern hts_md5_context *hts_md5_init(void);
extern void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size);
extern void hts_md5_final(unsigned char *result, hts_md5_context *ctx);
extern void hts_md5_destroy(hts_md5_context *ctx);


static unsigned int curl_version_num;

#define BUF_SZ 4096
#define MD5_LEN 32

enum {
    DL_OK       = 1,
    DL_CLENGTH  = 2,
    DL_WAITING  = 4,
    DL_ABANDON  = 8,
};

struct Multi_data;
struct Downstream;



typedef struct Download {
    char               hexmd5[MD5_LEN];
    hts_md5_context   *md5_ctx;
    struct Multi_data *mdata;
    struct Download   *next;
    struct Download   *waiting;
    struct Downstream *downstream;
    const char        *cache_dir;
    char              *file;
    char              *url;
    CURL              *curl;
    int                cmd_fd;
    int                curlid;
    int                flags;
    int                cache_fd;
    int                file_fd;
    off_t              size;
    off_t              received;
} Download;

typedef struct Downstream {
    struct Download   *download;
    struct Downstream *prev;
    struct Downstream *next;
    int                cmd_fd;
    unsigned int       id;
} Downstream;

/* Lookup table of downloads by md5 */
#define ACTIVE_SIZE 0x4000
#define ACTIVE_MASK 0x3fff

typedef struct Multi_data {
    CURLM         *multi;
    Poll_wrap     *pw;
    long           timeout;
    Download     **downloads;
    Download      *waiting;
    Download      *last_waiting;
    unsigned int   ncurls;
    unsigned int   free_curls;
    int            running;
    CURL         **curls;
} Multi_data;

static int handle_curl_socket(CURLM *multi, Multi_data *mdata,
                              int events, int fd);

int upstream_send_cmd(int cmd_fd, const char hexmd5[MD5_LEN], unsigned int id) {
    struct msghdr msg;
    struct iovec  iov[2];
    ssize_t res;

    memset(&msg, 0, sizeof(msg));
    msg.msg_name = NULL;
    msg.msg_control = NULL;

    /* The MD5 for the reference we want */
    iov[0].iov_base   = (void *) hexmd5;
    iov[0].iov_len    = MD5_LEN;
    iov[1].iov_base   = &id;
    iov[1].iov_len    = sizeof(id);
    msg.msg_iov    = iov;
    msg.msg_iovlen = 2;

    do {
        res = sendmsg(cmd_fd, &msg, 0);
    } while (res == -1
             && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
    return res == MD5_LEN + sizeof(int) ? 0 : -1;
}

static ssize_t recv_cmd_data(int cmd_fd, char hexmd5[MD5_LEN], unsigned int *id) {
    ssize_t res;
    struct msghdr msg;
    struct iovec  iov[2];

    memset(&msg, 0, sizeof(msg));
    msg.msg_name = NULL;
    msg.msg_control = NULL;

    iov[0].iov_base   = hexmd5;
    iov[0].iov_len    = MD5_LEN;
    iov[1].iov_base   = id;
    iov[1].iov_len    = sizeof(*id);
    msg.msg_iov    = iov;
    msg.msg_iovlen = 2;

    do {
        res = recvmsg(cmd_fd, &msg, 0);
    } while (res == -1
             && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
    if (res == 0) return 0; /* Shutdown */
    if (res != MD5_LEN + sizeof(*id)) return -1;

    return res;
}

static int upstream_send_msg(int cmd_fd, Upstream_msg *umsg, int fd) {
    ssize_t res;
    struct iovec iov[1];
    struct msghdr msg;
    char   buf[256]; // Must be bigger than CMSG_SPACE(sizeof(int))

    memset(&msg, 0, sizeof(msg));
    msg.msg_name = NULL;
    msg.msg_control = NULL;

    memset(&buf, 0, sizeof(buf));

    iov[0].iov_base = umsg;
    iov[0].iov_len  = sizeof(*umsg);
    msg.msg_iov     = iov;
    msg.msg_iovlen  = 1;

    if (umsg->code == US_START) {
        /* Add the descriptor for the downloaded file */
        if (make_scm_rights_cmsg(&msg, fd, buf, sizeof(buf)) < 0) {
            fprintf(stderr, "upstream_send_msg: cmsg buffer not big enough.\n");
            return -1;
        }
    }

    do {
        res = sendmsg(cmd_fd, &msg, 0);
    } while (res == -1
             && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
    return res == sizeof(*umsg) ? 0 : -1;
}

static int send_msg_all(Download *download,
                        Upstream_msg_code code, int64_t val) {
    Downstream *d;
    int res = 0;
    for (d = download->downstream; d != NULL; d = d->next) {
        Upstream_msg msg = { d->id, code, val };
        if (upstream_send_msg(d->cmd_fd, &msg, -1) < 0) res = -1;
    }
    return res;
}

int upstream_recv_msg(int cmd_fd, Upstream_msg *umsg, int *fd) {
    struct msghdr msg;
    struct iovec   iov[1];
    ssize_t res;
    char   buf[16384];

    memset(&msg, 0, sizeof(msg));
    msg.msg_name = NULL;
    msg.msg_control = NULL;

    iov[0].iov_base   = umsg;
    iov[0].iov_len    = sizeof(*umsg);
    msg.msg_iov    = iov;
    msg.msg_iovlen = 1;
    msg.msg_control = buf;
    msg.msg_controllen = sizeof(buf);
    memset(umsg, 0, sizeof(*umsg));

    do {
        res = recvmsg(cmd_fd, &msg, 0);
    } while (res == -1
             && (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK));
    if (res == 0) return 0; /* Shutdown */
    if (res != sizeof(*umsg)) return -1;

    if (umsg->code == US_START) {
        /* Message should include a file descriptor */
        *fd = get_scm_rights_fd(&msg);
        if (*fd < 0) {
            fprintf(stderr, "Failed to get file descriptor in upstream message\n");
            return -1;
        }
    }

    return 1;
}

static int make_subdir(Options *opts, char *hexmd5) {
    char path[6];

    memcpy(path, hexmd5, 2); path[2] = '\0';
    if (mkdirat(opts->cache_fd, path, 01755) != 0) {
        if (errno != EEXIST) {
            fprintf(stderr, "Couldn't make directory %s/%s : %s\n",
                    opts->cache_dir, path, strerror(errno));
            return -1;
        }
    }

    path[2] = '/'; memcpy(path + 3, hexmd5 + 2, 2); path[5] = '\0';
    if (mkdirat(opts->cache_fd, path, 01755) != 0) {
        if (errno != EEXIST) {
            fprintf(stderr, "Couldn't make directory %s/%s : %s\n",
                    opts->cache_dir, path, strerror(errno));
            return -1;
        }
    }
    return 0;
}



static int get_free_curl(Multi_data *mdata) {
    unsigned int i;

    if (mdata->free_curls == 0) return -1;
    for (i = 0; i < mdata->ncurls; i++) {
        if (mdata->free_curls & (1U << i)) {
            mdata->free_curls &= ~(1U << i);
            return (int) i;
        }
    }
    return -1; /* Should never happen */
}

static void release_curl(Download *download) {
    download->mdata->free_curls |= (1U << download->curlid);
    download->curlid = -1;
    download->curl   = NULL;
}

static inline Downstream * new_downstream(int cmd_fd,
                                          unsigned int downstream_id) {
    Downstream *downstream = calloc(1, sizeof(Downstream));

    if (downstream == NULL) { perror("new_downstream"); return NULL; }
    downstream->cmd_fd = cmd_fd;
    downstream->id     = downstream_id;
    return downstream;
}

Download *new_download(Options *opts, Multi_data *mdata, char hexmd5[MD5_LEN]) {
    Download *download = calloc(1, sizeof(Download));
    if (download == NULL) { perror("new_download"); return NULL; }

    download->mdata = mdata;
    download->cache_dir = opts->cache_dir;
    memcpy(download->hexmd5, hexmd5, MD5_LEN);
    download->cache_fd = opts->cache_fd;
    download->file_fd = -1;
    download->curlid = -1;
    return download;
}

static void remove_downstream(Downstream *downstream, Download *download) {
    /* Remove the downstream record from the list */
    if (downstream->prev == NULL) {
        download->downstream = downstream->next;
    } else {
        assert(downstream != download->downstream);
        downstream->prev->next = downstream->next;
    }
    if (downstream->next != NULL) downstream->next->prev = downstream->prev;
    free(downstream);
}

static void free_download(Download *download) {
    int hash;
    Multi_data *mdata = download->mdata;

    /* Remove downstream information */
    while (download->downstream != NULL) {
        remove_downstream(download->downstream, download);
    }

    /* Remove from hash table (if in there) */
    hash = (hexval(download->hexmd5[0]) << 12
            | hexval(download->hexmd5[1]) << 8
            | hexval(download->hexmd5[2]) << 4
            | hexval(download->hexmd5[3])) & ACTIVE_MASK;

    if (mdata->downloads[hash] == download) {
        mdata->downloads[hash] = download->next;
    } else {
        Download *d = mdata->downloads[hash];
        while (d != NULL && d->next != download) d = d->next;
        if (d != NULL) d->next = download->next;
    }

    if (download->flags & DL_WAITING) {
        /* Remove from waiting list */
        if (mdata->waiting == download) {
            if (mdata->last_waiting == download) {
                mdata->waiting = mdata->last_waiting = NULL;
            } else {
                mdata->waiting = download->waiting;
            }
        } else {
            Download *d = mdata->waiting;
            while (d != NULL && d->waiting != download) d = d->waiting;
            if (d != NULL) {
                if (mdata->last_waiting == download) mdata->last_waiting = d;
                d->waiting = download->waiting;
            }
        }
    }

    if (download->file_fd != -1) {
        close(download->file_fd);
        if ((download->flags & DL_OK) == 0) {
            unlinkat(download->cache_fd, download->file, 0);
        }
    }
    free(download->file);
    free(download->url);
    if (download->curlid != -1) release_curl(download);
    if (download->md5_ctx != NULL) hts_md5_destroy(download->md5_ctx);
    free(download);
}

static int start_new_download(CURLM *multi, Multi_data *mdata) {
    Download *download = mdata->waiting;
    CURLcode cc;
    CURLMcode mc;

    if (download == NULL) return 0;
    assert((download->flags & DL_WAITING) != 0);

    download->curlid = get_free_curl(mdata);
    if (download->curlid == -1) return 0;

    download->curl = mdata->curls[download->curlid];
    mdata->waiting = download->waiting;
    if (mdata->last_waiting == download) mdata->last_waiting = NULL;
    download->flags &= ~DL_WAITING;

    cc = curl_easy_setopt(download->curl, CURLOPT_URL, download->url);
    if (cc != CURLE_OK) {
        fprintf(stderr, "Couldn't set URL %s : %s\n", download->url,
                curl_easy_strerror(cc));
        free_download(download);
        return -1;
    }
    cc = curl_easy_setopt(download->curl, CURLOPT_WRITEDATA, download);
    if (cc != CURLE_OK) {
        fprintf(stderr, "Couldn't set user data in CURL handle : %s\n",
                curl_easy_strerror(cc));
        free_download(download);
        return -1;
    }
    cc = curl_easy_setopt(download->curl, CURLOPT_PRIVATE, download);
    if (cc != CURLE_OK) {
        fprintf(stderr, "Couldn't set private data in CURL handle : %s\n",
                curl_easy_strerror(cc));
        free_download(download);
        return -1;
    }
#if LIBCURL_VERSION_NUM >= 0x072000 // 7.32.0
    cc = curl_easy_setopt(download->curl, CURLOPT_XFERINFODATA, download);
#else
    cc = curl_easy_setopt(download->curl, CURLOPT_PROGRESSDATA, download);
#endif
    if (cc != CURLE_OK) {
        fprintf(stderr, "Couldn't set progress data in CURL handle : %s\n",
                curl_easy_strerror(cc));
        free_download(download);
        return -1;
    }

    mc = curl_multi_add_handle(multi, download->curl);
    if (mc != CURLM_OK) {
        fprintf(stderr, "Couldn't add handle to curl_multi : %s\n",
                curl_multi_strerror(mc));
        free_download(download);
        return -1;
    }
    mdata->running++;
    return 1;
}

static Download *get_download(Options *opts, CURLM *multi, Multi_data *mdata,
                              char hexmd5[MD5_LEN]) {
    unsigned int hash;
    Download *download;
    struct stat st;
    size_t url_len;
    int count, need_sep;

    /* Check for existing downloads */
    hash = (hexval(hexmd5[0]) << 12
            | hexval(hexmd5[1]) << 8
            | hexval(hexmd5[2]) << 4
            | hexval(hexmd5[3])) & ACTIVE_MASK;
    download = mdata->downloads[hash];

    while (download != NULL) {
        if (memcmp(hexmd5, download->hexmd5, MD5_LEN) == 0) return download;
        download = download->next;
    }

    /* Need to start a new one */
    download = new_download(opts, mdata, hexmd5);
    if (download == NULL) return NULL;

    /* Look for completed downloads */
    download->file = malloc(MD5_LEN + 16);
    if (download->file == NULL) {
        perror("Allocating download->file");
        free(download);
        return NULL;
    }

    snprintf(download->file, MD5_LEN + 16, "%.2s/%.2s/%.28s",
             download->hexmd5, download->hexmd5 + 2, download->hexmd5 + 4);

    download->file_fd = openat(opts->cache_fd, download->file, O_RDONLY);
    if (download->file_fd >= 0) {
        if (fstat(download->file_fd, &st) != 0) {
            fprintf(stderr, "Couldn't stat %s/%s : %s\n",
                    opts->cache_dir, download->file, strerror(errno));
            free_download(download);
            return NULL;
        }
        download->size     = st.st_size;
        download->received = st.st_size;
        download->flags = DL_OK | DL_CLENGTH;
        download->next = mdata->downloads[hash];
        mdata->downloads[hash] = download;
        return download;
    }

    if (errno != ENOENT) {
        fprintf(stderr, "Couldn't open %s/%s : %s\n",
                opts->cache_dir, download->file, strerror(errno));
        free_download(download);
        return NULL;
    }

    /* Need to start a new download */
    url_len = opts->upstream_url_len + MD5_LEN + 2;
    need_sep = (opts->upstream_url_len == 0
                || opts->upstream_url[opts->upstream_url_len - 1] != '/');
    download->url = malloc(url_len);
    if (download->url == NULL) {
        perror("Allocating download->url");
        free_download(download);
        return NULL;
    }
    snprintf(download->url, url_len,
             "%s%s%.32s", opts->upstream_url, need_sep ? "/" : "", hexmd5);

    /* Ensure the subdir exists */
    if (make_subdir(opts, hexmd5) != 0) {
        free_download(download);
        return NULL;
    }

    /* Open the destination file */
    for (count = 0; count < 1000; count++) {
        snprintf(download->file, MD5_LEN + 16, "%.2s/%.2s/%.28s.%03d",
                 hexmd5, hexmd5 + 2, hexmd5 + 4, count);
        do {
            download->file_fd = openat(opts->cache_fd, download->file,
                                       O_RDWR | O_CREAT | O_EXCL, 0644);
        } while (download->file_fd == -1 && errno == EINTR);
        if (download->file_fd >= 0) break;
        if (errno != EEXIST) break;
    }
    if (download->file_fd == -1) {
        fprintf(stderr, "Couldn't open %s/%s for writing: %s\n",
                opts->cache_dir, download->file, strerror(errno));
        free_download(download);
        return NULL;
    }

    /* Set up MD5 calculation */
    download->md5_ctx = hts_md5_init();
    if (download->md5_ctx == NULL) {
        free_download(download);
        return NULL;
    }

    /* Put on waiting list */
    download->waiting = NULL;
    if (mdata->last_waiting != NULL) mdata->last_waiting->waiting = download;
    if (mdata->waiting == NULL) mdata->waiting = download;
    mdata->last_waiting = download;
    download->flags |= DL_WAITING;

    /* Try to start it off */
    if (start_new_download(multi, mdata) < 0) return NULL;

    download->next = mdata->downloads[hash];
    mdata->downloads[hash] = download;

    return download;
}

static int get_cmd_multi(int cmd_fd, Options *opts,
                         CURLM *multi, Multi_data *mdata, int running) {
    char hexmd5[MD5_LEN] = { 0 };
    Downstream *downstream = NULL;
    Download *download;
    Upstream_msg msg;
    ssize_t clen;
    int res = -1;
    int downstream_fd;
    unsigned int downstream_id = 0;

    clen = recv_cmd_data(cmd_fd, hexmd5, &downstream_id);
    if (clen < 0)  return -1;
    if (clen == 0) return  0;

    if (!running) {
        /* Tell the other end we can't help */
        msg.id   = downstream_id;
        msg.code = US_RESULT;
        msg.val  = 503; /* Service unavailable */
        upstream_send_msg(cmd_fd, &msg, -1);
    }

    download = get_download(opts, multi, mdata, hexmd5);
    if (download == NULL) goto out;

    /* Check if we already know about this one */
    for (downstream = download->downstream;
         downstream != NULL && downstream->cmd_fd != cmd_fd;
         downstream = downstream->next) {}

    if (downstream != NULL) {
        return 1;
    }

    downstream = new_downstream(cmd_fd, downstream_id);
    if (downstream == NULL) goto out;

    downstream_fd = dup(download->file_fd);
    if (downstream_fd < 0) goto out;

    /* Send message to say the download is starting and pass on the file
       descriptor */
    msg.id   = downstream_id;
    msg.code = US_START;
    msg.val  = (download->flags & DL_CLENGTH) != 0 ? download->size : -1;
    if (upstream_send_msg(cmd_fd, &msg, downstream_fd) < 0) {
        close(downstream_fd);
        goto out;
    }
    close(downstream_fd);

    downstream->download = download;
    downstream->prev = NULL;
    downstream->next = download->downstream;
    if (download->downstream != NULL) download->downstream->prev = downstream;
    download->downstream = downstream;

    res = 1;

 out:
    if (res != 1) {
        /* Tell the other end we can't help */
        msg.id   = downstream_id;
        msg.code = US_RESULT;
        msg.val  = 500;
        upstream_send_msg(cmd_fd, &msg, -1);
        free(downstream);
    }
    return res;
}

static int send_result_code(Downstream *downstream, int code) {
    Upstream_msg msg = { downstream->id, US_RESULT, code };
    if (upstream_send_msg(downstream->cmd_fd, &msg, -1) != 0) {
        fprintf(stderr, "Error sending result code to downstream #%u : %s\n",
                downstream->id, strerror(errno));
        return -1;
    }

    return 0;
}

static int send_result_code_all(Download *download, int code) {
    int res = 0;
    Downstream *d;
    for (d = download->downstream; d != NULL; d = d->next) {
        res |= send_result_code(d, code);
    }
    return res;
}

#if LIBCURL_VERSION_NUM >= 0x072000 // 7.32.0
static int progress_callback(void *clientp,
                             curl_off_t dltotal HTS_UNUSED,
                             curl_off_t dlnow HTS_UNUSED,
                             curl_off_t ultotal HTS_UNUSED,
                             curl_off_t ulnow HTS_UNUSED) {
    Download *download = (Download *) clientp;

    return (download->flags & DL_ABANDON) == 0 ? 0 : 1;
}
#else
static int progress_callback(void *clientp,
                             double dltotal HTS_UNUSED, double dlnow HTS_UNUSED,
                             double ultotal HTS_UNUSED, double ulnow HTS_UNUSED) {
    Download *download = (Download *) clientp;

    return (download->flags & DL_ABANDON) == 0 ? 0 : 1;
}
#endif

static size_t multi_receive_data(void *buffer, size_t size, size_t nmemb,
                                 void *userp) {
    Download *download = (Download *) userp;
    size_t bytes = size * nmemb, in, out;
    unsigned char *ucb = (unsigned char *) buffer;

    if ((download->flags & DL_ABANDON) != 0) return bytes;

    if ((download->flags & DL_CLENGTH) == 0) {
#if LIBCURL_VERSION_NUM >= 0x073700 // 7.55.0
        curl_off_t clen;
#else
        double clen;
#endif
        long rcode = 500;

        /* Check the result code.  It's not worth downloading the file if
           it isn't 200 */
        if (curl_easy_getinfo(download->curl,
                              CURLINFO_RESPONSE_CODE, &rcode) != CURLE_OK) {
            return 0;
        }
        if (rcode != 200) {
            if (send_result_code_all(download, (int) rcode) != 0) return 0;
            download->flags |= DL_ABANDON;
            return bytes;
        }
#if LIBCURL_VERSION_NUM >= 0x073700 // 7.55.0
        if (curl_easy_getinfo(download->curl,
                              CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
                              &clen) != CURLE_OK) {
            return 0;
        }
#else
        if (curl_version_num >= 0x071904) {
            /* Get content-length if available */
            if (curl_easy_getinfo(download->curl,
                                  CURLINFO_CONTENT_LENGTH_DOWNLOAD,
                                  &clen) != CURLE_OK) {
                return 0;
            }
            if (clen < 0) clen = 0;
        } else {
            clen = 0; /* Not reliable before 7.19.4 */
        }
#endif
        download->flags |= DL_CLENGTH;
        download->size = (off_t) clen;

        /* Send on to upstream */
        if (send_msg_all(download, US_CONTENT_LENGTH, download->size) != 0) {
            return 0;
        }
    }

    /* Write the data to the file */
    if (do_write_all(download->file_fd, buffer, bytes) != 0) {
        fprintf(stderr, "Error writing to %s/%s : %s\n",
                download->cache_dir, download->file, strerror(errno));
        return 0;
    }

    /* MD5 calculation - assumes we can mess with the contents of buffer */
    for (in = out = 0; in < bytes; in++) {
        if (!isspace(ucb[in])) ucb[out++] = (unsigned char) toupper(ucb[in]);
    }
    if (out > 0) hts_md5_update(download->md5_ctx, ucb, out);

    download->received += (off_t) bytes;

    /* Tell upstream the new length */
    if (send_msg_all(download, US_PARTIAL_LENGTH, download->received) != 0) {
        return 0;
    }

    return bytes;
}

static int sock_func(CURL *easy HTS_UNUSED,
                     curl_socket_t s, int action,
                     void *userp, void *socketp) {
    Multi_data *mdata = (Multi_data *) userp;
    Pw_item *polled = (Pw_item *) socketp;
    int res;
    CURLMcode mc;

    if (action == CURL_POLL_REMOVE) {
        assert(polled != NULL);
        res = pw_remove(mdata->pw, polled, 0);
        if (res != 0) {
            /* Ugh! Libcurl may or may not have closed the file descriptor
               before this function gets called.  So the only thing we can
               do is ignore the error we get if it has been closed already. */
            if (errno == EBADF) return 0;
            perror("Removing file descriptor from poller");
            return res;
        }
        mc = curl_multi_assign(mdata->multi, s, NULL);
        if (mc != CURLM_OK) {
            fprintf(stderr, "curl_multi_assign failed : %s\n",
                    curl_multi_strerror(mc));
            return -1;
        }
        return 0;
    }

    if (polled == NULL) {
        /* Register the file descriptor */
        polled = pw_register(mdata->pw, s, US_CURL,
                             ((action & CURL_POLL_IN) ? PW_IN : 0
                              | (action & CURL_POLL_OUT) ? PW_OUT : 0), NULL);
        if (polled == NULL) return -1;
        mc = curl_multi_assign(mdata->multi, s, polled);
        if (mc != CURLM_OK) {
            fprintf(stderr, "curl_multi_assign failed : %s\n",
                    curl_multi_strerror(mc));
            return -1;
        }
        return 0;
    } else {
        return pw_mod(mdata->pw, polled,
                      (action & CURL_POLL_IN) ? PW_IN : 0
                      | (action & CURL_POLL_OUT) ? PW_OUT : 0);
    }
}

static int timer_func(CURLM *multi HTS_UNUSED, long timeout, void *userp) {
    Multi_data *mdata = (Multi_data *) userp;

    mdata->timeout = timeout;
    return 0;
}

static CURLM *get_multi(Multi_data *mdata) {
    CURLM *multi = curl_multi_init();
    CURLMcode mc;

    if (multi == NULL) {
        fprintf(stderr, "curl_multi_init() failed");
        return NULL;
    }

    mc = curl_multi_setopt(multi, CURLMOPT_SOCKETFUNCTION, sock_func);

    if (mc == CURLM_OK)
        mc = curl_multi_setopt(multi, CURLMOPT_SOCKETDATA, mdata);

    if (mc == CURLM_OK)
        mc = curl_multi_setopt(multi, CURLMOPT_TIMERFUNCTION, timer_func);

    if (mc == CURLM_OK)
        mc = curl_multi_setopt(multi, CURLMOPT_TIMERDATA, mdata);

    if (mc != CURLM_OK) {
        fprintf(stderr, "Failed to set options for curl multi handle : %s\n",
                curl_multi_strerror(mc));
        curl_multi_cleanup(multi);
        return NULL;
    }

    return multi;
}

static int init_multi_data(CURLM *multi, Multi_data *mdata) {
    unsigned int i;
    mdata->multi     = multi;
    mdata->pw        = pw_init(0);
    if (mdata->pw == NULL) {
        perror("Initalizing poller");
        return -1;
    }
    mdata->downloads = NULL;
    mdata->waiting = mdata->last_waiting = NULL;
    mdata->timeout   = -1;
    mdata->downloads = calloc(ACTIVE_SIZE, sizeof(Download *));
    if (mdata->downloads == NULL) { perror(""); return -1; }
    mdata->ncurls    = 4;
    mdata->free_curls = 0;
    mdata->running = 0;
    mdata->curls = calloc(mdata->ncurls, sizeof(CURL *));
    if (mdata->curls == NULL) { perror(""); goto fail; }

    for (i = 0; i < mdata->ncurls; i++) {
        CURLcode  cc;

        mdata->curls[i] = curl_easy_init();
        if (mdata->curls[i] == NULL) goto fail;
        cc = curl_easy_setopt(mdata->curls[i],
                              CURLOPT_WRITEFUNCTION, multi_receive_data);
        if (cc != CURLE_OK) {
            fprintf(stderr, "Couldn't set WRITEFUNCTION on curl handle : %s\n",
                    curl_easy_strerror(cc));
            goto fail;
        }
        cc = curl_easy_setopt(mdata->curls[i], CURLOPT_NOPROGRESS, 0L);
#if LIBCURL_VERSION_NUM >= 0x072000 // 7.32.0
        cc = curl_easy_setopt(mdata->curls[i],
                              CURLOPT_XFERINFOFUNCTION, progress_callback);
#else
        cc = curl_easy_setopt(mdata->curls[i],
                              CURLOPT_PROGRESSFUNCTION, progress_callback);
#endif
        if (cc != CURLE_OK) {
            fprintf(stderr, "Couldn't set progress callback in CURL handle : %s\n",
                    curl_easy_strerror(cc));
            goto fail;
        }
        mdata->free_curls |= 1U << i;
    }
    return 0;

 fail:
    if (mdata->curls != NULL) {
        for (i = 0; i < mdata->ncurls; i++) {
            if (mdata->curls[i] != NULL) curl_easy_cleanup(mdata->curls[i]);
        }
        free(mdata->curls);
    }
    free(mdata->downloads);
    return -1;
}

static int rename_download_file(Download *download) {
    if (download->flags & DL_OK) {
        char dest[MD5_LEN + 3];

        snprintf(dest, MD5_LEN + 3, "%.2s/%.2s/%.28s",
                 download->hexmd5, download->hexmd5 + 2, download->hexmd5 + 4);
        if (renameat(download->cache_fd, download->file,
                     download->cache_fd, dest) != 0) {
            fprintf(stderr, "Couldn't rename %s/%s to %s/%s: %s\n",
                    download->cache_dir, download->file,
                    download->cache_dir, dest, strerror(errno));
            download->flags &= ~DL_OK;
            return -1;
        }
        memcpy(download->file, dest, MD5_LEN + 3); /* Will always be shorter */
    }
    return 0;
}

static int finish_download(Download *download, CURLM *multi, CURLMsg *msg) {
    long rcode  = 500;
    CURLcode cc = msg->data.result;
    unsigned char md5[16];
    CURLMcode mc;
    int i;

    /* Check for curl errors */

    if (cc != CURLE_OK) {
        fprintf(stderr, "Download of %.32s failed: %s\n",
                download->hexmd5, curl_easy_strerror(cc));
        goto bad_download;
    }

    /* Check the response code */

    cc = curl_easy_getinfo(msg->easy_handle, CURLINFO_RESPONSE_CODE, &rcode);
    if (cc != CURLE_OK) {
        fprintf(stderr, "Couldn't get response code : %s\n",
                curl_easy_strerror(cc));
        return -1;
    }

    if (rcode != 200) goto bad_download;

    /* Check the content length was as expected */

    if (download->size == 0) download->size = download->received;
    if (download->received != download->size) {
        fprintf(stderr, "Downloading %.32s : Content-Length was a lie. "
                "Expected %ld, got %ld\n",
                download->hexmd5, (long) download->size, (long) download->received);
        rcode = 502;
        goto bad_download;
    }

    /* Check the MD5 */
    hts_md5_final(md5, download->md5_ctx);
    for (i = 0; i < 16; i++) {
        int byte = (hexval(download->hexmd5[i * 2]) << 4 |
                    hexval(download->hexmd5[i * 2 + 1]));
        if (byte != md5[i]) {
            fprintf(stderr, "Downloading %.32s : MD5 checksum didn't match.\n",
                    download->hexmd5);
            rcode = 502;
            goto bad_download;
        }
    }

    download->flags |= DL_OK;
    if (rename_download_file(download) != 0) return -1;

    mc = curl_multi_remove_handle(multi, download->curl);
    if (mc != CURLM_OK) {
        fprintf(stderr, "Couldn't remove easy handle from curl_multi : %s\n",
                curl_multi_strerror(mc));
        return -1;
    }

    /* Tell everyone that we have finished */
    if (send_result_code_all(download, 200) != 0) {
        return -1;
    }

    free_download(download);
    return 0;

 bad_download:
    if (send_result_code_all(download, (int) rcode) != 0) return -1;
    mc = curl_multi_remove_handle(multi, download->curl);
    if (mc != CURLM_OK) {
        fprintf(stderr, "Couldn't remove easy handle from curl_multi : %s\n",
                curl_multi_strerror(mc));
        return -1;
    }
    free_download(download);
    return 0;
}

static int handle_curl_socket(CURLM *multi, Multi_data *mdata,
                              int events, int fd) {
    CURLMcode mc;
    int running = 0;
    do {
        mc = curl_multi_socket_action(multi, fd, events, &running);
    } while (mc == CURLM_CALL_MULTI_PERFORM);

    if (mc != CURLM_OK) {
        fprintf(stderr, "Error from curl socket #%d: %s\n",
                fd, curl_multi_strerror(mc));
        return -1;
    }
    if (running < mdata->running) {
        CURLMsg *msg;
        int msgs_in_queue = 0;
        int nfinished = 0, started;

        mdata->running = running;

        while ((msg = curl_multi_info_read(multi, &msgs_in_queue)) != NULL) {
            if (msg->msg == CURLMSG_DONE) {
                CURLcode c;
                Download *download;

                c = curl_easy_getinfo(msg->easy_handle,
                                      CURLINFO_PRIVATE, (char **) &download);
                if (c != CURLE_OK) {
                    fprintf(stderr, "curl_easy_getinfo failed: %s\n",
                            curl_easy_strerror(c));
                    return -1;
                }
                if (finish_download(download, multi, msg) != 0) return -1;
                nfinished++;
            }
        }

        if (nfinished) {
            do {
                started = start_new_download(multi, mdata);
            } while (started > 0);
            if (started < 0) return -1;
        }
    } else {
        mdata->running = running;
        if (mdata->running == 0) mdata->timeout = -1;
    }
    return 0;
}

#define MAX_EVENTS 128
static void run_epoll_loop(Options *opts, unsigned int nfds, int *cmd_fds,
                           int liveness_fd, CURLM *multi, Multi_data *mdata) {
    Pw_events events[MAX_EVENTS];
    Pw_item **cmd_pollers;
    unsigned int npolled, i;
    int running = 1;

    cmd_pollers = malloc((nfds + 1) * sizeof(Pw_item *));
    if (cmd_pollers == NULL) return;

    for (npolled = 0; npolled < nfds; npolled++) {
        cmd_pollers[npolled] = pw_register(mdata->pw, cmd_fds[npolled],
                                           US_COMMAND, PW_IN, NULL);
        if (cmd_pollers[npolled] == NULL)
            goto cleanup;
    }

    cmd_pollers[nfds] = pw_register(mdata->pw, liveness_fd,
                                    US_LIVE, PW_HUP|PW_ERR, NULL);
    if (cmd_pollers[nfds] == NULL)
        goto cleanup;

    while (running || mdata->running > 0) {
        int n;
        int nevents = pw_wait(mdata->pw, events, MAX_EVENTS,
                              mdata->timeout < INT_MAX
                              ? (int) mdata->timeout : INT_MAX);
        if (nevents == -1) {
            if (errno == EINTR) continue;  /* Signal */
            perror("poll_wait");
            goto cleanup;
        }

        if (nevents == 0) { /* Timeout */
            if (handle_curl_socket(multi, mdata, 0, CURL_SOCKET_TIMEOUT) != 0) {
                goto cleanup;
            }
        }

        for (n = 0; n < nevents; n++) {
            unsigned int evts = PWE(events[n]);
            Pw_item *polled   = PWI(events[n]);
            switch (polled->fd_type) {
            case US_COMMAND:
                if (opts->verbosity > 2) {
                    fprintf(stderr, "Upstream received command\n");
                }
                if (get_cmd_multi(polled->fd, opts, multi, mdata, running) <= 0) {
                    goto cleanup;
                }
                break;
            case US_CURL: {
                int e = (( (evts  & PW_IN)  ? CURL_CSELECT_IN  : 0)
                         | ((evts & PW_OUT) ? CURL_CSELECT_OUT : 0)
                         | ((evts & PW_ERR) ? CURL_CSELECT_ERR : 0));
                if (opts->verbosity > 2) {
                    fprintf(stderr, "Upstream received curl event %d on fd #%d\n",
                            e, polled->fd);
                }
                if (handle_curl_socket(multi, mdata, e, polled->fd) != 0) {
                    goto cleanup;
                }
                break;
            }
            case US_LIVE:
                if (evts & (PW_HUP | PW_ERR)) {
                    running = 0;
                }
                break;
            default:
                break;
            }
        }
    }

 cleanup:
    for (i = 0; i < npolled; i++) {
        if (cmd_pollers[i] != NULL)
            pw_remove(mdata->pw, cmd_pollers[i], 0);
    }
    free(cmd_pollers);
    return;
}

static int run_multi_upstream_handler(Options *opts, int *cmd_fds,
                                      int liveness_fd) {
    int res = -1;
    Multi_data mdata;
    CURLM *multi = get_multi(&mdata);
    unsigned int i;

    if (multi == NULL) return -1;

    if (init_multi_data(multi, &mdata) != 0) {
        curl_multi_cleanup(multi);
        return -1;
    }

    run_epoll_loop(opts, opts->max_kids, cmd_fds, liveness_fd, multi, &mdata);

    for (i = 0; i < mdata.ncurls; i++) {
        curl_multi_remove_handle(multi, mdata.curls[i]);
        curl_easy_cleanup(mdata.curls[i]);
    }

    free(mdata.downloads);
    curl_multi_cleanup(multi);
    return res;
}

int run_upstream_handler(Options *opts, int *sockets, int liveness_fd) {
    int res;
    struct sigaction sa;
    struct sigaction old_sa;
    curl_version_info_data *cvi = NULL;

    sa.sa_handler = SIG_IGN;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;

    if (sigaction(SIGPIPE, &sa, &old_sa) != 0) {
        perror("sigaction(SIGPIPE)");
        return -1;
    }

    if (curl_global_init(CURL_GLOBAL_ALL) != 0) {
        fprintf(stderr, "Couldn't initialize libcurl\n");
        return -1;
    }

    cvi = curl_version_info(CURLVERSION_NOW);
    if (cvi == NULL) {
        fprintf(stderr, "Couldn't get curl version information\n");
        return -1;
    }
    curl_version_num = cvi->version_num;

    res = run_multi_upstream_handler(opts, sockets, liveness_fd);

    curl_global_cleanup();

    if (sigaction(SIGPIPE, &old_sa, NULL) != 0) {
        perror("sigaction(SIGPIPE)");
    }

    return res;
}
