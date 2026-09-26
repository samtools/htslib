/*  test/test_hfile_s3.c -- Test cases for the hfile_s3 backend.

    Copyright (C) 2026 Peter Dowdy.

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


/*
 * Needs a live S3-compatible endpoint, so this is run by "make check-s3"
 * rather than "make check", and skips unless HTS_S3_HOST is set.  Reading and
 * writing file formats via S3 is covered by test/s3/s3.tst.
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>

#include "../htslib/hfile.h"

#define MiB (1024 * 1024)

static int failures = 0;
static char bucket[256];
static char run_id[64];

#define PASS(name) fprintf(stderr, "  PASS: %s\n", (name))
#define FAIL(name, ...) do { \
    fprintf(stderr, "  FAIL: %s: ", (name)); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    failures++; \
} while (0)

static void s3_url(char *buf, size_t bufsz, const char *key)
{
    snprintf(buf, bufsz, "s3+http://%s/%s/%s", bucket, run_id, key);
}

static void generate_bytes(unsigned char *buf, size_t n, unsigned seed)
{
    size_t i;
    for (i = 0; i < n; i++)
        buf[i] = (unsigned char)((i * 7 + seed) & 0xFF);
}

// Write len bytes to url in hwrite() calls of at most chunk bytes
static int put(const char *url, const void *data, size_t len, size_t chunk)
{
    const char *p = data;
    size_t off, n;
    hFILE *fp = hopen(url, "w");

    if (!fp) return -1;
    for (off = 0; off < len; off += n) {
        n = len - off < chunk ? len - off : chunk;
        if (hwrite(fp, p + off, n) != (ssize_t) n) {
            hclose_abruptly(fp);
            return -1;
        }
    }
    return hclose(fp) == 0 ? 0 : -1;
}

// Read up to len bytes from url, returning the number read or -1
static ssize_t get(const char *url, void *buf, size_t len)
{
    char *p = buf;
    size_t got = 0;
    ssize_t n = 0;
    hFILE *fp = hopen(url, "r");

    if (!fp) return -1;
    while (got < len && (n = hread(fp, p + got, len - got)) > 0)
        got += n;
    if (hclose(fp) != 0 || n < 0) return -1;
    return got;
}

static void test_roundtrip(const char *name, const char *key,
                           size_t len, size_t chunk)
{
    char url[512];
    unsigned char *wbuf = malloc(len), *rbuf = malloc(len + 1);
    ssize_t got;

    if (!wbuf || !rbuf) {
        FAIL(name, "out of memory");
        goto out;
    }
    s3_url(url, sizeof(url), key);
    generate_bytes(wbuf, len, len & 0xff);

    if (put(url, wbuf, len, chunk) < 0)
        FAIL(name, "write failed: %s", strerror(errno));
    else if ((got = get(url, rbuf, len + 1)) < 0)
        FAIL(name, "read failed: %s", strerror(errno));
    else if ((size_t) got != len)
        FAIL(name, "read %zd bytes, expected %zu", got, len);
    else if (memcmp(wbuf, rbuf, len) != 0)
        FAIL(name, "data mismatch");
    else
        PASS(name);

 out:
    free(wbuf);
    free(rbuf);
}

// Check hopen(url, "r") fails with errno err (or any errno if err is 0)
static void test_open_fails(const char *name, const char *url, int err)
{
    hFILE *fp = hopen(url, "r");

    if (fp) {
        FAIL(name, "hopen unexpectedly succeeded");
        hclose_abruptly(fp);
    } else if (err && errno != err) {
        FAIL(name, "expected errno %d (%s), got %d (%s)",
             err, strerror(err), errno, strerror(errno));
    } else {
        PASS(name);
    }
}

static void test_missing(void)
{
    char url[512];

    s3_url(url, sizeof(url), "does/not/exist.bin");
    test_open_fails("missing key", url, ENOENT);
    test_open_fails("missing bucket",
                    "s3+http://htslib-test-no-such-bucket/probe", ENOENT);
}

static void test_wrong_credentials(void)
{
    const char *name = "wrong credentials";
    static const char data[] = "wrong-credentials-probe";
    char url[512];

    s3_url(url, sizeof(url), "wrong_creds.bin");
    if (put(url, data, sizeof(data), sizeof(data)) < 0) {
        FAIL(name, "setup write failed: %s", strerror(errno));
        return;
    }

    // Credentials in the URL override AWS_ACCESS_KEY_ID etc.
    snprintf(url, sizeof(url), "s3+http://wrongid:wrongsecret@%s/%s/%s",
             bucket, run_id, "wrong_creds.bin");
    test_open_fails(name, url, EACCES);
}

static void test_connection_refused(void)
{
    char url[512], *host = strdup(getenv("HTS_S3_HOST"));

    setenv("HTS_S3_HOST", "127.0.0.1:1", 1);
    s3_url(url, sizeof(url), "unreachable.bin");
    test_open_fails("connection refused", url, 0);
    if (host) setenv("HTS_S3_HOST", host, 1);
    free(host);
}

// hfile_s3 only supports SigV2 for reads
static void test_sigv2_read(void)
{
    const char *name = "SigV2 read";
    static const char data[] = "sigv2-probe";
    char url[512], rbuf[64];
    ssize_t got;

    s3_url(url, sizeof(url), "sigv2.bin");
    if (put(url, data, sizeof(data), sizeof(data)) < 0) {
        FAIL(name, "setup write failed: %s", strerror(errno));
        return;
    }

    setenv("HTS_S3_V2", "1", 1);
    got = get(url, rbuf, sizeof(rbuf));
    unsetenv("HTS_S3_V2");

    if (got < 0)
        FAIL(name, "read failed: %s", strerror(errno));
    else if (got != sizeof(data) || memcmp(rbuf, data, got) != 0)
        FAIL(name, "data mismatch");
    else
        PASS(name);
}

static void test_credentials_file(void)
{
    char path[] = "/tmp/htslib_s3_creds_XXXXXX";
    char *id = strdup(getenv("AWS_ACCESS_KEY_ID"));
    char *secret = strdup(getenv("AWS_SECRET_ACCESS_KEY"));
    FILE *fp;
    int fd;

    if (!id || !secret || (fd = mkstemp(path)) < 0) {
        FAIL("credentials file", "setup failed");
        goto out;
    }
    fp = fdopen(fd, "w");
    fprintf(fp, "[default]\naws_access_key_id = %s\n"
            "aws_secret_access_key = %s\n", id, secret);
    fclose(fp);

    unsetenv("AWS_ACCESS_KEY_ID");
    unsetenv("AWS_SECRET_ACCESS_KEY");
    setenv("AWS_SHARED_CREDENTIALS_FILE", path, 1);
    test_roundtrip("credentials file", "creds_file.bin", 100, 100);
    unsetenv("AWS_SHARED_CREDENTIALS_FILE");
    setenv("AWS_ACCESS_KEY_ID", id, 1);
    setenv("AWS_SECRET_ACCESS_KEY", secret, 1);
    unlink(path);

 out:
    free(id);
    free(secret);
}

int main(void)
{
    const char *host = getenv("HTS_S3_HOST");
    const char *b = getenv("HTSLIB_TEST_S3_BUCKET");

    if (!host || !*host) {
        fprintf(stderr, "HTS_S3_HOST not set, skipping hfile_s3 tests\n");
        return EXIT_SUCCESS;
    }

    snprintf(bucket, sizeof(bucket), "%s", (b && *b) ? b : "htslib-test");
    snprintf(run_id, sizeof(run_id), "test-hfile-s3-%ld-%ld",
             (long) getpid(), (long) time(NULL));

    fprintf(stderr, "test_hfile_s3: s3+http://%s/%s on %s\n",
            bucket, run_id, host);

    test_roundtrip("single part", "single.bin", 5000, 5000);
    // Writes of 6, 6 and 2 MiB force a multipart upload (the minimum part
    // size is 5 MiB), and reading back takes several 1 MiB ranged GETs.
    // Sizes are kept off read part boundaries, as reading to EOF currently
    // fails for those (see the known failure in test/s3/s3.tst).
    test_roundtrip("multipart", "multipart.bin", 14 * MiB + 1000, 6 * MiB);
    test_roundtrip("URL-escaped key",
                   "weird key/has space+plus#hash%percent.bin", 100, 100);

    setenv("HTS_S3_PART_SIZE", "6", 1);
    setenv("HTS_S3_READ_PART_SIZE", "2", 1);
    test_roundtrip("part size overrides", "part_size.bin",
                   8 * MiB + 1000, 8 * MiB);
    unsetenv("HTS_S3_PART_SIZE");
    unsetenv("HTS_S3_READ_PART_SIZE");

    test_sigv2_read();
    test_credentials_file();
    test_missing();
    test_wrong_credentials();
    test_connection_refused();

    if (failures > 0) {
        fprintf(stderr, "%d test(s) FAILED\n", failures);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "All tests passed.\n");
    return EXIT_SUCCESS;
}
