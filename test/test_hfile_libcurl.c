/*  test/test_hfile_libcurl.c -- Test cases for libcurl retry/resilience.

    Copyright (C) 2025 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>

#include "../htslib/hfile.h"
#include "../htslib/hts.h"
#include "../hts_internal.h"

// Test data file to serve
#define TEST_DATA_FILE "hfile_libcurl.tmp"
#define TEST_DATA_SIZE 16384

static int failures = 0;

#define PASS(name) fprintf(stderr, "  PASS: %s\n", (name))
#define FAIL(name, ...) do { \
    fprintf(stderr, "  FAIL: %s: ", (name)); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    failures++; \
} while (0)

// Start mock server, return its PID and port
static pid_t start_server(const char *mode, int fail_count, int *port_out)
{
    int pipefd[2];
    pid_t pid;
    char fail_count_str[32];
    char port_buf[32];
    ssize_t n;
    int i;

    if (pipe(pipefd) < 0) {
        perror("pipe");
        return -1;
    }

    snprintf(fail_count_str, sizeof(fail_count_str), "%d", fail_count);

    pid = fork();
    if (pid < 0) {
        perror("fork");
        close(pipefd[0]);
        close(pipefd[1]);
        return -1;
    }

    if (pid == 0) {
        // Child: redirect stdout to pipe, exec python server
        close(pipefd[0]);
        dup2(pipefd[1], STDOUT_FILENO);
        close(pipefd[1]);
        execlp("python3", "python3", "test/mock_http_server.py",
               "--mode", mode,
               "--file", "test/" TEST_DATA_FILE,
               "--fail-count", fail_count_str,
               "--port", "0",
               NULL);
        perror("execlp python3");
        _exit(127);
    }

    // Parent: read port from pipe
    close(pipefd[1]);

    // Read with timeout (give server up to 5 seconds to start)
    n = 0;
    for (i = 0; i < 50 && n == 0; i++) {
        n = read(pipefd[0], port_buf, sizeof(port_buf) - 1);
        if (n <= 0) {
            hts_usleep(100000); // 100ms
            n = 0;
        }
    }
    close(pipefd[0]);

    if (n <= 0) {
        fprintf(stderr, "Failed to read port from mock server\n");
        kill(pid, SIGTERM);
        waitpid(pid, NULL, 0);
        return -1;
    }

    port_buf[n] = '\0';
    *port_out = atoi(port_buf);
    if (*port_out <= 0) {
        fprintf(stderr, "Invalid port from mock server: %s\n", port_buf);
        kill(pid, SIGTERM);
        waitpid(pid, NULL, 0);
        return -1;
    }

    // Give the server a moment to be ready for connections
    hts_usleep(100000);
    return pid;
}

static void stop_server(pid_t pid)
{
    int status;
    pid_t ret;

    if (pid <= 0)
        return;

    kill(pid, SIGKILL);
    ret = waitpid(pid, &status, 0);
    if (ret < 0)
        perror("waitpid");
}

// Generate deterministic test data
static void generate_test_data(void)
{
    FILE *f;
    int i;
    char path[256];

    snprintf(path, sizeof(path), "test/%s", TEST_DATA_FILE);
    f = fopen(path, "wb");
    if (!f) {
        perror("fopen test data");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < TEST_DATA_SIZE; i++) {
        fputc((i * 7 + 13) & 0xFF, f);
    }
    fclose(f);
}

// Read expected data from test file
static unsigned char *read_expected_data(size_t *size_out)
{
    FILE *f;
    unsigned char *data;
    char path[256];

    snprintf(path, sizeof(path), "test/%s", TEST_DATA_FILE);
    f = fopen(path, "rb");
    if (!f) {
        perror("fopen expected data");
        return NULL;
    }
    data = malloc(TEST_DATA_SIZE);
    if (!data) {
        fclose(f);
        return NULL;
    }
    *size_out = fread(data, 1, TEST_DATA_SIZE, f);
    fclose(f);
    return data;
}

// Read entire file via hFILE
static unsigned char *hfile_read_all(const char *url, size_t *size_out)
{
    hFILE *fp;
    unsigned char *buf;
    size_t total = 0;
    ssize_t n;

    fp = hopen(url, "r");
    if (!fp)
        return NULL;

    buf = malloc(TEST_DATA_SIZE + 1024);
    if (!buf) {
        (void) hclose(fp);
        return NULL;
    }

    while ((n = hread(fp, buf + total, 4096)) > 0) {
        total += n;
        if (total > TEST_DATA_SIZE + 512)
            break;
    }

    if (n < 0) {
        free(buf);
        (void) hclose(fp);
        *size_out = 0;
        return NULL;
    }

    (void) hclose(fp);
    *size_out = total;
    return buf;
}

// Test 1: Normal transfer
static void test_normal_transfer(void)
{
    const char *name = "Normal transfer";
    pid_t pid;
    int port;
    char url[256];
    unsigned char *got, *expected;
    size_t got_size, exp_size;

    pid = start_server("normal", 0, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    got = hfile_read_all(url, &got_size);
    expected = read_expected_data(&exp_size);

    if (!got) {
        FAIL(name, "hfile_read_all returned NULL, errno=%d", errno);
    } else if (got_size != exp_size) {
        FAIL(name, "size mismatch: got %zu, expected %zu", got_size, exp_size);
    } else if (memcmp(got, expected, exp_size) != 0) {
        FAIL(name, "data mismatch");
    } else {
        PASS(name);
    }

    free(got);
    free(expected);
    stop_server(pid);
}

// Test 2: 503 retry succeeds
static void test_503_retry(void)
{
    const char *name = "503 retry succeeds";
    pid_t pid;
    int port;
    char url[256];
    unsigned char *got, *expected;
    size_t got_size, exp_size;

    pid = start_server("503_then_ok", 2, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "3", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    got = hfile_read_all(url, &got_size);
    expected = read_expected_data(&exp_size);

    if (!got) {
        FAIL(name, "hfile_read_all returned NULL, errno=%d", errno);
    } else if (got_size != exp_size) {
        FAIL(name, "size mismatch: got %zu, expected %zu", got_size, exp_size);
    } else if (memcmp(got, expected, exp_size) != 0) {
        FAIL(name, "data mismatch");
    } else {
        PASS(name);
    }

    free(got);
    free(expected);
    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

// Test 3: 429 retry succeeds
static void test_429_retry(void)
{
    const char *name = "429 retry succeeds";
    pid_t pid;
    int port;
    char url[256];
    unsigned char *got, *expected;
    size_t got_size, exp_size;

    pid = start_server("429_then_ok", 2, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "3", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    got = hfile_read_all(url, &got_size);
    expected = read_expected_data(&exp_size);

    if (!got) {
        FAIL(name, "hfile_read_all returned NULL, errno=%d", errno);
    } else if (got_size != exp_size) {
        FAIL(name, "size mismatch: got %zu, expected %zu", got_size, exp_size);
    } else if (memcmp(got, expected, exp_size) != 0) {
        FAIL(name, "data mismatch");
    } else {
        PASS(name);
    }

    free(got);
    free(expected);
    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

// Test 4: Connection drop retry
static void test_drop_mid_transfer(void)
{
    const char *name = "Connection drop retry";
    pid_t pid;
    int port;
    char url[256];
    unsigned char *got, *expected;
    size_t got_size, exp_size;

    pid = start_server("drop_mid_transfer", 2, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "5", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    got = hfile_read_all(url, &got_size);
    expected = read_expected_data(&exp_size);

    if (!got) {
        FAIL(name, "hfile_read_all returned NULL, errno=%d", errno);
    } else if (got_size != exp_size) {
        FAIL(name, "size mismatch: got %zu, expected %zu", got_size, exp_size);
    } else if (memcmp(got, expected, exp_size) != 0) {
        FAIL(name, "data mismatch");
    } else {
        PASS(name);
    }

    free(got);
    free(expected);
    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

// Test 5: 404 not retried
static void test_404_no_retry(void)
{
    const char *name = "404 not retried";
    pid_t pid;
    int port;
    char url[256];
    hFILE *fp;

    pid = start_server("404", 0, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "3", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    fp = hopen(url, "r");
    if (fp != NULL) {
        FAIL(name, "hopen should have failed for 404");
        (void) hclose(fp);
    } else if (errno != ENOENT) {
        FAIL(name, "expected ENOENT, got errno=%d (%s)", errno, strerror(errno));
    } else {
        PASS(name);
    }

    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

// Test 6: Retry exhaustion
static void test_retry_exhaustion(void)
{
    const char *name = "Retry exhaustion";
    pid_t pid;
    int port;
    char url[256];
    hFILE *fp;

    pid = start_server("503_then_ok", 999, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "2", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    fp = hopen(url, "r");
    if (fp != NULL) {
        FAIL(name, "hopen should have failed after retry exhaustion");
        (void) hclose(fp);
    } else {
        PASS(name);
    }

    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

// Test 7: Retry disabled
static void test_retry_disabled(void)
{
    const char *name = "Retry disabled";
    pid_t pid;
    int port;
    char url[256];
    hFILE *fp;

    pid = start_server("503_then_ok", 1, &port);
    if (pid < 0) { FAIL(name, "could not start server"); return; }

    snprintf(url, sizeof(url), "http://127.0.0.1:%d/data", port);
    setenv("HTS_RETRY_MAX", "0", 1);
    setenv("HTS_RETRY_DELAY", "50", 1);

    fp = hopen(url, "r");
    if (fp != NULL) {
        FAIL(name, "hopen should have failed with retries disabled");
        (void) hclose(fp);
    } else {
        PASS(name);
    }

    unsetenv("HTS_RETRY_MAX");
    unsetenv("HTS_RETRY_DELAY");
    stop_server(pid);
}

int main(void)
{
    // Check python3 is available
    if (system("python3 --version >/dev/null 2>&1") != 0) {
        fprintf(stderr, "python3 not found, skipping libcurl retry tests\n");
        return 0; // Skip rather than fail
    }

    // Check that HTTP URLs are supported (i.e. libcurl backend is loaded).
    // Try opening a URL that will fail to connect — if errno is ENOTSUP
    // or similar, the http:// scheme isn't registered.
    {
        hFILE *probe;
        setenv("HTS_RETRY_MAX", "0", 1);
        probe = hopen("http://0.0.0.0:1/probe", "r");
        unsetenv("HTS_RETRY_MAX");
        if (probe) {
            (void) hclose(probe);
        } else if (errno == ENOTSUP) {
            fprintf(stderr, "HTTP not supported, skipping libcurl retry tests\n");
            return 0;
        }
        // ECONNREFUSED or similar means libcurl is working, continue
    }

    generate_test_data();

    fprintf(stderr, "test_hfile_libcurl:\n");

    test_normal_transfer();
    test_503_retry();
    test_429_retry();
    test_drop_mid_transfer();
    test_404_no_retry();
    test_retry_exhaustion();
    test_retry_disabled();

    // Clean up test data
    unlink("test/" TEST_DATA_FILE);

    if (failures > 0) {
        fprintf(stderr, "%d test(s) FAILED\n", failures);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "All tests passed.\n");
    return EXIT_SUCCESS;
}
