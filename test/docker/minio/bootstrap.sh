#!/bin/sh
# test/docker/minio/bootstrap.sh -- entrypoint for the htslib S3 test MinIO image.
#
# Copyright (C) 2026 Broad Institute.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
# Starts MinIO in the background, waits for it to come up, provisions the
# test bucket, then waits on the server process so it remains PID 1.
#
# NOTE for callers: MinIO's own /minio/health/ready endpoint goes green as
# soon as the server itself is up, which is *before* the `mc mb` below has
# created the test bucket -- polling only /minio/health/ready races this
# script and can see "bucket does not exist" errors. Wait for the
# "bootstrap: bucket ... ready" line below (e.g. via `docker logs`) too, as
# .github/workflows/s3-integration.yml does.
set -eu

: "${MINIO_ROOT_USER:=minioadmin}"
: "${MINIO_ROOT_PASSWORD:=minioadmin}"
: "${MINIO_TEST_BUCKET:=htslib-test}"

minio server /data --console-address :9001 &
server_pid=$!

tries=0
until mc alias set local http://127.0.0.1:9000 \
        "$MINIO_ROOT_USER" "$MINIO_ROOT_PASSWORD" >/dev/null 2>&1; do
    tries=$((tries + 1))
    if [ "$tries" -ge 60 ]; then
        echo "bootstrap: minio did not become ready in time" >&2
        exit 1
    fi
    sleep 1
done

mc mb -p "local/$MINIO_TEST_BUCKET"

echo "bootstrap: bucket '$MINIO_TEST_BUCKET' ready"

wait "$server_pid"
