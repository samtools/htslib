#!/bin/sh
#
#    Copyright (C) 2026 Peter Dowdy.
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

# Runs s3.tst against the S3-compatible server given by HTS_S3_HOST.
# Needs credentials in AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY, and an
# existing bucket named by HTSLIB_TEST_S3_BUCKET (default htslib-test).

# Load in the test driver
. ../simple_test_driver.sh

if [ -z "$HTS_S3_HOST" ]
then
    echo "HTS_S3_HOST not set, skipping S3 tests"
    exit 0
fi

echo "Testing S3..."

s3="s3+http://${HTSLIB_TEST_S3_BUCKET:-htslib-test}/test-s3-$$-`date +%s`"
bgzip="../../bgzip"
tabix="../../tabix"
htsfile="../../htsfile"
test_view="../test_view"
test_faidx="../test_faidx"

test_driver $@

exit $?
