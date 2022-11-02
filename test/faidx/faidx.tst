#    Copyright (C) 2022 Genome Research Ltd.
#
#    Author: Robert Davies <rmd@sanger.ac.uk>
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

# First field:
#   INIT = initialisation, not counted in testing
#   P = expected to pass (zero return; expected output matches, if present)
#   N = expected to return non-zero
#   F = expected to fail
#
# Second field (P/N/F only):
#   Filename of expected output.  If '.', output is not checked
#
# Rest:
#   Command to execute.  $bgzip and $test_faidx are replaced with the path to
# bgzip and test_faidx.

# Index fasta
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -e faidx.fa.expected.fai

# Test various functions on the fasta index
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t fai_line_length -e 24 trailingblank2
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai	-t faidx_has_seq -e 1 foo
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t faidx_has_seq -e 0 absent
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t faidx_iseq -e trailingblank3 3
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t faidx_seq_len -e 33 trailingblank1
P . $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t faidx_seq_len64 -e 72 trailingblank2

# Index fastq
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -e fastqs.fq.expected.fai

# Test various functions on the fastq index
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -Q -t fai_line_length -e 63 FAKE0005_3
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -Q -t fai_line_length -e 144 SRR014849.203935_3
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -t faidx_has_seq -e 1 SRR014849.203935_3
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -t faidx_has_seq -e 0 absent
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -t faidx_iseq -e FAKE0005_1 0
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -t faidx_seq_len -e 453 FSRRS4401CM938_1
P . $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -t faidx_seq_len64 -e 309 FSRRS4401AOV6A_4

# Fasta retrieval tests
P faidx.1.expected.fa $test_faidx -i faidx.fa -f faidx.fa.tmp.fai trailingblank2:28-33 trailingblank3:4-5 bar:4-5
P faidx.1.expected.fa $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t fai_fetch trailingblank2:28-33 trailingblank3:4-5 bar:4-5
P faidx.1.expected.fa $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t faidx_fetch_seq64 trailingblank2:28-33 trailingblank3:4-5 bar:4-5
P faidx.1.expected.fa $test_faidx -i faidx.fa -f faidx.fa.tmp.fai -t fai_adjust_region trailingblank2:28-33 trailingblank3:4-5 bar:4-5

# Fastq retrieval tests
P fastqs.1.expected.fq $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -Q FAKE0006_1:4-12 FSRRS4401BE7HA_1:81-120 FAKE0010_2 SRR014849.50939_3:71-90
P fastqs.1.expected.fq $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -Q -t fai_fetch FAKE0006_1:4-12 FSRRS4401BE7HA_1:81-120 FAKE0010_2 SRR014849.50939_3:71-90
P fastqs.1.expected.fq $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai -Q -t faidx_fetch_seq64 FAKE0006_1:4-12 FSRRS4401BE7HA_1:81-120 FAKE0010_2 SRR014849.50939_3:71-90
P fastqs.2.expected.fa $test_faidx -i fastqs.fq -f fastqs.fq.tmp.fai FAKE0006_1:4-12 FSRRS4401BE7HA_1:81-120 FAKE0010_2 SRR014849.50939_3:71-90

# Indexing and retrieval on bgzip compressed fasta
INIT $bgzip -c < ../ce.fa > ce.fa.tmp.gz
P . $test_faidx -i ce.fa.tmp.gz -f ce.fa.tmp.gz.fai -g ce.fa.tmp.gz.gzi -e ../ce.fa.fai
P ce.1.expected.fa $test_faidx -i ce.fa.tmp.gz -f ce.fa.tmp.gz.fai -g ce.fa.tmp.gz.gzi CHROMOSOME_I:5001-5125 CHROMOSOME_X:101-225
