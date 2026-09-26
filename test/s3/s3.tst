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
#   Command to execute.  $s3 is the S3 URL prefix for this run, and $bgzip,
# $tabix, $htsfile, $test_view and $test_faidx are the paths to those programs.
#
# Expected outputs not already in test/ are made by running the same command
# on a local copy, so each test checks that going via S3 changes nothing.
# Indexes loaded for a remote file are cached in the current directory, so
# S3 keys are named *.tmp.* and those caches removed first.

INIT rm -f *.tmp.*

# Raw copy up and back down
P . $htsfile -C ../bgziptest.txt $s3/copy.tmp.txt
P ../bgziptest.txt $htsfile -C $s3/copy.tmp.txt -
N . $htsfile $s3/does-not-exist.tmp.bam

# Known bug: reading to EOF fails with EINVAL when the size is a multiple of
# the read part size (1 MiB by default), as hfile_s3 requests a range starting
# at EOF and gets a 416.  Change F to P once fixed.
INIT head -c 1048576 /dev/zero > aligned.tmp.bin
P . $htsfile -C aligned.tmp.bin $s3/aligned.tmp.bin
F . $htsfile -C $s3/aligned.tmp.bin aligned-copy.tmp.bin

# BGZF with a .gzi index, random access through the remote index
INIT seq 1 20000 > seq.tmp.txt
INIT $bgzip -i -o seq.tmp.txt.gz seq.tmp.txt
INIT $bgzip -b 71088 -s 60 seq.tmp.txt.gz > seq.tmp.out
P . $bgzip -i -o $s3/seq.tmp.txt.gz seq.tmp.txt
P seq.tmp.out $bgzip -b 71088 -s 60 $s3/seq.tmp.txt.gz

# BAM with CSI and BAI indexes, threaded read and write
INIT $test_view -b -m 14 -x local.tmp.bam.csi -p local.tmp.bam ../index.sam
INIT $test_view local.tmp.bam CHROMOSOME_II:2960-2970 > bam_region.tmp.out
INIT $test_view local.tmp.bam > bam_all.tmp.out
INIT $htsfile local.tmp.bam | cut -f 2 > htsfile.tmp.out
P . $test_view -b -m 14 -x $s3/csi.tmp.bam.csi -p $s3/csi.tmp.bam ../index.sam
P bam_region.tmp.out $test_view $s3/csi.tmp.bam CHROMOSOME_II:2960-2970
P htsfile.tmp.out $htsfile $s3/csi.tmp.bam | cut -f 2
P . $test_view -b -m 0 -x $s3/bai.tmp.bam.bai -p $s3/bai.tmp.bam ../index.sam
P bam_region.tmp.out $test_view $s3/bai.tmp.bam CHROMOSOME_II:2960-2970
P . $test_view -@ 4 -b -p $s3/threads.tmp.bam ../index.sam
P bam_all.tmp.out $test_view -@ 4 $s3/threads.tmp.bam

# CRAM with the reference and .crai index both on S3 (headers are skipped
# as @SQ UR: records where the reference was read from)
INIT $test_view -C -t ../ce.fa -x local.tmp.cram.crai -p local.tmp.cram ../index.sam
INIT $test_view -i reference=../ce.fa local.tmp.cram CHROMOSOME_II:2960-2970 | grep -v "^@" > cram_region.tmp.out
P . $htsfile -C ../ce.fa $s3/ce.tmp.fa
P . $htsfile -C ../ce.fa.fai $s3/ce.tmp.fa.fai
P . $test_view -C -t $s3/ce.tmp.fa -x $s3/ref.tmp.cram.crai -p $s3/ref.tmp.cram ../index.sam
P cram_region.tmp.out $test_view -i reference=$s3/ce.tmp.fa $s3/ref.tmp.cram CHROMOSOME_II:2960-2970 | grep -v "^@"

# VCF, BCF and BED with tabix indexes built and queried remotely
P . $bgzip -o $s3/vcf.tmp.vcf.gz ../tabix/vcf_file.vcf
P . $tabix -p vcf $s3/vcf.tmp.vcf.gz
P ../tabix/vcf_file.1.3000151.out $tabix -D $s3/vcf.tmp.vcf.gz 1:3000151-3000151
P . $test_view -b -m 14 -x $s3/bcf.tmp.bcf.csi -p $s3/bcf.tmp.bcf ../tabix/vcf_file.vcf
P ../tabix/vcf_file.2.3199812.out $tabix -D $s3/bcf.tmp.bcf 2:3199812-3199812
P . $bgzip -o $s3/bed.tmp.bed.gz ../tabix/bed_file.bed
P . $tabix -p bed $s3/bed.tmp.bed.gz
P ../tabix/bed_file.Y.100200.out $tabix -D $s3/bed.tmp.bed.gz Y:100200-100200

# Bgzipped FASTA: .fai built on S3, then loaded with the .gzi
INIT $test_faidx -c -i ../faidx/faidx.fa -f local.tmp.fai trailingblank1 > faidx.tmp.out
P . $bgzip -i -o $s3/faidx.tmp.fa.gz ../faidx/faidx.fa
P faidx.tmp.out $test_faidx -c -i $s3/faidx.tmp.fa.gz trailingblank1
P faidx.tmp.out $test_faidx -i $s3/faidx.tmp.fa.gz trailingblank1
