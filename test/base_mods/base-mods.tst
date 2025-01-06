#    Copyright (C) 2020, 2023 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
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
#   P = expected to pass
#   N = expected to return non-zero
#   F = expected to fail

# Second field:
#   Filename of expected output

# Third onwards; command to execute. $fmt is replaced by the current file
# format, ie sam, bam or cram. $samtools is a pointer to the desired
# samtools binary. This can be useful for testing older versions.

# Test files from SAM spec
P MM-chebi.out       $test_mod    MM-chebi.sam
P MM-double.out      $test_mod    MM-double.sam
P MM-multi.out       $test_mod    MM-multi.sam
P MM-explicit.out    $test_mod    MM-explicit.sam
P MM-explicit-x.out  $test_mod -x MM-explicit.sam

# Report bases outside the explicitly called ranges, so we could exclude
# these in any depth based consensus analysis and only gather statistics 
# for sites known to be have been scanned.
P MM-explicit-f.out  $test_mod -f 1 MM-explicit.sam

# Ensure state gets reset correctly between reads
P MM-not-all-modded.out	$test_mod MM-not-all-modded.sam

# Pileup testing
P MM-pileup.out $pileup_mod < MM-pileup.sam
P MM-pileup2.out $pileup_mod < MM-pileup2.sam

# Validation testing.  We just care about exit status here, but the
# test data is a copy of MM-pileup.sam so that suffices too.
P MM-pileup.out $pileup_mod < MM-MNp.sam
N MM-pileup.out $pileup_mod < MM-MNf1.sam
N MM-pileup.out $pileup_mod < MM-MNf2.sam
N MM-pileup.out $test_mod < MM-MNf1.sam
N MM-pileup.out $test_mod < MM-MNf2.sam
N MM-pileup.out $test_mod < MM-bounds+.sam
N MM-pileup.out $test_mod < MM-bounds-.sam
