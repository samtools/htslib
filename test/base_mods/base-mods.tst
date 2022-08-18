#    Copyright (C) 2020 Genome Research Ltd.
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

# Pileup testing
P MM-pileup.out $pileup_mod < MM-pileup.sam
P MM-pileup2.out $pileup_mod < MM-pileup2.sam
