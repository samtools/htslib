#    Copyright (C) 2017-2018 Genome Research Ltd.
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
#   Command to execute.  $pileup is replaced with the path to the pileup test
# program

# Deletions
P mp_D.out $pileup mp_D.sam
P mp_D.out $pileup -m mp_D.sam

# Deletions followed by insertions
P mp_DI.out $pileup mp_DI.sam
P mp_DI.out $pileup -m mp_DI.sam

# NB: pileup currently cannot return leading insertions.
# Test output reflects this.
# Insertions
P mp_I.out $pileup mp_I.sam
P mp_I.out $pileup -m mp_I.sam
P mp_P.out $pileup mp_P.sam
P mp_P.out $pileup -m mp_P.sam

# Insertions followed by deletions
P mp_ID.out $pileup mp_ID.sam
P mp_ID.out $pileup -m mp_ID.sam

# Ref skips
P mp_N.out $pileup mp_N.sam
P mp_N.out $pileup -m mp_N.sam

# Ref skips and deletions
P mp_N2.out $pileup mp_N2.sam
P mp_N2.out $pileup -m mp_N2.sam

# Various combinations of insertions, deletions and pads
P c1#pad1.out $pileup c1#pad1.sam
P c1#pad1.out $pileup -m c1#pad1.sam
P c1#pad2.out $pileup c1#pad2.sam
P c1#pad2.out $pileup -m c1#pad2.sam
P c1#pad3.out $pileup c1#pad3.sam
P c1#pad3.out $pileup -m c1#pad3.sam

# Issue #852.  Problem caused by alignments with entirely S/I ops in CIGAR.
P small.out $pileup -m small.bam

# Overlap removal and the effect on quality values
P mp_overlap1.out $pileup -m mp_overlap1.sam
P mp_overlap2.out $pileup -m mp_overlap2.sam

