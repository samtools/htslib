#    Copyright (C) 2025 Genome Research Ltd.
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
#   P = expected to pass (zero return; expected output matches, if present)
#   N = expected to return non-zero
#   F = expected to fail
#
# Second field (P/N/F only):
#   Filename of expected output.  If '.', output is not checked
#
# Rest:
#   Command to execute.
#   $test_view is replaced with the path to the test_view program

# Pairs with start and end differing.
P a7.sam  $test_view a7.cram
P a7b.sam $test_view a7b.cram
P a8.sam  $test_view a8.cram
P a8b.sam $test_view a8b.cram
P a9.sam  $test_view a9.cram
P a9b.sam $test_view a9b.cram

# Pairs with start matching and end differing
P b7.sam  $test_view b7.cram
P b7b.sam $test_view b7b.cram
P b8.sam  $test_view b8.cram
P b8b.sam $test_view b8b.cram

# Pairs with start differing and end matching
P c7.sam  $test_view c7.cram
P c7b.sam $test_view c7b.cram
P c8.sam  $test_view c8.cram
P c8b.sam $test_view c8b.cram

# Pairs with start and end both matching
P d7.sam  $test_view d7.cram
P d7b.sam $test_view d7b.cram

# Triplets with all start/ends matching
P d4.sam  $test_view d4.cram
P d4b.sam $test_view d4b.cram
P d4c.sam $test_view d4c.cram
P d4d.sam $test_view d4d.cram
P d4e.sam $test_view d4e.cram
P d4f.sam $test_view d4f.cram
P d5.sam  $test_view d5.cram
P d5b.sam $test_view d5b.cram
P d5c.sam $test_view d5c.cram
P d5d.sam $test_view d5d.cram
P d5e.sam $test_view d5e.cram
P d5f.sam $test_view d5f.cram

# Triplets with differing starts and ends
P a4.sam  $test_view a4.cram
P a5.sam  $test_view a5.cram
