#    Copyright (C) 2020, 2022 Genome Research Ltd.
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
#   Command to execute.  $tv is replaced with the path to test_view

# String matches
P string1.out $tv -i 'filter=qname =~ "\.1" && cigar =~ "D"' ../ce#1000.sam
P string2.out $tv -i 'filter=rname=="CHROMOSOME_II"' ../ce#5b.sam
P string3.out $tv -i 'filter=rname=~"CHROMOSOME_II"' ../ce#5b.sam
P string4.out $tv -i 'filter=cigar=~"D"' ../ce#1000.sam
P string5.out $tv -i 'filter=seq =~ "(AT){2}"' ../ce#1000.sam
P string6.out $tv -i 'filter=library=="x"' ../xx#rg.sam
P string7.out $tv -i 'filter=library!="x"' ../xx#rg.sam

# Integer ops
P int1.out    $tv -i 'filter=pos % 23 == 11' ../ce#1000.sam | grep -E -cv '^@'
P int2.out    $tv -i 'filter=qlen/(flag*mapq+pos)>5' ../ce#1000.sam | grep -E -cv '^@'

# Aux tags
P int3.out    $tv -i 'filter=[NM]>=10 || [MD]=~"A.*A.*A"' -t4 ../ce#1000.sam | grep -E -cv '^@'

# Functions.
P func1.out   $tv -i 'filter=length(seq) != qlen' ../ce#5b.sam | grep -E -cv '^@'
P func2.out   $tv -i 'filter=min(qual) >= 20' ../ce#1000.sam | grep -E -cv '^@'
P func3.out   $tv -i 'filter=max(qual) <= 20' ../ce#1000.sam | grep -E -cv '^@'
P func4.out   $tv -i 'filter=avg(qual) >= 20 && avg(qual) <= 30' ../ce#1000.sam | grep -E -cv '^@'
P func5.out   $tv -i 'filter=sclen>=20' ../realn02.sam | grep -E -v '^@'
P func6.out   $tv -i 'filter=rlen<50'   ../realn02.sam | grep -E -v '^@'
P func7.out   $tv -i 'filter=qlen>100'   ../realn02.sam | grep -E -v '^@'
P func8.out   $tv -i 'filter=hclen>=4'   ../c1#clip.sam | grep -E -v '^@'
