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

# Single file, unpaired data, with / without aux tags
P single_noaux.sam $tview single.fq
P single_aux.sam $tview -i fastq_aux single.fq

# Single file, interleaved paired data, no aux
P inter_noaux.sam $tview interleaved.fq

# Single file, interleaved paired data, with aux
P inter_aux.sam $tview -i fastq_aux interleaved.fq

# Single file, interleaved paired data, using CASAVA
P inter_casava.sam $tview -i fastq_casava interleaved_casava.fq

# CASAVA with filtering
P filter_casava.sam $tview -i fastq_casava filter_casava.fq


# Paired data is mainly tested by the Samtools test harness.
# Basically though it's just reading two files and relying on either
# this code or explicit overloading of READ1/READ2.
# We simply test here we can read r1 and r2 as separate files
P r1.sam $tview -i fastq_aux r1.fq
P r2.sam $tview -i fastq_aux r2.fq
