#! /usr/bin/env perl
#
#    Copyright (C) 2013 Genome Research Ltd.
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
use strict;
use warnings;

my $err_count = 0;
my $suc_count = 0;

sub test {
    my ($cmd) = @_;
    print "  $cmd\n";
    if (system("$cmd || exit 1") != 0) {
        print "FAIL $!\n";
        $err_count++;
    } else {
        $suc_count++;
    }
}

foreach my $sam (glob("*#*.sam")) {
    my ($base, $ref) = ($sam =~ /((.*)#.*)\.sam/);
    $ref .= ".fa";

    my $bam  = "$base.tmp.bam";
    my $cram = "$base.tmp.cram";

    print "\n=== Testing $sam, ref $ref ===\n";

    # SAM -> BAM -> SAM
    test "./test_view -S -b $sam > $bam";
    test "./test_view $bam > $bam.sam_";
    test "./compare_sam.pl $sam $bam.sam_";

    # SAM -> CRAM -> SAM
    test "./test_view -t $ref -S -C $sam > $cram";
    test "./test_view -D $cram > $cram.sam_";
    test "./compare_sam.pl -nomd $sam $cram.sam_";

    # BAM -> CRAM -> BAM -> SAM
    $cram = "$bam.cram";
    test "./test_view -t $ref -C $bam > $cram";
    test "./test_view -b -D $cram > $cram.bam";
    test "./test_view $cram.bam > $cram.bam.sam_";
    test "./compare_sam.pl -nomd $sam $cram.bam.sam_";

    # SAM -> CRAM3 -> SAM
    $cram = "$base.tmp.cram";
    test "./test_view -t $ref -S -C -o VERSION=3.0 $sam > $cram";
    test "./test_view -D $cram > $cram.sam_";
    test "./compare_sam.pl -nomd $sam $cram.sam_";

    # BAM -> CRAM3 -> BAM -> SAM
    $cram = "$bam.cram";
    test "./test_view -t $ref -C -o VERSION=3.0 $bam > $cram";
    test "./test_view -b -D $cram > $cram.bam";
    test "./test_view $cram.bam > $cram.bam.sam_";
    test "./compare_sam.pl -nomd $sam $cram.bam.sam_";

    # CRAM3 -> CRAM2
    $cram = "$base.tmp.cram";
    test "./test_view -t $ref -C -o VERSION=2.1 $cram > $cram.cram";

    # CRAM2 -> CRAM3
    test "./test_view -t $ref -C -o VERSION=3.0 $cram.cram > $cram";
    test "./test_view $cram > $cram.sam_";
    test "./compare_sam.pl -nomd $sam $cram.sam_";
}

print "\nSuccesses $suc_count\n";
print "\nFailures  $err_count\n";

exit ($err_count > 0);
