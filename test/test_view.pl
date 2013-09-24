#! /usr/bin/env perl
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
    test "./compare_sam.pl $sam $cram.sam_";

    # BAM -> CRAM -> BAM -> SAM
    $cram = "$bam.cram";
    test "./test_view -t $ref -C $bam > $cram";
    test "./test_view -b -D $cram > $cram.bam";
    test "./test_view $cram.bam > $cram.bam.sam_";
    test "./compare_sam.pl $sam $cram.bam.sam_";
}

print "\nSuccesses $suc_count\n";
print "\nFailures  $err_count\n";

exit ($err_count > 0);
