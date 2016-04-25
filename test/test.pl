#!/usr/bin/env perl
#
#    Copyright (C) 2012-2013 Genome Research Ltd.
#
#    Author: Petr Danecek <pd3@sanger.ac.uk>
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
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use IO::Handle;

my $opts = parse_params();

test_view($opts);
test_vcf_api($opts,out=>'test-vcf-api.out');
test_vcf_sweep($opts,out=>'test-vcf-sweep.out');
test_vcf_various($opts);
test_convert_padded_header($opts);

print "\nNumber of tests:\n";
printf "    total   .. %d\n", $$opts{nok}+$$opts{nfailed};
printf "    passed  .. %d\n", $$opts{nok};
printf "    failed  .. %d\n", $$opts{nfailed};
print "\n";

exit ($$opts{nfailed} > 0);

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print
        "About: samtools/htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit 1;
}
sub parse_params
{
    my $opts = { keep_files=>0, nok=>0, nfailed=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            't|temp-dir:s' => \$$opts{keep_files},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    return $opts;
}
sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    }
    else
    {
        # child
        exec('bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed [$ret]: $cmd\n", $out); }
    return $out;
}
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out) = _cmd("$args{cmd} 2>&1");
    if ( $ret ) { failed($opts,$test); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
        if ( !$ret && $out eq '' ) { unlink("$$opts{path}/$args{out}.old"); }
        else
        {
            print "\tthe expected output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{out}.old\n";
            print "\t  new .. $$opts{path}/$args{out}\n";
        }
    }
    my $exp = '';
    if ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}: $!"); return; }

    if ( $exp ne $out )
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        if ( !-e "$$opts{path}/$args{out}" )
        {
            rename("$$opts{path}/$args{out}.new","$$opts{path}/$args{out}") or error("rename $$opts{path}/$args{out}.new $$opts{path}/$args{out}: $!");
            print "\tthe file with expected output does not exist, creating new one:\n";
            print "\t\t$$opts{path}/$args{out}\n";
        }
        else
        {
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new");
        }
        return;
    }
    passed($opts,$test);
}
sub failed
{
    my ($opts,$test,$reason) = @_;
    $$opts{nfailed}++;
    print "\n";
    STDOUT->flush();
    if ( defined $reason ) { print STDERR "\t$reason\n"; }
    print STDERR ".. failed ...\n\n";
    STDERR->flush();
}
sub passed
{
    my ($opts,$test) = @_;
    $$opts{nok}++;
    print ".. ok\n\n";
}
sub is_file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or return 0;
    my (@bstat) = stat($bfile) or return 0;
    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}


# The tests --------------------------

my $test_view_failures;
sub testv {
    my ($cmd) = @_;
    print "  $cmd\n";
    my ($ret, $out) = _cmd($cmd);
    if ($ret != 0) {
        STDOUT->flush();
        print STDERR "FAILED\n\n";
        STDERR->flush();
        $test_view_failures++;
    }
}

sub test_view
{
    my ($opts, %args) = @_;

    foreach my $sam (glob("*#*.sam")) {
        my ($base, $ref) = ($sam =~ /((.*)#.*)\.sam/);
        $ref .= ".fa";

        my $bam  = "$base.tmp.bam";
        my $cram = "$base.tmp.cram";

        print "test_view testing $sam, ref $ref:\n";
        $test_view_failures = 0;

        # SAM -> BAM -> SAM
        testv "./test_view -S -b $sam > $bam";
        testv "./test_view $bam > $bam.sam_";
        testv "./compare_sam.pl $sam $bam.sam_";

        # SAM -> CRAM -> SAM
        testv "./test_view -t $ref -S -C $sam > $cram";
        testv "./test_view -D $cram > $cram.sam_";
        testv "./compare_sam.pl -nomd $sam $cram.sam_";

        # BAM -> CRAM -> BAM -> SAM
        $cram = "$bam.cram";
        testv "./test_view -t $ref -C $bam > $cram";
        testv "./test_view -b -D $cram > $cram.bam";
        testv "./test_view $cram.bam > $cram.bam.sam_";
        testv "./compare_sam.pl -nomd $sam $cram.bam.sam_";

        # SAM -> CRAM3 -> SAM
        $cram = "$base.tmp.cram";
        testv "./test_view -t $ref -S -C -o VERSION=3.0 $sam > $cram";
        testv "./test_view -D $cram > $cram.sam_";
        testv "./compare_sam.pl -nomd $sam $cram.sam_";

        # BAM -> CRAM3 -> BAM -> SAM
        $cram = "$bam.cram";
        testv "./test_view -t $ref -C -o VERSION=3.0 $bam > $cram";
        testv "./test_view -b -D $cram > $cram.bam";
        testv "./test_view $cram.bam > $cram.bam.sam_";
        testv "./compare_sam.pl -nomd $sam $cram.bam.sam_";

        # CRAM3 -> CRAM2
        $cram = "$base.tmp.cram";
        testv "./test_view -t $ref -C -o VERSION=2.1 $cram > $cram.cram";

        # CRAM2 -> CRAM3
        testv "./test_view -t $ref -C -o VERSION=3.0 $cram.cram > $cram";
        testv "./test_view $cram > $cram.sam_";
        testv "./compare_sam.pl -nomd $sam $cram.sam_";

        # Java pre-made CRAM -> SAM
        my $jcram = "${base}_java.cram";
        if (-e $jcram) {
            my $jsam = "${base}_java.tmp.sam_";
            testv "./test_view -i reference=$ref $jcram > $jsam";
            testv "./compare_sam.pl -nomd $sam $jsam";
        }

        if ($test_view_failures == 0)
        {
            passed($opts, "$sam conversions");
        }
        else
        {
            failed($opts, "$sam conversions", "$test_view_failures subtests failed");
        }
    }
}

sub test_vcf_api
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{path}/test-vcf-api $$opts{tmp}/test-vcf-api.bcf");
}

sub test_vcf_sweep
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{path}/test-vcf-sweep $$opts{tmp}/test-vcf-api.bcf");
}

sub test_vcf_various
{
    my ($opts, %args) = @_;

    # Excess spaces in header lines
    test_cmd($opts, %args, out => "test-vcf-hdr.out",
        cmd => "$$opts{bin}/htsfile -ch $$opts{path}/test-vcf-hdr-in.vcf");

    # Various VCF parsing issues
    test_cmd($opts, %args, out => "formatcols.vcf",
        cmd => "$$opts{bin}/htsfile -c $$opts{path}/formatcols.vcf");
    test_cmd($opts, %args, out => "noroundtrip-out.vcf",
        cmd => "$$opts{bin}/htsfile -c $$opts{path}/noroundtrip.vcf");
}

sub test_convert_padded_header
{
    my ($opts, %args) = @_;

    $args{out} = "headernul.tmp.cram";
    cmd("$$opts{path}/test_view -t ce.fa -C ce#1.sam > $args{out}");

    foreach my $nuls (0, 1, 678) {
        my $nulsbam = "$$opts{tmp}/headernul$nuls.bam";
        cmd("$$opts{path}/test_view -b -Z $nuls ce#1.sam > $nulsbam");
        test_cmd($opts, %args,
            cmd => "$$opts{path}/test_view -t ce.fa -C $nulsbam");
    }
}
