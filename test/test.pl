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

test_view($opts,0);
test_view($opts,4);

test_vcf_api($opts,out=>'test-vcf-api.out');
test_vcf_sweep($opts,out=>'test-vcf-sweep.out');
test_vcf_various($opts);
test_bcf_sr_sort($opts);
test_command($opts,cmd=>'test-bcf-translate -',out=>'test-bcf-translate.out');
test_convert_padded_header($opts);
test_rebgzip($opts);
test_logging($opts);
test_realn($opts);

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
        "   -f, --fail-fast                 Fail-fast mode: exit as soon as a test fails.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit 1;
}

sub cygpath {
    my ($path) = @_;
    $path = `cygpath -m $path`;
    $path =~ s/\r?\n//;
    return $path
}

sub safe_tempdir
{
    my $dir = tempdir(CLEANUP=>1);
    if ($^O =~ /^msys/) {
        $dir = cygpath($dir);
    }
    return $dir;
}

sub parse_params
{
    my $opts = { keep_files=>0, nok=>0, nfailed=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            't|temp-dir:s' => \$$opts{keep_files},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'f|fail-fast' => \$$opts{fail_fast},
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : safe_tempdir();
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    if ($^O =~ /^msys/) {
	$$opts{path} = cygpath($$opts{path});
	$$opts{bin}  = cygpath($$opts{bin});
    }

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
	# Example of how to embed Valgrind into the testing framework.
	# TEST_PRECMD="valgrind --leak-check=full --suppressions=$ENV{HOME}/valgrind.supp" make check
	$cmd = "$ENV{TEST_PRECMD} $cmd" if exists $ENV{TEST_PRECMD};

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

    my ($ret,$out) = _cmd("$args{cmd}");
    if ( $ret ) { failed($opts,$test); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("cmp $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
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
        $exp =~ s/\015?\012/\n/g;
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}: $!"); return; }

    (my $out_lf = $out) =~ s/\015?\012/\n/g;
    if ( $exp ne $out_lf )
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
    if ($$opts{fail_fast}) {
      die "\n";
    }
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
    my ($opts, $cmd) = @_;
    print "  $cmd\n";
    my ($ret, $out) = _cmd($cmd);
    if ($ret != 0) {
        STDOUT->flush();
        print STDERR "FAILED\n$out\n";
        STDERR->flush();
        $test_view_failures++;
        if ($$opts{fail_fast}) {
          die "\n";
        }
    }
}

sub test_view
{
    my ($opts, $nthreads) = @_;
    my $tv_args = $nthreads ? "-\@$nthreads" : "";

    foreach my $sam (glob("*#*.sam")) {
        my ($base, $ref) = ($sam =~ /((.*)#.*)\.sam/);
        $ref .= ".fa";

        my $bam  = "$base.tmp.bam";
        my $cram = "$base.tmp.cram";

        my $md = "-nomd";
        if ($sam =~ /^md/) {
            $md = "";
        }

        print "test_view testing $sam, ref $ref:\n";
        $test_view_failures = 0;

        # SAM -> BAM -> SAM
        testv $opts, "./test_view $tv_args -S -b $sam > $bam";
        testv $opts, "./test_view $tv_args $bam > $bam.sam_";
        testv $opts, "./compare_sam.pl $sam $bam.sam_";

        # SAM -> BAMu -> SAM
        testv $opts, "./test_view $tv_args -S -l0 -b $sam > $bam";
        testv $opts, "./test_view $tv_args $bam > $bam.sam_";
        testv $opts, "./compare_sam.pl $sam $bam.sam_";

        # SAM -> CRAM2 -> SAM
        testv $opts, "./test_view $tv_args -t $ref -S -C -o VERSION=2.1 $sam > $cram";
        testv $opts, "./test_view $tv_args -D $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl $md $sam $cram.sam_";

        # BAM -> CRAM2 -> BAM -> SAM
        $cram = "$bam.cram";
        testv $opts, "./test_view $tv_args -t $ref -C -o VERSION=2.1 $bam > $cram";
        testv $opts, "./test_view $tv_args -b -D $cram > $cram.bam";
        testv $opts, "./test_view $tv_args $cram.bam > $cram.bam.sam_";
        testv $opts, "./compare_sam.pl $md $sam $cram.bam.sam_";

        # SAM -> CRAM3u -> SAM
        $cram = "$base.tmp.cram";
        testv $opts, "./test_view $tv_args -t $ref -S -l0 -C -o VERSION=3.0 $sam > $cram";
        testv $opts, "./test_view $tv_args -D $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl $md $sam $cram.sam_";

        # BAM -> CRAM3 -> BAM -> SAM
        $cram = "$bam.cram";
        testv $opts, "./test_view $tv_args -t $ref -C -o VERSION=3.0 $bam > $cram";
        testv $opts, "./test_view $tv_args -b -D $cram > $cram.bam";
        testv $opts, "./test_view $tv_args $cram.bam > $cram.bam.sam_";
        testv $opts, "./compare_sam.pl $md $sam $cram.bam.sam_";

        # CRAM3 -> CRAM2
        $cram = "$base.tmp.cram";
        testv $opts, "./test_view $tv_args -t $ref -C -o VERSION=2.1 $cram > $cram.cram";

        # CRAM2 -> CRAM3
        testv $opts, "./test_view $tv_args -t $ref -C -o VERSION=3.0 $cram.cram > $cram";

	# CRAM3 -> CRAM3 + multi-slice
	testv $opts, "./test_view $tv_args -t $ref -C -o VERSION=3.0 -o seqs_per_slice=7 -o slices_per_container=5 $cram.cram > $cram";
        testv $opts, "./test_view $tv_args $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl $md $sam $cram.sam_";

        # Java pre-made CRAM -> SAM
        my $jcram = "${base}_java.cram";
        if (-e $jcram) {
            my $jsam = "${base}_java.tmp.sam_";
            testv $opts, "./test_view $tv_args -i reference=$ref $jcram > $jsam";
            testv $opts, "./compare_sam.pl -Baux $md $sam $jsam";
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

    # BAM and CRAM range queries on prebuilt BAM and CRAM
    # The cram file has @SQ UR: set to point to an invalid location to
    # force the reference to be reloaded from the one given on the
    # command line and nowhere else.  REF_PATH should also point to nowhere
    # (currently done by the Makefile).  This is to test the refseq reference
    # counting and reload (Issue #654).
    my $regions = "CHROMOSOME_II:2980-2980 CHROMOSOME_IV:1500-1500 CHROMOSOME_II:2980-2980 CHROMOSOME_I:1000-1100";
    testv $opts, "./test_view $tv_args -i reference=ce.fa range.cram $regions > range.tmp";
    testv $opts, "./compare_sam.pl range.tmp range.out";

    testv $opts, "./test_view $tv_args range.bam $regions > range.tmp";
    testv $opts, "./compare_sam.pl range.tmp range.out";
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
    test_cmd($opts, %args, out => "formatmissing-out.vcf",
        cmd => "$$opts{bin}/htsfile -c $$opts{path}/formatmissing.vcf");
}

sub write_multiblock_bgzf {
    my ($name, $frags) = @_;

    my $tmp = "$name.tmp";
    open(my $out, '>', $name) || die "Couldn't open $name $!\n";
    for (my $i = 0; $i < @$frags; $i++) {
	local $/;
	open(my $f, '>', $tmp) || die "Couldn't open $tmp : $!\n";
	print $f $frags->[$i];
	close($f) || die "Error writing to $tmp: $!\n";
	open(my $bgz, '-|', "$$opts{bin}/bgzip -c $tmp")
	    || die "Couldn't open pipe to bgzip: $!\n";
	my $compressed = <$bgz>;
	close($bgz) || die "Error running bgzip\n";
	if ($i < $#$frags) {
	    # Strip EOF block
	    $compressed =~ s/\x1f\x8b\x08\x04\x00{5}\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00{9}$//;
	}
	print $out $compressed;
    }
    close($out) || die "Error writing to $name: $!\n";
    unlink($tmp);
}

sub test_rebgzip
{
    my ($opts, %args) = @_;

    # Write a file that should match the one we ship
    my @frags = qw(1 22 333 4444 55555);
    my $mb = "$$opts{path}/bgziptest.txt.tmp.gz";
    write_multiblock_bgzf($mb, \@frags);

    # See if it really does match
    my ($ret, $out) = _cmd("cmp $mb $$opts{path}/bgziptest.txt.gz");

    if (!$ret && $out eq '') { # If it does, use the original
	test_cmd($opts, %args, out => "bgziptest.txt.gz",
		 cmd => "$$opts{bin}/bgzip -I $$opts{path}/bgziptest.txt.gz.gzi -c -g $$opts{path}/bgziptest.txt");
    } else {
	# Otherwise index the one we just made and test that
	print "test_rebgzip: Alternate zlib/deflate library detected\n";
	cmd("$$opts{bin}/bgzip -I $mb.gzi -r $mb");
	test_cmd($opts, %args, out => "bgziptest.txt.tmp.gz",
		 cmd => "$$opts{bin}/bgzip -I $mb.gzi -c -g $$opts{path}/bgziptest.txt");
    }
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

sub test_bcf_sr_sort
{
    my ($opts, %args) = @_;
    for (my $i=0; $i<10; $i++)
    {
        my $seed = int(rand(time));
        my $test = 'test-bcf-sr';
        my $cmd  = "$$opts{path}/test-bcf-sr.pl -t $$opts{tmp} -s $seed";
        print "$test:\n";
        print "\t$cmd\n";
        my ($ret,$out) = _cmd($cmd);
        if ( $ret ) { failed($opts,$test); }
        else { passed($opts,$test); }
    }
}

sub test_command
{
    my ($opts, %args) = @_;
    my $cmd  = "$$opts{path}/$args{cmd}";
    test_cmd($opts, %args, cmd=>$cmd);
}

sub test_logging
{
  my ($opts) = @_;
  my $test = 'test-logging';
  my $cmd  = "$$opts{path}/test-logging.pl";
  print "$test:\n";
  print "\t$cmd\n";
  my ($ret,$out) = _cmd($cmd);
  if ( $ret ) { failed($opts,$test); }
  else { passed($opts,$test); }
}

sub test_realn {
    my ($opts) = @_;

    my $test_realn = "$$opts{path}/test_realn";
    # Calculate BAQ
    test_cmd($opts, cmd => "$test_realn -f $$opts{path}/realn01.fa -i $$opts{path}/realn01.sam -o -", out => "realn01_exp.sam");
    test_cmd($opts, cmd => "$test_realn -f $$opts{path}/realn02.fa -i $$opts{path}/realn02.sam -o -", out => "realn02_exp.sam");

    # Calculate and apply BAQ
    test_cmd($opts, cmd => "$test_realn -a -f $$opts{path}/realn01.fa -i $$opts{path}/realn01.sam -o -", out => "realn01_exp-a.sam");
    test_cmd($opts, cmd => "$test_realn -a -f $$opts{path}/realn02.fa -i $$opts{path}/realn02.sam -o -", out => "realn02_exp-a.sam");

    # Calculate extended BAQ
    test_cmd($opts, cmd => "$test_realn -e -f $$opts{path}/realn01.fa -i $$opts{path}/realn01.sam -o -", out => "realn01_exp-e.sam");
    test_cmd($opts, cmd => "$test_realn -e -f $$opts{path}/realn02.fa -i $$opts{path}/realn02.sam -o -", out => "realn02_exp-e.sam");

    # Recalculate BAQ
    test_cmd($opts, cmd => "$test_realn -r -f $$opts{path}/realn02.fa -i $$opts{path}/realn02-r.sam -o -", out => "realn02_exp.sam");

    # Apply from existing BQ tags
    test_cmd($opts, cmd => "$test_realn -a -f $$opts{path}/realn02.fa -i $$opts{path}/realn02_exp.sam -o -", out => "realn02_exp-a.sam");

    # Revert quality values (using data in ZQ tags)
    test_cmd($opts, cmd => "$test_realn -f $$opts{path}/realn02.fa -i $$opts{path}/realn02_exp-a.sam -o -", out => "realn02_exp.sam");
}
