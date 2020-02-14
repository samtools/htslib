#!/usr/bin/env perl
#
#    Copyright (C) 2012-2019 Genome Research Ltd.
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

test_bgzip($opts, 0);
test_bgzip($opts, 4);

ce_fa_to_md5_cache($opts);
test_index($opts, 0);
test_index($opts, 4);

test_multi_ref($opts,0);
test_multi_ref($opts,4);

test_view($opts,0);
test_view($opts,4);

test_MD($opts);

test_vcf_api($opts,out=>'test-vcf-api.out');
test_bcf2vcf($opts);
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

# Run cmd, producing file out, and compare contents against exp
sub test_compare
{
    my ($opts,$cmd,$exp_fn,$out_fn, %args) = @_;
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n\t$cmd\n";

    my ($ret,$stdout) = _cmd($cmd);
    if ( $ret ) { failed($opts,$test); return; }

    local $/;
    my ($exp,$out) = ("","");
    if ( exists($args{"gz"}) ) {
        if ( open(my $fh,'-|',"$$opts{bin}/bgzip -d < $exp_fn") ) {
            $exp = <$fh>;
            close($fh);
        } else {
            failed($opts,$test,"bgzip -d < $exp_fn $!"); return;
        }
    } else {
        if ( open(my $fh,'<',$exp_fn) ) {
            $exp = <$fh>;
            close($fh);
        } else {
            failed($opts,$test,"$exp_fn $!"); return;
        }
    }

    if ( exists($args{"gz"}) ) {
        if ( open(my $fh,'-|',"$$opts{bin}/bgzip -d < $out_fn") ) {
            $out = <$fh>;
            close($fh);
        } else {
            failed($opts,$test,"bgzip -d < $out_fn $!"); return;
        }
    } else {
        if ( open(my $fh,'<',$out_fn) ) {
            $out = <$fh>;
            close($fh);
        } else {
            failed($opts,$test,"$out_fn $!"); return;
        }
    }

    if (exists($args{fix_newlines})) {
        $exp =~ s/\015\012/\n/g;
        $out =~ s/\015\012/\n/g;
    }

    if ( $exp ne $out )
    {
        failed($opts,$test,"The outputs differ:\n\t\t$exp_fn\n\t\t$out_fn");
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

sub ce_fa_to_md5_cache {
    my ($opts) = @_;

    # These should really be worked out from the file contents, but
    # pre-calculating them avoids a dependency on Digest::MD5
    my %csums = (CHROMOSOME_I     => '8ede36131e0dbf3417807e48f77f3ebd',
                 CHROMOSOME_II    => '8e7993f7a93158587ee897d7287948ec',
                 CHROMOSOME_III   => '3adcb065e1cf74fafdbba1e8c352b323',
                 CHROMOSOME_IV    => '251af66a69ee589c9f3757340ec2de6f',
                 CHROMOSOME_V     => 'cf200a65fb754836dcc56b24b3170ee8',
                 CHROMOSOME_X     => '6f9368fd2192c89c613718399d2d31fc',
                 CHROMOSOME_MtDNA => 'cd05857ece6411f40257a565ccfe15bb');

    my $m5_dir = "$$opts{tmp}/md5";
    if (!-d $m5_dir) {
        mkdir($m5_dir) || die "Couldn't make directory $m5_dir\n";
    }
    my $out;
    open(my $fa, '<', "$$opts{path}/ce.fa")
        || die "Couldn't open $$opts{path}/ce.fa : $!\n";
    my $name = '';
    while (<$fa>) {
        chomp;
        if (/^>(\S+)/) {
            if ($out) {
                close($out) || die "Error closing $m5_dir/$csums{$name} : $!\n";
            }
            $name = $1;
            if (!exists($csums{$name})) {
                die "Unexpected fasta entry : $name\n";
            }
            open($out, '>', "$m5_dir/$csums{$name}")
        } else {
            if (!$out) {
                die "$$opts{path}/ce.fa : Got data before fasta header\n";
            }
            $_ = uc($_);
            s/\s+//g;
            print $out $_;
        }
    }
    if ($out) {
        close($out) || die "Error closing $m5_dir/$csums{$name} : $!\n";
    }
    close($fa) || die "Error reading $$opts{path}/ce.fa : $!\n";
    $$opts{m5_dir} = $m5_dir;
}


# The tests --------------------------

sub test_bgzip {
    my ($opts, $threads) = @_;

    my $at = $threads ? "-@ $threads" : '';
    my $data = "$$opts{path}/ce.fa";
    my $compressed = "$$opts{tmp}/ce.fa.$threads.gz";
    my $compressed_copy = "$$opts{tmp}/ce.fa.$threads.copy.gz";
    my $uncompressed = "$$opts{tmp}/ce.fa.$threads.uncomp";
    my $offset = 1055584; # Start of MT in ce.fa
    my $uncompressed_part = "$$opts{tmp}/ce.fa.$threads.part";
    my $uncompressed_part2 = "$$opts{tmp}/ce.fa.$threads.part2";
    my $expected_part = "$$opts{tmp}/ce.fa.$threads.tail";
    my $index = "${compressed}.gzi";
    my $test = sprintf('%s %2s threads', 'bgzip round-trip',
                       $threads ? $threads : 'no');

    # Round-trip test
    print "$test: ";
    my $c = "$$opts{bin}/bgzip $at -i -I '$index' < '$data' > '$compressed'";
    my ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "$$opts{bin}/bgzip $at -d < '$compressed' > '$uncompressed'";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "cmp '$data' '$uncompressed'";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, $out ? $out : "'$data' '$uncompressed' differ");
        return;
    }
    passed($opts,$test);

    # Extract from an offset
    $test = sprintf('%s %2s threads', 'bgzip -b',
                    $threads ? $threads : 'no');
    print "$test: ";
    $c = sprintf("tail -c +%d '%s' > '%s'", $offset + 1, $data, $expected_part);
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "$$opts{bin}/bgzip $at -b $offset -d '$compressed' > $uncompressed_part";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "cmp '$expected_part' '$uncompressed_part'";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test,
               $out ? $out : "'$expected_part' '$uncompressed_part' differ");
        return;
    }
    passed($opts,$test);

    # Extract from an offset with named index
    $test = sprintf('%s %2s threads', 'bgzip -b -I',
                    $threads ? $threads : 'no');
    print "$test: ";
    $c = "cp '$compressed' '$compressed_copy'";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "$$opts{bin}/bgzip $at -b $offset -d -I '$index' '$compressed_copy' > $uncompressed_part2";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test, "non-zero exit from $c");
        return;
    }
    $c = "cmp '$expected_part' '$uncompressed_part2'";
    ($ret, $out) = _cmd($c);
    if ($ret) {
        failed($opts, $test,
               $out ? $out : "'$expected_part' '$uncompressed_part2' differ");
        return;
    }
    passed($opts,$test);
}

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

sub fake_multi_ref_data
{
    open(SAM, ">multi_ref.tmp.sam") || die;
    for (my $r=0;$r<1000;$r++) {
        print SAM "\@SQ\tSN:c$r\tLN:10000\n";
    }

    # Single ref
    my $rnum=0;
    for (my $p=1;$p<1000;$p++) {
        print SAM "X\t0\tc$rnum\t$p\t40\t10M\t*\t0\t0\tCCTAGCCCTA\tB?8B?BACCD\n";
    }

    # Multi ref; 1 seq per ref
    for (my $r=1;$r<300;$r++) {
        print SAM "X\t0\tc$rnum\t1\t40\t10M\t*\t0\t0\tCCTAGCCCTA\tB?8B?BACCD\n";
        $rnum++;
    }

    # Single ref again
    for (my $p=1;$p<1000;$p++) {
        print SAM "X\t0\tc$rnum\t$p\t40\t10M\t*\t0\t0\tCCTAGCCCTA\tB?8B?BACCD\n";
    }

    # Multi ref; 1 seq per ref
    for (my $r=1;$r<300;$r++) {
        print SAM "X\t0\tc$rnum\t1\t40\t10M\t*\t0\t0\tCCTAGCCCTA\tB?8B?BACCD\n";
        $rnum++;
    }
    close(SAM);
}

sub test_multi_ref
{
    my ($opts, $nthreads) = @_;
    my $tv_args = $nthreads ? "-\@$nthreads" : "";

    fake_multi_ref_data;
    print "test_view testing multi-ref CRAM modes:\n";
    $test_view_failures = 0;

    for (my $mf = -1; $mf <= 1; $mf++) {
        testv $opts, "./test_view $tv_args -o seqs_per_slice=100 -o no_ref=1 -o multi_seq_per_slice=$mf -S -C multi_ref.tmp.sam > multi_ref.tmp.cram";
        testv $opts, "./test_view $tv_args multi_ref.tmp.cram > multi_ref.tmp.sam_";
        testv $opts, "./compare_sam.pl multi_ref.tmp.sam multi_ref.tmp.sam_";
    }

    if ($test_view_failures == 0) {
        passed($opts, "multi-ref conversions");
    } else {
        failed($opts, "multi-ref conversions", "$test_view_failures subtests failed");
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
    print "test_view testing region queries:\n";
    $test_view_failures = 0;

    my $regions = "CHROMOSOME_II:2980-2980 CHROMOSOME_IV:1500-1500 CHROMOSOME_II:2980-2980 CHROMOSOME_I:1000-1100";
    testv $opts, "./test_view $tv_args -i reference=ce.fa range.cram $regions > range.tmp";
    testv $opts, "./compare_sam.pl range.tmp range.out";

    testv $opts, "./test_view $tv_args range.bam $regions > range.tmp";
    testv $opts, "./compare_sam.pl range.tmp range.out";

    if ($test_view_failures == 0) {
        passed($opts, "range.cram tests");
    } else {
        failed($opts, "range.cram tests", "$test_view_failures subtests failed");
    }

    # Test BAM files with references in targets list but no corresponding @SQ
    # lines in the text header.
    print "test_view testing BAM files with absent \@SQ lines:\n";
    $test_view_failures = 0;
    testv $opts, "./test_view $tv_args -p no_hdr_sq_1.tmp.sam no_hdr_sq_1.bam";
    testv $opts, "./compare_sam.pl no_hdr_sq_1.tmp.sam no_hdr_sq_1.expected.sam";

    # Try a range query to ensure id <-> name mapping works
    # Input only has reads from CHROMOSOME_I, so same "expected" file is used
    testv $opts, "./test_view $tv_args -p no_hdr_sq_1.chr1.tmp.sam no_hdr_sq_1.bam CHROMOSOME_I";
    testv $opts, "./compare_sam.pl no_hdr_sq_1.chr1.tmp.sam no_hdr_sq_1.expected.sam";
    if ($test_view_failures == 0) {
        passed($opts, "no_hdr_sq tests");
    } else {
        failed($opts, "no_hdr_sq tests", "$test_view_failures subtests failed");
    }

    # File with large (> 2Gbases) positions
    # Only works for SAM at the moment, but we can still round-trip it.
    print "test_view testing large (> 2Gbases) positions:\n";
    $test_view_failures = 0;
    testv $opts, "./test_view $tv_args -z -p longrefs/longref.tmp.sam.gz -x longrefs/longref.tmp.sam.gz.csi.otf -m 14 longrefs/longref.sam";
    testv $opts, "./test_view $tv_args -p longrefs/longref.tmp.sam_ longrefs/longref.tmp.sam.gz";
    testv $opts, "./compare_sam.pl longrefs/longref.sam longrefs/longref.tmp.sam_";

    # CRAM disabled for now as the positions cannot be 32-bit.  (These tests are useful for
    # checking SQ headers only.)
    # testv $opts, "./test_view $tv_args -C -o no_ref -p longrefs/longref.tmp.cram longrefs/longref.sam";
    # testv $opts, "./test_view $tv_args -p longrefs/longref.tmp.sam_ longrefs/longref.tmp.cram";
    # testv $opts, "./compare_sam.pl longrefs/longref.sam longrefs/longref.tmp.sam_";

    # Build index and compare with on-the-fly one made earlier.
    test_compare $opts, "$$opts{path}/test_index -c longrefs/longref.tmp.sam.gz", "longrefs/longref.tmp.sam.gz.csi.otf", "longrefs/longref.tmp.sam.gz.csi", gz=>1;

    # Large position iterator tests
    testv $opts, "./test_view $tv_args -p longrefs/longref_itr.tmp.sam longrefs/longref.tmp.sam.gz CHROMOSOME_I:10000000000-10000000003";
    testv $opts, "./compare_sam.pl longrefs/longref_itr.expected.sam longrefs/longref_itr.tmp.sam";
    testv $opts, "./test_view $tv_args -M -p longrefs/longref_multi.tmp.sam longrefs/longref.tmp.sam.gz CHROMOSOME_I:10000000000-10000000003 CHROMOSOME_I:10000000100-10000000110";
    testv $opts, "./compare_sam.pl longrefs/longref_multi.expected.sam longrefs/longref_multi.tmp.sam";

    # 64-bit positions are currently not compiled in by default for VCF
    #   # VCF round trip
    #   unlink("longrefs/index.tmp.vcf.gz.csi"); # To stop vcf_hdr_read from reading a stale index
    #   testv $opts, "./test_view $tv_args -z -p longrefs/index.tmp.vcf.gz -x longrefs/index.tmp.vcf.gz.csi.otf -m 14 longrefs/index.vcf";
    #   testv $opts, "./test_view $tv_args -p longrefs/index.tmp.vcf_ longrefs/index.tmp.vcf.gz";
    #   testv $opts, "cmp longrefs/index.vcf longrefs/index.tmp.vcf_";
    #
    #   # Build index and compare with on-the-fly one made earlier.
    #   test_compare $opts, "$$opts{path}/test_index -c longrefs/index.tmp.vcf.gz", "longrefs/index.tmp.vcf.gz.csi.otf", "longrefs/index.tmp.vcf.gz.csi", gz=>1;
    #
    #   # test_view can't do indexed look-ups on vcf, but we can use tabix
    #   test_compare $opts, "$$opts{bin}/tabix longrefs/index.tmp.vcf.gz 1:10010000100-10010000105 > longrefs/index.tmp.tabix1.vcf", "longrefs/index.expected1.vcf", "longrefs/index.tmp.tabix1.vcf", fix_newlines => 1;
    #   test_compare $opts, "$$opts{bin}/tabix longrefs/index.tmp.vcf.gz 1:10010000120-10010000130 > longrefs/index.tmp.tabix2.vcf", "longrefs/index.expected2.vcf", "longrefs/index.tmp.tabix2.vcf", fix_newlines => 1;

    if ($test_view_failures == 0) {
        passed($opts, "large position tests");
    } else {
        failed($opts, "large position tests", "$test_view_failures subtests failed");
    }
}

# Tests CRAM's ability to correctly preserve MD and NM, irrespective of whether
# they are correct.
sub test_MD
{
    my ($opts) = @_;

    foreach my $sam (glob("*#MD*.sam")) {
        my ($base, $ref) = ($sam =~ /((.*)#.*)\.sam/);
        $ref .= ".fa";

        my $bam  = "$base.tmp.bam";
        my $cram = "$base.tmp.cram";

        print "\ntest_MD testing $sam, ref $ref:\n";
        $test_view_failures = 0;
        $cram = "$base.tmp.cram";

        # Forcibly store MD and NM and don't auto-generate.
        # ALL NM/MD should match and be present only when originally present
        testv $opts, "./test_view -o store_nm=1 -o store_md=1 -t $ref -C $sam > $cram";
        testv $opts, "./test_view -i decode_md=0 -D $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl $sam $cram.sam_";

        # Skip auto-MD generation; check MD iff in output file.
        # (NB this does not check that all erroneous values are stored.)
        testv $opts, "./test_view -t $ref -C $sam > $cram";
        testv $opts, "./test_view -i decode_md=0 -D $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl -partialmd=2 $sam $cram.sam_";

        # Also check we haven't added NM or MD needlessly for xx#MD.sam.
        # This file has no errors so without auto-generation there must be
        # no NM or MD records.
        if ($sam eq "xx#MD.sam") {
            print "  Checking for MD/NM in $sam\n";
            open(my $fh, "<$cram.sam_") || die;
            while (<$fh>) {
                if (/(MD|NM):/) {
                    print STDERR "Failed\nLine contains MD/NM:\n$_";
                    $test_view_failures++;
                    last;
                }
            }
            close($fh);
        }

        # Force auto-MD generation; check MD iff in input file.
        # This will ensure any erroneous values have been round-tripped.
        testv $opts, "./test_view -t $ref -C $sam > $cram";
        testv $opts, "./test_view -i decode_md=1 -D $cram > $cram.sam_";
        testv $opts, "./compare_sam.pl -partialmd=1 $sam $cram.sam_";

        if ($test_view_failures == 0) {
            passed($opts, "$sam MD tests");
        } else {
            failed($opts, "$sam MD tests", "$test_view_failures subtests failed");
        }
    }
}

sub test_index
{
    my ($opts, $nthreads) = @_;
    $nthreads = $nthreads ? "-\@$nthreads" : "";

    # BAM
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -b -m 14 -x $$opts{tmp}/index.bam.csi $$opts{path}/index.sam > $$opts{tmp}/index.bam", "$$opts{tmp}/index.bam.csi", "$$opts{path}/index.bam.csi", gz=>1);
    unlink("$$opts{tmp}/index.bam.csi");
    test_compare($opts,"$$opts{path}/test_index -c $$opts{tmp}/index.bam", "$$opts{tmp}/index.bam.csi", "$$opts{path}/index.bam.csi", gz=>1);
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -b -m 0 -x $$opts{tmp}/index.bam.bai $$opts{path}/index.sam > $$opts{tmp}/index.bam", "$$opts{tmp}/index.bam.bai", "$$opts{path}/index.bam.bai");
    unlink("$$opts{tmp}/index.bam.bai");
    test_compare($opts,"$$opts{path}/test_index -b $$opts{tmp}/index.bam", "$$opts{tmp}/index.bam.bai", "$$opts{path}/index.bam.bai");

    # SAM
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -z -m 14 -x $$opts{tmp}/index.sam.gz.csi $$opts{path}/index.sam > $$opts{tmp}/index.sam.gz", "$$opts{tmp}/index.sam.gz.csi", "$$opts{path}/index.sam.gz.csi", gz=>1);
    unlink("$$opts{tmp}/index.bam.bai");
    test_compare($opts,"$$opts{path}/test_index -c $$opts{tmp}/index.sam.gz", "$$opts{tmp}/index.sam.gz.csi", "$$opts{path}/index.sam.gz.csi", gz=>1);
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -z -m 0 -x $$opts{tmp}/index.sam.gz.bai $$opts{path}/index.sam > $$opts{tmp}/index.sam.gz", "$$opts{tmp}/index.sam.gz.bai", "$$opts{path}/index.sam.gz.bai");
    unlink("$$opts{tmp}/index.sam.gz.bai");
    test_compare($opts,"$$opts{path}/test_index -b $$opts{tmp}/index.sam.gz", "$$opts{tmp}/index.sam.gz.bai", "$$opts{path}/index.sam.gz.bai");

    # CRAM
    local $ENV{REF_PATH} = $$opts{m5_dir};
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -C -x $$opts{tmp}/index.cram.crai $$opts{path}/index.sam > $$opts{tmp}/index.cram", "$$opts{tmp}/index.cram.crai", "$$opts{path}/index.cram.crai", gz=>1);
    unlink("$$opts{tmp}/index.cram.crai");
    test_compare($opts,"$$opts{path}/test_index $$opts{tmp}/index.cram", "$$opts{tmp}/index.cram.crai", "$$opts{path}/index.cram.crai", gz=>1);

    # BCF
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -b -m 14 -x $$opts{tmp}/index.bcf.csi $$opts{path}/index.vcf > $$opts{tmp}/index.bcf", "$$opts{tmp}/index.bcf.csi", "$$opts{path}/index.bcf.csi", gz=>1);
    unlink("$$opts{tmp}/index.bcf.csi");
    test_compare($opts,"$$opts{path}/test_index -c $$opts{tmp}/index.bcf", "$$opts{tmp}/index.bcf.csi", "$$opts{path}/index.bcf.csi", gz=>1);

    # VCF
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -z -m 14 -x $$opts{tmp}/index.vcf.gz.csi $$opts{path}/index.vcf > $$opts{tmp}/index.vcf.gz", "$$opts{tmp}/index.vcf.gz.csi", "$$opts{path}/index.vcf.gz.csi", gz=>1);
    unlink("$$opts{tmp}/index.vcf.gz.csi");
    test_compare($opts,"$$opts{path}/test_index -c $$opts{tmp}/index.vcf.gz", "$$opts{tmp}/index.vcf.gz.csi", "$$opts{path}/index.vcf.gz.csi", gz=>1);
    test_compare($opts,"$$opts{path}/test_view $nthreads -l 0 -z -m 0 -x $$opts{tmp}/index.vcf.gz.tbi $$opts{path}/index.vcf > $$opts{tmp}/index.vcf.gz", "$$opts{tmp}/index.vcf.gz.tbi", "$$opts{path}/index.vcf.gz.tbi", gz=>1);
    unlink("$$opts{tmp}/index.vcf.gz.tbi");
    test_compare($opts,"$$opts{path}/test_index -t $$opts{tmp}/index.vcf.gz", "$$opts{tmp}/index.vcf.gz.tbi", "$$opts{path}/index.vcf.gz.tbi", gz=>1);

    # Tabix and custom index names
    _cmd("$$opts{bin}/tabix -fp vcf $$opts{tmp}/index.vcf.gz");
    my $wtmp = $$opts{tmp};
    if ($^O =~ /^msys/) {
        $wtmp =~ s/\//\\\\/g;
    }
    test_cmd($opts,out=>'tabix.out',cmd=>"$$opts{bin}/tabix $wtmp/index.vcf.gz##idx##$wtmp/index.vcf.gz.tbi 1:10000060-10000060");
}

sub test_bcf2vcf
{
    my ($opts) = @_;
    test_cmd($opts,
             out => "tabix/vcf_file.vcf",
             cmd => "$$opts{path}/test_view $$opts{path}/tabix/vcf_file.bcf");
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
  if ( $ret ) {
      print $out;
      failed($opts,$test);
  }
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
