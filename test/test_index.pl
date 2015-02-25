#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Temp qw(tempdir);
use Getopt::Long;
use Sys::Hostname;
use Test::More tests => 34;

my $with_irods;
my $num_irods_tests = 16;

GetOptions('with-irods' => \$with_irods);

# Test programs
my $SAMGEN     = './samgen.pl';
my $TEST_INDEX = './test_index';
my $TEST_VIEW  = './test_view';

# iRODS clients
my $IPUT = 'iput';
my $IGET = 'iget';
my $IRM  = 'irm';

my $NUM_REFS    = 3;
my $REF_SIZE    = 10000;
my $READ_LENGTH = 120;
my $SAM_FILE    = 'test.sam';

sub setup {
  my ($dir, $sam_file) = @_;

  my $cmd = "$SAMGEN " .
    "--num-refs $NUM_REFS " .
    "--ref-size $REF_SIZE " .
    "--read-length $READ_LENGTH " .
    "--output $dir/$sam_file " .
    "|| exit 1";
  if (system($cmd) != 0) {
    die "system call '$cmd' failed\n";
  }
}

sub count_matching_lines {
  my ($fh, @regexes) = @_;

  my @counts = split '', 0 x scalar @regexes;
  while (<$fh>) {
    for (my $i = 0; $i < scalar @regexes; $i++) {
      m{$regexes[$i]} and $counts[$i]++;
    }
  }

  return @counts;
}

sub test_index_range {
  my ($cmd, $ref_num, $ref_range) = @_;

  open my $in, "$cmd |" or fail("$cmd");
  my ($num_seq_header, $num_reads) =
    count_matching_lines($in, qr{^\@SQ}, qr{^read\.\d+});
  close($in);

  cmp_ok($num_seq_header, '==', 3, "Num \@SQ lines in ref$ref_num");
  cmp_ok($num_reads,      '==', 1,
         "Num reads in range ref$ref_num:$ref_range");
}

sub test_bam_index_range {
  my $tmp_dir = tempdir(getcwd . '/test_bam_index_range_XXXX', CLEANUP => 1);
  my $sam_file = 'test.sam';
  my $bam_file = 'test.bam';
  setup($tmp_dir, $sam_file);

  my $make_bam = "$TEST_VIEW $tmp_dir/$sam_file -b > $tmp_dir/$bam_file";
  system($make_bam) == 0 or warn "Failed to convert SAM to BAM";

  my $make_bam_index = "$TEST_INDEX $tmp_dir/$bam_file";
  ok(system($make_bam_index) == 0, 'BAM index');

  my $ref_range = '1-1';
  foreach my $ref_num (0 .. $NUM_REFS -1) {
    my $use_bam_index = "$TEST_VIEW $tmp_dir/$bam_file ref$ref_num:$ref_range";
    test_index_range($use_bam_index, $ref_num, $ref_range);
  }

 TODO: {
    local $TODO = 'Known issue';
    ok(system("$TEST_VIEW $tmp_dir/$bam_file ref0:2-1 >/dev/null 2>&1") != 0,
       "Should fail on BAM invalid reference range 1");
  }

  ok(system("$TEST_VIEW $tmp_dir/$bam_file ref0:2-0 >/dev/null 2>&1") != 0,
     "Should fail on BAM invalid reference range 2");
}

sub test_cram_index_range {
  my $tmp_dir = tempdir(getcwd . '/test_cram_index_range_XXXX', CLEANUP => 1);
  my $sam_file  = 'test.sam';
  my $cram_file = 'test.cram';
  my $ref_file  = 'ref.fa';
  setup($tmp_dir, $sam_file);

  my $make_cram = "$TEST_VIEW $tmp_dir/$sam_file -t $tmp_dir/$ref_file -C " .
    "> $tmp_dir/$cram_file";
  system($make_cram) == 0 or warn "Failed to convert SAM to CRAM";

  my $make_cram_index = "$TEST_INDEX $tmp_dir/$cram_file";
  ok(system($make_cram_index) == 0, 'CRAM index');

  my $ref_range = '1-1';
  foreach my $ref_num (0 .. $NUM_REFS -1) {
    my $use_cram_index =
      "$TEST_VIEW $tmp_dir/$cram_file ref$ref_num:$ref_range";
    test_index_range($use_cram_index, $ref_num, $ref_range);
  }

 TODO: {
    local $TODO = 'Known issue';
    ok(system("$TEST_VIEW $tmp_dir/$cram_file ref0:2-1 >/dev/null 2>&1") != 0,
       "Should fail on CRAM invalid reference range 1");
  }

  ok(system("$TEST_VIEW $tmp_dir/$cram_file ref0:2-0 >/dev/null 2>&1") != 0,
     "Should fail on CRAM invalid reference range 2");
}

sub test_irods_local_bam_index_range {
  my $tmp_dir =
    tempdir(getcwd . '/test_irods_local_bam_index_range_XXXX', CLEANUP => 1);
  my $sam_file = 'test.sam';
  my $bam_file = sprintf "test.%s.%d.bam", hostname, $$;
  setup($tmp_dir, $sam_file);

  my $make_bam = "$TEST_VIEW $tmp_dir/$sam_file -b -O irods:$bam_file";
  system($make_bam) == 0 or warn "Failed to convert SAM to BAM in iRODS";

  my $make_bam_index = "$TEST_INDEX irods:$bam_file";
  ok(system($make_bam_index) == 0, 'BAM index iRODS');
  ok(-f "irods:$bam_file.bai", 'Local iRODS BAM index');

  my $ref_range = '1-1';
  foreach my $ref_num (0 .. $NUM_REFS -1) {
    my $use_bam_index =
      "$TEST_VIEW irods:$bam_file ref$ref_num:$ref_range";
    test_index_range($use_bam_index, $ref_num, $ref_range);
  }

  unlink "irods:$bam_file.bai" or
    warn "Failed to clean up local irods:$bam_file.bai";

  system("$IRM $bam_file") == 0 or
    warn "Failed to clean up remote $bam_file";
}

sub test_irods_remote_bam_index_range {
  my $tmp_dir =
    tempdir(getcwd . '/test_irods_remote_bam_index_range_XXXX', CLEANUP => 1);
  my $sam_file = 'test.sam';
  my $bam_file = sprintf "test.%s.%d.bam", hostname, $$;
  setup($tmp_dir, $sam_file);

  my $make_bam = "$TEST_VIEW $tmp_dir/$sam_file -b -O $tmp_dir/$bam_file";
  system($make_bam) == 0 or warn "Failed to convert SAM to BAM";
  system("$IPUT $tmp_dir/$bam_file") == 0 or
    warn 'Failed to uploaded a BAM file to iRODS';

  my $make_bam_index = "$TEST_INDEX $tmp_dir/$bam_file";
  ok(system($make_bam_index) == 0, 'BAM index');

  system("$IPUT $tmp_dir/$bam_file.bai") == 0 or
    warn 'Failed to uploaded a BAM index to iRODS';
  ok(! -f "irods:$bam_file.bai", 'No local iRODS BAM index');

  my $ref_range = '1-1';
  foreach my $ref_num (0 .. $NUM_REFS -1) {
    my $use_bam_index =
      "$TEST_VIEW irods:$bam_file ref$ref_num:$ref_range";
    test_index_range($use_bam_index, $ref_num, $ref_range);
  }

  unlink "irods:$bam_file.bai" or
    warn "Failed to clean up local irods:$bam_file.bai";

  system("$IRM $bam_file") == 0 or
    warn "Failed to clean up remote $bam_file";
  system("$IRM $bam_file.bai") == 0 or
    warn "Failed to clean up remote $bam_file.bai";
}

# Run the tests
test_bam_index_range;
test_cram_index_range;

SKIP: {
  skip 'testing without iRODS', $num_irods_tests unless $with_irods;

  test_irods_local_bam_index_range;
  test_irods_remote_bam_index_range;
}

__END__

=head1 NAME

 test_index.pl

=head1 SYNOPSIS

  test_index.pl [--with-irods]

Options:

  --with-irods  Enable iRODS tests (normally skipped).

=head1 DESCRIPTION

    As simple test harness for BAM/CRAM indexing. Currently it does
    some cursory tests for the creation and usability of indices. It
    does not test that index queries return the correct results, for
    example.

=head1 AUTHOR

Keith James <kdj@sanger.ac.uk>

=head1 COPYRIGHT AND DISCLAIMER

Copyright (C) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

=cut
