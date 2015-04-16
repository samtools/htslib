#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $num_refs;
my $output;
my $read_len;
my $ref_name;
my $ref_size;
my $ref_start;
my $ref_stride;
my $ref_end;

GetOptions('help'          => sub { pod2usage(-verbose => 2, -exitval => 0) },
           'num-refs=i'    => \$num_refs,
           'output=s'      => \$output,
           'read-length=i' => \$read_len,
           'ref-end=i'     => \$ref_end,
           'ref-name=s'    => \$ref_name,
           'ref-size=i'    => \$ref_size,
           'ref-start=i'   => \$ref_start,
           'ref-stride=i'  => \$ref_stride);

my $ref_name_regex = '^[A-Za-z0-9-_.]+$';
$ref_name ||= 'ref';
$ref_name =~ m{$ref_name_regex} or
  die "--ref-name must match $ref_name_regex, was '$ref_name'";

# Use this rather than ||= so that if a user explicity sets a value to
# 0, they get an error message rather than their value being
# overridden silently with a default.
!defined $num_refs   and $num_refs   = 1;
!defined $ref_size   and $ref_size   = 1;
!defined $ref_start  and $ref_start  = 1;
!defined $ref_stride and $ref_stride = 1;
!defined $ref_end    and $ref_end    = $ref_size;

$read_len ||= 0;

$num_refs > 0 or pod2usage
  (-msg     => "--num-refs must be > 0, was $num_refs\n", -exitval => 2);

$read_len >= 0 or pod2usage
  (-msg     => "--read-length must be >= 0, was $read_len\n", -exitval => 2);

$ref_size >= 1 or pod2usage
  (-msg     => "--ref-size must be >= 1, was $ref_size\n", -exitval => 2);

$ref_start >= 1 or pod2usage
  (-msg     => "--ref-start must be >= 1, was $ref_start\n", -exitval => 2);

$ref_end >= 1 or pod2usage
  (-msg     => "--ref-end must be >= 1, was $ref_end\n", -exitval => 2);

$ref_stride >= 1 or pod2usage
  (-msg     => "--ref-stride must be >= 1, was $ref_stride\n", -exitval => 2);

$ref_end <= $ref_size or pod2usage
  (-msg     => "--ref-end must be <= the reference size [$ref_size], " .
   "was $ref_end\n", -exitval => 2);

$ref_start <= $ref_end or pod2usage
  (-msg     => "--ref-start must be <= the reference range end [$ref_end], " .
   "was $ref_start\n", -exitval => 2);

my $sam;
my $ref_path = '.';
if ($output) {
  open $sam, '>', $output or die "Failed to open $output: $!\n";

  my ($name, $path, $suffix) = fileparse($output, '.sam');
  if ($path) {
    $ref_path = $path;
  }
}
else {
  $sam = \*STDOUT;
}

# Write reference Fasta files and SAM header. Write Fasta to the same
# directory as the SAM file if a file path is given.


my $fa_name = "$ref_path/$ref_name.fa";
open my $fa, '>', $fa_name or die "Failed to open $fa_name: $!\n";

my %refs;
my @ref_names = map { $ref_name . $_ } 0 .. ($num_refs -1);
foreach my $name (@ref_names) {
  my $seq = rand_seq($ref_size);
  $refs{$name} = $seq;

  write_fasta($fa, $name, $seq);
  write_sam_seq_header($sam, $name, $seq);
}

close $fa or warn "Failed to close $fa_name: $!\n";

# Write SAM records
my $read_count = 0;
foreach my $name (@ref_names) {
  for (my $i = $ref_start,
       my $j = $ref_start + $read_len -1;
       $i <= $ref_end && $j <= $ref_end;
       $i += $ref_stride,
       $j += $ref_stride) {

    my $qname = "read.$read_count";
    my $flag  = 16;
    my $rname = $name;
    my $pos   = $i;
    my $mapq  = 40;
    my $seq   = substr $refs{$name}, $pos - 1, $read_len;
    my $qual  = "I" x $read_len;
    my $cigar = $read_len . 'M';

    write_sam_record($sam, $qname, $flag, $rname, $pos, $mapq, $seq,
                     $qual, $cigar);
    $read_count++;
  }
}

if ($output) {
  close $sam or warn "Failed to close $output: $!\n";
}

sub write_sam_seq_header {
  my ($fh, $name, $seq) = @_;

  my $len = length $seq;
  print $fh q{@SQ} . "\tSN:$name\tLN:$len\n";
}

sub write_sam_record {
  my ($fh, $qname, $flag, $rname, $pos, $mapq,
      $seq, $qual, $cigar, $rnext, $pnext, $tlen) = @_;

  $seq   ||= '*';
  $qual  ||= '*';
  $cigar ||= '*';
  $rnext ||= '*';
  $pnext ||= '0';
  $tlen  ||= '0';

  print $fh join("\t", $qname, $flag, $rname, $pos, $mapq, $cigar,
                 $rnext, $pnext, $tlen, $seq, $qual), "\n";
}

sub write_fasta {
  my ($fh, $name, $seq) = @_;

  $seq =~ s/(.{1,60})/$1\n/g;

  print $fh ">$name\n$seq";
}

sub rand_seq {
  my ($len) = @_;

  my @bases = ('A', 'C', 'G', 'T');
  my $seq = q{};
  for (my $i = 0; $i < $len; $i++) {
    $seq .= $bases[int(rand(4))];
  }

  return $seq;
}

__END__

=head1 NAME

 samgen.pl - Generate simple, regular SAM files for tests.

=head1 SYNOPSIS

  samgen.pl [--help] [--num-refs <int>] [--output <filename>]
    [--read-length <int>]  [--ref-end <int>] [--ref-name <name>]
    [--ref-size <int>] [--ref-start <int>] [--ref-stride <int>]

Options:

  --help        Print this help.
  --num-refs    The number of reference sequences to generate.
                Optional, defaults to 1.
  --output      The SAM file to write. Optional, defaults to STDOUT.
  --read-length The read length for all reads. Optional, defaults to 0.
  --ref-end     The 1-based position on the reference beyond which no
                part of any read will map. Optional, defaults 1.
  --ref-name    The base name of the reference sequences. This must match
                the regex ^[A-Za-z0-9-_.]+$. Optional, defaults to 'ref';
  --ref-size    The size (length) of the reference sequence. Optional,
                defaults to 1.
  --ref-start   The 1-based position on the reference before which no
                part of any read will map. Optional, defaults 1.
  --ref-stride  The number of bases between the position of each
                successive read. Optional, defaults to 1.

=head1 DESCRIPTION

      This program generates simple, regularly structured SAM files
  and corresponding random reference sequences in Fasta format. A read
  is generated at a start position and then every n (stride) bases up
  to (but not including) the position where the read would extend
  beyond an end position.

  The reads it generates are unpaired and mapped with a score of 40.

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
