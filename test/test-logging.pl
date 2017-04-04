#!/usr/bin/env perl
#
#    Copyright (C) 2017 Genome Research Ltd.
#
#    Author: Anders Kaplan
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

my $log_statement_count = 0;
my $file_count = 0;
my $failure_count = 0;

sub check_log_statement
{
  my ($argsstr, $filename) = @_;
  $log_statement_count++;

  # collapse function calls to get rid of extra commas
  my $collapsed_argsstr = $argsstr;
  $collapsed_argsstr =~ s/\([^\)]*\)/()/g;

  my @args = split(/\s*,\s*/, $collapsed_argsstr);
  my $format = @args[0];

  unless ($format =~ /\\n\"$/)
  {
    print "$filename:\n";
    print "Missing newline in format string: $format.\n";
    $failure_count++;
  }

  my @placeholders = $format =~ /%[^%]/g;
  unless (scalar(@placeholders) == scalar(@args) - 1)
  {
    print "$filename:\n";
    print "Number of placeholders in the format string doesn't match the number of arguments: $argsstr\n";
    $failure_count++;
  }
}

sub check_file
{
  my ($filename) = @_;
  $file_count++;

  open(my $fh, '<', $filename) or die "Could not open $filename.";
  my @lines = <$fh>;

  # strip line comments
  foreach my $line (@lines)
  {
    $line =~ s/\/\/.*$//;
  }

  my $text = join(' ', @lines);

  # strip /* */ comments
  $text =~ s/\/\*([^\*]|\*[^\/])*\*\// /g;

  while ($text =~ /log_\w+\(\s*([^;]+)\s*\)\s*;/g)
  {
    check_log_statement($1, $filename);
  }
}

sub check_dir
{
  my ($path) = @_;
  foreach my $filename (glob("$path/*.c"))
  {
    check_file($filename);
  }
}

check_dir("..");
check_dir("../cram");

print "$file_count files scanned\n";
print "$log_statement_count log statements checked\n";
print "$failure_count errors found\n";
exit($failure_count > 0);
