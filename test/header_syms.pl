#!/usr/bin/env perl
#
#    Copyright (C) 2019 Genome Research Ltd.
#
#    Author: Rob Davies <rmd@sanger.ac.uk>
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

# Script to extract symbol names from HTSlib header files.  Used to
# check the shared library for missing exports.

# Instead of implementing a full C parser, this attempts to do the minimum
# amount it can get away with by scrubbing out most of the header text and
# then looking through the rest for function declarations.

# Roughly equivalent Exuberant-ctags command is:
# ctags -f - -n -I HTS_RESULT_USED -I HTS_DEPRECATED+ -I HTS_FORMAT+ \
#       -I KS_ATTR_PRINTF+ -I knet_win32_destroy+ -I knet_win32_init+
# Unfortunately this is not the default ctags on all platforms, hence this
# script.

use strict;
use warnings;
use Getopt::Long;

# Use this option to show the processed version of the header text
# instead of the function list.
my $show_processed = 0;

GetOptions('show-processed' => \$show_processed);

# List of functions to strip from the output
my %ignore = map { $_ => 1 } qw(knet_win32_init knet_win32_destroy);

foreach my $file (@ARGV) {
    extract_symbols($file, $show_processed, \%ignore);
}

sub extract_symbols {
    my ($file, $show_processed, $ignore) = @_;

    local $/ = undef;

    open(my $f, '<', $file) || die "Couldn't open $file : $!\n";
    my $text = <$f>;
    close($f) || die "Error reading $file : $!\n";

    # Get rid of comments
    $text =~ s#/\*.*?\*/##sg;
    $text =~ s#//.*$##mg;

    # Remove extern "C" brackets
    $text =~ s/#ifdef\s+__cplusplus.*?#endif//sg;

    # Remove #if 0 sections
    $text =~ s/^\s*#\s*if\s+0\s+.*?#\s*endif\s//msg;

    # Remove #defines
    $text =~ s/\n\s*?#\s*?define\s+(?:[^\n]+\\\n)*[^\n]+//sg;

    # Remove content inside curly braces
    $text =~ s/(\{(?:(?>[^{}]+)|(?1))*\})/{}/sg;

    # Get rid of typedefs
    $text =~ s/typedef\s+[^;]+;//sg;

    # Get rid of some macros
    $text =~ s/HTS_RESULT_USED//g;
    $text =~ s/HTSLIB_EXPORT//g;

    $text =~ s/HTS_DEPRECATED\s*?\(\"[^"]+\"?\)//g;
    $text =~ s/HTS_FORMAT\s*?\(.*?\)//g;
    $text =~ s/KS_ATTR_PRINTF\s*?\(.*?\)//g;

    # Get rid of static inline functions
    $text =~ s/static\s+inline\s+(?:\S+\s+)+?(\S+)\s*(\((?:(?>[^()]+)|(?-1))*\))\s*{}//g;

    if ($show_processed) {
        print $text;
        return;
    }

    # Find functions and print them
    while ($text =~ m/^\s+(?:\S+\s+)+?(?:\*+\s*)?(\S+)\s*(\((?:(?>[^()]+)|(?-1))*\))\s*;/msg) {
        next if (exists($ignore->{$1}));
        print "$1\n";
    }
}
