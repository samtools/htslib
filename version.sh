#!/bin/sh
# version.sh -- Script to build the htslib version string
#
#     Author : James Bonfield <jkb@sanger.ac.uk>
#
#     Copyright (C) 2017-2018, 2021 Genome Research Ltd.
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

# Master version, for use in tarballs or non-git source copies
VERSION=1.17

# If we have a git clone, then check against the current tag
srcdir=${0%/version.sh}
if [ -e $srcdir/.git ]
then
    # If we ever get to 10.x this will need to be more liberal
    v=`cd $srcdir && git describe --always --match '[0-9].[0-9]*' --dirty`
    case $v in
        [0-9]*.[0-9]*) VERSION="$v" ;;
        [0-9a-f][0-9a-f]*) VERSION="$VERSION-1-g$v" ;;
    esac
fi

# Numeric version is for use in .dylib or .so libraries
#
# Follows the same logic from the Makefile commit c2e93911
# as non-numeric versions get bumped to patch level 255 to indicate
# an unknown value.
if [ "$1" = "numeric" ]
then
    v1=`expr "$VERSION" : '\([0-9]*\)'`
    v2=`expr "$VERSION" : '[0-9]*.\([0-9]*\)'`
    v3=`expr "$VERSION" : '[0-9]*.[0-9]*.\([0-9]*\)'`
    if [ -z "`expr "$VERSION" : '^\([0-9.]*\)$'`" ]
    then
        VERSION="$v1.$v2.255"
    else
        VERSION="$v1.$v2${v3:+.}$v3"
    fi
fi

echo $VERSION
