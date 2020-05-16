#!/bin/sh -e
# test/with-shlib.sh -- make shared libhts available via $LD_LIBRARY_PATH etc.
#
#    Copyright (C) 2020 University of Glasgow.
#
#    Author: John Marshall <John.W.Marshall@glasgow.ac.uk>
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

libdir=${0%/*}/libdir-$$.tmp
case $libdir in
/*) abslibdir=$libdir ;;
*)  abslibdir=$PWD/$libdir ;;
esac

# Create a directory containing *only* the shared libhts, and add it
# to the platform-appropriate $LD_LIBRARY_PATH environment variable.

mkdir $libdir

case `uname -s` in
Darwin)
    (cd $libdir; ln -s ../../libhts.*.dylib .)
    export DYLD_LIBRARY_PATH=$abslibdir${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}
    ;;

*CYGWIN*)
    (cd $libdir; ln -s ../../cyghts-*.dll .)
    export PATH="$abslibdir${PATH:+;$PATH}"
    ;;

*MSYS*|*MINGW*)
    (cd $libdir; cp -p ../../hts-*.dll .)
    export PATH="$abslibdir${PATH:+;$PATH}"
    ;;

*)
    (cd $libdir; ln -s ../../libhts.so.* .)
    export LD_LIBRARY_PATH=$abslibdir${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
    ;;
esac

status=0
"$@" || status=$?

rm $libdir/*hts*
rmdir $libdir

exit $status
