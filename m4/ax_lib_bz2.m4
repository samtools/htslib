# SYNOPSIS
#
#   AX_LIB_BZ2([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of bz2 library.
#   It does this by checking for the header file bzlib.h and the bz2 library
#   object file. The location of these may be specified using the
#   --with-bz2=DIR command line option (eg --with-bz2=/usr/local),
#   using $DIR/include and $DIR/lib for the search path.
#
#   The following output variables are set using AC_SUBST:
#
#     BZ2_VERSION (if MINIMUM-VERSION is not "")
#
#   The C preprocessor symbol HAVE_LIBBZ2 will be also defined with
#   AC_DEFINE if a functioning libbz2 is available.
#
# LICENSE
#
#   Copyright (C) 2010,2016 Genome Research Ltd.
#   Author: James Bonfield <jkb@sanger.ac.uk>
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.



AC_DEFUN([AX_LIB_BZ2],
[
  AC_ARG_WITH(bz2,
              AC_HELP_STRING([--with-bz2=DIR],[look for libbz2 in DIR]),
              [_bz2_with=$withval],[_bz2_with=""])

  BZ2_ROOT=""
  if test "$_bz2_with" != "no"
  then
     if test -f "$_bz2_with/include/bzlib.h"
     then
         BZ2_ROOT=$_bz2_with
     fi
  
    # Check if it's a working library
    bz2_ok=no
    if test "x$BZ2_ROOT" != "x"
    then
      _cppflags=$CPPFLAGS
      CPPFLAGS="$CPPFLAGS -I${BZ2_ROOT}/include"
      _ldflags=$LDFLAGS
      LDFLAGS="$LFDLAGS -L${BZ2_ROOT}/lib"
      AC_LANG_PUSH([C])
      AC_CHECK_LIB(bz2, BZ2_bzBuffToBuffCompress,
          [AC_CHECK_HEADER(bzlib.h, bz2_ok=yes, bz2_ok=no)])
      AC_LANG_POP([C])
      if test "$bz2_ok" != "yes"
      then
          # Backout and whinge
          CPPFLAGS=$_cppflags
          LDFLAGS=$_ldflags
          AC_MSG_WARN(["--with-bz2 specified, but non functioning"])
      fi
  
    else
      # Maybe it works "out of the box"?
      AC_CHECK_LIB(bz2, BZ2_bzBuffToBuffCompress,
          [AC_CHECK_HEADER(bzlib.h, bz2_ok=yes, bz2_ok=no)])
    fi
  
    # Check version
    if test "x$1" != "x" && test "$bz2_ok" = "yes"
    then
        AC_MSG_CHECKING([if bz2 version >= $1])
  
        for i in "$BZ2_ROOT/include" "/usr/include" "/usr/share/include" "/usr/local/include"
        do
            if test -f "$i/bzlib.h"
            then
	        v1=`sed -n 's/.* version \([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\).*/\1/p' "$i/bzlib.h"`
	        v2=`sed -n 's/.* version \([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\).*/\2/p' "$i/bzlib.h"`
	        v3=`sed -n 's/.* version \([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\).*/\3/p' "$i/bzlib.h"`
                BZ2_VERSION=$v1.$v2.$v3
                break
            fi
        done
  
        have_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
        v1=`expr "$1" : '\([[0-9]]*\)'`
        v2=`expr "$1" : '[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$1" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        want_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
  
        if test `expr "$have_vers" ">=" "$want_vers"` = "1"
        then
            AC_MSG_RESULT([yes])
            AC_SUBST([BZ2_VERSION])
        else
	    AC_MSG_WARN([Version is too old])
            AC_MSG_RESULT([no])
            bz2_ok="no"
        fi
    fi
  
    # perform substitutions
    if test "$bz2_ok" = "yes"
    then
        AC_DEFINE(HAVE_LIBBZ2, 1,
           [Define to 1 if you have a functional libbz2.])
        if test "$BZ2_ROOT" != ""
        then
            CPPFLAGS="$CPPFLAGS -I${BZ2_ROOT}/include"
            LDFLAGS="$LDFLAGS -L${BZ2_ROOT}/lib"
            LIBS="-lbz2 $LIBS"
        else
            LIBS="-lbz2 $LIBS"
        fi
    else
      AC_MSG_WARN("No functioning bz2 found")
    fi
  
    # Not sure how many of these are needed, but it's belt-and-braces mode
    AH_TEMPLATE([HAVE_LIBBZ2], [Define if libbz2 is installed])
    AM_CONDITIONAL(HAVE_LIBBZ2, test "$bz2_ok" = "yes")
  
  
    # Execute the conditional expressions
    if test "$bz2_ok" = "yes"
    then
       # This is the IF-YES path
       ifelse([$2],,:,[$2])
    else
       # This is the IF-NO path
       ifelse([$3],,:,[$3])
    fi
  
    # Tidy up
    unset bz2_ok
    unset _cppflags
    unset _ldflags
  else
    AC_MSG_WARN([No bz2 support enabled.  Htslib will be unable to read/write some CRAM files.])
    AH_TEMPLATE([HAVE_LIBBZ2], [Define if libbz2 is installed])
    AM_CONDITIONAL(HAVE_LIBBZ2, false)
  fi
])
