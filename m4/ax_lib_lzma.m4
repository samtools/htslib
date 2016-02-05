# SYNOPSIS
#
#   AX_LIB_LZMA([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of lzma library.
#   It does this by checking for the header file lzma.h and the lzma library
#   object file. The location of these may be specified using the
#   --with-lzma=DIR command line option (eg --with-lzma=/usr/local),
#   using $DIR/include and $DIR/lib for the search path.
#
#   The following output variables are set using AC_SUBST:
#
#     LZMA_VERSION (if MINIMUM-VERSION is not "")
#
#   The C preprocessor symbol HAVE_LIBLZMA will be also defined with
#   AC_DEFINE if a functioning liblzma is available.
#
# LICENSE
#
#   Copyright (c) 2010,2015 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_LZMA],
[
  AC_ARG_WITH(lzma,
              AC_HELP_STRING([--with-lzma=DIR],[look for liblzma in DIR]),
              [_lzma_with=$withval],[_lzma_with=""])

  LZMA_ROOT=""
  if test "$_lzma_with" != "no"
  then
     if test -f "$_lzma_with/include/lzma.h"
     then
         LZMA_ROOT=$_lzma_with
     fi
  
    # Check if it's a working library
    lzma_ok=no
    if test "x$LZMA_ROOT" != "x"
    then
      _cppflags=$CPPFLAGS
      CPPFLAGS="$CPPFLAGS -I${LZMA_ROOT}/include"
      _ldflags=$LDFLAGS
      LDFLAGS="$LFDLAGS -L${LZMA_ROOT}/lib"
      AC_LANG_PUSH([C])
      AC_CHECK_LIB(lzma, lzma_easy_buffer_encode,
          [AC_CHECK_HEADER(lzma.h, lzma_ok=yes, lzma_ok=no)])
      AC_LANG_POP([C])
      if test "$lzma_ok" != "yes"
      then
          # Backout and whinge
          CPPFLAGS=$_cppflags
          LDFLAGS=$_ldflags
          AC_MSG_ERROR(["--with-lzma specified, but non functioning"])
      fi
  
    else
      # Maybe it works "out of the box"?
      AC_CHECK_LIB(lzma, lzma_easy_buffer_encode,
          [AC_CHECK_HEADER(lzma.h, lzma_ok=yes, lzma_ok=no)])
    fi
  
    # Check version
    if test "x$1" != "x" && test "$lzma_ok" = "yes"
    then
        AC_MSG_CHECKING([if lzma version >= $1])
  
        for i in "$LZMA_ROOT/include" "/usr/include" "/usr/share/include" "/usr/local/include"
        do
            if test -f "$i/lzma/version.h"
            then
                v1=`sed -n 's/.*#define *LZMA_VERSION_MAJOR *\(.*\)/\1/p' "$i/lzma/version.h"`
                v2=`sed -n 's/.*#define *LZMA_VERSION_MINOR *\(.*\)/\1/p' "$i/lzma/version.h"`
                v3=`sed -n 's/.*#define *LZMA_VERSION_PATCH *\(.*\)/\1/p' "$i/lzma/version.h"`
                LZMA_VERSION=$v1.$v2.$v3
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
            AC_SUBST([LZMA_VERSION])
        else
	    AC_MSG_WARN([Version is too old])
            AC_MSG_RESULT([no])
            lzma_ok="no"
        fi
    fi
  
    # perform substitutions
    if test "$lzma_ok" = "yes"
    then
        AC_DEFINE(HAVE_LIBLZMA, 1,
           [Define to 1 if you have a functional liblzma.])
        if test "$LZMA_ROOT" != ""
        then
            CPPFLAGS="$CPPFLAGS -I${LZMA_ROOT}/include"
            LDFLAGS="$LDFLAGS -L${LZMA_ROOT}/lib"
            LIBS="-llzma $LIBS"
        else
            LIBS="-llzma $LIBS"
        fi
    else
      AC_MSG_ERROR("No functioning lzma found")
    fi
  
    # Not sure how many of these are needed, but it's belt-and-braces mode
    AH_TEMPLATE([HAVE_LIBLZMA], [Define if liblzma is installed])
    AM_CONDITIONAL(HAVE_LIBLZMA, test "$lzma_ok" = "yes")
  
  
    # Execute the conditional expressions
    if test "$lzma_ok" = "yes"
    then
       # This is the IF-YES path
       ifelse([$2],,:,[$2])
    else
       # This is the IF-NO path
       ifelse([$3],,:,[$3])
    fi
  
    # Tidy up
    unset lzma_ok
    unset _cppflags
    unset _ldflags
  else
    AC_MSG_WARN([No lzma support enabled.  Htslib will be unable to read/write some CRAM files.])
    AH_TEMPLATE([HAVE_LIBLZMA], [Define if liblzma is installed])
    AM_CONDITIONAL(HAVE_LIBLZMA, false)
  fi
])
