# SYNOPSIS
#
#   AX_WITH_LIBDEFLATE([ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro checks whether Libdeflate is installed and adds a
# --with-libdeflate=DIR option to override the search path.
# See https://github.com/ebiggers/libdeflate for the library itself.
#
#   The following output variables are amended by this macro:
#
#     CPPFLAGS     Preprocessor flags for compiling
#     LDFLAGS      Linker flags for linking against the library
#     LIBS         Library list
#
#   It also sets LIBDEFLATE_LDFLAGS variable, to aid creation of
#   pkg-config files.
#
#   The HAVE_LIBDEFLATE cpp variable will be defined in a working
#   libdeflate was found.
#
# LICENSE
#
#   Copyright (C) 2018 Genome Research Ltd
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
AC_DEFUN([AX_LIBDEFLATE],
[
  AC_ARG_WITH(libdeflate,
	      AC_HELP_STRING([--with-libdeflate=DIR],[look for libdeflate in DIR]),
	      [_libdeflate_with=$withval],[_libdeflate_with="no"])

  # Check if it's a working library
  libdeflate_ok=no
  _cppflags=$CPPFLAGS
  _ldflags=$LDFLAGS
  if test "x$_libdeflate_with" != "xno"
  then
    if test "$_libdeflate_with" != "yes"
    then
      if test -f "${_libdeflate_with}/include/libdeflate.h"
      then
        CPPFLAGS="$CPPFLAGS -I${_libdeflate_with}/include"
      else
        CPPFLAGS="$CPPFLAGS -I${_libdeflate_with}"
      fi
      if test -f "${_libdeflate_with}/lib/libdeflate.a" -o -f "${_libdeflate_with}/lib/libdeflate.so"
      then
        LIBDEFLATE_LDFLAGS="-L${_libdeflate_with}/lib"
      else
        LIBDEFLATE_LDFLAGS="-L${_libdeflate_with}"
      fi
      LDFLAGS="$LDFLAGS ${LIBDEFLATE_LDFLAGS}"
    fi
    AC_SEARCH_LIBS([libdeflate_deflate_compress], [deflate],
	[AC_CHECK_HEADER(libdeflate.h, [libdeflate_ok=yes LIBS="$LIBS -ldeflate"], libdeflate_ok=no)])
    if test "$libdeflate_ok" != "yes"
    then
        AC_MSG_WARN("--with-libdeflate specified, but non functioning")
    fi

    # perform substitutions
    if test "$libdeflate_ok" = "yes"
    then
        AC_DEFINE(HAVE_LIBDEFLATE, 1,
           [Define to 1 if you have a functional libz.])
	LIBDEFLATE_LDFLAGS="$LIBDEFLATE_LDFLAGS $ac_cv_search_libdeflate_deflate_compress"
    else
      AC_MSG_WARN("No functioning libdeflate found")
      CPPFLAGS=$_cppflags
      LDFLAGS=$_ldflags
    fi
  fi

  AH_TEMPLATE([HAVE_LIBDEFLATE], [Define if libdeflate is installed])
  AM_CONDITIONAL(HAVE_LIBDEFLATE, test "$libdeflate_ok" = "yes")

  # Execute the conditional expressions
  if test "$libdeflate_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$1],,:,[$1])
  else
     # This is the IF-NO path
     ifelse([$2],,:,[$2])
  fi

  # Tidy up
  unset libdeflate_ok
  unset _cppflags
  unset _ldflags
])
