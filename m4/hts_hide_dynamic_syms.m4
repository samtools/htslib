dnl @synopsis HTS_HIDE_DYNAMIC_SYMBOLS
dnl
dnl Turn on compiler options that prevent unwanted symbols from being exported
dnl by shared libraries.
dnl
dnl @author Rob Davies <rmd@sanger.ac.uk>
dnl @license MIT/Expat
dnl
dnl Copyright (C) 2018 Genome Research Ltd.
dnl
dnl Permission is hereby granted, free of charge, to any person obtaining a copy
dnl of this software and associated documentation files (the "Software"), to
dnl deal in the Software without restriction, including without limitation the
dnl rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
dnl  sell copies of the Software, and to permit persons to whom the Software is
dnl furnished to do so, subject to the following conditions:
dnl
dnl The above copyright notice and this permission notice shall be included in
dnl all copies or substantial portions of the Software.
dnl
dnl THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
dnl IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
dnl FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
dnl THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
dnl LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
dnl FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
dnl DEALINGS IN THE SOFTWARE.

# SYNOPSIS
#
# HTS_TEST_CC_C_LD_FLAG(FLAG, FOUND_VAR)
#
# Test if FLAG can be used on both CFLAGS and LDFLAGS.  It it works,
# variable FOUND_VAR is set to FLAG.

AC_DEFUN([HTS_TEST_CC_C_LD_FLAG],
 [AS_VAR_PUSHDEF([hts_cv_check_flag],[hts_cv_check_$1])dnl
  AC_CACHE_CHECK([whether the compiler accepts $1],
   [hts_cv_check_flag],
   [ac_check_save_cflags=$CFLAGS
    ac_check_save_ldflags=$LDFLAGS
    CFLAGS="$CFLAGS $1"
    LDFLAGS="$LDFLAGS $1"
    AC_LINK_IFELSE([AC_LANG_PROGRAM()],
      [AS_VAR_SET([hts_cv_check_flag],[yes])
       AS_IF([test "x$2" != x],[eval AS_TR_SH([$2])="$1"])],
      [AS_VAR_SET([hts_cv_check_flag],[no])])
    CFLAGS=$ac_check_save_cflags
    LDFLAGS=$ac_check_save_ldflags])
  AS_VAR_POPDEF([hts_cv_check_flag])dnl
])

AC_DEFUN([HTS_HIDE_DYNAMIC_SYMBOLS], [
  # Test for flags to set default shared library visibility to hidden
  # -fvisibility=hidden : GCC compatible
  # -xldscope=hidden    : SunStudio
  ac_opt_found=no
  m4_foreach_w([ac_opt],[-fvisibility=hidden -xldscope=hidden],
   [AS_IF([test "x$ac_opt_found" = "xno"],
     [HTS_TEST_CC_C_LD_FLAG(ac_opt,[ac_opt_found])])
   ])
  AS_IF([test "x$ac_opt_found" != "xno"],
   [CFLAGS="$CFLAGS $ac_opt_found"
    LDFLAGS="$LDFLAGS $ac_opt_found"])
])
