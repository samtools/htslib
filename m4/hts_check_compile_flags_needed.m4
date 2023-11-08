# hts_check_compile_flags_needed.m4
#
# SYNOPSIS
#
#   HTS_CHECK_COMPILE_FLAGS_NEEDED(FEATURE, FLAGS, [INPUT], [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS])
#
# DESCRIPTION
#
#   Check whether the given FLAGS are required to build and link INPUT with
#   the current language's compiler.  Compilation and linking are first
#   tries without FLAGS.  If that fails it then tries to compile and
#   link again with FLAGS.
#
#   FEATURE describes the feature being tested, and is used when printing
#   messages and to name the cache entry (along with the tested flags).
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.  In ACTION-SUCCESS, $flags_needed will be set to
#   either an empty string or FLAGS depending on the test results.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   If omitted, INPUT defaults to AC_LANG_PROGRAM(), although that probably
#   isn't very useful.
#
#   NOTE: Implementation based on AX_CHECK_COMPILE_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#   Copyright (c) 2023 Robert Davies <rmd@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#   HTS_CHECK_COMPILE_FLAGS_NEEDED(FEATURE, FLAGS, [INPUT], [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS])

AC_DEFUN([HTS_CHECK_COMPILE_FLAGS_NEEDED],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[hts_cv_check_[]_AC_LANG_ABBREV[]flags_needed_$1_$6_$2])dnl
AC_CACHE_CHECK([_AC_LANG compiler flags needed for $1], CACHEVAR, [
  AC_LINK_IFELSE([m4_default([$3],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[none])],
    [ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
     _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $6 $2"
     AC_LINK_IFELSE([m4_default([$3],[AC_LANG_PROGRAM()])],
       [AS_VAR_SET(CACHEVAR,["$2"])],
       [AS_VAR_SET(CACHEVAR,[unsupported])])
     _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])])
AS_VAR_IF(CACHEVAR,unsupported, [
  m4_default([$5], :)
], [
  AS_VAR_IF(CACHEVAR,none,[flags_needed=""], [flags_needed="$CACHEVAR"])
  m4_default([$4], :)
])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl HTS_CHECK_COMPILE_FLAGS_NEEDED
