# Configure script for htslib, a C library for high-throughput sequencing data.
#
#    Copyright (C) 2015-2024 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
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

dnl Process this file with autoconf to produce a configure script
AC_INIT([HTSlib], m4_esyscmd_s([./version.sh 2>/dev/null]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ(2.63)  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR(hts.c)
AC_CONFIG_HEADERS(config.h)

m4_include([m4/hts_prog_cc_warnings.m4])
m4_include([m4/hts_check_compile_flags_needed.m4])
m4_include([m4/hts_hide_dynamic_syms.m4])
m4_include([m4/pkg.m4])

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2020-2024 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

dnl Notes to be copied (by autoheader) into the generated config.h.in
AH_TOP([/* If you use configure, this file provides @%:@defines reflecting your
   configuration choices.  If you have not run configure, suitable
   conservative defaults will be used.

   Autoheader adds a number of items to this template file that are not
   used by HTSlib: STDC_HEADERS and most HAVE_*_H header file defines
   are immaterial, as we assume standard ISO C headers and facilities;
   the PACKAGE_* defines are unused and are overridden by the more
   accurate PACKAGE_VERSION as computed by the Makefile.  */])

dnl Variant of AC_MSG_ERROR that ensures subsequent make(1) invocations fail
dnl until the configuration error is resolved and configure is run again.
AC_DEFUN([MSG_ERROR],
  [cat > config.mk <<'EOF'
ifneq ($(MAKECMDGOALS),distclean)
$(error Resolve configure error first)
endif
EOF
   AC_MSG_ERROR([$1], [$2])])

AC_PROG_CC
AC_PROG_RANLIB

dnl Turn on compiler warnings, if possible
HTS_PROG_CC_WARNINGS
dnl Flags to treat warnings as errors.  These need to be applied to CFLAGS
dnl later as they can interfere with some of the tests (notably AC_SEARCH_LIBS)
HTS_PROG_CC_WERROR(hts_late_cflags)

# HTSlib uses X/Open-only facilities (M_SQRT2 etc, drand48() etc), and
# various POSIX functions that are provided by various _POSIX_C_SOURCE values
# or by _XOPEN_SOURCE >= 500. It also uses usleep(), which is removed when
# _XOPEN_SOURCE >= 700. Additionally, some definitions may require
# _XOPEN_SOURCE >= 600 on some platforms (snprintf on MinGW,
# PTHREAD_MUTEX_RECURSIVE on some Linux distributions). Hence we set it to 600.

# Define _XOPEN_SOURCE unless the user has already done so via $CPPFLAGS etc.
AC_CHECK_DECL([_XOPEN_SOURCE], [],
  [AC_DEFINE([_XOPEN_SOURCE], [600], [Specify X/Open requirements])],
  [])

dnl Check that we have cpuid, and if so run the x86 SIMD checks
AC_CHECK_DECLS([__get_cpuid_max, __cpuid_count], [
   hts_have_cpuid=yes
], [
   hts_have_cpuid=no
], [[#include <cpuid.h>]])

AS_IF(test "x$hts_have_cpuid" = "xyes", [
dnl Options for rANS32x16 sse4.1 version - sse4.1
HTS_CHECK_COMPILE_FLAGS_NEEDED([sse4.1], [-msse4.1 -mssse3 -mpopcnt],
 [AC_LANG_PROGRAM([[
    #ifdef __x86_64__
    #include "x86intrin.h"
    #endif
  ]],[[
    #ifdef __x86_64__
    __m128i a = _mm_set_epi32(1, 2, 3, 4), b = _mm_set_epi32(4, 3, 2, 1);
    __m128i c = _mm_shuffle_epi8(_mm_max_epu32(a, b), b);
    return _mm_popcnt_u32(*((char *) &c));
    #endif
  ]])], [
  hts_cflags_sse4="$flags_needed"
  AC_DEFINE([HAVE_SSSE3],1,[Defined to 1 if rANS source using SSSE3 can be compiled.])
  AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
  AC_DEFINE([HAVE_SSE4_1],1,[Defined to 1 if rANS source using SSE4.1 can be compiled.
])

dnl Propagate HTSlib's unaligned access preference to htscodecs
  AH_VERBATIM([UBSAN],[
/* Prevent unaligned access in htscodecs SSE4 rANS codec */
#if defined(HTS_ALLOW_UNALIGNED) && HTS_ALLOW_UNALIGNED == 0
#undef UBSAN
#endif])
  AC_DEFINE([UBSAN],1,[])
])
AC_SUBST([hts_cflags_sse4])

dnl Options for rANS32x16 avx2 version
HTS_CHECK_COMPILE_FLAGS_NEEDED([avx2], [-mavx2 -mpopcnt], [AC_LANG_PROGRAM([[
    #ifdef __x86_64__
    #include "x86intrin.h"
    #endif
  ]],[[
    #ifdef __x86_64__
    __m256i a = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
    __m256i b = _mm256_add_epi32(a, a);
    long long c = _mm256_extract_epi64(b, 0);
    return _mm_popcnt_u32((int) c);
    #endif
  ]])], [
  hts_cflags_avx2="$flags_needed"
  AC_SUBST([hts_cflags_avx2])
  AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
  AC_DEFINE([HAVE_AVX2],1,[Defined to 1 if rANS source using AVX2 can be compiled.])
])

dnl Options for rANS32x16 avx512 version
HTS_CHECK_COMPILE_FLAGS_NEEDED([avx512f], [-mavx512f -mpopcnt],
 [AC_LANG_PROGRAM([[
    #ifdef __x86_64__
    #include "x86intrin.h"
    #endif
  ]],[[
    #ifdef __x86_64__
    __m512i a = _mm512_set1_epi32(1);
    __m512i b = _mm512_add_epi32(a, a);
    __m256i c = _mm512_castsi512_si256(b);
    __m256i d = _mm512_extracti64x4_epi64(a, 1);
    return _mm_popcnt_u32(*((char *) &c)) + (*(char *) &d);
    #endif
  ]])], [
  hts_cflags_avx512="$flags_needed"
  AC_SUBST([hts_cflags_avx512])
  AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
  AC_DEFINE([HAVE_AVX512],1,[Defined to 1 if rANS source using AVX512F can be compiled.])
])

dnl Check for working __builtin_cpu_supports (ssse3 is broken on some clangs)
AC_MSG_CHECKING([for working __builtin_cpu_supports("ssse3")])
AC_LINK_IFELSE([AC_LANG_PROGRAM([],[
  if (__builtin_cpu_supports("ssse3")) {
    return 0;
  }
])], [
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_BUILTIN_CPU_SUPPORT_SSSE3], 1,
            [Defined to 1 if __builtin_cpu_supports("ssse3") works])
], [
  AC_MSG_RESULT([no])
])

dnl Check for function attribute used in conjunction with __builtin_cpu_supports
AC_MSG_CHECKING([for __attribute__((target))])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  __attribute__((target("ssse3")))
  int zero(void) {
    return 0;
  }
]], [[zero();]])], [
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_ATTRIBUTE_TARGET], 1,
            [Define if __attribute__((target(...))) is available.])
], [
  AC_MSG_RESULT([no])
])

]) dnl End of AS_IF(hts_have_cpuid)

dnl Avoid chicken-and-egg problem where pkg-config supplies the
dnl PKG_PROG_PKG_CONFIG macro, but we want to use it to check
dnl for pkg-config...
m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG], [PKG_CONFIG=""])

need_crypto=no
pc_requires=
static_LDFLAGS=$LDFLAGS
static_LIBS='-lpthread -lz -lm'
private_LIBS=$LDFLAGS

AC_ARG_ENABLE([versioned-symbols],
  [AS_HELP_STRING([--disable-versioned-symbols],
                  [disable versioned symbols in shared library])],
  [], [enable_versioned_symbols=yes])

AC_ARG_ENABLE([bz2],
  [AS_HELP_STRING([--disable-bz2],
                  [omit support for BZ2-compressed CRAM files])],
  [], [enable_bz2=yes])

AC_ARG_ENABLE([gcs],
  [AS_HELP_STRING([--enable-gcs],
                  [support Google Cloud Storage URLs])],
  [], [enable_gcs=check])

AC_SYS_LARGEFILE

AC_ARG_ENABLE([libcurl],
  [AS_HELP_STRING([--enable-libcurl],
                  [enable libcurl-based support for http/https/etc URLs])],
  [], [enable_libcurl=check])

AC_ARG_ENABLE([lzma],
  [AS_HELP_STRING([--disable-lzma],
                  [omit support for LZMA-compressed CRAM files])],
  [], [enable_lzma=yes])

AC_ARG_ENABLE([plugins],
  [AS_HELP_STRING([--enable-plugins],
                  [enable separately-compiled plugins for file access])],
  [], [enable_plugins=no])
AC_SUBST(enable_plugins)

AC_ARG_WITH([external-htscodecs],
  [AS_HELP_STRING([--with-external-htscodecs],
                  [get htscodecs functions from a shared library])],
  [], [with_external_htscodecs=no])
AC_SUBST(with_external_htscodecs)

AC_ARG_WITH([libdeflate],
  [AS_HELP_STRING([--with-libdeflate],
                  [use libdeflate for faster crc and deflate algorithms])],
  [], [with_libdeflate=check])

AC_ARG_WITH([plugin-dir],
  [AS_HELP_STRING([--with-plugin-dir=DIR],
                  [plugin installation location [LIBEXECDIR/htslib]])],
  [case $withval in
     yes|no) MSG_ERROR([no directory specified for --with-plugin-dir]) ;;
   esac],
   [with_plugin_dir='$(libexecdir)/htslib'])
AC_SUBST([plugindir], $with_plugin_dir)

AC_ARG_WITH([plugin-path],
  [AS_HELP_STRING([--with-plugin-path=PATH],
                  [default HTS_PATH plugin search path [PLUGINDIR]])],
  [case $withval in
     yes) MSG_ERROR([no path specified for --with-plugin-path]) ;;
     no)  with_plugin_path= ;;
   esac],
  [with_plugin_path=$with_plugin_dir])
AC_SUBST([pluginpath], $with_plugin_path)

AC_ARG_ENABLE([s3],
  [AS_HELP_STRING([--enable-s3],
                  [support Amazon AWS S3 URLs])],
  [], [enable_s3=check])

basic_host=${host_alias:-unknown-`uname -s`}
AC_MSG_CHECKING([shared library type for $basic_host])
case $basic_host in
  *-cygwin* | *-CYGWIN*)
    host_result="Cygwin DLL"
    PLATFORM=CYGWIN
    PLUGIN_EXT=.cygdll
    ;;
  *-darwin* | *-Darwin*)
    host_result="Darwin dylib"
    PLATFORM=Darwin
    PLUGIN_EXT=.bundle
    ;;
  *-msys* | *-MSYS* | *-mingw* | *-MINGW*)
    host_result="MSYS dll"
    PLATFORM=MSYS
    PLUGIN_EXT=.dll
    # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
    # %lld and %z printf formats work.  It also enforces the snprintf to
    # be C99 compliant so it returns the correct values (in kstring.c).

    # Now set by default, so no need to do it here.
    # CPPFLAGS="$CPPFLAGS -D_XOPEN_SOURCE=600"
    ;;
  *)
    host_result="plain .so"
    PLATFORM=default
    PLUGIN_EXT=.so
    ;;
esac
AC_MSG_RESULT([$host_result])
AC_SUBST([PLATFORM])

dnl Check for versioned symbol support
dnl Only try for .so shared libraries as other types won't work
AS_IF([test x"$PLATFORM" = xdefault && test x"$enable_versioned_symbols" = xyes],
  [AC_CACHE_CHECK([whether the linker supports versioned symbols],
    [hts_cv_have_versioned_symbols], [
      save_LDFLAGS=$LDFLAGS
      LDFLAGS="-Wl,-version-script,$srcdir/htslib.map $LDFLAGS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM()],
                     [hts_cv_have_versioned_symbols=yes],
                     [hts_cv_have_versioned_symbols=no])
      LDFLAGS=$save_LDFLAGS
    ])
   AS_IF([test "x$hts_cv_have_versioned_symbols" = xyes],[
     VERSION_SCRIPT_LDFLAGS='-Wl,-version-script,$(srcprefix)htslib.map'
     AC_SUBST([VERSION_SCRIPT_LDFLAGS])
   ])
])

dnl Try to get more control over which symbols are exported in the shared
dnl library.
HTS_HIDE_DYNAMIC_SYMBOLS

dnl FIXME This pulls in dozens of standard header checks
AC_FUNC_MMAP
AC_CHECK_FUNCS([gmtime_r fsync drand48 srand48_deterministic])

# Darwin has a dubious fdatasync() symbol, but no declaration in <unistd.h>
AC_CHECK_DECL([fdatasync(int)], [AC_CHECK_FUNCS(fdatasync)])

AC_MSG_CHECKING([for __attribute__((constructor))])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
  static __attribute__((constructor)) void noop(void) {}
]], [])], [
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_ATTRIBUTE_CONSTRUCTOR], 1,
            [Define if __attribute__((constructor)) is available.])
], [AC_MSG_RESULT([no])])

AC_MSG_CHECKING([for clock_gettime with CLOCK_PROCESS_CPUTIME_ID])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <time.h>]], [[
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
]])], [
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_CLOCK_GETTIME_CPUTIME], 1,
            [Define if clock_gettime exists and accepts CLOCK_PROCESS_CPUTIME_ID.])
], [AC_MSG_RESULT([no])])

if test $enable_plugins != no; then
  AC_SEARCH_LIBS([dlsym], [dl], [],
    [MSG_ERROR([dlsym() not found

Plugin support requires dynamic linking facilities from the operating system.
Either configure with --disable-plugins or resolve this error to build HTSlib.])])
  # Check if the compiler understands -rdynamic 
  # TODO Test whether this is required and/or needs tweaking per-platform
  HTS_TEST_CC_C_LD_FLAG([-rdynamic],[rdynamic_flag])
  AS_IF([test x"$rdynamic_flag" != "xno"],
    [LDFLAGS="$LDFLAGS $rdynamic_flag"
     static_LDFLAGS="$static_LDFLAGS $rdynamic_flag"])
  case "$ac_cv_search_dlsym" in
    -l*) static_LIBS="$static_LIBS $ac_cv_search_dlsym" ;;
  esac
  AC_DEFINE([ENABLE_PLUGINS], 1, [Define if HTSlib should enable plugins.])
  AC_SUBST([PLUGIN_EXT])
  AC_DEFINE_UNQUOTED([PLUGIN_EXT], ["$PLUGIN_EXT"],
                     [Platform-dependent plugin filename extension.])
fi

AC_SEARCH_LIBS([log], [m], [],
  [MSG_ERROR([log() not found

HTSLIB requires a working floating-point math library.
FAILED.  This error must be resolved in order to build HTSlib successfully.])])

zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])

if test $zlib_devel != ok; then
  MSG_ERROR([zlib development files not found

HTSlib uses compression routines from the zlib library <http://zlib.net>.
Building HTSlib requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions or Cygwin)
is installed.

FAILED.  This error must be resolved in order to build HTSlib successfully.])
fi

dnl connect() etc. fns are in libc on linux, but libsocket on illumos/Solaris
AC_SEARCH_LIBS([recv], [socket ws2_32], [
if test "$ac_cv_search_recv" != "none required"
then
  static_LIBS="$static_LIBS $ac_cv_search_recv"
fi],
  dnl on MinGW-i686, checking recv() linking requires an annotated declaration
  [AC_MSG_CHECKING([for library containing recv using declaration])
   LIBS="-lws2_32 $LIBS"
   AC_LINK_IFELSE(
     [AC_LANG_PROGRAM([[#include <winsock2.h>]], [[recv(0, 0, 0, 0);]])],
     [AC_MSG_RESULT([-lws2_32])
      static_LIBS="$static_LIBS -lws2_32"],
     [AC_MSG_RESULT([no])
      MSG_ERROR([unable to find the recv() function])])])

if test "$enable_bz2" != no; then
  bz2_devel=ok
  AC_CHECK_HEADER([bzlib.h], [], [bz2_devel=missing], [;])
  AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [], [bz2_devel=missing])
  if test $bz2_devel != ok; then
    MSG_ERROR([libbzip2 development files not found

The CRAM format may use bzip2 compression, which is implemented in HTSlib
by using compression routines from libbzip2 <http://www.bzip.org/>.

Building HTSlib requires libbzip2 development files to be installed on the
build machine; you may need to ensure a package such as libbz2-dev (on Debian
or Ubuntu Linux) or bzip2-devel (on RPM-based Linux distributions or Cygwin)
is installed.

Either configure with --disable-bz2 (which will make some CRAM files
produced elsewhere unreadable) or resolve this error to build HTSlib.])
  fi
dnl Unfortunately the 'bzip2' package-cfg module is not standard.
dnl Redhat/Fedora has it; Debian/Ubuntu does not.
  if test -n "$PKG_CONFIG" && "$PKG_CONFIG" --exists bzip2; then
     pc_requires="$pc_requires bzip2"
  else
     private_LIBS="$private_LIBS -lbz2"
  fi
  static_LIBS="$static_LIBS -lbz2"
fi

if test "$enable_lzma" != no; then
  lzma_devel=ok
  AC_CHECK_HEADERS([lzma.h], [], [lzma_devel=header-missing], [;])
  AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [], [lzma_devel=missing])
  if test $lzma_devel = missing; then
    MSG_ERROR([liblzma development files not found

The CRAM format may use LZMA2 compression, which is implemented in HTSlib
by using compression routines from liblzma <http://tukaani.org/xz/>.

Building HTSlib requires liblzma development files to be installed on the
build machine; you may need to ensure a package such as liblzma-dev (on Debian
or Ubuntu Linux), xz-devel (on RPM-based Linux distributions or Cygwin), or
xz (via Homebrew on macOS) is installed; or build XZ Utils from source.

Either configure with --disable-lzma (which will make some CRAM files
produced elsewhere unreadable) or resolve this error to build HTSlib.])
  fi
  pc_requires="$pc_requires liblzma"
  static_LIBS="$static_LIBS -llzma"
fi

AS_IF([test "x$with_external_htscodecs" != "xno"],
  [libhtscodecs=ok
   AC_CHECK_HEADER([htscodecs/rANS_static4x16.h],[],
                   [libhtscodecs='missing header'],[;])
   AC_CHECK_LIB([htscodecs],[rans_compress_bound_4x16],
                [:],[libhtscodecs='missing library'])
   AS_IF([test "$libhtscodecs" = "ok"],
     [AC_DEFINE([HAVE_EXTERNAL_LIBHTSCODECS], 1, [Define if using an external libhtscodecs])
      LIBS="-lhtscodecs $LIBS"
      private_LIBS="-lhtscodecs $private_LIBS"
      static_LIBS="-lhtscodecs $static_LIBS"
      selected_htscodecs_mk="htscodecs_external.mk"],
      [MSG_ERROR([libhtscodecs development files not found: $libhtscodecs

You asked to use an external htscodecs library, but do not have the
required header / library files.  You either need to supply these and
if necessary set CPPFLAGS and LDFLAGS so the compiler can find them;
or configure using --without-external-htscodecs to build the required
functions from the htscodecs submodule.
])])],
  [AC_MSG_CHECKING([whether htscodecs files are present])
   AS_IF([test -e "$srcdir/htscodecs/htscodecs/rANS_static4x16.h"],
     [AC_MSG_RESULT([yes])
      selected_htscodecs_mk="htscodecs_bundled.mk"],
     [AC_MSG_RESULT([no])
      AS_IF([test -e "$srcdir/.git"],
        [MSG_ERROR([htscodecs submodule files not present.

HTSlib uses some functions from the htscodecs project, which is normally
included as a submodule.  Try running:

  git submodule update --init --recursive

in  the top-level htslib directory to update it, and then re-run configure.
])],
        [MSG_ERROR([htscodecs submodule files not present.

You have an incomplete distribution.  Please try downloading one of the
official releases from https://www.htslib.org
])])])])

AS_IF([test "x$with_libdeflate" != "xno"],
  [libdeflate=ok
   AC_CHECK_HEADER([libdeflate.h],[],[libdeflate='missing header'],[;])
   AC_CHECK_LIB([deflate], [libdeflate_deflate_compress],[:],[libdeflate='missing library'])
   AS_IF([test "$libdeflate" = "ok"],
    [AC_DEFINE([HAVE_LIBDEFLATE], 1, [Define if libdeflate is available.])
     LIBS="-ldeflate $LIBS"
     private_LIBS="$private_LIBS -ldeflate"
     static_LIBS="$static_LIBS -ldeflate"],
    [AS_IF([test "x$with_libdeflate" != "xcheck"],
       [MSG_ERROR([libdeflate development files not found: $libdeflate

You requested libdeflate, but do not have the required header / library
files.  The source for libdeflate is available from
<https://github.com/ebiggers/libdeflate>.  You may have to adjust
search paths in CPPFLAGS and/or LDFLAGS if the header and library
are not currently on them.

Either configure with --without-libdeflate or resolve this error to build
HTSlib.])])])])

libcurl=disabled
if test "$enable_libcurl" != no; then
  libcurl_devel=ok
  AC_CHECK_HEADER([curl/curl.h], [], [libcurl_devel="headers not found"], [;])
  AC_CHECK_LIB([curl], [curl_easy_pause], [:],
    [AC_CHECK_LIB([curl], [curl_easy_init],
       [libcurl_devel="library is too old (7.18+ required)"],
       [libcurl_devel="library not found"])])

  if test "$libcurl_devel" = ok; then
    AC_DEFINE([HAVE_LIBCURL], 1, [Define if libcurl file access is enabled.])
    libcurl=enabled
  elif test "$enable_libcurl" = check; then
    AC_MSG_WARN([libcurl not enabled: $libcurl_devel])
  else
    MSG_ERROR([libcurl $libcurl_devel

Support for HTTPS and other SSL-based URLs requires routines from the libcurl
library <http://curl.se/libcurl/>.  Building HTSlib with libcurl enabled
requires libcurl development files to be installed on the build machine; you
may need to ensure a package such as libcurl4-{gnutls,nss,openssl}-dev (on
Debian or Ubuntu Linux) or libcurl-devel (on RPM-based Linux distributions
or Cygwin) is installed.

Either configure with --disable-libcurl or resolve this error to build HTSlib.])
  fi

dnl -lcurl is only needed for static linking if hfile_libcurl is not a plugin
  if test "$libcurl" = enabled ; then
    if test "$enable_plugins" != yes ; then
      static_LIBS="$static_LIBS -lcurl"
    fi
  fi
fi
AC_SUBST([libcurl])

gcs=disabled
if test "$enable_gcs" != no; then
  if test $libcurl = enabled; then
    AC_DEFINE([ENABLE_GCS], 1, [Define if HTSlib should enable GCS support.])
    gcs=enabled
  else
    case "$enable_gcs" in
      check) AC_MSG_WARN([GCS support not enabled: requires libcurl support]) ;;
      *) MSG_ERROR([GCS support not enabled

Support for Google Cloud Storage URLs requires libcurl support to be enabled
in HTSlib.  Configure with --enable-libcurl in order to use GCS URLs.])
      ;;
    esac
  fi
fi
AC_SUBST([gcs])

s3=disabled
if test "$enable_s3" != no; then
  if test $libcurl = enabled; then
    s3=enabled
    need_crypto="$enable_s3"
  else
    case "$enable_s3" in
      check) AC_MSG_WARN([S3 support not enabled: requires libcurl support]) ;;
      *) MSG_ERROR([S3 support not enabled

Support for Amazon AWS S3 URLs requires libcurl support to be enabled
in HTSlib.  Configure with --enable-libcurl in order to use S3 URLs.])
      ;;
    esac
  fi
fi

CRYPTO_LIBS=
if test $need_crypto != no; then
  AC_CHECK_FUNC([CCHmac],
    [AC_DEFINE([HAVE_COMMONCRYPTO], 1,
               [Define if you have the Common Crypto library.])],
    [save_LIBS=$LIBS
     AC_SEARCH_LIBS([HMAC], [crypto],
       [AC_DEFINE([HAVE_HMAC], 1, [Define if you have libcrypto-style HMAC().])
        case "$ac_cv_search_HMAC" in
          -l*) CRYPTO_LIBS=$ac_cv_search_HMAC ;;
        esac],
     [case "$need_crypto" in
     check) AC_MSG_WARN([S3 support not enabled: requires SSL development files])
         s3=disabled ;;
     *) MSG_ERROR([SSL development files not found

Support for AWS S3 URLs requires routines from an SSL library.  Building
HTSlib with libcurl enabled requires SSL development files to be installed
on the build machine; you may need to ensure a package such as libgnutls-dev,
libnss3-dev, or libssl-dev (on Debian or Ubuntu Linux, corresponding to the
libcurl4-*-dev package installed), or openssl-devel (on RPM-based Linux
distributions or Cygwin) is installed.

Either configure with --disable-s3 or resolve this error to build HTSlib.]) ;;
       esac])
     LIBS=$save_LIBS])
dnl Only need to add to static_LIBS if not building as a plugin
  if test "$enable_plugins" != yes ; then
     static_LIBS="$static_LIBS $CRYPTO_LIBS"
  fi
fi

dnl Look for regcomp in various libraries (needed on windows/mingw).
AC_SEARCH_LIBS(regcomp, regex, [libregex=needed], [])

dnl Look for PTHREAD_MUTEX_RECURSIVE.
dnl This is normally in pthread.h except on some broken glibc implementations.
dnl Now set by default
dnl AC_CHECK_DECL(PTHREAD_MUTEX_RECURSIVE, [], [AC_DEFINE([_XOPEN_SOURCE],[600], [Needed for PTHREAD_MUTEX_RECURSIVE])], [[#include <pthread.h>]])

if test "$s3" = enabled ; then
   AC_DEFINE([ENABLE_S3], 1, [Define if HTSlib should enable S3 support.])
fi

dnl Apply value from HTS_PROG_CC_WERROR (if set)
AS_IF([test "x$hts_late_cflags" != x],[CFLAGS="$CFLAGS $hts_late_cflags"])

AC_SUBST([s3])
AC_SUBST([CRYPTO_LIBS])

AC_SUBST([pc_requires])
AC_SUBST([private_LIBS])
AC_SUBST([static_LDFLAGS])
AC_SUBST([static_LIBS])

AC_CONFIG_FILES([config.mk htslib.pc.tmp:htslib.pc.in])
AC_CONFIG_LINKS([htscodecs.mk:$selected_htscodecs_mk])

if test "$srcdir" != .; then
  # Set up for a separate build directory. As HTSlib uses a non-recursive
  # makefile, we need to create additional build subdirectories explicitly.
  AC_CONFIG_LINKS([Makefile:Makefile htslib.mk:htslib.mk])
  AC_CONFIG_FILES([htslib_vars.mk:builddir_vars.mk.in])
  AC_CONFIG_COMMANDS([mkdir],
    [AS_MKDIR_P([cram])
     AS_MKDIR_P([htscodecs/htscodecs])
     AS_MKDIR_P([htscodecs/tests])
     AS_MKDIR_P([test/fuzz])
     AS_MKDIR_P([test/longrefs])
     AS_MKDIR_P([test/tabix])])
fi

# @HTSDIRslash_if_relsrcdir@ will be empty when $srcdir is absolute
case "$srcdir" in
  /*) HTSDIRslash_if_relsrcdir= ;;
  *)  HTSDIRslash_if_relsrcdir='$(HTSDIR)/' ;;
esac
AC_SUBST([HTSDIRslash_if_relsrcdir])

AC_OUTPUT
