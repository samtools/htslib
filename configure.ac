#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([htslib], [1.0], [https://github.com/samtools/htslib/issues], [htslib], [https://github.com/samtools/htslib])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_DEFINE([_USE_KNETFILE], [], [TODO missing description])
AC_DEFINE([BGZF_CACHE], [], [TODO missing description])
AC_DEFINE([BGZF_MT], [], [TODO missing description])

# Checks for programs.
AC_PROG_CC
AM_PROG_AR

# Checks for libtool.
LT_PREREQ([2.2])
LT_INIT([dlopen])

# Checks for libraries.
AC_CHECK_LIB([m], [exp], [], [AC_MSG_ERROR(math not found.)])
AC_CHECK_LIB([pthread], [pthread_create], [], [AC_MSG_ERROR(pthread not found.)])
PKG_CHECK_MODULES(ZLIB, zlib, [], [AC_MSG_ERROR(zlib not found.)])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h inttypes.h limits.h netdb.h pthread.h stddef.h stdint.h stdlib.h string.h sys/socket.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([atexit fdatasync ftruncate getcwd gethostbyaddr gethostbyname getpagesize gettimeofday memmove memset mkdir select socket strcasecmp strchr strdup strerror strrchr strstr strtol strtoul])

# Finish configuring.
AC_CONFIG_FILES([Makefile htslib.pc])
AC_OUTPUT
