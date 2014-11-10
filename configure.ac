AC_INIT([htslib], [1.1], [jm18@sanger.ac.uk], [], [http://www.htslib.org])

AC_PREREQ(2.68)

AC_PROG_CC

AC_CHECK_LIB([pthread], [pthread_kill])
AC_CHECK_LIB([gssapi_krb5], [gss_acquire_cred])

dnl Moved here from the hand-rolled Makefile
AC_SUBST([LIBHTS_SOVERSION], [1])

dnl Save the current values to restore later
CPPFLAGS_RESTORE=$CPPFLAGS
LDFLAGS_RESTORE=$LDFLAGS

dnl Offer support for reading and writing directly to iRODS data objects.
with_irods=no
AC_ARG_WITH([irods],
            [AS_HELP_STRING([--with-irods],
              [Enable support for iRODS (default: disabled)])],
            [],
            [with_irods=no])

IRODS_HOME=
AS_IF([test "x$with_irods" != xno && test "x$with_irods" != xyes],
  [IRODS_HOME="$with_irods"
   AC_SUBST(IRODS_HOME, [$with_irods])])

dnl iRODS leaves its headers and shared libraries in non-standard
dnl places. However, these places may be determined relative to
dnl IRODS_HOME.
IRODS_CPPFLAGS=\
"-I${IRODS_HOME}/lib/api/include \
-I${IRODS_HOME}/lib/core/include \
-I${IRODS_HOME}/lib/md5/include \
-I${IRODS_HOME}/lib/sha1/include \
-I${IRODS_HOME}/server/core/include \
-I${IRODS_HOME}/server/drivers/include \
-I${IRODS_HOME}/server/icat/include \
-I${IRODS_HOME}/server/re/include"
IRODS_LDFLAGS="-L${IRODS_HOME}/lib/core/obj"

dnl Set these temporarily for AC_CHECK_LIB
CPPFLAGS="$CPPFLAGS ${IRODS_CPPFLAGS}"
LDFLAGS="$LDFLAGS ${IRODS_LDFLAGS}"

HAVE_LIBRODSAPIS="no"
AS_IF([test "x$with_irods" != xno],
  [AC_CHECK_LIB([RodsAPIs], [getRodsEnvFileName],
    [AC_SUBST([HAVE_LIBRODSAPIS],  ["yes"])
     AC_SUBST([IRODS_CPPFLAGS],    ["${IRODS_CPPFLAGS}"])
     AC_SUBST([IRODS_LDFLAGS],     ["${IRODS_LDFLAGS}"])
     AC_SUBST([IRODS_LDLIBS],      ["-lRodsAPIs -lgssapi_krb5"])
     AC_DEFINE([HAVE_LIBRODSAPIS], [1], [libRodsAPIs is present])],
    [AC_MSG_FAILURE([--with-irods was given, but test for libRodsAPIs failed])],
    [-lpthread -lgssapi_krb5])])

dnl Restore the original values
CPPFLAGS=$CPPFLAGS_RESTORE
LDFLAGS=$LDFLAGS_RESTORE

dnl Ensure that these defines from the original hand-rolled config.h
dnl are maintained.
AC_DEFINE([_USE_KNETFILE], [])
AC_DEFINE([BGZF_CACHE],    [])
AC_DEFINE([BGZF_MT],       [])

AH_TEMPLATE([_USE_KNETFILE],    [Use KNETFILE])
AH_TEMPLATE([BGZF_CACHE],       [Use BGZF cache])
AH_TEMPLATE([BGZF_MT],          [Use BGZF MT])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([hts.c])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
