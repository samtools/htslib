AC_INIT([htslib], [], [jm18@sanger.ac.uk], [], [http://www.htslib.org])

AC_PREREQ(2.68)

AC_PROG_CC

AC_CHECK_LIB([pthread], [pthread_kill])
AC_CHECK_LIB([gssapi_krb5], [gss_acquire_cred])

dnl Save the current values to restore later
CPPFLAGS_RESTORE=$CPPFLAGS
LDFLAGS_RESTORE=$LDFLAGS

dnl Offer support for reading and writing directly to iRODS data objects.
IRODS_HOME=/usr/local/lib/irods
AC_ARG_WITH([irods],
 [AS_HELP_STRING([--with-irods],
   [Enable support for iRODS (default: disabled)])],
 [AS_IF([test "x$with_irods" == xno],  [],
        [test "x$with_irods" == xyes],
        [AC_MSG_NOTICE(Using default iRODS location $IRODS_HOME)],
        [IRODS_HOME="$withval"])],
 [with_irods="no"])

AC_SUBST(IRODS_HOME, [${IRODS_HOME}])

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

AC_CONFIG_SRCDIR([hts.c])
AC_CONFIG_HEADERS([config_auto.h])
AC_CONFIG_FILES([config_auto.mk])
AC_OUTPUT
