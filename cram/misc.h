#ifndef _misc_h
#define _misc_h

#include "cram/os.h"

#include <stdio.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This informs gcc that crash() doesn't return, so it doesn't need to
 * concern itself that code paths going via crash could mean some variables
 * being undefined and then issuing uninitialised variable warnings.
 * This particularly affected convert.
 */
#ifdef __GNUC__
#    define __NORETURN__ __attribute__ ((__noreturn__))
#else
#    define __NORETURN__
#endif

/*
 * Used for printf style argument checking. We can request a function such
 * as vTcl_SetResult does argument checking, avoiding bugs with using
 * %d and passing in a 64-bit record.
 */
#ifdef __GNUC__
#    define __PRINTF_FORMAT__(a,b) __attribute__ ((format (printf, a, b)))
#else
#    define __PRINTF_FORMAT__(a,b)
#endif

extern int is_directory(char * fn);
extern int is_file(char * fn);
extern int file_size(char * fn);

#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )

#ifdef __cplusplus
}
#endif

#endif /*_misc_h*/
