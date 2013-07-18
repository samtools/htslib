/*
 * File: os.h
 *
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */

#ifndef _OS_H_
#define _OS_H_

#include <limits.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 * Detection of endianness. The main part of this is done in autoconf, but
 * for the case of MacOS FAT binaries we fall back on auto-sensing based on
 * processor type too.
 */

/* Set by autoconf */
#define SP_LITTLE_ENDIAN

/* Mac FAT binaries or unknown. Auto detect based on CPU type */
#if !defined(SP_BIG_ENDIAN) && !defined(SP_LITTLE_ENDIAN)

/*
 * x86 equivalents
 */
#if defined(__i386__) || defined(__i386) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686)
#  if defined(SP_BIG_ENDIAN)
#    undef SP_BIG_ENDIAN
#  endif
#  define SP_LITTLE_ENDIAN
#endif

/*
 * DEC Alpha
 */
#if defined(__alpha__) || defined(__alpha)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/*
 * SUN Sparc
 */
#if defined(__sparc__) || defined(__sparc)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/*
 * PowerPC
 */
#if defined(__ppc__) || defined(__ppc)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/* Some catch-alls */
#if defined(__LITTLE_ENDIAN__) || defined(__LITTLEENDIAN__)
#    define SP_LITTLE_ENDIAN
#endif

#if defined(__BIG_ENDIAN__) || defined(__BIGENDIAN__)
#    define SP_BIG_ENDIAN
#endif

#if defined(SP_BIG_ENDIAN) && defined(SP_LITTLE_ENDIAN)
#    error Both BIG and LITTLE endian defined. Fix os.h and/or Makefile
#endif

#if !defined(SP_BIG_ENDIAN) && !defined(SP_LITTLE_ENDIAN)
#    error Neither BIG nor LITTLE endian defined. Fix os.h and/or Makefile
#endif

#endif

/*-----------------------------------------------------------------------------
 * Allow for unaligned memory access. This is used in BAM code as the packed
 * structure has 4-byte cigar ints after the variable length name.
 *
 * Consider using AX_CHECK_ALIGNED_ACCESS_REQUIRED in autoconf.
 */
#if defined(__i386__) || defined(__i386) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686)
#  define ALLOW_UAC
#endif

/*-----------------------------------------------------------------------------
 * Byte swapping macros
 */

/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define iswap_int8(x) \
    (((x & 0x00000000000000ffLL) << 56) + \
     ((x & 0x000000000000ff00LL) << 40) + \
     ((x & 0x0000000000ff0000LL) << 24) + \
     ((x & 0x00000000ff000000LL) <<  8) + \
     ((x & 0x000000ff00000000LL) >>  8) + \
     ((x & 0x0000ff0000000000LL) >> 24) + \
     ((x & 0x00ff000000000000LL) >> 40) + \
     ((x & 0xff00000000000000LL) >> 56))

#define iswap_int4(x) \
    (((x & 0x000000ff) << 24) + \
     ((x & 0x0000ff00) <<  8) + \
     ((x & 0x00ff0000) >>  8) + \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x) \
    (((x & 0x00ff) << 8) + \
     ((x & 0xff00) >> 8))

/*
 * Linux systems may use byteswap.h to get assembly versions of byte-swap
 * on intel systems. This can be as trivial as the bswap opcode, which works
 * out at over 2-times faster than iswap_int4 above.
 */
#if 0
#if defined(__linux__)
#    include <byteswap.h>
#    undef iswap_int8
#    undef iswap_int4
#    undef iswap_int2
#    define iswap_int8 bswap_64
#    define iswap_int4 bswap_32
#    define iswap_int2 bswap_16
#endif
#endif


/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#ifdef SP_BIG_ENDIAN
#define le_int4(x) iswap_int4((x))
#endif

#ifdef SP_LITTLE_ENDIAN
#define le_int4(x) (x)
#endif

/*-----------------------------------------------------------------------------
 * <inttypes.h> definitions, incase they're not present
 */

#ifndef PRId64
#define __PRI64__ "l"
#define PRId64 __PRI64__ "d"
#define PRId32 "d"
#define PRId16 "d"
#define PRId8  "d"
#define PRIu64 __PRI64__ "u"
#define PRIu32 "u"
#define PRIu16 "u"
#define PRIu8  "u"
#endif

/*-----------------------------------------------------------------------------
 * Operating system specifics.
 * These ought to be done by autoconf, but are legacy code.
 */
/*
 * SunOS 4.x
 * Even though we use the ANSI gcc, we make use the the standard SunOS 4.x
 * libraries and include files, which are non-ansi
 */
#if defined(__sun__) && !defined(__svr4__)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

/*
 * Microsoft Visual C++
 * Windows
 */
#if defined(_MSC_VER)
#define popen _popen
#define pclose _pclose
#define ftruncate(fd,len) _chsize(fd,len)
#endif


/*
 * Microsoft Windows running MinGW
 */
#if defined(__MINGW32__)
/* #define mkdir(filename,mode) mkdir((filename)) */
#define sysconf(x) 512
#define ftruncate(fd,len) _chsize(fd,len)
#endif

/* Generic WIN32 API issues */
#ifdef _WIN32
#  ifndef HAVE_FSEEKO
#    if __MSVCRT_VERSION__ >= 0x800
       /* if you have MSVCR80 installed then you can use these definitions: */
#      define off_t __int64
#      define fseeko _fseeki64
#      define ftello _ftelli64
#    else
       /* otherwise we're stuck with 32-bit file support */
#      define off_t long
#      define fseeko fseek
#      define ftello ftell
#    endif
#  endif /* !HAVE_FSEEKO */
#endif /* _WIN32 */

#ifdef __cplusplus
}
#endif

#endif /*_OS_H_*/
