/* The MIT License

   Copyright (C) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef kroundup_size_t
#define kroundup_size_t(x) (--(x),                                       \
                            (x)|=(x)>>(sizeof(size_t)/8), /*  0 or  1 */ \
                            (x)|=(x)>>(sizeof(size_t)/4), /*  1 or  2 */ \
                            (x)|=(x)>>(sizeof(size_t)/2), /*  2 or  4 */ \
                            (x)|=(x)>>(sizeof(size_t)),   /*  4 or  8 */ \
                            (x)|=(x)>>(sizeof(size_t)*2), /*  8 or 16 */ \
                            (x)|=(x)>>(sizeof(size_t)*4), /* 16 or 32 */ \
                            ++(x))
#endif

#if defined __GNUC__ && (__GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4))
#define KS_ATTR_PRINTF(fmt, arg) __attribute__((__format__ (__printf__, fmt, arg)))
#else
#define KS_ATTR_PRINTF(fmt, arg)
#endif

#ifndef HAVE___BUILTIN_CLZ
#if defined __GNUC__ && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#define HAVE___BUILTIN_CLZ 1
#endif
#endif

/* kstring_t is a simple non-opaque type whose fields are likely to be
 * used directly by user code (but see also ks_str() and ks_len() below).
 * A kstring_t object is initialised by either of
 *       kstring_t str = { 0, 0, NULL };
 *       kstring_t str; ...; str.l = str.m = 0; str.s = NULL;
 * and either ownership of the underlying buffer should be given away before
 * the object disappears (see ks_release() below) or the kstring_t should be
 * destroyed with  free(str.s);  */
#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct {
	uint64_t tab[4];
	int sep, finished;
	const char *p; // end of the current token
} ks_tokaux_t;

#ifdef __cplusplus
extern "C" {
#endif

	int kvsprintf(kstring_t *s, const char *fmt, va_list ap) KS_ATTR_PRINTF(2,0);
	int ksprintf(kstring_t *s, const char *fmt, ...) KS_ATTR_PRINTF(2,3);
    int kputd(double d, kstring_t *s); // custom %g only handler
	int ksplit_core(char *s, int delimiter, int *_max, int **_offsets);
	char *kstrstr(const char *str, const char *pat, int **_prep);
	char *kstrnstr(const char *str, const char *pat, int n, int **_prep);
	void *kmemmem(const void *_str, int n, const void *_pat, int m, int **_prep);

	/* kstrtok() is similar to strtok_r() except that str is not
	 * modified and both str and sep can be NULL. For efficiency, it is
	 * actually recommended to set both to NULL in the subsequent calls
	 * if sep is not changed. */
	char *kstrtok(const char *str, const char *sep, ks_tokaux_t *aux);

	/* kgetline() uses the supplied fgets()-like function to read a "\n"-
	 * or "\r\n"-terminated line from fp.  The line read is appended to the
	 * kstring without its terminator and 0 is returned; EOF is returned at
	 * EOF or on error (determined by querying fp, as per fgets()). */
	typedef char *kgets_func(char *, int, void *);
	int kgetline(kstring_t *s, kgets_func *fgets, void *fp);
	typedef ssize_t kgets_func2(char *, int, void *);
	int kgetline2(kstring_t *s, kgets_func2 *fgets, void *fp);

#ifdef __cplusplus
}
#endif

static inline int ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		char *tmp;
		kroundup_size_t(size);
		tmp = (char*)realloc(s->s, size);
		if (!tmp)
			return -1;
		s->s = tmp;
		s->m = size;
	}
	return 0;
}

static inline char *ks_str(kstring_t *s)
{
	return s->s;
}

static inline size_t ks_len(kstring_t *s)
{
	return s->l;
}

// Give ownership of the underlying buffer away to something else (making
// that something else responsible for freeing it), leaving the kstring_t
// empty and ready to be used again, or ready to go out of scope without
// needing  free(str.s)  to prevent a memory leak.
static inline char *ks_release(kstring_t *s)
{
	char *ss = s->s;
	s->l = s->m = 0;
	s->s = NULL;
	return ss;
}

static inline int kputsn(const char *p, size_t l, kstring_t *s)
{
	size_t new_sz = s->l + l + 2;
	if (new_sz < s->l || ks_resize(s, new_sz) < 0)
		return EOF;
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (ks_resize(s, s->l + 2) < 0)
		return EOF;
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

static inline int kputc_(int c, kstring_t *s)
{
	if (ks_resize(s, s->l + 1) < 0)
		return EOF;
	s->s[s->l++] = c;
	return 1;
}

static inline int kputsn_(const void *p, size_t l, kstring_t *s)
{
	size_t new_sz = s->l + l;
	if (new_sz < s->l || ks_resize(s, new_sz) < 0)
		return EOF;
	memcpy(s->s + s->l, p, l);
	s->l += l;
	return l;
}

static inline int kputuw(unsigned x, kstring_t *s)
{
#if HAVE___BUILTIN_CLZ && UINT_MAX == 4294967295U
    static const unsigned int kputuw_num_digits[32] = {
        10, 10, 10,  9,  9,  9,  8,  8,
        8,   7,  7,  7,  7,  6,  6,  6,
        5,   5,  5,  4,  4,  4,  4,  3,
        3,   3,  2,  2,  2,  1,  1,  1
    };
    static const unsigned int kputuw_thresholds[32] = {
        0,        0, 1000000000U, 0,       0, 100000000U,   0,      0,
        10000000, 0,          0,  0, 1000000,         0,    0, 100000,
        0,        0,      10000,  0,       0,         0, 1000,      0,
        0,      100,          0,  0,      10,         0,    0,      0
    };
#else
    uint64_t m;
#endif
    static const char kputuw_dig2r[] =
        "00010203040506070809"
        "10111213141516171819"
        "20212223242526272829"
        "30313233343536373839"
        "40414243444546474849"
        "50515253545556575859"
        "60616263646566676869"
        "70717273747576777879"
        "80818283848586878889"
        "90919293949596979899";
    unsigned int l, j;
    char *cp;

    // Trivial case - also prevents __builtin_clz(0), which is undefined
    if (x < 10) {
        if (ks_resize(s, s->l + 2) < 0)
            return EOF;
        s->s[s->l++] = '0'+x;
        s->s[s->l] = 0;
        return 0;
    }

    // Find out how many digits are to be printed.
#if HAVE___BUILTIN_CLZ && UINT_MAX == 4294967295U
    /*
     * Table method - should be quick if clz can be done in hardware.
     * Find the most significant bit of the value to print and look
     * up in a table to find out how many decimal digits are needed.
     * This number needs to be adjusted by 1 for cases where the decimal
     * length could vary for a given number of bits (for example,
     * a four bit number could be between 8 and 15).
     */

    l = __builtin_clz(x);
    l = kputuw_num_digits[l] - (x < kputuw_thresholds[l]);
#else
    // Fallback for when clz is not available
    m = 1;
    l = 0;
    do {
        l++;
        m *= 10;
    } while (x >= m);
#endif

    if (ks_resize(s, s->l + l + 2) < 0)
        return EOF;

    // Add digits two at a time
    j = l;
    cp = s->s + s->l;
    while (x >= 10) {
        const char *d = &kputuw_dig2r[2*(x%100)];
        x /= 100;
        memcpy(&cp[j-=2], d, 2);
    }

    // Last one (if necessary).  We know that x < 10 by now.
    if (j == 1)
        cp[0] = x + '0';

    s->l += l;
    s->s[s->l] = 0;
    return 0;
}

static inline int kputw(int c, kstring_t *s)
{
    unsigned int x = c;
    if (c < 0) {
        x = -x;
        if (ks_resize(s, s->l + 3) < 0)
            return EOF;
        s->s[s->l++] = '-';
    }

    return kputuw(x, s);
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	int i, l = 0;
	unsigned long x = c;
	if (c < 0) x = -x;
	do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
	if (c < 0) buf[l++] = '-';
	if (ks_resize(s, s->l + l + 2) < 0)
		return EOF;
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

/*
 * Returns 's' split by delimiter, with *n being the number of components;
 *         NULL on failue.
 */
static inline int *ksplit(kstring_t *s, int delimiter, int *n)
{
	int max = 0, *offsets = 0;
	*n = ksplit_core(s->s, delimiter, &max, &offsets);
	return offsets;
}

#endif
