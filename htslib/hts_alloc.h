/*  hts_alloc.h -- Allocation functions.

    Copyright (C) 2026 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_HTS_ALLOC_H
#define HTSLIB_HTS_ALLOC_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Compute saturating addition of two size_t inputs
 *
 *  @return sum of @p a and @p b, or SIZE_MAX if the sum would overflow
 */

static inline size_t hts_add_sat2(size_t a, size_t b) {
    size_t r = a + b;
    return r >= a ? r : SIZE_MAX;
}

/**
 *  Compute saturating addition of three size_t inputs
 *
 *  @return sum of @p a, @p b and @p c, or SIZE_MAX if the sum would overflow
 */

static inline size_t hts_add_sat3(size_t a, size_t b, size_t c) {
    return hts_add_sat2(hts_add_sat2(a, b), c);
}

/**
 *  Compute saturating product of two size_t inputs
 *
 *  @return product of @p a, and @p b, or SIZE_MAX if the sum would overflow
 */

static inline size_t hts_prod_sat2(size_t a, size_t b) {
    const size_t safe_limit = ((size_t) 1) << (sizeof(size_t) * 8 / 2);
    if (a < safe_limit && b < safe_limit)
        return a * b;
    if (a == 0)
        return 0;
    if ((SIZE_MAX / a) < b)
        return SIZE_MAX;
    return a * b;
}

#define HTS_SIZE_LIMIT (PTRDIFF_MAX < SIZE_MAX ? PTRDIFF_MAX : SIZE_MAX)

/**
 *  Allocate memory with checks for overflow
 *
 *  @param size Number of bytes to allocate
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be @p size bytes.  NULL will be returned
 *  if the product is too big for a size_t or ptrdiff_t, or if the memory
 *  could not be allocated.
 *
 *  @note Using this function helps avoid -Walloc-size-larger-than warnings
 *  from gcc, which doesn't like attempts to malloc more than PTRDIFF_MAX.
 *
 */

static inline void * hts_malloc(size_t size) {
    if (size >= HTS_SIZE_LIMIT) {
        errno = ENOMEM;
        return NULL;
    }
    return malloc(size);
}

/**
 *  Allocate memory with checks for overflow, product of sizes
 *
 *  @param a First term of product
 *  @param b Secord term of product
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be @p a * @p b bytes.  NULL will be returned
 *  if the product is too big for a size_t or ptrdiff_t, or if the memory
 *  could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = malloc(number * sizeof(*ptr));
 *  @endcode
 *  as it will detect cases where the product could overflow.
 */

static inline void * hts_malloc_p(size_t a, size_t b) {
    size_t prod = hts_prod_sat2(a, b);
    return hts_malloc(prod);
}

/**
 *  Allocate memory with checks for overflow, product of sum of sizes
 *
 *  @param sz Size of item
 *  @param a  First term of sum
 *  @param b  Secord term of sum
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be @p sz * (@p a + @p b) bytes.  NULL will be
 *  returned if the product is too big for a size_t or ptrdiff_t, or if the
 *  memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = malloc((a + b) * sizeof(*ptr));
 *  @endcode
 *  as it will detect cases where the sum and/or product could overflow.
 */

static inline void * hts_malloc_ps(size_t sz, size_t a, size_t b) {
    size_t sum = hts_add_sat2(a, b);
    size_t prod = hts_prod_sat2(sz, sum);
    return hts_malloc(prod);
}

/**
 *  Allocate memory with checks for overflow, product of sum of sizes plus extra bytes
 *
 *  @param sz     Element size
 *  @param a      First term of sum
 *  @param b      Second term of sum
 *  @param extra  Extra space to allocate, in bytes
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be (@p a + @p b) * @p sz + @p extra bytes.
 *  NULL will be returned if the total is too big for a size_t or
 *  ptrdiff_t, or if the memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = malloc((a + b) * sizeof(*ptr) + extra);
 *  @endcode
 *  as it will detect cases where the sums and/or product could overflow.
 */

static inline void * hts_malloc_pse(size_t sz, size_t a, size_t b, size_t extra) {
    size_t sum = hts_add_sat2(a, b);
    size_t prod = hts_prod_sat2(sz, sum);
    size_t bytes = hts_add_sat2(prod, extra);
    return hts_malloc(bytes);
}

/**
 *  Allocate and clear memory with checks for overflow
 *
 *  @param sz  Size in bytes of each element
 *  @param num Number of elements
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be @p size bytes.  NULL will be returned
 *  if the product is too big for a size_t or ptrdiff_t, or if the memory
 *  could not be allocated.
 *
 *  @note Using this function helps avoid -Walloc-size-larger-than warnings
 *  from gcc, which doesn't like attempts to malloc more than PTRDIFF_MAX.
 *
 *  @note In this function the size parameter is listed first to match
 *  the other hts_alloc interfaces, which is the opposite way round to calloc().
 *  In practice this makes no difference as the values get multiplied together.
 */

static inline void * hts_calloc(size_t sz, size_t num) {
    if (num >= HTS_SIZE_LIMIT || sz >= HTS_SIZE_LIMIT) {
        errno = ENOMEM;
        return NULL;
    }
    return calloc(num, sz);
}

/**
 *  Allocate and clear memory with checks for overflow, product of sum of sizes
 *
 *  @param sz     Element size
 *  @param a      First term of sum
 *  @param b      Second term of sum
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be (@p a + @p b) * @p sz + @p extra bytes.
 *  NULL will be returned if the total is too big for a size_t or
 *  ptrdiff_t, or if the memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = calloc(a + b, sizeof(*ptr));
 *  @endcode
 *  as it will detect cases where the sum could overflow.
 */

static inline void * hts_calloc_ps(size_t sz, size_t a, size_t b) {
    size_t sum = hts_add_sat2(a, b);
    return hts_calloc(sz, sum);
}

/**
 *  Allocate and clear memory with checks for overflow
 *
 *  @param sz     Element size
 *  @param a      First term of sum
 *  @param b      Second term of sum
 *  @param extra  Extra space to allocate, in bytes
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be (@p a + @p b) * @p sz + @p extra bytes.
 *  NULL will be returned if the total is too big for a size_t or
 *  ptrdiff_t, or if the memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = calloc((a + b) * sizeof(*ptr) + extra, 1);
 *  @endcode
 *  as it will detect cases where the sums and/or product could overflow.
 */

static inline void * hts_calloc_pse(size_t sz, size_t a, size_t b, size_t extra) {
    size_t sum = hts_add_sat2(a, b);
    size_t prod = hts_prod_sat2(sz, sum);
    size_t bytes = hts_add_sat2(prod, extra);
    return hts_calloc(1, bytes);
}

/**
 *  Rellocate memory with checks for overflow
 *
 *  @return Pointer to the allocated memory, or NULL on error
 */

static inline void * hts_realloc(void *orig, size_t size) {
    if (size >= HTS_SIZE_LIMIT) {
        errno = ENOMEM;
        return NULL;
    }
    return realloc(orig, size);
}

/**
 *  Rellocate memory with checks for overflow
 *
 *  @param sz Size of item
 *  @param a  First term of sum
 *  @param b  Secord term of sum
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be @p sz * (@p a + @p b) bytes.  NULL will be
 *  returned if the product is too big for a size_t or ptrdiff_t, or if the
 *  memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = realloc(orig, (a + b) * sizeof(*ptr));
 *  @endcode
 *  as it will detect cases where the sum and/or product could overflow.
 */

static inline void * hts_realloc_p(void *orig, size_t a, size_t b) {
    size_t prod = hts_prod_sat2(a, b);
    return hts_realloc(orig, prod);
}

/**
 *  Rellocate memory with checks for overflow
 *
 *  @param sz     Element size
 *  @param a      First term of sum
 *  @param b      Second term of sum
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be (@p a + @p b) * @p sz bytes.
 *  NULL will be returned if the total is too big for a size_t or
 *  ptrdiff_t, or if the memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = realloc(orig, (a + b) * sizeof(*ptr));
 *  @endcode
 *  as it will detect cases where the sum and/or product could overflow.
 */

static inline void * hts_realloc_ps(void *orig, size_t sz, size_t a, size_t b) {
    size_t sum = hts_add_sat2(a, b);
    size_t prod = hts_prod_sat2(sz, sum);
    return hts_realloc(orig, prod);
}

/**
 *  Rellocate memory with checks for overflow
 *
 *  @param sz     Element size
 *  @param a      First term of sum
 *  @param b      Second term of sum
 *  @param extra  Extra space to allocate, in bytes
 *
 *  @return Pointer to the allocated memory, or NULL on error
 *
 *  The total allocated will be (@p a + @p b) * @p sz + @p extra bytes.
 *  NULL will be returned if the total is too big for a size_t or
 *  ptrdiff_t, or if the memory could not be allocated.
 *
 *  This function is a safer replacement for the idiom
 *  @code{c}
 *   ptr = realloc(orig, (a + b) * sizeof(*ptr) + extra);
 *  @endcode
 *  as it will detect cases where the sums and/or product could overflow.
 */

static inline void * hts_realloc_pse(void *orig, size_t sz, size_t a, size_t b,
                                     size_t extra) {
    size_t sum = hts_add_sat2(a, b);
    size_t prod = hts_prod_sat2(sz, sum);
    size_t bytes = hts_add_sat2(prod, extra);
    return hts_realloc(orig, bytes);
}

#undef HTS_SIZE_LIMIT

#ifdef __cplusplus
}
#endif

#endif /* HTSLIB_HTS_ALLOC_H */
