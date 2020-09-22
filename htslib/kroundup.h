/* The MIT License

   Copyright (C) 2020 Genome Research Ltd.

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

#ifndef KROUNDUP_H
#define KROUNDUP_H

// Value of this macro is 1 if x is a signed type; 0 if unsigned
#define k_signed_type(x) (!(-((x) * 0 + 1) > 0))

/*
  Macro with value 1 if the highest bit in x is set for any integer type

  This is written avoiding conditionals (?: operator) to reduce the likelihood
  of gcc attempting jump thread optimisations for code paths where (x) is
  large.  These optimisations can cause gcc to issue warnings about excessively
  large memory allocations when the kroundup64() macro below is used with
  malloc().  Such warnings can be misleading as they imply only the large
  allocation happens when it's actually working fine for normal values of (x).

  See https://developers.redhat.com/blog/2019/03/13/understanding-gcc-warnings-part-2/
*/
#define k_high_bit_set(x) ((((x) >> (sizeof(x) * 8 - 1 - k_signed_type(x))) & 1))

/*! @hideinitializer
  @abstract  Round up to next power of two
  @discussion
  This macro will work for unsigned types up to uint64_t.

  If the next power of two does not fit in the given type, it will set
  the largest value that does.
 */
#define kroundup64(x) ((x) > 0 ?                                        \
                       (--(x),                                          \
                        (x)|=(x)>>(sizeof(x)/8),                        \
                        (x)|=(x)>>(sizeof(x)/4),                        \
                        (x)|=(x)>>(sizeof(x)/2),                        \
                        (x)|=(x)>>(sizeof(x)),                          \
                        (x)|=(x)>>(sizeof(x)*2),                        \
                        (x)|=(x)>>(sizeof(x)*4),                        \
                        (x) += !k_high_bit_set(x),                      \
                        (x))                                            \
                       : 0)

// Historic interfaces for 32-bit and size_t values.  The macro above
// works for both (as long as size_t is no more than 64 bits).

#ifndef kroundup32
#define kroundup32(x) kroundup64(x)
#endif
#ifndef kroundup_size_t
#define kroundup_size_t(x) kroundup64(x)
#endif

#endif
