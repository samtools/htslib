/* The MIT License

   Copyright (c) 2016 Illumina Cambridge Ltd.

   Author: Peter Krusche <pkrusche@illumina.com>

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef HTSLIB_HFILE_MEM_H
#define HTSLIB_HFILE_MEM_H

#include "hfile.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Buffer lookup callback. Given a file name, returns a buffer and size.
 *
 * When hopen_mem is called to read a file with a name that starts with '@',
 * it will use such a function to obtain a buffer pointer. This allows us to
 * feed arbitrary memory blocks into htslib for decompression / parsing.
 *
 * @param name the file / internal handle name.
 * @param buffer void pointer that will receive the buffer
 * @param length size_t pointer that will receive the length of the data pointed to in buffer
 */
typedef int (*buffer_lookup_fn)(const char * name, void** buffer, size_t * length);

/**
 * Set buffer lookup function for memory files.
 * @param fn function of type buffer_lookup_fn
 */
extern void hfile_mem_set_lookup_function(buffer_lookup_fn fn);

/**
 * Get buffer for a hfile
 * @param file the file to use. This should be a hFILE that was opened using hfile_mem
 * @param buffer void pointer that will receive the buffer
 * @param length size_t pointer that will receive the length of the data pointed to in buffer
 *
 * @return 0 if successful an error code otherwise
 */
extern int hfile_mem_get_buffer(hFILE * file, void ** buffer, size_t * length);

#ifdef __cplusplus
};
#endif
#endif //HTSLIB_HFILE_MEM_H
