/**********************************************************************
The MIT License

Copyright (c) 2014 Intel Corporation

	Permission is hereby granted, free of charge, to any person
	obtaining a copy of this software and associated documentation
	files (the "Software"), to deal in the Software without
	restriction, including without limitation the rights to use,
	copy, modify, merge, publish, distribute, sublicense, and/or
	sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following
	conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the
	Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
	KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
	WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
	COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
	OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**********************************************************************/
#include "internal_state_size.h"
#include "types.h"

typedef struct {
    UINT8 opaque[INTERNAL_STATE_SIZE];
} LZ_State2;

typedef struct {
    UINT8 *next_in;   // Next input byte
    UINT32 avail_in;  // number of bytes available at next_in
    UINT32 total_in;  // total number of bytes read so far

    UINT8 *next_out;  // Next output byte
    UINT32 avail_out; // number of bytes available at next_out
    UINT32 total_out; // total number of bytes written so far
    UINT32 end_of_stream; // non-zero if this is the last input buffer

    LZ_State2 internal_state;
} LZ_Stream2;


void init_stream(LZ_Stream2 *stream);
void fast_lz(LZ_Stream2 *stream);
