/* test/test-ref.c -- ref unit tests

   Copyright (C) 2017 Genome Research Ltd

   Author: Thomas Hickman <th10@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE.
*/

#include "htslib/ref.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
    const char* m5_str = "bbf4de6d8497a119dda6e074521643dc";

    int error_code = EXIT_SUCCESS;

    hFILE* ref;

    if (!(ref = m5_to_ref(m5_str))){
        fprintf(stderr, "Error in m5_to_ref\n");
        return EXIT_FAILURE;
    }

    char buf[100];

    size_t size_read = hread(ref, buf, 100);
    if(size_read <= 0){
        fprintf(stderr, "Invalid hfile size read\n");
        return EXIT_FAILURE;
    }

    if(hclose(ref) != 0){
        fprintf(stderr, "Cannot close hfile\n");
        return EXIT_FAILURE;
    }

    return error_code;
}
