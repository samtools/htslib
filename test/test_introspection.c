/*  test/test_introspection.c -- demonstration of introspection function usage

    Copyright (C) 2020-2021 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

#include <config.h>
#include <stdio.h>

#include "../htslib/hts.h"
#include "../htslib/hfile.h"

int main(void) {
    printf("Version string: %s\n", hts_version());
    printf("Version number: %d\n", HTS_VERSION);
    printf("\nhtscodecs version: %s\n",
           hts_test_feature(HTS_FEATURE_HTSCODECS));

    printf("\nCC:             %s\n", hts_test_feature(HTS_FEATURE_CC));
    printf("CPPFLAGS:       %s\n", hts_test_feature(HTS_FEATURE_CPPFLAGS));
    printf("CFLAGS:         %s\n", hts_test_feature(HTS_FEATURE_CFLAGS));
    printf("LDFLAGS:        %s\n", hts_test_feature(HTS_FEATURE_LDFLAGS));

    unsigned int feat = hts_features();
    printf("\nFeature number: 0x%x\n", feat);
    if (feat & HTS_FEATURE_CONFIGURE)
        printf("                HTS_FEATURE_CONFIGURE\n");
    if (feat & HTS_FEATURE_PLUGINS)
        printf("                HTS_FEATURE_PLUGINS\n");
    if (feat & HTS_FEATURE_LIBCURL)
        printf("                HTS_FEATURE_LIBCURL\n");
    if (feat & HTS_FEATURE_S3)
        printf("                HTS_FEATURE_S3\n");
    if (feat & HTS_FEATURE_GCS)
        printf("                HTS_FEATURE_GCS\n");
    if (feat & HTS_FEATURE_LIBDEFLATE)
        printf("                HTS_FEATURE_LIBDEFLATE\n");
    if (feat & HTS_FEATURE_LZMA)
        printf("                HTS_FEATURE_LZMA\n");
    if (feat & HTS_FEATURE_BZIP2)
        printf("                HTS_FEATURE_BZIP2\n");
    if (feat & HTS_FEATURE_HTSCODECS)
        printf("                HTS_FEATURE_HTSCODECS\n");

    printf("\nFeature string: %s\n", hts_feature_string());


    // Plugins and schemes
    printf("\nPlugins present:\n");
    const char *plugins[100];
    int np = 100, i, j;

    if (hfile_list_plugins(plugins, &np) < 0)
        return 1;

    for (i = 0; i < np; i++) {
        const char *sc_list[100];
        int nschemes = 100;
        if (hfile_list_schemes(plugins[i], sc_list, &nschemes) < 0)
            return 1;

        printf("    %s:\n", plugins[i]);
        for (j = 0; j < nschemes; j++)
            printf("\t%s\n", sc_list[j]);
        puts("");
    }

    return 0;
}
