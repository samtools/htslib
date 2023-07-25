# Makefile fragment to add settings needed when bundling htscodecs functions
#
#    Copyright (C) 2021-2022 Genome Research Ltd.
#
#    Author: Rob Davies <rmd@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


HTSCODECS_SOURCES = $(HTSPREFIX)htscodecs/htscodecs/arith_dynamic.c \
        $(HTSPREFIX)htscodecs/htscodecs/fqzcomp_qual.c \
        $(HTSPREFIX)htscodecs/htscodecs/htscodecs.c \
        $(HTSPREFIX)htscodecs/htscodecs/pack.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static4x16pr.c \
	$(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_avx2.c \
	$(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_avx512.c \
	$(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_sse4.c \
	$(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr_neon.c \
	$(HTSPREFIX)htscodecs/htscodecs/rANS_static32x16pr.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static.c \
        $(HTSPREFIX)htscodecs/htscodecs/rle.c \
        $(HTSPREFIX)htscodecs/htscodecs/tokenise_name3.c \
	$(HTSPREFIX)htscodecs/htscodecs/utils.c


HTSCODECS_OBJS = $(HTSCODECS_SOURCES:.c=.o)

# htscodecs public headers
htscodecs_arith_dynamic_h = htscodecs/htscodecs/arith_dynamic.h
htscodecs_fqzcomp_qual_h = htscodecs/htscodecs/fqzcomp_qual.h
htscodecs_htscodecs_h = htscodecs/htscodecs/htscodecs.h $(htscodecs_version_h)
htscodecs_pack_h = htscodecs/htscodecs/pack.h
htscodecs_rANS_static_h = htscodecs/htscodecs/rANS_static.h
htscodecs_rANS_static4x16_h = htscodecs/htscodecs/rANS_static4x16.h
htscodecs_rle_h = htscodecs/htscodecs/rle.h
htscodecs_tokenise_name3_h = htscodecs/htscodecs/tokenise_name3.h
htscodecs_varint_h = htscodecs/htscodecs/varint.h

# htscodecs internal headers
htscodecs_htscodecs_endian_h = htscodecs/htscodecs/htscodecs_endian.h
htscodecs_c_range_coder_h = htscodecs/htscodecs/c_range_coder.h
htscodecs_c_simple_model_h = htscodecs/htscodecs/c_simple_model.h $(htscodecs_c_range_coder_h)
htscodecs_permute_h = htscodecs/htscodecs/permute.h
htscodecs_pooled_alloc_h = htscodecs/htscodecs/pooled_alloc.h
htscodecs_rANS_byte_h = htscodecs/htscodecs/rANS_byte.h
htscodecs_rANS_static16_int_h = htscodecs/htscodecs/rANS_static16_int.h $(htscodecs_varint_h) $(htscodecs_utils_h)
htscodecs_rANS_static32x16pr_h = htscodecs/htscodecs/rANS_static32x16pr.h
htscodecs_rANS_word_h = htscodecs/htscodecs/rANS_word.h $(htscodecs_htscodecs_endian_h)
htscodecs_utils_h = htscodecs/htscodecs/utils.h
htscodecs_version_h = htscodecs/htscodecs/version.h

# Add htscodecs tests into the HTSlib test framework

HTSCODECS_TEST_TARGETS = test_htscodecs_rans4x8 \
    test_htscodecs_rans4x16 test_htscodecs_arith test_htscodecs_tok3 \
    test_htscodecs_fqzcomp test_htscodecs_varint
