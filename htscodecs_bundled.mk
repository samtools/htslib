HTSCODECS_SOURCES = $(HTSPREFIX)htscodecs/htscodecs/arith_dynamic.c \
        $(HTSPREFIX)htscodecs/htscodecs/fqzcomp_qual.c \
        $(HTSPREFIX)htscodecs/htscodecs/pack.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static4x16pr.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static.c \
        $(HTSPREFIX)htscodecs/htscodecs/rle.c \
        $(HTSPREFIX)htscodecs/htscodecs/tokenise_name3.c

HTSCODECS_OBJS = $(HTSCODECS_SOURCES:.c=.o)

htscodecs_arith_dynamic_h = htscodecs/htscodecs/arith_dynamic.h
htscodecs_fqzcomp_qual_h = htscodecs/htscodecs/fqzcomp_qual.h
htscodecs_pack_h = htscodecs/htscodecs/pack.h
htscodecs_rANS_static_h = htscodecs/htscodecs/rANS_static.h
htscodecs_rANS_static4x16_h = htscodecs/htscodecs/rANS_static4x16.h
htscodecs_rle_h = htscodecs/htscodecs/rle.h
htscodecs_tokenise_name3_h = htscodecs/htscodecs/tokenise_name3.h
htscodecs_varint_h = htscodecs/htscodecs/varint.h

htscodecs_rANS_byte_h = htscodecs/htscodecs/rANS_byte.h
htscodecs_rANS_word_h = htscodecs/htscodecs/rANS_word.h
htscodecs_c_range_coder_h = htscodecs/htscodecs/c_range_coder.h
htscodecs_c_simple_model_h = htscodecs/htscodecs/c_simple_model.h $(htscodecs_c_range_coder_h)
htscodecs_pooled_alloc_h = htscodecs/htscodecs/pooled_alloc.h

# Add htscodecs tests into the HTSlib test framework

HTSCODECS_TEST_TARGETS = test_htscodecs_rans4x8 \
    test_htscodecs_rans4x16 test_htscodecs_arith test_htscodecs_tok3 \
    test_htscodecs_fqzcomp test_htscodecs_varint
