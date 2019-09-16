HTSCODECS_SOURCES = $(HTSPREFIX)htscodecs/htscodecs/arith_dynamic.c \
        $(HTSPREFIX)htscodecs/htscodecs/fqzcomp_qual.c \
        $(HTSPREFIX)htscodecs/htscodecs/pack.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static4x16pr.c \
        $(HTSPREFIX)htscodecs/htscodecs/rANS_static.c \
        $(HTSPREFIX)htscodecs/htscodecs/tokenise_name3.c

HTSCODECS_OBJS = $(HTSCODECS_SOURCES:.c=.o)
