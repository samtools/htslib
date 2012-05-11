CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
DFLAGS=
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o
INCLUDES=
LIBPATH=

.SUFFIXES:.c .o
.PHONY:lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

libhts.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bgzf.o:bgzf.c bgzf.h knetfile.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_USE_KNETFILE -DBGZF_MT -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

kstring.o:kstring.h
knetfile.o:knetfile.h
hts.o:hts.h khash.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
main.o:vcf.h sam.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
