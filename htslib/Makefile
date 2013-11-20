CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
DFLAGS=
LOBJS=		knetfile.o kurl.o bgzf.o hts.o vcf.o sam.o tbx.o
INCLUDES=
LIBPATH=

.SUFFIXES:.c .o
.PHONY:lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

libhts.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bgzf.o:bgzf.c bgzf.h knetfile.h khash.h kurl.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_USE_KURL -DBGZF_MT -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

hts.o:hts.c bgzf.h bgzf.h khash.h kurl.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_USE_KURL $(INCLUDES) hts.c -o $@

knetfile.o:knetfile.h
kurl.o:kurl.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
tbx.o:tbx.h bgzf.h kstring.h hts.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
