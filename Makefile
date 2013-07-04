CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
DFLAGS=
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o faidx.o razf.o synced_bcf_reader.o vcfutils.o
INCLUDES=
LIBPATH=

.SUFFIXES:.c .o
.PHONY:lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

libhts.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bgzf.o: bgzf.c
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_USE_KNETFILE -DBGZF_MT -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

kstring.o: kstring.c
knetfile.o: knetfile.c
hts.o: hts.c
vcf.o: vcf.c
sam.o: sam.c
tbx.o: tbx.c
faidx.o: faidx.c
razf.o: razf.c
synced_bcf_reader.o: synced_bcf_reader.c
vcfutils.o: vcfutils.c

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
