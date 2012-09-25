CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
DFLAGS=
LOBJS=		knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o synced_bcf_reader.o vcfutils.o
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

knetfile.o:knetfile.h
hts.o:hts.h bgzf.h khash.h kseq.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
tbx.o:tbx.h bgzf.h kstring.h hts.h
synced_bcf_reader.o:vcf.h bgzf.h kstring.h khash.h hts.h tbx.h
vcfutils.o:vcf.h bgzf.h kstring.h khash.h hts.h tbx.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
