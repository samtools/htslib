CC     = gcc
AR     = ar
RANLIB = ranlib

CFLAGS=		-g -Wall -O2 -Wc++-compat
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o faidx.o razf.o synced_bcf_reader.o vcfutils.o


all: lib-static

HTSPREFIX =
include htslib_vars.mk

lib-static: libhts.a


.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@


libhts.a: $(LOBJS)
	@-rm -f $@
	$(AR) -rc $@ $(LOBJS)
	-$(RANLIB) $@

bgzf.o: bgzf.c $(htslib_bgzf_h) htslib/knetfile.h htslib/khash.h
	$(CC) -c $(CFLAGS) -D_USE_KNETFILE -DBGZF_MT -DBGZF_CACHE bgzf.c -o $@

kstring.o: kstring.c htslib/kstring.h
knetfile.o: knetfile.c htslib/knetfile.h
hts.o: hts.c $(htslib_hts_h) $(htslib_bgzf_h) htslib/khash.h htslib/kseq.h htslib/ksort.h
vcf.o: vcf.c $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_tbx_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
sam.o: sam.c $(htslib_sam_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
tbx.o: tbx.c $(htslib_tbx_h) htslib/khash.h
faidx.o: faidx.c $(htslib_faidx_h) htslib/khash.h $(htslib_razf_h) htslib/knetfile.h
razf.o: razf.c $(htslib_razf_h)
synced_bcf_reader.o: synced_bcf_reader.c $(htslib_synced_bcf_reader_h) htslib/kseq.h
vcfutils.o: vcfutils.c $(htslib_vcfutils_h)


mostlyclean:
	-rm -f *.o *.dSYM

clean: mostlyclean
	-rm -f libhts.a

distclean: clean
	-rm -f TAGS


tags:
	ctags -f TAGS *.[ch] htslib/*.h


.PHONY: all clean distclean lib-static mostlyclean tags
