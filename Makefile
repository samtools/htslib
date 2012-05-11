CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
DFLAGS=
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o
AOBJS=		main.o
PROG=		bcf2ls
INCLUDES=
SUBDIRS=	.
LIBPATH=
LIBCURSES=	

.SUFFIXES:.c .o
.PHONY:all

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bcf2ls:$(LOBJS) $(AOBJS)
		$(CC) $(CFLAGS) $(LOBJS) $(AOBJS) -o $@ -lz

bgzf.o:bgzf.c bgzf.h knetfile.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_USE_KNETFILE $(INCLUDES) bgzf.c -o $@

BCFv2.pdf:BCFv2.tex
		pdflatex BCFv2

kstring.o:kstring.h
knetfile.o:knetfile.h
hts.o:hts.h khash.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
main.o:vcf.h sam.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM $(PROG) *~ *.a BCFv2.aux BCFv2.idx BCFv2.log BCFv2.pdf
