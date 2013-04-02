VERSION=
CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2 $(VERSION)
DFLAGS=
OBJS=		main.o samview.o vcfview.o bamidx.o bcfidx.o bamshuf.o bam2fq.o tabix.o \
			abreak.o bam2bed.o vcfcheck.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o \
            vcfnorm.o vcfgtcheck.o
INCLUDES=	-Ihtslib
PROG=		htscmd

.SUFFIXES:.c .o
.PHONY:all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

git-stamp:
		make VERSION="-DVERSION='\"`git describe --always --dirty`\"'"

lib:
		cd htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

htscmd:lib $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) -Lhtslib -lhts -lpthread -lz -lm

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG); cd htslib; $(MAKE) clean; cd ..
