CC=         gcc
CXX=        g++
CFLAGS=     -g -Wall -O2 #-m64 #-arch ppc
CXXFLAGS=   $(CFLAGS)
DFLAGS=     -D_FILE_OFFSET_BITS=64
OBJS=       
PROG=       test-vcf-api
INCLUDES=   -I../htslib
SUBDIRS=    .

HTSDIR = ../
HTSLIB = $(HTSDIR)/libhts.a


.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES:-I.%=-I../.%) $< -o $@

all:$(PROG)

test:$(PROG) run-tests

run-tests:
	./test.pl

test-vcf-api:test-vcf-api.o $(HTSLIB)
	$(CC) $(CFLAGS) -o $@ test-vcf-api.o -lm $(HTSLIB) -lpthread -lz

clean:
	rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a

