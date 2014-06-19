HTSLIB=../
INCFLAGS=-I$(HTSLIB) -g -g3
all : memcheck

memcheck:memcheck.c
	gcc -o memcheck $(INCFLAGS) memcheck.c -Wl,-Bstatic -L$(HTSLIB) -lhts -Wl,-Bdynamic -lz -lpthread

clean:
	rm -f memcheck
