CC=			gcc
#CFLAGS=		-g -Wall -O2 -Wc++-compat
CFLAGS=		-g -Wall
# SAMTOOLS=1 is used for cram sub-dir.
DFLAGS= -DSAMTOOLS=1
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o faidx.o razf.o synced_bcf_reader.o vcfutils.o $(CRAM_OBJS)

# Cram objects. TODO: merge data structures and migrate rest to this dir?
CRAM_OBJS= \
	cram/cram_codecs.o \
	cram/cram_encode.o \
	cram/cram_io.o \
	cram/cram_decode.o \
	cram/cram_index.o \
	cram/cram_stats.o \
	cram/cram_samtools.o \
	cram/sam_header.o \
	cram/vlen.o \
	cram/zfio.o \
	cram/mFILE.o \
	cram/md5.o \
	cram/open_trace_file.o \
	cram/pooled_alloc.o \
	cram/string_alloc.o \
	cram/files.o \
	cram/thread_pool.o

INCLUDES=-I.
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
hts.o:hts.h bgzf.h khash.h kseq.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
tbx.o:tbx.h bgzf.h kstring.h hts.h
faidx.o:faidx.h khash.h
razdf.o:razf.h knetfile.h zlib.h
synced_bcf_reader.o:vcf.h bgzf.h kstring.h khash.h hts.h tbx.h
vcfutils.o:vcf.h bgzf.h kstring.h khash.h hts.h tbx.h

clean:
	rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
	rm cram/*.o
