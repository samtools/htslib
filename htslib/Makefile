CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
# SAMTOOLS=1 is used for io_lib sub-dir.
DFLAGS= -DSAMTOOLS=1
LOBJS=		kstring.o knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o faidx.o razf.o synced_bcf_reader.o vcfutils.o $(CRAM_OBJS)

# Cram objects. TODO: merge data structures and migrate rest to this dir?
CRAM_OBJS= \
	io_lib/cram_codecs.o \
	io_lib/cram_encode.o \
	io_lib/cram_io.o \
	io_lib/cram_decode.o \
	io_lib/cram_index.o \
	io_lib/cram_stats.o \
	io_lib/cram_samtools.o \
	io_lib/sam_header.o \
	io_lib/hash_table.o \
	io_lib/jenkins_lookup3.o \
	io_lib/vlen.o \
	io_lib/zfio.o \
	io_lib/mFILE.o \
	io_lib/md5.o \
	io_lib/open_trace_file.o \
	io_lib/pooled_alloc.o \
	io_lib/string_alloc.o \
	io_lib/dstring.o \
	io_lib/files.o \
	io_lib/error.o \
	io_lib/xalloc.o \
	io_lib/thread_pool.o

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
	rm io_lib/*.o
