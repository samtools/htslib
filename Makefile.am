ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

lib_LTLIBRARIES = libhts.la

libhts_la_SOURCES = \
	bgzf.c \
	faidx.c \
	hfile.c \
	hfile_internal.h \
	hfile_net.c \
	hts.c \
	kfunc.c \
	knetfile.c \
	kstring.c \
	sam.c \
	synced_bcf_reader.c \
	tbx.c \
	vcf.c \
	vcf_sweep.c \
	vcfutils.c \
	cram/cram_codecs.c \
	cram/cram_codecs.h \
	cram/cram_decode.c \
	cram/cram_decode.h \
	cram/cram_encode.c \
	cram/cram_encode.h \
	cram/cram.h \
	cram/cram_index.c \
	cram/cram_index.h \
	cram/cram_io.c \
	cram/cram_io.h \
	cram/cram_samtools.c \
	cram/cram_samtools.h \
	cram/cram_stats.c \
	cram/cram_stats.h \
	cram/cram_structs.h \
	cram/files.c \
	cram/mFILE.c \
	cram/mFILE.h \
	cram/md5.c \
	cram/md5.h \
	cram/misc.h \
	cram/open_trace_file.c \
	cram/open_trace_file.h \
	cram/os.h \
	cram/pooled_alloc.c \
	cram/pooled_alloc.h \
	cram/sam_header.c \
	cram/sam_header.h \
	cram/string_alloc.c \
	cram/string_alloc.h \
	cram/thread_pool.c \
	cram/thread_pool.h \
	cram/vlen.c \
	cram/vlen.h \
	cram/zfio.c \
	cram/zfio.h

libhts_la_CFLAGS = -DSAMTOOLS=1
libhts_la_CPPFLAGS = -I$(top_srcdir) $(ZLIB_CPPFLAGS)
libhts_la_LIBADD = $(ZLIB_LIBS)

bin_PROGRAMS = bgzip tabix

bgzip_SOURCES = bgzip.c
bgzip_LDADD = libhts.la

tabix_SOURCES = tabix.c
tabix_LDADD = libhts.la

pkginclude_HEADERS = \
	htslib/bgzf.h \
	htslib/faidx.h \
	htslib/hfile.h \
	htslib/hts.h \
	htslib/hts_defs.h \
	htslib/kfunc.h \
	htslib/khash.h \
	htslib/khash_str2int.h \
	htslib/klist.h \
	htslib/knetfile.h \
	htslib/kseq.h \
	htslib/ksort.h \
	htslib/kstring.h \
	htslib/sam.h \
	htslib/synced_bcf_reader.h \
	htslib/tbx.h \
	htslib/vcf.h \
	htslib/vcf_sweep.h \
	htslib/vcfutils.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = htslib.pc

dist_man_MANS = tabix.1 faidx.5 sam.5 vcf.5
