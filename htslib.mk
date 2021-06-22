# Makefile rules useful for third-party code using htslib's public API.
#
#    Copyright (C) 2013-2017, 2019, 2021 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# The makefile fragment included below provides variables that can be used
# to express dependencies on headers supplied by an in-development htslib.
# If your source file foo.c #includes <htslib/hts.h> and <htslib/kstring.h>,
# you can write the correct prerequisites for foo.o as:
#
#	HTSDIR = <path to htslib top-level (build) directory>
#	include $(HTSDIR)/htslib.mk
#
#	foo.o: foo.c $(htslib_hts_h) $(htslib_kstring_h)

HTSSRCDIR = $(HTSDIR)
HTSPREFIX = $(HTSSRCDIR)/
include $(HTSDIR)/htslib_vars.mk

# This file provides the HTSCODECS_SOURCES variable.  It may not be present
# in a freshly checked-out htslib, so is only included if available.  The
# absence is unlikely to cause a problem as there will be plenty of other
# missing files that will trigger a build in htslib, and when that happens
# htslib's makefile will create it.
-include $(HTSDIR)/htscodecs.mk

# Rules for rebuilding an in-development htslib's static and shared libraries.
# If your program foo links with libhts, adding the appropriate prerequisite
# will cause the library to be rebuilt as necessary:
#
#	foo: foo.o $(HTSDIR)/libhts.a
#
# or similarly if your target requires any of the tools supplied:
#
#	bar.bed.bgz.tbi: bar.bed.bgz $(HTSDIR)/tabix
#		$(HTSDIR)/tabix -p bed bar.bed.bgz

HTSLIB_PUBLIC_HEADERS = \
	$(HTSSRCDIR)/htslib/bgzf.h \
	$(HTSSRCDIR)/htslib/cram.h \
	$(HTSSRCDIR)/htslib/faidx.h \
	$(HTSSRCDIR)/htslib/hfile.h \
	$(HTSSRCDIR)/htslib/hts.h \
	$(HTSSRCDIR)/htslib/hts_defs.h \
	$(HTSSRCDIR)/htslib/hts_endian.h \
	$(HTSSRCDIR)/htslib/hts_expr.h \
	$(HTSSRCDIR)/htslib/hts_log.h \
	$(HTSSRCDIR)/htslib/hts_os.h \
	$(HTSSRCDIR)/htslib/kbitset.h \
	$(HTSSRCDIR)/htslib/kfunc.h \
	$(HTSSRCDIR)/htslib/khash.h \
	$(HTSSRCDIR)/htslib/khash_str2int.h \
	$(HTSSRCDIR)/htslib/klist.h \
	$(HTSSRCDIR)/htslib/kseq.h \
	$(HTSSRCDIR)/htslib/ksort.h \
	$(HTSSRCDIR)/htslib/kstring.h \
	$(HTSSRCDIR)/htslib/regidx.h \
	$(HTSSRCDIR)/htslib/sam.h \
	$(HTSSRCDIR)/htslib/synced_bcf_reader.h \
	$(HTSSRCDIR)/htslib/tbx.h \
	$(HTSSRCDIR)/htslib/thread_pool.h \
	$(HTSSRCDIR)/htslib/vcf.h \
	$(HTSSRCDIR)/htslib/vcf_sweep.h \
	$(HTSSRCDIR)/htslib/vcfutils.h

HTSLIB_ALL = \
	$(HTSLIB_PUBLIC_HEADERS) \
	$(HTSSRCDIR)/bcf_sr_sort.c \
	$(HTSSRCDIR)/bcf_sr_sort.h \
	$(HTSSRCDIR)/bgzf.c \
	$(HTSDIR)/config.h \
	$(HTSSRCDIR)/errmod.c \
	$(HTSSRCDIR)/faidx.c \
	$(HTSSRCDIR)/header.c \
	$(HTSSRCDIR)/header.h \
	$(HTSSRCDIR)/hfile_internal.h \
	$(HTSSRCDIR)/hfile.c \
	$(HTSSRCDIR)/hfile_gcs.c \
	$(HTSSRCDIR)/hfile_libcurl.c \
	$(HTSSRCDIR)/hfile_s3.c \
	$(HTSSRCDIR)/hfile_s3_write.c \
	$(HTSSRCDIR)/hts.c \
	$(HTSSRCDIR)/hts_expr.c \
	$(HTSSRCDIR)/hts_internal.h \
	$(HTSSRCDIR)/hts_os.c \
	$(HTSSRCDIR)/kfunc.c \
	$(HTSSRCDIR)/kstring.c \
	$(HTSSRCDIR)/md5.c \
	$(HTSSRCDIR)/multipart.c \
	$(HTSSRCDIR)/plugin.c \
	$(HTSSRCDIR)/probaln.c \
	$(HTSSRCDIR)/realn.c \
	$(HTSSRCDIR)/regidx.c \
	$(HTSSRCDIR)/region.c \
	$(HTSSRCDIR)/sam.c \
	$(HTSSRCDIR)/sam_internal.h \
	$(HTSSRCDIR)/synced_bcf_reader.c \
	$(HTSSRCDIR)/tbx.c \
	$(HTSSRCDIR)/textutils.c \
	$(HTSSRCDIR)/textutils_internal.h \
	$(HTSSRCDIR)/thread_pool.c \
	$(HTSSRCDIR)/thread_pool_internal.h \
	$(HTSSRCDIR)/vcf.c \
	$(HTSSRCDIR)/vcf_sweep.c \
	$(HTSSRCDIR)/vcfutils.c \
	$(HTSSRCDIR)/cram/cram.h \
	$(HTSSRCDIR)/cram/cram_codecs.c \
	$(HTSSRCDIR)/cram/cram_codecs.h \
	$(HTSSRCDIR)/cram/cram_decode.c \
	$(HTSSRCDIR)/cram/cram_decode.h \
	$(HTSSRCDIR)/cram/cram_encode.c \
	$(HTSSRCDIR)/cram/cram_encode.h \
	$(HTSSRCDIR)/cram/cram_external.c \
	$(HTSSRCDIR)/cram/cram_index.c \
	$(HTSSRCDIR)/cram/cram_index.h \
	$(HTSSRCDIR)/cram/cram_io.c \
	$(HTSSRCDIR)/cram/cram_io.h \
	$(HTSSRCDIR)/cram/cram_samtools.h \
	$(HTSSRCDIR)/cram/cram_stats.c \
	$(HTSSRCDIR)/cram/cram_stats.h \
	$(HTSSRCDIR)/cram/cram_structs.h \
	$(HTSSRCDIR)/cram/mFILE.c \
	$(HTSSRCDIR)/cram/mFILE.h \
	$(HTSSRCDIR)/cram/misc.h \
	$(HTSSRCDIR)/cram/open_trace_file.c \
	$(HTSSRCDIR)/cram/open_trace_file.h \
	$(HTSSRCDIR)/cram/os.h \
	$(HTSSRCDIR)/cram/pooled_alloc.c \
	$(HTSSRCDIR)/cram/pooled_alloc.h \
	$(HTSSRCDIR)/cram/string_alloc.c \
	$(HTSSRCDIR)/cram/string_alloc.h \
	$(HTSSRCDIR)/os/lzma_stub.h \
	$(HTSSRCDIR)/os/rand.c \
	$(HTSCODECS_SOURCES)

$(HTSDIR)/config.h:
	+cd $(HTSDIR) && $(MAKE) config.h

$(HTSDIR)/hts-object-files : $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) hts-object-files

$(HTSDIR)/libhts.a: $(HTSDIR)/hts-object-files
	+cd $(HTSDIR) && $(MAKE) lib-static

$(HTSDIR)/libhts.so: $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) lib-shared

$(HTSDIR)/libhts.dylib $(HTSDIR)/libhts.dll.a $(HTSDIR)/hts.dll.a: $(HTSDIR)/hts-object-files
	+cd $(HTSDIR) && $(MAKE) lib-shared

$(HTSDIR)/bgzip: $(HTSSRCDIR)/bgzip.c $(HTSLIB_PUBLIC_HEADERS) $(HTSDIR)/libhts.a
	+cd $(HTSDIR) && $(MAKE) bgzip

$(HTSDIR)/htsfile: $(HTSSRCDIR)/htsfile.c $(HTSLIB_PUBLIC_HEADERS) $(HTSDIR)/libhts.a
	+cd $(HTSDIR) && $(MAKE) htsfile

$(HTSDIR)/tabix: $(HTSSRCDIR)/tabix.c $(HTSLIB_PUBLIC_HEADERS) $(HTSDIR)/libhts.a
	+cd $(HTSDIR) && $(MAKE) tabix

$(HTSDIR)/htslib_static.mk: $(HTSDIR)/htslib.pc.tmp
	+cd $(HTSDIR) && $(MAKE) htslib_static.mk

$(HTSDIR)/htslib.pc.tmp:
	+cd $(HTSDIR) && $(MAKE) htslib.pc.tmp

# Rules for phony targets.  You may wish to have your corresponding phony
# targets invoke these in addition to their own recipes:
#
#	clean: clean-htslib

all-htslib check-htslib clean-htslib distclean-htslib install-htslib mostlyclean-htslib plugins-htslib test-htslib testclean-htslib:
	+cd $(HTSDIR) && $(MAKE) $(@:-htslib=)

.PHONY: all-htslib check-htslib clean-htslib distclean-htslib install-htslib
.PHONY: mostlyclean-htslib plugins-htslib test-htslib testclean-htslib
