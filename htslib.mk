# Makefile rules useful for third-party code using htslib's public API.
#
#    Copyright (C) 2013 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>

# The makefile fragment included below provides variables that can be used
# to express dependencies on headers supplied by an in-development htslib.
# If your source file foo.c #includes <htslib/hts.h> and <htslib/kstring.h>,
# you can write the correct prerequisites for foo.o as:
#
#	HTSDIR = <path to htslib top-level directory>
#	include $(HTSDIR)/htslib.mk
#
#	foo.o: foo.c $(htslib_hts_h) $(HTSDIR)/htslib/kstring.h
#
# Variables are not provided for k*.h, as those never include other headers.

HTSPREFIX = $(HTSDIR)/
include $(HTSDIR)/htslib_vars.mk

# Rules for rebuilding an in-development htslib's static and shared libraries.
# If your program foo links with libhts, adding the appropriate prerequisite
# will cause the library to be rebuilt as necessary:
#
#	foo: foo.o $(HTSDIR)/libhts.a

HTSLIB_ALL = \
	$(HTSDIR)/htslib/bgzf.h \
	$(HTSDIR)/htslib/faidx.h \
	$(HTSDIR)/htslib/hts.h \
	$(HTSDIR)/htslib/khash.h \
	$(HTSDIR)/htslib/klist.h \
	$(HTSDIR)/htslib/knetfile.h \
	$(HTSDIR)/htslib/kseq.h \
	$(HTSDIR)/htslib/ksort.h \
	$(HTSDIR)/htslib/kstdint.h \
	$(HTSDIR)/htslib/kstring.h \
	$(HTSDIR)/htslib/razf.h \
	$(HTSDIR)/htslib/sam.h \
	$(HTSDIR)/htslib/synced_bcf_reader.h \
	$(HTSDIR)/htslib/tbx.h \
	$(HTSDIR)/htslib/vcf.h \
	$(HTSDIR)/htslib/vcfutils.h \
	$(HTSDIR)/bgzf.c \
	$(HTSDIR)/faidx.c \
	$(HTSDIR)/hts.c \
	$(HTSDIR)/knetfile.c \
	$(HTSDIR)/kstring.c \
	$(HTSDIR)/razf.c \
	$(HTSDIR)/sam.c \
	$(HTSDIR)/synced_bcf_reader.c \
	$(HTSDIR)/tbx.c \
	$(HTSDIR)/vcf.c \
	$(HTSDIR)/vcfutils.c

$(HTSDIR)/libhts.a: $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) lib-static

$(HTSDIR)/libhts.so $(HTSDIR)/libhts.dylib: $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) lib-shared

# Rules for phony targets.  You may wish to have your corresponding phony
# targets invoke these in addition to their own recipes:
#
#	clean: clean-htslib

clean-htslib install-htslib:
	+cd $(HTSDIR) && $(MAKE) $(@:-htslib=)

.PHONY: clean-htslib install-htslib
