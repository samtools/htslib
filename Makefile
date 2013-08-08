# Makefile for htslib, a C library for high-throughput sequencing data formats.
#
#    Copyright (C) 2013 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>

CC     = gcc
AR     = ar
RANLIB = ranlib

CPPFLAGS = -I.
CFLAGS   = -g -Wall -O2 -Wc++-compat
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS  =
LDLIBS   =

prefix      = /usr/local
exec_prefix = $(prefix)
includedir  = $(prefix)/include
libdir      = $(exec_prefix)/lib
mandir      = $(prefix)/share/man
man5dir     = $(mandir)/man5

INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644


all: lib-static lib-shared test/test-vcf-api

HTSPREFIX =
include htslib_vars.mk

lib-static: libhts.a

# $(shell), :=, and ifeq/.../endif are GNU Make-specific.  If you don't have
# GNU Make, comment out the parts of this conditional that don't apply.
PLATFORM := $(shell uname -s)
ifeq "$(PLATFORM)" "Darwin"
SHLIB_FLAVOUR = dylib
lib-shared: libhts.dylib
else
SHLIB_FLAVOUR = so
lib-shared: libhts.so
endif


PACKAGE_VERSION  = 0.0.1
LIBHTS_SOVERSION = 0


# $(NUMERIC_VERSION) is for items that must have a numeric X.Y.Z string
# even if this is a dirty or untagged Git working tree.
NUMERIC_VERSION = $(PACKAGE_VERSION)

# If building from a Git repository, replace $(PACKAGE_VERSION) with the Git
# description of the working tree: either a release tag with the same value
# as $(PACKAGE_VERSION) above, or an exact description likely based on a tag.
# Much of this is also GNU Make-specific.  If you don't have GNU Make and/or
# are not building from a Git repository, comment out this conditional.
ifneq "$(wildcard .git)" ""
original_version := $(PACKAGE_VERSION)
PACKAGE_VERSION := $(shell git describe --always --dirty)

# Unless the Git description matches /\d*\.\d*(\.\d*)?/, i.e., is exactly a tag
# with a numeric name, revert $(NUMERIC_VERSION) to the original version number
# written above, but with the patchlevel field bumped to 255.
ifneq "$(subst ..,.,$(subst 0,,$(subst 1,,$(subst 2,,$(subst 3,,$(subst 4,,$(subst 5,,$(subst 6,,$(subst 7,,$(subst 8,,$(subst 9,,$(PACKAGE_VERSION))))))))))))" "."
empty :=
NUMERIC_VERSION := $(subst $(empty) ,.,$(wordlist 1,2,$(subst ., ,$(original_version))) 255)
endif

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define HTS_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES: .c .o .pico

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

.c.pico:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EXTRA_CFLAGS_PIC) -c -o $@ $<


LIBHTS_OBJS = \
	knetfile.o \
	kstring.o \
	bgzf.o \
	faidx.o \
	hts.o \
	razf.o \
	sam.o \
	synced_bcf_reader.o \
	tbx.o \
	vcf.o \
	vcfutils.o


libhts.a: $(LIBHTS_OBJS)
	@-rm -f $@
	$(AR) -rc $@ $(LIBHTS_OBJS)
	-$(RANLIB) $@


# The target here is libhts.so, as that is the built file that other rules
# depend upon and that is used when -lhts appears in other program's recipes.
# As a byproduct invisible to make, libhts.so.NN is also created, as it is the
# file used at runtime (when $LD_LIBRARY_PATH includes the build directory).

libhts.so: $(LIBHTS_OBJS:.o=.pico)
	$(CC) -shared -Wl,-soname,libhts.so.$(LIBHTS_SOVERSION) -o $@ $(LIBHTS_OBJS:.o=.pico) -lz
	ln -sf $@ libhts.so.$(LIBHTS_SOVERSION)

# Similarly this also creates libhts.NN.dylib as a byproduct, so that programs
# when run can find this uninstalled shared library (when $DYLD_LIBRARY_PATH
# includes this project's build directory).

libhts.dylib: $(LIBHTS_OBJS)
	$(CC) -dynamiclib -install_name $(libdir)/libhts.$(LIBHTS_SOVERSION).dylib -current_version $(NUMERIC_VERSION) -compatibility_version $(LIBHTS_SOVERSION) -o $@ $(LIBHTS_OBJS) -lz
	ln -sf $@ libhts.$(LIBHTS_SOVERSION).dylib


bgzf.o bgzf.pico: bgzf.c config.h $(htslib_hts_h) $(htslib_bgzf_h) htslib/knetfile.h htslib/khash.h
kstring.o kstring.pico: kstring.c htslib/kstring.h
knetfile.o knetfile.pico: knetfile.c htslib/knetfile.h
hts.o hts.pico: hts.c version.h $(htslib_hts_h) $(htslib_bgzf_h) htslib/khash.h htslib/kseq.h htslib/ksort.h
vcf.o vcf.pico: vcf.c $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_tbx_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
sam.o sam.pico: sam.c $(htslib_sam_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
tbx.o tbx.pico: tbx.c $(htslib_tbx_h) htslib/khash.h
faidx.o faidx.pico: faidx.c config.h $(htslib_bgzf_h) $(htslib_faidx_h) htslib/khash.h htslib/knetfile.h
razf.o razf.pico: razf.c $(htslib_razf_h)
synced_bcf_reader.o synced_bcf_reader.pico: synced_bcf_reader.c $(htslib_synced_bcf_reader_h) htslib/kseq.h
vcfutils.o vcfutils.pico: vcfutils.c $(htslib_vcfutils_h)


check test: test/test-vcf-api
	cd test && ./test.pl

test/test-vcf-api: test/test-vcf-api.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/test-vcf-api.o libhts.a $(LDLIBS) -lz

test/test-vcf-api.o: test/test-vcf-api.c $(htslib_hts_h) $(htslib_vcf_h) htslib/kstring.h


install: installdirs install-$(SHLIB_FLAVOUR)
	$(INSTALL_DATA) htslib/*.h $(DESTDIR)$(includedir)/htslib
	$(INSTALL_DATA) libhts.a $(DESTDIR)$(libdir)/libhts.a
	$(INSTALL_DATA) *.5 $(DESTDIR)$(man5dir)

installdirs:
	mkdir -p $(DESTDIR)$(includedir)/htslib $(DESTDIR)$(libdir) $(DESTDIR)$(man5dir)

# After installation, the real file in $(libdir) will be libhts.so.X.Y.Z,
# with symlinks libhts.so (used via -lhts during linking of client programs)
# and libhts.so.NN (used by client executables at runtime).

install-so: libhts.so installdirs
	$(INSTALL_DATA) libhts.so $(DESTDIR)$(libdir)/libhts.so.$(PACKAGE_VERSION)
	ln -sf libhts.so.$(PACKAGE_VERSION) $(DESTDIR)$(libdir)/libhts.so
	ln -sf libhts.so.$(PACKAGE_VERSION) $(DESTDIR)$(libdir)/libhts.so.$(LIBHTS_SOVERSION)

install-dylib: libhts.dylib installdirs
	$(INSTALL_PROGRAM) libhts.dylib $(DESTDIR)$(libdir)/libhts.$(PACKAGE_VERSION).dylib
	ln -sf libhts.$(PACKAGE_VERSION).dylib $(DESTDIR)$(libdir)/libhts.dylib
	ln -sf libhts.$(PACKAGE_VERSION).dylib $(DESTDIR)$(libdir)/libhts.$(LIBHTS_SOVERSION).dylib


mostlyclean:
	-rm -f *.o *.pico test/*.o test/*.dSYM version.h

clean: mostlyclean clean-$(SHLIB_FLAVOUR)
	-rm -f libhts.a test/test-vcf-api

distclean: clean
	-rm -f TAGS

clean-so:
	-rm -f libhts.so libhts.so.*

clean-dylib:
	-rm -f libhts.dylib libhts.*.dylib


tags:
	ctags -f TAGS *.[ch] htslib/*.h


force:


.PHONY: all check clean distclean force install installdirs
.PHONY: lib-shared lib-static mostlyclean tags test
.PHONY: clean-so install-so
.PHONY: clean-dylib install-dylib
