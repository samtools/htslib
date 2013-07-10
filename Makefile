# Makefile for htslib, a C library for high-throughput sequencing data formats.
#
#    Copyright (C) 2013 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>

CC     = gcc
AR     = ar
RANLIB = ranlib

CPPFLAGS =
CFLAGS   = -g -Wall -O2 -Wc++-compat
EXTRA_CFLAGS_PIC = -fpic

prefix      = /usr/local
exec_prefix = $(prefix)
includedir  = $(prefix)/include
libdir      = $(exec_prefix)/lib

INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644


all: lib-static lib-shared

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


# TODO Unify infrastructure for $(PACKAGE_VERSION), #define HTS_VERSION, etc
PACKAGE_VERSION  = 0.0.1
LIBHTS_SOVERSION = 0

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
	$(CC) -dynamiclib -install_name $(libdir)/libhts.$(LIBHTS_SOVERSION).dylib -current_version $(PACKAGE_VERSION) -compatibility_version $(LIBHTS_SOVERSION) -o $@ $(LIBHTS_OBJS) -lz
	ln -sf $@ libhts.$(LIBHTS_SOVERSION).dylib


bgzf.o: bgzf.c $(htslib_bgzf_h) htslib/knetfile.h htslib/khash.h
	$(CC) -c $(CFLAGS) -D_USE_KNETFILE -DBGZF_MT -DBGZF_CACHE bgzf.c -o $@

kstring.o: kstring.c htslib/kstring.h
knetfile.o: knetfile.c htslib/knetfile.h
hts.o: hts.c $(htslib_hts_h) $(htslib_bgzf_h) htslib/khash.h htslib/kseq.h htslib/ksort.h
vcf.o: vcf.c $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_tbx_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
sam.o: sam.c $(htslib_sam_h) htslib/khash.h htslib/kseq.h htslib/kstring.h
tbx.o: tbx.c $(htslib_tbx_h) htslib/khash.h
faidx.o: faidx.c $(htslib_faidx_h) htslib/khash.h $(htslib_razf_h) htslib/knetfile.h
razf.o: razf.c $(htslib_razf_h)
synced_bcf_reader.o: synced_bcf_reader.c $(htslib_synced_bcf_reader_h) htslib/kseq.h
vcfutils.o: vcfutils.c $(htslib_vcfutils_h)


install: installdirs install-$(SHLIB_FLAVOUR)
	$(INSTALL_DATA) htslib/*.h $(DESTDIR)$(includedir)/htslib
	$(INSTALL_DATA) libhts.a $(DESTDIR)$(libdir)/libhts.a

installdirs:
	mkdir -p $(DESTDIR)$(includedir)/htslib $(DESTDIR)$(libdir)

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
	-rm -f *.o *.pico *.dSYM

clean: mostlyclean clean-$(SHLIB_FLAVOUR)
	-rm -f libhts.a

distclean: clean
	-rm -f TAGS

clean-so:
	-rm -f libhts.so libhts.so.*

clean-dylib:
	-rm -f libhts.dylib libhts.*.dylib


tags:
	ctags -f TAGS *.[ch] htslib/*.h


.PHONY: all clean distclean install installdirs
.PHONY: lib-shared lib-static mostlyclean tags
.PHONY: clean-so install-so
.PHONY: clean-dylib install-dylib
