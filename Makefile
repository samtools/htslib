# Makefile for htslib, a C library for high-throughput sequencing data formats.
#
#    Copyright (C) 2013-2015 Genome Research Ltd.
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

CC     = gcc
AR     = ar
RANLIB = ranlib

CPPFLAGS =
# TODO: probably update cram code to make it compile cleanly with -Wc++-compat
CFLAGS   = -g -Wall -O2
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS  =
LIBS     =

# For now these don't work too well as samtools also needs to know to
# add -lbz2 and -llzma if linking against the static libhts.a library.
# TODO This needs configury and adding to htslib.pc.in.
#
# # Bzip2 support; optionally used by CRAM.
# HAVE_LIBBZ2 := $(shell echo -e "\#include <bzlib.h>\012int main(void){return 0;}" > .test.c && $(CC) $(CFLAGS) $(CPPFLAGS) -o .test .test.c -lbz2 2>/dev/null && echo yes)
# ifeq "$(HAVE_LIBBZ2)" "yes"
# CPPFLAGS += -DHAVE_LIBBZ2
# LIBS     += -lbz2
# endif
#
# # Lzma support; optionally used by CRAM.
# HAVE_LIBLZMA := $(shell echo -e "\#include <lzma.h>\012int main(void){return 0;}" > .test.c && $(CC) $(CFLAGS) $(CPPFLAGS) -o .test .test.c -llzma 2>/dev/null && echo yes)
# ifeq "$(HAVE_LIBLZMA)" "yes"
# CPPFLAGS += -DHAVE_LIBLZMA
# LIBS     += -llzma
# endif

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
includedir  = $(prefix)/include
libdir      = $(exec_prefix)/lib
libexecdir  = $(exec_prefix)/libexec
datarootdir = $(prefix)/share
mandir      = $(datarootdir)/man
man1dir     = $(mandir)/man1
man5dir     = $(mandir)/man5
pkgconfigdir= $(libdir)/pkgconfig

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644
INSTALL_DIR     = $(MKDIR_P) -m 755

# Set by config.mk if plugins are enabled
plugindir =

BUILT_PROGRAMS = \
	bgzip \
	htsfile \
	tabix

BUILT_TEST_PROGRAMS = \
	test/fieldarith \
	test/hfile \
	test/sam \
	test/test-regidx \
	test/test_view \
	test/test-vcf-api \
	test/test-vcf-sweep

all: lib-static lib-shared $(BUILT_PROGRAMS) plugins $(BUILT_TEST_PROGRAMS)

HTSPREFIX =
include htslib_vars.mk


PACKAGE_VERSION  = 1.3
LIBHTS_SOVERSION = 1


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

print-version:
	@echo $(PACKAGE_VERSION)


.SUFFIXES: .bundle .c .o .pico .so

.c.o:
	$(CC) $(CFLAGS) -I. $(CPPFLAGS) -c -o $@ $<

.c.pico:
	$(CC) $(CFLAGS) -I. $(CPPFLAGS) $(EXTRA_CFLAGS_PIC) -c -o $@ $<


LIBHTS_OBJS = \
	kfunc.o \
	knetfile.o \
	kstring.o \
	bgzf.o \
	faidx.o \
	hfile.o \
	hfile_net.o \
	hts.o \
	md5.o \
	regidx.o \
	sam.o \
	synced_bcf_reader.o \
	vcf_sweep.o \
	tbx.o \
	vcf.o \
	vcfutils.o \
	cram/cram_codecs.o \
	cram/cram_decode.o \
	cram/cram_encode.o \
	cram/cram_external.o \
	cram/cram_index.o \
	cram/cram_io.o \
	cram/cram_samtools.o \
	cram/cram_stats.o \
	cram/files.o \
	cram/mFILE.o \
	cram/open_trace_file.o \
	cram/pooled_alloc.o \
	cram/rANS_static.o \
	cram/sam_header.o \
	cram/string_alloc.o \
	cram/thread_pool.o \
	cram/vlen.o \
	cram/zfio.o

PLUGIN_OBJS =

cram_h = cram/cram.h $(cram_samtools_h) $(cram_sam_header_h) $(cram_structs_h) $(cram_io_h) cram/cram_encode.h cram/cram_decode.h cram/cram_stats.h cram/cram_codecs.h cram/cram_index.h $(htslib_cram_h)
cram_io_h = cram/cram_io.h $(cram_misc_h)
cram_misc_h = cram/misc.h cram/os.h
cram_sam_header_h = cram/sam_header.h cram/string_alloc.h cram/pooled_alloc.h $(htslib_khash_h) $(htslib_kstring_h)
cram_samtools_h = cram/cram_samtools.h $(htslib_sam_h) $(cram_sam_header_h)
cram_structs_h = cram/cram_structs.h cram/thread_pool.h cram/string_alloc.h $(htslib_khash_h)
cram_open_trace_file_h = cram/open_trace_file.h cram/mFILE.h
hfile_internal_h = hfile_internal.h $(htslib_hfile_h)
hts_internal_h = hts_internal.h $(htslib_hts_h)


# To be effective, config.mk needs to appear after most Makefile variables are
# set but before most rules appear, so that it can both use previously-set
# variables in its own rules' prerequisites and also update variables for use
# in later rules' prerequisites.

# If your make doesn't accept -include, change this to 'include' if you are
# using the configure script or just comment the line out if you are not.
-include config.mk

# Usually config.h is generated by running configure or config.status,
# but if those aren't used create a default config.h here.
config.h:
	echo '/* Empty config.h generated by Makefile */' > $@


lib-static: libhts.a

# $(shell), :=, and ifeq/.../endif are GNU Make-specific.  If you don't have
# GNU Make, comment out the parts of this conditional that don't apply.
PLATFORM := $(shell uname -s)
ifeq "$(PLATFORM)" "Darwin"
SHLIB_FLAVOUR = dylib
lib-shared: libhts.dylib
BUILT_PLUGINS = $(PLUGIN_OBJS:.o=.bundle)
else
SHLIB_FLAVOUR = so
lib-shared: libhts.so
BUILT_PLUGINS = $(PLUGIN_OBJS:.o=.so)
endif

plugins: $(BUILT_PLUGINS)


libhts.a: $(LIBHTS_OBJS)
	@-rm -f $@
	$(AR) -rc $@ $(LIBHTS_OBJS)
	-$(RANLIB) $@


# The target here is libhts.so, as that is the built file that other rules
# depend upon and that is used when -lhts appears in other program's recipes.
# As a byproduct invisible to make, libhts.so.NN is also created, as it is the
# file used at runtime (when $LD_LIBRARY_PATH includes the build directory).

libhts.so: $(LIBHTS_OBJS:.o=.pico)
	$(CC) -shared -Wl,-soname,libhts.so.$(LIBHTS_SOVERSION) -pthread $(LDFLAGS) -o $@ $(LIBHTS_OBJS:.o=.pico) -lz -lm $(LIBS)
	ln -sf $@ libhts.so.$(LIBHTS_SOVERSION)

# Similarly this also creates libhts.NN.dylib as a byproduct, so that programs
# when run can find this uninstalled shared library (when $DYLD_LIBRARY_PATH
# includes this project's build directory).

libhts.dylib: $(LIBHTS_OBJS)
	$(CC) -dynamiclib -install_name $(libdir)/libhts.$(LIBHTS_SOVERSION).dylib -current_version $(NUMERIC_VERSION) -compatibility_version $(LIBHTS_SOVERSION) $(LDFLAGS) -o $@ $(LIBHTS_OBJS) -lz $(LIBS)
	ln -sf $@ libhts.$(LIBHTS_SOVERSION).dylib


.pico.so:
	$(CC) -shared -Wl,-E -pthread $(LDFLAGS) -o $@ $< $(LIBS)

.o.bundle:
	$(CC) -bundle -Wl,-undefined,dynamic_lookup $(LDFLAGS) -o $@ $< $(LIBS)


bgzf.o bgzf.pico: bgzf.c config.h $(htslib_hts_h) $(htslib_bgzf_h) $(htslib_hfile_h) $(htslib_khash_h)
kstring.o kstring.pico: kstring.c $(htslib_kstring_h)
knetfile.o knetfile.pico: knetfile.c $(htslib_knetfile_h)
hfile.o hfile.pico: hfile.c config.h $(htslib_hfile_h) $(hfile_internal_h) $(hts_internal_h) $(htslib_khash_h)
hfile_irods.o hfile_irods.pico: hfile_irods.c config.h $(hfile_internal_h) $(htslib_hts_h) $(htslib_kstring_h)
hfile_libcurl.o hfile_libcurl.pico: hfile_libcurl.c config.h $(hfile_internal_h) $(htslib_hts_h) $(htslib_kstring_h)
hfile_net.o hfile_net.pico: hfile_net.c config.h $(hfile_internal_h) $(htslib_knetfile_h)
hts.o hts.pico: hts.c config.h version.h $(htslib_hts_h) $(htslib_bgzf_h) $(cram_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_ksort_h) $(hts_internal_h)
vcf.o vcf.pico: vcf.c config.h $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_tbx_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_kstring_h) $(htslib_khash_str2int_h)
sam.o sam.pico: sam.c config.h $(htslib_sam_h) $(htslib_bgzf_h) $(cram_h) $(hts_internal_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_kstring_h)
tbx.o tbx.pico: tbx.c config.h $(htslib_tbx_h) $(htslib_bgzf_h) $(htslib_khash_h)
faidx.o faidx.pico: faidx.c config.h $(htslib_bgzf_h) $(htslib_faidx_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kstring_h)
synced_bcf_reader.o synced_bcf_reader.pico: synced_bcf_reader.c config.h $(htslib_synced_bcf_reader_h) $(htslib_kseq_h) $(htslib_khash_str2int_h) $(htslib_bgzf_h)
vcf_sweep.o vcf_sweep.pico: vcf_sweep.c config.h $(htslib_vcf_sweep_h) $(htslib_bgzf_h)
vcfutils.o vcfutils.pico: vcfutils.c config.h $(htslib_vcfutils_h) $(htslib_kbitset_h)
kfunc.o kfunc.pico: kfunc.c $(htslib_kfunc_h)
regidx.o regidx.pico: regidx.c config.h $(htslib_hts_h) $(htslib_kstring_h) $(htslib_kseq_h) $(htslib_khash_str2int_h) $(htslib_regidx_h)
md5.o md5.pico: md5.c config.h $(htslib_hts_h)
plugin.o plugin.pico: plugin.c config.h $(hts_internal_h) $(htslib_kstring_h)

cram/cram_codecs.o cram/cram_codecs.pico: cram/cram_codecs.c config.h $(cram_h)
cram/cram_decode.o cram/cram_decode.pico: cram/cram_decode.c config.h $(cram_h) cram/os.h $(htslib_hts_h)
cram/cram_encode.o cram/cram_encode.pico: cram/cram_encode.c config.h $(cram_h) cram/os.h $(htslib_hts_h)
cram/cram_external.o cram/cram_external.pico: cram/cram_external.c config.h $(htslib_hfile_h) $(cram_h)
cram/cram_index.o cram/cram_index.pico: cram/cram_index.c config.h $(htslib_hfile_h) $(cram_h) cram/os.h cram/zfio.h $(hts_internal_h)
cram/cram_io.o cram/cram_io.pico: cram/cram_io.c config.h $(cram_h) cram/os.h $(htslib_hts_h) $(cram_open_trace_file_h) cram/rANS_static.h $(htslib_hfile_h) $(htslib_bgzf_h) $(htslib_faidx_h)
cram/cram_samtools.o cram/cram_samtools.pico: cram/cram_samtools.c config.h $(cram_h) $(htslib_sam_h)
cram/cram_stats.o cram/cram_stats.pico: cram/cram_stats.c config.h $(cram_h) cram/os.h
cram/files.o cram/files.pico: cram/files.c config.h $(cram_misc_h)
cram/mFILE.o cram/mFILE.pico: cram/mFILE.c config.h cram/os.h cram/mFILE.h cram/vlen.h
cram/open_trace_file.o cram/open_trace_file.pico: cram/open_trace_file.c config.h cram/os.h $(cram_open_trace_file_h) $(cram_misc_h) $(htslib_hfile_h)
cram/pooled_alloc.o cram/pooled_alloc.pico: cram/pooled_alloc.c config.h cram/pooled_alloc.h
cram/rANS_static.o cram/rANS_static.pico: cram/rANS_static.c config.h cram/rANS_static.h cram/rANS_byte.h
cram/sam_header.o cram/sam_header.pico: cram/sam_header.c config.h $(cram_sam_header_h) cram/string_alloc.h
cram/string_alloc.o cram/string_alloc.pico: cram/string_alloc.c config.h cram/string_alloc.h
cram/thread_pool.o cram/thread_pool.pico: cram/thread_pool.c config.h cram/thread_pool.h
cram/vlen.o cram/vlen.pico: cram/vlen.c config.h cram/vlen.h cram/os.h
cram/zfio.o cram/zfio.pico: cram/zfio.c config.h cram/os.h cram/zfio.h


bgzip: bgzip.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ bgzip.o libhts.a -lz $(LIBS)

htsfile: htsfile.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ htsfile.o libhts.a -lz $(LIBS)

tabix: tabix.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ tabix.o libhts.a -lz $(LIBS)

bgzip.o: bgzip.c config.h $(htslib_bgzf_h) $(htslib_hts_h)
htsfile.o: htsfile.c config.h $(htslib_hfile_h) $(htslib_hts_h) $(htslib_sam_h) $(htslib_vcf_h)
tabix.o: tabix.c config.h $(htslib_tbx_h) $(htslib_sam_h) $(htslib_vcf_h) $(htslib_kseq_h) $(htslib_bgzf_h) $(htslib_hts_h) $(htslib_regidx_h)


# For tests that might use it, set $REF_PATH explicitly to use only reference
# areas within the test suite (or set it to ':' to use no reference areas).
check test: htsfile $(BUILT_TEST_PROGRAMS)
	test/fieldarith test/fieldarith.sam
	test/hfile
	test/sam test/ce.fa test/faidx.fa
	test/test-regidx
	cd test && REF_PATH=: ./test_view.pl
	cd test && REF_PATH=: ./test.pl

test/fieldarith: test/fieldarith.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/fieldarith.o libhts.a -lz $(LIBS)

test/hfile: test/hfile.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/hfile.o libhts.a -lz $(LIBS)

test/sam: test/sam.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/sam.o libhts.a -lz $(LIBS)

test/test-regidx: test/test-regidx.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/test-regidx.o libhts.a -lz $(LIBS)

test/test_view: test/test_view.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/test_view.o libhts.a -lz $(LIBS)

test/test-vcf-api: test/test-vcf-api.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/test-vcf-api.o libhts.a -lz $(LIBS)

test/test-vcf-sweep: test/test-vcf-sweep.o libhts.a
	$(CC) -pthread $(LDFLAGS) -o $@ test/test-vcf-sweep.o libhts.a -lz $(LIBS)

test/fieldarith.o: test/fieldarith.c $(htslib_sam_h)
test/hfile.o: test/hfile.c $(htslib_hfile_h) $(htslib_hts_defs_h)
test/sam.o: test/sam.c $(htslib_sam_h) $(htslib_faidx_h) $(htslib_kstring_h)
test/test-regidx.o: test/test-regidx.c $(htslib_regidx_h)
test/test_view.o: test/test_view.c $(cram_h) $(htslib_sam_h)
test/test-vcf-api.o: test/test-vcf-api.c $(htslib_hts_h) $(htslib_vcf_h) $(htslib_kstring_h) $(htslib_kseq_h)
test/test-vcf-sweep.o: test/test-vcf-sweep.c $(htslib_vcf_sweep_h)


install: libhts.a $(BUILT_PROGRAMS) $(BUILT_PLUGINS) installdirs install-$(SHLIB_FLAVOUR) install-pkgconfig
	$(INSTALL_PROGRAM) $(BUILT_PROGRAMS) $(DESTDIR)$(bindir)
	if test -n "$(BUILT_PLUGINS)"; then $(INSTALL_PROGRAM) $(BUILT_PLUGINS) $(DESTDIR)$(plugindir); fi
	$(INSTALL_DATA) htslib/*.h $(DESTDIR)$(includedir)/htslib
	$(INSTALL_DATA) libhts.a $(DESTDIR)$(libdir)/libhts.a
	$(INSTALL_DATA) htsfile.1 tabix.1 $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) faidx.5 sam.5 vcf.5 $(DESTDIR)$(man5dir)

installdirs:
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(includedir) $(DESTDIR)$(includedir)/htslib $(DESTDIR)$(libdir) $(DESTDIR)$(man1dir) $(DESTDIR)$(man5dir) $(DESTDIR)$(pkgconfigdir)
	if test -n "$(plugindir)"; then $(INSTALL_DIR) $(DESTDIR)$(plugindir); fi

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

# Substitute these pseudo-autoconf variables only at install time
# so that "make install prefix=/prefix/path" etc continue to work.
install-pkgconfig: installdirs
	sed -e 's#@includedir@#$(includedir)#g;s#@libdir@#$(libdir)#g;s#@PACKAGE_VERSION@#$(PACKAGE_VERSION)#g' htslib.pc.in > $(DESTDIR)$(pkgconfigdir)/htslib.pc
	chmod 644 $(DESTDIR)$(pkgconfigdir)/htslib.pc

# A pkg-config file (suitable for copying to $PKG_CONFIG_PATH) that provides
# flags for building against the uninstalled library in this build directory.
htslib-uninstalled.pc: htslib.pc.in
	sed -e 's#@includedir@#'`pwd`'#g;s#@libdir@#'`pwd`'#g' htslib.pc.in > $@


testclean:
	-rm -f test/*.tmp test/*.tmp.*

mostlyclean: testclean
	-rm -f *.o *.pico cram/*.o cram/*.pico test/*.o test/*.dSYM version.h

clean: mostlyclean clean-$(SHLIB_FLAVOUR)
	-rm -f libhts.a $(BUILT_PROGRAMS) $(BUILT_PLUGINS) $(BUILT_TEST_PROGRAMS)

distclean maintainer-clean: clean
	-rm -f config.cache config.h config.log config.mk config.status
	-rm -f TAGS *-uninstalled.pc

clean-so:
	-rm -f libhts.so libhts.so.*

clean-dylib:
	-rm -f libhts.dylib libhts.*.dylib


tags TAGS:
	ctags -f TAGS *.[ch] cram/*.[ch] htslib/*.h

# We recommend libhts-using programs be built against a separate htslib
# installation.  However if you feel that you must bundle htslib source
# code with your program, this hook enables Automake-style "make dist"
# for this subdirectory.  If you do bundle an htslib snapshot, please
# add identifying information to $(PACKAGE_VERSION) as appropriate.
# (The wildcards attempt to omit non-exported files (.git*, README.md,
# etc) and other detritus that might be in the top-level directory.)
distdir:
	tar -c *.[ch15] [ILMNRcht]*[ELSbcekmnt] | (cd $(distdir) && tar -x)
	+cd $(distdir) && $(MAKE) distclean

force:


.PHONY: all check clean distclean distdir force
.PHONY: install install-pkgconfig installdirs lib-shared lib-static
.PHONY: maintainer-clean mostlyclean plugins print-version tags test testclean
.PHONY: clean-so install-so
.PHONY: clean-dylib install-dylib
