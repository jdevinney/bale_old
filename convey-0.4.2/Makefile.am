# Copyright (c) 2018, Institute for Defense Analyses,
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
#
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.

## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
# (Really we should condition on the compiler, not the programming model)
if !WITH_UPC
  AM_CFLAGS = -Wall -Wno-parentheses -Wno-unused -Wwrite-strings
endif

## Things that will be built and/or installed:

lib_LIBRARIES = libconvey.a
include_HEADERS = convey.h convey_alc8r.h convey_impl.h
if !WITH_UPC
  check_PROGRAMS = alltoallv histogram indexgather
else
  # Because these programs will be linked differently, with $(UPCLINK),
  # we need to give them different names
  check_PROGRAMS = alltoallu histogupc indexgupc
endif
dist_check_SCRIPTS = test_simple test_twohop test_tensor test_elastic

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = convey.pc

## The pieces of the library:

convey_headers = bolite.h convey.h convey_alc8r.h convey_impl.h common.h \
        porter.h porter_impl.h private.h sorter.h tensor.h
convey_sources = circle.c common.c convey.c elastic.c porter.c simple.c \
        tensor.c trivial.c twohop.c
if WITH_SHMEM
  convey_sources += putport.c
endif
if WITH_MPI
  convey_sources += mpiport.c
endif
if WITH_UPC
  convey_sources += putport.c
  convey_headers += mpp2upc.h
  convey_sources += mpp2upc.c
else
if !WITH_MPP_UTIL
  convey_headers += mpp2shmem.h
  convey_sources += mpp2shmem.c
endif
endif

CLEANFILES =
libconvey_a_SOURCES = $(convey_headers)
if WITH_UPC
  nodist_libconvey_a_SOURCES = $(convey_sources:.c=__upc.upc)
  CLEANFILES += *__upc.upc
else
  libconvey_a_SOURCES += $(convey_sources)
endif

# Automatic linking to create UPC copies of C files.
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || $(MKDIR_P) $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; \
  then $(LN_S) $(notdir $<) $@ || cp -f $< $@; \
  else $(LN_S) $< $@ || cp -f $< $@; fi

## Benchmark and example programs:

if !WITH_UPC
  alltoallv_SOURCES = alltoallv.c bolite.h convey.h
  histogram_SOURCES = examples/histo.c examples/example.h bolite.h convey.h
  indexgather_SOURCES = examples/gather.c examples/example.h bolite.h convey.h
  alltoallv_LDADD = libconvey.a
  histogram_LDADD = libconvey.a
  indexgather_LDADD = libconvey.a
else
  alltoallu_SOURCES = bolite.h convey.h
  histogupc_SOURCES = examples/example.h bolite.h convey.h
  indexgupc_SOURCES = examples/example.h bolite.h convey.h
  nodist_alltoallu_SOURCES = alltoallv__upc.upc
  nodist_histogupc_SOURCES = examples/histo__upc.upc
  nodist_indexgupc_SOURCES = examples/gather__upc.upc 
  alltoallu_LDADD = libconvey.a
  histogupc_LDADD = libconvey.a
  indexgupc_LDADD = libconvey.a
endif

## Testing:

LAUNCH ?= $(srcdir)/launchinfo
TESTS_ENVIRONMENT = LAUNCH="$(LAUNCH)"
TESTS = test_simple test_tensor test_elastic
if !WITH_UPC
  TESTS += test_twohop
endif

## Distribution:

EXTRA_DIST = INSTALL README launchinfo bench_simple.awk bench_tensor.awk \
        include/bitops.h include/micro64.h narrative.h sorter.c
CLEANFILES += tests.out

## Generate documentation using doxygen, if possible:

include m4/doxygen.am
