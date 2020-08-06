# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the copyright holder nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.


## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = $(warning_flags)

## Things that will be built and/or installed:

lib_LIBRARIES = libconvey.a
include_HEADERS = convey.h convey_alc8r.h convey_codec.h biconvey.h
if !WITH_UPC
  check_PROGRAMS = alltoallv histogram indexgather bigather meld
  convey_tuner = alltoallv
else
  # Because these programs will be linked differently, with $(UPCLINK),
  # we need to give them different names
  check_PROGRAMS = alltoallu histogupc indexgupc
  convey_tuner = alltoallu
endif
check_PROGRAMS += wringer
if WITH_SHMEM
  check_PROGRAMS += getput
endif
if WITH_SHMEM_PTR
  check_PROGRAMS += shmemptr
endif
dist_check_SCRIPTS = test_simple test_twohop test_tensor test_elastic \
    test_compress test_examples test_shmem_ptr

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = convey.pc

## Intermediate files produced by cases.pl

codelet_sources = accel_x.c
EXTRA_DIST = accel.c cases.pl
CLEANFILES = $(codelet_sources)

PERL = @PERL@
CASES = @CASES@
%_x.c: %.c
	$(PERL) $(CASES) -E $(srcdir) -p $< > $@


## The pieces of the library:

convey_headers = bolite.h convey.h convey_alc8r.h convey_codec.h convey_impl.h \
    common.h pivot.h porter.h porter_impl.h private.h router.h sorter.h tensor.h
convey_sources = circle.c codata.c common.c convey.c elastic.c packer.c porter.c \
    simple.c squeeze.c tensor.c trivial.c twohop.c $(codelet_sources)
convey_headers += biconvey.h biconvey_impl.h
convey_sources += biconvey.c bitensor.c

if WITH_SHMEM
  convey_sources += putport.c
endif
if WITH_MPI
  convey_sources += mpiport.c
endif
if WITH_UPC
  convey_sources += putport.c
endif

if !WITH_MPP_UTIL
if WITH_MPI
  convey_headers += mpp2mpi.h
  convey_sources += mpp2mpi.c
endif
if WITH_SHMEM
  convey_headers += mpp2shmem.h
  convey_sources += mpp2shmem.c
endif
if WITH_UPC
  convey_headers += mpp2upc.h
  convey_sources += mpp2upc.c
endif
if !WITH_MPI
if !WITH_SHMEM
if !WITH_UPC
  convey_headers += mpp2nil.h
  convey_sources += mpp2nil.c
endif
endif
endif
endif

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
  alltoallv_CFLAGS = $(AM_CFLAGS) $(ignore_deprecation)
  alltoallv_SOURCES = alltoallv.c bolite.h convey.h
  histogram_SOURCES = examples/histo.c examples/example.h bolite.h convey.h
  indexgather_SOURCES = examples/gather.c examples/example.h bolite.h convey.h
  bigather_SOURCES = examples/bigather.c examples/example.h bolite.h convey.h biconvey.h
  meld_SOURCES = examples/meld.c examples/example.h bolite.h convey.h biconvey.h
  alltoallv_LDADD = libconvey.a
  histogram_LDADD = libconvey.a
  indexgather_LDADD = libconvey.a
  bigather_LDADD = libconvey.a
  meld_LDADD = libconvey.a
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

wringer_SOURCES = wringer.c
wringer_LDADD = libconvey.a
if WITH_SHMEM
  getput_SOURCES = getput.c
endif
if WITH_SHMEM_PTR
  shmemptr_SOURCES = shmemptr.c
endif

LAUNCH ?= $(srcdir)/launchinfo
TESTS_ENVIRONMENT = LAUNCH="$(LAUNCH)"
TESTS =
if WITH_SHMEM_PTR
  TESTS += test_shmem_ptr
endif
TESTS += test_simple test_tensor test_elastic test_compress test_examples
if !WITH_UPC
  TESTS += test_twohop
endif

## Tuning:

tune: $(convey_tuner)
	@echo "Please wait a few minutes... (watch tune.log if you're bored)"
	LAUNCH="$(LAUNCH)" $(srcdir)/tune_tensor > tune.log
	$(AWK) '/args =/{b=($$6=="-B"); c=substr($$(b?7:6),3)} /bandw/{print b,c,$$NF*1e-6}' \
          tune.log | sort -k1n -k2n > tune.dat
	$(PERL) $(srcdir)/tune.pl < tune.dat

.PHONY: tune

## Distribution:

EXTRA_DIST += INSTALL README launchinfo bench_simple.awk bench_tensor.awk \
    tune_tensor tune.pl narrative.h sorter.c
CLEANFILES += tests.out tune.log tune.dat

## Generate documentation using doxygen, if possible:

include m4/doxygen.am
