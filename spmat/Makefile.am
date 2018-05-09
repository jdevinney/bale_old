#/******************************************************************
# * Copyright 2014, Institute for Defense Analyses
# * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# * This material may be reproduced by or for the US Government
# * pursuant to the copyright license under the clauses at DFARS
# * 252.227-7013 and 252.227-7014.
# *
# * POC: Bale <bale@super.org>
# * Please contact the POC before disseminating this code.
# *****************************************************************/ 

common_sources  = spmat_agi.upc 
common_sources += spmat_exstack.upc
common_sources += spmat_exstack2.upc 
common_sources += spmat_conveyor.upc
common_sources += spmat_utils.upc
pkconfigdir=@libdir@/pkgconfig

include_HEADERS=spmat.h

lib_LIBRARIES = libspmat.a
libspmat_a_CPPFLAGS = @CPPFLAGS@ @libdeps_CFLAGS@
pkconfig_DATA = @PACKAGE_NAME@.pc

if BUILD_UPC
libspmat_a_SOURCES = $(common_sources:.c=__upc.upc)
endif

if BUILD_SHMEM
libspmat_a_SOURCES = $(common_sources:.upc=__c.c)
libspmat_a_CPPFLAGS += -DUSE_SHMEM
endif

# Automatic linking as needed to create C copies of UPC files.               #
%__c.c: %.upc
	@echo "INFO: Linking   C copy of UPC file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi


# Automatic linking as needed to create UPC copies of C files.               #
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi


