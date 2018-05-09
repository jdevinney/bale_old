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

bin_PROGRAMS = histo ig topo randperm permute_matrix #unionfind

###################################################################
histo_sources = histo_src/histo.h 
histo_sources += histo_src/histo.upc 
histo_sources += histo_src/histo_agi.upc 
histo_sources += histo_src/histo_exstack.upc 
histo_sources += histo_src/histo_conveyor.upc 
histo_sources += histo_src/histo_exstack2.upc
#histo_sources += histo_src/histo_exstack2_cyclic.upc
#histo_sources += histo_src/histo_exstack2_goto.upc
#histo_sources += histo_src/histo_exstack_function.upc
#histo_sources += histo_src/histo_exstack2_function.upc

histo_CPPFLAGS = @libdeps_CFLAGS@ -I$(top_srcdir)/histo_src 
histo_LDADD = @libdeps_LIBS@

if BUILD_UPC
histo_SOURCES = $(histo_sources)
histo_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
histo_SOURCES = $(histo_sources:.upc=__c.c)
histo_CPPFLAGS += -DUSE_SHMEM
histo_LINK = $(LINK)
endif

###################################################################
ig_sources = ig_src/ig.h 
ig_sources += ig_src/ig.upc 
ig_sources += ig_src/ig_agi.upc 
ig_sources += ig_src/ig_exstack.upc 
ig_sources += ig_src/ig_exstack2.upc 
ig_sources += ig_src/ig_conveyor.upc
#ig_sources += ig_src/ig_exstack2_cyclic.upc
#ig_sources += ig_src/ig_exstack2_goto.upc
#ig_sources += ig_src/ig_exstack_function.upc

ig_CPPFLAGS = @libdeps_CFLAGS@ -I$(top_srcdir)/ig_src
ig_LDADD = @libdeps_LIBS@

if BUILD_UPC
ig_SOURCES = $(ig_sources)
ig_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
ig_SOURCES = $(ig_sources:.upc=__c.c)
ig_CPPFLAGS += -DUSE_SHMEM
ig_LINK=$(LINK)
endif
###################################################################
topo_sources = topo_src/toposort.h
topo_sources += topo_src/toposort.upc
topo_sources += topo_src/toposort_agi.upc
topo_sources += topo_src/toposort_exstack.upc
topo_sources += topo_src/toposort_exstack2.upc
topo_sources += topo_src/toposort_conveyor.upc
topo_CPPFLAGS = @libdeps_CFLAGS@ -I$(top_srcdir)/topo_src
topo_LDADD = @libdeps_LIBS@

if BUILD_UPC
topo_SOURCES = $(topo_sources)
topo_LINK=$(UPCLINK)
endif

if BUILD_SHMEM
topo_SOURCES = $(topo_sources:.upc=__c.c)
topo_CPPFLAGS += -DUSE_SHMEM
topo_LINK=$(LINK)
endif

###################################################################
randperm_sources = randperm_src/randperm.upc
randperm_CPPFLAGS = @libdeps_CFLAGS@ 
randperm_LDADD = @libdeps_LIBS@

if BUILD_UPC
randperm_SOURCES = $(randperm_sources)
randperm_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
randperm_SOURCES = $(randperm_sources:.upc=__c.c)
randperm_CPPFLAGS += -DUSE_SHMEM
randperm_LINK=$(LINK)
endif

###################################################################
permute_matrix_sources = permute_matrix_src/permute_matrix.upc
permute_matrix_CPPFLAGS = @libdeps_CFLAGS@ 
permute_matrix_LDADD = @libdeps_LIBS@

if BUILD_UPC
permute_matrix_SOURCES = $(permute_matrix_sources)
permute_matrix_LINK=$(UPCLINK)
endif
if BUILD_SHMEM
permute_matrix_SOURCES = $(permute_matrix_sources:.upc=__c.c)
permute_matrix_CPPFLAGS += -DUSE_SHMEM
permute_matrix_LINK=$(LINK)
endif

###################################################################
#unionfind_sources = unionfind_src/unionfind.h
#unionfind_sources += unionfind_src/unionfind.upc
#unionfind_sources += unionfind_src/unionfind_agi.upc
#unionfind_sources += unionfind_src/unionfind_exstack.upc
#unionfind_sources += unionfind_src/unionfind_exstack2.upc
#unionfind_sources += unionfind_src/unionfind_conveyor.upc
#unionfind_CPPFLAGS = @libdeps_CFLAGS@ -I$(top_srcdir)/unionfind_src
#unionfind_LDADD = @libdeps_LIBS@

#if BUILD_UPC
#unionfind_SOURCES = $(unionfind_sources)
#unionfind_LINK=$(UPCLINK)
#endif
#if BUILD_SHMEM
#unionfind_SOURCES = $(unionfind_sources:.upc=__c.c)
#unionfind_CPPFLAGS += -DUSE_SHMEM
#unionfind_LINK=$(LINK)
#endif

###################################################################
# Automatic linking as needed to create C copies of UPC files.               #
%__c.c: %.upc
	@echo "INFO: Linking   C copy of UPC file: $@ $(top_builddir) $(top_srcdir)"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s ../$< $@ || cp -f $< $@; fi


# Automatic linking as needed to create UPC copies of C files.               #
%__upc.upc: %.c
	@echo "INFO: Linking UPC copy of   C file: $@"
	@rm -f $@
	@test -d $(dir $@) || mkdir -p $(dir $@)
	@if test "$(top_srcdir)" = "$(top_builddir)"; then ln -s $(notdir $<) $@ || cp -f $< $@; else ln -s $< $@ || cp -f $< $@; fi
