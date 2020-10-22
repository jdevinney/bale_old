#/******************************************************************
#
#
#  Copyright(C) 2018, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#  This material may be reproduced by or for the US Government
#  pursuant to the copyright license under the clauses at DFARS
#  252.227-7013 and 252.227-7014.
# 
#
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the copyright holder nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
# 
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
#  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
#  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#  OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# *****************************************************************/ 

common_sources = libgetput.upc
pkconfigdir=@libdir@/pkgconfig

include_HEADERS=libgetput.h
#lib_LIBRARIES = 
lib_LTLIBRARIES = liblibgetput.la

pkconfig_DATA = @PACKAGE_NAME@.pc
liblibgetput_la_CFLAGS = -g -O2
#liblibgetput_la_CFLAGS = -g -O0

if BUILD_UPC
liblibgetput_la_SOURCES = $(common_sources:.c=__upc.upc)
endif

if BUILD_SHMEM
#liblibgetput_la_LINK = $(LINK)
liblibgetput_la_SOURCES = $(common_sources:.upc=__c.c)
liblibgetput_la_CFLAGS += -DUSE_SHMEM
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

