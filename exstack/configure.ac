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

AC_PREREQ([2.60])     #Require Minimum Autotools version 2.60

AC_INIT([exstack],   [1.0.0],   [])
AM_INIT_AUTOMAKE([1.10 foreign no-define tar-ustar])

AC_CONFIG_HEADERS([config.h])

AC_ARG_WITH([upc],
            AS_HELP_STRING([--with-upc], [build a UPC version]),
            [], [with_upc=yes])
AC_ARG_WITH([shmem],
            AS_HELP_STRING([--with-shmem], [build a SHMEM version]),
            [], [with_shmem=no])


AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CC_C99
mpp_use_shmem=no
if test "x$with_shmem" = xyes; then
   with_upc=no
   mpp_use_shmem=yes
fi

if test "x$with_upc" = xyes; then
AM_PROG_UPC
else
  # must define this ugly conditional to keep automake happy
  AM_CONDITIONAL(am__fastdepUPC, [test])
fi

PKG_CHECK_MODULES([libdeps],
  [
        libgetput
  ]
)

AM_CONDITIONAL(BUILD_UPC, [test "x$with_upc" = xyes])
AM_CONDITIONAL(BUILD_SHMEM, [test "x$mpp_use_shmem" = xyes])

AC_CONFIG_FILES([
exstack.pc:exstack.pc.in
Makefile
])

AC_OUTPUT

