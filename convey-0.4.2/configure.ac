# Copyright (c) 2018, Institute for Defense Analyses,
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500.
#
# This material may be reproduced by or for the U.S. Government 
# pursuant to the copyright license under the clauses at DFARS 
# 252.227-7013 and 252.227-7014.

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
AC_INIT([convey], [0.4.2])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([simple.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])

# Configure options and variables
AC_ARG_WITH([mpp],
            AS_HELP_STRING([--with-mpp], [build atop mpp_utilV4]),
            [], [with_mpp=no])
AC_ARG_WITH([upc],
            AS_HELP_STRING([--with-upc], [build a UPC version]),
            [], [with_upc=no])
AC_ARG_ENABLE([profiling],
              AS_HELP_STRING([--enable-profiling],
                             [enable mpp_util profiling]),
              [], [enable_profiling=no])
if test "x$with_mpp" = xno; then
    enable_profiling=no
fi
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generation program])

REQUIREMENTS=

# Checks for programs.  (Before we can use AC_CHECK_HEADERS, we
# must use AC_PROG_CPP to find the C preprocessor.)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_RANLIB
if test "x$with_upc" = xyes; then
  AM_PROG_UPC
  AC_PROG_MKDIR_P
  AC_PROG_LN_S
else
  # must define this ugly conditional to keep automake happy
  AM_CONDITIONAL(am__fastdepUPC, [test])
fi
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
    AC_MSG_WARN([doxygen not found; continuing without it])
fi

if test "x$with_mpp" = xyes; then
    ACX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS"
                 CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
                 CC="$PTHREAD_CC"],
                [AC_MSG_ERROR([cannot determine how to compile with pthreads])])
fi

# Define how to push and pop the UPC language
AC_DEFUN([AX_LANG_PUSH_UPC], [AC_LANG_PUSH([C])
ac_compile='$UPC -c $UPCFLAGS $CPPFLAGS conftest.$ac_ext >&5'
ac_link='$UPC -o conftest$ac_exeext $UPCFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
ac_compiler_gnu=no])

AC_DEFUN([AX_LANG_POP_UPC], [AC_LANG_POP([C])])

# Consider finding mpp_util via pkg-config...
# PKG_CHECK_MODULES([MPP_UTIL], [mpp_utilV4], , )
# CPPFLAGS="$CPPFLAGS $MPP_UTIL_CFLAGS"
# LIBS="$MPP_UTIL_LIBS $LIBS"

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADERS([inttypes.h])
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
if test "x$with_upc" = xno; then
    AC_CHECK_HEADERS([stdatomic.h])
    AC_CHECK_TYPES([_Atomic uint64_t], [], [], [[#include <stdatomic.h>
#include <stdint.h>]])
    AC_MSG_CHECKING([whether atomics are usable])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdatomic.h>
#if ATOMIC_LLONG_LOCK_FREE
#endif]], [[_Atomic long long foo; atomic_store(&foo, 1);]])],
        [AC_MSG_RESULT([yes])
         AC_DEFINE([HAVE_USABLE_ATOMICS], [1],
                   [Define to 1 if ATOMIC_LLONG_LOCK_FREE and atomic_store() are usable.])],
        [AC_MSG_RESULT([no])])
fi

# Checks for header files
AC_CHECK_HEADERS([unistd.h], ,
    [AC_MSG_ERROR([Couldn't find unistd.h. This is bad.])])
## Set up mpp_utilV4 in the conventional way for now
if test "x$with_mpp" = xyes; then
    AC_CHECK_HEADERS([micro64.h bitops.h], ,
        [AC_MSG_WARN([using possibly obsolete versions of bitops.h and micro64.h])
CPPFLAGS="$CPPFLAGS -I$srcdir/include"])
    AC_CHECK_HEADERS([actimer.h], ,
        [AC_MSG_ERROR([Couldn't find actimer.h. Update CPPFLAGS and try again.])],
        [#define ACTIMER_MODULE_NAME conftest
#define ACTIMER_MAIN])
    AC_CHECK_HEADERS([mpp_utilV4.h], ,
        [AC_MSG_ERROR([Couldn't find mpp_utilV4.h. Update CPPFLAGS and try again.])])
fi

# Checks for dependencies
AC_CHECK_FUNCS([gethostname])
AC_CHECK_LIB([m], [sqrt], ,
    [AC_MSG_ERROR([Couldn't find the math library. Something is very wrong.])])
if test "x$with_mpp" = xyes; then
    AC_CHECK_LIB([rt], [clock_gettime], ,
        [AC_MSG_ERROR([Couldn't find the real-time library. Something is very wrong.])])
    AC_CHECK_LIB([actimer], [actimer_summarize], ,
        [AC_MSG_ERROR([Couldn't find libactimer. Update LDFLAGS and try again.])])
fi

mpp_use_mpi=no
mpp_use_shmem=no
if test "x$with_mpp" = xyes; then
    REQUIREMENTS="mpp_utilV4, actimer"
    if test "x$enable_profiling" = xyes; then
        AC_CHECK_LIB([z], [gzopen], ,
            [AC_MSG_ERROR([Couldn't find libz. Something is very wrong.])])
        AC_CHECK_LIB([mpp_utilV4_prof], [mpp_alltoall], ,
            [AC_MSG_ERROR([Couldn't find libmpp_utilV4_prof. Update LDFLAGS and try again.])])
        
    else
        AC_CHECK_LIB([mpp_utilV4], [mpp_alltoall], ,
            [AC_MSG_ERROR([Couldn't find libmpp_utilV4. Update LDFLAGS and try again.])])
    fi
    AC_LANG_PUSH([C])
    AC_MSG_CHECKING([for MPI mode])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpp_utilV4.h>
#if !MPP_USE_MPI
#error MPP_USE_MPI is false
#endif]], [return 0;])],
        [AC_MSG_RESULT([yes])
         mpp_use_mpi=yes],
        [AC_MSG_RESULT([no])])
    AC_MSG_CHECKING([for SHMEM mode])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpp_utilV4.h>
#if !MPP_USE_SHMEM
#error MPP_USE_SHMEM is false
#endif]], [return 0;])],
        [AC_MSG_RESULT([yes])
         mpp_use_shmem=yes],
        [AC_MSG_RESULT([no])])
    AC_LANG_POP([C])
elif test "x$with_upc" = xno; then
    mpp_use_shmem=yes
else
    AX_LANG_PUSH_UPC
    AC_CHECK_HEADERS([upc.h])
    AC_CHECK_HEADERS([upc_castable.h], [],
        [AC_DEFINE([HAVE_UPC_CASTABLE_H], [1], [Define if upc_castable.h exists.])],
        [#include <upc.h>])
    AX_LANG_POP_UPC
fi

AS_IF([test "x$mpp_use_shmem" = xyes],
[AC_CHECK_HEADERS([shmem.h])
AC_SEARCH_LIBS([shmemx_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_ALLTOALLV], [1], [Define if shmemx_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_team_alltoallv], [],
    [AC_DEFINE([HAVE_SHMEMX_TEAM_ALLTOALLV], [1], [Define if shmemx_team_alltoallv exists.])])
AC_SEARCH_LIBS([shmemx_putmem_signal], [],
    [AC_DEFINE([HAVE_SHMEMX_PUTMEM_SIGNAL], [1], [Define if shmemx_putmem_signal exists.])])
AC_SEARCH_LIBS([shmem_putmem_nbi], [],
    [AC_DEFINE([HAVE_SHMEM_PUTMEM_NBI], [1], [Define if shmem_putmem_nbi exists.])])
AC_SEARCH_LIBS([shmem_free], [],
    [AC_DEFINE([HAVE_SHMEM_FREE], [1], [Define if shmem_free exists.])])
AC_SEARCH_LIBS([shmem_align], [],
    [AC_DEFINE([HAVE_SHMEM_ALIGN], [1], [Define if shmem_align exists.])])
AC_SEARCH_LIBS([shmem_ptr], [],
    [AC_DEFINE([HAVE_SHMEM_PTR], [1], [Define if shmem_ptr exists.])])
AC_SEARCH_LIBS([shmem_global_exit], [],
    [AC_DEFINE([HAVE_SHMEM_GLOBAL_EXIT], [1], [Define if shmem_global_exit exists.])])
])
AS_IF([test "x$with_mpp" = xyes],
    [AC_DEFINE([HAVE_MPP_UTIL], [1], [Define if the build is based on mpp_utilV4.])])
AS_IF([test "x$with_upc" = xyes],
    [AC_DEFINE([MPP_USE_UPC], [1], [Define if the build is based on UPC.])])
AS_IF([test "x$enable_profiling" = xyes],
    [AC_DEFINE([ENABLE_PROFILING], [1], [Define to activate mpp profiling and actimer.])])

# Define conditions for automake
AM_CONDITIONAL(WITH_UPC, [test "x$with_upc" = xyes])
AM_CONDITIONAL(WITH_MPP_UTIL, [test "x$with_mpp" = xyes])
AM_CONDITIONAL(WITH_MPI, [test "x$mpp_use_mpi" = xyes])
AM_CONDITIONAL(WITH_SHMEM, [test "x$mpp_use_shmem" = xyes])
AM_CONDITIONAL(HAVE_DOXYGEN, [test -n "$DOXYGEN"])

# Produce output
AC_SUBST([VERSION])
AC_SUBST([REQUIREMENTS])
AC_CONFIG_FILES([Makefile convey.pc])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES(Doxyfile)])
AC_OUTPUT

# Report the results
if test "x$with_upc" = xyes; then
    prog_model=UPC
elif test "x$with_mpp" = xno; then
    prog_model=SHMEM
elif test "x$mpp_use_mpi" = xyes; then
    prog_model=MPP_UTIL/MPI
elif test "x$mpp_use_shmem" = xyes; then
    prog_model=MPP_UTIL/SHMEM
else
    prog_model=MPP_UTIL
fi

echo "
     $PACKAGE_NAME $PACKAGE_VERSION
     Prefix.............: $prefix
     Programming Model..: $prog_model
     MPP Profiling......: $enable_profiling
     C Compiler.........: $CC $CFLAGS $CPPFLAGS
     Linker.............: $CC $LDFLAGS $LIBS
"
