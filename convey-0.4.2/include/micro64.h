//
// new micro64.h - v2.x
//

// micro64.h v2 leverages c99 to provide most features
// 1. stdint.h for integer typedefs and macros for constants (e.g. UINT64_C)
// 2a. endianness has been modified to avoid namespace colisions with endian.h
//     M64_BIG_ENDIAN, M64_LITTLE_ENDIAN, m64_big_endian( ), m64_little_endian( )
// 2b. precise endianness behavior under control of -DM64_ENDIAN_STANDARD=x
//     0 -> old/legacy behavior, 1 -> new behavior
// 3. optional legacy types (under control of M64_LEGACY_TYPES which defaults to '1')
//    I64/UI64 ; int8, uint32, etc. ; LONGIS32 ; _M64_INT_REGS

#ifndef MICRO64_H
#define MICRO64_H

#if (defined(__STDC_VERSION__) && (__STDC_VERSION__ < 199901L))
#error "this header only supports C99 or better"
#endif

//
// 1. get our typedefs and macros for integers
//
#ifdef __cplusplus
#define __STDC_CONSTANT_MACROS
#endif

#include <stdbool.h>
#include <stdint.h>

//
// 2. sort out endianness
//

// 2a. endian.h (or emulation)
#if defined(__linux__) || defined(__GNU_LIBRARY__)
// Linux / GLIB standard macros for endianness
#include <endian.h>
#elif defined(BSD)
// Note: need to test this
#include <sys/endian.h>
#endif  // End: linux / glib / bsd

#if   defined(BYTE_ORDER) && defined(LITTLE_ENDIAN) && defined(BIG_ENDIAN)
#if   (BYTE_ORDER == LITTLE_ENDIAN)
#define M64_LITTLE_ENDIAN 1
#elif (BYTE_ORDER == BIG_ENDIAN)
#define M64_BIG_ENDIAN 1
#endif  // End: which BYTE_ORDER

#elif defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && defined(__BIG_ENDIAN)
#if   (__BYTE_ORDER == __LITTLE_ENDIAN)
#define M64_LITTLE_ENDIAN 1
#elif (__BYTE_ORDER == __BIG_ENDIAN)
#define M64_BIG_ENDIAN 1
#endif  // End: which __BYTE_ORDER
#endif  // End: BYTE_ORDER / __BYTE_ORDER cases

#if !defined(M64_LITTLE_ENDIAN) && !defined(M64_BIG_ENDIAN)

// big-endian OS's
#if defined(__sgi) || defined(__sparc) || defined(__hpux) || defined(_AIX)
#define M64_BIG_ENDIAN 1

#elif defined(__osf__) || defined(__x86_64) || defined(__i386)
#define M64_LITTLE_ENDIAN 1

#endif  // End: which OS

#endif  // End: ndef M64_LITTLE_ENDIAN && ndef M64_BIG_ENDIAN

// 2b. handle legacy micro64.h endian behavior
// the goal here is to migrate folks to endian.h behavior in three phases:
//    
//    1. must set -DM64_ENDIAN_STANDARD=x to choose between classic(x=0)
//       and new(x!=0) behavior.
//       BIG_ENDIAN / LITTLE_ENDIAN are "poisoned" if M64_ENDIAN_STANDARD not defined
//
//    2. must assert -DM64_ENDIAN_STANDARD=1 as classic behavior will fail with
//       error "<useful info>"
//
//    3. only new behavior allowed, -DM64_ENDIAN_STANDARD=0 fails with #error and
//       otherwise -DM64_ENDIAN_STANDARD=1 is asserted for you
//
//    Phase 1 - 2015-2016 ; Phase 2 - ~2016 ; Phase 3 - ~2017+

#ifdef M64_ENDIAN_STANDARD

#if (0 == M64_ENDIAN_STANDARD)

// classic behavior, but warn
#ifdef __GNUC__
#warning "classic micro64.h implementation of BIG_ENDIAN/LITTLE_ENDIAN - will phase out in ~2016"
#else
#pragma  "classic micro64.h implementation of BIG_ENDIAN/LITTLE_ENDIAN - will phase out in ~2016"
#endif

#undef BIG_ENDIAN
#undef LITTLE_ENDIAN

#if   defined(M64_LITTLE_ENDIAN)
#define LITTLE_ENDIAN
#elif defined(M64_BIG_ENDIAN)
#define BIG_ENDIAN
#else
#error "fatal error for 0 == M64_ENDIAN_STANDARD"
#endif  // End: classic little/big endian behavior

#endif  // End: 0 == M64_ENDIAN_STANDARD

#else  // ndef M64_ENDIAN_STANDARD

// Phase 1 - poison these macros:
#ifdef __GNUC__
#undef LITTLE_ENDIAN  // undef prior to poisoning to avoid _always_ getting warnings
#pragma GCC poison LITTLE_ENDIAN
#undef BIG_ENDIAN     // undef prior to poisoning to avoid _always_ getting warnings
#pragma GCC poison BIG_ENDIAN
#endif

#endif  // ndef / def M64_ENDIAN_STANDARD

#if   defined(M64_LITTLE_ENDIAN)
static inline bool m64_little_endian (void) { return true;  }
static inline bool m64_big_endian    (void) { return false; }

#elif defined(M64_BIG_ENDIAN)
static inline bool m64_little_endian (void) { return false; }
static inline bool m64_big_endian    (void) { return true;  }

#endif  // which endian

// 3. if requested, support some legacy micro64.h typedefs and macros:

// we default this behavior to on:
#ifndef M64_LEGACY_TYPES
#define M64_LEGACY_TYPES 1
#endif

#if M64_LEGACY_TYPES

// 3a. macros for constants
#define  I64  INT64_C
#define UI64  UINT64_C

// 3b. legacy typdefs
typedef   int8_t   int8;
typedef  uint8_t  uint8;
typedef  int16_t  int16;
typedef uint16_t uint16;
typedef  int32_t  int32;
typedef uint32_t uint32;
typedef  int64_t  int64;
typedef uint64_t uint64;

// 3c. LONGIS32
#ifndef __LP64__
#if ( defined(__i386)                                      || \
      (defined(__sparc) && (!defined(__sparcv9)))          || \
      (defined(__hpux)  && (!defined(__LP64__)))           || \
      (defined(__sgi)   && (_MIPS_SIM != _MIPS_SIM_ABI64)) || \
      (defined(_AIX)    && (!defined(__64BIT__)))         )
#define LONGIS32 1
#endif  // End: 32-bit arches
#endif  // End: ndef __LP64__

// 3d. _M64_INT_REGS == number of GP regs
#ifndef _M64_INT_REGS

#if   defined(__i386)  // old x86
#define _M64_INT_REGS 8

#elif defined(__x86_64) || defined(__x86_64__)  // x86_64
#define _M64_INT_REGS 16

#elif defined(__ia64) || defined(__ia64__)  // Itanium
#define _M64_INT_REGS 128

#else
// Default 
#define _M64_INT_REGS 32

#endif  // End: which CPU

#endif  // End: ndef _M64_INT_REGS

#endif  // End: M64_LEGACY_TYPES

#endif  // End: ndef MICRO64_H 

