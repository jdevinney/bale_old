//
// new bitops.h - v2.x
//
// macros/functions for bit-fiddling

#include <micro64.h>

#ifndef BITOPS_H
#define BITOPS_H

#if (defined(__STDC_VERSION__) && (__STDC_VERSION__ < 199901L))
#error "this header only supports C99 or better"
#endif

//
// get our typedefs and macros for integers
//
#ifdef __cplusplus
#define __STDC_CONSTANT_MACROS
#endif

#include <stdint.h>

// eliminates GCC unused warnings:
#ifdef __GNUC__
#define BITOPS_NOWARN __attribute__ ((unused))
#else
#define BITOPS_NOWARN
#endif


//                       SECTION 1:
//                                                          
//   Inlining and restricted pointers are c99 primitives
//                                                          

#ifndef INLINE
#define INLINE inline
#endif

#ifndef RESTRICT
#ifdef __NVCC__
#define RESTRICT __restrict__
#else
#define RESTRICT restrict
#endif
#endif

//
//  End: SECTION 1


//                       SECTION 2:                         
//                                                          
//        Define CRI instrinsics for all systems            
//                                                          

// most relevant case these days is x86 (Intel and AMD)
// for max performance target your specific host (-march= , -xHost)
// Note: Cray PrgEnv probably sets these for you since their cc driver
//       is like a cross-compiler.
//
// The other option is to use -DBITOPS_HW_POPCNT
// Using both the arch target and -DBITOPS_HW_POPCNT is fine

#ifdef _CRAYC

// Cray C only needs to add _trailz
#include <intrinsics.h>

static INLINE int64_t _trailz (uint64_t x)

{
  return _popcnt( (~x) & (x-1) );
}

#else  // ndef _CRAYC

//                       SECTION 2:
//                                                          
//               Fully-supported intrinsics
//

// Here we set some constants used below. 
#define _MASK1BIT  UINT64_C(0x5555555555555555)
#define _MASK2BITS UINT64_C(0x3333333333333333)
#define _MASK4BITS UINT64_C(0x0f0f0f0f0f0f0f0f)
#define _ZERO64    UINT64_C(0x0)
#define _ONE64     UINT64_C(0xffffffffffffffff)

// MASK macros 
#define _maskl(x)        (((x) == 0) ? _ZERO64   : (_ONE64 << (64-(x))))
#define _maskr(x)        (((x) == 0) ? _ZERO64   : (_ONE64 >> (64-(x))))
#define _mask(x)         (((x) < 64) ? _maskl(x) : _maskr(2*64 - (x)))

// Power/PowerPC 
#if   defined(__PPC__) || defined(__POWERPC__) || defined(__ppc__) || defined(__powerpc__)

#if defined(__IBMC__) || defined(__IBMCPP__)

#define _popcnt(X)  ((int64_t)__popcnt8(X))
#define _poppar(X)  ((int64_t)__poppar8(X))
#define _leadz(X)   ((int64_t)__cntlz8(X))
#define _trailz(X)  ((int64_t)__cnttz8(X))

#elif defined(__GNUC__)

static INLINE int64_t _popcnt (uint64_t x)

{
  int64_t ans;
  
  asm ("popcntd %0, %1" : "=r" (ans) : "r" (x));
  return ans;
}

static INLINE int64_t _poppar (uint64_t x)

{
  return _popcnt(x)&1;
}

static INLINE int64_t _leadz (uint64_t x)

{
  int64_t ans;
  
  asm ("cntlzd %0, %1" : "=r" (ans) : "r" (x));
  return ans;
}

static INLINE int64_t _trailz (uint64_t x)

{
  return _popcnt ((~x) & (x-1));
}

#endif // End: which PowerPC compiler   

// x86-64 with GCC, ICC, or PGI (also PathScale which asserts __GNUC__)
#elif (defined(__x86_64) || defined(__x86_64__)) && (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PGI))

// double shifts, for shift amounts 0..63

#define _dshiftl63 _dshiftl63_fcn
static INLINE uint64_t _dshiftl63_fcn (uint64_t x, uint64_t y, int64_t s)

{
  asm("shld %1,%0" : "+r" (x) : "r" (y), "c" (s));
  return x;
}

#define _dshiftr63 _dshiftr63_fcn
static INLINE uint64_t _dshiftr63_fcn (uint64_t x, uint64_t y, int64_t s)

{
  asm("shrd %1,%0" : "+r" (y) : "r" (x), "c" (s));
  return y;
}

// leadz is a bit tricky to deal with

// PLAN: define _fast_lzcnt and _safe_leadz
//               fast -> lzcnt   safe -> bsr/cmovz

// then define _leadz to be _fast_leadz if __amdfam10__ || __3dNOW__ || __AVX2__
//                    to be _safe_leadz otherwise

static INLINE int64_t _fast_leadz (uint64_t x)

{
  int64_t ans;

  // Nov '09: correct and fast on modern AMD, but not Intel as yet
  // Note: lzcnt uses the same opcode as bsr, but different prefix.
  //       So older chips will actually give incorrect results!
  asm( "lzcnt %1,%0" : "=r" (ans) : "r" (x) );
  return ans;
}

// Nov '04: microcode on Opteron so not that great...
// Apr '07: very good on Intel Core 2
static INLINE int64_t _safe_leadz (uint64_t x)

{
  int64_t ans;

  // use BSR
  asm( "bsr %1,%0" : "=r" (ans) : "r" (x) );
  ans = (x) ? ans : -1;

  return ( 63-ans );
}

// lzcnt works on all modern AMD and Intel >= Haswell/AVX-2
#if __amdfam10__ || __3dNOW__ || __AVX2__
#define _leadz _fast_leadz
#else
#define _leadz _safe_leadz
#endif

#if __AVX2__

// Note: users should _not_ use the function by this name directly:

static INLINE int64_t _bo_avx2_trailz (uint64_t x)

{
  int64_t ans;

  // Nov '09: correct and fast on modern AMD, but not Intel as yet
  // Note: lzcnt uses the same opcode as bsr, but different prefix.
  //       So older chips will actually give incorrect results!
  asm( "tzcnt %1,%0" : "=r" (ans) : "r" (x) );
  return ans;
}

#define _trailz _bo_avx2_trailz

#endif  // End: __AVX2__

// we now prefer the gcc/icc popcount builtin when available

#if __SSE4_2__ || __POPCNT__  // use builtin (can save a cycle vs asm; but s/w builtins are slower)

static INLINE int64_t _popcnt (uint64_t x)

{
  return __builtin_popcountl( x );
}

static INLINE int64_t _poppar (uint64_t x)

{
  return (__builtin_popcountl( x ) & 0x1);
}

#ifndef _trailz
static INLINE int64_t _trailz (uint64_t x)

{
  return _popcnt( (~x) & (x-1) );
}
#endif

#elif BITOPS_HW_POPCNT  // use inline assembly

static INLINE int64_t _popcnt (uint64_t x)

{
  int64_t ans;
  asm ("popcnt %1,%0" : "=r" (ans) : "r" (x));
  return ans;
}

static INLINE int64_t _poppar (uint64_t x)

{
  int64_t ans;
  asm ("popcnt %1,%0" : "=r" (ans) : "r" (x));
  return ( ans & 0x1 );
}

#ifndef _trailz
static INLINE int64_t _trailz (uint64_t x)

{
  return _popcnt( (~x) & (x-1) );
}
#endif

#else  // pure software implementation

static INLINE int64_t _popcnt (uint64_t x)

{
  x = x - ((x >> 1) & _MASK1BIT);
  x = ((x >> 2) & _MASK2BITS) + (x & _MASK2BITS);
  x = ((x >> 4) + x) & _MASK4BITS;
  // preferred on machines with fast 64-bit multiply 
  x = x * UINT64_C(0x0101010101010101);
  return (x >> 56);
}

static INLINE int64_t _poppar (uint64_t x)

{
#if 1
  // Note: gcc may vectorize this code, which is sort of amazing, but can actually _reduce_ performance
  x ^= x << 1;
  x ^= x << 2;
  x = (x & UINT64_C(0x8888888888888888)) * UINT64_C(0x1111111111111111);
  //return (x & _maskl(1)) ? 1 : 0;  // this will inhibit vectorization
  return x >> 63;
#elif 0
  return __builtin_parityl( x );  // ok but not great on gcc; poor on icc
#endif
}

#ifndef _trailz
// enabled as of April 2007
// Nov '04: microcode on Opteron so not that great...
// Apr '07: very good on Intel Core 2
static INLINE int64_t _trailz (uint64_t x)

{
  int64_t ans;

  asm( "bsf %1,%0" : "=r" (ans) : "r" (x) );
  ans = (x) ? ans : 64;

  return ans;
}
#endif

#endif  // End: popcnt / poppar (/ _trailz) implementations

#else // chips without hardware support for special instructions 
      // also includes x86 for compilers other than GCC / ICC

// workarounds for MIPSpro Compilers: Version 7.30 
#if defined(__sgi) && defined(_COMPILER_VERSION) && (_COMPILER_VERSION >= 730)
#define _leadz  _leadzsgi
#define _poppar _popparsgi
#define _popcnt _popcntsgi
#endif 

#ifndef M64_POP_USE_MUL
#define M64_POP_USE_MUL 1
#endif

static INLINE int64_t _popcnt (uint64_t x)

{
  x = x - ((x >> 1) & _MASK1BIT);
  x = ((x >> 2) & _MASK2BITS) + (x & _MASK2BITS);
  x = ((x >> 4) + x) & _MASK4BITS;
#if M64_POP_USE_MUL
  // preferred on machines with fast 64-bit multiply 
  x = x * UINT64_C(0x0101010101010101);
  return x >> 56;
#else
  x = (x >>  8) + x;
  x = (x >> 16) + x;
  x = (x >> 32) + x;
  return x & 0xff;
#endif
}

static INLINE int64_t _poppar (uint64_t x)

{
#if 0  // classic code
  x ^= x >> 32;
  x ^= x >> 16;
  x ^= x >>  8;
  x ^= x >>  4;
  // 4-bit lookup of poppar vals 
  x  = 0x6996 >> (x & 0xf);
  return x&1;
#else  // new code
  x ^= x << 1;
  x ^= x << 2;
  x = (x & UINT64_C(0x8888888888888888)) * UINT64_C(0x1111111111111111);
  return x >> 63;
#endif
}

static INLINE int64_t _leadz (uint64_t x)

{
  // propogate 1's 
  x |= x >>  1;
  x |= x >>  2;
  x |= x >>  4;
  x |= x >>  8;
  x |= x >> 16;
  x |= x >> 32;

  // popcnt(~x) == # of leading zeroes 
  x = ~x - ((~x >> 1) & _MASK1BIT);
  x = ((x >>  2) & _MASK2BITS) + (x & _MASK2BITS);
  x = ((x >>  4) + x) & _MASK4BITS;
#if M64_POP_USE_MUL
  // preferred on machines with fast 64-bit multiply 
  x = x * UINT64_C(0x0101010101010101);
  return x >> 56;
#else
  x = (x >>  8) + x;
  x = (x >> 16) + x;
  x = (x >> 32) + x;
  return x & 0xff;
#endif
}

static INLINE int64_t _trailz (uint64_t x)

{
  x = (~x) & (x-1);

  // popcnt(x) == # of trailing zeroes 
  x = x - ((x >> 1) & _MASK1BIT);
  x = ((x >>  2) & _MASK2BITS) + (x & _MASK2BITS);
  x = ((x >>  4) + x) & _MASK4BITS;
#if M64_POP_USE_MUL
  // preferred on machines with fast 64-bit multiply 
  x = x * UINT64_C(0x0101010101010101);
  return x >> 56;
#else
  x = (x >>  8) + x;
  x = (x >> 16) + x;
  x = (x >> 32) + x;
  return x & 0xff;
#endif
}

#endif // End: chips without hardware support for special instructions 

#endif  // def / ndef _CRAYC

// Additions to all (not part of original _CRI)

// byte-swapping
static INLINE uint64_t _byteswap64 (uint64_t x)

{
  x = ( x                                 << 32) | ( x >> 32                                );
  x = ((x & UINT64_C(0x0000ffff0000ffff)) << 16) | ((x >> 16) & UINT64_C(0x0000ffff0000ffff));
  x = ((x & UINT64_C(0x00ff00ff00ff00ff)) <<  8) | ((x >>  8) & UINT64_C(0x00ff00ff00ff00ff));

  return x;
}

// byte-swapping
static INLINE uint32_t _byteswap32 (uint32_t x)

{
  x = ( x                 << 16) | ( x >> 16                );
  x = ((x & 0x00ff00ffuL) <<  8) | ((x >>  8) & 0x00ff00ffuL);

  return x;
}

#ifndef _dshiftl63
#define _dshiftl63(x,y,n)  (((n) ==  0) ? (x) : ((((uint64_t)(x))<<(n))|(((uint64_t)(y))>>(64-(n)))))
#endif

#ifndef _dshiftl
#define _dshiftl(x,y,n)    (((n) == 64) ? (y) : _dshiftl63(x,y,n))
#endif

#ifndef _dshiftr63
#define _dshiftr63(x,y,n)  (((n) ==  0) ? (y) : ((((uint64_t)(y))>>(n))|(((uint64_t)(x))<<(64-(n)))))
#endif

#ifndef _dshiftr
#define _dshiftr(x,y,n)    (((n) == 64) ? (x) : _dshiftr63(x,y,n))
#endif


//
//  End: SECTION 2


//                       SECTION 3:                         
//                                                          
//             Architecture-dependent features              
//                                                          

//                       SECTION 3a:
//                                                            
//                 x86-64-based systems                  
//                                                            

#if defined(__x86_64) || defined(__x86_64__)

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PGI)

////                   multiply high support              

#define _int_mult_upper _int_mult_upper_fcn
static INLINE uint64_t _int_mult_upper_fcn (uint64_t x, uint64_t y)

{
  uint64_t ans;

#ifdef __PGI  // Why doesn't the more concise asm() work??
  asm ("mov %1, %%rax ; mul %2 ; mov %%rdx,%0" :  "=r" (ans) : "r" (x), "r" (y) : "rax", "rdx");
#else
  // better asm
  asm ("mul %2" : "=d" (ans), "+a" (x) : "r" (y));
#endif
  return ans;
}

#define _rtc _rtc_fcn
static INLINE uint64_t _rtc_fcn (void)

{
  uint64_t ans;
  uint32_t hi=0, lo=0;
    
  asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
  ans = hi;
  ans = (ans << 32) | lo;
  return ans;
}

#endif // End: GNUC/INTEL_C/PGI

#endif // x86-64   

//                       SECTION 3b:
//                                                            
//                Power and PowerPC systems                   
//                                                            

#if defined(__PPC__) || defined(__POWERPC__) || defined(__ppc__) || defined(__powerpc__)

#if defined(__IBMC__) || defined(__IBMCPP__)

////                   multiply high support              

#define _int_mult_upper __mulhdu

// Note (Mar '11): on Power 7 mftb seems to run at 512 MHz

#define _rtc _rtc_fcn
static INLINE uint64_t _rtc_fcn (void)
  
{
  uint64_t ret;
  
  __fence();
  ret = __mftb();
  __fence();

  return ret;
}

#elif defined(__GNUC__)

////                   multiply high support              

#define _int_mult_upper _int_mult_upper_fcn
static INLINE uint64_t _int_mult_upper_fcn (uint64_t x, uint64_t y)

{
  uint64_t ans;

  asm ("mulhdu %0,%1,%2" : "=r" (ans) : "r" (x) , "r" (y));
  return ans;
}

#define _rtc _rtc_fcn
static INLINE uint64_t _rtc_fcn (void)
  
{
  uint64_t ret;

  asm volatile ("mftb %0" : "=r" (ret));

  return ret;
}

#endif // End: which Power compiler   

#endif // Power, PowerPC   

//                       SECTION 3c:
//                                                            
//                          Mips
//                                                            

////                   multiply high support              
#if defined (__mips64)

#if (defined(__PATHCC__) || defined(__GNUC__))
#define _int_mult_upper _int_mult_upper_fcn
static INLINE uint64_t _int_mult_upper_fcn (uint64_t x, uint64_t y)

{
  uint64_t ans;

  asm("dmultu %1,%2; mfhi %0" : "=r" (ans) : "r" (x), "r" (y));
  return ans;
}
#endif // End: pathcc/gcc
  
#endif // End: Mips


//
//  End: SECTION 3


//                       SECTION 4:                           
//                                                            
//           fast 64-bit random number generator              
//             (reasonably free of properties)                

#define _BR_RUNUP_      128
#define _BR_LG_TABSZ_     7
#define _BR_TABSZ_      (INT64_C(1)<<_BR_LG_TABSZ_)

typedef struct {
  uint64_t hi, lo, ind;
  uint64_t tab[_BR_TABSZ_];
} brand_t;

#define _BR_64STEP_(H,L,A,B) {                  \
    uint64_t x;                                 \
    x = H ^ (H << A) ^ (L >> (64-A));           \
    H = L | (x >> (B-64));                      \
    L = x << (128 - B);                         \
  }

static INLINE uint64_t brand (brand_t *p) {
  uint64_t hi=p->hi, lo=p->lo, i=p->ind, ret;

  ret = p->tab[i];

  // 64-step a primitive trinomial LRS:  0, 45, 118   
  _BR_64STEP_(hi,lo,45,118);

  p->tab[i] = ret + hi;

  p->hi  = hi;
  p->lo  = lo;
  p->ind = hi & _maskr(_BR_LG_TABSZ_);

  return ret;
}

static INLINE double dbrand (brand_t *p)

{
  uint64_t x;

#if 1
  // new preferred version
  const double n = 0x1.0p-52;  // C99 hexadecimal floating-point
  // const double n = 2.2204460492503130808e-16;  // decimal floating-point equivalent
  x = brand(p) & _maskr(52);
  return (x * n);
#elif 0
  // !!! only works for IEEE !!!   
  union {uint64_t x; double d;} u;
  u.x = brand(p) & _maskr(52);
  u.x |= INT64_C(1023) << 52;
  return u.d-1.;
#endif
}

static BITOPS_NOWARN 
void brand_init (brand_t *p, uint64_t val)

{
  int64_t i;
  uint64_t hi, lo;

  hi = UINT64_C(0x9ccae22ed2c6e578) ^ val;
  lo = UINT64_C(0xce4db5d70739bd22) & _maskl(118-64);

  // we 64-step 0, 33, 118 in the initialization   
  for (i = 0; i < 64; i++)
    _BR_64STEP_(hi,lo,33,118);
  
  for (i = 0; i < _BR_TABSZ_; i++) {
    _BR_64STEP_(hi,lo,33,118);
    p->tab[i] = hi;
  }
  p->ind = _BR_TABSZ_/2;
  p->hi  = hi;
  p->lo  = lo;

  for (i = 0; i < _BR_RUNUP_; i++)
    brand(p);
}

//
//  End: SECTION 4


//                       SECTION 5:                           
//                                                            
//       portable set of prefetching macros
//  
// Syntax is _prefetch[_SUFFIX] where SUFFIX can be any
// combination of f, n, x in that order.
//   f -> allow page faults
//   n -> non-temporal
//   x -> exclusive/modify intent
//
// Not all machines support all variants, so we will perform replacements
// as necessary, e.g. _pretetch_fx -> _prefetch_x.
//
// Care must be taken with 'f' in particular since prefetching off
// the end of any array in this case could cause SEGV.
//
// Note: most evidence to date suggests exclusive has no effect,
//   or a negative effect on codes.  We have it here because most
//   machines support it, but also to allow for easier experimentation
//   with this type of prefetch.

// Knights Corner
#if   defined(__KNC__) && (defined(__GNUC__) || defined(__INTEL_COMPILER))

#define _prefetch(a)     __asm__ __volatile__ ("vprefetch0 (%0)"     :: "r" ((void *)(a)))
#define _prefetch_n(a)   __asm__ __volatile__ ("vprefetchnta (%0)"   :: "r" ((void *)(a)))
#define _prefetch_x(a)   __asm__ __volatile__ ("vprefetche0 (%0)"    :: "r" ((void *)(a)))
#define _prefetch_nx(a)  __asm__ __volatile__ ("vprefetchenta (%0)"  :: "r" ((void *)(a)))

//
// x86-64/GCC,Intel C
//
#elif (defined(__x86_64) || defined(__x86_64__)) && (defined(__GNUC__) || defined(__INTEL_COMPILER))
#define _prefetch(a)    __asm__ __volatile__ ("prefetcht0 (%0)"  :: "r" ((void *)(a)))
#define _prefetch_n(a)  __asm__ __volatile__ ("prefetchnta (%0)" :: "r" ((void *)(a)))  
// note: prefetchw is not available on all x86-64 CPUs
#if 0
#define _prefetch_x(a)  __asm__ __volatile__ ("prefetchw (%0)"   :: "r" ((void *)(a)))
#endif

// Power, PowerPC

#elif defined(__PPC__) || defined(__POWERPC__) || defined(__ppc__) || defined(__powerpc__)

#if defined(__IBMC__) || defined(__IBMCPP__)

#define _prefetch   __dcbt
#define _prefetch_x __dcbtst

#if defined(__IBM_GCC_ASM) && (__IBM_GCC_ASM > 99)  // Disabled: As of Mar '11 this breaks xlc/as badly

#define _prefetch_n(a)							\
  __asm__ __volatile__ ("dcbt 0,%0,16"    :: "r" ((void *)(a)))
#define _prefetch_nx(a)							\
  __asm__ __volatile__ ("dcbtst 0,%0,16"  :: "r" ((void *)(a)))

#endif  // End: xlc w/ GNU asm

#elif defined(__GNUC__)

// Notes on last argument to prefetches:
//   ,0  -> into L1
//   ,2  -> into L2
//   ,16 -> transient
//   plus other ones related to streams

#define _prefetch(a)							\
  __asm__ __volatile__ ("dcbt 0,%0,0"     :: "r" ((void *)(a)))
#define _prefetch_n(a)							\
  __asm__ __volatile__ ("dcbt 0,%0,16"    :: "r" ((void *)(a)))
#define _prefetch_x(a)							\
  __asm__ __volatile__ ("dcbtst 0,%0,0"   :: "r" ((void *)(a)))
#define _prefetch_nx(a)							\
  __asm__ __volatile__ ("dcbtst 0,%0,16"  :: "r" ((void *)(a)))

#endif // End: which Power compiler

#endif  // End: Arch/Compiler combinations

//
// Fill in gaps for missing variants.  I consider "fault"
// to be the most esoteric variant, so drop it first.
//
#ifndef _prefetch
#define _prefetch
#endif

#ifndef _prefetch_f
#define _prefetch_f _prefetch
#endif

#ifndef _prefetch_n
#define _prefetch_n _prefetch
#endif

#ifndef _prefetch_x
#define _prefetch_x _prefetch
#endif

#ifndef _prefetch_fn
#define _prefetch_fn _prefetch_n
#endif

#ifndef _prefetch_fx
#define _prefetch_fx _prefetch_x
#endif

#ifndef _prefetch_nx
#define _prefetch_nx _prefetch_n
#endif

#ifndef _prefetch_fnx
#define _prefetch_fnx _prefetch_nx
#endif

//
//  End: SECTION 5

//                       SECTION 6:                           
//
// routines for accelerating division with fixed divisor
// requires _int_mult_upper
// currently we only allow for 63-bit numerators but this
// restriction may be lifted in the future
//
// Usage:
//   first call _divbymul_prep() to precompute for division by d
//   then call _divbymul63 to quickly compute x/d where x <= 2^63

#ifdef _int_mult_upper  // requires multiply high

typedef struct {
  int64_t shift;
  uint64_t magic;
} divbymul_t;

static INLINE
uint64_t _divbymul63       (uint64_t numerator, divbymul_t div)

{
  // note: only works if numerator <= 2^63
  return _int_mult_upper (numerator, div.magic) >> div.shift;
}

static BITOPS_NOWARN
int64_t  _divbymul_prep    (uint64_t divisor, divbymul_t *div)

{
  int64_t k, shift;
  uint64_t magic;
  
  if (divisor <= 1) {
    // we set these just to eliminate potential compiler warnings
    // [also might help clue the user in that there is a problem...]
    div->magic = 0;
    div->shift = 0;
    return 1;  // invalid
  }
  
  shift = 63 - _leadz( divisor-1 );
  magic = (~UINT64_C(0)) / divisor;
  for (k = 1; k <= shift; k++) {
    magic <<= 1;
    magic += _int_mult_upper( magic+1, divisor ) < (UINT64_C(1) << k);
  }
  magic += _int_mult_upper( magic, divisor ) < (UINT64_C(1) << k);

  div->magic = magic;
  div->shift = shift;
  return 0;
}

#endif // def _int_mult_upper

//
//  End: SECTION 6

//                       SECTION 7:
//
// timing / timer routines based on gettimeofday()
// Note: we expect gettimeofday() to be more expensive than _rtc or actimer_*
//
// btimer_wall(void) - gettimeofday() expressed in seconds
// btimer_init(void) - init timer to now and return it
// btimer_elapsed(*) - return time since last timer reset/init
// btimer_update(*)  - compute elapsed time and return it, after resetting timer to now

#include <string.h>
#include <sys/time.h>

static INLINE
double btimer_wall       (void)
     
{
  struct timeval tp;
  memset( &tp, 0, sizeof(tp));

  gettimeofday( &tp, NULL );
  return tp.tv_sec + tp.tv_usec/1.0e6;
}

typedef struct { double epoch; } btimer_t;

static BITOPS_NOWARN
btimer_t btimer_init     (void)

{
  btimer_t bt;
  bt.epoch = btimer_wall( );
  return bt;
}

static BITOPS_NOWARN
double btimer_elapsed    (btimer_t *bt)

{
  return btimer_wall( ) - bt->epoch;
}

static BITOPS_NOWARN
double btimer_update     (btimer_t *bt)

{
  double then = bt->epoch;
  double now  = btimer_wall( );

  bt->epoch = now;
  return now - then;
}

//
//  End: SECTION 7


#endif // ndef BITOPS_H   


