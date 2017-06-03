#ifndef _BINFIELDS__H
#define _BINFIELDS__H

#ifdef _MSC_VER
 typedef   signed __int64  int64_t;
 typedef   signed __int32  int32_t;
 typedef unsigned __int64 uint64_t;
 typedef unsigned __int32 uint32_t;
 typedef unsigned __int16 uint16_t;
 typedef unsigned __int8  uint8_t;
#else
# include <stdint.h>
#endif

#if defined(_WIN64) || defined(_M_X64) || defined(_M_AMD64) || defined(_M_IA64) 
# define BITS64
#endif
#if __GNUC__
# if __x86_64__ || __ppc64__
#  define BITS64
# endif
#endif

#ifdef BITS64
 typedef  int64_t sw; 
 typedef uint64_t word;
#define WORDSIZE  0x40
#else
 typedef  int32_t sw; 
 typedef uint32_t word;
#define WORDSIZE  0x20
#endif 

typedef uint16_t hword;
typedef uint8_t  byte;
#define BYTESIZE 0x08
#define HWSIZE   (2*BYTESIZE)
#define HW_PER_W (WORDSIZE / HWSIZE)

#define F571        (((word)1)<<(WORDSIZE-5))
#define EMPTY_MASK  (F571-1)
#define SIZE_BITS   576
#define SIZE_BYTES  (SIZE_BITS/BYTESIZE)
#define SIZE_WORDS  (SIZE_BITS/WORDSIZE)
#define SIZE_WORDS2 (2*SIZE_WORDS)
#define SIZE_BYTES2 (2*SIZE_BYTES)

void     polyLShift(word *a, unsigned long size, unsigned long places);
void     polyRShift(word *a, unsigned long size, unsigned long places);

unsigned int
        polyDeg(const word *a);

int     polyIsZero(const word *a);

word*   polyAdd  (word* c, const word* a, const word* b);
word*   polyAddTo(word* a, const word* b);

word*   polySqr(word* c, const word* a);
word*   polyMul(word* c, const word *a, const word *b); 
word*   polyDiv(word* c, const word *a, const word *b);
word*   polyInv(word *c, const word *a);

word*   polyGen();
word*   polyRand(word *p);
word*   polyDup(word* a);

#endif