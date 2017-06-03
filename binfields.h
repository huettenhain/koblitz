#ifndef _BINFIELDS__H
#define _BINFIELDS__H

#ifdef _MSC_VER
 typedef unsigned __int64 uint64_t;
 typedef unsigned __int32 uint32_t;
 typedef unsigned __int16 uint16_t;
 typedef unsigned __int8  uint8_t;
#else
# include <stdint.h>
#endif

#if _WIN64
# define BITS64
#endif
#if __GNUC__
# if __x86_64__ || __ppc64__
#  define BITS64
# endif
#endif

#ifdef BITS64
 typedef  int64_t signed_word; 
 typedef uint64_t word;
#define WORDSIZE  0x40
#else
 typedef  int32_t signed_word; 
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

void iLShiftN(word *a, word s, word n);
void iRShiftN(word *a, word s, word n);

int      deg(const word *a);
int      polyIsZero(word *a);

word*    polyAdd  (word* c, const word* a, const word* b);
word*    polyAddTo(word* a, const word* b);

word*    polySqr(word* c, const word* a);
word*    polyMul(word* c, const word *a, const word *b); 
word*    polyDiv(word* c, const word *a, const word *b);
word*    polyInv(word *c, const word *a);

word*    polyGen();
word*    polyRand(word *p);
word*    polyDup(word* a);

#endif