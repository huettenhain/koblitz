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

typedef uint8_t  byte;
#define BYTESIZE 0x08

#ifdef BITS64
 typedef uint64_t    word;
 typedef uint32_t    hword;
# define HWSIZE      0x20
# define WORDSIZE    0x40
# define SIZE_WORDS  0x0009
# define SIZE_WORDS2 0x0012
#else
 typedef uint32_t    word;
 typedef uint16_t    hword;
# define HWSIZE      0x10
# define WORDSIZE    0x20
# define SIZE_WORDS  0x0012
# define SIZE_WORDS2 0x0024
#endif 

#define F571        (((word)1)<<(WORDSIZE-5))
#define EMPTY_MASK  (F571-1)

#define SIZE_BYTES  0x0048
#define SIZE_BYTES2 0x0090
#define SIZE_BITS   0x023a
#define SIZE_BITS2  0x0474

void iLShiftN(word *a, word s, word n);
void iRShiftN(word *a, word s, word n);
void iRShift (word *a, word s);
void iLShift (word *a, word s);

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