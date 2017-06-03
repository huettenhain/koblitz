#include "binfields.h"
#include "integers.h"

void iAdd(word *a, word *b) {
  word carry = 0;
  for (int i = 0; i < SIZE_WORDS; i++) {
    word tmp = a[i] + b[i],
         cry = (tmp < a[i]);
    a[i] = tmp;
    carry = (carry && !++a[i]) || cry;
  }
}

int iIsNegative(word *a) {
    return (a[SIZE_WORDS-1] >> (WORDSIZE - 1));
}

void iInc(word *a) {
  word j = 0;
  do if (++a[j++]) break; while (j < SIZE_WORDS);
}

void iNeg(word *a) {
  word i;
  for (i = 0; i < SIZE_WORDS; i++) 
    a[i] = ~a[i];
  iInc(a);
}

void iRShift(word *a) {
  word *p = a + SIZE_WORDS - 1;
  while (a < p) {
    *a >>= 1;
    a++;
    *(a - 1) |= (*a & 1) << (WORDSIZE - 1);
  }
  *((signed_word*)a) >>= 1;
}

void iLShift(word *a) {
  word *p = a + SIZE_WORDS - 1;
  while (p > a) {
    *p <<= 1;
    p--;
    *(p + 1) |= *p >> (WORDSIZE - 1);
  }
  *p <<= 1;
}