#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "curve.h"
#include "integers.h"

#define Y(__P) ((__P) + SIZE_WORDS)
#define X(__P) (__P)

uint32_t PGEN[] = {
    0xA01C8972, 0xE2945283, 0x4DCA88C7, 0x988B4717, 0x494776FB, 0xBBD1BA39,
    0xB4CEB08C, 0x47DA304D, 0x93B205E6, 0x43709584, 0x01841CA4, 0x60248048,
    0x0012D5D4, 0xAC9CA297, 0xF8103FE4, 0x82189631, 0x59923FBC, 0x026EB7A8,
    /***********************************************************************/
    0x3EF1C7A3, 0x01CD4C14, 0x591984F6, 0x320430C8, 0x7BA7AF1B, 0xB620B01A,
    0xF772AEDC, 0x4FBEBBB9, 0xAC44AEA7, 0x9D4979C0, 0x006D8A2C, 0xFFC61EFC,
    0x9F307A54, 0x4DD58CEC, 0x3BCA9531, 0x4F4AEADE, 0x7F4FBF37, 0x0349DC80};

int isValidPoint(const word *P) {
  /* y^2+xy=x^3+1 */
  word t[SIZE_WORDS] = {0}, l[SIZE_WORDS] = {0};

  polyAdd(t, X(P), Y(P));
  polyMul(t, t, Y(P));

  polySqr(l, X(P));
  polyMul(l, l, X(P));
  l[0] ^= 1;

  return !memcmp(t, l, SIZE_BYTES);
}

word *pointDup(const word *P) {
  return memcpy(malloc(SIZE_WORDS2 * sizeof(word)), P, SIZE_BYTES2);
}

word *pointGen() {
  word *P = malloc(SIZE_WORDS2 * sizeof(word));
  return memcpy(P, PGEN, SIZE_BYTES2);
}

word *pointGenZero() {
  word *pt = calloc(SIZE_WORDS2, sizeof(word));
  pt[SIZE_WORDS - 1] = ((word)1) << (WORDSIZE - 1);
  return pt;
}

int isZero(const word *P) { return (int)(P[SIZE_WORDS - 1] >> (WORDSIZE - 1)); }

word *pointZero(word *P) {
  memset(P, 0, SIZE_BYTES2);
  P[SIZE_WORDS - 1] = ((word)1) << (WORDSIZE - 1);
  return P;
}

word *pointDbl(word *R, const word *P) {
  word t1[SIZE_WORDS] = {0};
  word t2[SIZE_WORDS];

  if (isZero(P))
    return (R == P) ? R : memcpy(R, P, SIZE_BYTES2);
  else if (!memcmp(t1, P, SIZE_BYTES))
    return pointZero(R);

  polyDiv(Y(R), Y(P), X(P));
  polyAddTo(Y(R), X(P)); // Y(R) contains l = y/x + x 
  R[SIZE_WORDS] ^= 1; // Y(R) contains l+1 
  polySqr(t2, X(P)); // t2 contains x^2 
  polyInv(t1, t2); // t1 contains inverse of x^2
  polyAddTo(t1, t2); // t1 contains x'
  polyMul(Y(R), Y(R), t1); // Y(R) contains (l+1)*x'
  polyAddTo(Y(R), t2); // Y(R) contains y' = (l+1)*x' + x^2

  return memcpy(R, t1, SIZE_BYTES);
}

word *pointNeg(word *P) {
  if (isZero(P)) return P;
  polyAddTo(Y(P), X(P));
  return P;
}

word *pointAdd(word *R, const word *P, const word *Q) {
  word l[SIZE_WORDS], x[SIZE_WORDS], y[SIZE_WORDS];
  if (isZero(Q)) {
    if (R == P)
      return R;
    else
      return memcpy(R, P, SIZE_BYTES2);
  } else if (isZero(P)) {
    if (R == Q)
      return R;
    else
      return memcpy(R, Q, SIZE_BYTES2);
  } else {
    if (!memcmp(P, Q, SIZE_BYTES2)) {
      return pointDbl(R, P);
    } else {
      if (!memcmp(X(Q), X(P), SIZE_BYTES)) {
        word mQY[SIZE_WORDS];
        polyAdd(mQY, Y(Q), X(Q));
        if (!memcmp(mQY, Y(P), SIZE_BYTES))
          return pointZero(R);
        else
          memset(x, 0, SIZE_BYTES);
      } else
        polyAdd(x, X(P), X(Q));

      polyAdd(l, Y(P), Y(Q));
      polyDiv(l, l, x);
      polyAddTo(x, polySqr(y, l));
      polyAddTo(x, l);

      polyAdd(y, x, X(P));
      polyMul(y, y, l);
      polyAddTo(y, x);
      polyAdd(Y(R), y, Y(P));

      return memcpy(X(R), x, SIZE_BYTES);
    }
  }
}

word *pointSub(word *R, const word *P, const word *Q) {
  word N[SIZE_WORDS2];
  memcpy(N, Q, SIZE_BYTES2);
  return pointAdd(R, P, pointNeg(N));
}

word* pointSqr(word* P) {
  if (!isZero(P)) {
    polySqr(X(P),X(P));
    polySqr(Y(P),Y(P));
  }
  return P;
}

word *pointMul(word *R, const word *k, const word *P) {
  word m[SIZE_WORDS] = {0}, 
       n[SIZE_WORDS] = {0};
  word Q[SIZE_WORDS2];
  word *a = n, *b = m, *t;
  memcpy(Q, P, SIZE_BYTES2);
  memcpy(n, k, SIZE_BYTES);

  pointZero(R);

  while ( !polyIsZero(n) || !polyIsZero(m) ) {
    if (*a & 1) {
      word w[SIZE_WORDS];
      memcpy(w,b,SIZE_BYTES);
      iLShift(w);
      iNeg(w);
      iAdd(w,a);
      int v = 2 - (int)(*w % 4);
      if (v<0) pointSub(R, R, Q);
      if (v>0) pointAdd(R, R, Q);
      if (iIsNegative(a)) {
        iNeg(a);
        if (v<0) *a -= 1; 
        else iInc(a);
      } else {
        if (v>0) *a -= 1; 
        else iInc(a);
        iNeg(a);
      }
    } else {
      iNeg(a);
    }
    iRShift(a);
    iAdd(b,a);
    pointSqr(Q);
    t = a, a = b, b = t;
  }
  return R;
}

#undef X
#undef Y

word* pointMul_Naive(word* R, const word* k, const word* P) {
    word Q[SIZE_WORDS2];
    memcpy(Q,P,SIZE_BYTES2);
    pointZero(R);
    for (int i=0;i<deg(k);i++) {
        if ((k[i/WORDSIZE]>>(i%WORDSIZE))&1) 
            pointAdd(R,R,Q);
        pointDbl(Q,Q);
    }
    return R;
}
