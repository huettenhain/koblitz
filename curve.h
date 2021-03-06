#ifndef _CURVE_H
#define _CURVE_H

#include "binfields.h"

/* We will use, by default, the Koblitz Curve
   y^2 + xy = x^3 + 1                           */

word*  pointMul(word* R, const word* k, const word* P);
word*  pointMul_Naive(word* R, const word* k, const word* P);

word*  pointAdd(word* R, const word* P, const word* Q);
word*  pointDbl(word* R, const word* P);
word*  pointNeg(word* P);
word*  pointZero(word* P);
word*  pointSqr(word* P);

int    isZero(const word* P);
int    isValidPoint(const word* P);


word*  pointDup(const word* P);
word*  pointGen();
word*  pointGenZero();

#endif