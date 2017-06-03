
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "binfields.h"
#include "curve.h"

void validtest() {
  word *x = polyGen(), *y = polyGen();
  word *q = pointGen(), *p = pointGen();
  word *s = pointGenZero(), *t = pointGenZero();

  // Check if the generator is Valid
  assert(isValidPoint(q));

  memset(x, 0, SIZE_BYTES);
  x[0] = 1;
  pointMul(q, x, p);
  assert(!memcmp(p, q, SIZE_BYTES2));

  x[0] = 2;
  pointMul(t, x, p);
  pointDbl(s, p);
  assert(!memcmp(t, s, SIZE_BYTES2));

  polyRand(x);

  // Check if point multiplication works
  pointMul(q, x, p);
  assert(isValidPoint(q));

  pointMul_Naive(t, x, p);
  assert(isValidPoint(t));
  assert(!memcmp(q, t, SIZE_BYTES2));

  // Check Frobenius map
  pointSqr(q);                // q = p^2
  memcpy(s, q, SIZE_BYTES2);  // s = p^2
  pointSqr(q);                // q = p^4
  pointDbl(t, t);             // t = 2p
  pointAdd(q, q, s);          // q  = p^4 + p^2
  pointAdd(q, q, t);          // q  = p^4 + p^2 + 2p
  assert(isZero(q));

  polyRand(y);

  pointMul(q, x, p); // q =   x*p
  pointMul(q, y, q); // q = y*x*p
  pointMul(p, y, p); // p =   y*p
  pointMul(p, x, p); // p = x*y*p

  assert(!memcmp(p, q, SIZE_BYTES2));

  free(t);
  free(s);
  free(p);
  free(q);
  free(y);
  free(x);
}

void exchangetest() {
  word *x = polyGen(), *p = pointGen();

  assert(x && p);

  pointMul(p, polyRand(x), p);
  pointMul(p, x, p);

  free(p);
  free(x);
}

int main(int argc, char **argv) {
  double timing = 0;
  double d, delta = 0;
  unsigned int s = 0;
  unsigned int tests;
  srand(time(NULL));

  assert(argc == 2);
  tests = atoi(argv[1]);

  for (s = 0; s < tests; s++) {
    validtest();
  }

  for (s = 0; s < 5; s++) {
    timing = (double)clock();
    exchangetest();
    d = (clock() - timing) / CLOCKS_PER_SEC;
    delta += d;
    printf("Key Exchange in %0.3f seconds (delta %0.3f).\n", d, delta / s);
  }

  return 0;
}
