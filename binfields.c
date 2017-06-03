
#include "binfields.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Precomputed Array. For each byte (a1,a2,...,a8), this
   array contains the halfword (a1,0,a2,0,a3,0,...,0,a8).
   This is used for fast squaring. */

static const hword POW[1 << BYTESIZE] = {
    0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015, 0x0040,
    0x0041, 0x0044, 0x0045, 0x0050, 0x0051, 0x0054, 0x0055, 0x0100, 0x0101,
    0x0104, 0x0105, 0x0110, 0x0111, 0x0114, 0x0115, 0x0140, 0x0141, 0x0144,
    0x0145, 0x0150, 0x0151, 0x0154, 0x0155, 0x0400, 0x0401, 0x0404, 0x0405,
    0x0410, 0x0411, 0x0414, 0x0415, 0x0440, 0x0441, 0x0444, 0x0445, 0x0450,
    0x0451, 0x0454, 0x0455, 0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511,
    0x0514, 0x0515, 0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554,
    0x0555, 0x1000, 0x1001, 0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015,
    0x1040, 0x1041, 0x1044, 0x1045, 0x1050, 0x1051, 0x1054, 0x1055, 0x1100,
    0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115, 0x1140, 0x1141,
    0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155, 0x1400, 0x1401, 0x1404,
    0x1405, 0x1410, 0x1411, 0x1414, 0x1415, 0x1440, 0x1441, 0x1444, 0x1445,
    0x1450, 0x1451, 0x1454, 0x1455, 0x1500, 0x1501, 0x1504, 0x1505, 0x1510,
    0x1511, 0x1514, 0x1515, 0x1540, 0x1541, 0x1544, 0x1545, 0x1550, 0x1551,
    0x1554, 0x1555, 0x4000, 0x4001, 0x4004, 0x4005, 0x4010, 0x4011, 0x4014,
    0x4015, 0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054, 0x4055,
    0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115, 0x4140,
    0x4141, 0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155, 0x4400, 0x4401,
    0x4404, 0x4405, 0x4410, 0x4411, 0x4414, 0x4415, 0x4440, 0x4441, 0x4444,
    0x4445, 0x4450, 0x4451, 0x4454, 0x4455, 0x4500, 0x4501, 0x4504, 0x4505,
    0x4510, 0x4511, 0x4514, 0x4515, 0x4540, 0x4541, 0x4544, 0x4545, 0x4550,
    0x4551, 0x4554, 0x4555, 0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011,
    0x5014, 0x5015, 0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054,
    0x5055, 0x5100, 0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115,
    0x5140, 0x5141, 0x5144, 0x5145, 0x5150, 0x5151, 0x5154, 0x5155, 0x5400,
    0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414, 0x5415, 0x5440, 0x5441,
    0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455, 0x5500, 0x5501, 0x5504,
    0x5505, 0x5510, 0x5511, 0x5514, 0x5515, 0x5540, 0x5541, 0x5544, 0x5545,
    0x5550, 0x5551, 0x5554, 0x5555};

word *__reduce(word *c) {
  int i;
  word T;
#if 0
    if (deg(&c[SIZE_WORDS]) || deg(c)>SIZE_BITS) {
#endif
  for (i = SIZE_WORDS2 - 1; i >= SIZE_WORDS; i--) {
    T = c[i];
    c[i - SIZE_WORDS] ^= (T << 5) ^ (T << 7) ^ (T << 10) ^ (T << 15);
    c[i - SIZE_WORDS + 1] ^= (T >> (WORDSIZE -  5)) ^ (T >> (WORDSIZE -  7)) ^
                             (T >> (WORDSIZE - 10)) ^ (T >> (WORDSIZE - 15));
  }
  T = c[SIZE_WORDS - 1] >> (WORDSIZE - 5);
  c[0] = c[0] ^ T ^ (T << 2) ^ (T << 5) ^ (T << 10);
  c[SIZE_WORDS - 1] = c[SIZE_WORDS - 1] & EMPTY_MASK;
#if 0
    }
#endif
  return c;
}

word *polyMul(word *r, const word *a, const word *b) {
  int k, j;
  word c[SIZE_WORDS2] = {0};

  for (k = WORDSIZE - 1; k; k--) {
    for (j = 0; j < SIZE_WORDS; j++)
      if ((a[j] >> k) & 1) polyAddTo(&c[j], b);
    iLShiftN(c, SIZE_WORDS2, 1);
    // TODO: possibly speed this up?
  }

  for (j = 0; j < SIZE_WORDS; j++)
    if (a[j] & 1) polyAddTo(&c[j], b);

  return memcpy(r, __reduce(c), SIZE_BYTES);
}

word *polySqr(word *r, const word *a) {
  word i, c[SIZE_WORDS2];
  for (i = 0; i < SIZE_WORDS; i++) {
    word t = a[i], x;
    int j;
#define LOBYTE(__t) (__t & ((1 << BYTESIZE) - 1))
#define LOOPS (WORDSIZE / HWSIZE)
    for (x = 0, j = 0; j < LOOPS; j++) {
      x += (word)POW[LOBYTE(t)] << (j * HWSIZE);
      t >>= BYTESIZE;
    }
    c[2 * i] = x;
    for (x = 0, j = 0; j < LOOPS; j++) {
      x += (word)POW[LOBYTE(t)] << (j * HWSIZE);
      t >>= BYTESIZE;
    }
    c[2 * i + 1] = x;
#undef LOOPS
#undef LOBYTE
  }
  return memcpy(r, __reduce(c), SIZE_BYTES);
}

word *polyDiv(word *c, const word *a, const word *b) {
  word t[SIZE_WORDS];
  return polyMul(c, a, polyInv(t, b));
}

word *polyInv(word *r, const word *a) {
  word *t;
  int j;

  word x[5 * SIZE_WORDS] = {0}, *v = &x[1 * SIZE_WORDS],
             *u = &x[2 * SIZE_WORDS], *g = &x[3 * SIZE_WORDS],
             *f = &x[4 * SIZE_WORDS];

  memcpy(u, a, SIZE_BYTES);

  v[0] = 0x425;
  v[SIZE_WORDS - 1] = F571;
  g[0] = 1;

inv_loop:
  if (u[0] == 1 || u[0] == 0) {
    for (j = 1; j < SIZE_WORDS; j++)
      if (u[j]) goto inv_run;
    goto inv_done;
  }
inv_run:
  if ((j = deg(u) - deg(v)) < 0) {
    t = v;
    v = u;
    u = t; /* v <-> u */
    t = g;
    g = f;
    f = t; /* g <-> f */
    j = -j;
  }

  memcpy(x, v, SIZE_BYTES);
  iLShiftN(x, SIZE_WORDS, j);
  polyAddTo(u, x); /* u = u + v>>j */

  memcpy(x, f, SIZE_BYTES);
  iLShiftN(x, SIZE_WORDS, j);
  polyAddTo(g, x); /* g = g + f>>j */

  goto inv_loop;
inv_done:

  memcpy(r, g, SIZE_BYTES);
  return r;
}

word *polyAddTo(word *a, const word *b) {
  unsigned int i;
  for (i = 0; i < SIZE_WORDS; i++) a[i] ^= b[i];
  return a;
}

word *polyAdd(word *c, const word *a, const word *b) {
  unsigned int i;
  for (i = 0; i < SIZE_WORDS; i++) c[i] = a[i] ^ b[i];
  return c;
}

void iLShiftN(word *a, word s, word n) {
  word x, *p = a + s - 1;
  if ((x = n / WORDSIZE) != 0) {
    memmove(a + x, a, (s - x) * sizeof(word));
    memset(a, 0, x * sizeof(word));
    a += x;
  }
  if (n %= WORDSIZE) {
    for (x = WORDSIZE - n; p > a; p--) {
      *p <<= n;
      *p |= *(p - 1) >> x;
    }
    *p <<= n;
  }
}

void iRShiftN(word *a, word s, word n) {
  word x, *p = a + s - 1;
  if ((x = n / WORDSIZE) != 0) {
    memmove(a, a + x, (s - x) * sizeof(word));
    memset(a + x, 0, x * sizeof(word));
    p -= x;
  }
  if (n %= WORDSIZE) {
    for (x = (((word)1) << n) - 1; a < p; a++) {
      *a >>= n;
      *a |= *(a + 1) & x;
    }
    *a >>= n;
  }
}

word *polyGen() { return calloc(SIZE_WORDS, sizeof(word)); }

word *polyRand(word *p) {
  int i, j;
  for (i = SIZE_WORDS - 1; i >= 0; i--) {
    for (j = 0, p[i] = 1; j < WORDSIZE / BYTESIZE; j++) p[i] *= rand();
  }
  p[SIZE_WORDS - 1] &= EMPTY_MASK;
  return p;
}

int deg(const word *a) {
  word h;
  int d = (SIZE_WORDS - 1) * WORDSIZE;
  for (a += SIZE_WORDS - 1; d && !*a; a--) d -= WORDSIZE;
  for (h = *a; h; h >>= 1) d++;
  return d;
}

int polyIsZero(word *a) {
  int cntr;
  for (cntr = 0; cntr < SIZE_WORDS; cntr++)
    if (a[cntr]) return 0;
  return 1;
}

word *dup(word *a) { return memcpy(polyGen(), a, SIZE_BYTES); }
