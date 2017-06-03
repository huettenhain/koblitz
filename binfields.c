#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "binfields.h"

unsigned int __log2(word x);

#if defined(_MSC_VER) && (_MSC_VER >= 1300)
# include <intrin.h>
# pragma message("Using BitScanReverse Intrinsic for MSC.")
# ifdef BITS64
#  pragma intrinsic(_BitScanReverse64)
#  define BSR _BitScanReverse64
# else
#  pragma intrinsic(_BitScanReverse)
#  define BSR _BitScanReverse
# endif
  inline unsigned int __log2(word x) {
    word log2;
    BSR( &log2, (word)x );
    return log2;
  }
#elif defined(__GNUC__) && ((__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
  inline unsigned int __log2(word x) {
    return WORDSIZE-1-__builtin_clzll(x);
  }
#else 
// DeBruijn method: https://stackoverflow.com/q/11376288
//  thanks to Desmond Hume
# ifdef BITS64
  unsigned int __log2(word x) {
    static const unsigned char _bitpos[64] = {
      63,  0, 58,  1, 59, 47, 53,  2,
      60, 39, 48, 27, 54, 33, 42,  3,
      61, 51, 37, 40, 49, 18, 28, 20,
      55, 30, 34, 11, 43, 14, 22,  4,
      62, 57, 46, 52, 38, 26, 32, 41,
      50, 36, 17, 19, 29, 10, 13, 21,
      56, 45, 25, 31, 35, 16,  9, 12,
      44, 24, 15,  8, 23,  7,  6,  5 };
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return _bitpos[(((x - (x >> 1))*0x07EDD5E59A4E28C2ULL)) >> 58];
  }
# else 
  unsigned int __log2(word x) {
    static const unsigned char _bitpos[32] = {
      0,  9,  1, 10, 13, 21, 2, 29,
     11, 14, 16, 18, 22, 25, 3, 30,
      8, 12, 20, 28, 15, 17, 24, 7, 
     19, 27, 23,  6, 26,  5,  4, 31 };
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return _bitpos[(x * 0x07C4ACDDU) >> 27];
  }
# endif 
#endif

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
  for (i = SIZE_WORDS2 - 1; i >= SIZE_WORDS; i--) {
    T = c[i];
    c[i - SIZE_WORDS] ^= (T << 5) ^ (T << 7) ^ (T << 10) ^ (T << 15);
    c[i - SIZE_WORDS + 1] ^= (T >> (WORDSIZE -  5)) ^ (T >> (WORDSIZE -  7)) ^
                             (T >> (WORDSIZE - 10)) ^ (T >> (WORDSIZE - 15));
  }
  T = c[SIZE_WORDS - 1] >> (WORDSIZE - 5);
  c[0] = c[0] ^ T ^ (T << 2) ^ (T << 5) ^ (T << 10);
  c[SIZE_WORDS - 1] = c[SIZE_WORDS - 1] & EMPTY_MASK;
  return c;
}


inline static void __shiftLeft(word *c) {
  for (word *p=c+SIZE_WORDS2-1; p > c; p--)
      *p = ((*p) << 1) | (*(p-1) >> (WORDSIZE-1));
  *c <<= 1;
}
word *polyMul(word *r, const word *a, const word *b) {
  int k, j;
  word c[SIZE_WORDS2] = {0};

  for (k = WORDSIZE - 1; k; k--) {
    for (j = 0; j < SIZE_WORDS; j++)
      if ((a[j] >> k) & 1) polyAddTo(&c[j], b);
    __shiftLeft(c);
  }

  for (j = 0; j < SIZE_WORDS; j++)
    if (a[j] & 1) polyAddTo(&c[j], b);

  return memcpy(r, __reduce(c), SIZE_BYTES);
}

inline static word __low_byte_pow(word t) {
  return (word) POW[ t & ((1 << BYTESIZE) - 1) ];
}
word *polySqr(word *r, const word *a) {
  word i, c[SIZE_WORDS2];
  for (i = 0; i < SIZE_WORDS; i++) {
    word t = a[i], x;
    int j;
    for (x = 0, j = 0; j < HW_PER_W; j++) {
      x += __low_byte_pow(t) << (j * HWSIZE);
      t >>= BYTESIZE;
    }
    c[2 * i] = x;
    for (x = 0, j = 0; j < HW_PER_W; j++) {
      x += __low_byte_pow(t) << (j * HWSIZE);
      t >>= BYTESIZE;
    }
    c[2 * i + 1] = x;
  }
  return memcpy(r, __reduce(c), SIZE_BYTES);
}

word *polyDiv(word *c, const word *a, const word *b) {
  word t[SIZE_WORDS];
  return polyMul(c, a, polyInv(t, b));
}

word *polyInv(word *r, const word *a) {
  word *t;
  unsigned int j, du, dv, dt;

  word x[5 * SIZE_WORDS] = {0}, 
      *v = &x[1 * SIZE_WORDS],
      *u = &x[2 * SIZE_WORDS], 
      *g = &x[3 * SIZE_WORDS],
      *f = &x[4 * SIZE_WORDS];

  memcpy(u, a, SIZE_BYTES);
  du = polyDeg(u);

  v[0] = 0x425;
  v[SIZE_WORDS - 1] = F571;
  dv = polyDeg(v);

  g[0] = 1;

inv_loop:
  if (u[0] == 1 || u[0] == 0) {
    for (j = 1; j < SIZE_WORDS; j++)
      if (u[j]) goto inv_run;
    goto inv_done;
  }
inv_run:
  if (du < dv) {
    t = v; v = u; u = t; /* v <-> u */
    t = g; g = f; f = t; /* g <-> f */
    dt=du; du=dv; dv=dt;
  }
  j = du - dv; 
  memcpy(x, v, SIZE_BYTES);
  polyLShift(x, j);
  polyAddTo(u, x); /* u = u + v>>j */
  
  for (dt=du/WORDSIZE; dt && !u[dt]; dt--)
    ;
  du = dt*WORDSIZE + __log2(u[dt]) + 1;
# ifdef _DEBUG 
  assert( du == polyDeg(u) );
# endif 

  memcpy(x, f, SIZE_BYTES);
  polyLShift(x, j);
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

void polyLShift(word *a, unsigned long n) {
  word x, *p = a + SIZE_WORDS - 1;
  if ((x = n / WORDSIZE) != 0) {
    memmove(a + x, a, (SIZE_WORDS - x) * sizeof(word));
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

word *polyGen() { return calloc(SIZE_WORDS, sizeof(word)); }

word *polyRand(word *p) {
  int i, j;
  for (i = SIZE_WORDS - 1; i >= 0; i--) {
    for (j = 0, p[i] = 1; j < WORDSIZE / BYTESIZE; j++) p[i] *= rand();
  }
  p[SIZE_WORDS - 1] &= EMPTY_MASK;
  return p;
}

unsigned int polyDeg(const word *a) {
  unsigned int d;
  for (d=SIZE_WORDS-1; d && !a[d]; d--) 
    ;
#ifdef _DEBUG
  unsigned int  log2 = 0;
  unsigned int _log2 = __log2(a[d]);
  for (word h = a[d]; h; h >>= 1) log2++;  
  log2--;
  assert(log2==_log2);
#endif
  return (d*WORDSIZE) + __log2(a[d]) + 1;
}

int polyIsZero(const word *a) {
  int cntr;
  for (cntr = 0; cntr < SIZE_WORDS; cntr++)
    if (a[cntr]) return 0;
  return 1;
}

word *dup(word *a) { return memcpy(polyGen(), a, SIZE_BYTES); }
