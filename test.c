
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "curve.h"
#include "binfields.h"

int validtest() {
    int rv;
    word *x = polyRand(polyGen()),
         *y = polyRand(polyGen()),
         *q = pointGen(),
         *p = pointGen();

    assert(isValidPoint(q));

    pointMul(q,x,p);
    pointMul(q,y,q);
    pointMul(p,y,p);
    pointMul(p,x,p);

    rv = memcmp(p,q,SIZE_BYTES2);

    free(p); free(q);
    free(y); free(x);

    return !rv;
}

void exchangetest() {

    word *x = polyGen(),
         *p = pointGen();

    assert(x && p);

    pointMul(p,polyRand(x),p);
    pointMul(p,x,p);

    free(p);
    free(x);
}


int main(int argc, char **argv) {
	
    double timing = 0;
    double d,delta = 0;
    unsigned int s=0;
    unsigned int tests;
    srand(time(NULL));

    assert(argc == 2);
    tests = atoi(argv[1]);

    while (++s <= tests) {
        if (!validtest()) {
            printf("A Validity Test Failed.\n");
            getchar();
        }
        timing = (double) clock();
        exchangetest();
        d = (clock()-timing)/CLOCKS_PER_SEC;
        delta += d;
        printf("Key Exchange in %0.3f seconds (delta %0.3f).\n",
           d,delta/s);
    }

	return 0;
}
