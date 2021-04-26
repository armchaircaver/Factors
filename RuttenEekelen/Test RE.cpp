#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <intrin.h>
#include "RuttenEekelen.h"
#include <assert.h>

int main() {

	for (uint64_t i = 1; i < 10; i++){
		printf("i = %llu, ceillog2(i) = %d, floorlog2(i) = %d\n", i, ceillog2(i), floorlog2(i));
	}

	struct RE re;

	for (uint64_t n = 1; n < 10; n++){
		re = RE_gen(n);
		printf("RE(%llu), Mdash = %llu, (%llx hex), p = %d, t = %d\n", n, re.Mdash, re.Mdash, re.p, re.t);
	}

	struct RE re10 = RE_gen(10);
	for (uint64_t n = 1; n < 10; n ++){
		uint64_t mod = mulmodRE(11ULL, n, re10);
		printf("n=%llu, 11*n=%llu, mod=%llu\n",n, 11 * n, mod);
	}

	for (uint64_t m = 10Ull; m <= 10000000000000000000Ull; m *= 10){
		struct RE re = RE_gen(m);
		uint64_t total = 0Ull;
		for (uint64_t a = 1000; a < 100000000; a++){
			total += mulmodRE(a, a, re);
		}
		printf("m=%llu, total=%llu\n", m, total);
	}
}