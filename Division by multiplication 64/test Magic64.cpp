#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <intrin.h>
#include "Magic64.h"
#include <assert.h>


int main() {
	//struct mu magicu(uint64_t);
	//struct mu magicu2(uint64_t);
	struct mu mag;

	mag = magicu(3);
	printf("magicu(3), M = %llx, a = %d, s = %d\n", mag.M, mag.a, mag.s);
	if (mag.M + mag.a + mag.s != 0xAAAAAAAAAAAAAAAB + 0 + 1) return 1;
	mag = magicu(7);
	printf("magicu(7), M = %llx, a = %d, s = %d\n", mag.M, mag.a, mag.s);
	if (mag.M + mag.a + mag.s != 0x2492492492492493 + 1 + 3) return 2;

	mag = magicu2(3);
	printf("magicu2(3),  M = %llx a = %d s = %d\n", mag.M, mag.a, mag.s);
	if (mag.M + mag.a + mag.s != 0xAAAAAAAAAAAAAAAB + 0 + 1) return 11;
	mag = magicu2(7);
	printf("magicu2(7), M = %llx a = %d s = %d\n", mag.M, mag.a, mag.s);
	if (mag.M + mag.a + mag.s != 0x2492492492492493 + 1 + 3) return 12;

	if (0x7fffffffffffffff < 1ULL << 63){
		printf("0x7fffffffffffffff < 1ULL << 63\n");
		printf("%016llX\n", (0x7fffffffffffffff<<1) + 1);
	}

	uint64_t m = 1;
	for (int i = 1; i <= 63;i++){
		 m = (m<<1)+1 ;
		printf("m=%016llX hex  ", m,m);
		if (m == 1) exit(0);
		struct mu mag = magicu(m);
		if (mag.M > 0x7fffffffffffffff)
			printf("Bigger M=% 016llX, s=%d , a=%d\n", mag.M, mag.s, mag.a);
		else
			printf("Smaller M=%016llX, s=%d , a=%d\n", mag.M, mag.s, mag.a);
		
		for (uint64_t i = 10 * m; i <= 10 * m + 10000; i++){
			if (i / m != magicdiv(i, mag)){
				printf("M=%llu, s=%d , a=%d\n", mag.M, mag.s, mag.a);
				printf("magicdiv failed for i=%llu, m=%llu, i/m=%llu, magicdiv=%llu\n", i, m, i / m, magicdiv(i, mag));
				exit(1);
			}
		}
	}

	printf("All succeeded\n");

	int counter = 0;
	for (int i = 0; i < 20; i++){
		for (uint64_t m = 1 << i; m < (1 << i) + 100000; m++){
			struct mu mag = magicu(m);
			struct mu mag2 = magicu2(m);
			assert(mag.M == mag2.M && mag.s == mag2.s && mag.a == mag2.a);
			//printf("m=%llu, mag M=% 016llX, s=%d , a=%d\n", m, mag.M, mag.s, mag.a);
			//("m=%llu, mag2 M=% 016llX, s=%d , a=%d\n", m, mag2.M, mag2.s, mag.a);
			counter++;
		}
	}
	printf("magicu and magicu2 match for all %d cases\n",counter);

	return 0;
}
