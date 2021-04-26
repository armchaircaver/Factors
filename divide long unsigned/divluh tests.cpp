#include <stdio.h>
#include <stdlib.h>     // To define "exit", req'd by XLC.
#include <stdint.h>
#include "divulh.h"

int errors;
void error(unsigned x, unsigned y, unsigned z, unsigned r) {
	errors = errors + 1;
	printf("Error for x = %08x, y = %08x, z = %08x,"
		" got %08x (%d dec)\n", x, y, z, r, r);
}

int errors64;
void error64(uint64_t x, uint64_t y, uint64_t z, uint64_t r, uint64_t rem) {
	errors64 = errors64 + 1;
	printf("Error64 for x = %llu, y = %llu, z = %llu, got %llu, rem=%llu \n", x, y, z, r, rem);
}

int main() {
	int i, r, n;
	uint64_t r64;

	static unsigned test[] = { 0, 0, 1, 0, 0, 1, 1, 1, 0, 21, 7, 3,
		0, 0xffffffff, 1, 0xffffffff, 0, 0xffffffff, 10, 0x19999999,
		1, 1, 2, 0x80000000, 1, 1, 3, 0x55555555, 2, 4, 4, 0x80000001,
		0x32123456, 0x789abcde, 0x45600000, 0xb8c45446,
		0x80000000, 0x00000000, 0x80000001, 0xfffffffe,
		0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff,
		0x40000000, 0x00000000, 0xc0000000, 0x55555555,
	};

	n = sizeof(test) / 4;

	printf("divluh:\n");
	for (i = 0; i < n; i += 4) {
		r = divluh(test[i], test[i + 1], test[i + 2]);
		if (r != test[i + 3]) error(test[i], test[i + 1], test[i + 2], r);
	}

	printf("divluh1:\n");
	for (i = 0; i < n; i += 4) {
		r = divluh1(test[i], test[i + 1], test[i + 2]);
		if (r != test[i + 3]) error(test[i], test[i + 1], test[i + 2], r);
	}

	if (errors == 0)
		printf("Passed all %d cases.\n", n / 4);

	static uint64_t test64[] = {
		0, 0, 1, 0, 0,
		0, 1, 1, 1, 0,
		0, 21, 7, 3, 0,
		0, 0xffffffff, 1, 0xffffffff, 0,
		0, 0xffffffff, 10, 0x19999999, 5,
		1, 1, 2, 0x8000000000000000, 1,
		1, 1, 3, 0x5555555555555555, 2,
		2, 4, 4, 0x8000000000000001, 0,
		12345678, 87654321, 10000000000009, 22773756248222, 1426056935971,
		123456787654321, 87654321012345678, 10000000000000009, 227737576602156248, 791856907872182,
		2, 0, 3, 12297829382473034410, 2,
	};

	n = sizeof(test64) / 8;
	printf("%d 64 bit tests:\n", n / 5);

	printf("divluh64:\n");
	uint64_t rem;
	for (i = 0; i < n; i += 5) {
		//printf("64 bit test %d\n", i);
		r64 = divluh64(test64[i], test64[i + 1], test64[i + 2], rem);
		if (r64 != test64[i + 3] || rem != test64[i + 4]) error64(test64[i], test64[i + 1], test64[i + 2], r64, rem);
	}

	if (errors64 == 0)
		printf("64 bit passed all %d cases.\n", n / 5);

}
