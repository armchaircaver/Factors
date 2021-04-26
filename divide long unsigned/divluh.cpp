// Divide long unsigned by hardware "restoring" method.
#include <stdio.h>
#include <stdlib.h>     // To define "exit", req'd by XLC.
#include <stdint.h>

unsigned divluh(unsigned x, unsigned y, unsigned z) {
	// Divides (x || y) by z.
	int i;
	unsigned t;

	for (i = 1; i <= 32; i++) {
		t = (int)x >> 31;         // All 1's if x(31) = 1.
		x = (x << 1) | (y >> 31); // Shift x || y left
		y = y << 1;               // one bit.
		if ((x | t) >= z) {
			x = x - z;
			y = y + 1;
		}
	}
	return y;                    // Remainder is x.
}

uint64_t divluh64(uint64_t x, uint64_t y, uint64_t z, uint64_t &rem) {
	// Divides (x || y) by z.
	uint64_t t;

	for (int i = 1; i <= 64; i++) {
		t = (int64_t)x >> 63;         // All 1's if x(63) = 1.
		x = (x << 1) | (y >> 63); // Shift x || y left
		y <<= 1;               // one bit.
		if ((x | t) >= z) {
			x -= z;
			y++;
		}
	}
	rem = x;
	return y;                    // Remainder is x.
}

// Version based on nonrestoring hardware algorithm:
// This more-or-less verifies the algorithm given in HD,
// but we don't include this C version because it's not very efficient.
unsigned divluh1(unsigned x, unsigned y, unsigned z) {
	// Divides (x || y) by z.
	int i;
	unsigned c;

	c = 0;
	for (i = 1; i <= 32; i++) {
		if (c == 0) {
			c = x >> 31;
			x = (x << 1) | (y >> 31); // Shift x || y left
			y = y << 1;               // one bit.
			c = c ^ (x < z);
			x = x - z;
		}
		else {
			c = x >> 31;
			x = (x << 1) | (y >> 31); // Shift x || y left
			y = y << 1;               // one bit.
			x = x + z;
			c = c ^ (x < z);
		}
		y = y + (1 - c);
	}
	return y;                    // Remainder is x.
}

