#ifndef MONTGOMERY
#define MONTGOMERY

#include <stdint.h>

struct monty {
	uint64_t m;		// original number for reference
	uint64_t rinv;
	uint64_t mprime;
	uint64_t hr;
};

struct monty prepareMonty(uint64_t m);

uint64_t montymulmod(uint64_t a, uint64_t b, struct monty M);

#endif
