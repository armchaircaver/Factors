#ifndef MONTGOMERY
#define MONTGOMERY

#include <stdint.h>

struct monty_t {
	uint64_t m;		// original number for reference
	uint64_t rinv;
	uint64_t mprime;
	uint64_t hr;
};

monty_t prepareMonty(uint64_t m);

uint64_t montmul(uint64_t abar, uint64_t bbar, monty_t &M);

uint64_t modul64(uint64_t x, uint64_t y, uint64_t z);

uint64_t montymulmod(uint64_t a, uint64_t b, struct monty_t M);

uint64_t reverse(uint64_t p, monty_t M);

#endif
