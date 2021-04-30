#pragma once

#include <stdint.h>

struct RE {
	uint64_t Mdash;     // Magic number,
	int p;          
	int t;			
	uint64_t sdashmask;
	uint64_t M=0;		// original number for reference
};

int ceillog2(uint64_t n);

int floorlog2(uint64_t n);

struct RE RE_gen(uint64_t M);

uint64_t mulmodRE(uint64_t a, uint64_t b, struct RE re);

