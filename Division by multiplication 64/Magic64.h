#pragma once

#include <stdint.h>

struct mu {
	uint64_t M;     // Magic number,
	int a;          // "add" indicator,
	int s;			// and shift amount.
	uint64_t n;		// original number for reference
};

struct mu magicu(uint64_t d);

struct mu magicu2(uint64_t d);

uint64_t magicdiv(uint64_t a, struct mu mag);

uint64_t magicmod(uint64_t i, struct mu mag);

