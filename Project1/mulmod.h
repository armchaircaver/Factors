#include <stdint.h>
#include "../Division by multiplication 64/Magic64.h"

uint64_t mulmodSA(uint64_t a, uint64_t b, uint64_t m);

uint64_t mulmodRR(uint64_t a, uint64_t b, struct mu mag);

int msb(uint64_t value);

uint64_t mulmodAS(uint64_t a, uint64_t b, uint64_t m);

uint64_t mulmod(uint64_t a, uint64_t  b, uint64_t n, struct mu magicn, struct RE re);

uint64_t addMod(uint64_t a, uint64_t b, uint64_t m);

#define MULMODMETHOD 4