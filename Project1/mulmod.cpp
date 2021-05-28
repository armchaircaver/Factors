#include <intrin.h>  // _BitScanReverse64
#include <stdint.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <algorithm> // std::swap
#include "Magic64.h"
#include "RuttenEekelen.h"
#include "mulmod.h"


int msb(uint64_t value){ //return the most significant bit set
	unsigned long index;

	_BitScanReverse64(&index, value);
	return ++index;
}

uint64_t addMod(uint64_t a, uint64_t b, uint64_t m) {
	return (a < m - b) ? (a + b) : (a - (m - b));
}

// deals with 64 bit m
uint64_t mulmodSA(uint64_t a, uint64_t b, uint64_t m) {
	uint64_t res = 0;
	uint64_t temp_b;

	/* Only needed if b may be >= m */
	if (b >= m) {
		if (m > UINT64_MAX / 2u)
			b -= m;
		else
			b %= m;
	}

	if (a >= m) {
		if (m > UINT64_MAX / 2u)
			a -= m;
		else
			a %= m;
	}

	if (a > b)std::swap(a, b);

	while (a != 0) {
		if (a & 1) {
			/* Add b to res, modulo m, without overflow */
			if (b >= m - res) /* Equiv to if (res + b >= m), without overflow */
				res -= m;
			res += b;
		}
		a >>= 1;

		/* Double b, modulo m */
		temp_b = b;
		if (b >= m - b)       /* Equiv to if (2 * b >= m), without overflow */
			temp_b -= m;
		b += temp_b;
	}
	return res;
}

// assembler version
//extern "C" uint64_t  _mulmod_64(uint64_t a, uint64_t b, uint64_t m);

uint64_t mulmodAS(uint64_t a, uint64_t b, uint64_t m){
	//a = a >= m ? a %= m : a;
	//b = b >= m ? b %= m : b;
	//return _mulmod_64(a, b, m);
	//replace assembler version with  version using intrinsics, now that _udiv128 is available in VS19
	uint64_t high;
	uint64_t low = _umul128(a, b, &high);
	uint64_t rem;
	uint64_t quotient = _udiv128(high, low, m, &rem);
	return rem;


}
// repeated reduction of most significant 64 bits

uint64_t mulmodRR(uint64_t a, uint64_t b, struct mu mag){
	uint64_t n = mag.n;

	//if (n < 0x100000000)
		// return(((a%n)*(b%n) % n));
	//	return magicmod(magicmod(a, mag)*magicmod(b, mag), mag);

	if (n &  0x8000000000000000ULL) // mulmodRR doesn't work for n >=2^63
		return mulmodAS(a,b,n);

	uint64_t high;
	uint64_t low = _umul128(a, b, &high);
	while (high){
		int msbhi = msb(high);
		high <<= (64 - msbhi);
		high |= low >> msbhi;  // slightly faster than __shiftright128
		uint64_t highmod = magicmod(high, mag);
		low = (low & ((1ULL << msbhi) - 1ULL)) | (highmod << msbhi);
		high = highmod >> (64 - msbhi);
	}
	return magicmod(low, mag);
}



uint64_t mulmod(uint64_t a, uint64_t  b, uint64_t n, struct mu magicn, struct RE re){
#if MULMODMETHOD == 1
	return mulmodSA(a, b, n);
#elif MULMODMETHOD == 2
	return mulmodRR(a, b, magicn);
#elif MULMODMETHOD == 3
	return mulmodAS(a, b, n);
#elif  MULMODMETHOD == 4
	return mulmodRE(a, b, re);
#endif
}