// Computes the magic number for unsigned division.
// adapted for 64 bit from Hackers Delight 
#include <stdint.h>
#include <intrin.h>  // _umul128
#include "Magic64.h"


struct mu magicu(uint64_t d) {
	// Must have 1 <= d <= 2**64-1.
	int p;
	uint64_t nc, delta, q1, r1, q2, r2;
	struct mu magu;

	magu.n = d;

	magu.a = 0;             // Initialize "add" indicator.
	//nc = -1 - (-d) % d;       // Unsigned arithmetic here.
	nc = 0xFFFFFFFFFFFFFFFF - (0xFFFFFFFFFFFFFFFF - d) % d;
	p = 63;                 // Init. p.
	q1 = 0x8000000000000000 / nc;     // Init. q1 = 2**p/nc.
	r1 = 0x8000000000000000 - q1*nc;// Init. r1 = rem(2**p, nc).
	q2 = 0x7FFFFFFFFFFFFFFF / d;      // Init. q2 = (2**p - 1)/d.
	r2 = 0x7FFFFFFFFFFFFFFF - q2*d; // Init. r2 = rem(2**p - 1, d).
	do {
		p = p + 1;
		if (r1 >= nc - r1) {
			q1 = 2 * q1 + 1;            // Update q1.
			r1 = 2 * r1 - nc;
		}          // Update r1.
		else {
			q1 = 2 * q1;
			r1 = 2 * r1;
		}
		if (r2 + 1 >= d - r2) {
			if (q2 >= 0x7FFFFFFFFFFFFFFF) magu.a = 1;
			q2 = 2 * q2 + 1;            // Update q2.
			r2 = 2 * r2 + 1 - d;
		}       // Update r2.
		else {
			if (q2 >= 0x8000000000000000) magu.a = 1;
			q2 = 2 * q2;
			r2 = 2 * r2 + 1;
		}
		delta = d - 1 - r2;
	} while (p < 128 &&
		(q1 < delta || (q1 == delta && r1 == 0)));

	magu.M = q2 + 1;        // Magic number
	magu.s = p - 64;        // and shift amount to return
	return magu;            // (magu.a was set above).
}

struct mu magicu2(uint64_t d) {

	// Must have 1 <= d <= 2**64-1.
	int p;
	uint64_t p64, q, r, delta;
	struct mu magu;
	magu.n = d;
	magu.a = 0;             // Initialize "add" indicator.
	p = 63;                 // Initialize p.
	q = 0x7FFFFFFFFFFFFFFF / d;       // Initialize q = (2**p - 1)/d.
	r = 0x7FFFFFFFFFFFFFFF - q*d;   // Init. r = rem(2**p - 1, d).
	do {
		p = p + 1;
		if (p == 64) p64 = 1;     // Set p32 = 2**(p-32).
		else p64 = 2 * p64;
		if (r + 1 >= d - r) {
			if (q >= 0x7FFFFFFFFFFFFFFF) magu.a = 1;
			q = 2 * q + 1;           // Update q.
			r = 2 * r + 1 - d;       // Update r.
		}
		else {
			if (q >= 0x8000000000000000) magu.a = 1;
			q = 2 * q;
			r = 2 * r + 1;
		}
		delta = d - 1 - r;
	} while (p < 128 && p64 < delta);
	magu.M = q + 1;         // Magic number and
	magu.s = p - 64;        // shift amount to return
	return magu;            // (magu.a was set above).
}

uint64_t magicdiv(uint64_t i, struct mu mag){
	uint64_t high;
	uint64_t low = _umul128(i, mag.M, &high);

	if (mag.a == 0)
		return high >> mag.s;
	else {
		// return (i + high) >> mag.s in overflow-safe way, 
		// using the technique described in http://ridiculousfish.com/blog/posts/labor-of-division-episode-i.html
		return  ( ((i - high) >> 1) + high ) >> (mag.s - 1);
	}
}

uint64_t magicmod(uint64_t i, struct mu mag){
	if (i >= mag.n)
		return i - mag.n*magicdiv(i, mag);
	else
		return i;
}
