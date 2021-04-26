#include <stdint.h>
#include <intrin.h>  // _umul128, _BitScanReverse64
#include "RuttenEekelen.h"
#include <algorithm> // std::swap
#include "mulmod.h"
#include <cassert>


int ceillog2(uint64_t n){ // ceiling of log base 2
	unsigned long index;

	_BitScanReverse64(&index, n);
	if (n&(n - 1)) // if n is not a power of 2
		return ++index;
	else
		return index;
}

int floorlog2(uint64_t n){ // floor of log base 2
	unsigned long index;
	_BitScanReverse64(&index, n);
	return index;
}

// stolen from Ridiculous Fish libdivide
int32_t libdivide__count_leading_zeros64(uint64_t val) {
	unsigned long result;
	if (_BitScanReverse64(&result, val)) {
		return 63 - result;
	}
	return 0;
}


/* Code taken from Hacker's Delight, http://www.hackersdelight.org/HDcode/divlu.c .  
   License permits inclusion here per http://www.hackersdelight.org/permissions.htm
*/
uint64_t libdivide_128_div_64_to_64(uint64_t u1, uint64_t u0, uint64_t v, uint64_t *r) {
	const uint64_t b = (1ULL << 32); // Number base (16 bits).
	uint64_t un1, un0,        // Norm. dividend LSD's.
		vn1, vn0,        // Norm. divisor digits.
		q1, q0,          // Quotient digits.
		un64, un21, un10,// Dividend digit pairs.
		rhat;            // A remainder.
	int s;                  // Shift amount for norm.

	if (u1 >= v) {            // If overflow, set rem.
		if (r != NULL)         // to an impossible value,
			*r = (uint64_t)(-1);    // and return the largest
		return (uint64_t)(-1);
	}    // possible quotient.

	/* count leading zeros */
	s = libdivide__count_leading_zeros64(v); // 0 <= s <= 63.
	if (s > 0) {
		v = v << s;           // Normalize divisor.
		un64 = (u1 << s) | ((u0 >> (64 - s)) & (-s >> 31));
		un10 = u0 << s;       // Shift dividend left.
	}
	else {
		// Avoid undefined behavior.
		un64 = u1 | u0;
		un10 = u0;
	}

	vn1 = v >> 32;            // Break divisor up into
	vn0 = v & 0xFFFFFFFF;     // two 32-bit digits.

	un1 = un10 >> 32;         // Break right half of
	un0 = un10 & 0xFFFFFFFF;  // dividend into two digits.

	q1 = un64 / vn1;            // Compute the first
	rhat = un64 - q1*vn1;     // quotient digit, q1.
again1:
	if (q1 >= b || q1*vn0 > b*rhat + un1) {
		q1 = q1 - 1;
		rhat = rhat + vn1;
		if (rhat < b) goto again1;
	}

	un21 = un64*b + un1 - q1*v;  // Multiply and subtract.

	q0 = un21 / vn1;            // Compute the second
	rhat = un21 - q0*vn1;     // quotient digit, q0.
again2:
	if (q0 >= b || q0*vn0 > b*rhat + un0) {
		q0 = q0 - 1;
		rhat = rhat + vn1;
		if (rhat < b) goto again2;
	}

	if (r != NULL)            // If remainder is wanted,
		*r = (un21*b + un0 - q0*v) >> s;     // return it.
	return q1*b + q0;
}


struct RE RE_gen(uint64_t d){

	struct RE re;

	re.p = ceillog2(d);
	re.sdashmask = UINT64_MAX << re.p;
	re.t = 64 - re.p - (int)(1/d);
	re.M = d;

	if ((d & (d - 1)) == 0){
		re.Mdash = 0;
		return re;
	}
	
	uint64_t rem;
	uint64_t proposed_m = libdivide_128_div_64_to_64(1ULL << floorlog2(d), 0, d, &rem); //== (1 << (64 + floor_log_2_d)) / d

	// LIBDIVIDE_ASSERT(rem > 0 && rem < d);
	uint64_t e = d - rem;

	/* We have to use the general 65-bit algorithm.
		We need to compute(2 * *power) / d.However, we already have(2 * *(power - 1)) / d and its remainder.
		By doubling both, and then correcting the remainder, we can compute the larger division. * /
	*/ 
	proposed_m += proposed_m; //don't care about overflow here - in fact, we expect it
	uint64_t twice_rem = rem + rem;
	if (twice_rem >= d || twice_rem < rem)
		proposed_m++;

	re.Mdash = proposed_m ;
	return re;
}

// noinline used for profiling
/*__declspec(noinline) */ uint64_t mulmodRE(uint64_t a, uint64_t b, struct RE re){

	if (re.M & 0x8000000000000000)
		return mulmodAS(a, b, re.M); //RE cannot process numbers above 2^63

	int p = re.p;
	int t = re.t;
	uint64_t M = re.M;
	uint64_t Mdash = re.Mdash;

	uint64_t v;
	uint64_t u = _umul128(a, b, &v); // v,u are the high and low words of a*b

	uint64_t	s = u >> p;

	// combine h= v<<t and hdash=h+s ( h+s == h|s as s is low (64-p) bits and h is high p bits_
	uint64_t	hdash = (v << t) | s;

	// combine q=hdash*Mdash and qdash=q+hdash
	uint64_t	qdash;
	uint64_t	low = _umul128(hdash, Mdash, &qdash);
	qdash += hdash;

	uint64_t	y = qdash*M;

	// combine sdash = s<<p and d=sdash-y
	uint64_t d = (s << p) - y;

	uint64_t	r = u - y;

	uint64_t rdash;
	uint64_t rdashdash;
	if (d < M)
		rdash = r;
	else
		rdash = (r - M);

	if (rdash < M)
		rdashdash = rdash;
	else
		rdashdash = (rdash - M);

	if (rdashdash >= M)
	rdashdash -= M;

	/*
	if (rdashdash != mulmodSA(a, b, M)){
		printf("MISMATCH for a=%llu,b=%llu,M=%llu\n", a, b, M);
		printf("rdashdash=%llu, mulmodSA=%llu\n", rdashdash, mulmodSA(a, b, M));
		exit(1);
	}
	*/
	return rdashdash;
}


