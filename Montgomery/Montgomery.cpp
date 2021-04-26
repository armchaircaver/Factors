#include <stdio.h>
#include <cstdlib>

typedef unsigned long long uint64;
typedef long long int64;

void xbinGCD(uint64 a, uint64 b, uint64 *pu, uint64 *pv)
{
	uint64 alpha, beta, u, v;

	u = 1; v = 0;
	alpha = a; beta = b;         // Note that alpha is
	// even and beta is odd.

	/* The invariant maintained from here on is:
	a = u*2*alpha - v*beta. */

	while (a > 0) {
		a = a >> 1;
		if ((u & 1) == 0) {             // Delete a common
			u = u >> 1; v = v >> 1;      // factor of 2 in 
		}                               // u and v. 
		else {
			/* We want to set u = (u + beta) >> 1, but
			that can overflow, so we use Dietz's method. */
			u = ((u ^ beta) >> 1) + (u & beta);
			v = (v >> 1) + alpha;
		}
	}

	*pu = u;
	*pv = v;
	return;
}

uint64 modul64(uint64 x, uint64 y, uint64 z) {

	/* Divides (x || y) by z, for 64-bit integers x, y,
	and z, giving the remainder (modulus) as the result.
	Must have x < z (to get a 64-bit result). This is
	checked for. */

	int64 i, t;

	if (x >= z) {
		printf("Bad call to modul64, must have x < z.");
		exit(1);
	}
	for (i = 1; i <= 64; i++) {  // Do 64 times.
		t = (int64)x >> 63;       // All 1's if x(63) = 1.
		x = (x << 1) | (y >> 63); // Shift x || y left
		y = y << 1;               // one bit.
		if ((x | t) >= z) {
			x = x - z;
			y = y + 1;
		}
	}
	return x;                    // Quotient is y.
}

void mulul64(uint64 u, uint64 v, uint64 *whi, uint64 *wlo)
{
	uint64 u0, u1, v0, v1, k, t;
	uint64 w0, w1, w2;

	u1 = u >> 32; u0 = u & 0xFFFFFFFF;
	v1 = v >> 32; v0 = v & 0xFFFFFFFF;

	t = u0*v0;
	w0 = t & 0xFFFFFFFF;
	k = t >> 32;

	t = u1*v0 + k;
	w1 = t & 0xFFFFFFFF;
	w2 = t >> 32;

	t = u0*v1 + w1;
	k = t >> 32;

	*wlo = (t << 32) + w0;
	*whi = u1*v1 + w2 + k;

	return;
}

uint64 montmul(uint64 abar, uint64 bbar, uint64 m,
	uint64 mprime) {

	uint64 thi, tlo, tm, tmmhi, tmmlo, uhi, ulo, ov;

	mulul64(abar, bbar, &thi, &tlo);  // t = abar*bbar.
	/* Now compute u = (t + ((t*mprime) & mask)*m) >> 64.
	The mask is fixed at 2**64-1. Because it is a 64-bit
	quantity, it suffices to compute the low-order 64
	bits of t*mprime, which means we can ignore thi. */

	tm = tlo*mprime;

	mulul64(tm, m, &tmmhi, &tmmlo);   // tmm = tm*m.

	ulo = tlo + tmmlo;                // Add t to tmm
	uhi = thi + tmmhi;                // (128-bit add).
	if (ulo < tlo) uhi = uhi + 1;     // Allow for a carry.

	// The above addition can overflow. Detect that here.

	ov = (uhi < thi) | ((uhi == thi) & (ulo < tlo));

	ulo = uhi;                   // Shift u right 
	uhi = 0;                     // 64 bit positions. 

	if (ov > 0 || ulo >= m)      // If u >= m,
		ulo = ulo - m;            // subtract m from u.

	return ulo;
}

/*
uint64 monty
mulul64(p, rinv, &phi, &plo);
p = modul64(phi, plo, m);
*/