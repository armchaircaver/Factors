/* Computes a*b mod m using Montgomery multiplication (MM). a, b, and m
are unsigned numbers with a, b < m < 2**64, and m odd. The code does
some 128-bit arithmetic.
The variable r is fixed at 2**64, and its log base 2 at 64.
Works with gcc on Windows and Linux. */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "montgomery.h"
#include <intrin0.h>
#include <immintrin.h>

typedef unsigned long long uint64;
typedef long long int64;

/* ---------------------------- mulul64 ----------------------------- */

/* Multiply unsigned long 64-bit routine, i.e., 64 * 64 ==> 128.
Parameters u and v are multiplied and the 128-bit product is placed in
(*whi, *wlo). It is Knuth's Algorithm M from [Knu2] section 4.3.1.
Derived from muldwu.c in the Hacker's Delight collection. */
// superseded by intrinsic _umul128

void mulul64(uint64 u, uint64 v, uint64 *whi, uint64 *wlo)
{
	*wlo = _umul128(u, v, &*whi);
	return;	
}

/* ---------------------------- modul64 ----------------------------- */

uint64 modul64(uint64 x, uint64 y, uint64 z) {

	/* Divides (x || y) by z, for 64-bit integers x, y,
	and z, giving the remainder (modulus) as the result.
	Must have x < z (to get a 64-bit result). This is
	checked for. */


	//printf("In modul64, x = %016llx, y = %016llx, z = %016llx\n", x, y, z);
	if (x >= z) {
		printf("Bad call to modul64, must have x < z.");
		exit(1);
	}
	uint64_t r;
	_udiv128(x, y, z, &r); // don't need the quotient function result
	return r;	
}

/* ---------------------------- montmul ----------------------------- */

uint64 montmul(uint64 abar, uint64 bbar, monty_t & M) {

	uint64 thi, tlo, tm, tmmhi, tmmlo, uhi, ulo, ov;

	mulul64(abar, bbar, &thi, &tlo);  // t = abar*bbar.

	/* Now compute u = (t + ((t*mprime) & mask)*m) >> 64.
	The mask is fixed at 2**64-1. Because it is a 64-bit
	quantity, it suffices to compute the low-order 64
	bits of t*mprime, which means we can ignore thi. */

	tm = tlo*M.mprime;

	mulul64(tm, M.m, &tmmhi, &tmmlo);   // tmm = tm*m.

	ulo = tlo + tmmlo;                // Add t to tmm
	uhi = thi + tmmhi;                // (128-bit add).
	if (ulo < tlo) uhi = uhi + 1;     // Allow for a carry.

	// The above addition can overflow. Detect that here.

	ov = (uhi < thi) | ((uhi == thi) & (ulo < tlo));

	ulo = uhi;                   // Shift u right
	uhi = 0;                     // 64 bit positions.

	// if (ov > 0 || ulo >= m)      // If u >= m,
	//    ulo = ulo - M.m;            // subtract m from u.
	ulo = ulo - (M.m & ( (UINT64_MAX - (ov | (ulo >= M.m)))+1)  ); // Alternative
	// with no branching.

	if (ulo >= M.m){
		printf("ERROR in montmul, ulo = %016llx, m = %016llx\n", ulo, M.m);
		exit(1);
	}
	return ulo;
}

/* ---------------------------- xbinGCD ----------------------------- */

/* C program implementing the extended binary GCD algorithm. C.f.
http://www.ucl.ac.uk/~ucahcjm/combopt/ext_gcd_python_programs.pdf. This
is a modification of that routine in that we find s and t s.t.
gcd(a, b) = s*a - t*b,
rather than the same expression except with a + sign.
This routine has been greatly simplified to take advantage of the
facts that in the MM use, argument a is a power of 2, and b is odd. Thus
there are no common powers of 2 to eliminate in the beginning. The
parent routine has two loops. The first drives down argument a until it
is 1, modifying u and v in the process. The second loop modifies s and
t, but because a = 1 on entry to the second loop, it can be easily seen
that the second loop doesn't alter u or v. Hence the result we want is u
and v from the end of the first loop, and we can delete the second loop.
The intermediate and final results are always > 0, so there is no
trouble with negative quantities. Must have a either 0 or a power of 2
<= 2**63. A value of 0 for a is treated as 2**64. b can be any 64-bit
value.
Parameter a is half what it "should" be. In other words, this function
does not find u and v st. u*a - v*b = 1, but rather u*(2a) - v*b = 1. */

void xbinGCD(uint64 a, uint64 b, uint64 *pu, uint64 *pv)
{
	uint64 alpha, beta, u, v;
	//printf("Doing GCD(%llx hex, %llx hex)\n", a, b);

	u = 1; v = 0;
	alpha = a; beta = b;         // Note that alpha is even and beta is odd.

	/* The invariant maintained from here on is:
	a = u*2*alpha - v*beta. */

	// printf("Before, a u v = %016llx %016llx %016llx\n", a, u, v);
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
		//    printf("After,  a u v = %016llx %016llx %016llx\n", a, u, v);
	}

	// printf("At end,    a u v = %016llx %016llx %016llx\n", a, u, v);
	*pu = u;
	*pv = v;
	return;
}

/* ------------------------------ main ------------------------------ */

monty_t prepareMonty(uint64_t m){
	 monty_t M;

	M.m = m;
	M.hr = 0x8000000000000000LL;
	xbinGCD(M.hr, m, &M.rinv, &M.mprime);      // xbinGCD, in effect, doubles hr.
	if (2 * M.hr*M.rinv - m*M.mprime != 1) {
		printf("The Extended Euclidean algorithm failed.");
		printf("m=%llu, hr=%llu, rinv=%llu, mprime=%llu, test=%llu\n", m, M.hr, M.rinv, M.mprime, 2 * M.hr*M.rinv - m*M.mprime);
		exit(1);
	}
	return M;
}

//convert montgomery space number p to an ordinary number
uint64_t reverse(uint64_t p, monty_t M) {
	uint64 phi, plo;
	mulul64(p, M.rinv, &phi, &plo);
	return modul64(phi, plo, M.m);
}

uint64_t montymulmod(uint64_t a, uint64_t b, monty_t M){
	uint64_t abar = modul64(a, 0, M.m);
	uint64_t bbar = modul64(b, 0, M.m);

	uint64_t p = montmul(abar, bbar, M); /* Compute a*b (mod m). */

	/* Convert p back to a normal number by p = (p*rinv)%m. */
	return reverse(p, M);
}


#if 0  // ************ Main program removed **************
int main(int argc, char* argv[]) {
	char *q;
	uint64 a, b, m, hr, rinv, mprime, p1hi, p1lo, p1, p, abar, bbar;
	uint64 phi, plo;

	if (argc != 4) {
		printf("Need exactly three arguments, decimal numbers.");
		return 1;
	}

	a = strtoull(argv[1], &q, 0);
	if (*q != 0 || a == 0xFFFFFFFFFFFFFFFFLL) {
		printf("Invalid first argument.");
		return 1;
	}

	b = strtoull(argv[2], &q, 0);
	if (*q != 0 || b == 0xFFFFFFFFFFFFFFFFLL) {
		printf("Invalid second argument.");
		return 1;
	}

	m = strtoull(argv[3], &q, 0);         // The modulus, we are computing a*b (mod m).
	if (*q != 0) {
		printf("Invalid third argument.");
		return 1;
	}

	if ((m & 1) == 0) {
		printf("The modulus (third argument) must be odd.");
		return 1;
	}

	if (a >= m || b >= m) {
		printf("The first two args must be less than the modulus (third argument).");
		return 1;
	}

	printf("a, b, m = %016llx %016llx %016llx\n", a, b, m);

	/* The simple calculation: This computes (a*b)**4 (mod m) correctly for all a,
	b, m < 2**64. */

	mulul64(a, b, &p1hi, &p1lo);         // Compute a*b (mod m).
	p1 = modul64(p1hi, p1lo, m);
	mulul64(p1, p1, &p1hi, &p1lo);       // Compute (a*b)**2 (mod m).
	p1 = modul64(p1hi, p1lo, m);
	mulul64(p1, p1, &p1hi, &p1lo);       // Compute (a*b)**4 (mod m).
	p1 = modul64(p1hi, p1lo, m);
	printf("p1 = %016llx\n", p1);

	/* The MM method uses a quantity r that is the smallest power of 2
	that is larger than m, and hence also larger than a and b. Here we
	deal with a variable hr that is just half of r. This is because r can
	be as large as 2**64, which doesn't fit in one 64-bit word. So we
	deal with hr, where 2**63 <= hr <= 1, and make the appropriate
	adjustments wherever it is used.
	We fix r at 2**64, and its log base 2 at 64. It doesn't hurt if
	they are too big, it's just that some quantities (e.g., mprime) come
	out larger than they would otherwise be. */

	hr = 0x8000000000000000LL;

	/* Now, for the MM method, first compute the quantities that are
	functions of only r and m, and hence are relatively constant. These
	quantities can be used repeatedly, without change, when raising a
	number to a large power modulo m.
	First use the extended GCD algorithm to compute two numbers rinv
	and mprime, such that

	r*rinv - m*mprime = 1

	Reading this nodulo m, clearly r*rinv = 1 (mod m), i.e., rinv is the
	multiplicative inverse of r modulo m. It is needed to convert the
	result of MM back to a normal number. The other calculated number,
	mprime, is used in the MM algorithm. */

	xbinGCD(hr, m, &rinv, &mprime);      // xbinGCD, in effect, doubles hr.

	/* Do a partial check of the results. It is partial because the
	multiplications here give only the low-order half (64 bits) of the
	products. */

	printf("rinv = %016llx, mprime = %016llx\n", rinv, mprime);
	if (2 * hr*rinv - m*mprime != 1) {
		printf("The Extended Euclidean algorithm failed.");
		return 1;
	}

	/* Compute abar = a*r(mod m) and bbar = b*r(mod m). That is, abar =
	(a << 64)%m, and bbar = (b << 64)%m. */

	abar = modul64(a, 0, m);
	bbar = modul64(b, 0, m);

	p = montmul(abar, bbar, m, mprime); /* Compute a*b (mod m). */
	p = montmul(p, p, m, mprime);       /* Compute (a*b)**2 (mod m). */
	p = montmul(p, p, m, mprime);       /* Compute (a*b)**4 (mod m). */
	printf("p before converting back = %016llx\n", p);

	/* Convert p back to a normal number by p = (p*rinv)%m. */

	mulul64(p, rinv, &phi, &plo);
	p = modul64(phi, plo, m);
	printf("p = %016llx\n", p);
	if (p != p1) printf("ERROR, p != p1.\n");
	else         printf("Correct (p = p1).\n");

	return 0;
}
#endif
