#include <stdint.h>
#include <stdio.h>
#include <vector>
#include "../RuttenEekelen/RuttenEekelen.h"
#include "../miller rabin behaviour/FJ64_262K.h"



uint64_t gcd(uint64_t a, uint64_t b){
	uint64_t c;
	while (a != 0LL) {
		c = a; a = b%a;  b = c;
	}
	return b;
}


uint64_t pow_mod(uint64_t x, uint64_t y, uint64_t n, RE & re_mr){
	uint64_t number = 1LL;

	while (y){
		if (y & 1LL)
			number = mulmodRE(number, x, re_mr);
		y >>= 1;
		x = mulmodRE(x, x, re_mr);
	}
	return number;
}

bool miller_rabin_pass(uint64_t a, int s, uint64_t d, uint64_t n, uint64_t &factor, RE & re_mr){
	factor = 1ULL;
	if (n == a) return true;
	if (n%a == 0) return false;
	uint64_t x = pow_mod(a, d, n, re_mr);
	if ((x == 1ULL) || (x == n - 1LL)) return true;

	for (int i = 0; i<s; i++){
		uint64_t lastx = x;
		x = mulmodRE(x, x, re_mr);
		if (x == n - 1ULL) return true;
		if (x == 1LL) {
			factor = gcd(lastx-1ULL, n);
			return false;
		}
	}
	return false;
}

// old version kept for reference and testing newer versions
bool is_prime_ref(uint64_t n, uint64_t &factor){
	factor = 1uLL;
	if (n == 2uLL) return true;
	if ((!(n & 1uLL)) || (n < 2uLL)) return false;

	uint64_t d = n - 1LL;
	int s = 0;
	while (d % 2uLL == 0L){
		d >>= 1;
		s++;
	}


	struct RE re_mr = RE_gen(n);

	std::vector<uint64_t> witnesses = {};

	if (n < 2047LL) 
		witnesses = { 2 };
	else if (n<1373653LL) 
		 witnesses = { 2,3 };
	else if (n < 9080191LL) 
		 witnesses = { 31,73 }; 
	else if (n < 4759123141LL)  
		 witnesses = { 2, 7, 61 };
	else if (n < 1122004669633LL)  
		 witnesses = { 2, 13, 23, 1662803 };
	else if (n < 2152302898747LL)
		 witnesses = { 2, 3, 5, 7, 11 };
	else if (n < 3474749660383LL)
		witnesses = { 2, 3, 5, 7, 11, 13 };
	else if (n < 341550071728321LL)
		witnesses = { 2, 3, 5, 7, 11, 13, 17 };
	else if (n < 3825123056546413051LL)
		witnesses = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };
	else // from unverified note in https://oeis.org/A014233, 
		// this is supposed to work up to 2^64
		witnesses = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

	for (auto w = witnesses.begin(); w != witnesses.end(); ++w) 
		if (!miller_rabin_pass(*w, s, d, n, factor, re_mr))
			return false;

	return true; // passed tests for all witnesses
}


bool is_primeFJ(uint64_t n, uint64_t& factor) {
	factor = 1uLL;

	if (n == 2 || n == 3 || n == 5 || n == 7) return true;
	if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0) return false;
	if (n < 121) return (n > 1);

	if (n == 2uLL) return true;
	if ((!(n & 1uLL)) || (n < 2uLL)) return false;

	uint64_t d = n - 1LL;
	int s = 0;
	while (d % 2uLL == 0L) {
		d >>= 1;
		s++;
	}


	struct RE re_mr = RE_gen(n);

	return miller_rabin_pass(2, s, d, n, factor, re_mr) && miller_rabin_pass(bases[hashh(n)], s, d, n, factor, re_mr);
	
}

bool is_prime(uint64_t n, uint64_t& factor) {
	if (n < (1ull<<43) )
		return is_primeFJ(n, factor);

	// For larger numbers, use the montgomery mutiplication technique.
	// We don't find a potential factor using this technique
	factor = 1ull;
	return is_prime_2_64(n);
}
