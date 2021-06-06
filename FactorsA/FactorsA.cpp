#include <ctime>
#include <vector>
#include "../Factors/miller rabin.h"
#include "../Project1/mulmod.h"
#include <stdint.h>
#include "../Division by multiplication 64/Magic64.h"
#include <algorithm>  // sort, count
#include <unordered_set>
#include "../RuttenEekelen/RuttenEekelen.h"
#include "../Pollard Rho trials/Pollard Rho Montgomery.h"
#include "libdivide.h"
#include <random>

bool verbose = false;

// random number generator used in pollard rho initial values
std::linear_congruential_engine<
	uint64_t, 6364136223846793005U, 1442695040888963407U, 0U> mmix_rng(1234567ull);

void setVerbose(bool v) {
	verbose = v;
}

class Sieve {

public:
	std::vector<int> smallprimes;
	int NUMSMALLPRIMES;
	std::vector< libdivide::divider<uint64_t> > smallprimedividers;
	static const int PFSIZE = 1000L * 1000L;

	// primefactor[n] is 0 if n is prime, otherwise holds a prime factor of n.
	int primefactor[PFSIZE] = { 0 };


	void setfactor(int m) {
		if (primefactor[m] == 0)
			for (int j = m * m; j < PFSIZE; j += m) primefactor[j] = m;
	}

	int getprimefactor(int n) {
		return primefactor[n];
	}

	 Sieve() {
		setfactor(2);
		setfactor(3);
		for (int k = 5; k * k < PFSIZE; k += 6) {
			setfactor(k); setfactor(k + 2);
		}
		for (int p = 2; smallprimes.size() < 100; p++) {
			if (primefactor[p] == 0) {
				smallprimes.push_back(p);
				// see https://libdivide.com/documentation.html
				libdivide::divider<uint64_t> fast_d((uint64_t)p);
				smallprimedividers.push_back(fast_d);
			}
		}
		NUMSMALLPRIMES = 60;  // from trials, this seems an optimal value
	}

	 void setNumSmallPrimes(int i) {
		 if (i > smallprimes.size()) {
			 printf("ERROR attempting to set NUMSMALLPRIMES=%d, larger than smallprimes size = %zd\n", i, smallprimes.size());
			 exit(1);
		 }
		 NUMSMALLPRIMES = i;
	 }
};


// single instance of sieve to be used by all functions
Sieve sieve;

// for testing purposes
void setNumSmallPrimes(int i) {
	sieve.setNumSmallPrimes(i);
}


//factorise using pre-calculated array primefactor[], where primefactor[x] is a prime factor of x
void factorise_small(int n, uint64_t *primearray, int &pasize) {
	int i = sieve.primefactor[n];
	while (i != 0){
		primearray[pasize++]=(uint64_t)i;
		n /= i;
		i = sieve.primefactor[n];
	}
	primearray[pasize++] = (uint64_t)n;
}

uint64_t squareaddmodRE(uint64_t y, uint64_t  a, uint64_t n, RE &re) {
	y = mulmodRE(y, y, re);
	return addMod(y, a, n);
}


uint64_t squareaddmodMU(uint64_t y, uint64_t  a, uint64_t n, mu& magicn) {
	y = mulmodRR(y, y, magicn);
	return addMod(y, a, n);
}


uint64_t pollard_rhoMU(uint64_t n) {
	//pollard rho algorithm with brent variant, using magic numbers to divide by multiplication

	if (n % 2 == 0)
		return 2;

	struct mu magicn = magicu2(n);

	const uint64_t m = 200;
	uint64_t a, x, y, ys, r = 1, q = 1, g = 1;
	if (verbose)
		printf("\nCalling pollard_rhoMU for n=%llu", n);

	do
		a = mmix_rng() % n;
	while (a == 0 || a == n - 2);

	do
		y = mmix_rng() % n;
	while (y == 0 || y == n - 2);

	while (g == 1LL) {
		x = y;
		for (uint64_t i = 0LL; i < r; i++) {
			y = squareaddmodMU(y, a, n, magicn);
		}


		uint64_t k = 0;
		while (k < r && g == 1LL) {
			ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				y = squareaddmodMU(y, a,n, magicn);

				// q = q |x-y| mod n
				q = mulmodRR(q, (x > y) ? x - y : y - x, magicn);
			}
			g = gcd(q, n);
			k += m;
		}
		r *= 2uLL;
	}

	if (g == n) {
		do {
			// ys = ys² + a mod n
			ys = squareaddmodMU(ys, a, n, magicn);
			g = gcd((x > ys) ? x - ys : ys - x, n);
		} while (g == 1LL);
	}

	return g;
}


uint64_t pollard_rhoRE(uint64_t n) {
	// pollard rho algorithm with brent variation, using Rutten-Eekelin division

	if (n % 2 == 0)
		return 2;


	RE re = RE_gen(n);

	const uint64_t m = 200;
	uint64_t a, x, y, ys, r = 1, q = 1, g = 1;
	if (verbose)
		printf("\nCalling brent for n=%llu", n);

	do
		a = mmix_rng() % n;
	while (a == 0 || a == n - 2);

	do
		y = mmix_rng() % n;
	while (y == 0 || y == n - 2);

	while (g == 1LL) {
		x = y;
		for (uint64_t i = 0LL; i < r; i++) {
			y = squareaddmodRE(y, a, n, re);
		}


		uint64_t k = 0;
		while (k < r && g == 1LL) {
			ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				y = squareaddmodRE(y, a, n, re);

				// q = q |x-y| mod n
				q = mulmodRE(q, (x > y) ? x - y : y - x, re);
			}
			g = gcd(q, n);
			k += m;
		}
		r *= 2uLL;
	}

	if (g == n) {
		do {
			// ys = ys² + a mod n
			ys = squareaddmodRE(ys, a, n, re);
			g = gcd((x > ys) ? x - ys : ys - x, n);
		} while (g == 1LL);
	}

	return g;
}


void factorise_large(uint64_t n, uint64_t *primearray, int &pasize){
	//printf("calling factorise_large for n=%llu\n", n);
	uint64_t f; // might return a factor if is_prime returns false
	uint64_t p, q;

	if (is_prime(n, f)){
		primearray[pasize++] = (uint64_t)n;
		return;
	}
	if (verbose)
		printf("\nis_prime returned factor %llu for n=%llu", f, n);

	if (1 < f && f < n){
		if (verbose)
			printf("\nis_prime found factor %llu for n=%llu", f, n);
		p = f;
		q = n / p;
	}
	else {
		p = pollard_brent_montgomery_retry(n);
		q = n / p;
	}

	if (p > q){
		uint64_t t = p;
		p = q;
		q = t;
	}
	if (p == 1){
		// factorise_large going round again as it failed to find a factor
		factorise_large(q, primearray, pasize);
		return;
	}
	else{
		factorise_large(p, primearray, pasize);
		factorise_large(q, primearray, pasize);
		return;
	}

}

// to find out which method was used to factor a number
char method;

char getmethod(){
	return method;
}

void factorise(uint64_t n, uint64_t *primearray, int &pasize){

	//sieve();

	pasize = 0;

	if (n == 1LL){
		method = 'I';
		return;
	}
	
	if (n < sieve.PFSIZE){
		method = 'S';
		factorise_small(int(n), primearray, pasize);
		return;
	}
	
	// try trial division
	
	
	method = 'L';
	/*
	for (auto it = smallprimes.begin(); it != smallprimes.begin() + NUMSMALLPRIMES ; ++it) {
		uint64_t p = *it;
		if (p * p > n)
			break;
		while (n%p == 0){
			primearray[pasize++] = (uint64_t)p;
			n /= p;
		};
	};
	*/

	// trial division using libdivide
	// based on timing trials, libdivide halves the time spend in trial division	

	// handle factors of 2 using bit twiddling
	while ((n & 1) == 0) {
		primearray[pasize++] = 2ull;
		n >>= 1;
	}
	// trial division starting from 3, as we have already found factors of 2 above
	for (int i = 1; i < sieve.NUMSMALLPRIMES; ++i) {
		uint64_t p = sieve.smallprimes[i];
		if (n < sieve.PFSIZE)
			break;
		uint64_t q = n / sieve.smallprimedividers[i];
		while (q * (uint64_t)p == n) {
			primearray[pasize++] = (uint64_t)p;
			n = q;
			q = n / sieve.smallprimedividers[i];
		}
	}
	

	if (n == 1)
		return ;

	if (n < sieve.PFSIZE){
		factorise_small(int(n), primearray, pasize);
		method = 'P';
		return ;
	}

	factorise_large(n, primearray, pasize);
	method = 'Q';
	return ;
}

void allfactors(uint64_t n, uint64_t *factorsarray, int &factorssize){
	int pasize=0;
	uint64_t primearray[64];
	factorise(n, primearray, pasize);
	std::sort(primearray, primearray + pasize);
	factorssize = 1;
	factorsarray[0] = 1uLL;
	uint64_t lastprime = 0;
	uint64_t prev = 1;
	uint64_t multiplier = 1uLL;
	uint64_t p;
	for (int i = 0; i < pasize; i++){
		p = primearray[i];
		if (p == lastprime){
			multiplier *= p;
		}
		else{
			multiplier = p;
			prev = factorssize;
		}
		lastprime = p;
		for (int j = 0; j < prev; j++)
			factorsarray[factorssize++] = factorsarray[j] * multiplier;
	}
	std::sort(factorsarray, factorsarray + factorssize);
}

uint64_t totient(uint64_t n){
	int pasize = 0;
	uint64_t primearray[64];
	factorise(n, primearray, pasize);

	if (pasize == 1)
		return n - 1;

	std::unordered_set<uint64_t> distinctPrimes(primearray, primearray + pasize);
	uint64_t t = n;
	for (auto &p : distinctPrimes) 
		t -= t/p;
	return t;


}

uint64_t squarefreepart(uint64_t n){
	int pasize = 0;
	uint64_t primearray[64];
	factorise(n, primearray, pasize);
	std::sort(primearray, primearray + pasize);
	uint64_t sfp = 1ULL;
	std::unordered_set<uint64_t> distinctPrimes(primearray, primearray + pasize);
	for (auto &p : distinctPrimes) {
		uint64_t e = std::count(primearray, primearray + pasize, p);
		if (e%2){
			sfp *= p;
		}
	}
	return sfp;
}