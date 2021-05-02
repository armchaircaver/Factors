#include <ctime>
#include <vector>
#include "miller rabin.h"
#include "mulmod.h"
#include <stdint.h>
#include "Magic64.h"
#include <algorithm>  // sort, count
#include <unordered_set>
#include "RuttenEekelen.h"
#include "../Pollard Rho trials/Pollard Rho Montgomery.h"


const uint64_t PFSIZE = 1000L * 1000L;
int primefactor[PFSIZE] = { 0 };
int MAXSMALLPRIME = 300;
std::vector<int> smallprimes;
bool verbose = false;
bool sieved = false;

void setMaxSmallPrime(int i){
	MAXSMALLPRIME = i;
}

void setVerbose(bool v){
	verbose = v;
}
void setfactor(int m){
	if (primefactor[m] == 0)
		for (int j = m*m; j<PFSIZE; j += m) primefactor[j] = m;
}

void sieve(){
	if (sieved)
		return;
	setfactor(2);
	setfactor(3);
	for (int k = 5; k*k<PFSIZE; k += 6) {
		setfactor(k); setfactor(k + 2);
	}
	for (int p = 2; p<MAXSMALLPRIME; p++) {
		if (primefactor[p] == 0)
			smallprimes.push_back(p);
	}
	srand(int(time(NULL)));
	sieved = true;
}

void factorise_small(int n, uint64_t *primearray, int &pasize) {
	int i = primefactor[n];
	while (i != 0){
		primearray[pasize++]=(uint64_t)i;
		n /= i;
		i = primefactor[n];
	}
	primearray[pasize++] = (uint64_t)n;
}

uint64_t random_uint64_t(){
	uint64_t out = 0ULL;
	for (int i = 0; i < 8; i++)
		out |= ((uint64_t)(rand() & 255)) << (i * 8);
	return out;
}

//struct mu magicn;
//struct RE re;

/*
uint64_t squareaddmod(uint64_t y, uint64_t  a, uint64_t n){
	y = mulmod(y, y, n, magicn, re);
	y = addMod(y, a, n);
	return y;
}
*/

uint64_t squareaddmodRE(uint64_t y, uint64_t  a, uint64_t n, RE &re) {
	y = mulmodRE(y, y, re);
	return addMod(y, a, n);
}
/*
uint64_t brent(uint64_t n) {

	if (n % 2 == 0)
		return 2;


#if MULMODMETHOD == 2
	magicn = magicu2(n);
#elif  MULMODMETHOD == 4
	re = RE_gen(n);
#endif

	const uint64_t m = 200;
	uint64_t a, x, y, ys, r = 1, q = 1, g = 1;
	if (verbose)
		printf("\nCalling brent for n=%llu", n);

	do
		a = random_uint64_t() % n;
	while (a == 0 || a == n - 2);

	do
		y = random_uint64_t() % n;
	while (y == 0 || y == n - 2);

	while (g == 1LL) {
		x = y;
		for (uint64_t i = 0LL; i < r; i++) {
			y = squareaddmod(y, a, n);
		}


		uint64_t k = 0;
		while (k < r && g == 1LL) {
			ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				y = squareaddmod(y, a, n);

				// q = q |x-y| mod n
				q = mulmod(q, (x > y) ? x - y : y - x, n, magicn, re);
			}
			g = gcd(q, n);
			k += m;
		}
		r *= 2uLL;
	}

	if (g == n) {
		do {
			// ys = ys² + a mod n
			ys = squareaddmod(ys, a, n);
			g = gcd((x > ys) ? x - ys : ys - x, n);
		} while (g == 1LL);
	}

	return g;
}
*/

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
		a = random_uint64_t() % n;
	while (a == 0 || a == n - 2);

	do
		y = random_uint64_t() % n;
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
		a = random_uint64_t() % n;
	while (a == 0 || a == n - 2);

	do
		y = random_uint64_t() % n;
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

	if (!sieved)
		sieve();

	pasize = 0;

	if (n == 1LL){
		method = 'I';
		return;
	}
	
	if (n < PFSIZE){
		method = 'S';
		factorise_small(int(n), primearray, pasize);
		return;
	}
	
	method = 'L';
	int p = 1;
	for (auto pp = smallprimes.begin(); pp != smallprimes.end() && p*p <= n && p <= MAXSMALLPRIME; ++pp) {
		p = *pp;
		while (n%p == 0){
			primearray[pasize++] = (uint64_t)p;
			n /= p;
		};
	};

	if (n == 1)
		return ;

	if (n < PFSIZE){
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