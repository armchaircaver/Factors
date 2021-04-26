#include <ctime>
#include <vector>

#include "miller rabin.h"
#include "ipow.h"
#include "mulmod.h"
#include <stdint.h>

const int PFSIZE = 1000L * 1000L;
int primefactor[PFSIZE] = { 0 };
int MAXSMALLPRIME = 1000000;
std::vector<int> smallprimes;

void setMaxSmallPrime(int i){
	MAXSMALLPRIME = i;
}
void setfactor(int m){
	if (primefactor[m]==0)
		for (int j = m*m; j<PFSIZE; j += m) primefactor[j] = m;
}

void sieve(){
	setfactor(2);
	setfactor(3);
	for (int k = 5; k*k<PFSIZE; k += 6) {
		setfactor(k); setfactor(k + 2);
	}
	for (int p = 2; p<MAXSMALLPRIME;p++) {
		if (primefactor[p] == 0)
			smallprimes.push_back(p);
	}
	srand(int(time(NULL)));
}

std::vector<uint64_t> factorise_small(int n) {
	int i = primefactor[n];
	std::vector<uint64_t> f;
	while (i != 0){
		f.push_back(i);
		n /= i;
		i = primefactor[n];
	}
	f.push_back(n);
	return f;
}

const uint64_t max_uint64_t = 18446744073709551615ULL;

uint64_t random_uint64_t(){
	uint64_t out = 0ULL;
	for (int i = 0; i < 8; i++)
		out |= ((uint64_t)(rand() & 255)) << (i * 8);
	return out;
}


uint64_t squareaddmod(uint64_t y, uint64_t  a, uint64_t n){
	y = mulmod4(y, y, n);
	y += a;
	if (y < a)
		y += (max_uint64_t - n) + 1;
	return y %= n;
}
uint64_t brent(uint64_t n){
	const uint64_t m = 1000;
	uint64_t a, x, y, ys, r=1, q=1, g=1;
	if (n % 2 == 0)
		return 2;

	do
		a = random_uint64_t() % n;
	while (a == 0 || a == n - 2);

	do
		y = random_uint64_t() % n;
	while (y == 0 || y == n - 2);

	//printf("\nBrent,n=%llu,  a=%llu, y=%llu",n, a, y);
	r = 1LL;
	q = 1LL;

	while (g == 1LL) {
		x = y;
		for (uint64_t i = 0LL; i < r; i++) 
			y = squareaddmod(y, a, n);
		
		uint64_t k = 0;
		while (k < r && g == 1LL) {
			ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				y = squareaddmod(y, a, n);

				// q = q |x-y| mod n
		 	    q = mulmod(q, (x>y) ? x - y : y - x, n);
			}
			g = gcd(q, n);
			k += m;
		} ;

		r *= 2LL;
	};

	if (g == n) {
		do {
			// ys = ys² + a mod n
			ys = squareaddmod(ys, a, n);
			g = gcd((x>ys) ? x - ys : ys - x, n);
		} while (g == 1LL);
	}

	return g;
}


std::vector<uint64_t> factorise_large(uint64_t n){
	std::vector<uint64_t> factors;
	uint64_t f; // might return a factor if is_prime returns false
	uint64_t p, q;

	if (is_prime(n, f)){
		factors.push_back(n);
		return factors;
	}
	if (1 < f && f < n){
		p = f;
		q = n / p;
	}
	else {
		p = brent(n);
		q = n / p;
	}
	
	if (p > q){
		uint64_t t = p;
		p = q;
		q = t;
	}
	if (p == 1){
		//printf("\n factorise_large going round again, p=%llu, q=%llu, n=%llu", p, q, n);
		return factorise_large(q);
	}
	else{
		std::vector<uint64_t> fp = factorise_large(p);
		std::vector<uint64_t> fq = factorise_large(q);
		factors.reserve(factors.size() +fp.size() + fq.size());
		factors.insert(factors.end(), fp.begin(), fp.end());
		factors.insert(factors.end(), fq.begin(), fq.end());
		fp.clear();
		fq.clear();
		return factors;
	}

}

std::vector<uint64_t> factorise(uint64_t n){

	std::vector<uint64_t> factors;
	if (n == 1)
		return factors;
	if (n < PFSIZE)
		return factorise_small(int(n));
	//printf("\n too large for factorise_small");
	int p = 1;
	for (auto pp = smallprimes.begin(); pp != smallprimes.end() && p*p<=n && p<=MAXSMALLPRIME; ++pp) {
		p = *pp;
		while (n%p == 0){
			factors.push_back(p);
			n /= p;
		};
	};
	//printf("\n extracted small primes");

	if (n == 1)
		return factors;

	uint64_t f;
	if (is_prime(n,f)){
		//printf("\n remaining quotient=%llu is prime", n);
		factors.push_back(n);
		return factors;
	}

	if (n < PFSIZE){
		//printf("\n remaining quotient=%llu, now calling factorise_small",n);
		std::vector<uint64_t> f2 = factorise_small(int(n));
		factors.reserve(factors.size() + f2.size());
		factors.insert(factors.end(), f2.begin(), f2.end());
		f2.clear();
		return factors;
	}

	std::vector<uint64_t> f2 = factorise_large(n);
	factors.reserve(factors.size() + f2.size());
	factors.insert(factors.end(), f2.begin(), f2.end());
	return factors;

}
