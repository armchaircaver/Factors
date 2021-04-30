#include "miller rabin.h"
#include "mt64.h"
#include <cassert>
#include <iostream>
#include <ctime>
#include <stdint.h>
#include "FactorsA.h"
#include "FJ64_262K.h"
#include <RuttenEekelen.h>
#include <mulmod.h>

int countcomposites;
int factorcount;

void examine(uint64_t p, uint64_t mod) {
	uint64_t fac;
	uint64_t primearray[64];
	int pasize = 0;
	factorise(p, primearray, pasize);

	for (int it = 0; it != pasize; ++it) {
		if (primearray[it] < 300)
			p /= primearray[it];
	}

	if (p >= mod/10ULL){
		if (!is_prime(p, fac)){
			countcomposites++;
		}
		else
			assert(fac == 1ULL);


		if (fac > 1ULL){
			factorcount++;
			assert(p%fac == 0);
			printf("\t%llu, factors found=%llu,%llu\n", p, fac, p / fac);
		}
	}
}

// performance comparison of is_prime and   FJ64_262K is_prime_2_64
void performance_trial() {
	init_genrand64((uint64_t)time(NULL));
	uint64_t f;

	for (int index = 10; index < 64; index++) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);
		clock_t startTime = clock();
		int primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + count;
			if (is_prime(n, f))
				primecount++;
		}
		clock_t endTime = clock();
		printf("is_prime     : 2^%d\t %llu \t %8.2f\t", index, 1ull<<index,
			double(endTime - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);

		startTime = clock();
		primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + count;
			if (is_prime_2_64(n))
				primecount++;
		}
		printf("%8.2f\t", double(clock() - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);



		startTime = clock();
		primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + count;
			if (is_primeFJ(n,f))
				primecount++;
		}
		endTime = clock();
		printf("%8.2f\n", double(clock() - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);

	

	}
}

void verify() {
	init_genrand64((uint64_t)time(NULL));
	uint64_t f;

	for (int index = 60; index < 64; index++) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);
		int primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + genrand64_int64() % interval;
			bool n_is_prime = is_prime(n, f);
			if (n_is_prime != is_prime_2_64(n))
				printf("%llu, discrepancy is_prime_2_64\n", n); fflush(stdout);
			if (n_is_prime != is_primeFJ(n,f))
				printf("%llu, discrepancy is_primeFJ\n", n); fflush(stdout);
		}
		printf("10^%d verified\n", index);
	}
}

void longrun() {
	init_genrand64((uint64_t)time(NULL));
	uint64_t primecount = 0;
	for (int index = 63; index < 64; index++) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);
		uint64_t primecount = 0;
		for (int count = 0; count < 10000000000; count++) {
			uint64_t n = (1ULL << index) + genrand64_int64() % interval;
			if (is_prime_2_64(n)) {
				primecount++;
				if (primecount%10000==0)
					printf("%llu, primes\n", primecount); fflush(stdout);
			}
		}
		printf("10^%d verified\n", index);
	}
}

uint64_t pow_mod(uint64_t x, uint64_t y, RE re_mr) {
	uint64_t number = 1LL;

	while (y) {
		if (y & 1LL)
			number = mulmodRE(number, x,  re_mr);
		y >>= 1;
		x = mulmodRE(x, x,  re_mr);
	}
	return number;
}
bool miller_rabin_pass(uint64_t a, int s, uint64_t d, uint64_t n, uint64_t& factor, RE re_mr) {
	factor = 1ULL;
	if (n == a) return true;
	if (n % a == 0) return false;
	uint64_t x = pow_mod(a, d, re_mr);
	if ((x == 1ULL) || (x == n - 1LL)) return true;

	for (int i = 0; i < s; i++) {
		uint64_t lastx = x;
		x = mulmodRE(x, x, re_mr);
		if (x == n - 1ULL) return true;
		if (x == 1LL) {
			//factor = gcd(lastx - 1ULL, n);
			return false;
		}
	}
	return false;
}



void simple_sprp_compare() {
	for (auto w : { 2ull,3ull,5ull,7ull,101ull }) {
		printf("witness %llu\n", w);
		clock_t startTime = clock();
		int primes = 0;
		uint64_t start = (1ull << 62) + 1;
		for (uint64_t n = start; n < start + 1000000; n += 2) {
			if (is_SPRP(w, n))
				primes++;
		}
		clock_t endTime = clock();
		printf("is_SPRP     :  %8.2f, %d primes\n",
			double(endTime - startTime) / (double)CLOCKS_PER_SEC, primes); fflush(stdout);

		startTime = clock();
		struct RE re_mr;
		primes = 0;
		for (uint64_t n = start; n < start + 1000000; n += 2) {
			uint64_t factor = 1uLL;
			uint64_t d = n - 1LL;
			int s = 0;
			while (d % 2uLL == 0L) {
				d >>= 1;
				s++;
			}
			re_mr = RE_gen(n);
			if (miller_rabin_pass(w, s, d, n, factor, re_mr))
				primes++;
		}
		endTime = clock();
		printf("miller_rabin_pass:  %8.2f, %d primes\n",
			double(endTime - startTime) / (double)CLOCKS_PER_SEC, primes); fflush(stdout);

	}
}

void simple_is_prime_compare() {
	clock_t startTime = clock();
	for (int i = 40; i < 64; i += 2) {
		int primes = 0;
		printf("numbers starting 2**%d\n", i);
		uint64_t start = (1ull << i) + 1;
		for (uint64_t n = start; n < start + 10000000; n += 2) {

			if (n == 2 || n == 3 ||n == 5 || n == 7) {
				primes++; continue;
			}
			if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0) continue;
			if (n < 121) {
				if (n > 1) primes++;
				continue;
			}


			if (is_SPRP(2ull, n) && is_SPRP(bases[hashh(n)], n))
				primes++;
		}
		clock_t endTime = clock();
		printf("is_SPRP     :  %8.2f, %d primes\n",
			double(endTime - startTime) / (double)CLOCKS_PER_SEC, primes); fflush(stdout);


		startTime = clock();
		primes = 0;
		uint64_t factor = 1;
		for (uint64_t n = start; n < start + 10000000; n += 2) {
			if (is_primeFJ(n, factor))
				primes++;
		}
		endTime = clock();
		printf("is_primeFJ:  %8.2f, %d primes\n",
			double(endTime - startTime) / (double)CLOCKS_PER_SEC, primes); fflush(stdout);
	}

}


int main(){

	simple_is_prime_compare();
	simple_sprp_compare();
	performance_trial();
	verify();
	longrun();

	uint64_t seed = (uint64_t)time(NULL);

	uint64_t mod = 1000000LL;
	init_genrand64(seed);
	for (int index = 7; index<20; index++){
		countcomposites = 0;
		factorcount = 0;
		mod *= 10LL;
		clock_t startTime = clock();
		for (int i = 0; (i<100000) || (factorcount == 0); i++){
			uint64_t p = genrand64_int64() % mod;
			if (p >= mod / 10)
				examine(p,mod);
			}
		printf("10^%d: %d composites, %d factors, %8.2f sec\n", index, countcomposites, factorcount,
			   double(clock() - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);
	}
	return 0;
}