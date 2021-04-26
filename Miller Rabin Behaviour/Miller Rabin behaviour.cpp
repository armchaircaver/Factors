#include "miller rabin.h"
#include "mt64.h"
#include <cassert>
#include <iostream>
#include <ctime>
#include <stdint.h>
#include "FactorsA.h"
#include "FJ64_262K.h"

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
			uint64_t n = (1ULL << index) + genrand64_int64() % interval;
			if (is_prime(n, f))
				primecount++;
		}
		clock_t endTime = clock();
		printf("is_prime     : 10^%d\t %d \t %8.2f\t", index, primecount,
			double(endTime - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);

		startTime = clock();
		primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + genrand64_int64() % interval;
			if (is_prime_2_64(n))
				primecount++;
		}
		printf("%8.2f\t", double(clock() - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);



		startTime = clock();
		primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + genrand64_int64() % interval;
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
int main(){

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