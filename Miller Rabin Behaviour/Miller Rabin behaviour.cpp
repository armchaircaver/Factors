#include <cassert>
#include <iostream>
#include <ctime>
#include <stdint.h>
#include <random>

#include "../FactorsA/FactorsA.h"
#include "../mpzPollard/Timer.h"
#include "../Factors/miller rabin.h"
#include "../Project1/mulmod.h"
#include "../Division by multiplication 64/Magic64.h"
#include "../RuttenEekelen/RuttenEekelen.h"
#include "../Pollard Rho trials/Pollard Rho Montgomery.h"
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


uint64_t nextcoprime(uint64_t n) {
	static bool residues[210] = { false };
	static bool initialised = false;
	if (!initialised) {
		residues[0] = true;
		for (int p : {2, 3, 5, 7}) {
			for (int i = p; i < 210; i += p)
				residues[i] = true;
		}
		initialised = true;
	}
	int r = n % 210;
	while (residues[r])
		r++;
	return (n / 210) * 210 + r;
}
// performance comparison of is_prime and   FJ64_262K is_prime_2_64
void performance_trial() {

	uint64_t f;

	printf("\t\t\t\tis_prime  \t is_prime_2_64 \t is_primeFJ \t  is_prime_ref \n");
	for (int index = 11; index < 64; index+=4) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);

		// create a random set of numbers, adjusting to make them all coprime to 2.3.5.7=210
		std::mt19937_64 mt;
		mt.seed(std::random_device{}());
		mt.discard(700000);  // http://www.iro.umontreal.ca/~lecuyer/myftp/papers/lfsr04.pdf
		std::vector<uint64_t> samples;
		int trials = 3000000;
		for (int i = 0; i < trials; i++) {
			samples.push_back( nextcoprime ( (1ull << index) + mt() % interval) );
		}


		clock_t startTime = clock();
		int primecount1 = 0;
		for (auto n : samples) {
			if (is_prime(n, f))
				primecount1++;
		}
		clock_t endTime = clock();
		printf("is_prime     : 2^%d ... 2^%d\t %8.2f\t", index, index+1,
			double(endTime - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);

		startTime = clock();
		int primecount2 = 0;
		for(auto n : samples) {
			if (is_prime_2_64(n))
				primecount2++;
		}
		endTime = clock();
		printf("%8.2f\t", double(endTime - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);



		startTime = clock();
		int primecount3 = 0;
		for (auto n : samples) {
			if (is_primeFJ(n, f))
				primecount3++;
		}
		endTime = clock();
		printf("%8.2f\t", double(endTime - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);

		startTime = clock();
		int primecount4 = 0;
		for (auto n : samples) {
			if (is_prime_ref(n, f))
				primecount4++;
		}
		endTime = clock();
		printf("%8.2f\t%d\n", double(endTime - startTime) / (double)CLOCKS_PER_SEC, primecount1); fflush(stdout);


		if (primecount1 != primecount2 || primecount2 != primecount3 || primecount2 != primecount4)
		{
			printf("mismatch in number of primes found\n"); exit(1);
		}
	

	}
}

void verify() {
	std::mt19937_64 mt;
	mt.seed(std::random_device{}());
	uint64_t f;

	for (int index = 7; index < 64; index+=4) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);
		int primecount = 0;
		for (int count = 0; count < 1000000; count++) {
			uint64_t n = (1ULL << index) + mt() % interval;
			bool n_is_prime = is_prime_ref(n, f);
			if (n_is_prime != is_prime_2_64(n))
			{
				printf("%llu, discrepancy is_prime_2_64\n", n); exit(1);
			}
			if (n_is_prime != is_primeFJ(n, f))
			{
				printf("%llu, discrepancy is_primeFJ\n", n); exit(1);
			}
		}
		printf("2^%d verified\n", index);
	}
}

void longrun() {
	std::mt19937_64 mt;
	mt.seed(std::random_device{}());
	uint64_t primecount = 0;
	for (int index = 63; index < 64; index++) {
		uint64_t interval = (1ULL << (index + 1)) - (1ULL << index);
		uint64_t primecount = 0;
		for (int count = 0; count < 10000000000; count++) {
			uint64_t n = (1ULL << index) + mt() % interval;
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
	printf("\ncompare performance of is_SPRP (using Montgomery multiplication) and  is_primeFJ (using Rutten-Eekelen)");;
	for (int i = 39; i < 64; i += 4) {
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

			uint64_t pu,pv;
			xbinGCD(1ull << 63, n, &pu, &pv);

			if (is_SPRP(2ull, n,pv) && is_SPRP(bases[hashh(n)], n,pv))
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

void identify_RE_gt_MG() {
	Timer tim;
	//find numbers where Rutten Eekelen takes longer than Montgomery
	printf("\ncompare performance of is_SPRP (using Montgomery multiplication) and  is_primeFJ (using Rutten-Eekelen)");;
	for (int i = 39; i < 64; i += 4) {
		int primes = 0;
		printf("numbers starting 2**%d\n", i);
		uint64_t start = (1ull << i) + 1;
		for (uint64_t n = start; n < start + 10000000; n += 2) {

			if (n == 2 || n == 3 || n == 5 || n == 7) {
				primes++; continue;
			}
			if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0) continue;
			if (n < 121) {
				if (n > 1) primes++;
				continue;
			}


			tim.start();
			if (is_SPRP(2ull, n) && is_SPRP(bases[hashh(n)], n))
				primes++;
			tim.end();
			int64_t SPRP_time = tim.microsec();

			uint64_t factor = 0;
			tim.start();
			if (is_primeFJ(n, factor))
				primes++;
			tim.end();
			int64_t is_primeFJ_time = tim.microsec();

			if (is_primeFJ_time - SPRP_time > 100  || SPRP_time - is_primeFJ_time > 100)
				printf("%llu %lld %lld \n", n, SPRP_time, is_primeFJ_time);
		}
	}
}

int main(){

	//verify();
	//simple_is_prime_compare();
	//simple_sprp_compare();
	performance_trial();
	//longrun();
	//identify_RE_gt_MG();

	/*
	uint64_t mod = 1000000LL;
	std::mt19937_64 mt;
	mt.seed(std::random_device{}());
	for (int index = 7; index<20; index++){
		countcomposites = 0;
		factorcount = 0;
		mod *= 10LL;
		clock_t startTime = clock();
		for (int i = 0; (i<100000) || (factorcount == 0); i++){
			uint64_t p = mt() % mod;
			if (p >= mod / 10)
				examine(p,mod);
			}
		printf("10^%d: %d composites, %d factors, %8.2f sec\n", index, countcomposites, factorcount,
			   double(clock() - startTime) / (double)CLOCKS_PER_SEC); fflush(stdout);
	}
	*/
	return 0;
}