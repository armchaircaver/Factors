#include <vector>
#include <miller rabin.h>
#include "../FactorsA/FactorsA.h"
#include "Pollard Rho Montgomery.h"
#include <ctime>

std::vector<uint64_t>primes;

void construct_primes(int power, int numprimes=100) {
	uint64_t n = (1ULL << power) + 1ULL;
	uint64_t f;
	do {
		if (is_prime(n, f)) {
			primes.push_back(n);
			//printf("%llu\n", n);
		}
		n += 2;
	} while (primes.size() < numprimes);
}


void test_brent() {
	clock_t startTime = clock();

	uint64_t fac;
	for (int i = 0; i < primes.size(); i++) {
		uint64_t p = primes[i];
		for (int j = i + 1; j < primes.size(); j++) {
			uint64_t q = primes[j];
			uint64_t pq = p * q;
			fac = brent(pq);
			//printf("%llu: factor %llu\n", pq, fac);
			if (pq % fac) {
				printf("Not a factor\n");
				exit(1);
			}
		}
	}
	clock_t endTime = clock();
	printf("Brent tests completed  %8.2f sec\n", double(endTime - startTime) / (double)CLOCKS_PER_SEC);
}


void test_pollard_brent_mont() {
	clock_t startTime = clock();
	uint64_t fac;
	for (int i = 0; i < primes.size(); i++) {
		uint64_t p = primes[i];
		for (int j = i + 1; j < primes.size(); j++) {
			uint64_t q = primes[j];
			uint64_t pq = p * q;
			fac = pollard_brent_montgomery_retry(pq);
			//printf("%llu: factor %llu\n", pq, fac);
			if (pq % fac) {
				printf("%llu Not a factor of %llu\n", fac, pq);
			}
			if (fac == 1) {
				printf("factor returned 1 (pq=%llu)\n", pq);
			}
			if (fac == pq) {
				printf("factor returned pq=%llu\n", pq);
			}
		}
	}
	clock_t endTime = clock();
	printf("Pollard Brent Montgomery tests completed  %8.2f sec\n", double(endTime - startTime) / (double)CLOCKS_PER_SEC);
}


void lengthy_test_pollard_brent() {
	uint64_t fac;
	for (int k = 10; k < 32; k++) {

		// generate primes
		primes.clear();
		construct_primes(k);
		clock_t startTime = clock();

		for (int i = 0; i < primes.size(); i++) {
			uint64_t p = primes[i];
			for (int j = i + 1; j < primes.size(); j++) {
				uint64_t q = primes[j];
				uint64_t pq = p * q;
				fac = pollard_brent_montgomery_retry(pq);
				//printf("%llu: factor %llu\n", pq, fac);
				if (pq % fac) {
					printf("%llu Not a factor of %llu\n", fac, pq);
				}
				if (fac == 1) {
					printf("factor returned 1 (pq=%llu)\n", pq);
				}
				if (fac == pq) {
					printf("factor returned pq=%llu\n", pq);
				}
			}
		}	
		clock_t endTime = clock();
		double elapsed = double(endTime - startTime) / (double)CLOCKS_PER_SEC;
		printf("Pollard Brent test completed for 2**%d > p,q > 2**%d, %8.2f sec\n", k+1, k, elapsed);
	}
	printf("Pollard Brent tests completed\n");
}

void test_powers() {
	uint64_t f;
	for (int e = 2; e < 32; e++) {
		printf("Testing primes between 2**%d and 2**%d\n", e, e + 1);
		primes.clear();
		construct_primes(e, 10000);
		for (auto p : primes) {
			uint64_t n = p * p;
			printf("testing pollard for %llu**2 ... ", p);
			uint64_t fac = pollard_brent_montgomery_retry(n);
			printf("factor %llu found\n", fac);


			n = p * p * p;
			if (!is_prime(n, f)) {
				printf("testing pollard for %llu**3 = %llu ... ", p, n);
				fac = pollard_brent_montgomery_retry(n);
				printf("factor %llu found\n", fac);
			}

			n = p * p * p * p;
			if (!is_prime(n, f)) {
				printf("testing pollard for %llu**4 = %llu ... ", p, n);
				fac = pollard_brent_montgomery_retry(n);
				printf("factor %llu found\n", fac);
			}
		}
	}
	printf("Powers tests completed\n");
}

void test_powers2() {
	for (uint64_t m = 3; m < (1ULL << 32); m += 2ULL) {
		uint64_t n = m;

		// highest exponent of m to give a number < 2**64
		int max_e = (int)(64.0 * log(2) / log((double)m));
		
		if (m%100000==1)
			printf("Examining powers of %llu\n", m);
		
		int e = 5;
		while (e < max_e) {
			n *= m;
			e++;
			//printf("testing pollard_brent_montgomery_retry for %llu = %llu**%d ... ", n, m, e);
			uint64_t fac = pollard_brent_montgomery_retry(n);
			int tries = getTries();
			if (tries>1)
				printf("factor %llu of %llu**%d found, %d tries\n", fac,m,e, tries );
			if ( n % fac) {
				printf("Invalid factor\n");
				exit(1);
			}
		}
	}
	printf("Powers tests 2 completed\n");
}

int main() {
	//test_powers2();
	for (int i = 10; i < 32; i++) {
		printf("Numbers starting at 2**%d\n", i);
		primes.clear();
		construct_primes(i);
		test_brent();
		test_pollard_brent_mont();
	}
	lengthy_test_pollard_brent();
}