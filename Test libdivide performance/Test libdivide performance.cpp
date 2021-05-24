#include <stdint.h>
//#define LIBDIVIDE_SSE2 
#include "../FactorsA/libdivide.h"
#include <vector>
#include <ctime>

std::vector<uint64_t> primes = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
				73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
				179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,
				283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409 };

int main() {

#if defined(LIBDIVIDE_SSE2)
	printf("LIBDIVIDE_SSE2 defined\n");
#else
	printf("LIBDIVIDE_SSE2 not defined\n");

#endif

	printf("divisibility tests\n");

	int count = 0;
	clock_t begin = clock();

	for (auto p : primes) {
		uint64_t n = 1ull << 63;
		libdivide::divider<uint64_t> d(p);
		for (uint64_t m = 0; m < 1000000; m++) {
			n = (1ull << 63) + m;
			while (n) {
				uint64_t q = n / d;
				if (q*p==n) count++;
				n = q;
			}
		}
	}
	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC;

	printf("libdivide,    count=%d, %4.1f sec \n", count, elapsed );

	count = 0;
	begin = clock();
	uint64_t r = 0;
	for (auto p : primes) {
		uint64_t n = 1ull << 63;
		for (uint64_t m = 0; m < 1000000; m++) {
			n = (1ull << 63) + m;
			while (n) {
				n = _udiv128(0ull, n, p, &r);
				if (r==0) count++;
			}
		}
	}
	end = clock();
	elapsed = double(end - begin) / CLOCKS_PER_SEC;
	printf("conventional, count=%d, %4.1f sec \n", count, elapsed );


	printf("simple division tests\n");

	count = 0;
	begin = clock();

	for (auto p : primes) {
		uint64_t n = 1ull << 63;
		libdivide::divider<uint64_t> d(p);
		for (uint64_t m = 0; m < 1000000; m++) {
			uint64_t n = (1ull << 63) + m;
			while (n) {
				n /= d;
				count++;
			}
		}
	}
	end = clock();
	elapsed = double(end - begin) / CLOCKS_PER_SEC;

	printf("libdivide,    count=%d, %4.1f sec \n", count, elapsed);

	count = 0;
	begin = clock();
	for (auto p : primes) {
		uint64_t n = 1ull << 63;
		for (uint64_t m = 0; m < 1000000; m++) {
			n = (1ull << 63) + m;
			while (n) {
				n = n / p;
				count++;
			}
		}
	}
	end = clock();
	elapsed = double(end - begin) / CLOCKS_PER_SEC;
	printf("conventional, count=%d, %4.1f sec \n", count, elapsed);

	return 0;
}

	