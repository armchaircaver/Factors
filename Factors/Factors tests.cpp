#include <vector>
#include <stdio.h>
#include <stdint.h>
#include "Factors.h"
#include "ipow.h"
#include <ctime>	
#include <windows.h>

void printvector(std::vector<uint64_t> v){
	printf("[");
	for (auto it = v.begin(); it != v.end(); ++it) {
		if (it == v.end() - 1)
			printf("%llu", *it);
		else
			printf("%llu, ", *it);
	}
	printf("]");
}

void test(uint64_t n){
	printf("\nfactors of %llu: factorise: (", n);
	std::vector<uint64_t> f = factorise(n);
	printvector(f);
}

void verifyfactorisation(uint64_t n){
	std::vector<uint64_t> f = factorise(n);
	uint64_t c = 1;
	for (auto it = f.begin(); it != f.end(); ++it) {
		c *= *it;
	}
	if (c != n){
		printf("\n MISMATCH %llu", n);
		printvector(f);
		printf("\n");
		exit(1);
	}
}

void timetrial(int maxsp){
	printf("\nMAXSP\tbase\tnumbers\tmicrosec/number");
	setMaxSmallPrime(maxsp); 
	for (int e = 5; e <= 19; e++){
		uint64_t base = ipow(10, e);
		double elapsed = 0.0;
		clock_t begin = clock();
		int n = 0;
		for (; elapsed <= 2.0;){
			verifyfactorisation(base + n);
			clock_t end = clock();
			elapsed = double(end - begin) / CLOCKS_PER_SEC;
			n++;
		}
		printf("\n%d\t10^%d\t%d\t %10.1f", maxsp, e, n, elapsed*1000000.0 / double(n));
	}
}

int main(int argc, char **argv){
	clock_t begin = clock();
	sieve();
	clock_t end = clock();
	setMaxSmallPrime(10);
	printvector(factorise(100000000000000006));

	setMaxSmallPrime(10000);
	// some semi primes to test the algorithm
	std::vector<uint64_t> semis = { 18446743721522234449, 10000000000000000049, 8980935344490257, 100000000000000009 };
	for (uint32_t i = 0; i < semis.size(); i++){
		uint64_t n = semis[i];
		clock_t begin = clock();
		printf("\n %llu = ", n);
		printvector(factorise(n));
		clock_t end = clock();
		printf(" %10.6f s", double(end - begin) / CLOCKS_PER_SEC);

	}

	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double ElapsedMicroseconds;

	QueryPerformanceFrequency(&Frequency);

	// mersenne numbers test
	uint64_t p2 = 1;
	for (int n = 1; n < 65; ++n){
		p2 *= 2;
		uint64_t m = p2 - 1;
		printf("\n 2^%d-1\t", n);
		QueryPerformanceCounter(&StartingTime);
		std::vector<uint64_t> f = factorise(m);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds = double(EndingTime.QuadPart - StartingTime.QuadPart)*1000000.0 / double(Frequency.QuadPart);

		printvector(f);
		printf("\t%10.1f ", ElapsedMicroseconds);
	};

	timetrial(10);
	timetrial(100);
	timetrial(1000);
	timetrial(10000);
	timetrial(100000);
	timetrial(1000000);

	return 0L;
}