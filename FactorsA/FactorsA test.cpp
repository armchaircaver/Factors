#include <vector>
#include <stdio.h>
#include <stdint.h>
#include "FactorsA.h"
#include <ctime>	
#include <windows.h>
#include <algorithm>  // sort


void printarray(uint64_t *primearray, int pasize){
	printf("[");
	std::sort(primearray, primearray + pasize);
	for (int it = 0; it != pasize; ++it) {
		if (it == pasize - 1)
			printf("%llu", primearray[it]);
		else
			printf("%llu, ", primearray[it]);
	}
	printf("]");
}

void test(uint64_t n){
	int pasize;
	uint64_t primearray[64];
	printf("\nfactors of %llu: factorise: (", n);
	factorise(n, primearray, pasize);
	printarray(primearray, pasize);
}

void performfactorisation(uint64_t n){
	int pasize=0;
	uint64_t primearray[64];
	factorise(n, primearray, pasize);
}

uint64_t knuthrand() {
	// Random number Linear Congruential Generator MMIX from D.E. Knuth
	static uint64_t r = 0xaccedeafacade;
	r = r * 6364136223846793005ull + 1442695040888963407ull;
	return r;
}


void smallprimetrial() {
	//printf("\nMAXSP\tbase\tnumbers\tmicrosec/number");

	// construct a list of smooth numbers
	printf("constructing smooth numbers...");
	std::vector<uint64_t> smoothnumbers;
	uint64_t base = 1000000000000000000ull;
	for (int j = 0; j < 1000; j++)
		for (auto p : { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
						73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
						179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,
						283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409 })
		{
			uint64_t n = 1ull;
			int e = (int)(log(1.8446744e+19) / log((double)p));
			for (int i = 0; i < e; i++) {
				n *= p;
				smoothnumbers.push_back(n);
			}
		}
	printf("completed. %d smooth numbers.\n", (int)smoothnumbers.size() );

	// construct random numbers
	std::vector<uint64_t> randomnumbers;

	// keep doubling p2 until it wraps around to zero
	for (uint64_t p2 = 1024ull; p2; p2 <<= 1) 
		for (int i = 0; i < 10000; i++) 
			randomnumbers.push_back( knuthrand() % p2 );

	for (int i = 0; i < 20000; i++)
		randomnumbers.push_back( knuthrand() );
	printf("created %d random numbers\n", (int)randomnumbers.size());


	printf("smooth number performance:\n");
	for (auto numsp : { 3, 7, 10, 30, 60, 100 }) {
		setNumSmallPrimes(numsp);
		printf("NumSmallPrimes=%d: ", numsp);
		clock_t begin = clock();
		for (auto n : smoothnumbers)
			performfactorisation(n);
		clock_t end = clock();
		double elapsed = double(end - begin) / CLOCKS_PER_SEC;
		printf(" %4.1f sec\n", elapsed);
	}

	printf("random number performance:\n");
	for (auto numsp : { 3, 7, 10, 30, 60, 100 }) {
		setNumSmallPrimes(numsp);
		printf("NumSmallPrimes=%d: ", numsp);
		clock_t begin = clock();
		for (auto n : randomnumbers)
			performfactorisation(n);
		clock_t end = clock();
		double elapsed = double(end - begin) / CLOCKS_PER_SEC;
		printf(" %4.1f sec\n", elapsed);
	}

}


int main(int argc, char **argv){
	clock_t begin = clock();
	//sieve();
	clock_t end = clock();
	//printf("called sieve\n");
	int pasize = 0;
	uint64_t primearray[64];

	
	printf("\nTesting totient\n");
	for (uint64_t n = 1; n < 100; n += 7){
		printf("%llu, totient = %llu\n", n, totient(n));
	}

	printf("\nTesting square free part\n");
	for (uint64_t n = 1; n < 100; n += 7){
		printf("%llu, squarefreepart = %llu\n",n,squarefreepart(n));
	}


	printf("3,4,5, squarefreepart = %llu\n", squarefreepart(576ULL));
	printf("5,12,13, squarefreepart = %llu\n", squarefreepart(14400ULL));

	int allfacsize = 0;
	uint64_t* allfactorsarray = new uint64_t[193536] ; 

	printf("\nTesting large semiprimes\n");

	for (int i = 0; i < 10; i++){

		factorise(8333333333333333341uLL, primearray, pasize);
		printarray(primearray, pasize);
	}
	for (uint64_t n : {125, 72, 600}) {

		printf("\ncalling allfactors\n");
		allfactors(n, allfactorsarray, allfacsize);
		printf("allfactors for %llu: ", n);
		printf("allfacsize = %d\n", allfacsize);
		printarray(allfactorsarray, allfacsize); printf("\n");
	}

	uint64_t n = 18401055938125660800uLL;
	printf("\n %llu=", n);
	factorise(n, primearray, pasize);
	printarray(primearray, pasize); 
	allfactors(n, allfactorsarray, allfacsize);
	printf(", %d factors\n", allfacsize);
	if (allfacsize != 184320){
		printf("unexpected value for all factors size, expecting 184320");
		exit(1);
	}


	factorise(100000000000000006, primearray, pasize);
	printarray(primearray, pasize); printf("\n");

	factorise(720, primearray, pasize);
	printarray(primearray, pasize); printf("\n");

	// some semi primes to test the algorithm
	std::vector<uint64_t> semis = { 18446743721522234449, 10000000000000000049, 8980935344490257, 100000000000000009 };
	for (uint32_t i = 0; i < semis.size(); i++){
		uint64_t n = semis[i];
		clock_t begin = clock();
		printf("\n %llu = ", n);
		pasize=0;
		factorise(n, primearray, pasize);
		printarray(primearray, pasize);
		clock_t end = clock();
		printf(" %10.6f s", double(end - begin) / CLOCKS_PER_SEC);
		performfactorisation(n);

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
		pasize = 0;
		factorise(m, primearray, pasize);
		performfactorisation(m);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds = double(EndingTime.QuadPart - StartingTime.QuadPart)*1000000.0 / double(Frequency.QuadPart);

		printarray(primearray, pasize);
		printf("\t%10.1f ", ElapsedMicroseconds);
	};
	
	// carmichael numbers
	std::vector<uint64_t> A202562 = { 561, 84350561, 851703301, 2436691321, 34138047673, 60246018673, 63280622521,
		83946864769, 110296864801, 114919915021, 155999871721, 225593397919,
		342267565249, 534919693681, 660950414671, 733547013841, 1079942171239, 1301203515361,
		1333189866793 };
	for (uint32_t i = 0; i < A202562.size(); i++){
		uint64_t n = A202562[i];
		printf("\n %llu = ", n);
		pasize = 0;
		factorise(n, primearray, pasize);
		printarray(primearray, pasize);
		performfactorisation(n);
	}

	std::vector<uint64_t> A141768 = { 9, 25, 49, 91, 341, 481, 703, 1541, 1891, 2701, 5461, 6533,
		8911, 12403, 18721, 29341, 31621, 38503, 79003, 88831, 146611, 188191, 218791, 269011,
		286903, 385003, 497503, 597871, 736291, 765703, 954271, 1024651, 1056331, 1152271,
		1314631, 4234438493911, 4234787717131, 4234927410451, 4236953222503, 4240447139503, 4240866506311,
		4242893738203, 4244921454511, 4247159496751, 4247859005911, 4248488613403, 4249468095571, 4255137315703,
		4255207329691, 4258568678491, 4260460020391, 4261861285111, 4263543106903, 4279607175691 };
	for (uint32_t i = 0; i < A141768.size(); i++){
		uint64_t n = A141768[i];
		printf("\n %llu = ", n);
		pasize = 0;
		factorise(n, primearray, pasize);
		printarray(primearray, pasize);
		performfactorisation(n);
		if (pasize < 2){
			printf("pseudoprime not detected\n");
			exit(1);
		}
	}


	printf("\n");
	printf("\nfactorisation timing for 10000 numbers starting with 10^n n=5..19 for maxsp values\n");
	

	smallprimetrial();

	return 0L;
}