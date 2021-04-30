// FactorsDLL.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdint.h>
#include <FactorsA.h>

#define DLLEXPORT extern "C" __declspec(dllexport)

DLLEXPORT void factor(uint64_t n, uint64_t *res, int &reslen){
	reslen = 0;
	factorise(n, res, reslen);
}
	
DLLEXPORT void allfactor(uint64_t n, uint64_t *res, int &reslen){
	reslen = 0;
	allfactors(n, res, reslen);
}

DLLEXPORT bool isprime(uint64_t n, uint64_t &factor){
	return is_prime(n, factor);
}

/*
// shouldn't need this
DLLEXPORT void Sieve(){
	sieve();
}
*/

DLLEXPORT uint64_t SquareFreePart(uint64_t n){
	return squarefreepart(n);
}

DLLEXPORT uint64_t Totient(uint64_t n){
	return totient(n);
}

DLLEXPORT char GetMethod(){
	return getmethod();
}

// to test the overhead of using vectors
DLLEXPORT void nulltest(uint64_t n, uint64_t *res, int &reslen){
	std::vector<uint64_t> f(64);
	reslen = 63;
	for (int i = 0; i < reslen; ++i){
		res[i] = f[i];
	}
}

// to test the overhead of calling C++ from Python
DLLEXPORT void nothing(uint64_t n, uint64_t *res, int &reslen){
	reslen = 0;
}



