// DLL attempt.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <vector>

#define DLLEXPORT extern "C" __declspec(dllexport)

typedef unsigned long long u64;

DLLEXPORT uint64_t addone() {
	static uint64_t n = 0;
	n++;
	return n;
}

DLLEXPORT u64  knuthrand() {
	// Random number Linear Congruential Generator MMIX from D.E. Knuth
	static uint64_t r = 0xaccedeafacade;
	r = r * 6364136223846793005ull + 1442695040888963407ull;
	return r;
}


DLLEXPORT u64 ipow(u64 x, int p){
	if (p == 0) return 1;
	if (p == 1) return x;

	u64 tmp = ipow(x, p / 2);
	if (p % 2 == 0)return tmp * tmp;
	else return  x * tmp * tmp;
}

DLLEXPORT void ipow2(int x, int p, u64 &res){
	res = ipow(x, p);
}

DLLEXPORT void ipow3(int x, int p, u64 *res, int &reslen ){
	u64 y = ipow(x, p);
	reslen = p % 4;
	for (int i = 0; i < reslen; ++i){
		res[i] = y + i;
	}
}
