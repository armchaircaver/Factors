#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <ctime>
#include <string>
#include "mulmod.h"
#include "Magic64.h"
#include <windows.h>
#include "montgomery.h"
#include "libdivide.h"
#include "ipow.h"
#include "divluh.h"
#include <algorithm>    // std::swap
#include "RuttenEekelen.h"

uint64_t mulmod1(uint64_t x, uint64_t y, uint64_t n){
	if (n < 0x100000000)
		return(((x%n)*(y%n) % n));

	uint64_t ret = 0ULL;
	if (x>n) x %= n;
	if (y>n) y %= n;
	if (msb(x) + msb(y) < 63)
		return (x*y) % n;

	if (n >= 0x8000000000000000)
		return mulmodSA(x, y, n);

	while (y) {
		if (y & 1) {
			ret += x;
			if (ret >= n) ret -= n;
			//printf("mulmod ret=%llu\n", ret);
		}
		y >>= 1;
		x <<= 1;
		if (x >= n) x -= n;
	}
	return ret;
}

// "Russian peasant" multiplication, speeded up 
uint64_t mulmod2(uint64_t a, uint64_t b, uint64_t m) {
	if (a >= m) a %= m;
	if (b >= m) b %= m;
	uint64_t res = 0;
	while (a != 0) {
		if (a & 1) {
			res += b;
			if (res >= m) res -= m;
		}
		a >>= 1;
		b <<= 1;
		if (b >= m) b -= m;
	}
	return res;
}

// from http://stackoverflow.com/questions/12168348/ways-to-do-modulo-multiplication-with-primitive-types
// modified
uint64_t mulmod4(uint64_t a, uint64_t b, uint64_t n)
{
	if (n < 0x100000000)
		return(((a%n)*(b%n) % n));

	int ms = msb(n);

	if (ms == 64){
		//printf("mulmod4 calling mulmodSA\n");
		return mulmodSA(a, b, n);
	}

	uint64_t result = 0;
	int N = 64 - ms;
	uint64_t mask = (1ULL << N) - 1ULL;

	if (a >= n) a %= n;
	if (b >= n) b %= n;  // Make sure all values are originally in the proper range?

	if (a == 0)
		return 0;

	while (b){
		result += (b & mask)*a;
		if (result >= n) result %= n;

		b >>= N;
		a <<= N;
		if (a >= n) a %= n;
	}
	return result;
}

/*
uint64_t fastmod(uint64_t a, struct libdivide_u64_t fastn){
return a - libdivide_u64_do(a, &fastn)*a;
}


uint64_t mulmod4ld(uint64_t a, uint64_t b, uint64_t n, struct libdivide_u64_t fastn)
{
//printf("mulmod4ld fastn = %llu, %d\n", fastn.magic, fastn.more);
if (n < 0x100000000)
// return(((a%n)*(b%n) % n));
if ( ((a%n)*(b%n)) % n != fastmod(fastmod(a, fastn)*fastmod(b, fastn), fastn)){
printf("discrepancy %llu, %llu\n", ((a%n)*(b%n)) % n, fastmod(fastmod(a, fastn)*fastmod(b, fastn), fastn)) ;
printf("%llu, %llu\n", a%n, fastmod(a, fastn));
printf("%llu, %llu\n", b%n, fastmod(b, fastn));
uint64_t x = fastmod(a, fastn)*fastmod(b, fastn);
printf("%llu, %llu\n", fastmod(a, fastn)*fastmod(b, fastn), fastmod(fastmod(a, fastn)*fastmod(b, fastn), fastn));
printf("%llu, %llu\n", x, fastmod(x, fastn));
exit(1);
}
return fastmod(fastmod(a, fastn)*fastmod(b, fastn), fastn);

int ms = msb(n);

if (ms == 64)
return mulmodSA(a, b, n);

uint64_t result = 0;
int N = 64 - ms;
uint64_t mask = (1ULL << N) - 1ULL;

if (a >= n) a = fastmod(a, fastn);
if (b >= n) b = fastmod(b, fastn);   // Make sure all values are originally in the proper range?

if (a == 0)
return 0;

while (b){
result += (b & mask)*a;
if (result >= n) result = fastmod(result, fastn);

b >>= N;
a <<= N;
if (a >= n) a = fastmod(a, fastn);
}
return result;
}
*/
uint64_t mulmod_divluh(uint64_t p, uint64_t q, uint64_t m){
	if (p >= m) p %= m;
	if (q >= m) q %= m;
	uint64_t high;
	//printf("calling umul128 p=%llu, q=%llu, n=%llu \n", p, q, m);
	uint64_t low = _umul128(p, q, &high);
	//printf(" umul128 returns low=%llu, high=%llu \n", low, high);
	uint64_t remainder;
	uint64_t quotient = divluh64(high, low, m, remainder);

	return remainder;
}





// from http://apps.topcoder.com/forums/;jsessionid=4BF5C317BC7B4253C4D7F3FB9ACE2BB9?module=Thread&threadID=670443&start=0&mc=9#1221719
// doesn't appear to work for m>=2^63
// seems slower than the best algorithms
uint64_t modmup(uint64_t a, uint64_t b, uint64_t m)
{
	if (a>m)
		a = a%m;
	if (b>m)
		b = b%m;
	uint64_t ret = 0;
	if (a<b)
		std::swap(a, b);
	while (b)
	{
		while (a<m)
		{
			if (b & 1)
				ret += a;
			a <<= 1;
			b >>= 1;
		}
		a -= m;
		while (ret >= m)
			ret -= m;
		if (a<b)
			std::swap(a, b);
	}
	return ret;
}

// from http://stackoverflow.com/questions/12168348/ways-to-do-modulo-multiplication-with-primitive-types
// modified
uint64_t mulmod4m(uint64_t a, uint64_t b, struct mu mag)
{
	uint64_t n = mag.n;

	if (n < 0x100000000)
		// return(((a%n)*(b%n) % n));
		return magicmod(magicmod(a, mag)*magicmod(b, mag), mag);

	int ms = msb(n);

	if (ms == 64)
		return mulmodSA(a, b, n);

	uint64_t result = 0;
	int N = 64 - ms;
	uint64_t mask = (1ULL << N) - 1ULL;

	if (a >= n) a = magicmod(a, mag);
	if (b >= n) b = magicmod(b, mag);   // Make sure all values are originally in the proper range?

	if (a == 0)
		return 0;

	if (b>a) std::swap(a, b);

	while (b){
		//assert(result + (b & mask)*a >= result);
		result += (b & mask)*a;
		if (result >= n) result = magicmod(result, mag);

		b >>= N;
		a <<= N;
		if (a >= n) a = magicmod(a, mag);
	}
	return result;
}


void testfunc(uint64_t mbase, uint64_t(*f)(uint64_t, uint64_t, uint64_t), std::string fname){
	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double Elapsed;

	QueryPerformanceFrequency(&Frequency);

	uint64_t m = mbase;
	//clock_t begin = clock();
	QueryPerformanceCounter(&StartingTime);
	uint64_t x = 0;
	//uint64_t abase = mbase / 2;
	//uint64_t bbase = mbase / 2;
	for (uint64_t m = mbase; m < mbase + 1000; m+=2){
		for (uint64_t b = mbase - 1000; b < mbase; ++b) {
			for (uint64_t a = mbase - 1000; a < mbase; ++a) {
				x += f(a, b, m);
			}
		}
	}
	QueryPerformanceCounter(&EndingTime);
	Elapsed = double(EndingTime.QuadPart - StartingTime.QuadPart) / double(Frequency.QuadPart);
	//clock_t end = clock();
	//double elapsed = double(end - begin) / CLOCKS_PER_SEC;
	printf("%s x=%llu, time = %f \n", fname.c_str(), x, Elapsed);

}

void testfuncmagic(uint64_t mbase, uint64_t(*f)(uint64_t, uint64_t, struct mu), std::string fname){
	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double Elapsed;

	QueryPerformanceFrequency(&Frequency);

	uint64_t m = mbase;
	//clock_t begin = clock();
	QueryPerformanceCounter(&StartingTime);
	uint64_t x = 0;
	//uint64_t abase = mbase / 2;
	//uint64_t bbase = mbase / 2;
	for (uint64_t m = mbase; m < mbase + 1000; m += 2){
		struct mu mag = magicu(m);
		for (uint64_t b = mbase - 1000; b < mbase; ++b) {
			for (uint64_t a = mbase - 1000; a < mbase; ++a) {
				x += f(a, b, mag);
			}
		}
	}
	QueryPerformanceCounter(&EndingTime);
	Elapsed = double(EndingTime.QuadPart - StartingTime.QuadPart) / double(Frequency.QuadPart);
	//clock_t end = clock();
	//double elapsed = double(end - begin) / CLOCKS_PER_SEC;
	printf("%s x=%llu, time = %f \n", fname.c_str(), x, Elapsed);

}

void testfuncRE(uint64_t mbase, std::string fname) {
	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double Elapsed;

	QueryPerformanceFrequency(&Frequency);

	uint64_t m = mbase;
	//clock_t begin = clock();
	QueryPerformanceCounter(&StartingTime);
	uint64_t x = 0ull;

	for (uint64_t m = mbase; m < mbase + 2; m += 2) {
		struct RE re = RE_gen(m);
		uint64_t b = m - 5;
		for (uint64_t i = 1; i < 100000000; ++i) {
			b = mulmodRE(b, b, re);
		}
		x = b;
	}

	QueryPerformanceCounter(&EndingTime);
	Elapsed = double(EndingTime.QuadPart - StartingTime.QuadPart) / double(Frequency.QuadPart);
	printf("%s mbase=%llu, x=%llu, time = %f \n", fname.c_str(), mbase, x, Elapsed);


}


void testfuncMontgomery(uint64_t mbase, std::string fname){
	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double Elapsed;

	QueryPerformanceFrequency(&Frequency);

	uint64_t m = mbase;
	QueryPerformanceCounter(&StartingTime);
	uint64_t x = 0;

	for (uint64_t m = mbase; m < mbase + 2; m += 2) {
		monty_t M = prepareMonty(m);
		uint64_t b = m - 5;
		uint64_t bbar = modul64(b, 0, M.m);
		for (uint64_t i = 1; i < 100000000; ++i) {
			bbar = montmul(bbar, bbar, M);
		}
		x = reverse(bbar, M);
	}

	QueryPerformanceCounter(&EndingTime);
	Elapsed = double(EndingTime.QuadPart - StartingTime.QuadPart) / double(Frequency.QuadPart);
	printf("%s mbase=%llu, x=%llu, time = %f \n", fname.c_str(), mbase, x, Elapsed);

}
/*
void testfunclibdivide(uint64_t mbase, std::string fname){
	LARGE_INTEGER StartingTime, EndingTime, Frequency;
	double Elapsed;

	QueryPerformanceFrequency(&Frequency);

	uint64_t m = mbase;
	QueryPerformanceCounter(&StartingTime);
	uint64_t x = 0;
	uint64_t abase = mbase / 2;
	uint64_t bbase = mbase / 2;
	for (uint64_t m = mbase; m < mbase + 10; m += 2){
		struct libdivide_u64_t fast_m = libdivide_u64_gen(m);
		for (uint64_t b = bbase; b < bbase + 1000; ++b){
			for (uint64_t a = abase; a < abase + 1000; ++a){
				x += mulmod4ld(a, b, m, fast_m);
			}
		}
	}
	QueryPerformanceCounter(&EndingTime);
	Elapsed = double(EndingTime.QuadPart - StartingTime.QuadPart) / double(Frequency.QuadPart);
	printf("%s x=%llu, time = %f \n", fname.c_str(), x, Elapsed);

}
*/

//-------------------------------------------  An implementation of Schrage's method -----------------------------
// from https://groups.google.com/forum/#!topic/comp.programming/Qw4rUccdcIw

// schrageFast2 assumes a <= b < m-1
uint64_t schrageFast2(uint64_t a, uint64_t b, uint64_t m) {

	if (a < 2ULL) {
		if (a) return b;
		return 0ULL;
	}

	uint64_t q = m / a;
	uint64_t r = m % a;

	uint64_t result;
	if (r >= q) {
		result = schrageFast2(a >> 1, b, m);
		result = addMod(result, result, m);
		if (a & 1ULL) { result = addMod(result, b, m); }
		return result;

	}
	else {
		// Schrage method
		uint64_t res_a = a * (b % q);
		uint64_t res_b = r * (b / q);
		if (res_b <= res_a)
			return res_a - res_b;
		else
			return res_a - res_b + m;
	}
}

uint64_t schrageFast(uint64_t a, uint64_t b, uint64_t m) {
	if (a >= m) a /= m;
	if (b >= m) b /= m;
	if (a > b) std::swap(a, b);

	uint64_t result;
	if (a < 2ULL) {
		return  (a * b) % m;
	}

	if (b == m - 1ULL) {
		return (m - a);
	}

	uint64_t q = m / a;
	uint64_t r = m % a;

	if (r >= q) {
		result = schrageFast2(a >> 1, b, m);
		result = addMod(result, result, m);
		if (a & 1ULL) { result = addMod(result, b, m); }
		return result;

	}
	else {
		// Schrage method
		uint64_t res_a = a * (b % q);
		uint64_t res_b = r * (b / q);
		if (res_b <= res_a)
			return res_a - res_b;
		else
			return res_a - res_b + m;
	}

}

int main(int argc, char **argv){

	for (int i = 7; i < 20; i+=4) {
		uint64_t mbase = ipow(10, i) + 1;
		uint64_t m = mbase;
		uint64_t abase = mbase / 2;
		uint64_t bbase = mbase / 2;
		for (uint64_t m = mbase; m < mbase + 10; m += 2) {
			monty_t M = prepareMonty(m);
				for (uint64_t b = m - 1000; b < m; ++b) {
					for (uint64_t a = m - 1000; a < m; ++a) {
						uint64_t x = mulmodSA(a, b, m);
						uint64_t y = mulmodAS(a, b, m);
						uint64_t z = montymulmod(a, b, M);
						if (x != y) {
							printf("MISMATCH for a=%llu,b=%llu,m=%llu, x=%llu,y=%llu\n", a, b, m, x, y);
						}
						if (x != z) {
							printf("MISMATCH mulmodSA and montymulmod for a=%llu,b=%llu,m=%llu, x=%llu,z=%llu\n", a, b, m, x, z);
						}
					}
				}
		}
		printf("verified i=%d\n", i);
	}


	for (int i = 10; i < 64; i+=2){
		uint64_t mbase = ipow(2, i) + 1;
		//printf("\n numbers starting at mbase=2^%d+1 = %llu \n",i, mbase);
		  
		//testfunc(mbase, mulmod,         "mulmod       ");
		//testfunc(mbase, mulmodSA, "mulmodSA     ");
		//testfunc(mbase, mulmodAS, "mulmodAS      ");
		//testfunc(mbase, mulmod4, "mulmod4      ");
		//testfunc(mbase, schrageFast,      "schrageFast  ");
		//testfunc(mbase, mulmod_divluh,  "mulmod_divluh");
		//testfunc(mbase, modmup,         "modmup       ");
		//testfuncmagic(mbase, mulmod4m,    "mulmod4m     ");
		//testfuncmagic(mbase, mulmodRR,     "mulmodRR      ");
		testfuncRE(mbase,                 "RE           ");
		testfuncMontgomery(mbase,         "Montgomery   ");
		//testfunclibdivide(mbase, "mulmodlibdivide");

	}
}