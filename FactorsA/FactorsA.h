#include <stdint.h>

void factorise(uint64_t n, uint64_t *primearray, int &pasize);

void allfactors(uint64_t n, uint64_t *factorsarray, int &factorssize);

uint64_t totient(uint64_t n);

uint64_t squarefreepart(uint64_t n);

void setMaxSmallPrime(int i);

//uint64_t brent(uint64_t n);
uint64_t pollard_rhoRE(uint64_t n);
uint64_t pollard_rhoMU(uint64_t n);

//void sieve();

char getmethod();

void setVerbose(bool v);