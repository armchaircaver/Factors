#pragma once
#include <gmp.h>
#include <gmpxx.h>


mpz_class pollardcl(mpz_class n, bool verbose = false);
void pollard1(mpz_t& res, mpz_t n, bool verbose = false, float timeput=1.0);