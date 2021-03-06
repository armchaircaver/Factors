#include <stdint.h>

bool is_prime(uint64_t n, uint64_t &factor);

bool is_primeFJ(uint64_t n, uint64_t& factor);

uint64_t gcd(uint64_t a, uint64_t b);

uint64_t pow_mod(uint64_t x, uint64_t y, uint64_t n);

bool miller_rabin_pass(uint64_t a, int s, uint64_t d, uint64_t n, uint64_t& factor);

bool is_prime_ref(uint64_t n, uint64_t& factor);

