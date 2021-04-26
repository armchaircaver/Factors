#include <stdio.h>
#include <stdint.h>
#include "libdivide.h"

uint64_t sum_of_residues64(const uint64_t *numers, int count, uint64_t d) {
	uint64_t result = 0;
	for (int i = 0; i < count; i++)
		result += numers[i] % d; //this division is slow!
	return result;
}

//Here is how you would optimize it with libdivide in C++:
uint64_t sum_of_residues_opt64(const uint64_t *numers, int count, uint64_t d) {
	uint64_t result = 0;
	//libdivide::divider<uint64_t> fast_d(d); //constructs an instance of libdivide::divider
	struct libdivide_u64_t fast_d = libdivide_u64_gen(d);
	for (int i = 0; i < count; i++)
		//result += numers[i]- (numers[i] / fast_d)*d; //uses faster libdivide division
		result += numers[i] - d*libdivide_u64_do(numers[i], &fast_d); // performs faster libdivide division
	return result;
}

int main(int argc, char **argv){
	uint64_t numerators[] = { 300, 500, 700, 900,  25000000000000000 };
	auto len = sizeof(numerators) / sizeof(*numerators);
	uint64_t d = 101;
	printf("sum_of_residues64 %llu\n", sum_of_residues64(numerators, len, d));
	printf("sum_of_residues_opt64 %llu\n", sum_of_residues_opt64(numerators, len, d));

}
