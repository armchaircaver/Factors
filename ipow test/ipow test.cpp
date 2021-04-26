#include <stdio.h>
#include "ipow.h"

void ipow2(int x, int p, uint64_t &res){
	res = ipow(x, p);
}

int main(int argc, char **argv){
	for (int i = 0; i < 64; ++i){
		printf("\n %d, %llu", i, ipow(2, i));
	}
	uint64_t res;
	for (int i = 0; i < 64; ++i){
		ipow2(2, i, res);
		printf("\n %d, %llu", i, res);
	}
}