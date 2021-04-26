#include <stdint.h>

uint64_t ipow(int x, int p)
	{
		if (p == 0) return 1;
		if (p == 1) return x;

		uint64_t tmp = ipow(x, p / 2);
		if (p % 2 == 0) return tmp * tmp;
		else return x * tmp * tmp;
	}