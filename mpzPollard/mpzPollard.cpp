#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include "Timer.h"

void pollard1(mpz_t &res, mpz_t n, bool verbose = false, float timeout=1.0) {

	if (mpz_even_p (n)) {
		mpz_set_ui(res, 2);
		return;
	}

	Timer tim;
	tim.start();

	const uint64_t m = 1000;
	mpz_t  x, ys, q, g, ab;
	mpz_init(x);
	mpz_init(ys);
	mpz_init(q);
	mpz_init(g);
	mpz_init(ab);
	mpz_set_ui(q, 1);
	mpz_set_ui(g, 1);

	uint64_t r = 1;

	if (verbose)
		gmp_printf("Calling mpz pollard1 rho for n=%Zd\n", n);
	uint64_t squareadd_count = 0;

	gmp_randclass  ran(gmp_randinit_default);
	ran.seed(time(NULL));
	mpz_class y0 = ran.get_z_range(4) + 1;

	mpz_t a, y;
	mpz_init(a);
	mpz_init(y);
	mpz_set_ui( a, 1); // ran.get_z_range(n - 3) + 1;
	mpz_set( y , y0.get_mpz_t() );

	while (mpz_cmp_ui (g ,1) ==0 ){
		mpz_set(x, y); // x = y;
		for (uint64_t i = 0LL; i < r; i++) {
			//y = (y * y + a) % n;
			mpz_powm_ui(y, y, 2, n);
			mpz_add(y,y,a);
			if (mpz_cmp(y, n) >= 0) mpz_sub(y, y, n); // y -= n;
			squareadd_count++;
		}

		uint64_t k = 0;
		while (k < r && mpz_cmp_ui(g, 1) == 0) {
			mpz_set(ys, y); // ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				//y = (y * y + a) % n;
				mpz_powm_ui(y, y, 2, n);
				mpz_add(y, y, a);  //y += a;
				if (mpz_cmp(y, n) >= 0) mpz_sub(y, y, n); // y -= n;
				squareadd_count++;

				// q = q |x-y| mod n
				if (mpz_cmp(x, y) > 0)
					mpz_sub(ab, x, y);
				else
					mpz_sub(ab, y, x);

				mpz_mul(q, q, ab);
				mpz_mod(q, q, n); // q %= n;
			}
			if (mpz_cmp_ui(q, 0) == 0) {
				printf("q=0\n");
				exit(1);
			}
			mpz_gcd(g, q, n); // g = gcd(q, n);
			k += m;
		}
		r *= 2;
		tim.end();
		if ((float)tim.ms() > timeout*1000.0) {
			mpz_set_ui(res, 0);
			return;
		}
	}

	if (g == n) {
		do {
			//ys = (ys * ys + a) % n;
			mpz_powm_ui(ys, ys, 2, n);
			mpz_add(y, y, a); // y += a;
			if (mpz_cmp(y, n) >= 0) mpz_sub(y, y, n); // y -= n;
			squareadd_count++;
			//g = gcd((x > ys) ? x - ys : ys - x, n);
			if (mpz_cmp(x, ys) > 0)
				mpz_sub(ab, x, ys);
			else
				mpz_sub(ab, ys, x);
			mpz_gcd(g, ab, n);
		} while (mpz_cmp_ui(g, 1) == 0 /*g == 1*/);
	}

	if (verbose) printf("square add count = %llu\n", squareadd_count);
	mpz_set(res, g);
	return ;
}


mpz_class pollardcl(mpz_class n, bool verbose = false) {

	if (n % 2 == 0)
		return 2;

	const uint64_t m = 1000;
	mpz_class  x, ys, q = 1;
	uint64_t r = 1;
	mpz_class g = 1;

	if (verbose)
		gmp_printf("Calling mpz pollardcl rho for n=%Zd, ", n);
	uint64_t squareadd_count = 0;

	gmp_randclass  ran(gmp_randinit_default);
	ran.seed(time(NULL));

	mpz_class a = 1; // ran.get_z_range(n - 3) + 1;
	mpz_class y = ran.get_z_range(4) + 1;

	while (g == 1) {
		x = y;
		for (uint64_t i = 0LL; i < r; i++) {
			y = (y * y + a) % n;
			squareadd_count++;
		}

		uint64_t k = 0;
		while (k < r && g == 1) {
			ys = y;
			for (uint64_t i = 0LL; i < m && i < r - k; i++) {
				y = (y * y + a) % n;
				squareadd_count++;

				// q = q |x-y| mod n
				q *= (x > y) ? x - y : y - x;
				q %= n;
			}
			g = gcd(q, n);
			k += m;
		}
		r *= 2;
	}

	if (g == n) {
		do {
			ys = (ys * ys + a) % n;
			squareadd_count++;
			g = gcd((x > ys) ? x - ys : ys - x, n);
		} while (g == 1);
	}

	if (verbose) printf("square add count = %llu\n", squareadd_count);
	return g;
}


