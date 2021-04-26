#include <gmp.h>
#include <gmpxx.h>
#include "mpzPollard.h"
#include <vector>
#include "Timer.h"

int main() {

	// examples from https://stackoverflow.com/questions/19741967/how-can-i-improve-this-pollards-rho-algorithm-to-handle-products-of-semi-large
	mpz_class n;
	Timer t;

	std::vector<std::string> tests = { 
		/*
		"100000000000000001081", 
		"100000000000000001809",
		"100000000000000003049",	 
		"100000000000000003597" ,
		"100000000000000004363" , 
		"100000000000000004557" ,
		"100000000000000005067" ,
		"100000000000000008767", 
		"100000000000000006747",
		*/
	"10000000000000000000000183",
	"10000000000000000000000603",
	"10000000000000000000001389",
    "10000000000000000000001413",
    "10000000000000000000001603",
    "10000000000000000000001617",
    "10000000000000000000002301",
    "10000000000000000000002437"
	};

	printf("\npollardcl\n");
	uint64_t totaltime = 0;
	for (auto s : tests) {
		n = s;
		t.start();
		mpz_class f = pollardcl(n, false);
		t.end();
		if (n % f > 0) {
			gmp_printf("duff factor n=%Zd, factor=%Zd, %llu ms\n", n, f, t.ms());
			exit(1);
		}
		gmp_printf("n=%Zd, factor=%Zd, %llu ms\n", n, f, t.ms());
		totaltime += t.ms();
	}
	printf("totaltime %llu\n", totaltime);

	printf("\npollard1\n");
	totaltime = 0;
	mpz_t N;
	mpz_init(N);
	mpz_t F,Q;
	mpz_init(F);
	mpz_init(Q);
	for (auto s : tests) {
		mpz_set_str(N, s.c_str(), 10);
		t.start();
		pollard1(F,N, true);
		t.end();
		mpz_mod(Q, N, F);
		if (mpz_cmp_ui(Q,0) >0 ) {
			gmp_printf("duff factor N=%Zd, factor=%Zd, %llu ms\n", N, F, t.ms());
			exit(1);
		}
		gmp_printf("n=%Zd, factor=%Zd, %llu ms\n", N, F, t.ms());
		totaltime += t.ms();
	}
	printf("totaltime %llu\n", totaltime);

	printf("\npollard1 - four factors\n");
	for (int i = 10; i <= 15; i++) {
		mpz_class p = pow(10, i);
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		mpz_class q,r,s;
		mpz_nextprime(q.get_mpz_t(), p.get_mpz_t());
		mpz_nextprime(r.get_mpz_t(), q.get_mpz_t());
		mpz_nextprime(s.get_mpz_t(), r.get_mpz_t());
		mpz_class pqrs = p * q*r*s;
		int len_pqrs = (int)mpz_sizeinbase(pqrs.get_mpz_t(), 10);
		gmp_printf("pqrs=%Zd, (%d digits) ", pqrs, len_pqrs);
		t.start();
		pollard1(F, pqrs.get_mpz_t(), false);
		t.end();
		gmp_printf(", pollard1, factor=%Zd, %llu ms\n", F, t.ms());
	}

	printf("\npollard1 - unbalanced factors\n");
	for (int i = 10; i <= 15; i++) {
		mpz_class p = pow(10, i);
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		mpz_class q = pow(10,40);

		mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
		
		mpz_class pqrs = p * q;
		int len_pqrs = (int)mpz_sizeinbase(pqrs.get_mpz_t(), 10);
		gmp_printf("pq=%Zd, (%d digits) ", pqrs, len_pqrs);
		t.start();
		pollard1(F, pqrs.get_mpz_t(), false);
		t.end();
		gmp_printf(", pollard1, factor=%Zd, %llu ms\n", F, t.ms());
	}

	printf("\npollard1 - larger numbers\n");
	for (int i = 30; i <= 75; i++) {
		mpz_class p = pow(2, i);
		mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		mpz_class q;
		mpz_nextprime(q.get_mpz_t(), p.get_mpz_t());
		mpz_class pq = p * q;
		int len_pq = (int)mpz_sizeinbase(pq.get_mpz_t(), 10);
		gmp_printf("pq=%Zd, (%d digits), p=%Zd, q=%Zd ", pq, len_pq, p, q);
		t.start();
		pollard1(F, pq.get_mpz_t(), false);
		t.end();
		gmp_printf(", pollard1, factor=%Zd, %llu ms\n", F, t.ms());
	}
}
