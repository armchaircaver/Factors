// test whether the montgomery algorithm used in "pollard rho montgomery" 
// corresponds to the algorithm in montgomery.cpp

#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>
#include "../Project1/montgomery.h"
#include <intrin.h>
#pragma intrinsic(_umul128)
#pragma intrinsic(_addcarry_u64)

 std::pair<uint64_t, uint64_t> mulu64(uint64_t a, uint64_t b) {
    uint64_t high, low;
    low = _umul128(a, b, &high);
    return std::make_pair(high, low);
}

 void addcarry_tests() {
     unsigned char  carry=0;

     uint64_t x = 0ull, y = 0ull, res = 0;
     unsigned char carry2 = _addcarry_u64(carry, x, y, &res);
     printf("0+0 : carry2=%d\n", carry2);

     uint64_t a = 1ull << 63;
     for (uint64_t b = (1ull << 63) - 10; b < (1ull << 63) + 10; b++) {
         carry = _addcarry_u64(carry, a, b, &res);
         printf("carry test b=%llu, res=%llu, carry=%d\n", b, res, carry);
         carry = 0;
     }
     uint64_t Tl = 0xffffffffffffffff;
     uint64_t mNl = 2ull;
     uint64_t resl = 0;;
     uint32_t carryL = _addcarry_u64(0, Tl, mNl, &resl);
     uint64_t Th = 1;
     uint64_t mNh = 0xffffffffffffffff;
     uint64_t resh;
     uint32_t carryH = _addcarry_u64(carryL, Th, mNh, &resh);
     printf("carryL=%d, carryH=%d, Th=%llu, mNh=%llu, resh=%llu\n", carryL, carryH, Th, mNh, resh);

 }

// Computes (a + b) % m, assumes a < m and b < m.
inline uint64_t addmod64(uint64_t a, uint64_t b, uint64_t m) {
    if (b >= m - a) return a - m + b;
    return a + b;
}

// Computes aR * bR mod N with R = 2**64.
_declspec(noinline) uint64_t montmul64(uint64_t a, uint64_t b, uint64_t N, uint64_t Nneginv) {

    uint64_t Th,  mNh ;

    uint64_t Tl = _umul128(a, b, &Th);
    uint64_t m = Tl * Nneginv;
    uint64_t mNl = _umul128(m, N, &mNh);
    //printf("after multiplications: a=%llu, b=%llu, N=%llu, nneginv=%llu, Th=%llu, Tl=%llu, mNl=%llu, mNh=%llu\n", a,b, N, Nneginv,Th, Tl, mNl, mNh);
    
    /*
    bool lc = Tl + mNl < Tl;
    uint64_t th = Th + mNh + lc;
    bool hc = (th < Th) || (th == Th && lc);
    if (hc || (th >= N)) th = th - N;
    */

    // my attempt at speeding up- can't get intrinsics to work if adding 0+0, carryL is set to 2.
    
    uint64_t tl2, th2;
    uint32_t low_carry_in = 0;
 
    uint32_t carryL = _addcarry_u64(low_carry_in, Tl, mNl, &tl2);
 
    uint32_t carryH = _addcarry_u64(carryL, Th, mNh, &th2);
   
    if (carryH || (th2 >= N) ) th2 = th2 - N;

    return th2;
}

// Finds 2^-64 mod m and (-m)^-1 mod m for odd m (hacker's delight).
// equivalent to xbinGCD ?

std::pair<uint64_t, uint64_t> mont_modinv(uint64_t m) {
    uint64_t a = 1ull << 63;
    uint64_t u = 1;
    uint64_t v = 0;

    while (a > 0) {
        a = a >> 1;
        if ((u & 1) == 0) {
            u = u >> 1; v = v >> 1;
        }
        else {
            u = ((u ^ m) >> 1) + (u & m);
            v = (v >> 1) + (1ull << 63);
        }
    }
    return std::make_pair(u, v);
}

void largeprimetest() {
    std::vector<uint64_t> largeprimes = { 1073741827,1073741831,1073741833,1073741839,1073741843,268435459,268435463,268435493,268435523,268435537,268435561 };

    for (auto a : largeprimes)
        for (auto b : largeprimes) {
            uint64_t n = a * b;
            uint64_t nneginv = mont_modinv(n).second;
            printf("n=%llu = %llu*%llu\n", n, a, b);
            uint64_t y1 = 1000;
            for (int i = 0; i < 100000000; i++) {
                y1 = montmul64(y1, y1, n, nneginv);
            }
            printf("final y1 = % llu\n", y1);


        }
}

int main() {
    /* some primes to play with
1073741827 1073741831 1073741833 1073741839 1073741843 1073741857 1073741891 1073741909 1073741939 1073741953
268435459 268435463 268435493 268435523 268435537 268435561 268435577 268435579 268435597 268435631
1031 1033 1039 1049 1051 1061 1063 1069 1087 1091
4099 4111 4127 4129 4133 4139 4153 4157 4159 4177
*/

addcarry_tests();

    uint64_t y1 = 569396587663458304ull;
    uint64_t n =  1152921511049297929ull;
    uint64_t nneginv = mont_modinv(n).second;
    y1 = montmul64(y1, y1, n, nneginv);
    printf("y1=%llu\n", y1);
    //exit(0);

    std::vector<uint64_t> primes = { 1073741827, 1031,1033, 1073741831,4099,4111268435459, 268435463 };
    for (auto a : primes)
        for (auto b : primes) {
            // construct a semiprime
            uint64_t n = a * b;
            printf("n=%llu = %llu*%llu : ", n, a, b);

            uint64_t nneginv = mont_modinv(n).second;
            monty_t M = prepareMonty(n);

            uint64_t y1 = 1000, y2 = 1000;

            // verify that the two implementations of montgomery multiplication 
            // produce the same results
            // and "montgomery space" numbers are in the range [1,n-1]
            for (int i = 0; i < 100000000; i++) {
                uint64_t prev_y1 = y1;
                y1 = montmul64(y1, y1, n, nneginv);
                y2 = montmul(y2, y2, M);
                if (y1 != y2) {
                    printf("ERROR1 y1=%llu, y2=%llu, prev_y1=%llu\n", y1, y2, prev_y1); exit(1);
                }
                if (y1 >= n) {
                    printf("ERROR2 y1=%llu >= n=%llu\n", y1, n); exit(1);
                }
                if (y1 == 0ull) {
                    printf("ERROR3 y1=%llu, y2=%llu\n", y1, y2); exit(1);
                }
            }
            printf("final y1 = % llu, y2 = % llu, tests completed\n", y1, y2);
        }
}