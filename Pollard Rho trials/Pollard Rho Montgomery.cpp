// source : http://coliru.stacked-crooked.com/a/f57f11426d06acd8
// referenced in https://projecteuler.chat/viewtopic.php?t=3776

#include <cstdint>
#include <iostream>
#include <tuple>

std::pair<uint64_t, uint64_t> mulu64(uint64_t a, uint64_t b) {
    uint64_t high, low;
    low = _umul128(a, b, &high);
    return std::make_pair(high, low);
}


// Computes (a + b) % m, assumes a < m and b < m.
inline uint64_t addmod64(uint64_t a, uint64_t b, uint64_t m) {
    if (b >= m - a) return a - m + b;
    return a + b;
}

// Computes aR * bR mod N with R = 2**64.
uint64_t montmul64(uint64_t a, uint64_t b, uint64_t N, uint64_t Nneginv) {
    std::pair<uint64_t, uint64_t> T, mN;

    uint64_t Th, Tl, m, mNh, mNl, th;

    std::tie(Th, Tl) = mulu64(a, b);
    m = Tl * Nneginv;
    std::tie(mNh, mNl) = mulu64(m, N);

    bool lc = Tl + mNl < Tl;
    th = Th + mNh + lc;
    bool hc = (th < Th) || (th == Th && lc);

    if (hc  || (th >= N) ) th = th - N;

    return th;
}

// Finds 2^-64 mod m and (-m)^-1 mod m for odd m (hacker's delight).
// equivalent to xbinGCD ?

inline std::pair<uint64_t, uint64_t> mont_modinv(uint64_t m) {
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

uint64_t prgcd(uint64_t a, uint64_t b) {
    if (a == b) return a;
    while (b > 0) {
        uint64_t tmp = a;
        a = b;
        b = tmp % b;
    }

    return a;
}

// Returns a factor of n, assumes n is odd.
uint64_t pollard_brent_montgomery(uint64_t n) {

    // Random number generator from D.E. Knuth
    static uint64_t rng = 0xdeafbeef;
    uint64_t a = rng * 6364136223846793005ull + 1442695040888963407ull;
    uint64_t b = a * 6364136223846793005ull + 1442695040888963407ull;
    rng = (a + b) ^ (a * b);

    // y and c are "montgomery space" numbers
    uint64_t y = 1 + a % (n - 1);
    uint64_t c = 1 + b % (n - 3); // modified to avoid -2 mod n
    uint64_t m = 100;
    //printf("n= %llu, y=%llu, c=%llu \n", n, y, c);

    // nneginv is m' (mprime) in Warren
    uint64_t nneginv = mont_modinv(n).second;
    

    uint64_t g=1, r, q, x, ys;
    q = r = 1;

    do {
        x = y;
        for (uint64_t i = 0; i < r; ++i) {
            y = addmod64(montmul64(y, y, n, nneginv), c, n);
        }

        for (uint64_t k = 0; k < r && g == 1; k += m) {
            ys = y;
            for (uint64_t i = 0; i < std::min(m, r - k); ++i) {
                y = addmod64(montmul64(y, y, n, nneginv), c, n);
                q = montmul64(q, x < y ? y - x : x - y, n, nneginv);
            }

            g = prgcd(q, n);
        }

        r *= 2;
    } while (g == 1);

    if (g == n) {
        do {
            ys = addmod64(montmul64(ys, ys, n, nneginv), c, n);
            g = prgcd(x < ys ? ys - x : x - ys, n);
        } while (g == 1);
    }

    return g;
}

int tries = 0;
uint64_t pollard_brent_montgomery_retry(uint64_t n) {
    uint64_t f=0;
    tries = 0;
    do {
        f = pollard_brent_montgomery(n);
        tries++;
     } while (f == n);
    return f;
}

int getTries() {
    return tries;
}
