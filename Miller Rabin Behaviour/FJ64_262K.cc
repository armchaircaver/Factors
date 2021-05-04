#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <time.h> 
#include <random>
#include "FJ64_262K.h"
#include <immintrin.h>  // _div128

using namespace std;

typedef unsigned long long uint64; 
typedef long long int64; 
 


// Montgomery multiplication 
// http://www.hackersdelight.org/MontgomeryMultiplication.pdf
void xbinGCD(uint64 a, uint64 b, uint64 *pu, uint64 *pv)  { 
  uint64 alpha, beta, u, v; 
  u = 1; v = 0; 
  alpha = a; beta = b; // Note that alpha is 
  // even and beta is odd. 
 
  /* The invariant maintained from here on is: 
  a = u*2*alpha - v*beta. */ 

  while (a > 0) { 
    a = a >> 1; 
    if ((u & 1) == 0) { // Delete a common 
      u = u >> 1; v = v >> 1; // factor of 2 in 
    } // u and v. 
    else { 
      /* We want to set u = (u + beta) >> 1, but 
      that can overflow, so we use Dietz's method. */ 
      u = ((u ^ beta) >> 1) + (u & beta); 
      v = (v >> 1) + alpha; 
    }
  }
 
  *pu = u; 
  *pv = v; 
  return;
}



uint64 modul64(uint64 x, uint64 y, uint64 z) {


    /* Divides (x || y) by z, for 64-bit integers x, y,
    and z, giving the remainder (modulus) as the result.
    Must have x < z (to get a 64-bit result). This is
    checked for. */

    int64 i, t;

    if (x >= z) {
        printf("Bad call to modul64, must have x < z.");
        exit(1);
    }

    // can now use _udiv128 in VS19

    uint64_t r;
    uint64_t q = _udiv128(x, y, z, &r);
    return r;

    /*
    for (i = 1; i <= 64; i++) { // Do 64 times. 
        t = (int64)x >> 63; // All 1's if x(63) = 1. 
        x = (x << 1) | (y >> 63); // Shift x || y left 
        y = y << 1; // one bit. 
        if ((x | t) >= z) {
            x = x - z;
            y = y + 1;
        }
    }
    return x;
    */
}


void mulul64(uint64 u, uint64 v, uint64* whi, uint64* wlo) {
    *wlo = _umul128(u, v, &*whi);
    return;
    /*
  uint64 u0, u1, v0, v1, k, t;
  uint64 w0, w1, w2;

  u1 = u >> 32; u0 = u & 0xFFFFFFFF;
  v1 = v >> 32; v0 = v & 0xFFFFFFFF;

  t = u0*v0;
  w0 = t & 0xFFFFFFFF;
  k = t >> 32;

  t = u1*v0 + k;
  w1 = t & 0xFFFFFFFF;
  w2 = t >> 32;

  t = u0*v1 + w1;
  k = t >> 32;

  *wlo = (t << 32) + w0;
  *whi = u1*v1 + w2 + k;

  return;
  */
}


uint64 montmul(uint64 abar, uint64 bbar, uint64 m, uint64 mprime) { 
 
  uint64 thi, tlo, tm, tmmhi, tmmlo, uhi, ulo, ov; 
 
  mulul64(abar, bbar, &thi, &tlo); // t = abar*bbar. 

  /* Now compute u = (t + ((t*mprime) & mask)*m) >> 64. 
  The mask is fixed at 2**64-1. Because it is a 64-bit 
  quantity, it suffices to compute the low-order 64 
  bits of t*mprime, which means we can ignore thi. */ 
 
  tm = tlo*mprime; 
 
  mulul64(tm, m, &tmmhi, &tmmlo); // tmm = tm*m. 
 
  ulo = tlo + tmmlo; // Add t to tmm 
  uhi = thi + tmmhi; // (128-bit add). 
  if (ulo < tlo) uhi = uhi + 1; // Allow for a carry. 
 
  // The above addition can overflow. Detect that here. 
 
  ov = (uhi < thi) | ((uhi == thi) & (ulo < tlo)); 
 
  ulo = uhi; // Shift u right 
  uhi = 0; // 64 bit positions. 
 
  if (ov > 0 || ulo >= m) // If u >= m, 
  ulo = ulo - m; // subtract m from u. 
 
  return ulo; 
}


uint64 mulmodMont(uint64 baseM, uint64 e, uint64 modul, uint64 pv, uint64 oneM) {
    uint64 ans = oneM;
    while (e > 0) {
        if (e & 1) {
            ans = montmul(baseM, ans, modul, pv);
        }
        baseM = montmul(baseM, baseM, modul, pv);
        e >>= 1;
    }
    return ans;
}


bool is_SPRP(uint64 base, uint64 modul) {
    if (base >= modul) base = base % modul;
    uint64 pu, pv;
    xbinGCD(1ull << 63, modul, &pu, &pv);
    uint64 baseM = modul64(base, 0, modul);
    uint64 oneM = modul64(1, 0, modul);

    uint64 moneM = modul - oneM;
    uint64 e = modul - 1;
    while (e % 2 == 0) 
        e >>= 1;
    uint64 t = mulmodMont(baseM, e, modul, pv, oneM);
    if (t == oneM) return 1;
    while (e < modul - 1) {
        if (t == moneM) 
            return 1;
        t = montmul(t, t, modul, pv);
        if (t == oneM) 
            return 0;
        e <<= 1;
    }
    return 0;
}


int hashh(uint64 x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369;  // 0x3335b369
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3b;
    x = ((x >> 32) ^ x);
    return x&262143;
}


bool is_prime_2_64(uint64 a) {
    if (a==2 || a==3 || a==5 || a==7) return true;
    if (a%2==0 || a%3==0 || a%5==0 || a%7==0) return false;
    if (a<121) return (a>1);
     return is_SPRP(2ull,a) && is_SPRP(bases[hashh(a)],a);
}


//int main(void) {}
