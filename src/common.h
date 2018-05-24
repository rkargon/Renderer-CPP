//
//  common.h
//  Renderer
//
//  Utility functions and definitions.
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__common__
#define __Renderer__common__

#define USE_VECTOR

#ifdef USE_VECTOR
#include <immintrin.h>
#endif

#include <cstdarg>
#include <math.h>
#include <random>
#include <utility>

#define EPSILON 0.00000001

/* MATH UTILITY FUNCTIONS */

inline bool eps_equal(double a, double b) { return std::abs(a - b) < EPSILON; }
inline bool eps_zero(double a) { return eps_equal(a, 0); }

// Linearly interpolates between 'a' and 'b' using the ratio 'r'.
// r = 0 returns a, and r = 1 returns b.
inline double lerp(const double a, const double b, const double r) {
  return a + r * (b - a);
}

// interpolate between three values p1,2,3 using barycentric coordinates w1,2,3
inline double barycentric_lerp(const double p1, const double p2,
                               const double p3, const double w1,
                               const double w2, const double w3) {
  return p1 * w1 + p2 * w2 + p3 * w3;
}

// returns 1 if a>0,
//        0 if a==0, and
//       -1 is a<0
template <typename T> inline T signum(T a) { return ((a > 0) - (0 > a)); }

// vector functions
#ifdef USE_VECTOR
const __v4si zeroveci = __v4si{0, 0, 0, 0};
const __v4sf zerovecf = __v4sf{0, 0, 0, 0};

// multiplies two integer vectors
inline __v4si muli32(const __v4si &a, const __v4si &b) {
#ifdef __SSE4_1__ // modern CPU - use SSE 4.1
  return _mm_mullo_epi32(a, b);
#else // old CPU - use SSE 2
  __v4si tmp1 = _mm_mul_epu32(a, b); /* mul 2,0*/
  __v4si tmp2 =
      _mm_mul_epu32(_mm_srli_si128(a, 4), _mm_srli_si128(b, 4)); /* mul 3,1 */
  return _mm_unpacklo_epi32(
      _mm_shuffle_epi32(tmp1, _MM_SHUFFLE(0, 0, 2, 0)),
      _mm_shuffle_epi32(
          tmp2,
          _MM_SHUFFLE(0, 0, 2, 0))); /* shuffle results to [63..0] and pack */
#endif
}
#endif // USE_VECTOR

template <class T> inline void hash_combine(std::size_t &seed, const T &v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

#endif
