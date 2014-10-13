//
//  common.h
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__common__
#define __Renderer__common__

#include <immintrin.h>
#include <cstdarg>
#include <math.h>
#include <random>

#define HIPRECISION //whether to use floats or doubles
#define EPSILON 0.00000001

# ifdef HIPRECISION
#   define real double
#   define realsize long
#   define realbytes 8
#   define _abs fabs
#   define _acos acos
#   define _asin asin
#   define _atan atan
#   define _atan2 atan2
#   define _cos cos
#   define _exp exp
#   define _INFINITY HUGE_VAL
#   define _max fmax
#   define _min fmin
#   define _nan nan
#   define _pow pow
#   define _sin sin
#   define _sqrt sqrt
#   define _tan tan
# else
#   define real float
#   define realsize int
#   define realbytes 4
#   define _abs fabsf
#   define _acos acosf
#   define _asin asinf
#   define _atan atanf
#   define _atan2 atan2f
#   define _cos cosf
#   define _exp expf
#   define _INFINITY HUGE_VALF
#   define _max fmaxf
#   define _min fminf
#   define _nan nanf
#   define _pow powf
#   define _sin sinf
#   define _sqrt sqrtf
#   define _tan tanf
# endif

/* MATH UTILITY FUNCTIONS */

//these clamp the passed value 'a' to the range of min, max
inline real clamp(const real a, const real min, const real max){
    return a>max?max:(a<min?min:a);
}
inline void clamp(real * const a, const real min, const real max){
    *a = *a>max ? max : (*a<min?min:*a);
}

//Linearly interpolates between 'a' and 'b' using the ratio 'r'.
// r = 0 returns a, and r = 1 returns b.
inline real lerp(const real a, const real b, const real r){
    return a+r*(b-a);
}

//returns 1 if a>0,
//        0 if a==0, and
//       -1 is a<0
template <typename T>
inline T signum(T a){
    return ((a>0) - (0>a));
}

//vector functions
const __v4si zeroveci = __v4si{0,0,0,0};
const __v4sf zerovecf = __v4sf{0,0,0,0};

//multiplies two integer vectors 
inline __v4si muli32(const __v4si &a, const __v4si &b)
{
#ifdef __SSE4_1__  // modern CPU - use SSE 4.1
    return _mm_mullo_epi32(a, b);
#else               // old CPU - use SSE 2
    __v4si tmp1 = _mm_mul_epu32(a,b); /* mul 2,0*/
    __v4si tmp2 = _mm_mul_epu32( _mm_srli_si128(a,4), _mm_srli_si128(b,4)); /* mul 3,1 */
    return _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE (0,0,2,0)), _mm_shuffle_epi32(tmp2, _MM_SHUFFLE (0,0,2,0))); /* shuffle results to [63..0] and pack */
#endif
}

#endif
