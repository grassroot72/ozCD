
#ifndef _OZ_SSE_MATHFUN_H_
#define _OZ_SSE_MATHFUN_H_

#include <smmintrin.h>
#include "type.h"


# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))


/* SSE2 specific */
#define _MM_SET_PS(d,c,b,a)  _mm_setr_ps((a),(b),(c),(d))
#define _MM_SET1_PS(a)       _mm_set1_ps((a))
#define _MM_SETR_PS(a,b,c,d) _mm_setr_ps((a),(b),(c),(d))
#define _MM_SETZERO_PS()     _mm_setzero_ps()
#define _MM_SETZERO_SI128()  _mm_setzero_si128()
#define _MM_SETZERO_SI64()   _mm_setzero_si64()

#define VELEM(type,a,n)      (((type*)&a)[n])


/* declare some SSE constants -- why can't I figure a better way to do that? */
#define _PS_CONST(Name, Val)                                            \
  static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { (const float)(Val), (const float)(Val), (const float)(Val), (const float)(Val) }
#define _PI32_CONST(Name, Val)                                          \
  static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }



_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS_CONST_TYPE(sign_mask, int, 0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS_CONST(cephes_log_p0, 7.0376836292E-2);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS_CONST(cephes_log_p2, 1.1676998740E-1);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS_CONST(cephes_log_q1, -2.12194440e-4);
_PS_CONST(cephes_log_q2, 0.693359375);


/*
  natural logarithm computed for 4 simultaneous float
  return NaN for x <= 0
*/
static inline oz_sse4f_t log_ps(oz_sse4f_t x)
{
  oz_sse4f_t e;
  oz_sse4i_t emm0;

  oz_sse4f_t one = *(oz_sse4f_t *)_ps_1;
  oz_sse4f_t invalid_mask = _mm_cmple_ps(x, _MM_SETZERO_PS());

  /* cut off denormalized stuff */
  x = _mm_max_ps(x, *(oz_sse4f_t *)_ps_min_norm_pos);

  emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);

  /* keep only the fractional part */
  x = _mm_and_ps(x, *(oz_sse4f_t *)_ps_inv_mant_mask);
  x = _mm_or_ps(x, *(oz_sse4f_t *)_ps_0p5);

  emm0 = _mm_sub_epi32(emm0, *(oz_sse4i_t *)_pi32_0x7f);
  e = _mm_cvtepi32_ps(emm0);

  e = _mm_add_ps(e, one);

  /* part2:
    if ( x < SQRTHF ) {
      e -= 1;
      x = x + x - 1.0;
    } else { x = x - 1.0; }
  */

  oz_sse4f_t z, y;
  oz_sse4f_t mask = _mm_cmplt_ps(x, *(oz_sse4f_t *)_ps_cephes_SQRTHF);
  oz_sse4f_t tmp = _mm_and_ps(x, mask);
  x = _mm_sub_ps(x, one);
  e = _mm_sub_ps(e, _mm_and_ps(one, mask));
  x = _mm_add_ps(x, tmp);


  z = _mm_mul_ps(x,x);

  y = *(oz_sse4f_t *)_ps_cephes_log_p0;
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p1);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p2);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p3);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p4);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p5);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p6);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p7);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_log_p8);
  y = _mm_mul_ps(y, x);

  y = _mm_mul_ps(y, z);

  tmp = _mm_mul_ps(e, *(oz_sse4f_t *)_ps_cephes_log_q1);
  y = _mm_add_ps(y, tmp);


  tmp = _mm_mul_ps(z, *(oz_sse4f_t *)_ps_0p5);
  y = _mm_sub_ps(y, tmp);

  tmp = _mm_mul_ps(e, *(oz_sse4f_t *)_ps_cephes_log_q2);
  x = _mm_add_ps(x, y);
  x = _mm_add_ps(x, tmp);
  x = _mm_or_ps(x, invalid_mask); /* negative arg will be NAN */

  return x;
}

_PS_CONST(exp_hi,  88.3762626647949f);
_PS_CONST(exp_lo, -88.3762626647949f);

_PS_CONST(cephes_LOG2EF,  1.44269504088896341);
_PS_CONST(cephes_exp_C1,  0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1);

/*
  computes e**x of the 4 floats in x
*/
static inline oz_sse4f_t exp_ps(oz_sse4f_t x)
{
  oz_sse4f_t tmp = _MM_SETZERO_PS(), fx, mask, y, z;
  oz_sse4f_t pow2n;

  oz_sse4i_t emm0;
  oz_sse4f_t one = *(oz_sse4f_t *)_ps_1;

  x = _mm_min_ps(x, *(oz_sse4f_t *)_ps_exp_hi);
  x = _mm_max_ps(x, *(oz_sse4f_t *)_ps_exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm_mul_ps(x, *(oz_sse4f_t *)_ps_cephes_LOG2EF);
  fx = _mm_add_ps(fx, *(oz_sse4f_t *)_ps_0p5);

  /* how to perform a floorf with SSE: just below */
  emm0 = _mm_cvttps_epi32(fx);
  tmp  = _mm_cvtepi32_ps(emm0);

  /* if greater, substract 1 */
  mask = _mm_cmpgt_ps(tmp, fx);
  mask = _mm_and_ps(mask, one);
  fx = _mm_sub_ps(tmp, mask);

  tmp = _mm_mul_ps(fx, *(oz_sse4f_t *)_ps_cephes_exp_C1);
  z = _mm_mul_ps(fx, *(oz_sse4f_t *)_ps_cephes_exp_C2);
  x = _mm_sub_ps(x, tmp);
  x = _mm_sub_ps(x, z);

  z = _mm_mul_ps(x,x);

  y = *(oz_sse4f_t *)_ps_cephes_exp_p0;
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_exp_p1);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_exp_p2);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_exp_p3);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_exp_p4);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_cephes_exp_p5);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, x);
  y = _mm_add_ps(y, one);

  /* build 2^n */
  emm0 = _mm_cvttps_epi32(fx);
  emm0 = _mm_add_epi32(emm0, *(oz_sse4i_t *)_pi32_0x7f);
  emm0 = _mm_slli_epi32(emm0, 23);
  pow2n = _mm_castsi128_ps(emm0);

  y = _mm_mul_ps(y, pow2n);
  return y;
}


_PS_CONST(minus_cephes_DP1, -0.78515625);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8);
_PS_CONST(sincof_p0, -1.9515295891E-4);
_PS_CONST(sincof_p1,  8.3321608736E-3);
_PS_CONST(sincof_p2, -1.6666654611E-1);
_PS_CONST(coscof_p0,  2.443315711809948E-005);
_PS_CONST(coscof_p1, -1.388731625493765E-003);
_PS_CONST(coscof_p2,  4.166664568298827E-002);
_PS_CONST(cephes_FOPI, 1.27323954473516); /* 4 / M_PI */


static inline oz_sse4f_t sin_ps(oz_sse4f_t x)
{
  oz_sse4f_t xmm1, xmm2 = _MM_SETZERO_PS(), xmm3, sign_bit, y, y2, z, tmp;

  oz_sse4f_t swap_sign_bit, poly_mask;
  oz_sse4i_t emm0, emm2;

  sign_bit = x;
  /* take the absolute value */
  x = _mm_and_ps(x, *(oz_sse4f_t *)_ps_inv_sign_mask);
  /* extract the sign bit (upper one) */
  sign_bit = _mm_and_ps(sign_bit, *(oz_sse4f_t *)_ps_sign_mask);

  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(oz_sse4f_t *)_ps_cephes_FOPI);

  /* store the integer part of y in mm0 */
  emm2 = _mm_cvttps_epi32(y);
  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, *(oz_sse4i_t *)_pi32_1);
  emm2 = _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_inv1);
  y = _mm_cvtepi32_ps(emm2);
  /* get the swap sign flag */
  emm0 = _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_4);
  emm0 = _mm_slli_epi32(emm0, 29);
  /*
    get the polynom selection mask
    there is one polynom for 0 <= x <= Pi/4
    and another one for Pi/4<x<=Pi/2

    Both branches will be computed.
  */
  emm2 = _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_2);
  emm2 = _mm_cmpeq_epi32(emm2, _MM_SETZERO_SI128());

  swap_sign_bit = _mm_castsi128_ps(emm0);
  poly_mask = _mm_castsi128_ps(emm2);
  sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);


  /*
    The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3;
  */
  xmm1 = *(oz_sse4f_t *)_ps_minus_cephes_DP1;
  xmm2 = *(oz_sse4f_t *)_ps_minus_cephes_DP2;
  xmm3 = *(oz_sse4f_t *)_ps_minus_cephes_DP3;
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  y = *(oz_sse4f_t *)_ps_coscof_p0;
  z = _mm_mul_ps(x,x);

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  tmp = _mm_mul_ps(z, *(oz_sse4f_t *)_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
  y2 = *(oz_sse4f_t *)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(oz_sse4f_t *)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(oz_sse4f_t *)_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  y2 = _mm_and_ps(xmm3, y2); //, xmm3);
  y = _mm_andnot_ps(xmm3, y);
  y = _mm_add_ps(y, y2);
  /* update the sign */
  y = _mm_xor_ps(y, sign_bit);

  return y;
}

/*
  almost the same as sin_ps
*/
static inline oz_sse4f_t cos_ps(oz_sse4f_t x)
{
  oz_sse4f_t xmm1, xmm2 = _MM_SETZERO_PS(), xmm3, y, y2, z, sign_bit, poly_mask, tmp;
  oz_sse4i_t emm0, emm2;

  /* take the absolute value */
  x = _mm_and_ps(x, *(oz_sse4f_t *)_ps_inv_sign_mask);

  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(oz_sse4f_t *)_ps_cephes_FOPI);

  /* store the integer part of y in mm0 */
  emm2 = _mm_cvttps_epi32(y);
  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, *(oz_sse4i_t *)_pi32_1);
  emm2 = _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_inv1);
  y = _mm_cvtepi32_ps(emm2);

  emm2 = _mm_sub_epi32(emm2, *(oz_sse4i_t *)_pi32_2);

  /* get the swap sign flag */
  emm0 = _mm_andnot_si128(emm2, *(oz_sse4i_t *)_pi32_4);
  emm0 = _mm_slli_epi32(emm0, 29);
  /* get the polynom selection mask */
  emm2 = _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_2);
  emm2 = _mm_cmpeq_epi32(emm2, _MM_SETZERO_SI128());

  sign_bit = _mm_castsi128_ps(emm0);
  poly_mask = _mm_castsi128_ps(emm2);

  /*
    The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3;
  */
  xmm1 = *(oz_sse4f_t *)_ps_minus_cephes_DP1;
  xmm2 = *(oz_sse4f_t *)_ps_minus_cephes_DP2;
  xmm3 = *(oz_sse4f_t *)_ps_minus_cephes_DP3;
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  y = *(oz_sse4f_t *)_ps_coscof_p0;
  z = _mm_mul_ps(x,x);

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  tmp = _mm_mul_ps(z, *(oz_sse4f_t *)_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, *(oz_sse4f_t *)_ps_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
  y2 = *(oz_sse4f_t *)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(oz_sse4f_t *)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(oz_sse4f_t *)_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  y2 = _mm_and_ps(xmm3, y2); //, xmm3);
  y = _mm_andnot_ps(xmm3, y);
  y = _mm_add_ps(y, y2);
  /* update the sign */
  y = _mm_xor_ps(y, sign_bit);

  return y;
}

/*
  since sin_ps and cos_ps are almost identical, 
  sincos_ps could replace both of them..
  it is almost as fast, and gives you a free cosine with your sine
*/
static inline void sincos_ps(oz_sse4f_t x, oz_sse4f_t *s, oz_sse4f_t *c)
{
  oz_sse4f_t xmm1, xmm2, sign_bit_sin, y, y2, z, swap_sign_bit_sin, poly_mask;
  oz_sse4f_t sign_bit_cos;

  oz_sse4i_t emm2;

  sign_bit_sin = x;
  /* take the absolute value */
  x = _mm_and_ps(x, *(oz_sse4f_t *)_ps_inv_sign_mask);
  /* extract the sign bit (upper one) */
  sign_bit_sin = _mm_and_ps(sign_bit_sin, *(oz_sse4f_t *)_ps_sign_mask);

  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(oz_sse4f_t *)_ps_cephes_FOPI);

  /* store the integer part of y in emm2 */
  emm2 = _mm_cvttps_epi32(y);

  emm2 = _mm_and_si128( _mm_add_epi32( _mm_cvttps_epi32(y), *(oz_sse4i_t *)_pi32_1 ), *(oz_sse4i_t *)_pi32_inv1 );
  y = _mm_cvtepi32_ps(emm2);

  /* get the swap sign flag for the sine */
  swap_sign_bit_sin = _mm_castsi128_ps( _mm_slli_epi32( _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_4), 29) );

  /* get the polynom selection mask for the sine*/
  poly_mask = _mm_castsi128_ps( _mm_cmpeq_epi32( _mm_and_si128(emm2, *(oz_sse4i_t *)_pi32_2), _MM_SETZERO_SI128()) );


  /*
    The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3;
  */
  x += y * ( *(oz_sse4f_t *)_ps_minus_cephes_DP1 + *(oz_sse4f_t *)_ps_minus_cephes_DP2 + *(oz_sse4f_t *)_ps_minus_cephes_DP3 );

  sign_bit_cos = _mm_castsi128_ps( _mm_slli_epi32( _mm_andnot_si128( _mm_sub_epi32(emm2, *(oz_sse4i_t *)_pi32_2), *(oz_sse4i_t *)_pi32_4), 29) );
  sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  z = x * x;
  y = ( ( ( (*(oz_sse4f_t *)_ps_coscof_p0) * z + *(oz_sse4f_t *)_ps_coscof_p1 ) * z + *(oz_sse4f_t *)_ps_coscof_p2 ) * z
  		- *(oz_sse4f_t *)_ps_0p5 ) * z + *(oz_sse4f_t *)_ps_1;

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
  y2 = ( ( ( ( ((*(oz_sse4f_t *)_ps_sincof_p0) * z ) + *(oz_sse4f_t *)_ps_sincof_p1 ) * z ) + *(oz_sse4f_t *)_ps_sincof_p2 ) * z
      + *(oz_sse4f_t *)_ps_1 ) * x;

  /* select the correct result from the two polynoms */
  xmm1 = _mm_add_ps( _mm_andnot_ps( poly_mask, y), _mm_and_ps(poly_mask, y2) );
  xmm2 = _mm_sub_ps( _mm_add_ps( y, y2 ), xmm1 );
  /* update the sign */
  *s = _mm_xor_ps(xmm1, sign_bit_sin);
  *c = _mm_xor_ps(xmm2, sign_bit_cos);
}


#endif
