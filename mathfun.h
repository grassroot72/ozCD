/*
  mathfun.h

  Copyright (c) 2017 Nikos Tasios
  Copyright (C) 2019 by Edward LEI

  This code is licensed under the MIT license
*/

#ifndef _OZ_MATHFUN_H_
#define _OZ_MATHFUN_H_

#include <math.h>
#include <float.h>
#include "compiler.h"


#ifndef OZ_SINGLE
#  ifndef OZ_DOUBLE
#    error You must define OZ_SINGLE or OZ_DOUBLE
#  endif /* OZ_DOUBLE */
#endif /* OZ_SINGLE */

#ifdef OZ_SINGLE
#  ifdef OZ_DOUBLE
#    error You can define either OZ_SINGLE or OZ_DOUBLE, not both!
# endif /* OZ_DOUBLE */
#endif /* OZ_SINGLE */


#ifdef OZ_SINGLE
typedef float oz_real_t;

#  define OZ_EPS FLT_EPSILON
#  define OZ_REAL_MAX FLT_MAX

#  define OZ_REAL(x)        (x ## f)   /* form a constant */
#  define OZ_SQRT(x)        (sqrtf(x)) /* square root */
#  define OZ_FABS(x)        (fabsf(x)) /* absolute value */
#  define OZ_FMAX(x, y)     (fmaxf((x), (y))) /* maximum of two floats */
#  define OZ_FMIN(x, y)     (fminf((x), (y))) /* minimum of two floats */
#  define OZ_COPYSIGN(x, y) (copysignf((x), (y)))
#  define OZ_SIN(x)         (sinf(x))
#  define OZ_COS(x)         (cosf(x))
#endif  /* OZ_SINGLE */


#ifdef OZ_DOUBLE
typedef double oz_real_t;

#  define OZ_EPS DBL_EPSILON
#  define OZ_REAL_MAX DBL_MAX

#  define OZ_REAL(x)        (x)       /* form a constant */
#  define OZ_SQRT(x)        (sqrt(x)) /* square root */
#  define OZ_FABS(x)        (fabs(x)) /* absolute value */
#  define OZ_FMAX(x, y)     (fmax((x), (y))) /* maximum of two floats */
#  define OZ_FMIN(x, y)     (fmin((x), (y))) /* minimum of two floats */
#  define OZ_COPYSIGN(x, y) (copysign((x), (y)))
#  define OZ_SIN(x)         (sin(x))
#  define OZ_COS(x)         (cos(x))
#endif /* OZ_DOUBLE */


#define OZ_ONE OZ_REAL(1.)
#define OZ_ZERO OZ_REAL(0.)


#ifndef OZ_GJK_ERROR_TOLERANCE 
#define OZ_GJK_ERROR_TOLERANCE 10.0 * OZ_EPS
#endif

#ifndef OZ_GJK_RELATIVE_ERROR 
#define OZ_GJK_RELATIVE_ERROR 1.0e-4
#define OZ_GJK_RELATIVE_ERROR2 OZ_GJK_RELATIVE_ERROR * OZ_GJK_RELATIVE_ERROR
#endif


_oz_inline int
ozIsZero(oz_real_t val)
{
  return OZ_FABS(val) < OZ_EPS;
}

_oz_inline int
ozEq(const oz_real_t a,
     const oz_real_t b)
{
  oz_real_t ab;
  oz_real_t _a, _b;

  ab = OZ_FABS(a - b);
  if (OZ_FABS(ab) < OZ_EPS) return 1;

  _a = OZ_FABS(a);
  _b = OZ_FABS(b);
  if (_b > _a) {
    return ab < OZ_EPS * _b;
  }
  else {
    return ab < OZ_EPS * _a;
  }
}


#endif