/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#ifndef _OZ_VEC3_H_
#define _OZ_VEC3_H_

#include "mathfun.h"


_oz_inline oz_real_t
ozVec3Dot(const oz_real_t *a,
          const oz_real_t *b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

_oz_inline void
ozVec3Cross(oz_real_t *c,
            const oz_real_t *a,
            const oz_real_t *b)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/* Calculates d = (a x b) x c */
_oz_inline void
ozVec3TripleProduct(oz_real_t *d,
                    const oz_real_t *a,
                    const oz_real_t *b,
                    const oz_real_t *c)
{
  d[0] = b[0] * ozVec3Dot(a, c) - a[0] * ozVec3Dot(b, c);
  d[1] = b[1] * ozVec3Dot(a, c) - a[1] * ozVec3Dot(b, c);
  d[2] = b[2] * ozVec3Dot(a, c) - a[2] * ozVec3Dot(b, c);
}

_oz_inline void
ozVec3Add(oz_real_t *c,
          const oz_real_t *a,
          const oz_real_t *b)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

_oz_inline void
ozVec3Sub(oz_real_t *c,
          const oz_real_t *a,
          const oz_real_t *b)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

_oz_inline void
ozVec3SMul(oz_real_t *c,
           const oz_real_t f,
           const oz_real_t *a)
{
  c[0] = f * a[0];
  c[1] = f * a[1];
  c[2] = f * a[2];
}

_oz_inline void
ozVec3FMAdd(oz_real_t *c,
            const oz_real_t f,
            const oz_real_t *a,
            const oz_real_t *b)
{
  c[0] = f * a[0] + b[0];
  c[1] = f * a[1] + b[1];
  c[2] = f * a[2] + b[2];
}

_oz_inline oz_real_t
ozVec3Len2(const oz_real_t *vec)
{
  return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

_oz_inline oz_real_t
ozVec3Len(const oz_real_t *vec)
{
  return OZ_SQRT(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

_oz_inline int
ozVec3Eq(const oz_real_t *a,
         const oz_real_t *b)
{
  return ozEq(a[0], b[0]) && ozEq(a[1], b[1]) && ozEq(a[2], b[2]);
}


#endif
