/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#ifndef _OZ_QUAT_H_
#define _OZ_QUAT_H_

#include "vec3.h"


_oz_inline void
ozQuatInverse(oz_real_t *r,
              const oz_real_t *q)
{
  oz_real_t ilength2 = OZ_REAL(1.0)  / (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  r[0] = -q[0] * ilength2;
  r[1] = -q[1] * ilength2;
  r[2] = -q[2] * ilength2;
  r[3] =  q[3] * ilength2;
}

_oz_inline void
ozQuatVec3Rotate(oz_real_t *r,
                 const oz_real_t *q,
                 const oz_real_t *v)
{
  oz_real_t u[3];

  oz_real_t a[3], b[3];
  ozVec3Cross(a, q, v);
  ozVec3FMAdd(b, q[3], v, a);
  ozVec3Cross(u, q, b);

  ozVec3FMAdd(r, OZ_REAL(2.0), u, v);
}


#endif
