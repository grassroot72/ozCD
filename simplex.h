/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#ifndef _OZ_SIMPLEX_H_
#define _OZ_SIMPLEX_H_

#include <string.h>
#include "compiler.h"
#include "mathfun.h"


struct _oz_simplex_s {
  oz_uchar_t _bits;
  oz_uchar_t _last_sb;
  oz_uchar_t _size;
  oz_real_t _p[3 * 4];
  oz_real_t _a[3 * 4];
  oz_real_t _b[3 * 4];
  oz_real_t _max_vert2;
};
typedef struct _oz_simplex_s oz_simplex_t;


_oz_inline void
ozSimplexInit(oz_simplex_t *simplex)
{
  simplex->_bits = 0;
  simplex->_last_sb = 0;
  simplex->_size = 0;
  simplex->_max_vert2 = OZ_REAL(0.0);
}

_oz_inline void
ozSimplexRemovePoint(oz_simplex_t *simplex,
                     int p)
{
  simplex->_bits ^= (1 << p); /* Erase the bit at position p */
  --simplex->_size;
}


void
ozSimplexAddPoint(oz_simplex_t *simplex,
                  const oz_real_t *point);

void
ozSimplexAddPointWithInfo(oz_simplex_t *simplex,
                          const oz_real_t *point,
                          const oz_real_t *pa,
                          const oz_real_t *pb);

int
ozSimplexContains(const oz_simplex_t *simplex,
                  const oz_real_t *point);

void
ozSimplexTranslate(oz_simplex_t *simplex,
                   const oz_real_t *dr);

void
ozSimplexComputeClosestPoints(oz_simplex_t *simplex,
                              oz_real_t *pa,
                              oz_real_t *pb,
                              const oz_real_t *P);

void
ozSimplexClosest(oz_simplex_t *simplex,
                 oz_real_t *dir);

int
ozSimplexContainsOrigin(oz_simplex_t *simplex,
                        oz_real_t *dir);


#endif
