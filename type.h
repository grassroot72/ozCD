/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */


#ifndef _OZ_TYPE_H_
#define _OZ_TYPE_H_

#include <smmintrin.h>


typedef unsigned int oz_uint_t;
typedef unsigned char oz_uchar_t;

typedef __m128  oz_sse4f_t;   /* vector of 4 floats */
typedef __m128i oz_sse4i_t;

typedef __m128d oz_sse2d;     /* vector of 2 doubles (sse2) */


#endif
