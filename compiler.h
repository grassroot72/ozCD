/*
  compiler.h

  Copyright (c) 2017 Nikos Tasios
  Copyright (C) 2019 by Edward LEI

  This code is licensed under the MIT license
*/

#ifndef _OZ_COMPILER_H_
#define _OZ_COMPILER_H_


#include <stddef.h>

/* define oz inline keyword */
#ifdef __GNUC__
# define _oz_inline static inline __attribute__((always_inline))
#else /* __GNUC__ */
# define _oz_inline static __inline
#endif /* __GNUC__ */


typedef unsigned int oz_uint_t;
typedef unsigned char oz_uchar_t;


#endif