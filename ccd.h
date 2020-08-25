/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#ifndef _OZ_CCD_H_
#define _OZ_CCD_H_


/*
 * Check if two convex shapes ca and cb overlap
 * with their respective dynamic transform data ta and tb
 *  Return value:
 *  1  <- overlap
 *  0  <- not overlap
 */
int ozGJKOverlap (
  const oz_transform_t *ta, const void *ca,
  const oz_transform_t *tb, const void *cb
);

/*
 * Calculate the closest distance between convex shape ca and cb
 * with their respective dynamic transform data ta and tb
 */
void ozGJKDist (
  oz_real_t *dist,
  const oz_transform_t *ta, const void *ca,
  const oz_transform_t *tb, const void *cb
);

/*
 * Calculate the closest points,
 * p_on_a on convex shape ca, p_on_b on convex shape cb,
 * with their respective dynamic transform data ta and tb
 * If p_on_a and p_on_b are zeroed, that means, shape ca and cb overlap
 */
void ozGJKClosestPoints (
  oz_real_t *p_on_a, oz_real_t *p_on_b,
  const oz_transform_t *ta, const void *ca,
  const oz_transform_t *tb, const void *cb
);

/*
 * Calculate the distance between convex
 * shape ca with dynamic transform data ta and
 * shape cb with dynamic transform data tb
 * in a certain direction
 *  Return value:
 *  0  <- ca does not meet cb
 *  1  <- ca meets cb
 */
int ozGJKRaycast (
  oz_real_t *dist, oz_real_t *normal,
  const oz_transform_t *ta, const void *ca,
  const oz_transform_t *tb, const void *cb,
  const oz_real_t *raydir
);


#endif
