/*
  ccd.c

  Copyright (c) 2017 Nikos Tasios
  Copyright (C) 2019 by Edward LEI

  This code is licensed under the MIT license
*/

#include "vec3.h"
#include "quat.h"
#include "apex_memmove.h"
#include "simplex.h"
#include "shape.h"
#include "ccd.h"


int
ozGJKOverlap(const oz_transform_t *ta,
             const void *ca,
             const oz_transform_t *tb,
             const void *cb)
{
  oz_real_t dir[3];
  ozVec3Sub(dir, tb->pos, ta->pos);

  oz_simplex_t simplex;
  ozSimplexInit(&simplex);

  oz_uint_t fail_safe = 0;

  oz_real_t inv_rot_a[4], inv_rot_b[4];
  ozQuatInverse(inv_rot_a, ta->rot);
  ozQuatInverse(inv_rot_b, tb->rot);

  ozpfnSupport sa = *(const ozpfnSupport *)ca;
  ozpfnSupport sb = *(const ozpfnSupport *)cb;

  do {
    oz_real_t vertex_a[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozQuatVec3Rotate(inv_dir, inv_rot_a, dir);

      sa(support_point, ca, inv_dir);

      ozQuatVec3Rotate(vertex_a, ta->rot, support_point);
      ozVec3FMAdd(vertex_a, ta->size, vertex_a, ta->pos);
    }

    oz_real_t vertex_b[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozVec3SMul(inv_dir, -OZ_REAL(1.0), dir);
      ozQuatVec3Rotate(inv_dir, inv_rot_b, inv_dir);

      sb(support_point, cb, inv_dir);

      ozQuatVec3Rotate(vertex_b, tb->rot, support_point);
      ozVec3FMAdd(vertex_b, tb->size, vertex_b, tb->pos);
    }

    oz_real_t new_point[3];
    ozVec3Sub(new_point, vertex_a, vertex_b);

    oz_real_t dn = ozVec3Dot(dir, new_point);

    if (dn < OZ_REAL(0.0) || ozSimplexContains(&simplex, new_point)) {
      return 0;
    }
    ozSimplexAddPoint(&simplex, new_point);

    if (ozSimplexContainsOrigin(&simplex, dir) || ozVec3Len2(dir) == OZ_REAL(0.0)) {
      return 1;
    }
  } while (fail_safe++ < 2000);

  return 1;
}

void
ozGJKDist(oz_real_t *dist,
          const oz_transform_t *ta,
          const void *ca,
          const oz_transform_t *tb,
          const void *cb)
{
  oz_real_t inv_rot_a[4], inv_rot_b[4];
  ozQuatInverse(inv_rot_a, ta->rot);
  ozQuatInverse(inv_rot_b, tb->rot);

  oz_real_t dir[3] = {OZ_REAL(0.0)};
  oz_simplex_t simplex;
  ozSimplexInit(&simplex);

  oz_uint_t fail_safe = 0;

  oz_real_t dist2 = OZ_REAL_MAX;

  ozpfnSupport sa = *(const ozpfnSupport *)ca;
  ozpfnSupport sb = *(const ozpfnSupport *)cb;

  do {
    oz_real_t vertex_a[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozVec3SMul(inv_dir, -OZ_REAL(1.0), dir);
      ozQuatVec3Rotate(inv_dir, inv_rot_a, inv_dir);

      sa(support_point, ca, inv_dir);

      ozQuatVec3Rotate(vertex_a, ta->rot, support_point);
      ozVec3FMAdd(vertex_a, ta->size, vertex_a, ta->pos);
    }

    oz_real_t vertex_b[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozQuatVec3Rotate(inv_dir, inv_rot_b, dir);

      sb(support_point, cb, inv_dir);

      ozQuatVec3Rotate(vertex_b, tb->rot, support_point);
      ozVec3FMAdd(vertex_b, tb->size, vertex_b, tb->pos);
    }

    oz_real_t new_point[3];
    ozVec3Sub(new_point, vertex_a, vertex_b);

    if (ozSimplexContains(&simplex, new_point) ||
        dist2 - ozVec3Dot(dir, new_point) <= dist2 * OZ_GJK_RELATIVE_ERROR2) {
      apex_memcpy(dist, dir, 3 * sizeof(*dist));
      return;
    }
    ozSimplexAddPoint(&simplex, new_point);

    ozSimplexClosest(&simplex, dir);

    dist2 = ozVec3Len2(dir);

    if (simplex._size == 4 ||
        dist2 < OZ_GJK_ERROR_TOLERANCE * simplex._max_vert2) {
      memset(dist, 0, 3 * sizeof(*dist));
      return;
    }

  } while (++fail_safe < 200);

  apex_memcpy(dist, dir, 3 * sizeof(*dist));
}

void
ozGJKClosestPoints(oz_real_t *p_on_a,
                   oz_real_t *p_on_b,
                   const oz_transform_t *ta,
                   const void *ca,
                   const oz_transform_t *tb,
                   const void *cb)
{
  oz_real_t inv_rot_a[4], inv_rot_b[4];
  ozQuatInverse(inv_rot_a, ta->rot);
  ozQuatInverse(inv_rot_b, tb->rot);

  oz_real_t dir[3] = {OZ_REAL(0.0)};
  oz_simplex_t simplex;
  ozSimplexInit(&simplex);

  oz_uint_t fail_safe = 0;

  oz_real_t dist2 = OZ_REAL_MAX;

  ozpfnSupport sa = *(const ozpfnSupport *)ca;
  ozpfnSupport sb = *(const ozpfnSupport *)cb;

  do {
    oz_real_t vertex_a[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozVec3SMul(inv_dir, -OZ_REAL(1.0), dir);
      ozQuatVec3Rotate(inv_dir, inv_rot_a, inv_dir);

      sa(support_point, ca, inv_dir);

      ozQuatVec3Rotate(vertex_a, ta->rot, support_point);
      ozVec3FMAdd(vertex_a, ta->size, vertex_a, ta->pos);
    }

    oz_real_t vertex_b[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozQuatVec3Rotate(inv_dir, inv_rot_b, dir);

      sb(support_point, cb, inv_dir);

      ozQuatVec3Rotate(vertex_b, tb->rot, support_point);
      ozVec3FMAdd(vertex_b, tb->size, vertex_b, tb->pos);
    }

    oz_real_t new_point[3];
    ozVec3Sub(new_point, vertex_a, vertex_b);

    if (ozSimplexContains(&simplex, new_point) ||
        dist2 - ozVec3Dot(dir, new_point) <= dist2 * OZ_GJK_RELATIVE_ERROR2) {
      ozSimplexComputeClosestPoints(&simplex, p_on_a, p_on_b, dir);
      return;
    }

    ozSimplexAddPointWithInfo(&simplex, new_point, vertex_a, vertex_b);

    ozSimplexClosest(&simplex, dir);

    dist2 = ozVec3Len2(dir);

    if (simplex._size == 4 ||
        dist2 < OZ_GJK_ERROR_TOLERANCE * simplex._max_vert2) {
      memset(p_on_a, 0, 3 * sizeof(*p_on_a));
      memset(p_on_b, 0, 3 * sizeof(*p_on_b));
      return;
    }

  } while (++fail_safe < 2000);

  ozSimplexComputeClosestPoints(&simplex, p_on_a, p_on_b, dir);
}

int
ozGJKRaycast(oz_real_t *dist,
             oz_real_t *normal,
             const oz_transform_t *ta,
             const void *ca,
             const oz_transform_t *tb,
             const void *cb,
             const oz_real_t *raydir)
{
  oz_real_t inv_rot_a[4], inv_rot_b[4];
  ozQuatInverse(inv_rot_a, ta->rot);
  ozQuatInverse(inv_rot_b, tb->rot);

  oz_real_t dir[3];
  ozVec3Sub(dir, tb->pos, ta->pos);
  oz_simplex_t simplex;
  ozSimplexInit(&simplex);

  oz_uint_t fail_safe = 0;

  oz_real_t x[3] = {OZ_REAL(0.0)};
  oz_real_t lambda = OZ_REAL(0.0);

  oz_real_t dist2 = OZ_REAL_MAX;

  ozpfnSupport sa = *(const ozpfnSupport *)ca;
  ozpfnSupport sb = *(const ozpfnSupport *)cb;

  do {
    oz_real_t vertex_a[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozVec3SMul(inv_dir, -1.0, dir);
      ozQuatVec3Rotate(inv_dir, inv_rot_a, inv_dir);

      sa(support_point, ca, inv_dir);

      ozQuatVec3Rotate(vertex_a, ta->rot, support_point);
      ozVec3FMAdd(vertex_a, ta->size, vertex_a, ta->pos);
    }

    oz_real_t vertex_b[3];
    {
      oz_real_t inv_dir[3], support_point[3];
      ozQuatVec3Rotate(inv_dir, inv_rot_b, dir);

      sb(support_point, cb, inv_dir);

      ozQuatVec3Rotate(vertex_b, tb->rot, support_point);
      ozVec3FMAdd(vertex_b, tb->size, vertex_b, tb->pos);
    }
    oz_real_t new_point[3];
    ozVec3Sub(new_point, vertex_a, vertex_b);

    oz_real_t new_point_trans[3];
    ozVec3Sub(new_point_trans, new_point, x);

    if (ozVec3Dot(dir, new_point_trans) > OZ_REAL(0.0)) {
      if (ozVec3Dot(dir, raydir) >= OZ_REAL(0.0)) {
        return 0;
      }

      oz_real_t delta = ozVec3Dot(dir, new_point_trans) / ozVec3Dot(dir, raydir);
      lambda -= delta;
      if (lambda > *dist) {
        return 0;
      }

      ozVec3SMul(x, -lambda, raydir);
      ozVec3SMul(normal, -OZ_REAL(1.0) / ozVec3Len(dir), dir);
      oz_real_t dr[3];
      ozVec3SMul(dr, -delta, raydir);
      ozSimplexTranslate(&simplex, dr);
    }

    ozSimplexAddPoint(&simplex, new_point_trans);
    ozSimplexClosest(&simplex, dir);

    dist2 = ozVec3Len2(dir);

    if (simplex._size == 4 ||
        dist2 < OZ_GJK_ERROR_TOLERANCE * simplex._max_vert2) {
      *dist = lambda;
      return 1;
    }
  } while (fail_safe++ < 1000);

  *dist = lambda;

  return 0;
}
