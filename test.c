/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#include <stdio.h>
#include "shape.h"
#include "ccd.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

int main(int argc, char* argv[])
{
  oz_cylinder_t cyl;
  ozCylinderInit(&cyl, OZ_REAL(1.0), OZ_REAL(1.0));

  /* Two cylinders of height 1, rotate by 45 degrees around the z-axis,
     should overlap at x = 1 / sin(45).
  */

  oz_real_t d = OZ_REAL(1.0) / sinf(M_PI / OZ_REAL(4.0));

  oz_transform_t ta = {
    {OZ_REAL(0.0), OZ_REAL(0.0), OZ_REAL(0.0)},
    {OZ_REAL(0.0), OZ_REAL(0.0), OZ_SIN(M_PI / OZ_REAL(8.0)), OZ_COS(M_PI / OZ_REAL(8.0))},
    OZ_REAL(1.0)
  };

  oz_transform_t tb = {
    {d - OZ_REAL(0.01), OZ_REAL(0.0), OZ_REAL(0.0)},
    {OZ_REAL(0.0), OZ_REAL(0.0), OZ_SIN(M_PI / OZ_REAL(8.0)), OZ_COS(M_PI / OZ_REAL(8.0))},
    OZ_REAL(1.0)
  };

  oz_transform_t tc = {
    {d + OZ_REAL(0.01), OZ_REAL(0.0), OZ_REAL(0.0)},
    {OZ_REAL(0.0), OZ_REAL(0.0), OZ_SIN(M_PI / OZ_REAL(8.0)), OZ_COS(M_PI / OZ_REAL(8.0))},
    OZ_REAL(1.0)
  };

  int overlap_a = ozGJKOverlap(&ta, &cyl, &tb, &cyl);
  fprintf(stdout, "A & B overlap? %d\n", overlap_a);

  int overlap_b = ozGJKOverlap(&ta, &cyl, &tc, &cyl);
  fprintf(stdout, "A & B overlap? %d\n", overlap_b);

  oz_real_t dist[3];
  ozGJKDist(dist, &ta, &cyl, &tc, &cyl);
  oz_real_t len = OZ_SQRT(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
  fprintf(stdout, "Distance between A & B is %.3f\n", len);

  return 0;
}
