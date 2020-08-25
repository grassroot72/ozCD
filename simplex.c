/*
 * Copyright (c) 2017 Nikos Tasios
 * Copyright (C) 2019 Edward LEI <edward_lei72@hotmail.com>
 *
 * The code is licensed under the MIT license
 */

#include "type.h"
#include "vec3.h"
#include "apex_memmove.h"
#include "simplex.h"


/* Lookup table for single enabled bit position:
   Calculates the log2 of a 4-bit integer.
*/
static const oz_uchar_t s_pos[] = {0, 0, 1, 0, 2, 0, 0, 0, 3};
/*____________________________________^__^_____^___________^*/

/* Lookup table which tells us at which positions in the simplex array
   to get our points a, b, c, d. ex. if our bits are 0111 -> 7 -> {0, 1, 2}
*/
static const oz_uchar_t p_pos[16][3] =
{
  {0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
  {2, 0, 0}, {0, 2, 0}, {1, 2, 0}, {0, 1, 2},
  {3, 0, 0}, {0, 3, 0}, {1, 3, 0}, {0, 1, 3},
  {2, 3, 0}, {0, 2, 3}, {1, 2, 3}, {0, 0, 0}
};


#if defined(BARY_ERICSON)
static oz_real_t
_ozTriangleArea2D(const oz_real_t x1,
                  const oz_real_t y1,
                  const oz_real_t x2,
                  const oz_real_t y2,
                  const oz_real_t x3,
                  const oz_real_t y3)
{
  return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

/* Algorithm for calculating barycentric coordinates from
  'Real-time collision detection' by Christer Ericson.
*/
static void
_ozBarycentricCoordinates(oz_real_t *R,
                          const oz_real_t *P,
                          const oz_real_t *A,
                          const oz_real_t *B,
                          const oz_real_t *C)
{
  c u, v, w;

  oz_real_t m[3];
  {
    oz_real_t AB[3], AC[3];
    ozVec3Sub(AB, B, A);
    ozVec3Sub(AC, C, A);
    ozVec3Cross(m, AB, AC);
  }

  oz_real_t nu, nv, ood;
  oz_real_t x = OZ_FABS(m[0]), y = OZ_FABS(m[1]), z = OZ_FABS(m[2]);

  if (x >= y && x >= z) {
    nu = _ozTriangleArea2D(P[1], P[2], B[1], B[2], C[1], C[2]);
    nv = _ozTriangleArea2D(P[1], P[2], C[1], C[2], A[1], A[2]);
    ood = OZ_REAL(1.0) / m[0];
  }
  else if (y >= x && y >= z) {
    nu = _ozTriangleArea2D(P[0], P[2], B[0], B[2], C[0], C[2]);
    nv = _ozTriangleArea2D(P[0], P[2], C[0], C[2], A[0], A[2]);
    ood = OZ_REAL(1.0) / -m[1];
  }
  else {
    nu = _ozTriangleArea2D(P[0], P[1], B[0], B[1], C[0], C[1]);
    nv = _ozTriangleArea2D(P[0], P[1], C[0], C[1], A[0], A[1]);
    ood = 1.0 / m[2];
  }

  R[0] = nu * ood;
  R[1] = nv * ood;
  R[2] = OZ_REAL(1.0) - R[0] - R[1];
}

#elif defined(BARY_CRAMER)
static void
_ozBarycentricCoordinates(oz_real_t *R,
                          const oz_real_t *P,
                          const oz_real_t *A,
                          const oz_real_t *B,
                          const oz_real_t *C)
{
  oz_real_t v0[3], v1[3], v2[3];
  ozVec3Sub(v0, B, A);
  ozVec3Sub(v1, C, A);
  ozVec3Sub(v2, P, A);

  oz_real_t d00 = ozVec3Dot(v0, v0);
  oz_real_t d01 = ozVec3Dot(v0, v1);
  oz_real_t d02 = ozVec3Dot(v0, v2);
  oz_real_t d11 = ozVec3Dot(v1, v1);
  oz_real_t d12 = ozVec3Dot(v1, v2);
  oz_real_t denom = d00 * d11 - d01 * d01;

  R[0] = (d11 * d02 - d01 * d12) / denom;
  R[1] = (d00 * d12 - d01 * d02) / denom;
  R[2] = OZ_REAL(1.0) - R[0] - R[1];
}

#elif defined(BARY_GEPP)
static void
_ozBarycentricCoordinates(oz_real_t *R,
                          const oz_real_t *P,
                          const oz_real_t *A,
                          const oz_real_t *B,
                          const oz_real_t *C)
{
  oz_real_t v0[3], v1[3], v2[3];
  ozVec3Sub(v0, B, A);
  ozVec3Sub(v1, C, A);
  ozVec3Sub(v2, P, A);

  oz_real_t d00 = ozVec3Dot(v0, v0);
  oz_real_t d01 = ozVec3Dot(v0, v1);
  oz_real_t d02 = ozVec3Dot(v0, v2);
  oz_real_t d11 = ozVec3Dot(v1, v1);
  oz_real_t d12 = ozVec3Dot(v1, v2);

  R[0] = (d00 * d12 - d01 * d02) / (d00 * d11 - d01 * d01);
  R[1] = (d00 >= d01)? (d02 - d01 * R[0]) / d00: (d12 - d11 * R[0]) / d01;
  R[2] = OZ_REAL(1.0) - R[0] - R[1];
}
#endif


/*----- Simplex functions ---------------------------------------------------*/

void
ozSimplexAddPoint(oz_simplex_t *simplex,
                  const oz_real_t *point)
{
  oz_uchar_t b = ~simplex->_bits;  /* Flip bits */
  b &= -b; /* Last set (available) bit */

  oz_uchar_t pos = s_pos[b];  /* Get the bit position from the lookup table */
  simplex->_last_sb = pos;
  simplex->_bits |= b;        /* Insert the new bit */
  ++simplex->_size;
  apex_memcpy(simplex->_p + 3 * pos, point, 3 * sizeof(*point));
  oz_real_t l2 = ozVec3Len2(point);
  if (l2 > simplex->_max_vert2) {
    simplex->_max_vert2 = l2;
  }
}

void
ozSimplexAddPointWithInfo(oz_simplex_t *simplex,
                          const oz_real_t *point,
                          const oz_real_t *pa,
                          const oz_real_t *pb)
{
  ozSimplexAddPoint(simplex, point);
  apex_memcpy(simplex->_a + 3 * simplex->_last_sb, pa, 3 * sizeof(*pa));
  apex_memcpy(simplex->_b + 3 * simplex->_last_sb, pb, 3 * sizeof(*pb));
}

int
ozSimplexContains(const oz_simplex_t *simplex,
                  const oz_real_t *point)
{
  oz_uchar_t bits = simplex->_bits;

  int i;
  for (i = 0; i < 4; ++i, bits >>= 1) {
    if ((bits & 1) && ozVec3Eq(simplex->_p + 3 * i, point)) {
      return 1;
    }
  }

  return 0;
}

void
ozSimplexTranslate(oz_simplex_t *simplex,
                   const oz_real_t *dr)
{
  /* for(int k = 0; k < 4; ++k) p_[k] += dr; */
  simplex->_max_vert2 = OZ_REAL(0.0);
  oz_uchar_t bits = simplex->_bits;

  int i;
  for (i = 0; i < 4; ++i, bits >>= 1) {
    if (bits & 1) {
      ozVec3Add(simplex->_p + 3 * i, dr, simplex->_p + 3 * i);

      if(ozVec3Len2(simplex->_p + 3 * i) > simplex->_max_vert2) {
        simplex->_max_vert2 = ozVec3Len2(simplex->_p + 3 * i);
      }
    }
  }
}

void
ozSimplexComputeClosestPoints(oz_simplex_t *simplex,
                              oz_real_t *pa,
                              oz_real_t *pb,
                              const oz_real_t *P)
{
  switch (simplex->_size) {
    /* IMPORTANT: We are having accuracy problems with this projection. */
  case 3: {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];
    const oz_real_t *aA = simplex->_a + 3 * simplex->_last_sb;
    const oz_real_t *aB = simplex->_a + 3 * pos[0];
    const oz_real_t *aC = simplex->_a + 3 * pos[1];
    const oz_real_t *bA = simplex->_b + 3 * simplex->_last_sb;
    const oz_real_t *bB = simplex->_b + 3 * pos[0];
    const oz_real_t *bC = simplex->_b + 3 * pos[1];

    const oz_real_t *A = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *B = simplex->_p + 3 * pos[0];
    const oz_real_t *C = simplex->_p + 3 * pos[1];

    oz_real_t bary[3];
    _ozBarycentricCoordinates(bary, P, A, B, C);

    int i;
    for (i = 0; i < 3; ++i) {
      pa[i] = aA[i] * bary[0] + aB[i] * bary[1] + aC[i] * bary[2];
      pb[i] = bA[i] * bary[0] + bB[i] * bary[1] + bC[i] * bary[2];
    }

    break;
  }

  case 2: {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];
    const oz_real_t *aA = simplex->_a + 3 * simplex->_last_sb;
    const oz_real_t *aB = simplex->_a + 3 * pos[0];
    const oz_real_t *bA = simplex->_b + 3 * simplex->_last_sb;
    const oz_real_t *bB = simplex->_b + 3 * pos[0];
    double u, v;
    {
      const oz_real_t* A = simplex->_p + 3 * simplex->_last_sb;
      const oz_real_t* B = simplex->_p + 3 * pos[0];

      {
        oz_real_t AB[3], AP[3];
        ozVec3Sub(AB, B, A);
        ozVec3Sub(AP, P, A);
        v = ozVec3Dot(AB, AP) / ozVec3Len2(AB);
        u = OZ_REAL(1.0) - v;
      }
    }

    int i;
    for (i = 0; i < 3; ++i) {
      pa[i] = aA[i] * u + aB[i] * v;
      pb[i] = bA[i] * u + bB[i] * v;
    }

    break;
  }

  case 1: {
    apex_memcpy(pa, simplex->_a + 3 * simplex->_last_sb, 3 * sizeof(*pa));
    apex_memcpy(pb, simplex->_b + 3 * simplex->_last_sb, 3 * sizeof(*pb));
    break;
  }

  default:
    break;
  }
}

/* Reduce the simplex to its feature closest to the origin and
   update the search direction, dir.
*/
void
ozSimplexClosest(oz_simplex_t *simplex,
                 oz_real_t *dir)
{
  switch (simplex->_size) {
  case 4:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];
    const oz_real_t *c = simplex->_p + 3 * pos[1];
    const oz_real_t *d = simplex->_p + 3 * pos[2];

    oz_real_t ab[3], ac[3], ad[3];
    ozVec3Sub(ab, b, a);
    ozVec3Sub(ac, c, a);
    ozVec3Sub(ad, d, a);

    /*----- Vertex Case -----------------------------*/

    oz_real_t dot_aba = ozVec3Dot(ab, a);
    oz_real_t dot_aca = ozVec3Dot(ac, a);
    oz_real_t dot_ada = ozVec3Dot(ad, a);

    if (dot_aba >= OZ_REAL(0.0) &&
        dot_aca >= OZ_REAL(0.0) &&
        dot_ada >= OZ_REAL(0.0)) {
      /* Take direction passing through origin */
      apex_memcpy(dir, a, 3 * sizeof(*dir));
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[1]);
      ozSimplexRemovePoint(simplex, pos[2]);
      break;
    }

    /*----- Edge Cases ------------------------------*/

    /* ab Edge case */
    oz_real_t dot_abb = ozVec3Dot(ab, b);
    oz_real_t dot_abPerp1 = dot_aba * ozVec3Dot(ac, b) - dot_abb * dot_aca;
    oz_real_t dot_abPerp2 = dot_aba * ozVec3Dot(ad, b) - dot_abb * dot_ada;
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces abc, abd
    */
    if (dot_abPerp1 <= OZ_REAL(0.0) &&
        dot_abPerp2 <= OZ_REAL(0.0) &&
        dot_aba <= OZ_REAL(0.0)) {
      oz_real_t f = dot_aba / (dot_aba - dot_abb);
      ozVec3FMAdd(dir, f, ab, a);
      ozSimplexRemovePoint(simplex, pos[1]);
      ozSimplexRemovePoint(simplex, pos[2]);
      break;
    }

    /* ac Edge case */
    oz_real_t dot_acc = ozVec3Dot(ac, c);
    oz_real_t dot_acPerp1 = dot_aca * ozVec3Dot(ad, c) - dot_acc * dot_ada;
    oz_real_t dot_acPerp2 = dot_aca * ozVec3Dot(ab, c) - dot_acc * dot_aba;
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces abc, acd
    */
    if (dot_acPerp1 <= OZ_REAL(0.0) &&
        dot_acPerp2 <= OZ_REAL(0.0) &&
        dot_aca <= OZ_REAL(0.0)) {
      oz_real_t f = dot_aca / (dot_aca - dot_acc);
      ozVec3FMAdd(dir, f, ac, a);
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[2]);
      break;
    }

    /* ad Edge case */
    oz_real_t dot_add = ozVec3Dot(ad, d);
    oz_real_t dot_adPerp1 = dot_ada * ozVec3Dot(ab, d) - dot_add * dot_aba;
    oz_real_t dot_adPerp2 = dot_ada * ozVec3Dot(ac, d) - dot_add * dot_aca;
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces acd, abd
    */
    if (dot_adPerp1 <= OZ_REAL(0.0) &&
        dot_adPerp2 <= OZ_REAL(0.0) &&
        dot_ada <= OZ_REAL(0.0)) {
      oz_real_t f = dot_ada / (dot_ada - dot_add);
      ozVec3FMAdd(dir, f, ad, a);
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[1]);
      break;
    }

    /*----- Face Cases ------------------------------*/

    /* On abc side */
    /* The origin should be on abc's side and
       between the half-spaces defined by ac and ab (normal to abc)
    */
    {
      oz_real_t abxac[3];
      ozVec3Cross(abxac, ab, ac);

      if (ozVec3Dot(ad, abxac) * ozVec3Dot(a, abxac) > OZ_REAL(0.0) &&
          dot_abPerp1 >= OZ_REAL(0.0) && dot_acPerp2 >= OZ_REAL(0.0)) {
        /* Remove point d */
        ozSimplexRemovePoint(simplex, pos[2]);
        ozVec3SMul(dir, ozVec3Dot(abxac, a) / ozVec3Len2(abxac), abxac);
        break;
      }
    }

    /* On abd side */
    /* The origin should be on abd's side and
       between the half-spaces defined by ab and ad (normal to abd)
    */
    {
      oz_real_t abxad[3];
      ozVec3Cross(abxad, ab, ad);

      if (ozVec3Dot(ac, abxad) * ozVec3Dot(a, abxad) > OZ_REAL(0.0) &&
          dot_abPerp2 >= OZ_REAL(0.0) && dot_adPerp1 >= OZ_REAL(0.0)) {
        /* Remove point c */
        ozSimplexRemovePoint(simplex, pos[1]);
        ozVec3SMul(dir, ozVec3Dot(abxad, a) / ozVec3Len2(abxad), abxad);
        break;
      }
    }

    /* On acd side */
    /* The origin should be on acd's side and
       between the half-spaces defined by ac and ad (normal to acd)
    */
    {
      oz_real_t acxad[3];
      ozVec3Cross(acxad, ac, ad);

      if (ozVec3Dot(ab, acxad) * ozVec3Dot(a, acxad) > OZ_REAL(0.0) &&
          dot_acPerp1 >= OZ_REAL(0.0) && dot_adPerp2 >= OZ_REAL(0.0)) {
        /* Remove point b */
        ozSimplexRemovePoint(simplex, pos[0]);
        ozVec3SMul(dir, ozVec3Dot(acxad, a) / ozVec3Len2(acxad), acxad);
        break;
      }
    }

    /* On bcd side */
    /* The origin should be on bcd's side */
    {
      oz_real_t bc[3], bd[3], bcxbd[3];
      ozVec3Sub(bc, c, b);
      ozVec3Sub(bd, d, b);
      ozVec3Cross(bcxbd, bc, bd);

      if (ozVec3Dot(bcxbd, ab) * ozVec3Dot(bcxbd, b) < OZ_REAL(0.0)) {
        /* Remove point a */
        ozSimplexRemovePoint(simplex, simplex->_last_sb);
        simplex->_last_sb = pos[0];
        ozVec3SMul(dir, ozVec3Dot(bcxbd, b) / ozVec3Len2(bcxbd), bcxbd);
        break;
      }
    }

    /* 'else' should only be when the origin is inside the tetrahedron */
    break;
  }

  case 3:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];
    const oz_real_t *c = simplex->_p + 3 * pos[1];

    oz_real_t ab[3], ac[3];
    ozVec3Sub(ab, b, a);
    ozVec3Sub(ac, c, a);

    /* Check if O in vertex region A */
    oz_real_t dot_aba = -ozVec3Dot(ab, a);
    oz_real_t dot_aca = -ozVec3Dot(ac, a);

    if (dot_aba <= OZ_REAL(0.0) && dot_aca <= OZ_REAL(0.0)) {
      /* Take direction passing through origin */
      apex_memcpy(dir, a, 3 * sizeof(*dir));
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[1]);
      break;
    }

    /* Check if O in edge region AB */
    oz_real_t dot_abb = -ozVec3Dot(ab, b);
    oz_real_t dot_acb = -ozVec3Dot(ac, b);
    oz_real_t vc = dot_aba * dot_acb - dot_abb * dot_aca;

    if (vc <= OZ_REAL(0.0) && dot_aba >= OZ_REAL(0.0) && dot_abb <= OZ_REAL(0.0)) {
      oz_real_t f = dot_aba / (dot_aba - dot_abb);
      ozVec3FMAdd(dir, f, ab, a);
      /* Remove Point c */
      ozSimplexRemovePoint(simplex, pos[1]);
      break;
    }

    /* Check if O in edge region AC */
    oz_real_t dot_abc = -ozVec3Dot(ab, c);
    oz_real_t dot_acc = -ozVec3Dot(ac, c);
    oz_real_t vb = dot_abc * dot_aca - dot_aba * dot_acc;

    if (vb <= OZ_REAL(0.0) && dot_aca >= OZ_REAL(0.0) && dot_acc <= OZ_REAL(0.0)) {
      oz_real_t f = dot_aca / (dot_aca - dot_acc);
      ozVec3FMAdd(dir, f, ac, a);
      /* Remove Point b */
      ozSimplexRemovePoint(simplex, pos[0]);
      break;
    }

    oz_real_t va = dot_abb * dot_acc - dot_abc * dot_acb;
    oz_real_t w = OZ_REAL(1.0) / (va + vb + vc);

    int i;
    for (i = 0; i < 3; ++i) {
      dir[i] = a[i] + (ab[i] * vb + ac[i] * vc) * w;
    }
    break;
  }

  case 2:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];

    oz_real_t  ab[3];
    ozVec3Sub(ab, b, a);

    oz_real_t t = -ozVec3Dot(ab, a);

    if (t <= OZ_REAL(0.0)) {
      /* Take direction passing through origin */
      apex_memcpy(dir, a, 3 * sizeof(*dir));
      ozSimplexRemovePoint(simplex, pos[0]);
      break;
    }

    double denom = ozVec3Len2(ab);
    if (t >= denom) {
      ozSimplexRemovePoint(simplex, simplex->_last_sb);
      simplex->_last_sb = pos[0];
      apex_memcpy(dir, b, 3 * sizeof(*dir));
      break;
    }

    ozVec3FMAdd(dir, t / denom, ab, a);
    break;
  }

  case 1:
  {
    apex_memcpy(dir, simplex->_p + 3 * simplex->_last_sb, 3 * sizeof(*dir));
    break;
  }

  default: break;
  }
}

/* Check if the origin is contained in the simplex and
   update the search direction, dir.
*/
int
ozSimplexContainsOrigin(oz_simplex_t *simplex,
                        oz_real_t *dir)
{
  switch (simplex->_size) {
  case 4:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];
    const oz_real_t *c = simplex->_p + 3 * pos[1];
    const oz_real_t *d = simplex->_p + 3 * pos[2];

    oz_real_t ab[3], ac[3], ad[3];
    ozVec3Sub(ab, b, a);
    ozVec3Sub(ac, c, a);
    ozVec3Sub(ad, d, a);

    /*----- Face Cases --------------------------------*/

    /* On abc side */
    oz_real_t abxac[3];
    ozVec3Cross(abxac, ab, ac);

    int abPerp1Pos, acPerp2Pos;
    {
      oz_real_t cross_temp_a[3], cross_temp_b[3];
      ozVec3Cross(cross_temp_a, abxac, ab);
      ozVec3Cross(cross_temp_b, ac, abxac);
      abPerp1Pos = (ozVec3Dot(cross_temp_a, a) > OZ_REAL(0.0));
      acPerp2Pos = (ozVec3Dot(cross_temp_b, a) > OZ_REAL(0.0));
    }

    /* The origin should be on abc's side and
       between the half-spaces defined by ac and ab (normal to abc)
    */
    {
      oz_real_t abcPerp[3];

      if (ozVec3Dot(abxac, ad) > OZ_REAL(0.0)) {
        ozVec3SMul(abcPerp, -OZ_REAL(1.0), abxac);
      }
      else {
        apex_memcpy(abcPerp, abxac, 3 * sizeof(*abcPerp));
      }

      if ((ozVec3Dot(abcPerp, a) < OZ_REAL(0.0)) && !abPerp1Pos && !acPerp2Pos) {
        /* Remove point d */
        ozSimplexRemovePoint(simplex, pos[2]);
        apex_memcpy(dir, abcPerp, 3 * sizeof(*dir));
        break;
      }
    }

    /* On abd side */
    oz_real_t abxad[3];
    ozVec3Cross(abxad, ab, ad);

    int abPerp2Pos, adPerp1Pos;
    {
      oz_real_t cross_temp_a[3], cross_temp_b[3];
      ozVec3Cross(cross_temp_a, abxad, ab);
      ozVec3Cross(cross_temp_b, ad, abxad);
      abPerp2Pos = (ozVec3Dot(cross_temp_a, a) > OZ_REAL(0.0));
      adPerp1Pos = (ozVec3Dot(cross_temp_b, a) > OZ_REAL(0.0));
    }

    /* The origin should be on abd's side and
       between the half-spaces defined by ab and ad (normal to abd)
    */
    {
      oz_real_t abdPerp[3];

      if (ozVec3Dot(abxad, ac) > OZ_REAL(0.0)) {
        ozVec3SMul(abdPerp, -OZ_REAL(1.0), abxad);
      }
      else {
        apex_memcpy(abdPerp, abxad, 3 * sizeof(*abdPerp));
      }

      if ((ozVec3Dot(abdPerp, a) < OZ_REAL(0.0)) && !abPerp2Pos && !adPerp1Pos) {
        /* Remove point c */
        ozSimplexRemovePoint(simplex, pos[1]);
        apex_memcpy(dir, abdPerp, 3 * sizeof(*dir));
        break;
      }
    }

    /* On acd side */
    oz_real_t acxad[3];
    ozVec3Cross(acxad, ac, ad);

    int acPerp1Pos, adPerp2Pos;
    {
      oz_real_t cross_temp_a[3], cross_temp_b[3];
      ozVec3Cross(cross_temp_a, acxad, ac);
      ozVec3Cross(cross_temp_b, ad, acxad);
      acPerp1Pos = (ozVec3Dot(cross_temp_a, a) > OZ_REAL(0.0));
      adPerp2Pos = (ozVec3Dot(cross_temp_b, a) > OZ_REAL(0.0));
    }

    /* The origin should be on acd's side and
       between the half-spaces defined by ac and ad (normal to acd)
    */
    {
      oz_real_t acdPerp[3];

      if (ozVec3Dot(acxad, ab) > OZ_REAL(0.0)) {
        ozVec3SMul(acdPerp, -OZ_REAL(1.0), acxad);
      }
      else {
        apex_memcpy(acdPerp, acxad, 3 * sizeof(*acdPerp));
      }

      if ((ozVec3Dot(acdPerp, a) < OZ_REAL(0.0)) && !acPerp1Pos && !adPerp2Pos) {
        /* Remove point b */
        ozSimplexRemovePoint(simplex, pos[0]);
        apex_memcpy(dir, acdPerp, 3 * sizeof(*dir));
        break;
      }
    }

    /*----- Edge Cases --------------------------------*/

    /* ab Edge case */
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces abc, abd
    */
    if (abPerp1Pos && abPerp2Pos) {
      ozVec3TripleProduct(dir, a, ab, ab);
      ozSimplexRemovePoint(simplex, pos[1]);
      ozSimplexRemovePoint(simplex, pos[2]);
      break;
    }

    /* ac Edge case */
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces abc, acd
    */
    if (acPerp1Pos && acPerp2Pos) {
      ozVec3TripleProduct(dir, a, ac, ac);
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[2]);
      break;
    }

    /* ad Edge case */
    /* The origin must be inside the space defined by the intersection
       of two half-space normal to the adjacent faces acd, abd
    */
    if (adPerp1Pos && adPerp2Pos) {
      ozVec3TripleProduct(dir, a, ad, ad);
      ozSimplexRemovePoint(simplex, pos[0]);
      ozSimplexRemovePoint(simplex, pos[1]);
      break;
    }

    /* 'else' should only be when the origin is inside the tetrahedron */
    return 1;
  }

  case 3:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];
    const oz_real_t *c = simplex->_p + 3 * pos[1];

    oz_real_t ab[3], ac[3], abxac[3];
    ozVec3Sub(ab, b, a);
    ozVec3Sub(ac, c, a);
    ozVec3Cross(abxac, ab, ac);

    /*----- Edge Cases --------------------------------*/

    /* Origin on the outside of triangle and close to ab */
    oz_real_t cross_temp[3];
    ozVec3Cross(cross_temp, ab, abxac);

    if (ozVec3Dot(cross_temp, a) < OZ_REAL(0.0)) {
      ozVec3TripleProduct(dir, a, ab, ab);
      /* Remove Point c */
      ozSimplexRemovePoint(simplex, pos[1]);
      break;
    }

    /* Origin on the outside of triangle and close to ac */
    ozVec3Cross(cross_temp, abxac, ac);

    if (ozVec3Dot(cross_temp, a) < OZ_REAL(0.0)) {
      ozVec3TripleProduct(dir, a, ac, ac);
      /* Remove Point b */
      ozSimplexRemovePoint(simplex, pos[0]);
      break;
    }

    /*----- Face Case --------------------------------*/
    if (ozVec3Dot(abxac, a) > OZ_REAL(0.0)) {
      ozVec3SMul(dir, -OZ_REAL(1.0), abxac);
    }
    else {
      apex_memcpy(dir, abxac, 3 * sizeof(*dir));
    }

    break;
  }

  case 2:
  {
    const oz_uchar_t *pos = p_pos[(simplex->_bits ^ (1 << simplex->_last_sb))];

    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    const oz_real_t *b = simplex->_p + 3 * pos[0];

    oz_real_t ab[3];
    ozVec3Sub(ab, b, a);

    ozVec3TripleProduct(dir, a, ab, ab);
    break;
  }

  case 1:
  {
    const oz_real_t *a = simplex->_p + 3 * simplex->_last_sb;
    ozVec3SMul(dir, -OZ_REAL(1.0), a);
    break;
  }

  default: break;
  }

  return 0;
}
