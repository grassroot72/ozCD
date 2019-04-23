/*
  shape.c

  Copyright (c) 2017 Nikos Tasios
  Copyright (C) 2019 by Edward LEI

  This code is licensed under the MIT license
*/

#include <stdlib.h>
#include <string.h>

#include "type.h"
#include "vec3.h"
#include "quat.h"
#include "apex_memmove.h"
#include "shape.h"


/*----- static functions ----------------------------------------------------*/

static void
_ozSphereSupport(oz_real_t *support_point,
                 const void *shape,
                 const oz_real_t *dir)
{
  oz_real_t norm = ozVec3Len(dir);

  if (norm > OZ_REAL(0.0)) {
    support_point[0] = dir[0] / norm;
    support_point[1] = dir[1] / norm;
    support_point[2] = dir[2] / norm;
  }
  else {
    support_point[0] = OZ_REAL(1.0);
    support_point[1] = OZ_REAL(0.0);
    support_point[2] = OZ_REAL(0.0);
  }
}

static void
_ozPointSupport(oz_real_t *support_point,
                const void *shape,
                const oz_real_t *dir)
{
  memset(support_point, 0, 3 * sizeof(*dir));
}

static void
_ozLineSupport(oz_real_t *support_point,
               const void *shape,
               const oz_real_t *dir)
{
  support_point[0] = OZ_REAL(0.0);
  support_point[1] = OZ_COPYSIGN(OZ_REAL(1.0), dir[1]);
  support_point[2] = OZ_REAL(0.0);
}

static void
_ozDiskSupport(oz_real_t *support_point,
               const void *shape,
               const oz_real_t *dir)
{
  oz_real_t length2 = dir[0] * dir[0] + dir[2] * dir[2];

  if (ozIsZero(length2)) {
    support_point[0] = OZ_REAL(1.0);
    support_point[1] = OZ_REAL(0.0);
    support_point[2] = OZ_REAL(0.0);
  }
  else {
    oz_real_t length = OZ_REAL(1.0) / OZ_SQRT(length2);
    support_point[0] = dir[0] * length;
    support_point[1] = OZ_REAL(0.0);
    support_point[2] = dir[2] * length;
  }
}

static void
_ozMeshSupport(oz_real_t *support_point,
               const void *shape,
               const oz_real_t *dir)
{
  oz_mesh_t mesh = *(const oz_mesh_t *)shape;

  oz_uint_t next = 0, last = 0, curr = 0;
  oz_real_t p = OZ_REAL(0.0);
  oz_real_t max = ozVec3Dot(mesh._vertices, dir);

  for(;;) {
    for (oz_uint_t vid = 0; vid < mesh._n_vert_neighbours[curr]; ++vid) {
      next = mesh._vert_neighbours[curr][vid];

      if (next != last) {
        p = ozVec3Dot(mesh._vertices + 3 * next, dir);

        if (p > max) {
          max = p;
          last = curr;
          curr = next;
          break;
        }
      }

      if (vid == mesh._n_vert_neighbours[curr] - 1) {
        apex_memcpy(support_point, mesh._vertices + 3 * curr, 3 * sizeof(*support_point));
        return;
      }
    }
  }
}

static void
_ozCylinderSupport(oz_real_t *support_point,
                   const void *shape,
                   const oz_real_t *dir)
{
  oz_cylinder_t cyl = *(const oz_cylinder_t *)shape;

  oz_real_t length = OZ_SQRT(dir[0] * dir[0] + dir[2] * dir[2]);

  if (ozIsZero(length)) {
    support_point[0] = cyl._base_radius;
    support_point[1] = OZ_COPYSIGN(cyl._half_height, dir[1]);
    support_point[2] = OZ_REAL(0.0);
  }
  else {
    oz_real_t d = cyl._base_radius / length;
    support_point[0] = d * dir[0];
    support_point[1] = OZ_COPYSIGN(cyl._half_height, dir[1]);
    support_point[2] = d * dir[2];
  }
}

static void
_ozBoxSupport(oz_real_t *support_point,
              const void *shape,
              const oz_real_t *dir)
{
  const oz_real_t *extent = ((const oz_box_t *)shape)->_extent;
  support_point[0] = OZ_COPYSIGN(extent[0], dir[0]);
  support_point[1] = OZ_COPYSIGN(extent[1], dir[1]);
  support_point[2] = OZ_COPYSIGN(extent[2], dir[2]);
}

static void
_ozConeSupport(oz_real_t *support_point,
               const void *shape,
               const oz_real_t *dir)
{
  oz_cone_t cone = *(const oz_cone_t *)shape;

  oz_real_t test = dir[1] / ozVec3Len(dir);
  if (test >= cone._sintheta) {
    support_point[0] = OZ_REAL(0.0);
    support_point[1] = cone._half_height;
    support_point[2] = OZ_REAL(0.0);
  }
  else if (test < cone._sintheta && (dir[0] != OZ_REAL(0.0) || dir[2] != OZ_REAL(0.0))) {
    oz_real_t length = OZ_SQRT(dir[0] * dir[0] + dir[2] * dir[2]);
    support_point[0] = cone._base_radius * dir[0] / length;
    support_point[1] = -cone._half_height;
    support_point[2] = cone._base_radius * dir[2] / length;
  }
  else {
    support_point[0] = OZ_REAL(0.0);
    support_point[1] = -cone._half_height;
    support_point[2] = OZ_REAL(0.0);
  }
}

static void
_ozBiconeSupport(oz_real_t *support_point,
                 const void *shape,
                 const oz_real_t *dir)
{
  oz_bicone_t bicone = *(const oz_bicone_t *)shape;

  oz_real_t test = dir[1] / ozVec3Len(dir);
  if (test >= bicone._sintheta) {
    support_point[0] = OZ_REAL(0.0);
    support_point[1] = bicone._half_height;
    support_point[2] = OZ_REAL(0.0);
  }
  else if (test <= -bicone._sintheta) {
    support_point[0] = OZ_REAL(0.0);
    support_point[1] = -bicone._half_height;
    support_point[2] = OZ_REAL(0.0);
  }
  else {
    oz_real_t factor = bicone._base_radius / OZ_SQRT(dir[0] * dir[0] + dir[2] * dir[2]);
    support_point[0] = factor * dir[0];
    support_point[1] = OZ_REAL(0.0);
    support_point[2] = factor * dir[2];
  }
}

static void
_ozLeafCylinderSupport(oz_real_t *support_point,
                       const void *shape,
                       const oz_real_t *dir)
{
  oz_leaf_cylinder_t leaf = *(const oz_leaf_cylinder_t *)shape;

  oz_real_t x = OZ_REAL(0.0), z = OZ_REAL(0.0);
  if (dir[0] != OZ_REAL(0.0) || dir[2] != OZ_REAL(0.0)) {
    oz_real_t l = OZ_SQRT(dir[0] * dir[0] + dir[2] * dir[2]);
    oz_real_t test = dir[2] / l;

    if (test >= leaf._half_length / leaf._circle_radius) z = leaf._half_length;
    else if (test <= -leaf._half_length / leaf._circle_radius) z = -leaf._half_length;
    else {
      x = leaf._circle_radius * dir[0] / l - OZ_COPYSIGN(leaf._circle_distance, dir[0]);
      z = leaf._circle_radius * dir[2] / l;
    }
  }

  support_point[0] = x;
  support_point[1] = OZ_COPYSIGN(leaf._half_height, dir[1]);
  support_point[2] = z;
}

static void
_ozSphereSweptSupport(oz_real_t *support_point,
                      const void *shape,
                      const oz_real_t *dir)
{
  oz_sphere_swept_t swept = *(const oz_sphere_swept_t *)shape;
  ozpfnSupport support = *(const ozpfnSupport*)swept._shape;

  oz_real_t sphere_point[3];
  _ozSphereSupport(sphere_point, NULL, dir);
  support(support_point, swept._shape, dir);
  ozVec3FMAdd(support_point, swept._radius, sphere_point, support_point);
}

static void
_ozHullSupport(oz_real_t *support_point,
               const void *shape,
               const oz_real_t *dir)
{
  oz_hull_t hull = *(const oz_hull_t *)shape;
  ozpfnSupport sa = *(const ozpfnSupport *)hull._ca;
  ozpfnSupport sb = *(const ozpfnSupport *)hull._cb;

  oz_real_t inv_dir_a[3];
  ozQuatVec3Rotate(inv_dir_a, hull._inv_rot_a, dir);
  sa(support_point, hull._ca, inv_dir_a);
  ozQuatVec3Rotate(support_point, hull._ta.rot, support_point);
  ozVec3FMAdd(support_point, hull._ta.size, support_point, hull._ta.pos);

  oz_real_t inv_dir_b[3], temp_point[3];
  ozQuatVec3Rotate(inv_dir_b, hull._inv_rot_b, dir);
  sb(temp_point, hull._cb, inv_dir_b);
  ozQuatVec3Rotate(temp_point, hull._tb.rot, temp_point);
  ozVec3FMAdd(temp_point, hull._tb.size, temp_point, hull._tb.pos);

  if (ozVec3Dot(temp_point, dir) > ozVec3Dot(support_point, dir)) {
    apex_memcpy(support_point, temp_point, 3 * sizeof(*support_point));
  }
}

static void
_ozMinkowskiSumSupport(oz_real_t *support_point,
                       const void *shape,
                       const oz_real_t *dir)
{
  oz_minkowski_sum_t msum = *(const oz_minkowski_sum_t *)shape;
  ozpfnSupport sa = *(const ozpfnSupport*)msum._ca;
  ozpfnSupport sb = *(const ozpfnSupport*)msum._cb;

  sa(support_point, msum._ca, dir);

  oz_real_t inv_dir[3], support_point_b[3];
  ozQuatVec3Rotate(inv_dir, msum._inv_rot, dir);
  sb(support_point_b, msum._cb, inv_dir);
  ozQuatVec3Rotate(support_point_b, msum._t.rot, support_point_b);
  ozVec3FMAdd(support_point_b, msum._t.size, support_point_b, msum._t.pos);

  ozVec3Add(support_point, support_point, support_point_b);
}

/*----- end: static functions -----------------------------------------------*/


void
ozSphereInit(oz_sphere_t *sphere)
{
  sphere->support = _ozSphereSupport;
}

void
ozPointInit(oz_point_t *point)
{
  point->support = _ozPointSupport;
}

void
ozLineInit(oz_line_t *line)
{
  line->support = _ozLineSupport;
}

void
ozDiskInit(oz_disk_t *disk)
{
  disk->support = _ozDiskSupport;
}

void
ozMeshInit(oz_mesh_t *mesh,
           const oz_uint_t n_vertices,
           const oz_real_t *vertices,
           const oz_uint_t n_faces,
           const oz_uint_t *face_start,
           const oz_uint_t *faces)
{
  mesh->support = _ozMeshSupport;
  mesh->_n_vertices = n_vertices;
  mesh->_vertices = (oz_real_t *)malloc(3 * n_vertices * sizeof(*vertices));
  apex_memcpy(mesh->_vertices, vertices, 3 * n_vertices * sizeof(*vertices));

  /** Find Edges **/
  /* Euler's formula */
  oz_uint_t n_edges = n_vertices + n_faces - 2;
  /* Each edge is represented by its two vertex indices */
  oz_uint_t *edges = (oz_uint_t *)malloc(n_edges * 2 * sizeof(*edges));

  /* Iterate over each face and for each next face in the list, check if they
     share two vertices, this defines an edge
  */
  for (oz_uint_t fi = 0, eid = 0; fi < n_faces; ++fi) {
    const oz_uint_t *face = faces + face_start[fi];
    oz_uint_t face_size = ((fi < n_faces - 1)? face_start[fi + 1]: n_vertices) - face_start[fi];

    for (oz_uint_t fj = fi + 1; fj < n_faces; ++fj) {
      oz_uint_t fcount = 0;
      oz_uint_t edge[3];

      for (oz_uint_t i = 0; i < face_size; ++i) {
        oz_uint_t vid_fi = face[i];

        for (oz_uint_t j = 0; j < face_size; ++j) {
          if (vid_fi == faces[face_start[fj] + j]) {
            edge[fcount++] = vid_fi;
          }
        }

        if (fcount == 2) {
          apex_memcpy(edges + 2 * eid, edge, 2 * sizeof(*edge));
          ++eid;
          fcount = 0;
        }
      }
    }
  }

  /** Find Vertex Neighbours **/
  /* For all vertices, check if two edges share this vertex. If they do and it
     isn't vertex 0, append the other vertices of these edge to the neighbor list
  */
  mesh->_vert_neighbours = (oz_uint_t **)malloc(n_vertices * sizeof(*mesh->_vert_neighbours));
  mesh->_n_vert_neighbours = (oz_uint_t *)calloc(n_vertices, sizeof(*mesh->_n_vert_neighbours));

  for (oz_uint_t vid = 0; vid < n_vertices; ++vid) {
    oz_uint_t capacity = 5;
    oz_uint_t n_neighbours = 0;
    mesh->_vert_neighbours[vid] = (oz_uint_t *)malloc(capacity * sizeof(**mesh->_vert_neighbours));
    oz_uint_t *neighbours = mesh->_vert_neighbours[vid];

    for (oz_uint_t ei = 0; ei < n_edges; ++ei) {
      for (oz_uint_t i = 0; i < 2; ++i) {
        if (edges[2 * ei + i] == vid && edges[2 * ei + (i + 1) % 2] != 0) {
          neighbours[n_neighbours++] = edges[2 * ei + (i + 1) % 2];

          if (n_neighbours == capacity) {
            capacity *= 2;
            oz_uint_t *temp = (oz_uint_t *)realloc(mesh->_vert_neighbours[vid], capacity * sizeof(**mesh->_vert_neighbours));
            mesh->_vert_neighbours[vid] = temp;
            neighbours = temp;
          }
        }
      }
    }

    if (n_neighbours) mesh->_n_vert_neighbours[vid] = n_neighbours;
  }

  free(edges);
}

void
ozMeshTerminate(oz_mesh_t *mesh)
{
  free(mesh->_vertices);
  free(mesh->_n_vert_neighbours);

  for (oz_uint_t vid = 0; vid < mesh->_n_vertices; ++vid) {
    free(mesh->_vert_neighbours[vid]);
  }

  free(mesh->_vert_neighbours);
}

void
ozCylinderInit(oz_cylinder_t *cylinder,
               oz_real_t base_radius,
               oz_real_t height)
{
  cylinder->support = _ozCylinderSupport;
  cylinder->_base_radius = base_radius;
  cylinder->_half_height = OZ_REAL(0.5) * height;
}

void
ozBoxInit(oz_box_t *box,
          const oz_real_t *extent)
{
  box->support = _ozBoxSupport;
  apex_memcpy(box->_extent, extent, 3 * sizeof(*extent));
}

void
ozConeInit(oz_cone_t *cone,
           oz_real_t base_radius,
           oz_real_t height)
{
  cone->support = _ozConeSupport;
  cone->_base_radius = base_radius;
  cone->_half_height = OZ_REAL(0.5) * height;
  cone->_sintheta = base_radius / OZ_SQRT(base_radius * base_radius + height * height);
}

void
ozBiconeInit(oz_cone_t *bicone,
             oz_real_t base_radius,
             oz_real_t height)
{
  bicone->support = _ozBiconeSupport;
  bicone->_base_radius = base_radius;
  bicone->_half_height = OZ_REAL(0.5) * height;
  bicone->_sintheta = base_radius / OZ_SQRT(base_radius * base_radius + bicone->_half_height * bicone->_half_height);
}

void
ozLeafCylinderInit(oz_leaf_cylinder_t *leaf,
                   oz_real_t width,
                   oz_real_t length,
                   oz_real_t height)
{
  leaf->support = _ozLeafCylinderSupport;
  leaf->_half_width = OZ_REAL(0.5) * width;
  leaf->_half_length = OZ_REAL(0.5) * length;
  leaf->_half_height = OZ_REAL(0.5) * height;
  leaf->_circle_radius = OZ_REAL(0.25) * (length * length + width * width) / width;
  leaf->_circle_distance = OZ_REAL(0.25) * (length * length - width * width) / width;
}

void
ozSphereSweptInit(oz_sphere_swept_t *swept,
                  oz_real_t radius,
                  const void *shape)
{
  swept->support = _ozSphereSweptSupport;
  swept->_shape = shape;
  swept->_radius = radius;
}

void
ozHullInit(oz_hull_t *hull,
           const oz_transform_t *ta,
           const void *ca,
           const oz_transform_t *tb,
           const void *cb)
{
  hull->support = _ozHullSupport;
  hull->_ca = ca;
  hull->_cb = cb;
  hull->_ta = *ta;
  hull->_tb = *tb;

  ozQuatInverse(hull->_inv_rot_a, hull->_ta.rot);
  ozQuatInverse(hull->_inv_rot_b, hull->_tb.rot);
}

void
ozMinkowskiSumInit(oz_minkowski_sum_t *msum,
                   const oz_transform_t *t,
                   const void *ca,
                   const void *cb)
{
  msum->support = _ozMinkowskiSumSupport;
  msum->_ca = ca;
  msum->_cb = cb;
  msum->_t = *t;

  ozQuatInverse(msum->_inv_rot, msum->_t.rot);
}
