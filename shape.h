/*
  shape.h

  Copyright (c) 2017 Nikos Tasios
  Copyright (C) 2019 by Edward LEI

  This code is licensed under the MIT license
*/

#ifndef _OZ_SHAPE_H_
#define _OZ_SHAPE_H_

#include "mathfun.h"


/* The dynamic transform data of a 3d object:
   - pos[3]  <- position
   - rot[4]  <- rotation
   - scale   <- scale
*/
struct _oz_transform_s {
  oz_real_t pos[3];
  oz_real_t rot[4];    /* quaternion */
  oz_real_t size;      /* scale of the 3d object */
};
typedef struct _oz_transform_s oz_transform_t;


/* support function pointer used by GJK */
typedef void (*ozpfnSupport)(
  oz_real_t *support_point,
  const void *shape,
  const oz_real_t *dir
);


/* Sphere: radius = 1 */
struct _oz_sphere_s {
  ozpfnSupport support;
};
typedef struct _oz_sphere_s oz_sphere_t;

/* Point */
struct _oz_point_s {
  ozpfnSupport support;
};
typedef struct _oz_point_s oz_point_t;

/* Line: length = 1, parallel to the y-axis */
struct _oz_line_s {
  ozpfnSupport support;
};
typedef struct _oz_line_s oz_line_t;

/* Disk: radius = 1, its axis is parallel to the y-axis */
struct _oz_disk_s {
  ozpfnSupport support;
};
typedef struct _oz_disk_s oz_disk_t;

/* Mesh */
struct _oz_mesh_s {
  ozpfnSupport support;
  oz_uint_t _n_vertices;
  oz_real_t *_vertices;
  oz_uint_t *_n_vert_neighbours;
  oz_uint_t **_vert_neighbours;
};
typedef struct _oz_mesh_s oz_mesh_t;

/* Cylinder: its axis is parallel to the y-axis */
struct _oz_cylinder_s {
  ozpfnSupport support;
  oz_real_t _base_radius, _half_height;
};
typedef struct _oz_cylinder_s oz_cylinder_t;

/* Box: length is given by the extend array */
struct _oz_box_s {
  ozpfnSupport support;
  oz_real_t _extent[3];
};
typedef struct _oz_box_s oz_box_t;

/* Cone: its axis is parallel to the y-axis */
struct _oz_cone_s {
  ozpfnSupport support;
  oz_real_t _base_radius, _half_height;
  oz_real_t _sintheta;
};
typedef struct _oz_cone_s oz_cone_t;

/* Bicone: its axis is parallel to the y-axis */
struct _oz_bicone_s {
  ozpfnSupport support;
  oz_real_t _base_radius, _half_height;
  oz_real_t _sintheta;
};
typedef struct _oz_bicone_s oz_bicone_t;

/* Leaf Cylinder */
struct _oz_leaf_cylinder_s {
  ozpfnSupport support;
  oz_real_t _half_width, _half_length, _half_height;
  oz_real_t _circle_radius, _circle_distance;
};
typedef struct _oz_leaf_cylinder_s oz_leaf_cylinder_t;

/* Sphere-swept */
struct _oz_sphere_swept_s {
  ozpfnSupport support;
  const void *_shape;
  oz_real_t _radius;
};
typedef struct _oz_sphere_swept_s oz_sphere_swept_t;

/* Hull - the convex hull of two shapes */
struct _oz_hull_s {
  ozpfnSupport support;
  const void *_ca;     /* convex a's support function */
  const void *_cb;     /* convex b's support function */
  oz_transform_t _ta;  /* transform data of a */
  oz_transform_t _tb;  /* transform data of b */
  /* Cache inverse rotations */
  oz_real_t _inv_rot_a[4];
  oz_real_t _inv_rot_b[4];
};
typedef struct _oz_hull_s oz_hull_t;

/* Minkowski */
struct _oz_minkowski_sum_s {
  ozpfnSupport support;
  const void *_ca;  /* convex a's support function */
  const void *_cb;  /* convex b's support function */
  oz_transform_t _t;  /* transform data */
  /* Cache inverse rotation */
  oz_real_t _inv_rot[4];
};
typedef struct _oz_minkowski_sum_s oz_minkowski_sum_t;


void
ozShpereInit(oz_sphere_t *sphere);

void
ozPointInit(oz_point_t *point);

void
ozLineInit(oz_line_t *line);

void
ozDiskInit(oz_disk_t *disk);

void
ozMeshInit(oz_mesh_t *mesh,
           const oz_uint_t n_vertices,
           const oz_real_t *vertices,
           const oz_uint_t n_faces,
           const oz_uint_t *face_start,
           const oz_uint_t *faces);

void
ozMeshTerminate(oz_mesh_t *mesh);

void
ozCylinderInit(oz_cylinder_t *cylinder,
               oz_real_t base_radius,
               oz_real_t height);

void
ozBoxInit(oz_box_t *box,
          const oz_real_t *extent);

void
ozConeInit(oz_cone_t *cone,
           oz_real_t base_radius,
           oz_real_t height);

void
ozBiconeInit(oz_cone_t *bicone,
             oz_real_t base_radius,
             oz_real_t height);

void
ozLeafCylinderInit(oz_leaf_cylinder_t *leaf,
                   oz_real_t width,
                   oz_real_t length,
                   oz_real_t height);

void
ozSphereSweptInit(oz_sphere_swept_t *swept,
                  oz_real_t radius,
                  const void *shape);

void
ozHullInit(oz_hull_t *hull,
           const oz_transform_t *ta,
           const void *ca,
           const oz_transform_t *tb,
           const void *cb);

void
ozMinkowskiSumInit(oz_minkowski_sum_t *msum,
                   const oz_transform_t *t,
                   const void *ca,
                   const void *cb);


#endif
