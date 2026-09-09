/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 * This is an internal header that must not be installed.
 *
 * A block forest ("blockforest") is a regular grid of uniform blocks,
 * each block being a p8est tree refined to a level such that the block's
 * grid cells coincide with the octants of the tree.  The whole forest
 * covers an axis-aligned bounding box (AABB) given by the configuration.
 * This example is 3D only (p8est / octree).
 */

#ifndef P8EST_BLOCKFOREST_GLOBAL_H
#define P8EST_BLOCKFOREST_GLOBAL_H

#include <sc_options.h>
#include <p8est_geometry.h>

/* one value of "w" is associated with every block */
typedef struct blockforest_block_data {
    double w;
} blockforest_block_data_t;

/* one value of "w" is stored with every octant of the forest */
typedef struct blockforest_quadrant_data {
    double w;
} blockforest_quadrant_data_t;

/* we maintain application data in a global structure and pass it around */
typedef struct p8est_blockforest_global {
    /* global objects */
    sc_MPI_Comm           mpicomm;
    sc_options_t         *options;
    const char           *config_file;
    int                   help;
    int                   mpisize;    /* expected process count (config) */
    int                   mpirank;
    double                dx;         /* grid spacing (config) */
    const char           *min_str;    /* "x0 y0 z0" (config) */
    const char           *max_str;    /* "x1 y1 z1" (config) */
    const char           *block_size_str;     /* "bx by bz" (config) */
    double                min[3];     /* AABB lower corner (parsed) */
    double                max[3];     /* AABB upper corner (parsed) */
    int                   block_size[3];      /* grid cells per block (parsed) */
    double                w;          /* scalar value carried by every block */
    int                   novtk;      /* do not write VTK files */

    /* derived quantities */
    int                   n[3];       /* global number of grid cells */
    int                   block_count[3];     /* number of blocks per axis */
    int                   tree_level;         /* octants per block edge = 2^tree_level */
    int                   min_level;          /* forest refinement level (>= tree_level) */

    /* forest objects */
    p8est_geometry_t     *geom;
    p8est_connectivity_t *conn;
    p8est_t              *p8est;

    /* temporary per-octant w array for VTK output and its counter */
    p4est_locidx_t qcount;
    sc_array_t     *qarray;
} p8est_blockforest_global_t;

/* this function contains the core of the demonstration program */
int p8est_blockforest_run (p8est_blockforest_global_t * g);

#endif /* !P8EST_BLOCKFOREST_GLOBAL_H */
