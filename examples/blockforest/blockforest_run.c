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

/* the core functionality of the p8est blockforest example program */
#include <p4est_base.h>
#include "blockforest_global.h"
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

/***************** mapping between blocks, octants, and the AABB **************/

/* smallest level such that 2^level >= a */
static int
blockforest_ceillog2 (int a)
{
    int level = 0;
    int two_power = 1;

    P4EST_ASSERT (a > 0);
    while (two_power < a) {
        two_power *= 2;
        ++level;
    }
    return level;
}

/* copy the per-octant value w into a contiguous array for VTK output */
static void
blockforest_vtk_volume (p8est_iter_volume_info_t * v, void *user_data)
{
    p8est_blockforest_global_t *g =
        (p8est_blockforest_global_t *) user_data;
    blockforest_quadrant_data_t *qdat;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (v != NULL);
    P4EST_ASSERT (v->p8est == g->p8est);

    qdat = (blockforest_quadrant_data_t *) v->quad->p.user_data;
    P4EST_ASSERT (qdat != NULL);
    P4EST_ASSERT (g->qcount < g->p8est->local_num_quadrants);

    /* write octant value into the output array and advance counter */
    *(double *) sc_array_index (g->qarray, g->qcount++) = qdat->w;
}

/* write a VTK file with the mesh and the per-octant w field */
static int
blockforest_vtk_general (p8est_blockforest_global_t * g, const char *filename)
{
    const char          *fnames[1] = {"w"};
    sc_array_t          *fvalues[1] = {NULL};
    p8est_vtk_context_t *vtk, *rvtk;
    p4est_locidx_t       qcount;
    p4est_locidx_t       it;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->p8est != NULL);
    P4EST_ASSERT (g->qarray != NULL);
    P4EST_ASSERT (g->qarray->elem_size == sizeof (double));

    /* the cell data array must hold one value per local octant */
    qcount = g->p8est->local_num_quadrants;
    sc_array_resize (g->qarray, (size_t) qcount);
    for (it = 0; it < qcount; ++it) {
        *(double *) sc_array_index (g->qarray, it) = 0.;
    }

    /* fill the array from the octant user data */
    P4EST_ASSERT (g->qcount == 0);
    p8est_iterate (g->p8est, NULL, g, blockforest_vtk_volume, NULL,
                   NULL, NULL);
    P4EST_ASSERT (g->qcount == g->p8est->local_num_quadrants);
    g->qcount = 0;

    fvalues[0] = g->qarray;

    /* write file in multiple steps */
    P4EST_ASSERT (g->p8est != NULL);
    vtk = p8est_vtk_context_new (g->p8est, filename);
    p8est_vtk_context_set_geom (vtk, g->geom);
    if ((rvtk = p8est_vtk_write_header (vtk)) == NULL) {
        P4EST_GLOBAL_LERRORF ("ERROR: write VTK header for %s\n", filename);
        return -1;
    }
    P4EST_ASSERT (rvtk == vtk);

    /* the cell data is written from a contiguous array of values */
    if ((rvtk = p8est_vtk_write_cell_data (vtk, 1, 1, 1, 0, 1, 0,
                                           fnames, fvalues)) == NULL) {
        P4EST_GLOBAL_LERRORF ("ERROR: write VTK data for %s\n", filename);
        return -1;
    }
    P4EST_ASSERT (rvtk == vtk);

    /* finalize the output files and deallocate context */
    if (p8est_vtk_write_footer (vtk)) {
        P4EST_GLOBAL_LERRORF ("ERROR: write VTK footer for %s\n", filename);
        return -1;
    }

    /* return success */
    return 0;
}

/* write a VTK file with the per-octant w field */
static int
blockforest_vtk (p8est_blockforest_global_t * g, const char *filename)
{
    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->p8est != NULL);
    P4EST_ASSERT (g->qarray == NULL);
    if (g->novtk) {
        /* output disabled */
        return 0;
    }

    /* allocate the temporary per-octant w array */
    g->qarray = sc_array_new_count (sizeof (double),
                                    (size_t) g->p8est->local_num_quadrants);
    if (blockforest_vtk_general (g, filename)) {
        P4EST_GLOBAL_LERRORF ("ERROR: write VTK file %s\n", filename);
        sc_array_destroy_null (&g->qarray);
        return -1;
    }

    /* return memory neutral */
    sc_array_destroy_null (&g->qarray);
    return 0;
}

/* provide function for consistent deallocation */
static int
blockforest_run_return (int retval, p8est_blockforest_global_t * g)
{
    P4EST_ASSERT (g != NULL);

    /* clean up what has been allocated in the same function */
    if (g->p8est != NULL) {
        /* destroy forest */
        p8est_destroy (g->p8est);
        g->p8est = NULL;
    }
    return retval;
}

/* callback to initialize per-octant user data */
static void
blockforest_init (p8est_t * p8est, p4est_topidx_t which_tree,
                  p8est_quadrant_t * quadrant)
{
    p8est_blockforest_global_t *g =
        (p8est_blockforest_global_t *) p8est->user_pointer;
    blockforest_quadrant_data_t *qdat;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->conn != NULL);

    /* p8est is agnostic to the quadrant user data */
    qdat = (blockforest_quadrant_data_t *) quadrant->p.user_data;
    P4EST_ASSERT (qdat != NULL);

    /* every octant of a block carries the w value of its block, which is
       stored as a per-tree attribute of the connectivity */
    if (g->conn->tree_attr_bytes == sizeof (blockforest_block_data_t)) {
        blockforest_block_data_t *bdat =
            (blockforest_block_data_t *) (g->conn->tree_to_attr +
                                          which_tree * g->conn->tree_attr_bytes);
        qdat->w = bdat->w;
    } else {
        /* fallback in case attributes are absent */
        qdat->w = g->w;
    }
}

/* custom geometry mapping the reference cube of every tree to its AABB cell.
   Tree t of the brick connectivity is the t-th block in Morton order; we
   derive its (bx, by, bz) index, then map the reference cube [0,1]^3 to
   the block's cell inside the AABB. */
static void
blockforest_aabb_X (p8est_geometry_t * geom, p4est_topidx_t which_tree,
                    const double abc[3], double xyz[3])
{
    p8est_blockforest_global_t *g = (p8est_blockforest_global_t *) geom->user;
    int                          bx, by, bz;
    int                          nx, ny, nz;
    p4est_topidx_t               num_blocks;
    double                       x0, y0, z0;

    P4EST_ASSERT (g != NULL);

    num_blocks = (p4est_topidx_t) ((size_t) g->block_count[0] *
                                   (size_t) g->block_count[1] *
                                   (size_t) g->block_count[2]);
    P4EST_ASSERT (which_tree < num_blocks);

    /* block index within the Morton-ordered brick of trees */
    bz = (int) (which_tree % g->block_count[2]);
    by = (int) ((which_tree / g->block_count[2]) % g->block_count[1]);
    bx = (int) (which_tree / ((p4est_topidx_t) g->block_count[1] *
                              (p4est_topidx_t) g->block_count[2]));

    /* grid cells per block; the tree spans block_size grid cells */
    nx = g->block_size[0];
    ny = g->block_size[1];
    nz = g->block_size[2];

    x0 = g->min[0] + ((double) (bx * nx)) * g->dx;
    y0 = g->min[1] + ((double) (by * ny)) * g->dx;
    z0 = g->min[2] + ((double) (bz * nz)) * g->dx;

    /* linear mapping from reference cube to the block's cell */
    xyz[0] = x0 + abc[0] * ((double) nx) * g->dx;
    xyz[1] = y0 + abc[1] * ((double) ny) * g->dx;
    xyz[2] = z0 + abc[2] * ((double) nz) * g->dx;
}

/* attach the w value to every tree of the connectivity as a block attribute */
static void
blockforest_set_block_w (p8est_blockforest_global_t * g)
{
    p4est_topidx_t it;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->conn != NULL);

    /* if the connectivity already has attributes, make sure they have the
       requested size; otherwise allocate them */
    if (g->conn->tree_attr_bytes == 0) {
        p8est_connectivity_set_attr (g->conn, sizeof (blockforest_block_data_t));
    } else {
        P4EST_ASSERT (g->conn->tree_attr_bytes == sizeof (blockforest_block_data_t));
    }

    /* assign one w value per block */
    for (it = 0; it < g->conn->num_trees; ++it) {
        blockforest_block_data_t *bdat =
            (blockforest_block_data_t *) (g->conn->tree_to_attr +
                                          it * g->conn->tree_attr_bytes);
        bdat->w = g->w;
    }
}

/* compute derived quantities and create the forest */
int
p8est_blockforest_run (p8est_blockforest_global_t * g)
{
    int erres = 0;
    int d;
    int nx, ny, nz;
    p4est_locidx_t min_quadrants;
    int min_level;
    int tree_level;
    p8est_geometry_t *geom;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->conn == NULL);
    P4EST_ASSERT (g->geom == NULL);
    P4EST_ASSERT (g->p8est == NULL);

    /* validate config: AABB and dx determine the global grid */
    P4EST_ASSERT (g->dx > 0.);
    for (d = 0; d < 3; ++d) {
        P4EST_ASSERT (g->min[d] < g->max[d]);
        P4EST_ASSERT (g->block_size[d] > 0);
        P4EST_ASSERT ((g->max[d] - g->min[d]) > 0.);

        /* number of global grid cells per axis */
        g->n[d] = (int) (floor ((g->max[d] - g->min[d]) / g->dx + .5));
        P4EST_ASSERT (g->n[d] > 0);

        /* block count per axis; the domain must be covered by whole blocks */
        P4EST_ASSERT (g->n[d] % g->block_size[d] == 0);
        g->block_count[d] = g->n[d] / g->block_size[d];
    }

    nx = g->n[0];
    ny = g->n[1];
    nz = g->n[2];

    /* the block edge length in grid cells is block_size[d] = 2^tree_level */
    tree_level = blockforest_ceillog2 (g->block_size[0]);
    P4EST_ASSERT (tree_level >= 0);
    P4EST_ASSERT (blockforest_ceillog2 (g->block_size[0]) ==
                  blockforest_ceillog2 (g->block_size[1]));
    P4EST_ASSERT (blockforest_ceillog2 (g->block_size[0]) ==
                  blockforest_ceillog2 (g->block_size[2]));
    P4EST_ASSERT (tree_level <= P8EST_MAXLEVEL);
    /* one octant must coincide with exactly one grid cell, so the block edge
       length must be a power of two */
    P4EST_ASSERT ((1 << tree_level) == g->block_size[0]);
    P4EST_ASSERT ((1 << tree_level) == g->block_size[1]);
    P4EST_ASSERT ((1 << tree_level) == g->block_size[2]);
    g->tree_level = tree_level;

    /* each tree spans exactly block_size grid cells, so refining it uniformly
       to tree_level makes its octants coincide one-to-one with the grid cells;
       the global grid is tiled by block_count trees, not resolved inside one */
    min_level = tree_level;
    P4EST_ASSERT (min_level <= P8EST_MAXLEVEL);
    g->min_level = min_level;

    /* print the derived configuration on rank 0 */
    P4EST_GLOBAL_PRODUCTIONF
        ("blockforest: domain [%g, %g]x[%g, %g]x[%g, %g], dx %g\n",
         g->min[0], g->max[0], g->min[1], g->max[1], g->min[2], g->max[2],
         g->dx);
    P4EST_GLOBAL_PRODUCTIONF
        ("blockforest: global grid %dx%dx%d, block grid %dx%dx%d, "
         "block size %dx%dx%d, tree level %d, min level %d\n",
         nx, ny, nz, g->block_count[0], g->block_count[1],
         g->block_count[2], g->block_size[0], g->block_size[1],
         g->block_size[2], tree_level, min_level);
    P4EST_GLOBAL_PRODUCTIONF
        ("blockforest: running with %d processes, w = %g\n",
         g->mpisize, g->w);

    /* create the connectivity: one tree per block */
    g->conn = p8est_connectivity_new_brick (g->block_count[0],
                                            g->block_count[1],
                                            g->block_count[2], 0, 0, 0);
    P4EST_ASSERT (g->conn != NULL);

    /* attach the block attribute w to every tree */
    blockforest_set_block_w (g);

    /* create a geometry that maps the reference cube of every tree to the
       corresponding cell of the AABB.  The geometry holds the global config
       in its user pointer; it is destroyed by p8est_geometry_destroy. */
    geom = P4EST_ALLOC (p8est_geometry_t, 1);
    geom->name = "blockforest_aabb";
    geom->user = g;
    geom->X = blockforest_aabb_X;
    geom->destroy = NULL;
    g->geom = geom;

    /* create the forest: one uniform tree per block, carrying w as user data.
       min_quadrants is set to a safe minimum so every process obtains at
       least a share of the octants. */
    min_quadrants = 8;
    g->p8est = p8est_new_ext (g->mpicomm, g->conn, min_quadrants,
                              min_level, 1, sizeof (blockforest_quadrant_data_t),
                              blockforest_init, g);
    P4EST_ASSERT (g->p8est != NULL);

    /* distribute the blocks (trees) across the processes */
    p8est_partition (g->p8est, 1, NULL);

    /* write VTK files to visualize the block forest */
    if (blockforest_vtk (g, "blockforest")) {
        P4EST_GLOBAL_LERROR ("ERROR: write VTK output after forest_new\n");
        return blockforest_run_return (-1, g);
    }

    /* return memory neutral */
    return blockforest_run_return (erres, g);
}
