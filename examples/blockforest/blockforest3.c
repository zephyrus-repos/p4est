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
 * This example program generates a block forest from a configuration file.
 *
 *   p8est_blockforest [<OPTIONS>]
 *
 * The following options are recognized:
 *   --help          Display a usage and help message and exit successfully.
 *   --config FILE   Read the configuration from the INI file FILE.
 *                   A configuration file contains the AABB domain (min/max),
 *                   the grid spacing dx, the number of grid cells per block,
 *                   the expected MPI process count, and a scalar value w
 *                   that is assigned to every block.
 *
 * Invalid options or arguments result in an error message and exit status.
 */
static const char *p8est_blockforest_usage =
    "The configuration file (INI) is read via sc_options_load_ini.\n"
    "  The following options are legal (default values in parentheses):\n"
    "  o dx = <double>             Grid spacing; this determines the global\n"
    "                              number of grid cells together with the AABB.\n"
    "  o min = <x0> <y0> <z0>      Lower corner of the AABB (0 0 0).\n"
    "  o max = <x1> <y1> <z1>      Upper corner of the AABB (1 1 1).\n"
    "  o block_size = <bx> <by> <bz>  Grid cells per block (4 4 4).\n"
    "  o num_procs = <int>         Expected MPI process count (must equal\n"
    "                              the actual number of processes).\n"
    "  o w = <double>              Scalar value carried by every block (1.0).\n"
    "  o no_vtk = <bool>           Do not write VTK files (false).\n";

#include <p4est_base.h>
#include "blockforest_global.h"

/* process the command line and the configuration file */
static int
p8est_blockforest_process (p8est_blockforest_global_t * g)
{
    int erres = 0;
    int mpiret;
    int mpisize = 0;

    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->options != NULL);

    /* parse the three space-separated coordinate triples */
    if (sscanf (g->min_str, "%lf %lf %lf", &g->min[0], &g->min[1],
                &g->min[2]) != 3) {
        P4EST_GLOBAL_LERROR ("ERROR: min must be three numbers x0 y0 z0\n");
        erres = -1;
    }
    if (!erres && sscanf (g->max_str, "%lf %lf %lf", &g->max[0], &g->max[1],
                          &g->max[2]) != 3) {
        P4EST_GLOBAL_LERROR ("ERROR: max must be three numbers x1 y1 z1\n");
        erres = -1;
    }
    if (!erres && sscanf (g->block_size_str, "%d %d %d", &g->block_size[0],
                          &g->block_size[1], &g->block_size[2]) != 3) {
        P4EST_GLOBAL_LERROR
            ("ERROR: block_size must be three integers bx by bz\n");
        erres = -1;
    }

    /* check the expected process count against the actual one */
    if (!erres) {
        mpiret = sc_MPI_Comm_size (g->mpicomm, &mpisize);
        SC_CHECK_MPI (mpiret);
        if (g->mpisize != mpisize) {
            P4EST_GLOBAL_LERRORF
                ("ERROR: expected %d processes but %d are running\n",
                 g->mpisize, mpisize);
            erres = -1;
        }
    }

    /* the config values must satisfy the block forest construction */
    if (!erres && !(g->dx > 0.)) {
        P4EST_GLOBAL_LERROR ("ERROR: dx must be positive\n");
        erres = -1;
    }

    /* return the error status */
    return erres;
}

/* process the command line options */
static int
p8est_blockforest_options (int argc, char **argv,
                           p8est_blockforest_global_t * g)
{
    int           erres;
    int           firstarg;
    sc_options_t *o;

    P4EST_ASSERT (argc >= 1);
    P4EST_ASSERT (argv != NULL);
    P4EST_ASSERT (g != NULL);
    P4EST_ASSERT (g->options == NULL);

    /* allocate new options processor */
    o = g->options = sc_options_new (argv[0]);

    /* register the configuration variables */
    sc_options_add_switch (o, 'h', "help", &g->help,
                           "Print help message and exit cleanly");
    sc_options_add_string (o, 'c', "config", &g->config_file, NULL,
                           "Configuration file in INI format");
    sc_options_add_double (o, 'd', "dx", &g->dx, .125,
                           "Grid spacing (determines the global grid)");
    sc_options_add_string (o, 'm', "min", &g->min_str, "0 0 0",
                           "Lower corner x0 y0 z0 of the AABB");
    sc_options_add_string (o, 'M', "max", &g->max_str, "1 1 1",
                           "Upper corner x1 y1 z1 of the AABB");
    sc_options_add_string (o, 'b', "block_size", &g->block_size_str, "4 4 4",
                           "Grid cells per block bx by bz");
    sc_options_add_int (o, 'p', "num_procs", &g->mpisize, 1,
                        "Expected MPI process count");
    sc_options_add_double (o, 'w', "w", &g->w, 1.,
                           "Scalar value carried by every block");
    sc_options_add_switch (o, 'n', "no-vtk", &g->novtk,
                           "Do not write VTK graphics files");

    /* process the command line options */
    erres = 0;
    if (!erres && (erres = ((firstarg =
                             sc_options_parse (p4est_get_package_id (),
                                               SC_LP_ERROR, o, argc, argv)) < 0))) {
        P4EST_GLOBAL_LERROR ("ERROR: processing options\n");
    }
    P4EST_ASSERT (erres || (1 <= firstarg && firstarg <= argc));

    /* load the configuration file if given; it overwrites the options */
    if (!erres && g->config_file != NULL) {
        if (sc_options_load_ini (p4est_get_package_id (), SC_LP_ERROR, o,
                                 g->config_file, NULL)) {
            P4EST_GLOBAL_LERRORF ("ERROR: loading configuration file %s\n",
                                  g->config_file);
            erres = -1;
        }
    }

    /* validate the configuration */
    if (!erres && (erres = p8est_blockforest_process (g))) {
        P4EST_GLOBAL_LERROR ("ERROR: processing the configuration\n");
    }

    if (g->help || erres) {
        /* print a usage message to explain the command line */
        sc_options_print_usage (p4est_get_package_id (), SC_LP_PRODUCTION, o,
                                p8est_blockforest_usage);
    } else {
        /* on normal operation print options for posteriority */
        sc_options_print_summary (p4est_get_package_id (), SC_LP_PRODUCTION,
                                  o);
    }

    /* this function has processed the command line */
    return erres;
}

/* free allocated application memory */
static void
p8est_blockforest_cleanup (p8est_blockforest_global_t * g)
{
    P4EST_ASSERT (g != NULL);

    /* this data may or may not have been initialized */
    if (g->geom != NULL) {
        p8est_geometry_destroy (g->geom);
        g->geom = NULL;
    }
    if (g->conn != NULL) {
        p8est_connectivity_destroy (g->conn);
        g->conn = NULL;
    }

    /* options are always defined at this point */
    if (g->options != NULL) {
        sc_options_destroy (g->options);
        g->options = NULL;
    }
}

/* the main function of the program */
int
main (int argc, char **argv)
{
    int                         erres;
    int                         mpiret;
    p8est_blockforest_global_t  sglobal, *global = &sglobal;

    /* initialize MPI subsystem */
    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);

    /* initialize global application data context */
    memset (global, 0, sizeof (*global));
    global->mpicomm = sc_MPI_COMM_WORLD;

    /* this sets an abort handler and initializes the log system */
    sc_init (global->mpicomm, 1, 1, NULL, SC_LP_APPLICATION);
    p4est_init (NULL, SC_LP_APPLICATION);

    /* we identify an error status of the program with a nonzero value */
    erres = 0;

    /* process command line options */
    if (!erres && (erres = p8est_blockforest_options (argc, argv, global))) {
        P4EST_GLOBAL_LERROR ("ERROR: Usage/options\n");
    }

    /*
     * Run the actual demo (except when a help message has been requested).
     * We have moved the code for this function into a separate file.
     * The reason is that the present file is an excellent template
     * for your own p4est application.  Just copy it and hack away.
     */
    if (!erres && !global->help && (erres = p8est_blockforest_run (global))) {
        P4EST_GLOBAL_LERROR ("ERROR: running the program\n");
    }

    /* free all data regardless of the error condition */
    p8est_blockforest_cleanup (global);

    /* check memory balance and clean up internal registrations */
    sc_finalize ();

    /* release the MPI subsystem */
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    /* return failure or success to the calling shell */
    return erres ? EXIT_FAILURE : EXIT_SUCCESS;
}
