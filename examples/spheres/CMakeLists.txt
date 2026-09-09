# Mirrors example/spheres/Makefile.am.
#
# Note: upstream compiles p4est_spheres.c / p8est_spheres.c into the p4est
# library itself; since this repo links against the prebuilt static
# /opt/p4est, they are compiled into the example executables instead.
p4est_example(p4est_spheres spheres_global.h p4est_spheres.c spheres2.c)
p4est_example(p8est_spheres spheres_global.h p8est_spheres.c spheres3.c)
