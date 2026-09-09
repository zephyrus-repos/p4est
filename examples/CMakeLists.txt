# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(p4est_examples C)

# export compile_commands.json for clang tooling
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# MPI and threads
find_package(MPI REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

# p4est installation
set(P4EST_PREFIX "/opt/p4est")

# All example binaries are placed into ${CMAKE_BINARY_DIR}/bin and named
# after the upstream example names (see each subdir's Makefile.am), e.g.
# p4est_simple, p8est_mesh, p4est_step1.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Create one example executable from the given source file(s).
#   p4est_example(<target> <source1> [<source2> ...])
function(p4est_example target)
    add_executable(${target} ${ARGN})
    target_include_directories(${target} PRIVATE
        "${P4EST_PREFIX}/include")
    target_link_libraries(${target} PRIVATE
        "${P4EST_PREFIX}/lib/libp4est.a"
        "${P4EST_PREFIX}/lib/libsc.a"
        jansson
        z
        m
        MPI::MPI_C
        Threads::Threads)
endfunction()

# Each example directory declares its own binaries. The targets and source
# sets mirror the upstream example/<dir>/Makefile.am files: a *_2.c source
# builds the p4est (2D quadtree) program, a *_3.c source the p8est (3D
# octree) program.
add_subdirectory(balance)
add_subdirectory(blockforest)
add_subdirectory(mesh)
add_subdirectory(particles)
add_subdirectory(points)
add_subdirectory(search)
add_subdirectory(simple)
add_subdirectory(spheres)
add_subdirectory(steps)
add_subdirectory(tetgen)
add_subdirectory(timings)
add_subdirectory(userdata)
