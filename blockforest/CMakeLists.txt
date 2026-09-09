# CMakeLists.txt — top-level build for the blockforest program.
#
# The blockforest program lives at the repository root (not under examples/)
# as its own standalone build. It is a 3D-only p8est program, built against
# the preinstalled p4est at /opt/p4est (static libs, see notes/install.md).
cmake_minimum_required(VERSION 3.20)
project(p4est_blockforest C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# MPI and threads
find_package(MPI REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

# p4est installation
set(P4EST_PREFIX "/opt/p4est")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(p8est_blockforest
    blockforest3.c
    blockforest_run.c
    blockforest_global.h
    blockforest_config.ini)
target_include_directories(p8est_blockforest PRIVATE
    "${P4EST_PREFIX}/include")
target_link_libraries(p8est_blockforest PRIVATE
    "${P4EST_PREFIX}/lib/libp4est.a"
    "${P4EST_PREFIX}/lib/libsc.a"
    jansson
    z
    m
    MPI::MPI_C
    Threads::Threads)
