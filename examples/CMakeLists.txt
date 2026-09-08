cmake_minimum_required(VERSION 3.12)
project(p4est_examples C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# MPI and threads
find_package(MPI REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

# p4est installation
set(P4EST_PREFIX "/opt/p4est")

# Find all example sources
file(GLOB SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*/*.c"
)

# Exclude auxiliary sources
list(REMOVE_ITEM SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/spheres/p4est_spheres.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/spheres/p8est_spheres.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/userdata/userdata_run2.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/userdata/userdata_run3.c"
)

# Build one executable per source
foreach(src IN LISTS SOURCES)
    get_filename_component(name "${src}" NAME_WE)
    set(exe_sources "${src}")

    # Add userdata helpers
    if(name STREQUAL "userdata2")
        list(APPEND exe_sources
            "${CMAKE_CURRENT_SOURCE_DIR}/userdata/userdata_run2.c")
    elseif(name STREQUAL "userdata3")
        list(APPEND exe_sources
            "${CMAKE_CURRENT_SOURCE_DIR}/userdata/userdata_run3.c")
    endif()

    add_executable(${name} ${exe_sources})

    target_include_directories(${name} PRIVATE
        "${P4EST_PREFIX}/include"
    )

    # Link p4est and its dependencies
    target_link_libraries(${name} PRIVATE
        "${P4EST_PREFIX}/lib/libp4est.a"
        "${P4EST_PREFIX}/lib/libsc.a"
        jansson
        z
        m
        MPI::MPI_C
        Threads::Threads
    )

endforeach()