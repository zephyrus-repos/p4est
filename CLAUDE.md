# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

This is a personal study/scratch repository for the [p4est](https://www.p4est.org/) library. The source of
p4est (v2.8.7, tarball under `p4est-2.8.7.tar.gz`, unpacked into `tmp/`) is the upstream release as
shipped from p4est.org, and `examples/` is a copy of the upstream `example/` sources (plus a custom
CMake build) used as a testbed for learning p4est. There is no product code here.

## Build & run

p4est is prebuilt and installed at `/opt/p4est` (static libs `libp4est.a`, `libsc.a`, headers under
`/opt/p4est/include`). MPI is OpenMPI 4.1.8 at `/opt/openmpi-4.1.8`.

Configure and build the examples with CMake (out-of-source, from the repo root):

```bash
cmake -S examples -B examples/build
cmake --build examples/build -j$(nproc)
```

- `examples/CMakeLists.txt` is **organised by subdirectory**: it defines one reusable
  `p4est_example(<target> <sources...>)` function and then `add_subdirectory()`s every example
  directory (`balance/`, `mesh/`, ... `userdata/`). Each directory's `CMakeLists.txt` mirrors that
  directory's upstream `Makefile.am` and lists its binaries explicitly.
- Target names are the upstream program names (e.g. `p4est_simple`, `p8est_mesh`, `p4est_step1`),
  **not** source-file names. Binaries land in `examples/build/bin/`. The `p4est_example` function
  handles includes and linking, so a new program is added by calling it from the right directory.
- Helper sources are listed explicitly as part of their program's target rather than being
  auto-excluded from a glob: `spheres/p4est_spheres.c`+`p8est_spheres.c` are compiled into
  `p4est_spheres`/`p8est_spheres`, and `userdata/userdata_run2.c`+`run3.c` into
  `p4est_userdata`/`p8est_userdata` (upstream instead compiles the spheres helpers into `libp4est`;
  here they are folded into the examples because `/opt/p4est` is a prebuilt static lib).
- The hardcoded `P4EST_PREFIX` in `examples/CMakeLists.txt` is `/opt/p4est`. If p4est is ever
  installed elsewhere, that line (and `compile_flags.txt`, used by editor tooling / clangd) must be
  updated.

Run any example under `mpirun`:

```bash
mpirun -np 4 ./examples/build/bin/p4est_simple unit 3
```

There are no unit tests in this repo. Some examples print a forest checksum and self-check it
against a hardcoded regression table (see `simple/simple2.c`); these fail with `SC_CHECK_ABORT`
if the checksum diverges. To verify a rebuild is correct, run `p4est_simple`/`p8est_simple` and
check the `Checksum regression OK` line. (Upstream runs `make check` inside the p4est source build
dir; that is a test of the library, not this repo — see `notes/install.md`.)

## Layout & architecture

### Repo top level

- `examples/` — the code actually developed here. Source is in subdirectories matching upstream
  (`balance/`, `mesh/`, `particles/`, `points/`, `search/`, `simple/`, `spheres/`, `steps/`,
  `tetgen/`, `timings/`, `userdata/`), one `.c` per program with a `2`/`3` suffix for
  2D-quadtree (p4est) vs 3D-octree (p8est) variants. `build/` is the CMake build dir (gitignored).
- `notes/install.md` — how p4est was built and installed from the tarball (Autotools flow:
  `../configure --enable-mpi --disable-shared CFLAGS=... LIBS="-lm -lz -ljansson"`, `make`, `make
  install`). Rebuild/install the library, not the examples, per this doc.
- `tmp/p4est-2.8.7/` — unpacked upstream p4est source. Its `src/` holds the library and its
  `example/` dir is the unmodified origin of `examples/`. Useful for reading p4est internals and
  the upstream CMake files.
- `p4est-2.8.7.tar.gz` — the official release tarball, committed as-is.

### How every example program is structured

All examples follow one skeleton, regardless of subdirectory:

1. **Init MPI & p4est** — `sc_MPI_Init`, get `mpisize`/`mpirank`, then
   `sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT)` and `p4est_init (NULL, SC_LP_DEFAULT)`.
2. **Parse args** — a config/domain keyword plus a refinement level (see the `Usage:` comment at the
   top of each file); several take optional radius/dimension args.
3. **Build connectivity** — `p4est_connectivity_new_*()` (brick, moebius, cubed, disk, icosahedron,
   etc.). Geometry (for non-planar domains) via `p4est_geometry_new_*()`, attached with
   `p4est_new_ext (..., geom)`.
4. **Run forest operations** — the point of each example, e.g. `p4est_refine`/`p4est_coarsen`,
   `p4est_balance`, `p4est_partition`, `p4est_ghost_new`, `p4est_search_*`, mesh building, point
   finding, VTK output (`p4est_vtk_write_file`). Refinement/coarsening/init callbacks are static
   functions of signature `(p4est_t *, p4est_topidx_t, p4est_quadrant_t *, ...)`.
5. **Teardown** — `p4est_destroy`, `p4est_geometry_destroy`, `p4est_connectivity_destroy`,
   `sc_finalize`, `sc_MPI_Finalize`.

The `steps/` directory is a tutorial sequence: `p4est_step1.c` … `p4est_step5.c` build on each
other (step N includes the code of the previous steps), with `hw32.h` a generated header embedding a
bitmap used as the step-1 refinement criterion. Good starting point for understanding the API.

### Key conventions

- **2D vs 3D duality**: sources named `*2` use the p4est (quadtree) API, `*3` use p8est (octree).
  p8est programs are typically generated from the 2D ones by including `p4est_to_p8est.h`, which
  `#define`s p4est symbols to their p8est counterparts (see `spheres/p4est_to_p8est_spheres.h` and
  `userdata/userdata_global.h`). When modifying an example, check whether a matching `*2`/`*3` pair
  should change together.
- **User data**: per-quadrant data is attached via `quadrant->p.user_data` (set by the init
  callback and allocated through `p4est_new_ext`'s `data_size` argument); the `userdata/` examples
  demonstrate this, including internal vs external storage modes.
- **Error handling**: p4est/sc code uses `SC_CHECK_*` macros (abort on failure) and
  `P4EST_GLOBAL_*` / `P4EST_GLOBAL_PRODUCTIONF` for rank-0-printed messages; new code should follow
  the same macros. `p4estindent` (the upstream indentation script) enforces the GNU-style format
  visible in all sources.
- **Debug builds** of the library enable `P4EST_ENABLE_DEBUG`, under which some examples run extra
  assertions (e.g. `simple2` rebalances and asserts an unchanged checksum). The installed `/opt/p4est`
  build is a Release build (`-O2`), so those paths are compiled out.
