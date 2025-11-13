<!--
Project-specific Copilot instructions for the schism-esmf repository.
Keep this short and actionable: what an automated coding agent needs to know to be
productive here (architecture, build/test/dev flows, conventions, and integration points).
-->
# schism-esmf — Copilot / automated agent instructions

Quick orientation (big picture)
- This repository provides ESMF/NUOPC "caps" that wrap a pre-built SCHISM hydrodynamics
  library so it can be used as an ESMF/NUOPC component. Key outputs are libraries
  (eg. `libesmf_schism.a`, `libschism_nuopc.a`) and executables `main_esmf` / `main_nuopc`.
- Major directories: `src/` (Fortran sources and subcomponents), `example/` (setups/tests),
  `scripts/` (helpers to create setups), and `cmake/` (custom CMake helpers, notably ESMF discovery).

What to assume when editing/building
- SCHISM is expected to be pre-built. The build system requires `SCHISM_BUILD_DIR` to point
  to a SCHISM build (contains `include/` and `lib/` with `libhydro.a` etc.). See `CMakeLists.txt`.
- ESMF is discovered via an `esmf.mk` file; provide it via the environment variable `ESMFMKFILE`
  (or let `cmake/FindESMF.cmake` locate it). The CMake helper `cmake/FindESMF.cmake` parses the
  ESMF makefile and exposes include/lib variables and an `ESMF::ESMF` target.

Build and test commands (concrete)
- Fast: set env vars, configure, build
  ```bash
  export ESMFMKFILE=/path/to/esmf.mk
  export SCHISM_BUILD_DIR=/path/to/schism/build
  mkdir -p build && cd build
  cmake ..
  cmake --build . -- -j$(nproc)
  ```
- **macOS-specific**: use Homebrew-installed gfortran and explicitly set the Fortran compiler
  ```bash
  export ESMFMKFILE=/path/to/esmf.mk
  export SCHISM_BUILD_DIR=/path/to/schism/build
  mkdir -p build && cd build
  cmake .. \
    -DCMAKE_Fortran_COMPILER=$(which gfortran) \
    -DCMAKE_C_COMPILER=$(which gcc)
  cmake --build . -- -j$(sysctl -n hw.ncpu)
  ```
- Legacy Make: The repo also supports its historical top-level `Makefile` workflow; some docs/examples
  mention using `make all` in the repository root or `make` inside `example/` to exercise runs.
- To run example setups: go to `example/` and run `make` as described by `scripts/ReadMe.md` and `example/Makefile`.

CI-ready checklist (for automated builds)
1. Ensure `ESMFMKFILE` points to a valid `esmf.mk` (check existence: `test -f "$ESMFMKFILE"`).
2. Ensure `SCHISM_BUILD_DIR` points to a pre-built SCHISM (check lib: `test -f "$SCHISM_BUILD_DIR/lib/libhydro.a"`).
3. Set `CMAKE_Fortran_COMPILER` to an MPI-aware wrapper (e.g., `mpif90`, `esmpifort`) — CMake will warn if you use a bare compiler.
4. Run CMake configure with explicit compiler flags and capture stderr for "FATAL_ERROR" or "WARNING" messages.
5. Run `cmake --build . -- VERBOSE=1` to see full compile/link commands and catch missing symbols early.
6. Optional: run a quick smoke test by executing `./main_esmf --help` or similar (if applicable) to verify linkage.

**Automated validation**: Run `.github/test-cmake-instructions.sh` to validate your environment against this checklist.
See `.github/TEST_README.md` for details.

**Build troubleshooting**: If build fails, see `.github/BUILD_TROUBLESHOOTING.md` for detailed solutions to common issues:
- MPI/ParMETIS symbol mismatches
- Missing module files
- MOSSCO dependencies
- Compile option parsing errors

Important patterns & conventions
- Fortran modules are written in `.F90` files and CMake emits them into a centralized
  modules directory (`CMAKE_Fortran_MODULE_DIRECTORY` configured to `build/modules`). When adding
  new modules, ensure consumers see the `target_include_directories(...)` or that targets are linked
  in the correct CMake subdirectory (see `src/CMakeLists.txt`).
- SCHISM libraries expected in `${SCHISM_BUILD_DIR}/lib` with names found by `find_library` for `core` and `hydro`.
- Optional integrations: PDAF (data assimilation) and MOSSCO are optional — guard changes so they compile
  when absent. See `scripts/Readme.md` and `src/` subdirs referencing PDAF and MOSSCO.

Where to make changes
- High-level code: `src/` (subdirs `driver/`, `model/`, `schism/`). Add tests or example inputs under `example/`.
- Build logic and ESMF wiring: `cmake/` and top-level `CMakeLists.txt`. Small edits to linking are usually
  localized in `src/CMakeLists.txt` or `src/schism/CMakeLists.txt`.

Debugging and verification tips
- Build with verbose output: `cmake --build . -- VERBOSE=1` or `make VERBOSE=1` to see compiler/linker commands.
- If ESMF or SCHISM symbols are missing at link time, verify `ESMFMKFILE` and `SCHISM_BUILD_DIR` and inspect
  `cmake/FindESMF.cmake` behavior (it exports `ESMF_INCLUDE_DIRS` and `ESMF_LIBRARIES`).
- Module (.mod) problems: ensure `CMAKE_Fortran_MODULE_DIRECTORY` is consistent and that targets depend on the target
  that creates the module (use `target_link_libraries` and `target_include_directories` instead of global include hacks).

Examples from the codebase
- See `ReadMe.md` at project root for env var examples (ESMFMKFILE and SCHISM_BUILD_DIR).  
- See `cmake/FindESMF.cmake` for how ESMF variables and the `ESMF::ESMF` target are produced.
- See `src/CMakeLists.txt` for how `main_esmf` and `main_nuopc` are declared and linked.

Do not assume
- That SCHISM or ESMF will be built by this repository — they are external prerequisites.  
- That all example workflows are unit-tested; many are manual and exercised via `example/` and `scripts/`.

When in doubt
- Prefer minimal, scoped changes: update or add a CMake target in the specific `src/` subdirectory rather than global edits.
- Add small runnable examples under `example/` when introducing behavior changes so CI or local runs can validate them.

If anything here is unclear or you'd like me to expand a section (build matrix, tests, or common failure modes), tell me which part and I'll iterate.

### Style to consider

Always use camelcase for naming the ReadMe.md files in each directory.