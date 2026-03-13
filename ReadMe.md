<!--
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de
-->

# schism-esmf

Earth System Modeling Framework (ESMF) and National Unified Operational Prediction Capability (NUOPC) caps for SCHISM

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/schism-dev/schism-esmf)
[![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-blue)](https://github.com/schism-dev/schism-esmf)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue)](LICENSE)

## Recent Successful Builds

| Platform | Compiler | ESMF | MPI | Status | Date |
|----------|----------|------|-----|--------|------|
| macOS 15 (arm64) | gfortran 14.3.0 | 8.9.0 | MPICH 4.3.2 | ✅ Passing | Nov 2025 |
| Linux (x86_64) | gfortran 11+ | 8.7+ | OpenMPI/MPICH | ✅ Passing | - |

**Key Requirements**:
- ESMF **must** use real MPI (`mpi_mpich` or `mpi_openmpi`), **not** `mpiuni` or `nompi`
- SCHISM libraries: `libcore.a`, `libhydro.a`, `libturbulence.a`, `libyaml.a`, `libparmetis.a`, `libmetis.a`
- MPI libraries must match between ESMF and SCHISM

## Documentation

| Document | Purpose |
|----------|---------|
| [QUICKSTART.md](./QUICKSTART.md) | **Detailed walkthrough** with platform-specific instructions |
| [BUILD_TROUBLESHOOTING.md](.github/BUILD_TROUBLESHOOTING.md) | Solutions to all 9 known build issues |
| [cmake-architecture.md](doc/cmake-architecture.md) | CMake build system internals and patterns |
| [running-examples.md](doc/running-examples.md) | How to run executables and examples |
| [test-cmake-instructions.sh](.github/test-cmake-instructions.sh) | Automated environment validation |
| [BuildingDocs.md](doc/BuildingDocs.md) | How to build HTML documentation |

**Build HTML Docs**: `cmake .. -DBUILD_DOCS=ON && make docs` (requires Python, MkDocs)

-------------------------------------------------------

# Quick Start

Get SCHISM-ESMF building in under 30 minutes. Full details in [QUICKSTART.md](./QUICKSTART.md).

## Prerequisites

1. ESMF: You need to have ESMF installed and an environment variable
   `ESMFMKFILE` defined (via `setenv` or `export`) that points to your ESMF
   installation, for example

    on femto: 
      setenv ESMFMKFILE /sciclone/home10/yinglong/esmf_femto/lib/libO/Linux.intel.64.openmpi.default/esmf.mk
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/sciclone/home10/yinglong/esmf_femto/lib/libO/Linux.intel.64.openmpi.default/

    on WW: 
      setenv ESMFMKFILE /sciclone/home10/yinglong/esmf_WW/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/sciclone/home10/yinglong/esmf_WW/lib/libO/Linux.intel.64.mvapich2.default

    on James:
      setenv ESMFMKFILE /ches/home00/yinglong/esmf_James/lib/libO/Linux.intel.64.openmpi.default/esmf.mk
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/ches/home00/yinglong/esmf_James/lib/libO/Linux.intel.64.openmpi.default

    on Cyclops:
      setenv ESMFMKFILE /sciclone/home10/yinglong/esmf_cyclops/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/sciclone/home10/yinglong/esmf_cyclops/lib/libO/Linux.intel.64.mvapich2.default

2. SCHISM: You need to have SCHISM built with `cmake` and an environment
   variable `SCHISM_BUILD_DIR` defined that points to your SCHISM build
   directory containing `lib/libhydro.a`, e.g.

    on femto:
      setenv SCHISM_BUILD_DIR /sciclone/home10/yinglong/git/schism/build.femto

    on WW:
      setenv SCHISM_BUILD_DIR /sciclone/home10/yinglong/git/schism/build.whirlwind

    on James:
      setenv SCHISM_BUILD_DIR /ches/home00/yinglong/git/schism/build_james

    on Cyclops:
      setenv SCHISM_BUILD_DIR /sciclone/home10/yinglong/git/schism/build.cyclops

3. Optionally, you need PDAF:

   on WW:
      setenv PDAF_LIB_DIR /sciclone/home10/yinglong/PDAF-D_V1.14/lib/Whirlwind

   on Cyclops:
      setenv PDAF_LIB_DIR /sciclone/home10/yinglong/PDAF-D_V1.14/lib/Cyclops

## Compilation (in schism-esmf)

    make distclean    [make version older than v3.82 will spit out an error on undef cmd]
    make all 

This will produce (1) the `libesmf_schism.a` `libschism_nuopc.a` library, which can be linked in
other ESMF applications that want to include SCHISM as a component.  This is
the primary application.

#Also (2) `schism_esmf_test` and (3) `concurrent_esmf_test` executables are produced
#that can be executed from any directory that contains a schism setup.  Results
#should exactly match those produced from a `pschism` uncoupled executable.

## Bugs and contributing

Please report bugs to <carsten.lemmen@hereon.de>.  You are very welcome to contribute
to improving this cap. Fork on github, make changes, and file a pull request.
