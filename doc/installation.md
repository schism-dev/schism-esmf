---
title: Installation 
summary: Instructions to install the SCHISM/ESMF cap
SPDX-FileCopyrightText: 2022 Helmholtz-Zentrum Hereon
SPDX-License-Identifier: CC0-1.0
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de
---

# Installation
Earth System Modeling Framework (ESMF) and National Unified Operational Prediction Capability (NUOPC) caps for SCHISM

## Prerequisites

1. A development system including a Fortran and a C++ compiler, GNU make (>= 3.82) and CMake, and an installation of a message passing (MPI) implementation and the NetCDF library consistent with your Fortran/C++ toolchain.

2. Earth System Modeling Framework (ESMF): You need to have ESMF installed and an environment variable
   `ESMFMKFILE` defined (via `setenv` or `export`) that points to your ESMF
  installation.  For some system, you'd also need to amend the search path for libraries `LD_LIBRARY_PATH` to 
  include the ESMF library directory.   for example by issuing

  ```{csh}
    setenv ESMFMKFILE /my/path/to/esmf/lib/libO/Linux.intel.64.openmpi.default/esmf.mk
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/my/path/to/esmf/lib/libO/Linux.intel.64.openmpi.default/
  ```

3. SCHISM: You need to have SCHISM built with `cmake` and an environment
   variable `SCHISM_BUILD_DIR` defined that points to your SCHISM build
   directory containing `lib/libhydro.a`, e.g.

  ```
  setenv SCHISM_BUILD_DIR /my/path/to/schism/build
  ```

4. Optionally, you can have PDAF:

  ```
  setenv PDAF_LIB_DIR /my/path/to/pdaf/lib
  ````

## Compilation (in schism-esmf)

  ```bash
  make install-nuopc install-esmf 
  ```

This will produce (1) the `libesmf_schism.a` `libschism_nuopc.a` library, which can be linked in
other ESMF applications that want to include SCHISM as a component.  This is
the primary application.

## Bugs and wishes

Please report bugs and any other issues or wishes on GitHub's issue tracker https://github.com/schism-dev/schism-esmf/issues (preferred) or by email to <carsten.lemmen@hereon.de>.  

## Contributing

You are very welcome to contribute to improving this cap. Fork on github, make changes, and file a pull request. 

As the primary license of SCHISM-ESMF is Apache 2.0, you're bound by its section 5 stating that your contribution also falls under Apache 2.0.  