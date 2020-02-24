# schism-esmf
Earth System Modeling Framework cap for SCHISM

-------------------------------------------------------

# Using it

## Prerequisites

1. ESMF: You need to have ESMF installed and an environment variable
   `ESMFMKFILE` defined (via `setenv` or `export`) that points to your ESMF
   installation, for example

    setenv ESMFMKFILE /sciclone/home10/yinglong/esmf/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk

At the moment also need to:
     setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/sciclone/home10/yinglong/esmf/lib/libO/Linux.intel.64.mvapich2.default


2. SCHISM: You need to have SCHISM built with `cmake` and an environment
   variable `SCHISM_BUILD_DIR` defined that points to your SCHISM build
   directory containing `lib/libhydro.a`, e.g.

    setenv SCHISM_BUILD_DIR /sciclone/home10/yinglong/git/schism/build.whirlwind

## Compilation

    make distclean
    make all

This will produce (1) the `libesmf_schism.a` library, which can be linked in
other ESMF applications that want to include SCHISM as a component.  This is
the primary application.

Also (2) `schism_esmf_test` and (3) `concurrent_esmf_test` executables are produced
that can be executed from any directory that contains a schism setup.  Results
should exactly match those produced from a `pschism` uncoupled executable.

## Bugs and contributing

Please report bugs to <carsten.lemmen@hzg.de>.  You are very welcome to contribute
to improving this cap. Fork on github, make changes, and file a pull request.
