# schism-esmf
Earth System Modeling Framework (ESMF) and National Unified Operational Prediction Capability (NUOPC) caps for SCHISM

-------------------------------------------------------

# Using it

## Prerequisites

1. ESMF: You need to have ESMF installed and an environment variable
   `ESMFMKFILE` defined (via `setenv` or `export`) that points to your ESMF
   installation, for example

    on femto: 
      setenv ESMFMKFILE /sciclone/home10/yinglong/esmf_femto/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk
      setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/sciclone/home10/yinglong/esmf_femto/lib/libO/Linux.intel.64.intelmpi.default

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
      setenv PDAF_BUILD_DIR /sciclone/home10/yinglong/PDAF-D_V1.14/lib/Whirlwind

   on Cyclops:
      setenv PDAF_BUILD_DIR /sciclone/home10/yinglong/PDAF-D_V1.14/lib/Cyclops

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

Please report bugs to <carsten.lemmen@hereon.de>.  You are very welcome to contribute
to improving this cap. Fork on github, make changes, and file a pull request.

---------------------------------------------------------------------------------
Joseph's notes:
1) param.nml: set large nspool, ihfskip to bypass outputs. rnday, dt must match global.nml 
