# schism-esmf
Earth System Modeling Framework cap for SCHISM

-------------------------------------------------------
To make:
need to first setenv (using whirlwind as an e.g. below):

##Set esmf.mk from ESMF to get envars etc; $SCHISM_BUILD_DIR is where your lib/libhydro.a is located

setenv ESMFMKFILE /sciclone/home10/yinglong/esmf/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk

setenv SCHISM_BUILD_DIR /sciclone/home10/yinglong/git/schism/build.whirlwind/

make distclean
make all

