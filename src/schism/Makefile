# This code is part of the SCHISM-ESMF interface.
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

include ../include/Rules.mk

ifdef PDAF_LIB_DIR
USE_PDAF=ON
endif


# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lhydro -lcore -lparmetis -lmetis
F90FLAGS+= -I $(SCHISM_BUILD_DIR)/include
ifdef USE_PDAF
F90FLAGS +=  -I ../PDAF_bindings #-g -traceback -mcmodel=medium #-shared-intel
endif
LDFLAGS+= -L $(SCHISM_BUILD_DIR)/lib

.SUFFIXES: .F90

%.o : %.F90
ifdef USE_PDAF
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) -DUSE_PDAF  $<
else
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) $<
endif

.PHONY:  clean default all esmf nuopc common

default: all
all: nuopc esmf
esmf: dep-esmf dep-schism schism_esmf_cap.o
nuopc: dep-esmf dep-schism schism_nuopc_cap.o
common: dep-esmf dep-schism schism_bmi.o schism_esmf_util.o

schism_esmf_util.o:  dep-esmf dep-schism schism_bmi.o
schism_esmf_cap.o:   common
schism_nuopc_util.o: common
schism_nuopc_cap.o:  dep-esmf dep-schism schism_nuopc_util.o common

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
	$(RM) nohup.out
