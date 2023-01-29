# This code is part of the SCHISM-ESMF interface.
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

include include/Rules.mk

EXTRA_DIST=Readme.md

SUBDIRS=model schism driver

MOSSCO_LIB=
ifeq ($(origin MOSSCO_DIR), environment)
 SUBDIRS+=mossco
 MOSSCO_LIB+=-L mossco -lmossco
endif

ifdef USE_PDAF
  SUBDIRS+=PDAF_bindings
endif

MODELS=atmosphere_cmi_nuopc.o

SCHISM_NUOPC=schism_nuopc_cap.o schism_esmf_util.o schism_nuopc_util.o schism_bmi.o
SCHISM_ESMF=schism_nuopc_cap.o schism_esmf_util.o schism_bmi.o
MODEL_OBJS:=$(addprefix ../model/,$(MODELS))
SCHISM_NUOPC_OBJS:=$(addprefix ../schism/,$(SCHISM_NUOPC))
SCHISM_ESMF_OBJS:=$(addprefix ../schism/,$(SCHISM_ESMF))

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
NO_PARMETIS := $(shell echo ${SCHISM_NO_PARMETIS} | tr '[:lower:]' '[:upper:]')
ifeq ($(NO_PARMETIS),OFF)
  METIS_LDFLAGS =
  ifneq ($(PARMETISHOME),)
    METIS_LDFLAGS := -L$(PARMETISHOME)/lib
  endif
endif

LIBS+= -lhydro -lcore $(METIS_LDFLAGS) -lparmetis -lmetis
ifdef SCHISM_BUILD_DIR
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include -I./model -I./schism -I./driver -I./mossco
LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib
endif 

.SUFFIXES: .F90

default: all

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) -o $@ 	$<

main_esmf: main_esmf.o driver/toplevel_schism_netcdf.o model/atmosphere_cmi_nuopc.o \
	$(addprefix schism/,schism_esmf_cap.o schism_esmf_util.o schism_bmi.o)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -I. -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(F90FLAGS) $(MOSSCO_LIB) $(LDFLAGS) $(LIBS)

main_nuopc: main_nuopc.o driver/driver_schism_atm.o model/atmosphere_cmi_nuopc.o \
	$(addprefix schism/,schism_nuopc_cap.o schism_esmf_util.o schism_nuopc_util.o schism_bmi.o)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -I. -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(F90FLAGS) $(LDFLAGS) $(LIBS)

main_esmf.o: dep-esmf dep-schism driver/toplevel_schism_atm.o main_esmf.F90 Makefile
main_nuopc.o: dep-esmf dep-schism driver/driver_schism_atm.o driver/driver_schism.o main_nuopc.F90 Makefile

driver/driver_schism.o: driver/driver_schism.F90 driver/Makefile Makefile
	make -C driver driver_schism.o

driver/driver_schism_atm.o: driver/driver_schism_atm.F90 driver/Makefile Makefile
	make -C driver driver_schism_atm.o

driver/toplevel_schism_atm.o: driver/toplevel_schism_atm.F90 driver/Makefile Makefile
	make -C driver toplevel_schism_atm.o

.PHONY:  clean default all nuopc esmf


all: nuopc esmf
nuopc: main_nuopc
esmf: main_esmf
clean: 
	@for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done
	$(RM) PET*.ESMF_LogFile *.nc *.stdout *.stderr
	$(RM) main_nuopc main_esmf *.o *.mod */*.mod */*.o
