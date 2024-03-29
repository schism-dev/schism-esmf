# This code is part of the SCHISM-ESMF interface.
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de
# SPDX-FileContributor: Richard Hofmeister

include ../include/Rules.mk

EXTRA_DIST=Readme.md
SUBDIRS=

CMI_MOSSCO=
CMI_NUOPC=atmosphere_cmi_nuopc.o
CMI_ESMF=atmosphere_cmi_esmf.o
ifeq ($(origin MOSSCO_DIR), environment)
CMI_MOSSCO+=netcdf_component.o
endif

CMI_NUOPC_OBJS:=$(addprefix ../model/,$(CMI_NUOPC))
CMI_ESMF_OBJS:=$(addprefix ../model/,$(CMI_ESMF))
CMI_MOSSCO_OBJS:=$(addprefix ../mossco/,$(CMI_MOSSCO))
CMI_NUOPC_SCHISM_OBJ:=$(addprefix ../schism/,schism_nuopc_cap.o)
CMI_ESMF_SCHISM_OBJ:=$(addprefix ../schism/,schism_esmf_cap.o)

F90FLAGS += -I../model -I ../schism
ifeq ($(origin MOSSCO_DIR), environment)
F90FLAGS += -I ../mossco
endif

.SUFFIXES: .F90

default: all

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS)  $<

.PHONY:  clean default all

$(CMI_NUOPC_SCHISM_OBJ):
	make -C ../schism nuopc

$(CMI_NUOPC_OBJS):
	make -C ../model nuopc

$(CMI_ESMF_SCHISM_OBJ):
	make -C ../schism esmf

$(CMI_ESMF_OBJS):
	make -C ../model esmf

$(CMI_MOSSCO_OBJS):
	make -C ../mossco all

driver_schism.o: dep-esmf dep-schism $(CMI_NUOPC_SCHISM_OBJ) driver_schism.F90
driver_schism_atm.o: dep-esmf dep-schism $(CMI_NUOPC_SCHISM_OBJ) $(CMI_NUOPC_OBJS) driver_schism_atm.F90
toplevel_schism_atm.o: dep-esmf dep-schism $(CMI_ESMF_SCHISM_OBJ) $(CMI_ESMF_OBJS) toplevel_schism_atm.F90
toplevel_schism_netcdf.o: dep-esmf dep-schism $(CMI_ESMF_SCHISM_OBJ) $(CMI_MOSSCO_OBJS) toplevel_schism_netcdf.F90

all: subdirs nuopc esmf
esmf: dep-esmf dep-schism toplevel_schism_atm.o toplevel_schism_netcdf.o
nuopc: dep-esmf dep-schism driver_schism_atm.o driver_schism.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod

subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@
