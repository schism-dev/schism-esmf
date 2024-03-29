# This code is part of the SCHISM-ESMF interface.
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

include ../include/Rules.mk

ifneq ($(origin MOSSCO_DIR), environment)
$(error MOSSCO_DIR has to be set in environment.)
endif

ifeq ($(wildcard $(MOSSCO_DIR)/src/components/netcdf_component.F90),)
$(error MOSSCO_DIR=$(MOSSCO_DIR) does not seem to be correct)
endif

ifdef ESMFMKFILE
ifeq ($(FORTRAN_COMPILER),XLF)
CPPFLAGS += -WF,-DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) -WF,-DESMF_VERSION_MINOR=$(ESMF_VERSION_MINOR)
else
CPPFLAGS += -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) -DESMF_VERSION_MINOR=$(ESMF_VERSION_MINOR)
endif
endif

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(F90FLAGS) $<

default: all

time_aggregation_component.o: dep-esmf time_aggregation_component.F90 mossco_strings.o \
	mossco_component.o mossco_field.o mossco_state.o mossco_attribute.o \
	mossco_config.o mossco_time.o
netcdf_component.o: dep-esmf netcdf_component.F90 mossco_netcdf.o git-sha.h \
	mossco_component.o mossco_config.o
mossco_netcdf.o: dep-esmf mossco_netcdf.F90 git-sha.h mossco_strings.o mossco_field.o \
	mossco_state.o mossco_time.o mossco_gridspec.o mossco_variable_types.o \
	mossco_grid.o mossco_mesh.o

mossco_strings.o: dep-esmf mossco_strings.F90 mossco_memory.o
mossco_memory.o: dep-esmf mossco_memory.F90
mossco_field.o: dep-esmf mossco_field.F90 mossco_attribute.o mossco_grid.o mossco_geom.o \
	mossco_loc.o
mossco_attribute.o: dep-esmf mossco_attribute.F90
mossco_geom.o: dep-esmf mossco_geom.F90 mossco_mesh.o mossco_locstream.o
mossco_state.o: dep-esmf mossco_state.F90

time_aggregation_component.F90 netcdf_component.F90: dep-esmf 
	ln -sf $(MOSSCO_DIR)/src/components/$@  $@

mossco_netcdf.F90 mossco_strings.F90 mossco_field.F90 mossco-state.F90 \
mossco_time.F90 mossco_gridspec.F90 mossco_variable_types.F90 \
mossco_grid.F90 mossco_geom.F90 mossco_mesh.F90 mossco_memory.F90 \
mossco_loc.F90 mossco_attribute.F90 mossco_locstream.F90 mossco_state.F90 \
mossco_component.F90 mossco_config.F90:
	ln -sf $(MOSSCO_DIR)/src/utilities/$@  $@

git-sha.h:
	ln -sf $(MOSSCO_DIR)/src/include/$@  $@

.PHONY:  clean default all libmossco esmf

all: esmf libmossco

libmossco: netcdf_component.o time_aggregation_component.o
	$(AR) crvs libmossco.a *.o

distclean: clean
	$(RM) *.F90

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
