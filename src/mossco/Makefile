# This code is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2021 Helmholtz-Zentrum Hereon
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
# @author Carsten Lemmen carsten.lemmen@hereon.de
#
# @license Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# 		http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

ifneq ($(origin ESMFMKFILE), environment)
$(error ESMFMKFILE has to be set in environment.)
endif
include $(ESMFMKFILE)

ifneq ($(origin MOSSCO_DIR), environment)
$(error MOSSCO_DIR has to be set in environment.)
endif

MOSSCO_DIR:= $(shell readlink --canonicalize ${MOSSCO_DIR})

ifeq ($(wildcard $(MOSSCO_DIR)/src/components/netcdf_component.F90),)
$(error MOSSCO_DIR=$MOSSCO_DIR does not seem to be correct)
endif

ifeq ($(FORTRAN_COMPILER),XLF)
CPPFLAGS += -WF,-DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) -WF,-DESMF_VERSION_MINOR=$(ESMF_VERSION_MINOR)
else
CPPFLAGS += -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) -DESMF_VERSION_MINOR=$(ESMF_VERSION_MINOR)
endif

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(CPPFLAGS) $(F90FLAGS) $<

default: all

netcdf_component.o: netcdf_component.F90 mossco_netcdf.o git-sha.h \
	mossco_component.o mossco_config.o
mossco_netcdf.o: mossco_netcdf.F90 git-sha.h mossco_strings.o mossco_field.o \
	mossco_state.o mossco_time.o mossco_gridspec.o mossco_variable_types.o \
	mossco_grid.o mossco_mesh.o

mossco_strings.o: mossco_strings.F90 mossco_memory.o
mossco_memory.o: mossco_memory.F90
mossco_field.o: mossco_field.F90 mossco_attribute.o mossco_grid.o mossco_geom.o \
	mossco_loc.o
mossco_attribute.o: mossco_attribute.F90
mossco_geom.o: mossco_geom.F90 mossco_mesh.o mossco_locstream.o
mossco_state.o: mossco_state.F90

netcdf_component.F90:
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

libmossco: netcdf_component.o
	$(AR) crvs libmossco.a *.o

distclean: clean
	$(RM) *.F90

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
