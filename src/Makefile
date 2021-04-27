# This code is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2021 Helmholtz-Zentrum Hereon
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
# @author Carsten Lemmen carsten.lemmen@hereon.de
#
# @license Apache License, Version 2.0 (the "License")
# You may not use this file except in compliance with the License.
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

ifneq ($(origin SCHISM_BUILD_DIR), environment)
$(error SCHISM_BUILD_DIR has to be set in environment.)
endif

SCHISM_BUILD_DIR:= $(shell readlink --canonicalize ${SCHISM_BUILD_DIR})

ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

EXTRA_DIST=Readme.md

SUBDIRS=model schism driver

MOSSCO_LIB=
ifeq ($(origin MOSSCO_DIR), environment)
 SUBDIRS+=mossco
 MOSSCO_LIB+=-L mossco -lmossco
endif

MODELS=atmosphere_cmi_nuopc.o

SCHISM_NUOPC=schism_cmi_nuopc.o schism_esmf_util.o schism_bmi.o
SCHISM_ESMF=schism_cmi_nuopc.o schism_esmf_util.o schism_bmi.o
MODEL_OBJS:=$(addprefix ../model/,$(MODELS))
SCHISM_NUOPC_OBJS:=$(addprefix ../schism/,$(SCHISM_NUOPC))
SCHISM_ESMF_OBJS:=$(addprefix ../schism/,$(SCHISM_ESMF))

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lhydro -lcore -lparmetis -lmetis
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include -I./model -I./schism -I./driver -I./mossco
LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib

.SUFFIXES: .F90

default: all

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) -o $@ 	$<

main_esmf: main_esmf.o driver/toplevel_schism_netcdf.o model/atmosphere_cmi_nuopc.o \
	$(addprefix schism/,schism_cmi_esmf.o schism_esmf_util.o schism_bmi.o)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -I. -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(F90FLAGS) $(MOSSCO_LIB) $(LDFLAGS) $(LIBS)

main_nuopc: main_nuopc.o driver/driver_schism_atm.o model/atmosphere_cmi_nuopc.o \
	$(addprefix schism/,schism_cmi_nuopc.o schism_esmf_util.o schism_bmi.o)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -I. -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(F90FLAGS) $(LDFLAGS) $(LIBS)

main_esmf.o: driver/toplevel_schism_atm.o main_esmf.F90 Makefile
main_nuopc.o: driver/driver_schism_atm.o driver/driver_schism.o main_nuopc.F90 Makefile

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
	$(RM) PET*.ESMF_LogFile *.nc *.stdout *.stderr
	$(RM) main_nuopc main_esmf *.o *.mod */*.mod */*.o
