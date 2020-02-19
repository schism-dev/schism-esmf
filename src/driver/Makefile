# This code is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2020 Helmholtz-Zentrum Geesthacht
# @author Carsten Lemmen carsten.lemmen@hzg.de
#
# @license under the Apache License, Version 2.0 (the "License");
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

EXTRA_DIST=Readme.md
SUBDIRS=

MODELS=atmosphere.o
MODEL_OBJS:=$(addprefix ../model/,$(MODELS))
SCHISM_OBJS:=$(addprefix ../schism/,schism_cmi_nuopc.o)

F90FLAGS += -I../model -I ../schism

.SUFFIXES: .F90

default: all

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS)  $<

.PHONY:  clean default all nuopc esmf

$(MODEL_OBJS):
	make -C ../model nuopc

$(SCHISM_OBJS):
	make -C ../schism nuopc

driver_schism.o: $(SCHISM_OBJS)
driver_schism--atm.o: $(MODEL_OBJS) $(SCHISM_OBJS)
toplevel_schism--atm.o: $(MODEL_OBJS) $(SCHISM_OBJS)

all: subdirs nuopc #esmf
esmf: toplevel_schism--atm.o
nuopc: driver_schism--atm.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod

subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@
