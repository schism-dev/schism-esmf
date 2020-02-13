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

ifneq ($(origin SCHISM_BUILD_DIR), environment)
$(error SCHISM_BUILD_DIR has to be set in environment.)
endif

SCHISM_BUILD_DIR:= $(shell readlink --canonicalize ${SCHISM_BUILD_DIR})

ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif


.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $<

main_nuopc: main_nuopc.o toplevel_component.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)

main_nuopc.o: toplevel_component.o
toplevel_component.o: atm.o ocn.o

.PHONY:  clean

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) main_nuopc *.o *.mod
