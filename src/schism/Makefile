# This code is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2021-2022 Helmholtz-Zentrum Hereon
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
# 
# @author Carsten Lemmen <carsten.lemmen@hereon.de>
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
