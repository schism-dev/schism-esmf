# This code is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2021 Helmholtz-Zentrum Hereon
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
# 
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

ifneq ($(origin SCHISM_BUILD_DIR), environment)
$(error SCHISM_BUILD_DIR has to be set in environment.)
endif

SCHISM_BUILD_DIR:= $(shell readlink --canonicalize ${SCHISM_BUILD_DIR})

ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

ifdef PDAF_BUILD_DIR
PDAF_BUILD_DIR:= $(shell readlink --canonicalize ${PDAF_BUILD_DIR})
ifeq ($(wildcard $(PDAF_BUILD_DIR)/lib/libpdaf-d.a),)
$(warning PDAF has to be compiled before ESMF-SCHISM, continuing without PDAF)
else
USE_PDAF=ON
endif
endif


# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lhydro -lcore -lparmetis -lmetis
F90FLAGS+= -I $(SCHISM_BUILD_DIR)/include
ifdef USE_PDAF
F90FLAGS +=  -I ../PDAF_bindings
endif
LDFLAGS+= -L $(SCHISM_BUILD_DIR)/lib

.SUFFIXES: .F90

%.o : %.F90
ifdef USE_PDAF
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) -DUSE_PDAF  $<
else
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) $<
endif

schism_cmi_nuopc.o: schism_esmf_util.o schism_bmi.o
schism_cmi_esmf.o: schism_esmf_util.o schism_bmi.o
schism_esmf_util.o: schism_bmi.o

.PHONY:  clean default all

default: all
all: nuopc esmf
esmf: schism_cmi_esmf.o
nuopc: schism_cmi_nuopc.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
