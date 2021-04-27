# This code is part of the SCHISM-ESMF interface.
#
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

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -I ../schism $(F90FLAGS) $<

.PHONY:  clean default all

default: all
all: nuopc esmf
esmf: atmosphere_cmi_esmf.o
nuopc: atmosphere_cmi_nuopc.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
