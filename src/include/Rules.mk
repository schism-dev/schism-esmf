# This Makefile snippet is part of the SCHISM-ESMF interface
#
# @copyright (C) 2022 Helmholtz-Zentrum Hereon
# @author Carsten Lemmen <carsten.lemmen@hereon.de>
# @license Apache License, Version 2.0 (the "License")
#
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

ifdef ESMFMKFILE
# prevent multiple inclusion by checking for ESMF_OPENMP
ifndef ESMF_OPENMP 
$(info Found ESMF Makefile fragment $(ESMFMKFILE))
include $(ESMFMKFILE)
F90=$(ESMF_F90COMPILER)
LIBS=$(ESMF_F90ESMFLINKLIBS)
CPPFLAGS=$(ESMF_F90COMPILEOPTS)
F90FLAGS=$(ESMF_F90COMPILEPATHS)
LDFLAGS+=$(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS)
# Find out whether we have OPENMP (ist this needed for PDAF?), then the relevant
# compiler flag is already set
ESMF_OPENMP := $(strip $(shell grep "\# ESMF_OPENMP:" $(ESMFMKFILE) | cut -d':' -f2-))
endif 
endif 

ifdef SCHISM_BUILD_DIR
$(info Found SCHISM build directory $(SCHISM_BUILD_DIR))
ifneq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(info Found compiled SCHISM libraries)
ifneq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libfabm.a),)
$(info Found compiled FABM library)
endif
endif
endif

.PHONY: dep-esmf dep-fabm dep-schism dep-pdaf default all

default: all

# Upstream requirements to be met
dep-esmf:
ifndef ESMFMKFILE
	$(error ESMF Makefile snippet variable not defined)
endif 
ifeq ($(wildcard $(ESMFMKFILE)),)
	$(error ESMF Makefile snippet not defined)
endif 

dep-schism: 
ifndef SCHISM_BUILD_DIR
	$(error SCHISM_BUILD_DIR has to be set)
endif 
ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
	$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

dep-pdaf:
ifndef PDAF_BUILD_DIR
	$(error PDAF_BUILD_DIR has to be set)
endif 
ifeq ($(wildcard $(PDAF_BUILD_DIR)/libpdaf-d.a),)
	$(error PDAF has to be compiled before ESMF-SCHISM.)
endif

dep-fabm: dep-schism 
ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libfabm.a),)
	$(error FABM has to be compiled before ESMF-SCHISM.)
endif

%.o: %.F90
	echo "==============="
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<

%.mod: %.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<
