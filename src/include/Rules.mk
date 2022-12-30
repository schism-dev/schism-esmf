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

# We require at least ESMF 8.1
ESMF_GT_8_1 := $(shell [ $(ESMF_VERSION_MAJOR) -gt 8 -o \( $(ESMF_VERSION_MAJOR) -eq 8 -a $(ESMF_VERSION_MINOR) -ge 1 \) ] && echo true)

ifneq ($(ESMF_GT_8_1),true)
$(error Your ESMF version $(ESMF_VERSION_MAJOR).$(ESMF_VERSION_MINOR) is too old. At least ESMF 8.1 is required)
endif

F90=$(ESMF_F90COMPILER)
LIBS=$(ESMF_F90ESMFLINKLIBS)
CPPFLAGS=$(ESMF_F90COMPILEOPTS)
F90FLAGS=$(ESMF_F90COMPILEPATHS)
LDFLAGS+=$(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS)

ESMF_COMM = $(strip $(shell grep "^. ESMF_COMM:" $(ESMFMKFILE) | cut -d':' -f2-))
ESMF_COMPILER = $(strip $(shell grep "^. ESMF_COMPILER:" $(ESMFMKFILE) | cut -d':' -f2-))

ifeq ("x$(ESMF_COMM)","xmpiuni")
	USE_MPI ?= false
else
	USE_MPI ?= true
endif

# Determine the original compilers (fortran and c++) used for the combination of compiler and device

ifeq ($(ESMF_FC),)

# We only implemented for a subset of MPI implementations
ifeq (,$(filter $(ESMF_COMM),mpich mpich2 mpich3 mvapich2 openmpi intelmpi))
$(error The communicator $(ESMF_COMM) is not implemented yet, please file a bug report)
endif

# OpenMPI section
ifeq ($(ESMF_COMM),openmpi)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) --showme:command 2> /dev/null | rev | cut -d'/' -f1 | rev)
	ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) --showme:command 2> /dev/null | cut -d'-' -f1)
ifeq ($(ESMF_FC),)
ifeq ($(ESMF_F90COMPILER),mpifort)
	ESMF_FC:=$(shell mpif90 --showme:command 2> /dev/null | awk -F/ '{print $NF}')
endif
endif
endif

# IntelMPI section
ifeq ($(ESMF_COMM),intelmpi)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -show 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
	ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) -show 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
endif

# mpich mpich2, mpich3, mvapich2 sections
#ifeq ($(ESMF_COMM),mvapich2)
ifneq (,$(filter $(ESMF_COMM),mpich mpich2 mpich3 mvapich2))
  ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 )
  ESMF_FC:=$(shell basename $(ESMF_FC))
  ESMF_FC:=$(shell echo $(ESMF_FC) | cut -d'-' -f1)
  ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 )
  ESMF_CC:=$(shell basename $(ESMF_CC))
  ESMF_CC:=$(shell echo $(ESMF_CC) | cut -d'-' -f1)
endif

# Make a correction on quadruplets that end in the compiler name and start with x86_64
ifeq ($(ESMF_FC),x86_64)
  ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 )
  ESMF_FC:=$(shell basename $(ESMF_FC))
  ESMF_FC:=$(shell echo $(ESMF_FC) | cut -d'-' -f4)
  ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 )
  ESMF_CC:=$(shell basename $(ESMF_CC))
  ESMF_CC:=$(shell echo $(ESMF_CC) | cut -d'-' -f4)
endif

# Finally exit if none of the above produced a valid ESMF_FC or EMSF_CC, pointing
# to a possible inconsistency between ESMFMKFILE and PATH
ifeq ($(ESMF_FC),)
$(error $(ESMF_F90COMPILER) is *not* based on $(ESMF_COMM)!)
endif
ifeq ($(ESMF_FC),)
$(error $(ESMF_CXXCOMPILER) is *not* based on $(ESMF_COMM)!)
endif

export ESMF_FC
export ESMF_CC
export ESMF_COMPILER
export ESMF_COMM
$(info Fortran/C++ compilers are ${ESMF_FC} and ${ESMF_CC} based on ${ESMF_COMPILER}/${ESMF_COMM} device)
endif

ESMF_OPENMP = $(strip $(shell grep "^. ESMF_OPENMP:" $(ESMFMKFILE) | cut -d':' -f2-))
ifeq ("$(ESMF_OPENMP)","OFF")
	USE_OMP ?= false
else
	USE_OMP ?= true
endif

export USE_OMP
export USE_MPI

# git section 
CPPFLAGS+=-DSCHISM_ESMF_GIT_COMMIT_HASH="\"$(shell git rev-parse --short HEAD || unknown)\""

endif # ESMF_OPENMP, ensuring only 1 iteration of this part
endif # ESMFMKFILE

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
	$(error ESMF Makefile snippet does not exist as $(ESMFMKFILE))
endif 

dep-schism: 
ifndef SCHISM_BUILD_DIR
	$(error SCHISM_BUILD_DIR has to be set)
endif 
ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
	$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

dep-pdaf:
ifndef PDAF_LIB_DIR
	$(error PDAF_LIB_DIR has to be set)
endif 
ifeq ($(wildcard $(PDAF_LIB_DIR)/libpdaf-d.a),)
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
