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

ESMF_COMM = $(strip $(shell grep "\# ESMF_COMM:" $(ESMFMKFILE) | cut -d':' -f2-))
ESMF_COMPILER = $(strip $(shell grep "\# ESMF_COMPILER:" $(ESMFMKFILE) | cut -d':' -f2-))

ifeq ("x$(ESMF_COMM)","xmpiuni")
	USE_MPI ?= false
else
	USE_MPI ?= true
endif

# OpenMPI section
ifeq ($(ESMF_COMM),openmpi)
	ESMF_FC ?= $(shell $(ESMF_F90COMPILER) --showme:command 2> /dev/null)
ifeq ($(ESMF_FC),)
ifeq ($(ESMF_F90COMPILER),mpifort)
	ESMF_FC:=$(shell mpif90 --showme:command 2> /dev/null)
endif
endif

ifeq ($(ESMF_FC),)
	$(error $(ESMF_F90COMPILER) is *not* based on $(ESMF_COMM)!)
endif

ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) --showme:command 2> /dev/null)
endif

# IntelMPI section
ifeq ($(ESMF_COMM),intelmpi)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -show 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
ifeq ($(ESMF_FC),)
	$(error $(ESMF_F90COMPILER) is *not* based on $(ESMF_COMM)!)
endif
endif

# MPIch 2 section
ifeq ($(ESMF_COMM),mpich2)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
ifeq ($(ESMF_FC),x86_64)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f4)
endif
ifeq ($(ESMF_FC),)
	$(error $(ESMF_F90COMPILER) is *not* based on $(ESMF_COMM)!)
endif
ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
endif

# MPIch 3 section
ifeq ($(ESMF_COMM),mpich3)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
ifeq ($(ESMF_FC),x86_64)
	ESMF_FC:=$(shell $(ESMF_F90COMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f4)
endif
ifeq ($(ESMF_FC),)
	$(error $(ESMF_F90COMPILER) is *not* based on $(ESMF_COMM)!)
endif
ESMF_CC:=$(shell $(ESMF_CXXCOMPILER) -compile_info 2> /dev/null | cut -d' ' -f1 | cut -d'-' -f1)
endif

ESMF_OPENMP = $(strip $(shell grep "\# ESMF_OPENMP:" $(ESMFMKFILE) | cut -d':' -f2-))
ifeq ("$(ESMF_OPENMP)","OFF")
	USE_OMP ?= false
else
	USE_OMP ?= true
endif

export USE_OMP
export USE_MPI
export ESMF_COMPILER
export ESMF_COMM
export ESMF_FC
export ESMF_CC

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
