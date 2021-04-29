# This Makefile snippet is part of the SCHISM-ESMF interface.
#
# @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
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
#
# Find the communicator and determine whether this is parallel device, this
# is still buggy with mpiifort and needs improvement

include $(ESMFMKFILE)

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

compiler-info:
	@-echo USE_OMP=$(USE_OMP)
	@-echo USE_MPI=$(USE_MPI)
	@-echo ESMF_COMPILER=$(ESMF_COMPILER)
	@-echo ESMF_COMM=$(ESMF_COMM)
	@-echo ESMF_F90COMPILER=$(ESMF_F90COMPILER)
	@-echo ESMF_FC=$(ESMF_FC)
	@-echo ESMF_CC=$(ESMF_CC)
