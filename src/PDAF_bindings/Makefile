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

ifneq ($(origin SCHISM_BUILD_DIR), environment)
$(error SCHISM_BUILD_DIR has to be set in environment.)
endif

SCHISM_BUILD_DIR:= $(shell readlink --canonicalize ${SCHISM_BUILD_DIR})

ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

ifeq ($(wildcard ../compiler.mk),)
$(warning cannot find compiler Makefile snippet "../compiler.mk")
else
include ../compiler.mk
endif

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lhydro -lcore -lparmetis -lmetis
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include #-I $(PDAF_BUILD_DIR)/src

# Add compiler-specific flags to CFLAGS and F90FLAGS

ifneq (,$(findstring g95,$(ESMF_FC)))
  F90FLAGS += -O3 -fdefault-real-8 # double all reals
else
ifneq (,$(findstring gfortran, $(ESMF_FC)))
  #F90FLAGS += -freal-4-real-8 # make only real double precision
  F90FLAGS += -O3 -fdefault-real-8 # double all reals
else
ifneq (,$(findstring ifort,$(ESMF_FC)))
  F90FLAGS += -r8 -O3 -xHOST #Change to PDAF setting #-real-size 64 # make only real double precision
else
ifneq (,$(findstring nagf,$(ESMF_FC)))
  F90FLAGS += -double # modern version of -r8
  CFLAGS += -fno-strict-overflow
else
ifneq (,$(findstring pgf,$(ESMF_FC)))
  F90FLAGS += -O3 -fast -Mr8 -Mr8intrinsics -Mbyteswapio -Mpreprocess # modern version of -r8
else
$(error "Promotion of real to 8 byte not implemented for $(ESMF_FC)")
endif # pg
endif # nag
endif # intel
endif # gfortran
endif # g95

LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L$(PDAF_BUILD_DIR) -lpdaf-d

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) $<

#schism_cmi_nuopc.o: schism_esmf_util.o schism_bmi.o
#schism_cmi_esmf.o: schism_esmf_util.o schism_bmi.o
#schism_esmf_util.o: schism_bmi.o

.PHONY:  clean default all

default: all
all: nuopc esmf
esmf: parser_mpi.o mod_parallel_pdaf.o mod_assimilation.o init_parallel_pdaf.o init_pdaf.o init_pdaf_info.o finalize_pdaf.o \
      init_ens_pdaf.o next_observation_pdaf.o distribute_state_pdaf.o output_netcdf_pdaf.o prepoststep_ens.o prepoststep_pdaf.o prepoststep_seek.o init_enkf.o init_seek.o init_seik.o \
      collect_state_pdaf.o init_dim_obs_pdaf.o obs_op_pdaf.o init_obs_pdaf.o prodrinva_pdaf.o init_obsvar_pdaf.o assimilate_pdaf.o \
      init_dim_obs_f_pdaf.o init_dim_obs_l_pdaf.o obs_op_f_pdaf.o init_obs_f_pdaf.o init_obs_l_pdaf.o prodrinva_l_pdaf.o init_obsvar_l_pdaf.o \
      init_n_domains_pdaf.o init_dim_l_pdaf.o g2l_state_pdaf.o l2g_state_pdaf.o g2l_obs_pdaf.o 
nuopc: schism_cmi_nuopc.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
