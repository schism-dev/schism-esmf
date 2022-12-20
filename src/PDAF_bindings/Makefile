# This code is part of the SCHISM-ESMF interface.
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

include ../include/Rules.mk

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lhydro -lcore -lparmetis -lmetis
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include -I$(PDAF_LIB_DIR)/../include

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

LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L$(PDAF_LIB_DIR) -lpdaf-d 

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(F90FLAGS) $<

#schism_nuopc_cap.o: schism_esmf_util.o schism_bmi.o
#schism_esmf_cap.o: schism_esmf_util.o schism_bmi.o
#schism_esmf_util.o: schism_bmi.o

.PHONY:  clean default all

default: all
all: nuopc esmf
esmf: dep-esmf dep-pdaf mod_parallel_pdaf.o mod_assimilation.o init_parallel_pdaf.o init_pdaf_info.o finalize_pdaf.o \
      init_ens_pdaf.o next_observation_pdaf.o distribute_state_pdaf.o output_netcdf_pdaf.o prepoststep_pdaf.o \
      collect_state_pdaf.o init_dim_obs_all.o PDAFomi_obs_op_schism.o assimilate_pdaf.o \
      init_n_domains_pdaf.o init_dim_l_pdaf.o g2l_state_pdaf.o l2g_state_pdaf.o \
      obs_A_pdafomi.o obs_Z_pdafomi.o obs_S_pdafomi.o obs_T_pdafomi.o obs_U_pdafomi.o obs_V_pdafomi.o callback_obs_pdafomi.o init_pdaf.o
nuopc: schism_nuopc_cap.o

clean:
	$(RM) PET*.ESMF_LogFile traceout *.nc *.stdout
	$(RM) *.o *.mod
