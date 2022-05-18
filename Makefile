# This Makefile is part of the SCHISM-ESMF interface
#
# @copyright (C) 2021-2022 Helmholtz-Zentrum Hereon
# @copyright (C) 2018-2021 Helmholtz-Zentrum Geesthacht
# #
# @author Carsten Lemmen <carsten.lemmen@hereon.de>
# @author Richard Hofmeister
#
# @license Apache License, Version 2.0 (the "License");
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

include src/include/Rules.mk

DESTDIR?=./lib

ifdef PDAF_LIB_DIR
CPPFLAGS+= -DUSE_PDAF
USE_PDAF=ON
endif

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lschism_esmf -lparmetis -lmetis -lesmf
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include -I src/schism #-r8  ###-I src/model -I src/schism
##PDAF requires MKL (BLAS, LAPACK), this should already be provided by ESMF_FLAGS ...

ifdef USE_PDAF
LDFLAGS+= -L$(PDAF_LIB_DIR) -lpdaf-d -mkl -lpthread -lm -ldl
endif
ifeq ($(ESMF_COMPILER), intel)
LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L. -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm
else
ifeq ($(ESMF_COMPILER), gfortran)
# @todo still some lapack routines missing, so we need to link with either
# OpenBLAS or vecLibFort (osx), this should be configured automatically ... we
# really need to move to CMake
#LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L. -lpthread -lm -lvecLibFort
LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L. -lpthread -lm -lOpenBLAS
endif
endif

EXPAND_TARGETS= expand_schismlibs

ifneq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libfabm.a),)
  $(info Include fabm libraries from $(SCHISM_BUILD_DIR)/lib/libfabm*.a)
  EXPAND_TARGETS+= expand_fabmlibs
  F90FLAGS += -DUSE_FABM
endif

.SUFFIXES:
.PHONY: all lib test schism_nuopc_lib schism_esmf_lib install install-esmf install-nuopc pdaf
default: all

# User-callable make targets

all: lib test schism_nuopc_lib

lib: schism_esmf_lib schism_nuopc_lib

install: install-esmf install-nuopc

install-esmf:  schism_esmf_lib
	mkdir -p $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libhydro.a $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libcore.a $(DESTDIR)
	cp libschism_esmf.a $(DESTDIR)
	cp $(SCHISM_MODS) $(DESTDIR)

install-nuopc:  schism_nuopc_lib
	mkdir -p $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libhydro.a $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libmetis.a $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libparmetis.a $(DESTDIR)
	cp $(SCHISM_BUILD_DIR)/lib/libcore.a $(DESTDIR)
	cp libschism_cap.a $(DESTDIR)
	cp libschism_cap.a $(SCHISM_BUILD_DIR)/lib
	#cp $(SCHISM_NUOPC_MODS) $(DESTDIR)
	cp $(SCHISM_NUOPC_MODS) $(SCHISM_BUILD_DIR)/include/
	sed 's#@@SCHISM_BUILD_DIR@@#'$(SCHISM_BUILD_DIR)'#g' ./src/schism/schism_nuopc_cap.mk.in > $(DESTDIR)/schism.mk
	#sed 's#@@SCHISM_BUILD_DIR@@#'$(SCHISM_BUILD_DIR)'#g' ./src/schism/schism_nuopc_cap.mk.in > $(SCHISM_BUILD_DIR)/include/schism.mk

##test: concurrent_esmf_test triple_schism multi_schism schism_pdaf
test: pdaf 
pdaf: schism_pdaf

# Internal make targets for final linking
SCHISM_NUOPC_MODS=$(addprefix src/schism/,schism_nuopc_util.mod schism_nuopc_cap.mod)
SCHISM_NUOPC_OBJS=$(addprefix src/schism/,schism_nuopc_util.o schism_nuopc_cap.o)
SCHISM_ESMF_MODS=$(addprefix src/schism/,schism_esmf_cap.mod)
SCHISM_ESMF_OBJS=$(addprefix src/schism/,schism_esmf_cap.o)
SCHISM_MODS=$(addprefix src/schism/,schism_bmi.mod schism_esmf_util.mod)
SCHISM_OBJS=$(addprefix src/schism/,schism_bmi.o schism_esmf_util.o)
PDAF_OBJS=$(addprefix src/PDAF_bindings/,parser_mpi.o mod_parallel_pdaf.o mod_assimilation.o init_parallel_pdaf.o \
            init_pdaf.o init_pdaf_info.o finalize_pdaf.o init_ens_pdaf.o next_observation_pdaf.o \
            distribute_state_pdaf.o prepoststep_ens.o prepoststep_pdaf.o prepoststep_seek.o init_enkf.o init_seek.o init_seik.o \
            collect_state_pdaf.o init_dim_obs_pdaf.o obs_op_pdaf.o init_obs_pdaf.o prodrinva_pdaf.o init_obsvar_pdaf.o assimilate_pdaf.o \
            init_dim_obs_f_pdaf.o init_dim_obs_l_pdaf.o obs_op_f_pdaf.o init_obs_f_pdaf.o init_obs_l_pdaf.o prodrinva_l_pdaf.o init_obsvar_l_pdaf.o \
            init_n_domains_pdaf.o init_dim_l_pdaf.o g2l_state_pdaf.o l2g_state_pdaf.o g2l_obs_pdaf.o output_netcdf_pdaf.o)
#MODEL_OBJS=$(addprefix src/model/,atmosphere_cmi_esmf.o)

#concurrent_esmf_test: $(SCHISM_OBJS) $(MODEL_OBJS) concurrent_esmf_test.o
#	$(F90) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

ifdef USE_PDAF
schism_pdaf: install-esmf dep-pdaf $(PDAF_OBJS) $(SCHISM_OBJS) $(SCHISM_ESMF_OBJS) schism_pdaf.o
	$(F90) $(CPPFLAGS) $(PDAF_OBJS) $(SCHISM_OBJS) $(SCHISM_ESMF_OBJS) schism_pdaf.o -o $@ $(LDFLAGS)  -L./lib $(LIBS)
endif

schism_esmf_lib: dep-esmf dep-schism $(SCHISM_OBJS)  $(SCHISM_ESMF_OBJS) $(EXPAND_TARGETS)
	$(AR) crs libschism_esmf.a  $(SCHISM_OBJS) .objs/*/*.o

schism_nuopc_lib: dep-esmf dep-schism $(SCHISM_OBJS) $(SCHISM_NUOPC_OBJS) $(EXPAND_TARGETS)
	$(AR) crs libschism_cap.a  $(SCHISM_NUOPC_OBJS) $(SCHISM_OBJS) .objs/*/*.o

expand_schismlibs: dep-schism
	$(shell mkdir -p .objs/d; cd .objs/d; \
	$(AR) x $(SCHISM_BUILD_DIR)/lib/libcore.a ; \
		$(AR) x $(SCHISM_BUILD_DIR)/lib/libhydro.a ; \
		$(AR) x $(SCHISM_BUILD_DIR)/lib/libparmetis.a ; \
		$(AR) x $(SCHISM_BUILD_DIR)/lib/libmetis.a ; \
	)

# @todo the fabm lib symbols should be renamed, e.g., prefixed with schism_ to
# avoid duplicate symbols when coupling to other systems that also contain fabm
# A possible solution is provided by www.mossco.de/code in their
# scripts/rename_fabm_symbols.py
expand_fabmlibs: dep-fabm
	$(shell mkdir -p .objs/sf; cd .objs/sf; for L in $(SCHISM_BUILD_DIR)/lib/lib*fabm_schism.a ; do $(AR) x $$L; done)
	$(shell mkdir -p .objs/f; cd .objs/f; $(AR) x $(SCHISM_BUILD_DIR)/lib/libfabm.a )

$(PDAF_OBJS):
	make -C src/PDAF_bindings esmf

ifdef USE_PDAF
$(SCHISM_ESMF_OBJS): $(PDAF_OBJS) $(SCHISM_OBJS)
else
$(SCHISM_ESMF_OBJS): $(SCHISM_OBJS)
endif
	make -C src/schism esmf

ifdef USE_PDAF
$(SCHISM_NUOPC_OBJS): $(PDAF_OBJS) $(SCHISM_OBJS)
else
$(SCHISM_NUOPC_OBJS): $(SCHISM_OBJS)
endif
	make -C src/schism nuopc

ifdef USE_PDAF
$(SCHISM_OBJS): $(PDAF_OBJS)
else
$(SCHISM_OBJS):
endif
	make -C src/schism common

#$(MODEL_OBJS):
#	make -C src/model esmf

clean:
	$(MAKE) -C src clean
	$(RM) *.o *.mod
	$(RM) $(SCHISM_OBJS) $(MODEL_OBJS) $(PDAF_OBJS)

distclean: clean
	$(RM) -rf .objs
	$(RM) -f fort.* flux.dat param.out.nml total.dat total_TR.dat mirror.out
	$(RM) -f concurrent_esmf_test triple_schism multi_schism schism_pdaf libschism_esmf.a
	$(RM) -f outputs/*nc
	$(RM) -f outputs/nonfatal*nc
	$(RM) -f PET*

# $(shell cd .objs/sf; nm fabm_schism.F90.o|grep 'fabm_mp\|fabm_types_mp' | awk '{printf $$2 " "; gsub("fabm_mp","s_fabm_mp",$$2); gsub("fabm_types_mp","s_fabm_types_mp",$$2); print $$2}'>replace.tsv; objcopy --redefine-syms=replace.tsv fabm_schism.F90.o)
# $(shell cd .objs/f; for O in *.o ; do rm -f replace.tsv; nm -f posix $$O |grep 'fabm_mp\|fabm_types_mp' | awk '{printf $$1 " "; gsub("fabm_mp","s_fabm_mp",$$1); gsub("fabm_types_mp","s_fabm_types_mp",$$1); print $$1}'>replace.tsv; printf "fabm._ s_fabm._\nfabm_types._ s_fabm_types._\n">> replace.tsv; objcopy --redefine-syms=replace.tsv $$O; done)
