# This Makefile is part of the SCHISM-ESMF interface
#
# @copyright (C) 2018, 2019, 2020 Helmholtz-Zentrum Geesthacht
# @author Carsten Lemmen carsten.lemmen@hzg.de
# @author Richard Hofmeister richard.hofmeister@hzg.de
#
# @license under the Apache License, Version 2.0 (the "License");
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

# include your ESMF Makefile (esmf.mk), and obtain compiler settings
# from this file
ifndef ESMFMKFILE
$(error ESMFMKFILE has to be set in environment)
endif
include $(ESMFMKFILE)

F90=$(ESMF_F90COMPILER)
LIBS=$(ESMF_F90ESMFLINKLIBS)
CPPFLAGS=$(ESMF_F90COMPILEOPTS)
F90FLAGS=$(ESMF_F90COMPILEPATHS)
LDFLAGS+=$(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS)

# add SCHISM libraries
ifndef SCHISM_BUILD_DIR
$(error SCHISM_BUILD_DIR has to be set in environment.)
endif

SCHISM_BUILD_DIR:= $(shell readlink --canonicalize ${SCHISM_BUILD_DIR})

ifeq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libhydro.a),)
$(error SCHISM has to be compiled before ESMF-SCHISM.)
endif

# @todo parmetis should have been included in lschism_esmf, but
# that does not seem to work cross-platform ...
LIBS+= -lschism_esmf -lparmetis -lmetis
F90FLAGS+= -I$(SCHISM_BUILD_DIR)/include -I src/model -I src/schism
LDFLAGS+= -L$(SCHISM_BUILD_DIR)/lib -L.

EXPAND_TARGETS= expand_schismlibs

ifneq ($(wildcard $(SCHISM_BUILD_DIR)/lib/libfabm.a),)
  $(info Include fabm libraries from $(SCHISM_BUILD_DIR)/lib/libfabm*.a)
  EXPAND_TARGETS+= expand_fabmlibs
  F90FLAGS += -DUSE_FABM
endif

.PHONY: all lib test schism_esmf_test concurrent_esmf_test schism_esmf_lib triple_schism
default: all

# User-callable make targets

all: lib test

lib: schism_esmf_lib

test: concurrent_esmf_test triple_schism

# Internal make targets
SCHISM_OBJS=$(addprefix src/schism/,schism_cmi_esmf.o schism_esmf_util.o schism_bmi.o)
MODEL_OBJS=$(addprefix src/model/,atmosphere_cmi_esmf.o)

concurrent_esmf_test: $(SCHISM_OBJS) $(MODEL_OBJS) concurrent_esmf_test.o
	$(F90) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

triple_schism: $(SCHISM_OBJS) $(MODEL_OBJS) triple_schism.o
	$(F90) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

schism_esmf_lib: $(SCHISM_OBJS) $(MODEL_OBJS) $(EXPAND_TARGETS)
	$(AR) crs libschism_esmf.a  $(SCHISM_OBJS) $(MODEL_OBJS) .objs/*/*.o

expand_schismlibs:
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
expand_fabmlibs:
	$(shell mkdir -p .objs/sf; cd .objs/sf; for L in $(SCHISM_BUILD_DIR)/lib/lib*fabm_schism.a ; do $(AR) x $$L; done)
	$(shell mkdir -p .objs/f; cd .objs/f; $(AR) x $(SCHISM_BUILD_DIR)/lib/libfabm.a )

$(SCHISM_OBJS):
	make -C src/schism esmf

$(MODEL_OBJS):
	make -C src/model esmf

%.o: %.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<

%.mod: %.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<

clean:
	$(RM) *.o *.mod
	$(RM) $(SCHISM_OBJS) $(MODEL_OBJS)

distclean: clean
	$(RM) -rf .objs
	$(RM) -f fort.* flux.dat param.out.nml total.dat total_TR.dat mirror.out
	$(RM) -f concurrent_esmf_test triple_schism libschism_esmf.a
	$(RM) -f outputs/*nc
	$(RM) -f outputs/nonfatal*nc
	$(RM) -f PET*

# $(shell cd .objs/sf; nm fabm_schism.F90.o|grep 'fabm_mp\|fabm_types_mp' | awk '{printf $$2 " "; gsub("fabm_mp","s_fabm_mp",$$2); gsub("fabm_types_mp","s_fabm_types_mp",$$2); print $$2}'>replace.tsv; objcopy --redefine-syms=replace.tsv fabm_schism.F90.o)
# $(shell cd .objs/f; for O in *.o ; do rm -f replace.tsv; nm -f posix $$O |grep 'fabm_mp\|fabm_types_mp' | awk '{printf $$1 " "; gsub("fabm_mp","s_fabm_mp",$$1); gsub("fabm_types_mp","s_fabm_types_mp",$$1); print $$1}'>replace.tsv; printf "fabm._ s_fabm._\nfabm_types._ s_fabm_types._\n">> replace.tsv; objcopy --redefine-syms=replace.tsv $$O; done)
