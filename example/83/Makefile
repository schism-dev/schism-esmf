# This Makefile is part of the SCHISM-ESMF interface
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

ifndef SCHISM_DIR
$(error SCHISM_DIR has to be set in environment)
endif

ifeq ($(wildcard $(SCHISM_DIR)/sample_inputs/param.nml),)
$(error Please point SCHISM_DIR to the source directory of SCHISM.)
endif

ifndef MOSSCO_DIR
$(error MOSSCO_DIR has to be set in environment to generate the
example files here)
endif

D=${MOSSCO_DIR}/examples/esmf/schism

.PHONY: default all subdirs $(SUBDIRS)

default: all
all: subdirs setups
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

distclean: clean
	@$(RM) *.in *.nml *.ic *.prop *.gr3

clean:
	@$(RM) *.out PET* 
	@$(RM) outputs/*nc outputs/nonfatal* outputs/fatal* outputs/local_to* outputs/global_to*

bctides.in: copy-from-mossco

copy-from-mossco:
	ln -sf $D/bctides.in  $D/*.gr3 $D/*ic $D/*prop .

setups: copy-from-mossco sample_param

sample_param: param.nml vgrid.in

param.nml:
	ln -sf $(SCHISM_DIR)/sample_inputs/param.nml $@

vgrid.in:
	ln -sf $(SCHISM_DIR)/sample_inputs/vgrid.in.SZ $@

global.nml: param.nml
	python ../../scripts/param2global.py
