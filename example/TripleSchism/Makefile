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

.PHONY: default all subdirs $(SUBDIRS)

SUBDIRS=sim_1 sim_2 sim_3

SCHISM_COUNT?=8

default: all
all: subdirs setups
subdirs: $(SUBDIRS)

$(SUBDIRS):
	mkdir -p $@/outputs
	make -C ../83 setups
	cp ../83/*.in ../83/*.nml ../83/*.ic ../83/*.prop ../83/*.gr3 $@

distclean: clean
	@$(RM) *.in *.nml *.ic *.prop *.gr3
	@$(RM) -r $(SUBDIRS)

clean:
	@$(RM) *.out PET*

setups: global.nml

global.nml: ../83/global.nml
	cp $< $@

../83/global.nml:
	make -C ../83 global.nml

multi_schism:
	@for i in $(shell seq 1 $(SCHISM_COUNT)); do \
	  S=`printf 'schism_%03d' $$i`; \
          mkdir -p $$S/outputs ;  \
	  cp ../83/*.in ../83/*.nml ../83/*.ic ../83/*.prop ../83/*.gr3 $$S ; \
	done
