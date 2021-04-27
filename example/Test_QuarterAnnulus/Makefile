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

ifndef SCHISM_VERIFICATION_DIR
$(error SCHISM_VERIFICATION_DIR has to be set in environment)
endif

ifeq ($(wildcard $(SCHISM_VERIFICATION_DIR)/Test_QuarterAnnulus/param.nml),)
$(error Please point SCHISM_VERIFICATION_DIR to the directory containing Test_QuarterAnnulus.)
endif

D=$(SCHISM_VERIFICATION_DIR)/Test_QuarterAnnulus

COPY_FILES=bctides.in drag.gr3 hgrid.gr3 param.nml vgrid.in

.PHONY: default all subdirs $(SUBDIRS)

default: all
all: subdirs setups
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

distclean: clean
	@$(RM) *.in *.nml *.ic *.prop *.gr3

setups: $(COPY_FILES) global.nml

$(COPY_FILES):
	cp $(D)/$@ .

clean:
	@$(RM) *.out PET*
	@$(RM) outputs/*nc outputs/nonfatal* outputs/fatal* outputs/local_to* outputs/global_to*

global.nml: param.nml
	python ../../scripts/param2global.py
