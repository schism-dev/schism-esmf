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

SUBDIRS=cfl 83 Test_QuarterAnnulus TripleSchism Test_HeatPool

.PHONY: default all subdirs $(SUBDIRS)

default: all
all: subdirs
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
