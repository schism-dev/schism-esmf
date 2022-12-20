# This Makefile is part of the SCHISM-ESMF interface
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

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
ifndef SCHISM_DIR
	$(error SCHISM_DIR has to be set in environment)
endif
ifeq ($(wildcard $(SCHISM_DIR)/sample_inputs/param.nml),)
	$(error Please point SCHISM_DIR to the source directory of SCHISM.)
endif
	ln -sf $(SCHISM_DIR)/sample_inputs/param.nml $@

vgrid.in:
ifndef SCHISM_DIR
	$(error SCHISM_DIR has to be set in environment)
endif
ifeq ($(wildcard $(SCHISM_DIR)/sample_inputs/param.nml),)
	$(error Please point SCHISM_DIR to the source directory of SCHISM.)
endif
	ln -sf $(SCHISM_DIR)/sample_inputs/vgrid.in.SZ $@

global.nml: param.nml
	python ../../scripts/param2global.py
