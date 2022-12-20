# This Makefile is part of the SCHISM-ESMF interface
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

ifndef SCHISM_VERIFICATION_DIR
$(error SCHISM_VERIFICATION_DIR has to be set in environment)
endif

ifeq ($(wildcard $(SCHISM_VERIFICATION_DIR)/Test_HeatPool/param.nml),)
$(error Please point SCHISM_VERIFICATION_DIR to the directory containing Test_HeatPool)
endif

D=$(SCHISM_VERIFICATION_DIR)/Test_HeatPool

COPY_FILES=bctides.in drag.gr3 vgrid.in hgrid.ll

.PHONY: default all subdirs $(SUBDIRS)

default: all
all: subdirs setups
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

distclean: clean
	@$(RM) *.in *.nml *.ic *.prop *.gr3

setups: $(COPY_FILES) global.nml param.nml hgrid.gr3

$(COPY_FILES):
	cp $(D)/$@ .

clean:
	@$(RM) *.out PET*
	@$(RM) outputs/*nc outputs/nonfatal* outputs/fatal* outputs/local_to* outputs/global_to*

global.nml: param.nml
	python ../../scripts/param2global.py

hgrid.gr3: hgrid.ll
	ln -sf $< $@

param.nml:
	sed 's/ics = 1/ics = 2/g' $(D)/$@ > $@
