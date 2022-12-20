# This Makefile is part of the SCHISM-ESMF interface
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

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
