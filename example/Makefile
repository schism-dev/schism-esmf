# This Makefile is part of the SCHISM-ESMF interface
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

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
