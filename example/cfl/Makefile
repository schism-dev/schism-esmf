# This Makefile is part of the SCHISM-ESMF interface
#
# SPDX-FileCopyrightText: 2021-2022 Helmholtz-Zentrum Hereon
# SPDX-FileCopyrightText: 2020-2021 Helmholtz-Zentrum Geesthacht
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

export SCRIPTS_DIR=$(subst /example/cfl,/scripts,$(CURDIR))

SCRIPTS=${SCRIPTS_DIR}/run_timesteps.py ${SCRIPTS_DIR}/create_grid.py ${SCRIPTS_DIR}/make_setups.bash ${SCRIPTS_DIR}/run_setups.py

default: all setups

all: sample_param

clean:
	@$(RM) *.in param.nml
	@$(RM) -r dt_*

setups: timestep_setups resolution_setups

timestep_setups: $(SCRIPTS)
	python $<

resolution_setups: $(SCRIPTS)
	python $<

sample_param: param.nml vgrid.in

param.nml:
ifndef SCHISM_DIR
	$(error SCHISM_DIR has to be set in environment)
endif
	ln -s $(SCHISM_DIR)/sample_inputs/param.nml $@

vgrid.in:
ifndef SCHISM_DIR
	$(error SCHISM_DIR has to be set in environment)
endif
	ln -s $(SCHISM_DIR)/sample_inputs/vgrid.in.SZ $@
