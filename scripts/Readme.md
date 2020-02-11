# Utility scripts for SCHISM-ESMF

This directory contains utility scripts to test the
SCHISM-ESMF component. 

> Do not execute the scripts in this directory, but 
  run `make` in `../example/` to see these in action.

run_timesteps.py
: Run a simulation ensemble with different timesteps
  resulting in different CFL relationships.

run_setups.py
: Run a simulation ensemble with different grid
  resolution.

make_setups.bash
: Script called by `run_*.py` to prepare a complete setup.

create_grid.py
: Creates a sample grid to be used with SCHISM. Is called
  by make_setup.bash
