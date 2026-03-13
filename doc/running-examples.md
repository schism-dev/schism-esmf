# Running SCHISM-ESMF Examples

This guide explains how to run the SCHISM-ESMF executables and set up example cases.

## Executables Overview

After building, two executables are created in `build/src/`:

### `main_esmf`
- **Purpose**: ESMF-based SCHISM driver
- **Use case**: Standalone SCHISM runs with ESMF infrastructure
- **Features**: ESMF time management, ESMF State containers
- **Configuration**: Uses `esmf.rc` for ESMF-specific settings

### `main_nuopc`
- **Purpose**: NUOPC-based SCHISM driver with coupling capabilities
- **Use case**: Coupled earth system models (ocean-atmosphere, etc.)
- **Features**: NUOPC interoperability layer, multi-model coupling
- **Configuration**: Uses NUOPC field dictionaries and run sequences

## Required Input Files

Both executables require a standard SCHISM setup with these files:

### Core SCHISM Files
- `param.nml` - Main SCHISM configuration
- `vgrid.in` - Vertical grid definition
- `hgrid.gr3` - Horizontal grid
- `bctides.in` - Boundary conditions (if applicable)
- `drag.gr3`, `manning.gr3` - Bottom friction parameters

### ESMF/NUOPC-Specific Files
- **For `main_esmf`**:
  - `esmf.rc` - ESMF configuration (optional)
  - `schism_component.cfg` - SCHISM component settings

- **For `main_nuopc`**:
  - `field_dictionary.yaml` - NUOPC field mappings
  - `schism_component.cfg` - Component configuration

## Running Examples

### Example 1: Quarter Annulus Test

The Test_QuarterAnnulus example demonstrates basic SCHISM-ESMF functionality.

```bash
cd example/Test_QuarterAnnulus

# Generate grid and setup
make setup

# Run with ESMF driver
../../build/src/main_esmf

# Or run with NUOPC driver
../../build/src/main_nuopc
```

**What it does**:
- Creates a quarter-annulus domain
- Sets up tidal forcing
- Runs a short simulation (typically hours to days)
- Outputs to `outputs/` directory

**Expected output files**:
- `outputs/schout_*.nc` - Model output (elevation, velocity, salinity, etc.)
- `outputs/mirror.out` - Runtime log
- `outputs/nonfatal_*` - Warning messages (if any)

### Example 2: Heat Pool Test

Tests temperature advection and diffusion.

```bash
cd example/Test_HeatPool

# Run the example
make

# Check results
ls outputs/
```

### Example 3: CFL Test

Tests Courant-Friedrichs-Lewy condition handling.

```bash
cd example/cfl

# Run the test
make

# Analyze timestep adaptation
grep "dt=" outputs/mirror.out
```

## Configuration Files Explained

### `param.nml`

Standard SCHISM namelist. Key parameters for ESMF runs:

```fortran
&CORE
  dt = 100.0          ! Time step (seconds)
  rnday = 1.0         ! Run duration (days)
  start_year = 2000
  start_month = 1
  start_day = 1
  start_hour = 0
/

&OPT
  ipre = 1            ! Pre-processing option
  ibc = 1             ! Boundary condition type
  ibtp = 1            ! Baroclinic/barotropic option
/

&SCHOUT
  nhot_write = 8640   ! Hotstart output frequency
  iout_sta = 1        ! Station output
  nspool_sta = 10     ! Station output frequency
/
```

### `schism_component.cfg`

ESMF/NUOPC component configuration:

```yaml
component_name: SCHISM
mesh_file: hgrid.gr3
start_time: 2000-01-01T00:00:00
stop_time: 2000-01-02T00:00:00
timestep: 100s
output_interval: 3600s
```

### `field_dictionary.yaml` (NUOPC only)

Defines coupling fields:

```yaml
sea_surface_height_above_geoid:
  standard_name: sea_surface_height_above_geoid
  canonical_units: m
  description: Sea surface elevation relative to geoid

sea_water_temperature:
  standard_name: sea_water_temperature
  canonical_units: K
  description: Water temperature

sea_water_velocity:
  standard_name: sea_water_velocity
  canonical_units: m s-1
  description: Water velocity components
```

## Running with MPI

Both executables support parallel execution:

```bash
# Using mpirun
mpirun -np 4 ../../build/src/main_esmf

# Using mpiexec
mpiexec -n 4 ../../build/src/main_nuopc

# SLURM on HPC
srun -n 48 ../../build/src/main_esmf
```

**MPI considerations**:
- Number of processors must be compatible with domain decomposition
- SCHISM uses ParMETIS for optimal load balancing
- See `param.nml` → `nproc_nd` and `nproc_ev` for processor distribution

## Creating Your Own Setup

### 1. Prepare Grid Files

```bash
# Create horizontal grid (use SMS, QGIS, or Python tools)
# → hgrid.gr3

# Create vertical grid
# → vgrid.in

# Generate boundary conditions
# → bctides.in
```

### 2. Create SCHISM Configuration

```bash
# Copy template
cp example/Test_QuarterAnnulus/param.nml .

# Edit for your domain
vi param.nml
```

### 3. Create ESMF Configuration

```bash
# Copy ESMF config
cp example/Test_QuarterAnnulus/schism_component.cfg .

# For NUOPC coupling
cp example/TripleSchism/field_dictionary.yaml .
```

### 4. Run

```bash
# Link executable
ln -s /path/to/build/src/main_esmf .

# Execute
./main_esmf
```

## Debugging Runs

### Check for errors

```bash
# View runtime log
tail -f outputs/mirror.out

# Check for fatal errors
grep -i "fatal\|error" outputs/mirror.out

# Check for warnings
ls outputs/nonfatal_*
```

### Common issues

**"Cannot open file"**:
- Check that all required input files exist
- Verify file paths in `param.nml`

**"MPI error"**:
- Ensure number of processors matches configuration
- Check MPI is properly initialized

**"Grid error"**:
- Validate grid files with SCHISM's `utility/Grid_Scripts/`
- Check for negative depths, overlapping elements

**"Timestep too large"**:
- Reduce `dt` in `param.nml`
- Check CFL condition: `dt < dx/c` where `c` is wave speed

## Output Analysis

### NetCDF Output

SCHISM-ESMF produces CF-compliant NetCDF files:

```python
import xarray as xr

# Load output
ds = xr.open_dataset('outputs/schout_0001.nc')

# View variables
print(ds.data_vars)

# Plot elevation
ds.elevation.isel(time=0).plot()
```

### Visualize with VisIT

```bash
# Launch VisIT
visit

# Open schout_*.nc files
# Add plots: Pseudocolor, Vector, etc.
```

### Use SCHISM's Post-processing Tools

```bash
cd $SCHISM_DIR/utility/Post-Processing-Fortran/

# Combine outputs
./combine_output11 -i ../../outputs

# Extract time series
./read_output10_allnodes.f90
```

## Performance Tips

1. **Optimal processor count**: Use powers of 2 or domain-compatible numbers
2. **I/O frequency**: Balance between data needs and performance (`nspool` settings)
3. **Timestep**: Largest stable `dt` improves performance
4. **Vertical layers**: More layers = higher accuracy but slower
5. **Output variables**: Only write needed variables (`iof_*` flags in `param.nml`)

## Advanced: Coupled Runs

For coupling SCHISM with other models (atmosphere, wave, etc.):

1. **Use `main_nuopc`** - provides NUOPC coupling layer
2. **Configure field dictionary** - map fields between components
3. **Set up run sequence** - define coupling timestepping
4. **Use ESMF_App** - or create custom driver

See `example/TripleSchism/` for multi-SCHISM coupling example.

## Troubleshooting Guide

| Issue | Solution |
|-------|----------|
| Segmentation fault | Check array bounds, increase stack size (`ulimit -s unlimited`) |
| NaN values | Reduce timestep, check boundary conditions |
| Slow performance | Optimize processor count, reduce output frequency |
| Missing output | Check `iof_*` flags in `param.nml`, verify disk space |
| Coupling errors | Validate field dictionary, check unit conversions |

## Getting Help

- **SCHISM documentation**: https://ccrm.vims.edu/schismweb/
- **ESMF user guide**: https://earthsystemmodeling.org/docs/
- **GitHub issues**: https://github.com/schism-dev/schism-esmf/issues
- **SCHISM forum**: https://groups.google.com/g/schism-model-users

## Next Steps

- Explore `scripts/` directory for setup automation
- Review `Plot/` directory for visualization scripts
- Check `doc/` for detailed documentation
- Experiment with different model configurations

---

**Last updated**: November 2025  
**Tested with**: SCHISM 5.10+, ESMF 8.7-8.9
