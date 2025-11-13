# API Documentation Setup Summary

This document summarizes the FORD-based API documentation infrastructure now available in the SCHISM-ESMF project.

## What Was Added

### 1. FORD Configuration (`ford.yml`)

Main configuration file for generating Fortran API documentation:

- **Source directory**: `./src` (scans all `.F90` files)
- **Excluded**: `PDAF_bindings` subdirectory
- **Output**: `./build/ford_docs/`
- **Features enabled**:
  - Call graphs
  - Source code browser
  - Search functionality
  - External links to ESMF and NetCDF documentation

### 2. CMake Integration

Updated `doc/CMakeLists.txt` with new targets:

```bash
# Build user guide only (MkDocs)
make docs

# Build API reference only (FORD)
make docs-api

# Build both
make docs-all
```

FORD detection is automatic - CMake checks for the `ford` module and warns if not found but doesn't fail the build.

### 3. Documentation Requirements

Updated `doc/requirements.txt`:

```
mkdocs>=1.5.0
mkdocs-material>=9.4.0
ford>=7.0.0
markdown>=3.4.0
```

All dependencies are now installable via pip or Read the Docs.

### 4. Custom Styling

Added `doc/ford-custom.css` with Material Design colors matching the MkDocs theme for visual consistency.

### 5. Documentation Guides

Created comprehensive guides:

- **`doc/api-documentation-guide.md`**: Complete FORD syntax reference with examples
- **`doc/ford-example.md`**: Before/after documentation example using real code
- Updated **`doc/BuildingDocs.md`**: Added FORD build instructions

## Quick Start

### Install FORD

```bash
pip install ford
```

### Generate API Documentation

#### Direct Method

```bash
ford ford.yml
open build/ford_docs/index.html
```

#### Via CMake

```bash
cd build
cmake .. -DDOCS_ONLY=ON
cmake --build . --target docs-api
open ford_docs/index.html
```

## Documentation Workflow

### 1. Write FORD Comments

Add documentation to Fortran source files:

```fortran
!> @brief Brief one-line description
!!
!! Detailed multi-line description with full context.
!!
!! @param[in] input_var Description of input
!! @param[out] output_var Description of output
!! @author Your Name
subroutine my_subroutine(input_var, output_var)
  integer, intent(in) :: input_var   !< Inline description
  integer, intent(out) :: output_var !< Inline description
  
  ! Implementation...
end subroutine my_subroutine
```

### 2. Test Locally

```bash
ford ford.yml 2>&1 | tee ford_build.log
```

Check the log for warnings about:
- Undocumented procedures
- Missing parameters
- Broken cross-references

### 3. Verify Output

Open `build/ford_docs/index.html` and check:
- Module pages render correctly
- Call graphs display
- Search works
- Cross-links function

### 4. Commit Changes

```bash
git add src/**/*.F90
git commit -m "Add FORD documentation to module X"
```

## Current Documentation Coverage

Based on the source tree, priority modules to document:

### High Priority (Core Interfaces)

- [ ] `src/schism/schism_bmi.F90` - Basic Model Interface
- [ ] `src/schism/schism_esmf_cap.F90` - Main ESMF cap
- [ ] `src/schism/schism_nuopc_cap.F90` - NUOPC cap
- [ ] `src/driver/driver_schism.F90` - Main driver

### Medium Priority (Utilities)

- [ ] `src/schism/schism_esmf_util.F90` - ESMF utilities
- [ ] `src/schism/toplevel_schism_netcdf.F90` - NetCDF I/O
- [ ] `src/model/atmosphere_cmi_esmf.F90` - Atmosphere coupling

### Lower Priority (Optional Features)

- [ ] `src/PDAF_bindings/*.F90` - Data assimilation (already excluded from FORD)
- [ ] `src/mossco/*.F90` - MOSSCO integration

## Integration with Read the Docs

The `.readthedocs.yml` configuration is already set up to build both MkDocs and FORD documentation:

```yaml
build:
  os: ubuntu-22.04
  tools:
    python: "3.11"

mkdocs:
  configuration: mkdocs.yml

python:
  install:
    - requirements: doc/requirements.txt
```

When you push to the repository, Read the Docs will:
1. Install dependencies (including FORD)
2. Build MkDocs user guide
3. Build FORD API reference
4. Publish both to `https://schism-esmf.readthedocs.io/`

## Next Steps

### Immediate (Do Now)

1. **Test FORD locally**:
   ```bash
   pip install ford
   ford ford.yml
   ```

2. **Review generated output**:
   - Check build/ford_docs/index.html
   - Verify graphs render correctly
   - Test search functionality

3. **Fix any build warnings**:
   - Look for syntax errors in ford.yml
   - Check for missing module dependencies

### Short-term (This Week)

1. **Document key modules**:
   - Start with `schism_bmi.F90` (use ford-example.md as template)
   - Add docs to `schism_esmf_cap.F90`
   - Document main driver routines

2. **Update CI/CD**:
   - Add FORD build to GitHub Actions
   - Update GitLab CI to publish FORD docs
   - Test Read the Docs deployment

3. **Link documentation**:
   - Add "API Reference" section to index.md
   - Link from user guide to API docs
   - Cross-link between modules

### Long-term (This Month)

1. **Complete documentation coverage**:
   - Document all public interfaces
   - Add usage examples
   - Document data structures

2. **Generate documentation metrics**:
   - Track documentation coverage
   - Identify undocumented areas
   - Set coverage targets

3. **Establish documentation standards**:
   - Create templates for common patterns
   - Define required sections for modules
   - Set up automated checks

## Resources

### Documentation

- [API Documentation Guide](api-documentation-guide.md) - Complete FORD syntax reference
- [FORD Example](ford-example.md) - Before/after example
- [Building Documentation](BuildingDocs.md) - Build instructions

### External Links

- [FORD GitHub](https://github.com/Fortran-FOSS-Programmers/ford)
- [FORD Documentation](https://forddocs.readthedocs.io/)
- [FORD Wiki](https://github.com/Fortran-FOSS-Programmers/ford/wiki)

## Troubleshooting

### "ford: command not found"

```bash
pip install --user ford
# or
pip install ford
```

### FORD warnings about missing modules

The `ford.yml` is configured with preprocessor definitions:
```yaml
preprocess: true
macro: ESMF_CONTEXT
macro: _SCHISM_LOG_AND_FINALIZE_ON_ERROR_
```

If you see additional preprocessor macros in error messages, add them to `ford.yml`.

### Graphs not rendering

Ensure GraphViz is installed:
```bash
# macOS
brew install graphviz

# Linux
sudo apt-get install graphviz
```

### Output directory errors

FORD creates the output directory automatically. If you see permission errors:
```bash
mkdir -p build/ford_docs
chmod 755 build/ford_docs
```

## Questions?

- **User documentation**: See [Building Documentation](BuildingDocs.md)
- **Build issues**: See [Build Troubleshooting](build-troubleshooting.md)
- **Report problems**: Open an issue on GitHub
