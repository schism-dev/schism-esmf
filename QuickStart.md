# SCHISM-ESMF Quick Start Guide

**Complete walkthrough**: Get SCHISM-ESMF building and running in under 30 minutes.

> **Quick Reference**: See [Readme.md](./Readme.md) for build summary. This guide provides **detailed explanations** and **platform-specific tips**.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Step 1: Install ESMF with MPI](#step-1-install-esmf-with-mpi-support)
3. [Step 2: Build SCHISM](#step-2-build-schism)
4. [Step 3: Build SCHISM-ESMF](#step-3-build-schism-esmf)
5. [Step 4: Verify Build](#step-4-verify-the-build)
6. [Step 5: Run Example](#step-5-run-an-example-optional)
7. [Troubleshooting](#troubleshooting)
8. [Platform Notes](#platform-specific-notes)

---

## Prerequisites

- **Fortran compiler**: gfortran 11+ or Intel Fortran
- **C compiler**: gcc or clang
- **CMake**: 3.16 or later
- **MPI**: OpenMPI or MPICH (not MPIUNI stub)
- **Python**: 3.8+ (for setup scripts)

## Step 1: Install ESMF with MPI Support

### Using Conda/Mamba (Recommended)

```bash
# Create environment
mamba create -n esmf-env

# Activate environment
conda activate esmf-env

# Install ESMF with MPI (choose one)
mamba install -c conda-forge "esmf=8.9.0=mpi_mpich*"
# OR
mamba install -c conda-forge "esmf=8.9.0=mpi_openmpi*"

# Verify MPI support (should NOT show "mpiuni")
grep ESMF_COMM $CONDA_PREFIX/lib/esmf.mk
```

### Using Spack

```bash
spack install esmf+mpi
spack load esmf
```

### From Source

```bash
export ESMF_DIR=/path/to/esmf
export ESMF_COMM=mpich  # or openmpi
cd $ESMF_DIR
make
make install
```

## Step 2: Build SCHISM

SCHISM must be built with CMake (not the legacy Makefile).

```bash
# Clone SCHISM
git clone https://github.com/schism-dev/schism.git
cd schism

# Configure with CMake
mkdir -p build && cd build
cmake ../src \
  -DCMAKE_Fortran_COMPILER=$(which gfortran) \
  -DUSE_PARMETIS=ON \
  -DUSE_WWM=OFF

# Build (adjust -j for your CPU count)
make -j8

# Note the build directory path
export SCHISM_BUILD_DIR=$(pwd)
echo "SCHISM built in: $SCHISM_BUILD_DIR"
```

**Verify SCHISM libraries**:
```bash
ls -lh $SCHISM_BUILD_DIR/lib/
# Should see: libcore.a, libhydro.a, libturbulence.a, libyaml.a, libparmetis.a, libmetis.a
```

## Step 3: Build SCHISM-ESMF

```bash
# Clone this repository
git clone https://github.com/schism-dev/schism-esmf.git
cd schism-esmf

# Set environment variables
export ESMFMKFILE=/path/to/esmf.mk  # e.g., $CONDA_PREFIX/lib/esmf.mk
export SCHISM_BUILD_DIR=/path/to/schism/build

# Configure
mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=$(which gfortran)

# Build
cmake --build . -- -j8
```

**macOS specific**:
```bash
cmake .. \
  -DCMAKE_Fortran_COMPILER=$(which gfortran) \
  -DCMAKE_C_COMPILER=$(which gcc)
cmake --build . -- -j$(sysctl -n hw.ncpu)
```

**Expected output**:
```
[ 15%] Built target schism_interface_common
[ 31%] Built target schism_nuopc_interface
[ 42%] Built target schism_esmf_interface
[ 52%] Built target schism_model_libs
[ 84%] Built target schism_driver_libs
[ 92%] Linking Fortran executable main_esmf
[100%] Built target main_esmf
[100%] Built target main_nuopc
```

## Step 4: Verify the Build

```bash
# Check executables exist
ls -lh build/src/main_esmf build/src/main_nuopc

# Check what libraries are linked (macOS)
otool -L build/src/main_esmf | grep -E "esmf|mpi|schism"

# Check what libraries are linked (Linux)
ldd build/src/main_esmf | grep -E "esmf|mpi|schism"
```

## Step 5: Run an Example (Optional)

```bash
cd example/Test_QuarterAnnulus
make

# This will:
# 1. Create grid and boundary conditions
# 2. Set up ESMF configuration files
# 3. Run the coupled model
```

## Troubleshooting

### "ESMF_COMM is mpiuni"

Your ESMF was built without real MPI. Reinstall:
```bash
mamba remove esmf
mamba install -c conda-forge "esmf=*=mpi_mpich*"
```

### "Cannot find libhydro.a"

SCHISM wasn't built or `SCHISM_BUILD_DIR` is wrong:
```bash
# Rebuild SCHISM
cd /path/to/schism
mkdir -p build && cd build
cmake ../src
make
export SCHISM_BUILD_DIR=$(pwd)
```

### "Undefined symbols: mpi_*"

MPI library mismatch. Ensure ESMF and SCHISM use the same MPI:
```bash
# Check ESMF MPI
grep ESMF_COMM $ESMFMKFILE

# Rebuild SCHISM with matching MPI if needed
```

### "Module file not found"

Parallel build race condition. Try sequential build:
```bash
cmake --build . -- -j1
```

Or see `.github/BUILD_TROUBLESHOOTING.md` for the fix.

### Build still fails?

1. **Run validation script**:
   ```bash
   .github/test-cmake-instructions.sh
   ```

2. **Check detailed troubleshooting**:
   - `.github/BUILD_TROUBLESHOOTING.md` - All known issues and solutions
   - `.github/TEST_README.md` - Test framework documentation
   - `.github/copilot-instructions.md` - Build system details

3. **Enable verbose output**:
   ```bash
   cmake --build . -- VERBOSE=1 2>&1 | tee build.log
   ```

## Environment Summary

After successful build, your environment should have:

```bash
# ESMF
echo $ESMFMKFILE
grep ESMF_COMM $ESMFMKFILE  # Should show: mpich or openmpi

# SCHISM
echo $SCHISM_BUILD_DIR
ls $SCHISM_BUILD_DIR/lib/*.a  # Should list 6+ libraries

# Compilers
which gfortran
gfortran --version  # Should be 11+
```

## Next Steps

- **Read the documentation**: See `doc/` directory
- **Explore examples**: See `example/` directory
- **Understand the architecture**: See `.github/copilot-instructions.md`
- **Contribute**: Fork, branch, and submit PRs!

## Platform-Specific Notes

### macOS (Apple Silicon)

- Use Homebrew or conda gfortran (not Xcode's)
- MPI from conda-forge works well
- Expect ~5 minute build time on M1/M2

### Linux (x86_64)

- Most distro compilers work (gcc 11+)
- OpenMPI from package manager is fine
- Expect ~3 minute build time on modern CPU

### HPC Clusters

- Use module system for compilers/MPI
- May need to specify MPI wrappers explicitly
- Consider using Spack for dependency management

## Getting Help

- **Issues**: https://github.com/schism-dev/schism-esmf/issues
- **Discussions**: https://github.com/schism-dev/schism-esmf/discussions
- **SCHISM Forum**: https://ccrm.vims.edu/schismweb/

---

**Build time**: ~10 minutes (after dependencies installed)  
**Skill level**: Intermediate (Fortran/MPI experience helpful)  
**Last updated**: November 2025
