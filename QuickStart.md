# SCHISM-ESMF Quick Start Guide

**Complete walkthrough**: Get SCHISM-ESMF building and running in under 30 minutes.

> **Quick Reference**: See [ReadMe.md](./ReadMe.md) for build summary. This guide provides **detailed explanations** and **platform-specific tips**.

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
# Create  and activate environment
mamba create -n esmf "esmf=8.9.1=mpi_mpich*"
mamba activate esmf

# Verify MPI support (should NOT show "mpiuni")
grep ESMF_COMM $CONDA_PREFIX/lib/esmf.mk
```

You may also use the `mpi_openmpi*` variant if preferred.

### Using Homebrew

```bash
brew install gcc open-mpi esmf

# Point to the installed esmf.mk
export ESMFMKFILE=$(brew --prefix esmf)/lib/esmf.mk

# Verify MPI support (should NOT show "mpiuni")
grep ESMF_COMM $ESMFMKFILE
```

If you prefer MPICH, swap `open-mpi` for `mpich` in the install step.

### Using Spack

```bash
spack install esmf+mpi
spack load esmf
```


### From Source

You need to manually install in your operating system the dependencies of ESMF, typically `netcdf`, `xerces` and `mpi`; how these are installed, depends heavily on your package manager.

```bash
export ESMF_DIR=/devel/esmf/esmf # or any other
export ESMF_COMM=mpich           # or openmpi
cd $ESMF_DIR
make
make install
```

Then set an environment variable to point to the `esmf.mk` file:

```bash
export ESMFMKFILE=$ESMF_DIR/lib/esmf.mk
```

## Step 2: Build SCHISM

From the installation before, set the environment variables for the C and Fortran compiler, for a `mamba` system they would be:

```bash
export FC=$CONDA_PREFIX/bin/gfortran CC=$CONDA_PREFIX/bin/clang
```

Choose a directory where your sources reside and reference them with environment variables; you may choose different ones than those suggested here.

```bash
export SCHISM_BASE=$HOME/devel/schism/schism
export SCHISM_BUILD_DIR=$SCHISM_BASE/build
export SCHISM_ESMF_BASE=$SCHISM_BASE/../schism-esmf
export SCHISM_ESMF_BUILD_DIR=$SCHISM_ESMF_BASE/build
````

Then, clone and build SCHISM with CMake (not the legacy Makefile).

```bash
mkdir -p $SCHISM_BASE/.. # make sure the parent directory exists
git clone --recurse-submodules --depth=1 https://github.com/schism-dev/schism.git $SCHISM_BASE

mkdir -p $SCHISM_BUILD_DIR # make sure the build directory exists
cd $SCHISM_BUILD_DIR

# Configure with CMake
cmake -S $SCHISM_BASE/src -B $SCHISM_BUILD_DIR \
  -DCMAKE_Fortran_COMPILER=$FC \
  -DCMAKE_C_COMPILER=$CC \
  -DUSE_PARMETIS=ON -DBLD_STANDALONE=ON \
  -DUSE_WWM=OFF

# Build (adjust -j for your CPU count)
cmake --build $SCHISM_BUILD_DIR --parallel 8 --target pschism
```

**Verify SCHISM libraries**:
```bash
ls -lh $SCHISM_BUILD_DIR/lib/
# Should see: libcore.a, libhydro.a, libparmetis.a, libmetis.a
```

## Step 3: Build SCHISM-ESMF

```bash
# Clone this repository
git clone --depth=1 --recurse-submodules https://github.com/schism-dev/schism-esmf.git $SCHISM_ESMF_BASE
cd $SCHISM_ESMF_BASE

# Build the cap
mkdir -p $SCHISM_ESMF_BUILD_DIR
cmake -S $SCHISM_ESMF_BASE -B $SCHISM_ESMF_BUILD_DIR -DCMAKE_Fortran_COMPILER=mpifort -DSCHISM_REQUIRE_turbulence=OFF
cmake --build $SCHISM_ESMF_BUILD_DIR -- -j8
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
