#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: GitHub Copilot
#
# Test framework for validating CMake build instructions from copilot-instructions.md
# This script implements the CI-ready checklist and provides detailed diagnostics.

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test result counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# Function to print test results
print_test() {
    local status=$1
    local test_name=$2
    local message=${3:-}
    
    case $status in
        PASS)
            echo -e "${GREEN}✓ PASS${NC}: $test_name"
            ((TESTS_PASSED++))
            ;;
        FAIL)
            echo -e "${RED}✗ FAIL${NC}: $test_name"
            [[ -n "$message" ]] && echo -e "  ${RED}→${NC} $message"
            ((TESTS_FAILED++))
            ;;
        SKIP)
            echo -e "${YELLOW}⊘ SKIP${NC}: $test_name"
            [[ -n "$message" ]] && echo -e "  ${YELLOW}→${NC} $message"
            ((TESTS_SKIPPED++))
            ;;
        INFO)
            echo -e "${BLUE}ℹ INFO${NC}: $test_name"
            [[ -n "$message" ]] && echo -e "  ${BLUE}→${NC} $message"
            ;;
    esac
}

# Function to detect OS
detect_os() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    else
        echo "unknown"
    fi
}

# Function to get CPU count
get_cpu_count() {
    local os=$(detect_os)
    if [[ "$os" == "macos" ]]; then
        sysctl -n hw.ncpu 2>/dev/null || echo "1"
    else
        nproc 2>/dev/null || echo "1"
    fi
}

echo "========================================"
echo "schism-esmf CMake Instructions Test"
echo "========================================"
echo ""

OS_TYPE=$(detect_os)
print_test INFO "Operating System" "$OS_TYPE"
print_test INFO "CPU Count" "$(get_cpu_count)"
echo ""

# ========================================
# Test 1: Check ESMFMKFILE environment variable
# ========================================
echo "Test 1: ESMFMKFILE Environment Variable"
echo "----------------------------------------"
if [[ -z "${ESMFMKFILE:-}" ]]; then
    print_test FAIL "ESMFMKFILE is set" "Environment variable not defined"
    print_test INFO "Fix" "export ESMFMKFILE=/path/to/esmf.mk"
elif [[ ! -f "$ESMFMKFILE" ]]; then
    print_test FAIL "ESMFMKFILE points to valid file" "File does not exist: $ESMFMKFILE"
else
    print_test PASS "ESMFMKFILE is set" "$ESMFMKFILE"
    
    # Validate esmf.mk content
    if grep -q "ESMF_VERSION_MAJOR" "$ESMFMKFILE" 2>/dev/null; then
        print_test PASS "esmf.mk contains ESMF variables"
        
        # Extract and display ESMF version if available
        if grep -q "ESMF_VERSION_STRING" "$ESMFMKFILE"; then
            version=$(grep "ESMF_VERSION_STRING" "$ESMFMKFILE" | head -1 | cut -d= -f2 | tr -d ' ')
            print_test INFO "ESMF Version" "$version"
        fi
        
        # Check ESMF_COMM (critical for build success)
        if grep -q "ESMF_COMM" "$ESMFMKFILE"; then
            esmf_comm=$(grep "^ESMF_COMM=" "$ESMFMKFILE" | head -1 | cut -d= -f2 | tr -d ' ')
            print_test INFO "ESMF MPI Implementation" "$esmf_comm"
            
            if [[ "$esmf_comm" == "mpiuni" ]] || [[ "$esmf_comm" =~ "nompi" ]]; then
                print_test FAIL "ESMF has real MPI support" "ESMF_COMM=$esmf_comm (stub MPI, not compatible with SCHISM)"
                print_test INFO "Fix (conda)" "mamba install -c conda-forge 'esmf=*=mpi_mpich*'"
                print_test INFO "Fix (spack)" "spack install esmf+mpi"
            elif [[ "$esmf_comm" == "mpich" ]] || [[ "$esmf_comm" == "openmpi" ]]; then
                print_test PASS "ESMF has real MPI support" "$esmf_comm"
            else
                print_test WARN "ESMF MPI implementation unknown" "$esmf_comm"
            fi
        else
            print_test FAIL "ESMF_COMM variable found" "Cannot determine MPI implementation"
        fi
    else
        print_test FAIL "esmf.mk appears invalid" "Missing ESMF_VERSION_MAJOR"
    fi
fi
echo ""

# ========================================
# Test 2: Check SCHISM_BUILD_DIR environment variable
# ========================================
echo "Test 2: SCHISM_BUILD_DIR Environment Variable"
echo "----------------------------------------------"
if [[ -z "${SCHISM_BUILD_DIR:-}" ]]; then
    print_test FAIL "SCHISM_BUILD_DIR is set" "Environment variable not defined"
    print_test INFO "Fix" "export SCHISM_BUILD_DIR=/path/to/schism/build"
elif [[ ! -d "$SCHISM_BUILD_DIR" ]]; then
    print_test FAIL "SCHISM_BUILD_DIR points to valid directory" "Directory does not exist: $SCHISM_BUILD_DIR"
else
    print_test PASS "SCHISM_BUILD_DIR is set" "$SCHISM_BUILD_DIR"
    
    # Check for required subdirectories
    if [[ -d "$SCHISM_BUILD_DIR/include" ]]; then
        print_test PASS "SCHISM include directory exists"
    else
        print_test FAIL "SCHISM include directory exists" "$SCHISM_BUILD_DIR/include not found"
    fi
    
    if [[ -d "$SCHISM_BUILD_DIR/lib" ]]; then
        print_test PASS "SCHISM lib directory exists"
    else
        print_test FAIL "SCHISM lib directory exists" "$SCHISM_BUILD_DIR/lib not found"
    fi
    
    # Check for required libraries
    if [[ -f "$SCHISM_BUILD_DIR/lib/libhydro.a" ]]; then
        print_test PASS "libhydro.a exists"
        print_test INFO "Library size" "$(du -h "$SCHISM_BUILD_DIR/lib/libhydro.a" | cut -f1)"
    else
        print_test FAIL "libhydro.a exists" "$SCHISM_BUILD_DIR/lib/libhydro.a not found"
    fi
    
    if [[ -f "$SCHISM_BUILD_DIR/lib/libcore.a" ]]; then
        print_test PASS "libcore.a exists"
    else
        print_test FAIL "libcore.a exists" "Core library not found"
    fi
    
    # Check for additional required dependencies
    required_libs=("libturbulence.a" "libyaml.a")
    optional_libs=("libparmetis.a" "libmetis.a")
    
    missing_required=()
    for lib in "${required_libs[@]}"; do
        if [[ -f "$SCHISM_BUILD_DIR/lib/$lib" ]]; then
            print_test PASS "$lib exists"
        else
            print_test FAIL "$lib exists" "Required dependency not found"
            missing_required+=("$lib")
        fi
    done
    
    missing_optional=()
    for lib in "${optional_libs[@]}"; do
        if [[ -f "$SCHISM_BUILD_DIR/lib/$lib" ]]; then
            print_test PASS "$lib exists"
        else
            print_test WARN "$lib exists" "Optional library not found (may cause link errors)"
            missing_optional+=("$lib")
        fi
    done
    
    if [[ ${#missing_required[@]} -gt 0 ]] || [[ ${#missing_optional[@]} -gt 0 ]]; then
        print_test INFO "Available SCHISM libraries" "$(ls -1 $SCHISM_BUILD_DIR/lib/*.a 2>/dev/null | xargs -n1 basename || echo 'none')"
        if [[ ${#missing_required[@]} -gt 0 ]]; then
            print_test INFO "Fix" "Rebuild SCHISM with: cmake ../src && make"
        fi
    fi
fi
echo ""

# ========================================
# Test 2.5: Check MPI Libraries
# ========================================
echo "Test 2.5: MPI Library Availability"
echo "-----------------------------------"

# Try to find MPI libraries
mpi_found=false
mpi_impl=""

if command -v mpirun &> /dev/null || command -v mpiexec &> /dev/null; then
    print_test PASS "MPI runtime found"
    mpi_found=true
    
    # Detect MPI implementation
    if command -v mpichversion &> /dev/null; then
        mpi_impl="MPICH"
        mpi_version=$(mpichversion 2>/dev/null || echo "unknown")
        print_test INFO "MPI Implementation" "$mpi_impl $mpi_version"
    elif command -v ompi_info &> /dev/null; then
        mpi_impl="OpenMPI"
        mpi_version=$(ompi_info --version 2>/dev/null | head -1 || echo "unknown")
        print_test INFO "MPI Implementation" "$mpi_impl $mpi_version"
    else
        print_test INFO "MPI Implementation" "Unknown (mpirun found)"
    fi
else
    print_test WARN "MPI runtime found" "mpirun/mpiexec not in PATH"
fi

# Check for MPI libraries in common locations
mpi_lib_found=false
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    if ls "$CONDA_PREFIX/lib"/libmpi*.dylib &>/dev/null || ls "$CONDA_PREFIX/lib"/libmpi*.so &>/dev/null; then
        print_test PASS "MPI libraries found in conda environment"
        mpi_lib_found=true
        print_test INFO "MPI library path" "$CONDA_PREFIX/lib"
    fi
elif [[ -n "${SPACK_ROOT:-}" ]]; then
    # Spack environment - MPI should be in PATH
    print_test INFO "Spack environment detected" "MPI from spack load"
    mpi_lib_found=true
else
    # Check system locations
    for libdir in /usr/lib /usr/local/lib /opt/local/lib /opt/homebrew/lib; do
        if [[ -d "$libdir" ]] && (ls "$libdir"/libmpi*.so &>/dev/null || ls "$libdir"/libmpi*.dylib &>/dev/null); then
            print_test PASS "MPI libraries found" "$libdir"
            mpi_lib_found=true
            break
        fi
    done
fi

if ! $mpi_lib_found && ! $mpi_found; then
    print_test FAIL "MPI libraries available" "No MPI installation detected"
    print_test INFO "Fix (conda)" "mamba install -c conda-forge mpich or openmpi"
    print_test INFO "Fix (system)" "Install OpenMPI or MPICH via package manager"
fi

echo ""

# ========================================
# Test 3: Check Fortran compiler
# ========================================
echo "Test 3: Fortran Compiler Configuration"
echo "---------------------------------------"

if [[ -n "${CMAKE_Fortran_COMPILER:-}" ]]; then
    FC_COMPILER="$CMAKE_Fortran_COMPILER"
    print_test INFO "CMAKE_Fortran_COMPILER is set" "$FC_COMPILER"
else
    # Try to detect common Fortran compilers
    if command -v mpif90 &> /dev/null; then
        FC_COMPILER="mpif90"
        print_test INFO "Detected MPI Fortran compiler" "mpif90"
    elif command -v gfortran &> /dev/null; then
        FC_COMPILER="gfortran"
        print_test INFO "Detected Fortran compiler" "gfortran"
    else
        FC_COMPILER=""
        print_test FAIL "Fortran compiler found" "No compiler detected"
    fi
fi

if [[ -n "$FC_COMPILER" ]] && command -v "$FC_COMPILER" &> /dev/null; then
    print_test PASS "Fortran compiler available" "$FC_COMPILER"
    
    # Check if it's an MPI wrapper
    if "$FC_COMPILER" --version 2>&1 | grep -iq "mpi\|openmpi\|mpich"; then
        print_test PASS "Compiler is MPI-aware"
    elif "$FC_COMPILER" -show 2>&1 | grep -iq "mpi"; then
        print_test PASS "Compiler is MPI wrapper"
    else
        print_test FAIL "Compiler is MPI-aware" "Use mpif90, esmpifort, or similar"
        print_test INFO "Fix" "export CMAKE_Fortran_COMPILER=\$(which mpif90)"
    fi
    
    # Get compiler version
    if version=$("$FC_COMPILER" --version 2>/dev/null | head -1); then
        print_test INFO "Compiler version" "$version"
    fi
fi
echo ""

# ========================================
# Test 4: CMake configuration test (dry run)
# ========================================
echo "Test 4: CMake Configuration (Dry Run)"
echo "--------------------------------------"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-test-$$"

if [[ -z "${ESMFMKFILE:-}" ]] || [[ -z "${SCHISM_BUILD_DIR:-}" ]]; then
    print_test SKIP "CMake configure" "Prerequisites not met (ESMFMKFILE or SCHISM_BUILD_DIR)"
else
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    print_test INFO "Build directory" "$BUILD_DIR"
    
    # Build CMake command
    CMAKE_CMD="cmake .."
    if [[ "$OS_TYPE" == "macos" ]] && [[ -z "${CMAKE_Fortran_COMPILER:-}" ]]; then
        if command -v gfortran &> /dev/null; then
            CMAKE_CMD="cmake .. -DCMAKE_Fortran_COMPILER=$(which gfortran)"
            if command -v gcc &> /dev/null; then
                CMAKE_CMD="$CMAKE_CMD -DCMAKE_C_COMPILER=$(which gcc)"
            fi
            print_test INFO "Using macOS-specific CMake flags"
        fi
    fi
    
    print_test INFO "CMake command" "$CMAKE_CMD"
    
    # Run CMake configure
    if CMAKE_OUTPUT=$(eval "$CMAKE_CMD" 2>&1); then
        print_test PASS "CMake configure succeeded"
        
        # Check for warnings
        if echo "$CMAKE_OUTPUT" | grep -i "warning" > /dev/null; then
            warning_count=$(echo "$CMAKE_OUTPUT" | grep -ic "warning")
            print_test INFO "CMake warnings" "$warning_count warning(s) found"
            echo "$CMAKE_OUTPUT" | grep -i "warning" | head -5 | while read -r line; do
                echo "  $line"
            done
        else
            print_test PASS "No CMake warnings"
        fi
        
        # Check for ESMF target
        if [[ -f "CMakeCache.txt" ]] && grep -q "ESMF_FOUND:BOOL=TRUE" CMakeCache.txt; then
            print_test PASS "ESMF found by CMake"
        else
            print_test FAIL "ESMF found by CMake" "Check ESMFMKFILE"
        fi
        
        # Check for SCHISM libraries
        if grep -q "SCHISM_BUILD_DIR" CMakeCache.txt; then
            print_test PASS "SCHISM_BUILD_DIR set in CMake"
        fi
        
    else
        print_test FAIL "CMake configure succeeded"
        echo "$CMAKE_OUTPUT" | tail -20
    fi
    
    # Cleanup test build directory
    cd "$REPO_ROOT"
    rm -rf "$BUILD_DIR"
fi
echo ""

# ========================================
# Test 5: Check repository structure
# ========================================
echo "Test 5: Repository Structure"
echo "-----------------------------"

required_dirs=("src" "cmake" "example" "scripts")
for dir in "${required_dirs[@]}"; do
    if [[ -d "$REPO_ROOT/$dir" ]]; then
        print_test PASS "Directory exists: $dir"
    else
        print_test FAIL "Directory exists: $dir"
    fi
done

required_files=("CMakeLists.txt" "cmake/FindESMF.cmake" "src/main_esmf.F90" "src/main_nuopc.F90")
for file in "${required_files[@]}"; do
    if [[ -f "$REPO_ROOT/$file" ]]; then
        print_test PASS "File exists: $file"
    else
        print_test FAIL "File exists: $file"
    fi
done
echo ""

# ========================================
# Summary
# ========================================
echo "========================================"
echo "Test Summary"
echo "========================================"
echo -e "${GREEN}Passed:${NC}  $TESTS_PASSED"
echo -e "${RED}Failed:${NC}  $TESTS_FAILED"
echo -e "${YELLOW}Skipped:${NC} $TESTS_SKIPPED"
echo "========================================"

if [[ $TESTS_FAILED -eq 0 ]]; then
    echo -e "${GREEN}All tests passed!${NC}"
    echo ""
    echo "You can now build with:"
    echo "  mkdir -p build && cd build"
    if [[ "$OS_TYPE" == "macos" ]]; then
        echo "  cmake .. -DCMAKE_Fortran_COMPILER=\$(which gfortran) -DCMAKE_C_COMPILER=\$(which gcc)"
        echo "  cmake --build . -- -j\$(sysctl -n hw.ncpu)"
    else
        echo "  cmake .."
        echo "  cmake --build . -- -j\$(nproc)"
    fi
    exit 0
else
    echo -e "${RED}Some tests failed. Please fix the issues above.${NC}"
    exit 1
fi
