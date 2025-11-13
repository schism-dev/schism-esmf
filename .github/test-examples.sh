#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
#
# Example: How to use the CMake test framework

# This script demonstrates typical usage patterns for the test framework

echo "=========================================="
echo "CMake Test Framework - Usage Examples"
echo "=========================================="
echo ""

# Example 1: Basic test (assumes env vars are set)
echo "Example 1: Basic Test"
echo "---------------------"
echo "Prerequisites: ESMFMKFILE and SCHISM_BUILD_DIR must be set"
echo ""
echo "Command:"
echo "  .github/test-cmake-instructions.sh"
echo ""

# Example 2: Set env vars inline
echo "Example 2: Set Environment Variables Inline"
echo "--------------------------------------------"
echo "Command:"
echo "  ESMFMKFILE=/opt/esmf/lib/esmf.mk \\"
echo "  SCHISM_BUILD_DIR=/path/to/schism/build \\"
echo "  .github/test-cmake-instructions.sh"
echo ""

# Example 3: macOS with explicit compilers
echo "Example 3: macOS with Homebrew Compilers"
echo "-----------------------------------------"
echo "Command:"
echo "  export CMAKE_Fortran_COMPILER=\$(brew --prefix gcc)/bin/gfortran"
echo "  export CMAKE_C_COMPILER=\$(brew --prefix gcc)/bin/gcc"
echo "  export ESMFMKFILE=/opt/esmf/lib/esmf.mk"
echo "  export SCHISM_BUILD_DIR=~/schism/build"
echo "  .github/test-cmake-instructions.sh"
echo ""

# Example 4: Run only environment checks (skip CMake config)
echo "Example 4: Quick Environment Check"
echo "-----------------------------------"
echo "To check only env vars without running CMake config,"
echo "temporarily unset one of the required variables:"
echo ""
echo "Command:"
echo "  SCHISM_BUILD_DIR=/nonexistent \\"
echo "  .github/test-cmake-instructions.sh"
echo ""

# Example 5: Integration with build script
echo "Example 5: Integration with Build Script"
echo "-----------------------------------------"
echo "#!/usr/bin/env bash"
echo "set -e"
echo ""
echo "# Validate environment before building"
echo ".github/test-cmake-instructions.sh"
echo ""
echo "# If tests pass, proceed with build"
echo "mkdir -p build && cd build"
echo "cmake .."
echo "cmake --build . -- -j\$(nproc)"
echo ""

# Example 6: CI/CD usage
echo "Example 6: CI/CD Integration"
echo "-----------------------------"
echo "In your CI config file:"
echo ""
echo "GitHub Actions:"
echo "  - name: Validate build environment"
echo "    run: .github/test-cmake-instructions.sh"
echo "    env:"
echo "      ESMFMKFILE: \${{ env.ESMFMKFILE }}"
echo "      SCHISM_BUILD_DIR: \${{ env.SCHISM_BUILD_DIR }}"
echo ""
echo "GitLab CI:"
echo "  test:"
echo "    script:"
echo "      - .github/test-cmake-instructions.sh"
echo ""

echo "=========================================="
echo "For more details, see:"
echo "  .github/TEST_README.md"
echo "=========================================="
