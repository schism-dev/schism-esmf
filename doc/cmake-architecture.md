# CMake Build System Architecture

This document explains the CMake build system design for schism-esmf, focusing on the patterns and conventions used to handle complex Fortran module dependencies.

## Overview

The schism-esmf build system creates:
- **Libraries**: `libschism_interface_common.a`, `libschism_esmf_interface.a`, `libschism_nuopc_interface.a`, `libschism_model_libs.a`, `libschism_driver_libs.a`
- **Executables**: `main_esmf`, `main_nuopc`

## Directory Structure

```
CMakeLists.txt           # Top-level: ESMF/SCHISM discovery, installation
cmake/
  FindESMF.cmake         # Custom module to parse esmf.mk
  ConfigureESMF.cmake    # ESMF validation and warnings
  compilers.cmake        # Compiler-specific flags (optional)
src/
  CMakeLists.txt         # Main executables, module directory setup
  schism/
    CMakeLists.txt       # Interface libraries (ESMF/NUOPC caps)
  model/
    CMakeLists.txt       # Atmosphere CMI implementations
  driver/
    CMakeLists.txt       # Toplevel drivers
  PDAF_bindings/         # Optional: data assimilation
  mossco/                # Optional: MOSSCO integration
```

## Key Architectural Patterns

### 1. Shared Source Pattern (`schism_interface_common`)

**Problem**: Both `schism_esmf_interface` and `schism_nuopc_interface` need to compile the same source files (`schism_bmi.F90`, `schism_esmf_util.F90`). When building in parallel, both targets try to create identical `.mod` files simultaneously, causing race conditions:

```
Error copying Fortran module "schism_bmi.mod"
Cannot rename module file 'schism_bmi.mod0' to 'schism_bmi.mod'
```

**Solution**: Create a separate library for common sources that both interface libraries depend on:

```cmake
# src/schism/CMakeLists.txt

# Compile common sources exactly once
add_library(schism_interface_common STATIC
    schism_bmi.F90
    schism_esmf_util.F90
)

set_target_properties(schism_interface_common PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# ESMF interface depends on common library
add_library(schism_esmf_interface STATIC
    schism_esmf_cap.F90
)
target_link_libraries(schism_esmf_interface PUBLIC schism_interface_common)

# NUOPC interface depends on common library
add_library(schism_nuopc_interface STATIC
    schism_nuopc_util.F90
    schism_nuopc_cap.F90
)
target_link_libraries(schism_nuopc_interface PUBLIC schism_interface_common)
```

**Benefits**:
- ✅ Common sources compiled exactly once
- ✅ Modules generated once and shared
- ✅ No parallel build race conditions
- ✅ Clear dependency graph
- ✅ Faster builds (less compilation)

### 2. Centralized Fortran Module Directory

**Problem**: Fortran `.mod` files must be found by consumers. By default, CMake generates them in each target's build directory, requiring complex include path management.

**Solution**: Set `CMAKE_Fortran_MODULE_DIRECTORY` globally and per-target:

```cmake
# src/CMakeLists.txt

# Set centralized module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Ensure all targets use it
add_subdirectory(schism)
add_subdirectory(model)
add_subdirectory(driver)

# Make modules visible to executables
target_include_directories(main_esmf PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(main_nuopc PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
```

**Per-target enforcement**:
```cmake
set_target_properties(schism_interface_common PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
```

**Benefits**:
- ✅ Single location for all `.mod` files
- ✅ Simplified include paths
- ✅ Easier debugging (one place to check)
- ✅ Consistent module handling

### 3. Dependency Chain

The build system links libraries in a specific order to resolve symbols:

```
main_esmf/main_nuopc
  ↓
schism_driver_libs (toplevel_schism_atm, driver_schism)
  ↓
schism_model_libs (atmosphere_cmi_esmf, atmosphere_cmi_nuopc)
  ↓
schism_esmf_interface/schism_nuopc_interface (ESMF caps)
  ↓
schism_interface_common (BMI, utilities)
  ↓
SCHISM_LIBRARIES (core, hydro, turbulence, yaml, parmetis, metis)
  ↓
MPI_Fortran_LIBRARIES (mpi_mpich, pmpi_mpich)
  ↓
ESMF::ESMF (libesmf, pio, netcdf)
```

**Implementation in `src/CMakeLists.txt`**:
```cmake
# Add subdirectories first (creates targets)
add_subdirectory(schism)
add_subdirectory(model)
add_subdirectory(driver)

# Link executables to libraries (order matters for static libs)
target_link_libraries(main_esmf PRIVATE
    schism_driver_libs
    schism_model_libs
    schism_esmf_interface
    ESMF::ESMF
)

target_link_libraries(main_nuopc PRIVATE
    schism_driver_libs
    schism_model_libs
    schism_nuopc_interface
    ESMF::ESMF
)
```

**Why order matters**: With static libraries, the linker processes libraries left-to-right. Symbols are only resolved from libraries that come *after* the object files/libraries that need them. This is why `ESMF::ESMF` comes last — it provides symbols needed by all the schism libraries.

### 4. SCHISM Dependency Discovery

**Problem**: SCHISM is built separately and provides multiple libraries. CMake needs to discover and link all of them.

**Solution**: Loop through known library names and accumulate found libraries:

```cmake
# Top-level CMakeLists.txt

set(SCHISM_LIBRARY_DIR "${SCHISM_BUILD_DIR}/lib")
set(SCHISM_INCLUDE_DIR "${SCHISM_BUILD_DIR}/include")

set(SCHISM_LIBRARIES "")
foreach(LIB core hydro turbulence yaml parmetis metis)
  find_library(SCHISM_${LIB}_PATH
    NAMES ${LIB}
    HINTS ${SCHISM_LIBRARY_DIR}
    NO_DEFAULT_PATH
  )
  if(SCHISM_${LIB}_PATH)
    list(APPEND SCHISM_LIBRARIES ${SCHISM_${LIB}_PATH})
    message(STATUS "Found SCHISM library: ${SCHISM_${LIB}_PATH}")
  else()
    message(WARNING "SCHISM library ${LIB} not found in ${SCHISM_LIBRARY_DIR}")
  endif()
endforeach()

# Also need MPI libraries
find_package(MPI REQUIRED COMPONENTS Fortran)
if(MPI_Fortran_FOUND)
  list(APPEND SCHISM_LIBRARIES ${MPI_Fortran_LIBRARIES})
endif()
```

**Benefits**:
- ✅ Automatic discovery of available libraries
- ✅ Clear warning for missing libraries
- ✅ Flexible (works with different SCHISM builds)
- ✅ MPI properly linked

### 5. ESMF Configuration

**Problem**: ESMF is installed externally and provides configuration via `esmf.mk` makefile fragment, not CMake config.

**Solution**: Custom `FindESMF.cmake` module parses `esmf.mk`:

```cmake
# cmake/FindESMF.cmake

# Find esmf.mk
find_file(ESMF_MK_FILE
  NAMES esmf.mk
  HINTS ${ESMFMKFILE} ENV ESMFMKFILE
  PATH_SUFFIXES lib
)

# Parse variables from esmf.mk
file(STRINGS "${ESMF_MK_FILE}" _esmf_version_line REGEX "^ESMF_VERSION_STRING=")
string(REGEX REPLACE "^ESMF_VERSION_STRING=([0-9.]+).*" "\\1" ESMF_VERSION "${_esmf_version_line}")

# Parse compile options (must use separate_arguments!)
set(_esmf_opts_string "${ESMF_F90COMPILEOPTS} ${ESMF_F90COMPILEFREECPP}")
separate_arguments(_esmf_opts_list UNIX_COMMAND "${_esmf_opts_string}")
set(ESMF_COMPILE_OPTIONS ${_esmf_opts_list})

# Create imported target
add_library(ESMF::ESMF INTERFACE IMPORTED)
set_target_properties(ESMF::ESMF PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${ESMF_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${ESMF_LIBRARIES}"
  INTERFACE_COMPILE_OPTIONS "${ESMF_COMPILE_OPTIONS}"
)
```

**Critical detail**: Use `separate_arguments()` to properly parse space-separated compiler flags. Passing them as a single string causes compilation errors.

## Build Process Flow

1. **Configuration** (`cmake ..`):
   - Find ESMF via `$ESMFMKFILE`
   - Parse `esmf.mk` and create `ESMF::ESMF` target
   - Find SCHISM libraries in `$SCHISM_BUILD_DIR/lib`
   - Discover MPI libraries
   - Set up module directory
   - Configure subdirectories

2. **Generation** (`cmake --build .`):
   - Generate build files (Makefiles or Ninja)
   - Create `build/modules/` directory

3. **Compilation**:
   - Build `schism_interface_common` first (15%)
   - Build `schism_nuopc_interface` and `schism_esmf_interface` in parallel (31%, 42%)
   - Build `schism_model_libs` (52%)
   - Build `schism_driver_libs` (84%)
   - Compile `main_esmf.F90` and `main_nuopc.F90` (92%)
   - Link executables (100%)

4. **Linking**:
   - Linker processes libraries in dependency order
   - Resolves symbols from SCHISM, MPI, and ESMF
   - Creates final executables

## Common Pitfalls and Solutions

### Pitfall: "Module file not found"

**Cause**: Target doesn't depend on library that creates the module.

**Fix**: Use `target_link_libraries()` to create dependency:
```cmake
target_link_libraries(schism_model_libs PUBLIC schism_esmf_interface)
```

### Pitfall: "Undefined symbols" at link time

**Cause**: Library not linked or wrong link order.

**Fix**: Add missing library to `SCHISM_LIBRARIES` or `target_link_libraries()`.

### Pitfall: Parallel build fails with module errors

**Cause**: Multiple targets compiling same source files.

**Fix**: Extract common sources into separate library (see Pattern #1).

### Pitfall: "Unrecognized command-line option"

**Cause**: ESMF compile options passed as single quoted string.

**Fix**: Use `separate_arguments()` in FindESMF.cmake.

## Testing the Build System

Run the validation script:
```bash
.github/test-cmake-instructions.sh
```

This checks:
- ESMF configuration and MPI support
- SCHISM library availability
- Compiler configuration
- CMake configuration validity

## Extending the Build System

### Adding a New Library Target

```cmake
add_library(my_new_lib STATIC
    source1.F90
    source2.F90
)

# Set module directory
set_target_properties(my_new_lib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Include modules from other targets
target_include_directories(my_new_lib PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Link dependencies
target_link_libraries(my_new_lib PUBLIC
    schism_interface_common
    ESMF::ESMF
)
```

### Adding an Optional Dependency

```cmake
if(USE_OPTIONAL_LIB)
    find_package(OptionalLib)
    if(OptionalLib_FOUND)
        target_link_libraries(my_target PUBLIC OptionalLib::OptionalLib)
        target_compile_definitions(my_target PRIVATE USE_OPTIONAL_LIB)
    endif()
endif()
```

## References

- **CMake Fortran Support**: https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#fortran-modules
- **ESMF Documentation**: https://earthsystemmodeling.org/docs/
- **SCHISM Manual**: https://ccrm.vims.edu/schismweb/

## Revision History

- **2025-11**: Initial architecture documentation
- Added shared source pattern for parallel build fixes
- Documented SCHISM dependency chain
- Added MPI library discovery
