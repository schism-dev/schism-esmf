cmake_minimum_required(VERSION 3.16)

# Option to configure a documentation-only build (no Fortran/C targets)
option(DOCS_ONLY "Configure to build documentation only (skip Fortran/C targets)" OFF)

if(DOCS_ONLY)
  # Avoid enabling compiled languages in docs-only mode
  project(SCHISM_ESMF_Interface LANGUAGES NONE)
else()
  project(SCHISM_ESMF_Interface LANGUAGES Fortran C)
  enable_language(Fortran)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Build options
option(USE_PDAF "Enable PDAF data assimilation support" OFF)
option(USE_MOSSCO "Enable MOSSCO coupling framework support" OFF)
option(SCHISM_REQUIRE_PARMETIS "Require ParMETIS library (fail if not found)" OFF)
option(SCHISM_REQUIRE_TURBULENCE "Require turbulence library (fail if not found)" ON)
option(SCHISM_REQUIRE_YAML "Require YAML library (fail if not found)" ON)
option(BUILD_DOCS "Build HTML documentation with MkDocs (requires Python, mkdocs, mkdocs-material)" OFF)

# Force docs build in DOCS_ONLY mode
if(DOCS_ONLY)
  set(BUILD_DOCS ON CACHE BOOL "Build HTML documentation with MkDocs" FORCE)
endif()

message(STATUS "Build options:")
message(STATUS "  USE_PDAF: ${USE_PDAF}")
message(STATUS "  USE_MOSSCO: ${USE_MOSSCO}")
message(STATUS "  SCHISM_REQUIRE_PARMETIS: ${SCHISM_REQUIRE_PARMETIS}")
message(STATUS "  SCHISM_REQUIRE_TURBULENCE: ${SCHISM_REQUIRE_TURBULENCE}")
message(STATUS "  SCHISM_REQUIRE_YAML: ${SCHISM_REQUIRE_YAML}")
message(STATUS "  BUILD_DOCS: ${BUILD_DOCS}")
message(STATUS "  DOCS_ONLY: ${DOCS_ONLY}")
message(STATUS "")

if(NOT DOCS_ONLY)
  # Include GNUInstallDirs for standard installation directories
  include(GNUInstallDirs)
endif()

if(NOT DOCS_ONLY)
  # Setup ESMF dependency
  include(cmake/ConfigureESMF.cmake)
endif()

if(NOT DOCS_ONLY)
  # Find MPI (required for ParMETIS and SCHISM)
  find_package(MPI REQUIRED COMPONENTS Fortran C)
  if(MPI_FOUND)
    message(STATUS "MPI Fortran libraries: ${MPI_Fortran_LIBRARIES}")
    message(STATUS "MPI C libraries: ${MPI_C_LIBRARIES}")
  endif()
endif()


# SCHISM Dependency
# SCHISM is expected to be pre-built. Its location is specified by SCHISM_BUILD_DIR.
# MPI, NetCDF, LAPACK, and BLAS are expected to be handled by the ESMF dependency.
# If ESMF is not found or doesn't correctly export them, issues might arise.
# Users might need to ensure their ESMF build correctly wraps these or install them separately
# and adjust CMake configurations if direct linking becomes necessary.
# Allow both environment variable and CMake variable override
if(NOT DOCS_ONLY)
  set(SCHISM_BUILD_DIR "$ENV{SCHISM_BUILD_DIR}" CACHE PATH "Path to pre-built SCHISM installation")

  if(NOT SCHISM_BUILD_DIR OR NOT EXISTS "${SCHISM_BUILD_DIR}")
    message(FATAL_ERROR
      "SCHISM_BUILD_DIR must be set to a valid SCHISM build directory. "
      "Set it via environment variable or CMake -DSCHISM_BUILD_DIR=/path/to/schism"
    )
  endif()

  message(STATUS "Using SCHISM build directory: ${SCHISM_BUILD_DIR}")
endif()

if(NOT DOCS_ONLY)
  # Define SCHISM components
  set(SCHISM_INCLUDE_DIR "${SCHISM_BUILD_DIR}/include")
  set(SCHISM_LIBRARY_DIR "${SCHISM_BUILD_DIR}/lib")

  # Verify required directories exist
  foreach(DIR ${SCHISM_INCLUDE_DIR} ${SCHISM_LIBRARY_DIR})
    if(NOT EXISTS "${DIR}")
      message(FATAL_ERROR "Required SCHISM directory not found: ${DIR}")
    endif()
  endforeach()
endif()

# Import SCHISM libraries
# ===========================================================================
# Core libraries (required: core, hydro; recommended: turbulence, yaml; optional: parmetis, metis, util)
if(NOT DOCS_ONLY)
  set(SCHISM_LIBRARIES "")
  set(SCHISM_REQUIRED_LIBS core hydro)
  set(SCHISM_RECOMMENDED_LIBS)
  set(SCHISM_OPTIONAL_LIBS util)
endif()

# Configure required/recommended based on options
if(NOT DOCS_ONLY AND SCHISM_REQUIRE_TURBULENCE)
  list(APPEND SCHISM_REQUIRED_LIBS turbulence)
elseif(NOT DOCS_ONLY)
  list(APPEND SCHISM_RECOMMENDED_LIBS turbulence)
endif()

if(NOT DOCS_ONLY AND SCHISM_REQUIRE_YAML)
  list(APPEND SCHISM_REQUIRED_LIBS yaml)
elseif(NOT DOCS_ONLY)
  list(APPEND SCHISM_RECOMMENDED_LIBS yaml)
endif()

if(NOT DOCS_ONLY AND SCHISM_REQUIRE_PARMETIS)
  list(APPEND SCHISM_REQUIRED_LIBS parmetis metis)
elseif(NOT DOCS_ONLY)
  list(APPEND SCHISM_OPTIONAL_LIBS parmetis metis)
endif()

if(NOT DOCS_ONLY)
  message(STATUS "Discovering SCHISM libraries in ${SCHISM_LIBRARY_DIR}")
endif()

foreach(LIB ${SCHISM_REQUIRED_LIBS} ${SCHISM_RECOMMENDED_LIBS} ${SCHISM_OPTIONAL_LIBS})
  if(DOCS_ONLY)
    break()
  endif()
  find_library(SCHISM_${LIB}_PATH
    NAMES ${LIB}
    HINTS ${SCHISM_LIBRARY_DIR}
    NO_DEFAULT_PATH
  )
  if(SCHISM_${LIB}_PATH)
    list(APPEND SCHISM_LIBRARIES ${SCHISM_${LIB}_PATH})
    # Get library size for informational purposes
    file(SIZE ${SCHISM_${LIB}_PATH} LIB_SIZE)
    math(EXPR LIB_SIZE_KB "${LIB_SIZE} / 1024")
    message(STATUS "  ✓ Found SCHISM library: lib${LIB}.a (${LIB_SIZE_KB} KB)")
  else()
    # Check if library is required, recommended, or optional
    if(LIB IN_LIST SCHISM_REQUIRED_LIBS)
      message(FATAL_ERROR "  ✗ SCHISM ${LIB} library not found (REQUIRED) in ${SCHISM_LIBRARY_DIR}\n"
                          "     Hint: Rebuild SCHISM or disable requirement with -DSCHISM_REQUIRE_${LIB}=OFF")
    elseif(LIB IN_LIST SCHISM_RECOMMENDED_LIBS)
      message(WARNING "  ⚠ SCHISM ${LIB} library not found (RECOMMENDED) - may cause link errors\n"
                      "     Hint: Rebuild SCHISM or mark as optional with -DSCHISM_REQUIRE_${LIB}=OFF")
    else()
      message(STATUS "  ○ SCHISM ${LIB} library not found (optional)")
    endif()
  endif()
endforeach()

if(NOT DOCS_ONLY)
  # Find and link MPI libraries (required for SCHISM)
  message(STATUS "Discovering MPI libraries")
  find_package(MPI REQUIRED COMPONENTS Fortran)
  if(MPI_Fortran_FOUND)
    list(APPEND SCHISM_LIBRARIES ${MPI_Fortran_LIBRARIES})
    message(STATUS "  ✓ MPI Fortran libraries: ${MPI_Fortran_LIBRARIES}")
    message(STATUS "  ✓ MPI include directories: ${MPI_Fortran_INCLUDE_DIRS}")
    
    # Try to detect MPI implementation
    if(MPI_Fortran_LIBRARIES MATCHES "mpich")
      message(STATUS "  → MPI implementation: MPICH")
    elseif(MPI_Fortran_LIBRARIES MATCHES "openmpi")
      message(STATUS "  → MPI implementation: OpenMPI")
    else()
      message(STATUS "  → MPI implementation: Unknown")
    endif()
  else()
    message(FATAL_ERROR "  ✗ MPI Fortran libraries not found (REQUIRED)")
  endif()
endif()

if(NOT DOCS_ONLY)
  message(STATUS "Total SCHISM libraries found: ${SCHISM_LIBRARIES}")
  message(STATUS "")
endif()

# # Handle FABM dependency
# include(cmake/fabm.cmake)

if(NOT DOCS_ONLY)
  # Add src subdirectory
  add_subdirectory(src)
endif()

# Add documentation subdirectory
add_subdirectory(doc)

# add_subdirectory(example)

if(NOT DOCS_ONLY)
  install(TARGETS main_esmf main_nuopc
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
endif()

# Install the configured schism.mk file
if(NOT DOCS_ONLY)
  set(SCHISM_MK_CONFIGURED_PATH "${CMAKE_BINARY_DIR}/schism.mk")
  # Check for existence at configure time for the message.
  # install(FILES) will check again at install time.
  if(EXISTS "${SCHISM_MK_CONFIGURED_PATH}")
    install(FILES "${SCHISM_MK_CONFIGURED_PATH}"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}" # Matches original Makefile's DESTDIR for schism.mk
      RENAME schism.mk
    )
    message(STATUS "Installation rule added for schism.mk to ${CMAKE_INSTALL_LIBDIR}")
  else()
    # This message is for configure time. If schism.mk is generated by a build step
    # that hasn't run yet during configure, it might show this warning but still be installable.
    message(STATUS "Configured schism.mk not found at ${SCHISM_MK_CONFIGURED_PATH} during CMake configure phase. If it's generated by a build command, it might still be installed correctly.")
  endif()
endif()
