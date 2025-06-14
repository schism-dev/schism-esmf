cmake_minimum_required(VERSION 3.16)
project(SCHISM_ESMF_Interface LANGUAGES Fortran C)

enable_language(Fortran)

# Add more CMake code in the subsequent steps.

# Find MPI
find_package(MPI REQUIRED)

# Find NetCDF
find_package(NetCDF REQUIRED)

# Find LAPACK
find_package(LAPACK REQUIRED)

# Find BLAS
find_package(BLAS REQUIRED)

# Include directories for NetCDF
if(NetCDF_FOUND)
  include_directories(${NetCDF_INCLUDE_DIRS})
  link_libraries(${NetCDF_LIBRARIES})
endif()

# Include directories and libraries for LAPACK and BLAS
if(LAPACK_FOUND)
  include_directories(${LAPACK_INCLUDE_DIRS}) # LAPACK_INCLUDE_DIRS may not always be set, but good practice
  link_libraries(${LAPACK_LIBRARIES})
endif()

if(BLAS_FOUND)
  include_directories(${BLAS_INCLUDE_DIRS}) # BLAS_INCLUDE_DIRS may not always be set
  link_libraries(${BLAS_LIBRARIES})
endif()

# ESMF often requires MPI Fortran support
if (MPI_Fortran_FOUND)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_libraries(${MPI_Fortran_LIBRARIES})
else()
  message(WARNING "MPI Fortran libraries not found, this might cause issues with ESMF.")
endif()

# SCHISM Dependency
# SCHISM is expected to be pre-built. Its location is specified by SCHISM_BUILD_DIR.
set(SCHISM_BUILD_DIR_ENV $ENV{SCHISM_BUILD_DIR})

if(NOT SCHISM_BUILD_DIR_ENV)
  message(FATAL_ERROR "SCHISM_BUILD_DIR environment variable is not set. Please set it to point to your SCHISM build directory.")
else()
  message(STATUS "SCHISM_BUILD_DIR found: ${SCHISM_BUILD_DIR_ENV}")
endif()

set(SCHISM_INCLUDE_DIR "${SCHISM_BUILD_DIR_ENV}/include")
set(SCHISM_LIBRARY_DIR "${SCHISM_BUILD_DIR_ENV}/lib")

# Check if include and lib directories exist
if(NOT EXISTS "${SCHISM_INCLUDE_DIR}")
  message(FATAL_ERROR "SCHISM include directory not found: ${SCHISM_INCLUDE_DIR}")
endif()
if(NOT EXISTS "${SCHISM_LIBRARY_DIR}")
  message(FATAL_ERROR "SCHISM library directory not found: ${SCHISM_LIBRARY_DIR}")
endif()

include_directories(${SCHISM_INCLUDE_DIR})

# Define SCHISM libraries
# Using find_library to confirm existence and get full path
find_library(SCHISM_LIB_CORE_PATH NAMES core NAMES_PER_DIR PATHS "${SCHISM_LIBRARY_DIR}" NO_DEFAULT_PATH)
find_library(SCHISM_LIB_HYDRO_PATH NAMES hydro NAMES_PER_DIR PATHS "${SCHISM_LIBRARY_DIR}" NO_DEFAULT_PATH)
# ParMetis and Metis are often linked with SCHISM, but their presence might be conditional
# We will handle these as separate dependencies if needed, or assume they are part of the SCHISM libs provided
# For now, let's assume they are part of the general SCHISM linking process or handled by ESMF/NUOPC dependencies.
# The Makefile does "AR x" for these from SCHISM_BUILD_DIR/lib into .objs/d, suggesting they are expected there.

if(NOT SCHISM_LIB_CORE_PATH)
  message(FATAL_ERROR "SCHISM library 'core' not found in ${SCHISM_LIBRARY_DIR}")
endif()
if(NOT SCHISM_LIB_HYDRO_PATH)
  message(FATAL_ERROR "SCHISM library 'hydro' not found in ${SCHISM_LIBRARY_DIR}")
endif()

add_library(SCHISM_core STATIC IMPORTED GLOBAL)
set_target_properties(SCHISM_core PROPERTIES IMPORTED_LOCATION "${SCHISM_LIB_CORE_PATH}")

add_library(SCHISM_hydro STATIC IMPORTED GLOBAL)
set_target_properties(SCHISM_hydro PROPERTIES IMPORTED_LOCATION "${SCHISM_LIB_HYDRO_PATH}")

# Group SCHISM libraries for easier linking
set(SCHISM_LIBRARIES SCHISM_core SCHISM_hydro)
message(STATUS "SCHISM libraries found: ${SCHISM_LIBRARIES}")

# Metis and ParMetis (often bundled or direct dependencies of SCHISM)
# These are listed in the root Makefile's LDFLAGS and also extracted from SCHISM_BUILD_DIR/lib
# Let's try to find them in SCHISM_LIBRARY_DIR first.
find_library(SCHISM_LIB_PARMETIS_PATH NAMES parmetis NAMES_PER_DIR PATHS "${SCHISM_LIBRARY_DIR}" NO_DEFAULT_PATH)
find_library(SCHISM_LIB_METIS_PATH NAMES metis NAMES_PER_DIR PATHS "${SCHISM_LIBRARY_DIR}" NO_DEFAULT_PATH)

if(SCHISM_LIB_PARMETIS_PATH)
  add_library(SCHISM_parmetis STATIC IMPORTED GLOBAL)
  set_target_properties(SCHISM_parmetis PROPERTIES IMPORTED_LOCATION "${SCHISM_LIB_PARMETIS_PATH}")
  list(APPEND SCHISM_LIBRARIES SCHISM_parmetis)
  message(STATUS "SCHISM ParMetis library found: ${SCHISM_LIB_PARMETIS_PATH}")
else()
  message(WARNING "SCHISM ParMetis library not found in ${SCHISM_LIBRARY_DIR}. This might be optional or found elsewhere.")
endif()

if(SCHISM_LIB_METIS_PATH)
  add_library(SCHISM_metis STATIC IMPORTED GLOBAL)
  set_target_properties(SCHISM_metis PROPERTIES IMPORTED_LOCATION "${SCHISM_LIB_METIS_PATH}")
  list(APPEND SCHISM_LIBRARIES SCHISM_metis)
  message(STATUS "SCHISM Metis library found: ${SCHISM_LIB_METIS_PATH}")
else()
  message(WARNING "SCHISM Metis library not found in ${SCHISM_LIBRARY_DIR}. This might be optional or found elsewhere.")
endif()

# ESMF Dependency
# Attempt to find ESMF using find_package first
find_package(ESMF QUIET)

if(ESMF_FOUND)
  message(STATUS "ESMF found using find_package.")
  # If ESMF_FOUND, it should have set up targets like ESMF::ESMF or variables
  # Assume it provides ESMF_INCLUDE_DIRS and ESMF_LIBRARIES
  if(TARGET ESMF::ESMF)
    message(STATUS "ESMF target ESMF::ESMF found.")
    set(ESMF_LIBRARIES ESMF::ESMF)
    # Include directories should be handled by the target properties
  else()
    # Fallback if specific target isn't there but ESMF_FOUND is true
    if(DEFINED ESMF_INCLUDE_DIRS)
      message(STATUS "Using ESMF_INCLUDE_DIRS: ${ESMF_INCLUDE_DIRS}")
      include_directories(${ESMF_INCLUDE_DIRS})
    endif()
    if(NOT DEFINED ESMF_LIBRARIES)
      message(WARNING "find_package(ESMF) succeeded but ESMF_LIBRARIES is not defined. Manual configuration might be needed.")
    else()
      message(STATUS "Using ESMF_LIBRARIES: ${ESMF_LIBRARIES}")
    endif()
  endif()
else()
  message(STATUS "ESMF not found via find_package. Attempting to use ESMF environment variables (often set by sourcing esmf.mk).")
  set(ESMFMKFILE_ENV $ENV{ESMFMKFILE})

  if(NOT ESMFMKFILE_ENV)
    message(WARNING "ESMFMKFILE environment variable is not set. This is a strong indicator that ESMF environment variables might be missing.")
  else()
    message(STATUS "ESMFMKFILE found: ${ESMFMKFILE_ENV}. Relying on environment variables set by sourcing this or a similar script.")
  endif()

  # Check for essential ESMF environment variables
  set(ESMF_LIBS_ENV $ENV{ESMF_LIBS})
  set(ESMF_INCLUDE_OPTS_ENV $ENV{ESMF_INCLUDE_OPTS})
  set(ESMF_LINK_OPTS_ENV $ENV{ESMF_LINK_OPTS})
  set(ESMF_F90ESMFLINKLIBS_ENV $ENV{ESMF_F90ESMFLINKLIBS})
  set(ESMF_F90COMPILEOPTS_ENV $ENV{ESMF_F90COMPILEOPTS})
  set(ESMF_F90COMPILECPPFLAGS_ENV $ENV{ESMF_F90COMPILECPPFLAGS})
  set(ESMF_F90COMPILEFREECPP_ENV $ENV{ESMF_F90COMPILEFREECPP})

  set(ESMF_LIBRARIES_TEMP "") # Initialize a temporary list/string for libraries

  if(ESMF_LIBS_ENV)
    message(STATUS "Using ESMF_LIBS from environment: ${ESMF_LIBS_ENV}")
    list(APPEND ESMF_LIBRARIES_TEMP ${ESMF_LIBS_ENV})
  else()
    message(WARNING "ESMF_LIBS environment variable not found.")
  endif()

  if(ESMF_F90ESMFLINKLIBS_ENV)
    message(STATUS "Using ESMF_F90ESMFLINKLIBS from environment: ${ESMF_F90ESMFLINKLIBS_ENV}")
    list(APPEND ESMF_LIBRARIES_TEMP ${ESMF_F90ESMFLINKLIBS_ENV})
  else()
    message(WARNING "ESMF_F90ESMFLINKLIBS environment variable not found. Linking ESMF applications might fail.")
  endif()

  set(ESMF_LIBRARIES ${ESMF_LIBRARIES_TEMP}) # Assign collected libraries to ESMF_LIBRARIES

  if(ESMF_INCLUDE_OPTS_ENV)
    message(STATUS "Using ESMF_INCLUDE_OPTS from environment: ${ESMF_INCLUDE_OPTS_ENV}")
    # Parse -I options and add them to include directories
    string(REGEX MATCHALL "-I[^ ]+" _esmf_includes "${ESMF_INCLUDE_OPTS_ENV}")
    foreach(_inc ${_esmf_includes})
      string(REPLACE "-I" "" _inc_path "${_inc}")
      include_directories(${_inc_path})
      message(STATUS "Added ESMF include directory: ${_inc_path}")
    endforeach()
  else()
    message(WARNING "ESMF_INCLUDE_OPTS environment variable not found. ESMF includes might be missing.")
  endif()

  if(ESMF_F90COMPILEOPTS_ENV)
      message(STATUS "Adding ESMF compile options from ESMF_F90COMPILEOPTS: ${ESMF_F90COMPILEOPTS_ENV}")
      add_compile_options(${ESMF_F90COMPILEOPTS_ENV})
  endif()

  if(ESMF_LINK_OPTS_ENV)
      message(STATUS "Adding ESMF link options from ESMF_LINK_OPTS: ${ESMF_LINK_OPTS_ENV}")
      # These are linker flags like -L paths or rpaths.
      # Add them to linker flags.
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ESMF_LINK_OPTS_ENV}")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ESMF_LINK_OPTS_ENV}")
      set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ESMF_LINK_OPTS_ENV}")
  endif()

  if(ESMF_F90COMPILECPPFLAGS_ENV)
      message(STATUS "Adding ESMF CPP definitions from ESMF_F90COMPILECPPFLAGS: ${ESMF_F90COMPILECPPFLAGS_ENV}")
      # Parse -D options and add them as compile definitions
      string(REGEX MATCHALL "-D[^ ]+" _esmf_defs "${ESMF_F90COMPILECPPFLAGS_ENV}")
      foreach(_def ${_esmf_defs})
          string(REPLACE "-D" "" _def_val "${_def}")
          add_compile_definitions(${_def_val})
          message(STATUS "Added ESMF compile definition: ${_def_val}")
      endforeach()
  endif()

  if(ESMF_F90COMPILEFREECPP_ENV)
      message(STATUS "Adding ESMF free/fixed form flags from ESMF_F90COMPILEFREECPP: ${ESMF_F90COMPILEFREECPP_ENV}")
      add_compile_options(${ESMF_F90COMPILEFREECPP_ENV})
  endif()

  if(NOT ESMF_LIBS_ENV AND NOT ESMF_F90ESMFLINKLIBS_ENV AND NOT ESMF_INCLUDE_OPTS_ENV)
    message(WARNING "Key ESMF environment variables (ESMF_LIBS, ESMF_F90ESMFLINKLIBS, ESMF_INCLUDE_OPTS) are missing. ESMF configuration may be incomplete.")
  endif()

endif()

# Check CMAKE_Fortran_COMPILER, as ESMF usually expects an MPI-aware Fortran compiler.
# Users should typically set CMAKE_Fortran_COMPILER to the ESMF compiler wrapper (e.g., esmpifort).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  message(WARNING "CMAKE_Fortran_COMPILER is a standard compiler (${CMAKE_Fortran_COMPILER}). ESMF typically requires an MPI wrapper (e.g., esmpifort, mpif90). Ensure this compiler is MPI-aware or set CMAKE_Fortran_COMPILER to the appropriate ESMF wrapper.")
else()
  message(STATUS "CMAKE_Fortran_COMPILER is ${CMAKE_Fortran_COMPILER}. Assuming it is an MPI wrapper suitable for ESMF.")
endif()

# PDAF Dependency
set(PDAF_LIB_DIR_ENV $ENV{PDAF_LIB_DIR})
set(PDAF_INC_DIR_ENV $ENV{PDAF_INC_DIR}) # Often PDAF_LIB_DIR/../include or similar

option(USE_PDAF "Enable PDAF features" OFF)

if(PDAF_LIB_DIR_ENV)
  # If PDAF_LIB_DIR is set, we likely want to use PDAF by default,
  # but allow user to override with -DUSE_PDAF=OFF
  if(NOT DEFINED CMAKE_USER_SET_USE_PDAF) # Check if user explicitly set USE_PDAF via -D
    set(USE_PDAF ON CACHE BOOL "Enable PDAF features (auto-enabled due to PDAF_LIB_DIR)" FORCE)
  endif()
  message(STATUS "PDAF_LIB_DIR environment variable is set to: ${PDAF_LIB_DIR_ENV}")
else()
  if(USE_PDAF)
    message(WARNING "USE_PDAF is ON but PDAF_LIB_DIR environment variable is not set. PDAF configuration might fail.")
  endif()
endif()

if(USE_PDAF)
  if(NOT PDAF_LIB_DIR_ENV)
    message(FATAL_ERROR "USE_PDAF is ON but PDAF_LIB_DIR environment variable is not set. Please set it to point to your PDAF library directory.")
  endif()

  set(PDAF_INCLUDE_DIR ${PDAF_LIB_DIR_ENV}/../include) # Common pattern for PDAF
  if(PDAF_INC_DIR_ENV)
    set(PDAF_INCLUDE_DIR ${PDAF_INC_DIR_ENV}) # User-specified include dir takes precedence
  endif()

  if(EXISTS "${PDAF_INCLUDE_DIR}")
    include_directories(${PDAF_INCLUDE_DIR})
    message(STATUS "PDAF include directory added: ${PDAF_INCLUDE_DIR}")
  else()
    message(WARNING "PDAF include directory not found: ${PDAF_INCLUDE_DIR}. This might cause compilation errors for PDAF features.")
  endif()

  # Find PDAF library (typically libpdaf-d)
  # The Makefile uses -lpdaf-d, suggesting dynamic or static linking.
  # We look for both .a and .so/.dylib
  find_library(PDAF_LIBRARY NAMES pdaf-d NAMES_PER_DIR PATHS "${PDAF_LIB_DIR_ENV}" NO_DEFAULT_PATH)

  if(PDAF_LIBRARY)
    message(STATUS "PDAF library found: ${PDAF_LIBRARY}")
    set(PDAF_LIBRARIES ${PDAF_LIBRARY})
  else()
    message(FATAL_ERROR "PDAF library 'pdaf-d' not found in ${PDAF_LIB_DIR_ENV}. Please ensure it exists (e.g., libpdaf-d.a or libpdaf-d.so).")
  endif()

  add_compile_definitions(USE_PDAF)
  message(STATUS "PDAF is enabled. Compile definition USE_PDAF added.")
else()
  message(STATUS "PDAF is disabled.")
endif()

# Store CMAKE_USER_SET_USE_PDAF to respect user's -D option in subsequent runs
if(DEFINED USE_PDAF_INIT)
    # USE_PDAF_INIT is a temporary variable to capture the initial state provided by the user or default
    if(NOT USE_PDAF STREQUAL USE_PDAF_INIT)
        set(CMAKE_USER_SET_USE_PDAF ON CACHE INTERNAL "User has explicitly set USE_PDAF")
    endif()
else()
    # Initialize USE_PDAF_INIT. This only happens on the very first configure.
    set(USE_PDAF_INIT ${USE_PDAF} CACHE INTERNAL "Initial value of USE_PDAF for change detection")
    if(USE_PDAF AND PDAF_LIB_DIR_ENV) # If default was ON due to ENV var
        # Don't mark as user-set unless they explicitly set it
    elseif(DEFINED CMAKE_CACHE_ARGS_USE_PDAF) # A way to check if -D was used, though not standard
        set(CMAKE_USER_SET_USE_PDAF ON CACHE INTERNAL "User has explicitly set USE_PDAF")
    endif()
endif()

# FABM Dependency
# Check for FABM libraries within the SCHISM_LIBRARY_DIR
set(FABM_core_lib_path "${SCHISM_LIBRARY_DIR}/libfabm.a")
set(FABM_schism_lib_path "${SCHISM_LIBRARY_DIR}/libfabm_schism.a")

option(USE_FABM "Enable FABM features" OFF)

if(EXISTS ${FABM_core_lib_path})
  if(NOT DEFINED CMAKE_USER_SET_USE_FABM) # Check if user explicitly set USE_FABM via -D
    set(USE_FABM ON CACHE BOOL "Enable FABM features (auto-enabled due to found libfabm.a)" FORCE)
    message(STATUS "Found libfabm.a in SCHISM_LIBRARY_DIR. USE_FABM automatically enabled.")
  elseif(USE_FABM)
    message(STATUS "Found libfabm.a in SCHISM_LIBRARY_DIR. USE_FABM is ON as per user request or prior auto-enable.")
  else()
    message(STATUS "Found libfabm.a in SCHISM_LIBRARY_DIR, but USE_FABM is OFF as per user request.")
  endif()
else()
  if(USE_FABM)
    message(WARNING "USE_FABM is ON, but libfabm.a not found in ${SCHISM_LIBRARY_DIR}. FABM features might not compile/link correctly.")
  else()
    message(STATUS "libfabm.a not found in ${SCHISM_LIBRARY_DIR}. USE_FABM is OFF.")
  endif()
endif()

if(USE_FABM)
  add_compile_definitions(USE_FABM)
  message(STATUS "FABM is enabled. Compile definition USE_FABM added.")

  # Define FABM libraries as imported targets
  if(EXISTS ${FABM_core_lib_path})
    add_library(FABM_core STATIC IMPORTED GLOBAL)
    set_target_properties(FABM_core PROPERTIES IMPORTED_LOCATION "${FABM_core_lib_path}")
    list(APPEND SCHISM_LIBRARIES FABM_core) # Adding to SCHISM_LIBRARIES for now, or create a new list
    message(STATUS "FABM core library (libfabm.a) found and added as target.")
  else()
    message(FATAL_ERROR "USE_FABM is ON but FABM core library (libfabm.a) not found at ${FABM_core_lib_path}.")
  endif()

  if(EXISTS ${FABM_schism_lib_path})
    add_library(FABM_schism STATIC IMPORTED GLOBAL)
    set_target_properties(FABM_schism PROPERTIES IMPORTED_LOCATION "${FABM_schism_lib_path}")
    list(APPEND SCHISM_LIBRARIES FABM_schism) # Adding to SCHISM_LIBRARIES
    message(STATUS "FABM SCHISM library (libfabm_schism.a) found and added as target.")
  else()
    # This might be a critical error or a warning depending on how FABM is structured.
    # The Makefile implies both are needed if USE_FABM is active.
    message(FATAL_ERROR "USE_FABM is ON but FABM SCHISM library (libfabm_schism.a) not found at ${FABM_schism_lib_path}.")
  endif()

  set(FABM_LIBRARIES FABM_core FABM_schism) # Specific list for FABM libs
else()
  message(STATUS "FABM is disabled.")
endif()

# Logic to remember if user explicitly set USE_FABM
if(DEFINED USE_FABM_INIT)
    if(NOT USE_FABM STREQUAL USE_FABM_INIT)
        set(CMAKE_USER_SET_USE_FABM ON CACHE INTERNAL "User has explicitly set USE_FABM")
    endif()
else()
    set(USE_FABM_INIT ${USE_FABM} CACHE INTERNAL "Initial value of USE_FABM for change detection")
    if(USE_FABM AND EXISTS ${FABM_core_lib_path}) # If default was ON due to libfabm.a
        # Don't mark as user-set unless they explicitly set it
    elseif(DEFINED CMAKE_CACHE_ARGS_USE_FABM)
        set(CMAKE_USER_SET_USE_FABM ON CACHE INTERNAL "User has explicitly set USE_FABM")
    endif()
endif()

# Add src subdirectory
add_subdirectory(src)
