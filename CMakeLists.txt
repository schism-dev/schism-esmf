cmake_minimum_required(VERSION 3.16)
project(SCHISM_ESMF_Interface LANGUAGES Fortran C)

enable_language(Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ESMF Dependency
# Use the FindESMF.cmake module to locate ESMF and its settings.
find_package(ESMF MODULE REQUIRED) # Assumes FindESMF.cmake is in CMAKE_MODULE_PATH

# ESMF_FOUND is implicitly true if REQUIRED is used and no error occurred.
message(STATUS "ESMF found via FindESMF.cmake module.")
message(STATUS "ESMF Version: ${ESMF_VERSION}")
message(STATUS "ESMF Include Dirs: ${ESMF_INCLUDE_DIRS}")
message(STATUS "ESMF Libraries: ${ESMF_LIBRARIES}") # This might be the path to libesmf.so/a or an IMPORTED target
message(STATUS "ESMF Compile Definitions: ${ESMF_COMPILE_DEFINITIONS}")
message(STATUS "ESMF Compile Options: ${ESMF_COMPILE_OPTIONS}")
message(STATUS "ESMF Link Flags: ${ESMF_LINK_FLAGS}")

# Apply ESMF settings
if(ESMF_INCLUDE_DIRS)
  include_directories(${ESMF_INCLUDE_DIRS})
endif()

if(ESMF_COMPILE_DEFINITIONS)
  add_compile_definitions(${ESMF_COMPILE_DEFINITIONS})
endif()

if(ESMF_COMPILE_OPTIONS)
  # These are general compile options. Add them globally or per target.
  # For global options:
  add_compile_options(${ESMF_COMPILE_OPTIONS})
  # For target-specific options, use target_compile_options(myTarget PRIVATE ${ESMF_COMPILE_OPTIONS})
  # Global is likely fine here as these usually affect all Fortran code using ESMF.
endif()

# Linking:
# Prefer using the ESMF::ESMF target if it's robustly set up by FindESMF.cmake
# This target should encapsulate include directories, link libraries, and other properties.
# For now, we will explicitly add link flags and assume ESMF_LIBRARIES refers to the actual library files or names.
# If ESMF::ESMF is the intended way, target_link_libraries(your_target PUBLIC ESMF::ESMF) would be used.
# The variable ESMF_LIBRARIES from FindESMF.cmake is currently set to ESMF_LIBRARY_LOCATION.
# The variable ESMF_INTERFACE_LINK_LIBRARIES from FindESMF.cmake has more comprehensive link information.
# Let's rely on ESMF_LINK_FLAGS for additional linker flags, and ESMF_LIBRARIES for the core library.
# The ESMF::ESMF target's INTERFACE_LINK_LIBRARIES is set to ESMF_INTERFACE_LINK_LIBRARIES.
# So, linking against ESMF::ESMF should be preferred for targets.

if(ESMF_LINK_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ESMF_LINK_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ESMF_LINK_FLAGS}")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${ESMF_LINK_FLAGS}")
endif()

# Note: Actual linking of ESMF libraries to your targets (e.g., main_esmf, main_nuopc)
# will happen in their respective target_link_libraries() calls, typically with ESMF::ESMF or ${ESMF_LIBRARIES}.
# For example: target_link_libraries(my_target PRIVATE ESMF::ESMF)


# SCHISM Dependency
# SCHISM is expected to be pre-built. Its location is specified by SCHISM_BUILD_DIR.
# MPI, NetCDF, LAPACK, and BLAS are expected to be handled by the ESMF dependency.
# If ESMF is not found or doesn't correctly export them, issues might arise.
# Users might need to ensure their ESMF build correctly wraps these or install them separately
# and adjust CMake configurations if direct linking becomes necessary.
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

# # Handle FABM dependency
# include(cmake/fabm.cmake)

# Add src subdirectory
add_subdirectory(src)

# add_subdirectory(example)

# install(TARGETS main_esmf main_nuopc
#     RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
# )

# Install the configured schism.mk file
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
