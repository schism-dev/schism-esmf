# This script handles the setup for the ESMF (Earth System Modeling Framework) dependency.
# It finds the ESMF package, applies necessary include directories, compile options,
# link flags, and checks the Fortran compiler compatibility.

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

# Check CMAKE_Fortran_COMPILER, as ESMF usually expects an MPI-aware Fortran compiler.
# Users should typically set CMAKE_Fortran_COMPILER to the ESMF compiler wrapper (e.g., esmpifort, mpifort, mpif90).
# Check if the compiler name contains common MPI wrapper patterns.
get_filename_component(COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)
if(COMPILER_NAME MATCHES "^(mpi|esmpi|sxmpi)")
  message(STATUS "CMAKE_Fortran_COMPILER is ${CMAKE_Fortran_COMPILER}. Detected as an MPI wrapper.")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  message(WARNING "CMAKE_Fortran_COMPILER is a standard compiler (${CMAKE_Fortran_COMPILER}). ESMF typically requires an MPI wrapper (e.g., esmpifort, mpifort, mpif90). Ensure this compiler is MPI-aware or set CMAKE_Fortran_COMPILER to the appropriate ESMF wrapper.")
else()
  message(STATUS "CMAKE_Fortran_COMPILER is ${CMAKE_Fortran_COMPILER}. Assuming it is an MPI wrapper suitable for ESMF.")
endif()
