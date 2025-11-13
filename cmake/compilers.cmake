# This code is part of the SCHISM-ESMF interface.
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

# ESMF is expected to be found by the main CMakeLists.txt before this script is included.
# The necessary ESMF variables (e.g., ESMF_F90COMPILER) should already be in global scope.

set(CMAKE_Fortran_COMPILER ${ESMF_F90COMPILER}  CACHE FILEPATH "ESMF Fortran compiler" FORCE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEOPTS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ESMF_F90LINKOPTS}")
message(STATUS "Using ESMF Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Using ESMF Fortran compile options: ${ESMF_F90COMPILEOPTS}")
message(STATUS "Using ESMF Fortran linker options: ${ESMF_F90LINKOPTS}")

# Try --show first
execute_process(
  COMMAND ${CMAKE_Fortran_COMPILER} --show
  RESULT_VARIABLE _result
  OUTPUT_VARIABLE _compiler_cmd
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(_result EQUAL 0)
  set(MPI_SHOW_CMD "--show")
else()
  # Fallback to -show
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -show
    OUTPUT_VARIABLE _compiler_cmd
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(MPI_SHOW_CMD "-show")
endif()

message(STATUS "Using MPI show command: ${MPI_SHOW_CMD}")
string(REGEX MATCH "^[^ ]+" _underlying_compiler "${_compiler_cmd}")
message(STATUS "Detected underlying compiler: ${_underlying_compiler}")

set(CMAKE_Fortran_COMPILER ${_underlying_compiler}  CACHE FILEPATH "ESMF Fortran compiler" FORCE)
message(STATUS "Underlying compiler command: ${_underlying_compiler}")

if(_underlying_compiler MATCHES "gfortran$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

set(CMAKE_Fortran_STANDARD 2008)
enable_language(Fortran)

# -------------- C++ section --------------------
# if(NOT CMAKE_CXX_COMPILER)
#     set(CMAKE_CXX_COMPILER ${ESMF_CXXCOMPILER}  CACHE FILEPATH "Custom C++ compiler" FORCE)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ESMF_CXXCOMPILEOPTS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ESMF_CXXLINKOPTS}")
#     message(STATUS "Using ESMF C++ compiler: ${CMAKE_CXX_COMPILER}")
#     message(STATUS "Using ESMF C++ compile options: ${ESMF_CXXCOMPILEOPTS}")
#     message(STATUS "Using ESMF C++ linker options: ${ESMF_CXXLINKOPTS}")
# endif()

# execute_process(
#   COMMAND ${CMAKE_CXX_COMPILER} ${MPI_SHOW_CMD}
#   OUTPUT_VARIABLE _cxx_cmd
#   ERROR_QUIET
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# string(REGEX MATCH "^[^ ]+" _underlying_cxx_compiler "${_cxx_cmd}")
# set(CMAKE_CXX_COMPILER ${_underlying_cxx_compiler}  CACHE FILEPATH "Custom C++ compiler" FORCE)
# message(STATUS "Underlying C++ compiler: ${CMAKE_CXX_COMPILER}")

# set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)

# set(CMAKE_CXX_STANDARD 17)
# enable_language(CXX)

# -------------- C section --------------------

# if(NOT CMAKE_C_COMPILER)
#   get_filename_component(CXX_PATH "${CMAKE_CXX_COMPILER}" PROGRAM REALPATH)
#   get_filename_component(CXX_BIN "${CXX_PATH}" DIRECTORY)
#   get_filename_component(CXX_NAME "${CXX_PATH}" NAME_WE)

#   if(CXX_NAME MATCHES "clang\\+\\+$")
#     string(REPLACE "++" "" C_NAME "${CXX_NAME}")
#   elseif(CXX_NAME MATCHES "g\\+\\+$")
#     string(REPLACE "++" "" C_NAME "${CXX_NAME}")
#   elseif(CXX_NAME MATCHES "icpc$")
#     set(C_NAME "icc")
#   elseif(CXX_NAME MATCHES "mpiicpc$")
#     set(C_NAME "mpiicc")
#   elseif(CXX_NAME MATCHES "mpicxx$" OR CXX_NAME MATCHES "mpiCC$")
#     set(C_NAME "mpicc")
#   else()
#     message(FATAL_ERROR "Unknown C++ compiler: ${CXX_NAME}")
#   endif()

#   find_program(DERIVED_C_COMPILER NAMES ${C_NAME} PATHS ${CXX_BIN} NO_DEFAULT_PATH)
#   if(NOT DERIVED_C_COMPILER)
#     message(FATAL_ERROR "Could not find matching C compiler for ${CXX_NAME} ${CXX_PATH}")
#   endif()

#   message(STATUS "Derived C compiler: ${DERIVED_C_COMPILER}")
#   set(CMAKE_C_COMPILER ${DERIVED_C_COMPILER} CACHE FILEPATH "Derived C compiler" FORCE)
# endif()



