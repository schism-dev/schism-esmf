# src/CMakeLists.txt

# Add subdirectories. Order might matter for dependencies between them.
# Based on src/Makefile, typical order might be:
# (PDAF_bindings if USE_PDAF)
# schism (provides modules needed by model and driver)
# model (provides modules needed by driver)
# driver
# (mossco if MOSSCO_DIR is set - will handle MOSSCO later if needed)

if(USE_PDAF)
  add_subdirectory(PDAF_bindings)
endif()

add_subdirectory(schism)

# Ensure executables can find modules from schism interface libraries
target_include_directories(main_esmf PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(main_nuopc PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

# The actual linking of schism_esmf_interface and schism_nuopc_interface
# to main_esmf and main_nuopc respectively, is now handled inside src/schism/CMakeLists.txt
# using target_link_libraries on the main_esmf and main_nuopc targets directly.
# This is okay as long as those targets are already known (which they are, from src/CMakeLists.txt).
add_subdirectory(model)
add_subdirectory(driver)

# Compilation of main_esmf.F90 and main_nuopc.F90
# These are at the src/ level according to src/Makefile
# Their dependencies are objects from driver, model, schism.
# For now, we'll declare them and link them later once the library targets from subdirs are defined.

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules) # Central place for .mod files
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# ESMF executable
add_executable(main_esmf main_esmf.F90)
target_include_directories(main_esmf PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} # For .mod files from this project
    # Other necessary include paths from dependencies like ESMF, SCHISM etc.
    # These are often handled by linking to imported targets or via include_directories in parent scope.
)
# Target linking will be specified in subdirectories or after all targets are known.
# For now, ensure it links with ESMF libraries if found.
if(ESMF_LIBRARIES)
    target_link_libraries(main_esmf PRIVATE ${ESMF_LIBRARIES})
endif()
if(MPI_Fortran_FOUND) # ESMF usually needs MPI
    target_link_libraries(main_esmf PRIVATE MPI::MPI_Fortran)
endif()
 # Link SCHISM libraries
if(SCHISM_LIBRARIES)
    target_link_libraries(main_esmf PRIVATE ${SCHISM_LIBRARIES})
endif()
# Link core dependencies (NetCDF, LAPACK, BLAS)
if(NetCDF_FOUND)
    target_link_libraries(main_esmf PRIVATE NetCDF::NetCDF_Fortran) # Or however NetCDF exports its Fortran libs
endif()
if(LAPACK_FOUND)
    target_link_libraries(main_esmf PRIVATE ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
    target_link_libraries(main_esmf PRIVATE ${BLAS_LIBRARIES})
endif()


# NUOPC executable
add_executable(main_nuopc main_nuopc.F90)
target_include_directories(main_nuopc PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)
if(ESMF_LIBRARIES) # NUOPC also builds on ESMF
    target_link_libraries(main_nuopc PRIVATE ${ESMF_LIBRARIES})
endif()
if(MPI_Fortran_FOUND)
    target_link_libraries(main_nuopc PRIVATE MPI::MPI_Fortran)
endif()
if(SCHISM_LIBRARIES)
    target_link_libraries(main_nuopc PRIVATE ${SCHISM_LIBRARIES})
endif()
if(NetCDF_FOUND)
    target_link_libraries(main_nuopc PRIVATE NetCDF::NetCDF_Fortran)
endif()
if(LAPACK_FOUND)
    target_link_libraries(main_nuopc PRIVATE ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
    target_link_libraries(main_nuopc PRIVATE ${BLAS_LIBRARIES})
endif()


# The MOSSCO part from src/Makefile:
# MOSSCO_LIB=
# ifeq ($(origin MOSSCO_DIR), environment)
#  SUBDIRS+=mossco
#  MOSSCO_LIB+=-L mossco -lmossco
# endif
# This implies a mossco subdirectory and a libmossco.a library from it.
# We will handle this when/if we convert the mossco subdirectory.
set(MOSSCO_DIR_ENV $ENV{MOSSCO_DIR})
if(MOSSCO_DIR_ENV)
    message(STATUS "MOSSCO_DIR environment variable found: ${MOSSCO_DIR_ENV}")
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mossco")
        message(STATUS "mossco subdirectory exists. Will add it.")
        add_subdirectory(mossco)
        # Link main_esmf with libmossco if it's built by mossco/CMakeLists.txt
        # Assuming mossco will create a target, e.g., mossco_lib
        # if(TARGET mossco_lib)
        #   target_link_libraries(main_esmf PRIVATE mossco_lib)
        # endif()
    else()
        message(WARNING "MOSSCO_DIR is set, but src/mossco subdirectory does not exist. Skipping.")
    endif()
else()
    message(STATUS "MOSSCO_DIR environment variable not set. Skipping mossco subdirectory.")
endif()
