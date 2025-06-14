add_executable(main_esmf main_esmf.F90)
add_executable(main_nuopc main_nuopc.F90)

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
 # Link SCHISM libraries
if(SCHISM_LIBRARIES)
    target_link_libraries(main_esmf PRIVATE ${SCHISM_LIBRARIES})
endif()


# NUOPC executable
target_include_directories(main_nuopc PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)
if(ESMF_LIBRARIES) # NUOPC also builds on ESMF
    target_link_libraries(main_nuopc PRIVATE ${ESMF_LIBRARIES})
endif()
if(SCHISM_LIBRARIES)
    target_link_libraries(main_nuopc PRIVATE ${SCHISM_LIBRARIES})
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
