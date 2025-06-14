# src/driver/CMakeLists.txt
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules CACHE PATH "Directory for Fortran modules" FORCE)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# List all .F90 source files in the current directory
file(GLOB DRIVER_SOURCES *.F90)

if(NOT DRIVER_SOURCES)
    message(WARNING "No .F90 files found in src/driver. The schism_driver_libs library will be empty.")
else()
    message(STATUS "Driver sources: ${DRIVER_SOURCES}")
endif()

add_library(schism_driver_libs STATIC ${DRIVER_SOURCES})

install(TARGETS schism_driver_libs
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/modules/driver"
)

# Add include directories required by these driver components
target_include_directories(schism_driver_libs PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} # For its own modules and those from other components
    ${SCHISM_INCLUDE_DIR}             # For SCHISM base modules
    # Model and Schism interface modules should be found via CMAKE_Fortran_MODULE_DIRECTORY
)

# Link against libraries from schism and model subdirectories
# These are dependencies for compiling the driver sources
if(TARGET schism_esmf_interface)
    target_link_libraries(schism_driver_libs PUBLIC schism_esmf_interface)
else()
    message(WARNING "Target schism_esmf_interface not found. Linking schism_driver_libs might fail or be incomplete.")
endif()

if(TARGET schism_nuopc_interface)
    target_link_libraries(schism_driver_libs PUBLIC schism_nuopc_interface)
else()
    message(WARNING "Target schism_nuopc_interface not found. Linking schism_driver_libs might fail or be incomplete.")
endif()

if(TARGET schism_model_libs)
    target_link_libraries(schism_driver_libs PUBLIC schism_model_libs)
else()
    message(WARNING "Target schism_model_libs not found. Linking schism_driver_libs might fail or be incomplete.")
endif()

# Link against ESMF as driver components use ESMF functionalities
if(ESMF_FOUND AND TARGET ESMF::ESMF)
    target_link_libraries(schism_driver_libs PUBLIC ESMF::ESMF)
elseif(ESMF_LIBRARIES) # Fallback
    target_link_libraries(schism_driver_libs PUBLIC ${ESMF_LIBRARIES})
endif()

# Note: MOSSCO dependency for toplevel_schism_netcdf.F90 will be handled
# by linking the main_esmf executable to the mossco library if MOSSCO is enabled.
# So, schism_driver_libs itself doesn't need to link mossco_lib directly here.
