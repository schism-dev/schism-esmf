
set(MODEL_SOURCES
    atmosphere_cmi_esmf.F90
    atmosphere_cmi_nuopc.F90
)

add_library(schism_model_libs STATIC ${MODEL_SOURCES})

install(TARGETS schism_model_libs
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/modules/model"
)

# Add include directories required by these model components
# They likely need ESMF includes and modules from src/schism
target_include_directories(schism_model_libs PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} # For its own modules and those from schism_interface
    ${SCHISM_INCLUDE_DIR}             # For SCHISM base modules
)

# Link against ESMF as these are ESMF components (CMI = Custom Model Implementation)
if(ESMF_FOUND AND TARGET ESMF::ESMF)
    target_link_libraries(schism_model_libs PUBLIC ESMF::ESMF)
elseif(ESMF_LIBRARIES) # Fallback
    target_link_libraries(schism_model_libs PUBLIC ${ESMF_LIBRARIES})
endif()

# Link against schism interface libraries (for modules like schism_esmf_util)
if(TARGET schism_esmf_interface)
    target_link_libraries(schism_model_libs PUBLIC schism_esmf_interface)
endif()

if(TARGET schism_nuopc_interface)
    target_link_libraries(schism_model_libs PUBLIC schism_nuopc_interface)
endif()