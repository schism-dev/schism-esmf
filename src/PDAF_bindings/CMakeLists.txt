# src/PDAF_bindings/CMakeLists.txt
if(NOT USE_PDAF)
    message(STATUS "USE_PDAF is OFF, skipping PDAF_bindings.")
    return()
endif()

message(STATUS "Processing PDAF_bindings (USE_PDAF is ON)")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules CACHE PATH "Directory for Fortran modules" FORCE)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# List all .F90 source files in the current directory
file(GLOB PDAF_BINDING_SOURCES *.F90)

if(NOT PDAF_BINDING_SOURCES)
    message(WARNING "No .F90 files found in src/PDAF_bindings. The pdaf_schism_bindings library will be empty.")
else()
    message(STATUS "PDAF Binding sources: ${PDAF_BINDING_SOURCES}")
endif()

add_library(pdaf_schism_bindings STATIC ${PDAF_BINDING_SOURCES})

target_include_directories(pdaf_schism_bindings PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} # For its own modules
    ${SCHISM_INCLUDE_DIR}             # For SCHISM base modules
    ${PDAF_INCLUDE_DIR}               # For PDAF modules
)

# Link against the main PDAF library
if(PDAF_LIBRARIES)
    target_link_libraries(pdaf_schism_bindings PUBLIC ${PDAF_LIBRARIES})
else()
    message(WARNING "PDAF_LIBRARIES not defined. Linking pdaf_schism_bindings might fail.")
endif()

# Link against SCHISM libraries if there are direct dependencies
# This is often needed for utility functions or types from SCHISM
if(SCHISM_LIBRARIES)
    target_link_libraries(pdaf_schism_bindings PUBLIC ${SCHISM_LIBRARIES})
endif()

# Add USE_PDAF compile definition, as these bindings are PDAF specific
target_compile_definitions(pdaf_schism_bindings PRIVATE USE_PDAF)
