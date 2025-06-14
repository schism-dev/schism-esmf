# src/schism/CMakeLists.txt
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules CACHE PATH "Directory for Fortran modules" FORCE)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Common source files for both ESMF and NUOPC interface variants
set(SCHISM_INTERFACE_COMMON_SOURCES
    schism_bmi.F90
    schism_esmf_util.F90
)

# Source files specific to the ESMF interface variant
set(SCHISM_INTERFACE_ESMF_SOURCES
    schism_esmf_cap.F90
)

# Source files specific to the NUOPC interface variant
set(SCHISM_INTERFACE_NUOPC_SOURCES
    schism_nuopc_util.F90
    schism_nuopc_cap.F90
)

# --- schism_esmf_interface library ---
add_library(schism_esmf_interface STATIC
    ${SCHISM_INTERFACE_COMMON_SOURCES}
    ${SCHISM_INTERFACE_ESMF_SOURCES}
)

target_include_directories(schism_esmf_interface PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY} # For its own modules
    ${SCHISM_INCLUDE_DIR}             # For SCHISM base modules
)
if(ESMF_FOUND AND TARGET ESMF::ESMF)
    target_link_libraries(schism_esmf_interface PUBLIC ESMF::ESMF) # ESMF modules might be needed
elseif(ESMF_LIBRARIES)
    target_link_libraries(schism_esmf_interface PUBLIC ${ESMF_LIBRARIES}) # Fallback
endif()

# Link common dependencies (SCHISM core, PDAF if enabled)
target_link_libraries(schism_esmf_interface PUBLIC ${SCHISM_LIBRARIES})
if(USE_PDAF)
    # schism_esmf_interface depends on pdaf_schism_bindings, which links to actual PDAF libs
    target_link_libraries(schism_esmf_interface PUBLIC pdaf_schism_bindings)
    target_compile_definitions(schism_esmf_interface PRIVATE USE_PDAF)
endif()

# --- schism_nuopc_interface library ---
add_library(schism_nuopc_interface STATIC
    ${SCHISM_INTERFACE_COMMON_SOURCES}
    ${SCHISM_INTERFACE_NUOPC_SOURCES}
)

target_include_directories(schism_nuopc_interface PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${SCHISM_INCLUDE_DIR}
)
if(ESMF_FOUND AND TARGET ESMF::ESMF) # NUOPC also uses ESMF
    target_link_libraries(schism_nuopc_interface PUBLIC ESMF::ESMF)
elseif(ESMF_LIBRARIES)
    target_link_libraries(schism_nuopc_interface PUBLIC ${ESMF_LIBRARIES})
endif()

target_link_libraries(schism_nuopc_interface PUBLIC ${SCHISM_LIBRARIES})
if(USE_PDAF)
    target_link_libraries(schism_nuopc_interface PUBLIC pdaf_schism_bindings)
    target_compile_definitions(schism_nuopc_interface PRIVATE USE_PDAF)
endif()


# Add these libraries as dependencies for the main executables defined in src/CMakeLists.txt
# This creates a global property that the parent CMakeLists.txt can read.
# A more modern approach is to link them directly in src/CMakeLists.txt
# by making these library targets have ALIAS global scope or just linking directly.
# For now, we assume direct linking in parent scope after these targets are defined.

# Modify the main executables in the parent scope (src/CMakeLists.txt) to link these.
# This can be done by getting the targets and adding link libraries.
get_target_property(MAIN_ESMF_TARGET main_esmf TYPE)
if(MAIN_ESMF_TARGET)
    target_link_libraries(main_esmf PRIVATE schism_esmf_interface)
else()
    message(WARNING "main_esmf target not found in parent scope when linking from src/schism/CMakeLists.txt")
endif()

get_target_property(MAIN_NUOPC_TARGET main_nuopc TYPE)
if(MAIN_NUOPC_TARGET)
    target_link_libraries(main_nuopc PRIVATE schism_nuopc_interface)
else()
    message(WARNING "main_nuopc target not found in parent scope when linking from src/schism/CMakeLists.txt")
endif()

# Handle schism_nuopc_cap.mk.in -> schism.mk generation
# The root Makefile has:
# sed 's#@@SCHISM_BUILD_DIR@@#'$(SCHISM_BUILD_DIR)'#g' ./src/schism/schism_nuopc_cap.mk.in > $(DESTDIR)/schism.mk
set(SCHISM_MK_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/schism_nuopc_cap.mk.in")
set(SCHISM_MK_OUTPUT "${CMAKE_BINARY_DIR}/schism.mk") # Staging it in binary dir first

if(EXISTS ${SCHISM_MK_TEMPLATE})
    message(STATUS "Configuring ${SCHISM_MK_TEMPLATE} -> ${SCHISM_MK_OUTPUT}")
    configure_file(${SCHISM_MK_TEMPLATE} ${SCHISM_MK_OUTPUT} @ONLY)
    # Installation of this file will be handled in the root CMakeLists.txt's install rules.
else()
    message(WARNING "Template file schism_nuopc_cap.mk.in not found at ${SCHISM_MK_TEMPLATE}")
endif()
