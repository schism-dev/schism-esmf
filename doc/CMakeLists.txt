# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>

# Documentation build configuration for SCHISM-ESMF
# This CMakeLists.txt provides targets to build HTML documentation using MkDocs

# Check if BUILD_DOCS is enabled
if(NOT BUILD_DOCS)
  return()
endif()

message(STATUS "")
message(STATUS "=== Documentation Build Configuration ===")

# Find Python3
find_package(Python3 COMPONENTS Interpreter)

if(NOT Python3_FOUND)
  message(WARNING "Python3 not found. Documentation build targets will not be available.")
  message(WARNING "Install Python 3.8+ to enable documentation building.")
  return()
endif()

message(STATUS "✓ Python found: ${Python3_EXECUTABLE} (version ${Python3_VERSION})")

# Check for mkdocs and mkdocs-material
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import mkdocs; print(mkdocs.__version__)"
  RESULT_VARIABLE MKDOCS_CHECK
  OUTPUT_VARIABLE MKDOCS_VERSION
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import material; print(material.__version__)"
  RESULT_VARIABLE MKDOCS_MATERIAL_CHECK
  OUTPUT_VARIABLE MKDOCS_MATERIAL_VERSION
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(MKDOCS_AVAILABLE TRUE)

if(MKDOCS_CHECK EQUAL 0)
  message(STATUS "✓ MkDocs found: version ${MKDOCS_VERSION}")
else()
  message(WARNING "⚠ MkDocs not found. Install with: pip install mkdocs")
  set(MKDOCS_AVAILABLE FALSE)
endif()

if(MKDOCS_MATERIAL_CHECK EQUAL 0)
  message(STATUS "✓ MkDocs Material theme found: version ${MKDOCS_MATERIAL_VERSION}")
else()
  message(WARNING "⚠ MkDocs Material theme not found. Install with: pip install mkdocs-material")
  set(MKDOCS_AVAILABLE FALSE)
endif()

if(NOT MKDOCS_AVAILABLE)
  message(WARNING "")
  message(WARNING "Documentation dependencies missing. Install with:")
  message(WARNING "  pip install mkdocs mkdocs-material")
  message(WARNING "  OR")
  message(WARNING "  conda install -c conda-forge mkdocs mkdocs-material")
  message(WARNING "")
  return()
endif()

# Ensure install directory variables are available even in DOCS_ONLY mode
include(GNUInstallDirs)
if(NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "share/doc/${CMAKE_PROJECT_NAME}" CACHE PATH "" FORCE)
endif()

# Set documentation paths
set(DOCS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(DOCS_BUILD_DIR "${CMAKE_BINARY_DIR}/docs")
set(MKDOCS_CONFIG "${CMAKE_SOURCE_DIR}/mkdocs.yml")
set(FORD_CONFIG "${CMAKE_SOURCE_DIR}/ford.yml")
set(FORD_BUILD_DIR "${CMAKE_BINARY_DIR}/ford_docs")

# Create build directory for docs
file(MAKE_DIRECTORY "${DOCS_BUILD_DIR}")

# Custom target: Build HTML documentation
add_custom_target(docs
  COMMAND ${Python3_EXECUTABLE} -m mkdocs build
    --config-file "${MKDOCS_CONFIG}"
    --site-dir "${DOCS_BUILD_DIR}"
    --strict
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  COMMENT "Building HTML documentation with MkDocs..."
  VERBATIM
)

# Custom target: Build API reference with FORD
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import ford; print(ford.__version__)"
  RESULT_VARIABLE FORD_CHECK
  OUTPUT_VARIABLE FORD_VERSION
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(FORD_CHECK EQUAL 0)
  message(STATUS "✓ FORD found: version ${FORD_VERSION}")
  
  add_custom_target(docs-api
    COMMAND ${Python3_EXECUTABLE} -m ford
      "${FORD_CONFIG}"
      -o "${FORD_BUILD_DIR}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Building API reference with FORD..."
    VERBATIM
  )
  
  add_custom_target(docs-all
    DEPENDS docs docs-api
    COMMENT "Building all documentation (user guide + API reference)"
  )
  
  message(STATUS "✓ API documentation targets available:")
  message(STATUS "    make docs-api   - Build Fortran API reference")
  message(STATUS "    make docs-all   - Build user guide + API reference")
else()
  message(STATUS "○ FORD not found (API reference generation not available)")
  message(STATUS "    Install with: pip install ford")
endif()

# Custom target: Serve documentation locally for development
add_custom_target(docs-serve
  COMMAND ${Python3_EXECUTABLE} -m mkdocs serve
    --config-file "${MKDOCS_CONFIG}"
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  COMMENT "Serving documentation at http://127.0.0.1:8000 (Press Ctrl+C to stop)"
  VERBATIM
)

# Custom target: Clean documentation build artifacts
add_custom_target(docs-clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${DOCS_BUILD_DIR}"
  COMMENT "Cleaning documentation build directory"
  VERBATIM
)

# Installation: Install built HTML documentation
install(
  DIRECTORY "${DOCS_BUILD_DIR}/"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}/html"
  OPTIONAL
  MESSAGE_NEVER
)

# Also install source markdown files for reference
install(
  DIRECTORY "${DOCS_SOURCE_DIR}/"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}/markdown"
  FILES_MATCHING PATTERN "*.md"
  PATTERN "_config.yml" EXCLUDE
)

# Install top-level documentation files
install(
  FILES
    "${CMAKE_SOURCE_DIR}/Readme.md"
    "${CMAKE_SOURCE_DIR}/QUICKSTART.md"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)

message(STATUS "")
message(STATUS "Documentation targets available:")
message(STATUS "  make docs        - Build HTML documentation")
message(STATUS "  make docs-serve  - Serve docs locally at http://127.0.0.1:8000")
message(STATUS "  make docs-clean  - Clean documentation build artifacts")
message(STATUS "")
message(STATUS "Built documentation will be at: ${DOCS_BUILD_DIR}")
message(STATUS "Install destination: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR}")
message(STATUS "")
