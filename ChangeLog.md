# ChangeLog

All notable changes to schism-esmf will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - v0.3.0

### Added
- Comprehensive documentation infrastructure with MkDocs and Material theme
- FORD-based Fortran API documentation generation
- Consolidated documentation configuration in `doc/config/`
- CMake DOCS_ONLY build mode for documentation-only builds
- GitHub Actions workflow for automated documentation deployment
- GitLab CI pipeline for documentation publishing
- Read the Docs integration
- Build troubleshooting guides and architecture documentation
- QuickStart guide for new users

### Changed
- Reorganized documentation configuration files to `doc/config/` for consistency
- Updated all CI/CD references to new config paths
- Simplified MkDocs emoji extension configuration

## [0.2.0] - 2025-11-13

### Added
- Complete CMake build system conversion (Makefiles retained for compatibility)
- CMake discovery for ESMF via `FindESMF.cmake` module
- Automated CI/CD workflows for build testing
- Documentation integration in CMake build system
- CMake subdirectory structure for modular builds
- Support for DOCS_ONLY build mode

### Changed
- Primary build system migrated from Make to CMake
- Build configuration standardized across platforms
- ESMF detection improved with esmf.mk parsing

### Fixed
- Build issues with ESMF library discovery
- Module file generation and dependency tracking

## [0.1.0] - 2020-02-12

### Added
- Initial ESMF/NUOPC cap implementation for SCHISM
- Basic Makefile-based build system
- ESMF cap with parallel communication support
- NUOPC cap for standardized model coupling
- Example configurations (83, CFL, QuarterAnnulus)
- Script utilities for setup creation
- Support for FABM biogeochemical coupling
- Mesh creation utilities for ESMF grids
- Field export/import for ocean coupling (SST, SSH, currents, salinity)
- Integration with SCHISM hydrodynamics core

### Changed
- Forked from SCHISM SVN repository revision 5235
- Adapted for git-based development workflow

## [0.0.1] - 2019-09-11

### Added
- Initial commit and project setup
- Basic repository structure
- License files (GPL-3.0, Apache-2.0, CC-BY-SA-4.0, CC0-1.0)

---

## Version History Notes

- **v0.0.1**: Initial repository creation from SVN fork
- **v0.1.0**: First functional ESMF/NUOPC implementation with Makefile build
- **v0.2.0**: Major refactoring with CMake build system and CI/CD
- **v0.3.0**: Documentation infrastructure and enhanced developer experience (unreleased)

## Notable Intermediate Features (between v0.1.0 and v0.2.0)

- 2024-03-02: NUOPC cap restructuring (PR #30)
- 2024-02-17: NUOPC cap attribute fixes (PR #29)
- 2024-01-26: Coastal application support (PR #27)
- 2023-12-04: Integration with HiFIS framework
- 2023-04-13: Coastal app enhancements (PR #24)
- 2022-08-31: UGRID support (PR #76)
- 2022-01-29: Coupling model improvements (PR #21)

## Contributors

Major contributors (commits):
- Joseph Zhang (907 commits)
- Carsten Lemmen (490 commits)
- Zhengui Wang (155 commits)
- Fei Ye (122 commits)
- Nicole Cai (79 commits)
- Dan Yu (72 commits)
- Wei Huang (39 commits)
- And many others (see git shortlog for complete list)
