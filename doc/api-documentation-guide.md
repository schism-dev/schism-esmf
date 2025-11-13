# Documenting SCHISM-ESMF Fortran Code

This guide explains how to document the Fortran code for automatic API reference generation using FORD.

## FORD Documentation Syntax

FORD uses special comment markers to extract documentation:

- `!` or `!!` - Documentation comments (FORD will extract these)
- `!>` - Brief description of item that follows
- `!<` - Brief description of item on same line

### Module Documentation

```fortran
!> @brief Basic Model Interface (BMI) for SCHISM
!!
!! This module provides a generic BMI interface to SCHISM that allows
!! the model to be controlled externally through standardized function calls.
!!
!! @author Carsten Lemmen <carsten.lemmen@hereon.de>
!! @author Richard Hofmeister
!! @copyright 2021-2023 Helmholtz-Zentrum Hereon
!! @license Apache-2.0
module schism_bmi

  use esmf
  
  implicit none
  private
  
  public :: schism_parallel_init
  public :: schismTimeStep
  
contains
```

### Subroutine/Function Documentation

```fortran
!> @brief Initialize SCHISM parallel environment
!!
!! Initializes the MPI communicator for SCHISM using the provided
!! communicator from ESMF or another coupling framework.
!!
!! @param[in] communicator MPI communicator to use for SCHISM
!!
!! @note This must be called before any other SCHISM routines
!!
!! @author Carsten Lemmen
!! @date 2021-2023
subroutine schism_parallel_init(communicator)

  use schism_msgp, only: parallel_init
  
  implicit none
  integer, intent(in) :: communicator !< MPI communicator handle

  call parallel_init(communicator)

end subroutine schism_parallel_init
```

### Variable Documentation

```fortran
integer :: nsteps             !< Number of time steps to run
real(8) :: dt                 !< Time step size in seconds
character(len=256) :: indir   !< Input directory path

type :: SchismState
  real(8), pointer :: temp(:)  !< Temperature field [K]
  real(8), pointer :: salt(:)  !< Salinity field [PSU]
  integer :: nNodes            !< Number of grid nodes
end type SchismState
```

### Type Documentation

```fortran
!> @brief SCHISM model state container
!!
!! This derived type holds pointers to SCHISM's internal state variables
!! for access by coupling frameworks.
!!
!! @note All pointer members must be associated before use
type :: SchismState
  real(8), pointer :: temp(:)  !< Temperature field [K]
  real(8), pointer :: salt(:)  !< Salinity field [PSU]
  integer :: nNodes            !< Number of grid nodes
  integer :: nLevels           !< Number of vertical levels
end type SchismState
```

## FORD Tags Reference

Common FORD documentation tags:

| Tag | Purpose | Example |
|-----|---------|---------|
| `@brief` | Short one-line description | `@brief Initialize the model` |
| `@param[in]` | Input parameter | `@param[in] dt Time step size` |
| `@param[out]` | Output parameter | `@param[out] rc Return code` |
| `@param[inout]` | Input/output parameter | `@param[inout] state Model state` |
| `@return` | Function return value | `@return Pointer to array` |
| `@author` | Code author | `@author John Doe` |
| `@date` | Date range | `@date 2021-2023` |
| `@note` | Important note | `@note Thread-safe` |
| `@warning` | Warning message | `@warning Not validated` |
| `@bug` | Known bug | `@bug Fails on empty input` |
| `@todo` | Future work | `@todo Add error handling` |
| `@see` | Cross-reference | `@see schism_init` |

## Example: Fully Documented Module

```fortran
!> @brief SCHISM ESMF coupling utilities
!!
!! This module provides utility functions for coupling SCHISM with ESMF,
!! including grid creation, field registration, and data exchange helpers.
!!
!! The module handles:
!! - Creation of ESMF meshes from SCHISM grids
!! - Field bundle management
!! - Coordinate system transformations
!!
!! @author Carsten Lemmen <carsten.lemmen@hereon.de>
!! @copyright 2021-2025 Helmholtz-Zentrum Hereon
!! @license Apache-2.0
!!
!! @see schism_bmi
!! @see schism_esmf_cap
module schism_esmf_util

  use esmf
  use schism_glbl, only: np, ne
  
  implicit none
  private
  
  !> Maximum string length for field names
  integer, parameter :: MAXSTRLEN = 256
  
  !> Default fill value for missing data
  real(ESMF_KIND_R8), parameter :: MISSING_VALUE = -999.0_ESMF_KIND_R8
  
  ! Public interfaces
  public :: create_schism_mesh
  public :: register_schism_fields
  
contains

  !> @brief Create ESMF mesh from SCHISM grid
  !!
  !! Constructs an ESMF mesh object from SCHISM's internal unstructured grid
  !! representation. The mesh includes both node and element connectivity.
  !!
  !! @param[out] mesh      ESMF mesh object (allocated by this routine)
  !! @param[in]  coordSys  Coordinate system (ESMF_COORDSYS_CART or ESMF_COORDSYS_SPH_DEG)
  !! @param[out] rc        Return code (ESMF_SUCCESS on success)
  !!
  !! @note The mesh must be destroyed with ESMF_MeshDestroy when no longer needed
  !!
  !! @warning Only supports 2D horizontal meshes; vertical structure is handled separately
  !!
  !! @author Carsten Lemmen
  !! @date 2021-2023
  subroutine create_schism_mesh(mesh, coordSys, rc)
  
    type(ESMF_Mesh), intent(out) :: mesh
    type(ESMF_CoordSys_Flag), intent(in), optional :: coordSys
    integer, intent(out), optional :: rc
    
    ! Local variables
    integer :: localrc
    integer :: numNodes    !< Local number of nodes
    integer :: numElems    !< Local number of elements
    
    ! Implementation...
    
    if (present(rc)) rc = ESMF_SUCCESS
    
  end subroutine create_schism_mesh

end module schism_esmf_util
```

## Building API Documentation

### Local Build

```bash
# Install FORD
pip install ford

# Generate API docs
ford ford.yml

# Output in build/ford_docs/index.html
```

### With CMake

```bash
# Configure with docs enabled
cmake .. -DBUILD_DOCS=ON -DDOCS_ONLY=ON

# Build API reference
make docs-api

# Build both user guide and API
make docs-all

# Output:
# - User guide: build/docs/index.html
# - API reference: build/ford_docs/index.html
```

## Integration with MkDocs

You can link from the MkDocs user guide to the FORD API reference:

In `doc/index.md`:
```markdown
## API Reference

For detailed Fortran API documentation, see the [API Reference](../ford_docs/index.html).

Key modules:
- [schism_bmi](../ford_docs/module/schism_bmi.html) - Basic Model Interface
- [schism_esmf_cap](../ford_docs/module/schism_esmf_cap.html) - ESMF Cap
```

## Best Practices

1. **Document every public interface**
   - All public modules, subroutines, functions
   - Include parameter descriptions with units
   - Explain intent (in/out/inout)

2. **Use consistent formatting**
   - Start with `@brief` for one-line summary
   - Add detailed description in body
   - Use `@param` for all arguments
   - Add `@author` and `@date` to modules

3. **Cross-reference related items**
   - Use `@see` to link related routines
   - Reference types, modules in descriptions

4. **Include examples**
   - Show typical usage in module header
   - Document preconditions and postconditions

5. **Mark work-in-progress**
   - Use `@todo` for planned features
   - Use `@bug` for known issues
   - Use `@warning` for limitations

## Documentation Coverage

Check which files need documentation:

```bash
# Find undocumented modules
grep -L "^!>" src/**/*.F90

# Find modules without @brief
grep -L "@brief" src/**/*.F90
```

## See Also

- FORD documentation: <https://github.com/Fortran-FOSS-Programmers/ford>
- Example projects: <https://github.com/Fortran-FOSS-Programmers/ford/wiki/Projects-Using-FORD>
