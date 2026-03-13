# FORD Documentation Example

This document shows a before-and-after example of adding FORD documentation to `schism_bmi.F90`.

## Original Code (Undocumented)

```fortran
module schism_bmi

  use esmf

  interface
    subroutine parallel_init(communicator)
      implicit none
      integer, optional :: communicator
    end subroutine parallel_init

    subroutine schism_init(iorder,indir,iths,ntime)
      implicit none
      integer, intent(in)          :: iorder
      character(len=*), intent(in) :: indir
      integer, intent(out)         :: iths,ntime
    end subroutine schism_init

    subroutine schism_step(it)
      implicit none
      integer, intent(in) :: it
    end subroutine schism_step
  end interface

contains

subroutine schism_parallel_init(communicator)
  use schism_msgp, only: parallel_init
  implicit none
  integer :: communicator
  call parallel_init(communicator)
end subroutine schism_parallel_init

end module schism_bmi
```

## Documented Version (FORD-ready)

```fortran
!> @brief Basic Model Interface (BMI) for SCHISM ocean model
!!
!! This module provides a generic Basic Model Interface (BMI) to SCHISM,
!! enabling external control through standardized function calls. The interface
!! wraps SCHISM's internal parallel initialization, time-stepping, and
!! finalization routines.
!!
!! The BMI allows SCHISM to be coupled with other models through frameworks
!! like ESMF, NUOPC, or standalone coupling applications.
!!
!! @note This module accounts for halo (ghost) zones when used with ESMF,
!!       as ESMF partitions among nodes rather than elements by default.
!!
!! @author Carsten Lemmen <carsten.lemmen@hereon.de>
!! @author Richard Hofmeister
!! @copyright 2021-2023 Helmholtz-Zentrum Hereon
!! @copyright 2018-2021 Helmholtz-Zentrum Geesthacht
!! @license Apache-2.0
!!
!! @see https://csdms.colorado.edu/wiki/BMI for BMI specification
module schism_bmi

  use esmf

  !> External interface to SCHISM's MPI initialization routine
  !!
  !! @note The actual implementation is in SCHISM's libcore
  interface
    !> @brief Initialize parallel (MPI) environment for SCHISM
    !!
    !! @param[in] communicator Optional MPI communicator handle
    subroutine parallel_init(communicator)
      implicit none
      integer, optional :: communicator
    end subroutine parallel_init

    !> @brief Initialize SCHISM model state and read input files
    !!
    !! Reads the model configuration, grid, initial conditions, and boundary
    !! conditions from the specified input directory.
    !!
    !! @param[in]  iorder  Time integration order (1 or 2)
    !! @param[in]  indir   Input directory path containing param.nml and grid files
    !! @param[out] iths    Number of tracers in the model
    !! @param[out] ntime   Total number of time steps to run
    subroutine schism_init(iorder,indir,iths,ntime)
      implicit none
      integer, intent(in)          :: iorder
      character(len=*), intent(in) :: indir
      integer, intent(out)         :: iths,ntime
    end subroutine schism_init

    !> @brief Advance SCHISM model by one time step
    !!
    !! Executes one complete time step of the SCHISM model, including:
    !! - Solving hydrodynamic equations
    !! - Advancing tracers
    !! - Computing fluxes
    !! - Writing output if scheduled
    !!
    !! @param[in] it Time step number (1-based index)
    subroutine schism_step(it)
      implicit none
      integer, intent(in) :: it
    end subroutine schism_step
  end interface

contains

!> @brief Initialize SCHISM's parallel environment with given communicator
!!
!! This wrapper ensures that SCHISM uses the same MPI communicator as the
!! calling framework (e.g., ESMF), enabling proper multi-component coupling.
!!
!! @param[in] communicator MPI communicator handle from ESMF or other framework
!!
!! @note Must be called before schism_init()
!! @warning The communicator must remain valid for the lifetime of the SCHISM model
!!
!! @author Carsten Lemmen
subroutine schism_parallel_init(communicator)

  use schism_msgp, only: parallel_init
  use schism_msgp, only: schism_mpi_comm=>comm

  implicit none
  
  integer :: communicator  !< MPI communicator provided by coupling framework

  call parallel_init(communicator)

end subroutine schism_parallel_init

end module schism_bmi
```

## Key Documentation Elements Added

### Module-level Documentation

- **Purpose statement**: Explains what the module does
- **Context**: Where it fits in the overall system
- **Authors and copyright**: Proper attribution
- **Cross-references**: Links to related documentation (`@see`)

### Interface Documentation

Each external subroutine now has:
- **Brief description**: One-line summary with `@brief`
- **Detailed explanation**: Multi-paragraph description of behavior
- **Parameter documentation**: Each parameter with `@param[in/out]`
- **Important notes**: Using `@note` and `@warning`

### Subroutine Documentation

The wrapper subroutine includes:
- **Purpose**: What it does and why
- **Parameters**: Inline documentation with `!<` for each variable
- **Preconditions**: Using `@note`
- **Warnings**: Using `@warning`
- **Attribution**: Original author with `@author`

## How FORD Uses These Comments

When you run `ford ford.yml`, it will:

1. **Parse** all `!>` and `!!` comments before declarations
2. **Extract** `@param`, `@return`, `@author`, etc. tags
3. **Generate**:
   - Module documentation page
   - Call graphs showing relationships
   - Cross-linked source code browser
   - Search index

## Testing the Documentation

```bash
# Generate API docs
ford ford.yml

# Open in browser
open build/ford_docs/index.html

# You should see:
# - List of modules (schism_bmi, etc.)
# - Module page with full description
# - Procedure list with brief summaries
# - Call graphs if enabled
```

## Quick Reference Card

| Purpose | Syntax |
|---------|--------|
| Document next item | `!>` |
| Continue documentation | `!!` |
| Document same-line item | `!<` |
| Brief summary | `@brief` |
| Input parameter | `@param[in] name Description` |
| Output parameter | `@param[out] name Description` |
| In/out parameter | `@param[inout] name Description` |
| Return value | `@return Description` |
| Author | `@author Name <email>` |
| Copyright | `@copyright Year Entity` |
| License | `@license License-Name` |
| Important note | `@note Text` |
| Warning | `@warning Text` |
| Known bug | `@bug Description` |
| Future work | `@todo Description` |
| Cross-reference | `@see other_module` |
| External link | `@see http://example.com` |

## Next Steps

1. Copy this pattern to other modules in `src/schism/`
2. Document the main ESMF cap in `src/schism/schism_esmf_cap.F90`
3. Add documentation to driver routines in `src/driver/`
4. Run `ford ford.yml` to verify syntax
5. Fix any warnings FORD reports
6. Integrate into CI/CD pipeline

See the full [API Documentation Guide](api-documentation-guide.md) for more examples and best practices.
