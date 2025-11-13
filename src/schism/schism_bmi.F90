!> @brief Basic Model Interface (BMI) for SCHISM ocean model
!!
!! This module provides a generic Basic Model Interface (BMI) to SCHISM,
!! enabling external control through standardized function calls. The interface
!! wraps SCHISM's internal parallel initialization, time-stepping, and
!! finalization routines, as well as providing access to model state variables.
!!
!! The BMI allows SCHISM to be coupled with other models through frameworks
!! like ESMF, NUOPC, or standalone coupling applications.
!!
!! ## Key Features
!!
!! - Parallel (MPI) initialization with custom communicator
!! - Access to model time step size
!! - Pointer access to selected internal variables (e.g., air temperature)
!! - Accounts for halo (ghost) zones when used with ESMF
!!
!! @note This module accounts for halo (ghost) zones when used with ESMF,
!!       as ESMF partitions among nodes rather than elements by default.
!!
!! @note While this module uses ESMF for logging only, it should ideally
!!       not depend on ESMF to maintain BMI portability.
!!
!! @author Carsten Lemmen
!! @author Richard Hofmeister
!! @license Apache-2.0
!! @see https://csdms.colorado.edu/wiki/BMI
!! @see module schism_esmf_cap


#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_bmi.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_bmi

  use esmf  !< ESMF library (used for logging only)

  !> External interfaces to SCHISM core routines
  !!
  !! These interfaces define the signatures of SCHISM's internal routines
  !! that are implemented in the SCHISM core library (libcore.a).
  !!
  !! @note The actual implementations are in SCHISM's Fortran source code,
  !!       not in this module.
  interface

    !> @brief Initialize parallel (MPI) environment for SCHISM
    !!
    !! Initializes the MPI environment for SCHISM using an optional
    !! communicator provided by the calling framework.
    !!
    !! @param[in] communicator Optional MPI communicator handle. If not provided,
    !!                         MPI_COMM_WORLD is used by default.
    subroutine parallel_init(communicator)
      implicit none
      integer, optional :: communicator
    end subroutine parallel_init

    !> @brief Finalize parallel (MPI) environment for SCHISM
    !!
    !! Cleans up MPI resources used by SCHISM. Should be called before
    !! MPI_Finalize in the main program.
    subroutine parallel_finalize
      implicit none
    end subroutine parallel_finalize

    !> @brief Initialize SCHISM model state and read input files
    !!
    !! Reads the model configuration, grid, initial conditions, and boundary
    !! conditions from the specified input directory. This must be called after
    !! parallel_init() and before any time-stepping.
    !!
    !! @param[in]  iorder  Time integration order (1=Euler, 2=Runge-Kutta)
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

    !> @brief Finalize SCHISM model and clean up resources
    !!
    !! Closes open files, deallocates memory, and performs cleanup operations.
    !! Should be called before parallel_finalize().
    subroutine schism_finalize()
      implicit none
    end subroutine schism_finalize

  end interface

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "schism_parallel_init"
!> @brief Initialize SCHISM's parallel environment with given communicator
!!
!! This wrapper ensures that SCHISM uses the same MPI communicator as the
!! calling framework (e.g., ESMF), enabling proper multi-component coupling.
!! The communicator is passed to SCHISM's internal parallel initialization
!! routine.
!!
!! @param[in] communicator MPI communicator handle from ESMF or other framework
!!
!! @note Must be called before schism_init()
!! @warning The communicator must remain valid for the lifetime of the SCHISM model
subroutine schism_parallel_init(communicator)

  use schism_msgp, only: parallel_init
  use schism_msgp, only: schism_mpi_comm=>comm

  implicit none
  
  integer :: communicator  !< MPI communicator provided by coupling framework

  call parallel_init(communicator)

end subroutine schism_parallel_init

#undef  ESMF_METHOD
#define ESMF_METHOD "schismTimeStep"
!> @brief Get SCHISM model time step size in seconds
!!
!! Returns the current time step size used by SCHISM. This value is read
!! from the model configuration (param.nml) during initialization.
!!
!! @param[out] seconds Time step size in seconds
!!
!! @note This should be called after schism_init() when the time step
!!       is known.
subroutine schismTimeStep(seconds)

  use schism_glbl, only: dt

  implicit none
  
  double precision, intent(out) :: seconds  !< Time step size in seconds

  seconds = dt

end subroutine schismTimeStep

#undef ESMF_METHOD
#define ESMF_METHOD 'schismPtr1'
!> @brief Get pointer to SCHISM internal variable by name
!!
!! Returns a pointer to one of SCHISM's internal state variables, allowing
!! external code to read or modify model state. This is useful for:
!! - Coupling with other models
!! - Data assimilation
!! - Custom diagnostics
!!
!! Currently supported variables:
!! - `airt2`: Air temperature at 2m height [K]
!!
!! @param[in] varname Name of the variable to access
!!
!! @return Pointer to the requested variable array, or null pointer if not found
!!
!! @warning Modifying internal variables directly can lead to inconsistent
!!          model state. Use with caution.
!!
!! @note More variables can be added to the select case as needed
function schismPtr1(varname) result(farrayPtr)

  use schism_glbl, only: airt2

  implicit none
  
  character(len=*), intent(in) :: varname      !< Name of variable to retrieve
  double precision, pointer :: farrayPtr(:)    !< Pointer to variable array

  select case(trim(varname))
  case ('airt2')
    farrayPtr => airt2  !< Air temperature at 2m [K]
  case default
    nullify(farrayPtr)  !< Variable not found, return null pointer
  end select

end function schismPtr1

end module schism_bmi
