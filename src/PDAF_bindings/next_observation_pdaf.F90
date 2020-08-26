!$Id: next_observation_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The subroutine is called before each forecast phase
! by PDAF\_get\_state. It has to initialize the number 
! of time steps until the next available observation 
! (nsteps) and the current model time (time). In 
! addition the exit flag (exit) has to be initialized.
! It indicates if the data assimilation process is 
! completed such that the ensemble loop in the model 
! routine can be exited.
!
! The routine is called by all filter processes. 
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use schism_glbl, only : time_stamp
  use mod_assimilation, only : delt_obs
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  ! Template reminder - delete when implementing functionality
!  WRITE (*,*) 'TEMPLATE next_observation_pdaf.F90: Set number of time step in forecast!'

!  write(*,*) 'In next_observation_pdaf, check!',mype_world,task_id,filterpe

   nsteps = delt_obs

! *********************************
! *** Set current physical time ***
! *********************************

   time = real(time_stamp) ! time_stamp is real*8

! *********************
! *** Set exit flag ***
! *********************

   doexit = 0

END SUBROUTINE next_observation_pdaf
