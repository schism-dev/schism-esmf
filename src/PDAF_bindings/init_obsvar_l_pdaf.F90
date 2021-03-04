!$Id: init_obsvar_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
! !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_l_pdaf(domain, step, dim_obs_l, obs_l, meanvar_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! This routine will only be called, if the 
! local adaptive forgetting factor feature 
! is used. Please note that this is an 
! experimental feature.
!
! The routine is called in the loop over all
! local analysis domains during each analysis
! by the routine PDAF\_set\_forget\_local that 
! estimates a local adaptive forgetting factor.
! The routine has to initialize the mean observation 
! error variance for the current local analysis 
! domain.  (See init_obsvar() for a global variant.)
!
! The routine is executed by all filter processes.
!
! !REVISION HISTORY:
! 2007-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, ONLY: rms_obs,rms_type,obs_coords_f,obs_index_l
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain        ! Current local analysis domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  REAL, INTENT(out)   :: meanvar_l     ! Mean local observation error variance

! Local vars
  integer :: i,idx
! !CALLING SEQUENCE:
! Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
!EOP


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Template reminder - delete when implementing functionality
! WRITE (*,*) 'TEMPLATE init_obsvar_l_pdaf.F90: Set mean observation variance here!'

  if (rms_type==1) then
     meanvar_l = rms_obs**2
  else
     meanvar_l=0.d0
     do i=1,dim_obs_l
        idx=obs_index_l(i)
        meanvar_l=meanvar_l+obs_coords_f(4,idx)**2
     end do
     meanvar_l=meanvar_l/dim_obs_l
  end if

END SUBROUTINE init_obsvar_l_pdaf
