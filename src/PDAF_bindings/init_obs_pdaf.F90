!$Id: init_obs_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_obs_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_pdaf(step, dim_obs_p, observation_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step. 
! It has to provide the PE-local observation vector 
! for the current time step.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, ONLY: obs_p
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step             ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p        ! PE-local dimension of obs. vector
  REAL, INTENT(out)   :: observation_p(dim_obs_p) ! PE-local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_obs_ensemble
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis
! Called by: PDAF_netf_analysis
!EOP
!  write(*,*) 'In init_obs_pdaf, check!',mype_world,task_id,filterpe


! ***************************************************************
! *** Initialize observation vector for PE-local model domain ***
! ***************************************************************
  
     ! This is generic because it is just using
     ! the vector obs filled in init_dim_obs_f

  observation_p(:) = obs_p(:)


END SUBROUTINE init_obs_pdaf

