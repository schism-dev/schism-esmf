!$Id: assimilate_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

! !DESCRIPTION:
! This routine is called during the model integrations at each time 
! step. It check whether the forecast phase is completed. If so, 
! PDAF_put_state_X is called to perform the analysis step.
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code for NEMO
! Later revisions - see svn log
!
! !USES:
  USE pdaf_interfaces_module, ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
         PDAFomi_assimilate_lenkf,&
         PDAF_get_localfilter

  use schism_glbl, only: errmsg
  use schism_msgp, only: parallel_abort
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, task_id,filterpe
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP
! Interface between model and PDAF, and prepoststep
  EXTERNAL :: next_observation_pdaf,  & ! Provide time step, model time of next observation
              distribute_state_pdaf,  & ! Routine to distribute a state vector to model fields
              collect_state_pdaf,     & ! Collect a state vector from model fields
              prepoststep_ens           ! User supplied pre/poststep routine
! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf,    & ! Provide number of local analysis domains
              init_dim_l_pdaf,        & ! Initialize state dimension for local ana. domain
              g2l_state_pdaf,         & ! Get state on local ana. domain from global state
              l2g_state_pdaf            ! Init global state from state on local analysis domain
! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi,   & ! Get dimension of full obs. vector for PE-local domain
              obs_op_pdafomi,         & ! Obs. operator for full obs. vector for PE-local domain
              init_dim_obs_l_pdafomi, & ! Get dimension of obs. vector for local analysis domain
              localize_covar_pdafomi    ! Apply localization to covariance matrix in LEnKF


! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag 
  INTEGER localfilter          ! Flag whether the chosen filter is localized
  integer :: nsteps,doexit
  real :: timenow

! Using simplified interface with standard routine name

! *********************************
! *** Call assimilation routine ***
! *********************************
! write(*,*) 'In assimilate_pdaf, check!',mype_world, task_id,filterpe
! Update state
! write(*,*) 'Before get_state',nsteps,timenow, doexit
! call PDAF_get_state(nsteps,timenow, doexit, next_observation_pdaf, distribute_state_pdaf, prepoststep_ens, status_pdaf)
! write(*,*) 'after get_state',nsteps,timenow, doexit

! Check whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)
! write(*,*) 'in assimilate_pdaf, localfilter=',localfilter

  IF (localfilter == 1) THEN
     CALL PDAFomi_put_state_local_si(status_pdaf)
!    write(*,*) 'in assimilate_pdaf, afterlocal=',localfilter
  ELSE
     IF (filtertype/=8) THEN
!        All global filters, except LEnKF
         CALL PDAFomi_put_state_global_si(status_pdaf)
     ELSE
!        LEnKF has its own OMI interface routine
         CALL PDAFomi_put_state_lenkf_si(status_pdaf)
     END IF
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (errmsg,*) &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state(assimilate_pdaf) - stopping! (PE ', mype_world,')'
     CALL parallel_abort(errmsg)
  END IF

! Update state
! write(*,*) 'Before get_state',nsteps,timenow, doexit
! call PDAF_get_state(nsteps,timenow, doexit, next_observation_pdaf, distribute_state_pdaf, prepoststep_ens, status_pdaf)
! write(*,*) 'after get_state',nsteps,timenow, doexit

END SUBROUTINE assimilate_pdaf
