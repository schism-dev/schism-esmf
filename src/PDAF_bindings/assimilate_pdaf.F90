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
  use schism_glbl, only: errmsg
  use schism_msgp, only: parallel_abort
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP

! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag !pass back to schism_cmi, be careful with int4

! Using simplified interface with standard routine name

! *********************************
! *** Call assimilation routine ***
! *********************************
! write(*,*) 'In assimilate_pdaf, check!'

! Disable local filter for dev
  IF (filtertype == 4) THEN
     CALL PDAF_put_state_etkf_si(status_pdaf)
! ELSEIF (filtertype == 5) THEN
!    CALL PDAF_put_state_letkf_si(status_pdaf)
  ELSEIF (filtertype == 6) THEN
     CALL PDAF_put_state_estkf_si(status_pdaf)
! ELSEIF (filtertype == 7) THEN
!    CALL PDAF_put_state_lestkf_si(status_pdaf)
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (errmsg,*) &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state(assimilate_pdaf) - stopping! (PE ', mype_world,')'
     CALL parallel_abort(errmsg)
  END IF

END SUBROUTINE assimilate_pdaf
