!$Id: init_ens_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens_pdaf --- Initialize ensemble for filter
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! If only a single filter algorithm is used, the 
! ensemble initialization can be performed directly
! in this routine. If a single filter is implemented,
! one can perform the initialization directly here.
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the particular
! ensemble initialization routine for the selected filter.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Check only
  use mod_parallel_pdaf, only: mype_model,task_id,filterpe
! new28
! The following is only use for lock-exchange experiment, delete them after done
  use mod_assimilation, only: offset_field_p
  use schism_glbl, only: npa,nvrt
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! Local variables
  INTEGER ::  member,i,is,ie  ! Counters

! !CALLING SEQUENCE:
! Called by: PDAF_init       (as U_init_ens)
! Calls: init_seik
! Calls: init_seek
! Calls: init_enkf
!EOP

! write(*,*) 'In init_ens_pdaf, check!',mype_world,task_id,filterpe

! *** Read ensemble 
  CALL collect_state_pdaf(dim_p, state_p)

  is=offset_field_p(2)+1
  ie=offset_field_p(2)+nvrt*npa
! Debug
! write(*,'(a,6f8.2,i4)') 'state_p in ens',state_p(is-1:is+1),state_p(ie-1:ie+1),mype_model
! write(*,'(a,f6.2,2i4,l2)') 'In init_ens_pdaf, state_p(max)',maxval(state_p),kind(state_p),mype_world,task_id,filterpe
! *** Initialize ens_p
  DO member = 1,dim_ens
        ens_p(:,member) = state_p(:) !+ float(member-1)*1.d-2 ! temporary, add fesom example later
!       do i=is,ie ! temperature
           state_p(is:ie)=state_p(is:ie)+0.2d0 !this is just for lock-exchange test, will remove after test done!
!       end do
  END DO
! write(*,*) 'in init_ens_pdaf',state_p(1:3)
! FESOM example is followed by init_seik, need to gen_cov first, then use its
! output to add pertubation into ensembles



! Disable, just collect
! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

! IF (filtertype == 0) THEN
!    ! EOF initialization for SEEK
!    CALL init_seek(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! ELSE IF (filtertype == 2) THEN
!    ! Use random sampling initialization
!    CALL init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! ELSE
!    ! Use 2nd-order exact sampling
!    CALL init_seik(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! END IF

END SUBROUTINE init_ens_pdaf
