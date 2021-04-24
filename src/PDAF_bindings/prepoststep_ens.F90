!$Id: prepoststep_ens.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains. Also it is called once at the initial time
! before any forecasts are computed.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analized, e.g. by 
! computing the estimated variances. In addition, 
! the estimates can be written to disk. If a user 
! considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Check only
  use schism_glbl, only: npa,nvrt
  use mod_parallel_pdaf, only: mype_model,mype_filter,task_id,filterpe,mpierr,COMM_filter,MPI_REAL8, MPI_SUM
  use mod_assimilation, only: offset_field_p
  use output_schism_pdaf, only: write_netcdf_pdaf

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! Local vars
  integer i,j,member,field,istep
  INTEGER :: offset_field           ! Offset of a field in the state vector
  INTEGER :: dim_field              ! Dimension of a field
  REAL :: rmse_p(6)                ! PE-local estimated rms errors
  REAL :: rmse(6)                  ! Global estimated rms errors
  REAL, ALLOCATABLE :: var_p(:)     ! Estimated local model state variances
  CHARACTER(len=1) :: typestr       ! Character indicating call type

! !CALLING SEQUENCE:
! Called by: PDAF_get_state       (as U_prepoststep)
! Called by: PDAF_seik_update     (as U_prepoststep)
! Called by: PDAF_enkf_analysis   (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Called by: PDAF_estkf_analysis  (as U_prepoststep)
! Called by: PDAF_lestkf_analysis (as U_prepoststep)
! Called by: PDAF_etkf_analysis   (as U_prepoststep)
! Called by: PDAF_letkf_analysis  (as U_prepoststep)
! Called by: PDAF_lenkf_analysis  (as U_prepoststep)
! Called by: PDAF_netf_analysis   (as U_prepoststep)
! Called by: PDAF_lnetf_analysis  (as U_prepoststep)
!EOP

!write(*,*) 'In prepoststep_ens, check!',step,mype_world,task_id,filterpe

! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! Template reminder - delete when implementing functionality
!  WRITE (*,*) 'TEMPLATE prepoststep_ens.F90: Implement prepoststep here!'


! **********************
! *** INITIALIZATION ***
! **********************
  ! Allocate array for variances
  ALLOCATE(var_p(dim_p))

  ! Initialize numbers
  rmse_p = 0.0d0
  rmse   = 0.0d0

  IF (mype_filter==0) THEN
     IF (step==0) THEN
        !open rmse file
        open(70,file='rmse.dat')
        WRITE (*,'(a, i7,3x,a)') 'SCHISM-PDAF', step,'Analyze initial state ensemble'
        WRITE (typestr,'(a1)') 'i'
     ELSE IF (step>0) THEN
        WRITE (*,'(a, 8x,a)') 'SCHISM-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE IF (step<0) THEN
        WRITE (*,'(a, 8x,a)') 'SCHISM-PDAF', 'Analyze forecast state ensemble'
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF

!     *** Compute mean state ***
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'SCHISM-PDAF', '--- compute ensemble mean'
  state_p = 0.0d0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = state_p(:)/real(dim_ens,8)

! *** Compute local sampled variances of state vector ***
  var_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        var_p(j) = var_p(j) + &
             (ens_p(j, member) - state_p(j))* &
             (ens_p(j, member) - state_p(j))
     END DO
  END DO
  var_p(:) = var_p(:)/real(dim_ens-1,8)

! *** Compute different field RMSE
  offset_field = 0
  DO field = 1, 6
     ! Specify dimension of field
     ! SSH
     IF (field == 1) THEN
        dim_field = npa
     ! T,S,u,v,w
     ELSE
        dim_field = npa*nvrt
     END IF

     DO i = 1, dim_field
        rmse_p(field) = rmse_p(field) + var_p(i + offset_field)
     ENDDO
     rmse_p(field) = rmse_p(field) / real(dim_field,8)

     ! Set offset for next field
     offset_field = offset_field + dim_field

  END DO

  ! Global sum of RMS errors
  CALL MPI_Allreduce (rmse_p, rmse, 6, MPI_REAL8, MPI_SUM, &
       COMM_filter, MPIerr)

  rmse = SQRT(rmse)

  ! Display RMS errors
  IF (mype_filter==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'SCHISM-PDAF', 'RMS error according to sampled covariance'
     WRITE (*,'(a,7x,a9,1x,a14,a14,a14,a14,a14,/a, 10x,81a)') &
          'SCHISM-PDAF', 'ssh','temp','salt','u','v','w',&
          'SCHISM-PDAF', ('-',i=1,81)
     WRITE (*,'(a,10x,es11.4,5es14.4,1x,a5,a1,/a, 10x,81a)') &
          'SCHISM-PDAF', rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6), 'RMSe-', typestr,&
          'SCHISM-PDAF', ('-',i=1,81)
!    Write RMSE to file
     if (step==0) then
        WRITE (70,'(a, 10x,a)') &
              'SCHISM-PDAF', 'RMS error according to sampled covariance'
        WRITE (70,'(a,7x,a9,1x,a14,a14,a14,a14,a14,/a, 10x,81a)') &
              '   step    ', 'ssh','temp','salt','u','v','w',&
              '           ', ('-',i=1,81)
     end if
     WRITE (70,'(i,10x,es11.4,5es14.4,1x,a5,a1)') &
           step, rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6), 'RMSe-', typestr
  END IF

! **************************
! *** Write output nc files ***
! Here we only output ens-avg rank nc, follow schism output
! **************************
  istep=step
  CALL write_netcdf_pdaf(istep, dim_p, state_p)

! *************************
! *** Write output state_p for restart
! *************************
! Here we only output analysis for restart
! if (step>0) call write_ens_state(istep, dim_ens, ens_p) 

! Debug
! i=offset_field_p(2)+1
! if (mype_model.eq.2) then
! write(*,'(a,3f8.3,2i3,l)') 'In prepoststep_ens, check!',state_p(i:i+2),mype_model,task_id,filterpe
! end if

! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(var_p)


END SUBROUTINE prepoststep_ens
