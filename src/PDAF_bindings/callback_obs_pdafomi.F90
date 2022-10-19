!$Id: callback_obs_pdafomi.F90 883 2021-11-27 14:16:40Z lnerger $
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  use schism_glbl, only: errmsg
  use schism_msgp, only: parallel_abort
! use mod_assimilation, only: dim_obs_A,dim_obs_Z,dim_obs_S,dim_obs_T,dim_obs_U,dim_obs_V
  USE obs_A_pdafomi, ONLY: assim_A, init_dim_obs_A
  USE obs_Z_pdafomi, ONLY: assim_Z, init_dim_obs_Z
  USE obs_S_pdafomi, ONLY: assim_S, init_dim_obs_S
  USE obs_T_pdafomi, ONLY: assim_T, init_dim_obs_T
  USE obs_U_pdafomi, ONLY: assim_U, init_dim_obs_U
  USE obs_V_pdafomi, ONLY: assim_V, init_dim_obs_V
  use mod_parallel_pdaf, only: mype_world !Debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_p_all ! dimensions for all types of observation
  INTEGER :: dim_obs_A,dim_obs_Z,dim_obs_S,dim_obs_T,dim_obs_U,dim_obs_V


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

! WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_pdafomi: complete interface to observation modules'

  ! Initialize number of observations
   dim_obs_A=0
   dim_obs_Z=0
   dim_obs_S=0
   dim_obs_T=0
   dim_obs_U=0
   dim_obs_V=0

  !Get all observations and its index routine 
   call init_dim_obs_all(step,dim_obs_p_all)
!  write(*,*) 'step, dim_obs_p_all=',step,dim_obs_p_all,mype_world
! dim_obs_A=dim_obstype_mod(1)
! dim_obs_Z=dim_obstype_mod(2)
! dim_obs_S=dim_obstype_mod(3)
! dim_obs_T=dim_obstype_mod(4)
! dim_obs_U=dim_obstype_mod(5)
! dim_obs_V=dim_obstype_mod(6)
! if (dim_obs_A.eq.0) assim_A=.false.
! if (dim_obs_Z.eq.0) assim_Z=.false.
! if (dim_obs_S.eq.0) assim_S=.false.
! if (dim_obs_T.eq.0) assim_T=.false.
! if (dim_obs_U.eq.0) assim_U=.false.
! if (dim_obs_V.eq.0) assim_V=.false.

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A) CALL init_dim_obs_A(step, dim_obs_A)
  IF (assim_Z) CALL init_dim_obs_Z(step, dim_obs_Z)
  IF (assim_S) CALL init_dim_obs_S(step, dim_obs_S)
  IF (assim_T) CALL init_dim_obs_T(step, dim_obs_T)
  IF (assim_U) CALL init_dim_obs_U(step, dim_obs_U)
  IF (assim_V) CALL init_dim_obs_V(step, dim_obs_V)

  dim_obs = dim_obs_A + dim_obs_Z + dim_obs_S + dim_obs_T + dim_obs_U + dim_obs_V 
! if (dim_obs_p_all.ne.dim_obs) then
!    WRITE (errmsg,*)'obs count number not match error ', dim_obs, dim_obs_p_all &
!         ' in init_dim_obs_pdafomi - stopping! '
!    CALL parallel_abort(errmsg)
! end if

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_OBSTYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: obs_op_A
  USE obs_Z_pdafomi, ONLY: obs_op_Z
  USE obs_S_pdafomi, ONLY: obs_op_S
  USE obs_T_pdafomi, ONLY: obs_op_T
  USE obs_U_pdafomi, ONLY: obs_op_U
  USE obs_V_pdafomi, ONLY: obs_op_V
! use mod_assimilation, only: dim_obs_A,dim_obs_Z,dim_obs_S,dim_obs_T,dim_obs_U,dim_obs_V

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state

! Local vars
! integer :: ic1,ic2 ! counter


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

! WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/obs_op_pdafomi: complete interface to observation modules'

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi
! ic1=0
! ic2=0
! if (dim_obs_A.gt.0) then
!    ic1=1
!    ic2=dim_obs_A
! end if
  CALL obs_op_A(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_Z(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_S(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_T(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_U(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_V(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: init_dim_obs_l_A
  USE obs_Z_pdafomi, ONLY: init_dim_obs_l_Z
  USE obs_S_pdafomi, ONLY: init_dim_obs_l_S
  USE obs_T_pdafomi, ONLY: init_dim_obs_l_T
  USE obs_U_pdafomi, ONLY: init_dim_obs_l_U
  USE obs_V_pdafomi, ONLY: init_dim_obs_l_V

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! Local vars
! integer :: dim_local_sum


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

! WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_l_pdafomi: complete interface to observation modules'

  ! Call init_dim_obs_l specific for each observation
  CALL init_dim_obs_l_A(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_Z(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_S(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_T(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_U(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_V(domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type to apply covariance
!! localization in the LEnKF.
!!
SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: localize_covar_A
  USE obs_Z_pdafomi, ONLY: localize_covar_Z
  USE obs_S_pdafomi, ONLY: localize_covar_S
  USE obs_T_pdafomi, ONLY: localize_covar_T
  USE obs_U_pdafomi, ONLY: localize_covar_U
  USE obs_V_pdafomi, ONLY: localize_covar_V
  use schism_glbl,only : xnd,ynd,xlon,ylat,ics,pi,npa,nvrt,ntracers,errmsg
  use schism_msgp, only: parallel_abort

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs               !< number of observations
  REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

! *** local variables ***
  REAL, ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries
  integer :: i,j,k,ic


! **********************
! *** INITIALIZATION ***
! **********************

! WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/localize_covar_pdafomi: complete interface to observation modules'

  ! Initialize coordinate array

  ! One needs to provide the array COORDS_P holding the coordinates of each
  ! element of the process-local state vector. Each column of the array holds
  ! the information for one element. The array can be initialized here using
  ! information on the model grid.

  ALLOCATE(coords_p(2, dim_p))
  ic=0
  do i=1,npa
     ic=ic+1
     if (ics==2) then
         coords_p(1,ic)=xlon(i)/pi*180.d0
         coords_p(2,ic)=ylat(i)/pi*180.d0
     else
         coords_p(1,ic)=xnd(i)
         coords_p(2,ic)=ynd(i)
     end if
  end do

  do j=1,ntracers+3 ! + 3 = U,V,W
     do i=1,npa
        do k=1,nvrt
           ic=ic+1
           if (ics==2) then
              coords_p(1,ic)=xlon(i)/pi*180.d0
              coords_p(2,ic)=ylat(i)/pi*180.d0
           else
              coords_p(1,ic)=xnd(i)
              coords_p(2,ic)=ynd(i)
           end if
        end do
     end do
  end do
 
  if (ic.ne.dim_p) then
     WRITE (errmsg,*)'localize_covar dim_p not match, ', dim_p, ic
     CALL parallel_abort(errmsg)
  end if

! *************************************
! *** Apply covariance localization ***
! *************************************

  ! Call localize_covar specific for each observation
  CALL localize_covar_A(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_Z(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_S(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_T(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_U(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_V(dim_p, dim_obs, HP_p, HPH, coords_p)


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(coords_p)

END SUBROUTINE localize_covar_pdafomi

!-------------------------------------------------------------------------------
!> Call-back routine for deallocate_obs
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!! for each observation type
!!
SUBROUTINE deallocate_obs_pdafomi()

  ! Include observation types
  USE obs_A_pdafomi, ONLY: deallocate_obs_A
  USE obs_Z_pdafomi, ONLY: deallocate_obs_Z
  USE obs_S_pdafomi, ONLY: deallocate_obs_S
  USE obs_T_pdafomi, ONLY: deallocate_obs_T
  USE obs_U_pdafomi, ONLY: deallocate_obs_U
  USE obs_V_pdafomi, ONLY: deallocate_obs_V

  IMPLICIT NONE


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL deallocate_obs_A()
  CALL deallocate_obs_Z()
  CALL deallocate_obs_S()
  CALL deallocate_obs_T()
  CALL deallocate_obs_U()
  CALL deallocate_obs_V()

END SUBROUTINE deallocate_obs_pdafomi

