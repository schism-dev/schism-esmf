!$Id: init_dim_obs_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analysis domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY: local_range, coords_obs_f, obs_index_l, distance_l
  use schism_glbl, only: xnd,ynd
  use mod_assimilation, only: distance_l,local_range,obs_index_l,obs_coords_f
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! Local vars
  real distance
  integer cnt,i

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
! Called by: PDAF_lnetf_update   (as U_init_dim_l)
!EOP


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************


! Count observations within local_range
  dim_obs_l = 0
  DO i = 1, dim_obs_f
      distance = SQRT((xnd(domain_p) - obs_coords_f(1,i))**2 + &
                      (ynd(domain_p) - obs_coords_f(2,i))**2)
     IF (distance <= local_range) dim_obs_l = dim_obs_l + 1
  END DO

! Initialize index array for local observations in full observed vector & array
! of distances of local observations
  IF (ALLOCATED(distance_l)) DEALLOCATE(distance_l)
  IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
  ALLOCATE(distance_l(dim_obs_l))
  ALLOCATE(obs_index_l(dim_obs_l))

  cnt = 0
  DO i = 1, dim_obs_f
     distance = SQRT((xnd(domain_p) - obs_coords_f(1,i))**2 + &
                     (ynd(domain_p) - obs_coords_f(2,i))**2)
     IF (distance <= local_range) THEN
         cnt = cnt + 1
         distance_l(cnt) = distance
         obs_index_l(cnt) = i
     END IF
  END DO
  

END SUBROUTINE init_dim_obs_l_pdaf

