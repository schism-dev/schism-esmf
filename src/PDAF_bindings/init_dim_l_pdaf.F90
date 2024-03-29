!$Id: init_dim_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the loop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use schism_glbl,only : nvrt,ntracers,xnd,ynd,xlon,ylat,ics,pi
  use mod_assimilation, only: coords_l
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
! Called by: PDAF_lnetf_update   (as U_init_dim_l)
!EOP


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Template reminder - delete when implementing functionality
! WRITE (*,*) 'TEMPLATE init_dim_l_pdaf.F90: Set local state dimension here!'

! ssh + t,s,u,v,w

  dim_l = 1+(ntracers+3)*nvrt

  if (ics==2) then
     coords_l(1)=xlon(domain_p)/pi*180.d0
     coords_l(2)=ylat(domain_p)/pi*180.d0
  else
     coords_l(1)=xnd(domain_p)
     coords_l(2)=ynd(domain_p)
  end if


END SUBROUTINE init_dim_l_pdaf
