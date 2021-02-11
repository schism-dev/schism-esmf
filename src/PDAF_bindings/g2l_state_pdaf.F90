!$Id: g2l_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, ONLY: offset_field_p
  use schism_glbl,only : nvrt,ntracers,npa
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

! Local vars
  integer iid,cnt,nfield,k,ii,ic,offset_p(6)
! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
!EOP


! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! Template reminder - delete when implementing functionality
! WRITE (*,*) 'TEMPLATE g2l_state_pdaf.F90: Initialize local state vector here!'

!   state_l = ??

! field, start with z/tracer(t/s)/u/v/w
! local offset, because we only have 5 fields in offset_field_p
  offset_p(1) = offset_field_p(1) !elev
  offset_p(2) = offset_field_p(2) !tr_nd(1),temp
  offset_p(3) = offset_field_p(2)+npa*nvrt !tr_nd(2), salt
  offset_p(4) = offset_field_p(3) !uu
  offset_p(5) = offset_field_p(4) !vv
  offset_p(6) = offset_field_p(5) !ww


  nfield=5 ! z+tracer(t/s)+u+v+w

! elev
  state_l(1)=state_p(domain_p)
! 3D field
  cnt=2
  do iid=2,nfield+1
      do k=1,nvrt
         state_l(cnt)= state_p( (domain_p-1)*nvrt + offset_p(iid) + k  )
         cnt=cnt+1
      end do !k
  end do !iid
! do iid=2,nfield
!    ic=1 !iid=3~5
!    if (iid.eq.2) ic=ntracers
!       do ii=1,ic
!          do k=1,nvrt
!            state_l(cnt)= state_p(domain_p + offset_field_p(iid) + (k-1)*npa + nvrt*npa*(ii-1))
!            cnt=cnt+1
!          end do !k
!       end do !ii
! end do !iid
! write(*,'(a, i,106f10.3)') 'g2l_state:', ntracers,state_l

END SUBROUTINE g2l_state_pdaf
