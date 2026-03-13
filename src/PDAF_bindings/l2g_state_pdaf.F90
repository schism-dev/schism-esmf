!$Id: l2g_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: l2g_state_pdaf --- Initialize full state from local analysis
!
! !INTERFACE:
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_X\_update 
! after the analysis and ensemble transformation 
! on a single local analysis domain. It has to 
! initialize elements of the PE-local full state 
! vector from the provided analysis state vector 
! on the local analysis domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, ONLY: offset_field_p,idx_domain_p,Vlocal_opt,vidx_domain_p
  use schism_glbl,only : nvrt,ntracers,npa

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p         ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

! Local vars
  integer iid,cnt,nfield,k,ii,ic,offset_p(6)
  integer i,id,ilayer
  real tmp,wt

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
! Called by: PDAF_lestkf_update   (as U_l2g_state)
! Called by: PDAF_letkf_update    (as U_l2g_state)
! Called by: PDAF_lnetf_update    (as U_l2g_state)
!EOP


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Template reminder - delete when implementing functionality
!  WRITE (*,*) 'TEMPLATE l2g_state_pdaf.F90: Set part of global state vector here!'

!  state_p = ?
  
  if (Vlocal_opt.eq.1) then ! 2D+1D
     id=idx_domain_p(2,domain_p)
     state_p(id)=state_l(dim_l)

  else !2D
     ! local offset, because we only have 5 fields in offset_field_p
     offset_p(1) = offset_field_p(1) !elev
     offset_p(2) = offset_field_p(2) !tr_nd(1),temp
     offset_p(3) = offset_field_p(2)+npa*nvrt !tr_nd(2), salt
     offset_p(4) = offset_field_p(3) !uu
     offset_p(5) = offset_field_p(4) !vv
     offset_p(6) = offset_field_p(5) !ww

     nfield=5 ! z+tracer(t/s)+u+v+w

     ! elev
     state_p(domain_p)=state_l(1)
     ! 3D field
     cnt=2
     do iid=2,nfield+1
       if (vidx_domain_p(domain_p).ne.1) then
         ilayer=nvrt-vidx_domain_p(domain_p)
         do k=1,vidx_domain_p(domain_p)-1 !from specified bottom to surface
            !state_p( (domain_p-1)*nvrt + offset_p(iid) + vidx_domain_p(domain_p) ) = state_l(cnt) 
            cnt=cnt+1 !loop cnt is necessary
         end do !k

         do k=vidx_domain_p(domain_p),nvrt !from specified bottom to surface
            if ((ilayer.gt.10).and.((k-vidx_domain_p(domain_p).le.7))) then !7 layers to reduced sudden jump 
               wt=real(k-vidx_domain_p(domain_p))/7.d0
            else
               wt=1.d0
            end if
            tmp=state_p( (domain_p-1)*nvrt + offset_p(iid) + k )*(1.d0-wt)+state_l(cnt)*wt
            !state_p( (domain_p-1)*nvrt + offset_p(iid) + k ) = state_l(cnt) 
            state_p( (domain_p-1)*nvrt + offset_p(iid) + k ) = tmp
            cnt=cnt+1
         end do !k

       else
         do k=1,nvrt
            state_p( (domain_p-1)*nvrt + offset_p(iid) + k ) = state_l(cnt) 
            cnt=cnt+1
         end do !k
       end if
     end do !iid
  end if

END SUBROUTINE l2g_state_pdaf
