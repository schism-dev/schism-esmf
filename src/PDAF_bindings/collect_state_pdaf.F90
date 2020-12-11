!$Id: collect_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
! after the propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use schism_glbl, only: nea,nsa,npa,nvrt,ntracers,idry_e,we,tr_el, &
 &idry_s,su2,sv2, idry,eta2,tr_nd,uu2,vv2,ww2,kbp
! Check only
  use mod_parallel_pdaf, only: mype_model,task_id,filterpe
  use mod_assimilation, only: offset_field_p

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

  integer :: i,j,k,itot

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP
  
! write(*,*) 'In collect_state_pdaf, check!',mype_world,task_id,filterpe

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

!   Assign state vars to state_p in the resident domains

   !init count
   itot=0

   !node-base
   ! Elev
   do i=1,npa
     itot=itot+1
     state_p(itot)=eta2(i)
   enddo !i
   ! Tracer
   do j=1,ntracers
     do i=1,npa
       do k=1,nvrt
         itot=itot+1
         if (k.lt.kbp(i)) then
            state_p(itot)= -9999.d0 ! try fill -9999 at below bottom 
         else
            state_p(itot)=tr_nd(j,k,i)
         end if
       enddo !k
     enddo !i
   enddo !j
   ! U
   do i=1,npa
     do k=1,nvrt
        itot=itot+1
        if (k.lt.kbp(i)) then
           state_p(itot)= -9999.d0 ! try fill -9999 at below bottom 
        else
           state_p(itot)=uu2(k,i)
        end if
     enddo !k
   enddo !i
   ! V
   do i=1,npa
     do k=1,nvrt
        itot=itot+1
        if (k.lt.kbp(i)) then
           state_p(itot)= -9999.d0 ! try fill -9999 at below bottom 
        else
           state_p(itot)=vv2(k,i)
        end if
     enddo !k
   enddo !i
   ! W
   do i=1,npa
     do k=1,nvrt
        itot=itot+1
        if (k.lt.kbp(i)) then
           state_p(itot)= -9999.d0 ! try fill -9999 at below bottom 
        else
           state_p(itot)=ww2(k,i)
        end if
     enddo !k
   enddo !i
 
!  Debug
!  i=offset_field_p(2)+1
!  if (mype_model.eq.2) then
!  write(*,'(a,3f8.3,2i3,l)') 'In collect_state_pdaf, check!',state_p(i:i+2),mype_model,task_id,filterpe
!  end if

END SUBROUTINE collect_state_pdaf
