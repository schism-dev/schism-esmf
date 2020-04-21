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
 &idry_s,su2,sv2, idry,eta2,tr_nd,tr_nd0,q2,xl,dfv,dfh,dfq1,dfq2
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

  integer :: i,j,k,itot

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  !WRITE (*,*) 'TEMPLATE collect_state_pdaf.F90: Implement initialization of state vector here!'

!   Assign state vars to state_p in the resident domains

   !Elem
   itot=0
   do i=1,nea
     itot=itot+1
     state_p(itot)=idry_e(i)
   enddo !i
   do i=1,nea
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=we(k,i)
     enddo !k
   enddo !i
   do i=1,nea
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         state_p(itot)=tr_el(j,k,i) 
       enddo !j
     enddo !k
   enddo !i
  
   !Side
   do i=1,nsa
     itot=itot+1
     state_p(itot)=idry_s(i)
   enddo !i
   do i=1,nsa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=su2(k,i)
     enddo !k
   enddo !i
   do i=1,nsa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=sv2(k,i)
     enddo !k
   enddo !i

   !node
   do i=1,npa
     itot=itot+1
     state_p(itot)=idry(i)
   enddo !i
   do i=1,npa
     itot=itot+1
     state_p(itot)=eta2(i)
   enddo !i
   do i=1,npa
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         state_p(itot)=tr_nd(j,k,i)
       enddo !j
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         state_p(itot)=tr_nd0(j,k,i)
       enddo !j
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=q2(k,i)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=xl(k,i)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=dfv(k,i)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=dfh(k,i)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=dfq1(k,i)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       state_p(itot)=dfq2(k,i)
     enddo !k
   enddo !i

END SUBROUTINE collect_state_pdaf
