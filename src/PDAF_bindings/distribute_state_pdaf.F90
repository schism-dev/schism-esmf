!$Id: distribute_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use schism_glbl, only: nea,nsa,npa,nvrt,ntracers,idry_e,we,tr_el, &
 &idry_s,su2,sv2, idry,eta2,tr_nd,tr_nd0,q2,xl,dfv,dfh,dfq1,dfq2
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

  integer :: i,j,k,itot

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_dist_state)
! Called by: PDAF_assimilate_X   (as U_dist_state)
!EOP


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  ! Template reminder - delete when implementing functionality
  !WRITE (*,*) 'TEMPLATE distribute_state_pdaf.F90: Implement initialization of model fields here!'

!  Assign state_p to eta2, su2 etc; exchange ghost
!  ? = state_p
   
   !Elem
   itot=0
   do i=1,nea
     itot=itot+1
     idry_e(i)=nint(state_p(itot))
   enddo !i
   do i=1,nea
     do k=1,nvrt
       itot=itot+1
       we(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,nea
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         tr_el(j,k,i)=state_p(itot)
       enddo !j
     enddo !k
   enddo !i
  
   !Side
   do i=1,nsa
     itot=itot+1
     idry_s(i)=nint(state_p(itot))
   enddo !i
   do i=1,nsa
     do k=1,nvrt
       itot=itot+1
       su2(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,nsa
     do k=1,nvrt
       itot=itot+1
       sv2(k,i)=state_p(itot)
     enddo !k
   enddo !i

   !node
   do i=1,npa
     itot=itot+1
     idry(i)=nint(state_p(itot))
   enddo !i
   do i=1,npa
     itot=itot+1
     eta2(i)=state_p(itot)
   enddo !i
   do i=1,npa
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         tr_nd(j,k,i)=state_p(itot)
       enddo !j
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         tr_nd0(j,k,i)=state_p(itot)
       enddo !j
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       q2(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       xl(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       dfv(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       dfh(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       dfq1(k,i)=state_p(itot)
     enddo !k
   enddo !i
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       dfq2(k,i)=state_p(itot)
     enddo !k
   enddo !i

END SUBROUTINE distribute_state_pdaf
