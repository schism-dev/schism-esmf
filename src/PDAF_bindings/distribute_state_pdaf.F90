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
    & idry_s,su2,sv2,idry,eta2,tr_nd,uu2,vv2,ww2, &
    & elnode,i34,rkind,kbe,kbs,isidenode
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe

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
! write(*,*) 'In distribute_state_pdaf, check!',mype_world,task_id,filterpe
   
   !init count
   itot=0
   !node-base
   ! Elev
   do i=1,npa
     itot=itot+1
     eta2(i)=state_p(itot)
   enddo !i
   ! Tracer
   do i=1,npa
     do k=1,nvrt
       do j=1,ntracers
         itot=itot+1
         tr_nd(j,k,i)=state_p(itot)
       enddo !j
     enddo !k
   enddo !i
   ! U
   do i=1,npa
     do k=1,nvrt
        itot=itot+1
        uu2(k,i)=state_p(itot)
     enddo !k
   enddo !i
   ! V
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       vv2(k,i)=state_p(itot)
     enddo !k
   enddo !i
   ! W
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       ww2(k,i)=state_p(itot)
     enddo !k
   enddo !i

!  new28!!
!  Here we need to update tr_el, su2,sv2,we by node-base vars, do this later
!  Update tr_el
   do i=1,nea    
      if (idry_e(i).eq.1) cycle !dry
      do j=1,ntracers
         do k=kbe(i),nvrt
            tr_el(j,k,i)=sum(tr_nd(j,k,elnode(1:i34(i),i)))/real(i34(i),rkind)
         end do !k
         do k=1,kbe(i)-1
            tr_el(j,k,i)=tr_el(j,kbe(i),i) !extrapolation
         end do !k
      end do !j
   end do !i
!  Update su2,sv2
   do j=1,nsa
      if(idry_s(j)==1) cycle
      do k=kbs(j),nvrt
         su2(k,j)=sum(uu2(k,isidenode(:,j)))/2.d0
         sv2(k,j)=sum(vv2(k,isidenode(:,j)))/2.d0
      end do
      do k=1,kbs(j)-1
         su2(k,j)=0.d0  !zero-out
         sv2(k,j)=0.d0  !zero-out
      end do
   end do
!  Update we
   do i=1,nea
      if(idry_e(j)==1) cycle
      do k=kbe(i),nvrt
         we(k,i)=sum(ww2(k,elnode(1:i34(i),i)))/real(i34(i),rkind)
      end do
      do k=1,kbe(i)-1
         we(k,i)=0.d0  !zero-out
      end do
   end do
   

END SUBROUTINE distribute_state_pdaf
