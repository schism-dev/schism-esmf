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
    & elnode,i34,rkind,kbe,kbs,isidenode,kbp
! Check only
  use mod_parallel_pdaf, only: mype_model,task_id,filterpe
  use mod_assimilation, only: offset_field_p

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

  integer :: i,j,k,itot,ifill,ifill2

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
!  Debug
!  i=offset_field_p(2)+1
!  if (mype_model.eq.2) then
!  write(*,'(a,3f8.3,2i3,l)') 'In distribute_state_pdaf, check!',state_p(i:i+2),mype_model,task_id,filterpe
!  end if
   
!  ifill=0, original
!  ifill=1, overite -9999 as bottom value
   ifill=1
   ifill2=1 ! extra-control for side,element base vars
   !init count
   itot=0
   !node-base
   ! Elev
   do i=1,npa
     itot=itot+1
     eta2(i)=state_p(itot)
   enddo !i
   ! Tracer
   do j=1,ntracers
     do i=1,npa
       do k=1,nvrt
         itot=itot+1
         tr_nd(j,k,i)=state_p(itot)
       enddo !k
       if (ifill.eq.1) then
!      Fill -9999 to see difference of analysis
       do k=1,kbp(i)-1
!         write(*,*) 'check if -9999 in tracer',tr_nd(j,k,i)
          tr_nd(j,k,i)=tr_nd(j,kbp(i),i)
       end do
       end if
     enddo !i
   enddo !j
   ! U
   do i=1,npa
     do k=1,nvrt
        itot=itot+1
        uu2(k,i)=state_p(itot)
     enddo !k
     if (ifill.eq.1) then
!    Fill -9999 to see difference of analysis
     do k=1,kbp(i)-1
!       write(*,*) 'check if -9999 in U',uu2(k,i)
        uu2(k,i)= 0.
     end do
     end if
   enddo !i
   ! V
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       vv2(k,i)=state_p(itot)
     enddo !k
     if (ifill.eq.1) then
!    Fill -9999 to see difference of analysis
     do k=1,kbp(i)-1
!       write(*,*) 'check if -9999 in V',vv2(k,i)
        vv2(k,i)= 0.
     end do
     end if
   enddo !i
   ! W
   do i=1,npa
     do k=1,nvrt
       itot=itot+1
       ww2(k,i)=state_p(itot)
     enddo !k
     if (ifill.eq.1) then
!    Fill -9999 to see difference of analysis
     do k=1,kbp(i)-1
!       write(*,*) 'check if -9999 in V',ww2(k,i)
        ww2(k,i)= 0.
     end do
     end if
   enddo !i

!  new28!!
!  Here we need to update tr_el, su2,sv2,we by node-base vars
   if (ifill2.eq.1) then
!  Update tr_el
   do i=1,nea    
      if (idry_e(i).eq.1) cycle !dry

!     Here we need to assign first, then extrapolation 
!     do k=kbe(i),nvrt
!        do j=1,ntracers
!           tr_el(j,k,i)=sum(tr_nd(j,k,elnode(1:i34(i),i)))/real(i34(i),rkind)
!        end do !j  
!     end do !k
      do k=kbe(i)+1,nvrt
         do j=1,ntracers
            tr_el(j,k,i)=sum(tr_nd(j,k,elnode(1:i34(i),i))+tr_nd(j,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
         end do !j  
      end do !k

!     Specify 1st index tr_el
!     do j=1,ntracers
!        tr_el(j,1,i)=tr_el(j,2,i)
!     end do !j
!     Do extrapolation
      do k=1,kbe(i)!-1
         do j=1,ntracers
            tr_el(j,k,i)=tr_el(j,kbe(i),i) !extrapolation
         end do !j
      end do !k
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
      do k=kbe(i)+1,nvrt
         we(k,i)=sum(ww2(k,elnode(1:i34(i),i)))/real(i34(i),rkind)
      end do
!     do k=2,nvrt
!        we(k,i)=sum(ww2(k,elnode(1:i34(i),i))+ww2(k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
!     end do
!     we(1,i)=we(2,i)
      do k=1,kbe(i)!-1
         we(k,i)=0.d0  !zero-out
      end do
   end do
   
   end if !ifill

END SUBROUTINE distribute_state_pdaf
