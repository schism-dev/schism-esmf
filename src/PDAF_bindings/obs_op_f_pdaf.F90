!$Id: obs_op_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! SCHISM module
  use schism_glbl,only : elnode,i34,nvrt,idry_e,kbp,znl,npa,nsa,ntracers,errmsg,ics !,&
!                       &in_dir,out_dir,len_in_dir,len_out_dir,eta2,tr_nd,uu2,vv2               
  use schism_msgp, only: parallel_abort,myrank
! PDAF module
  use mod_assimilation, only: iep_obs_mod,obstype_mod,arco_obs_mod,obs_coords_p,dim_obs_p
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe,mype_filter
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! Local vars
  integer i,m,nd,k0,ie,k,ibad,istat,itot,iz
  real swild(max(100,nsa+nvrt+12+ntracers)),swild2(nvrt,12) 
  real zrat
  character(len=72) :: fdb
  integer :: lfdb
  real,allocatable :: m_state_p(:) ! PE-local observed state
  real,allocatable :: elev(:),salt(:,:),temp(:,:),uu(:,:),vv(:,:),ww(:,:) !  assign from state_p

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! write(*,*) 'In obs_op_pdaf, check!',mype_world,task_id,filterpe

! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

!   m_state_p = ??

  iz=3 ! ics=1, use zobs
! Assign state_p to local vars
!  Allocate var
  if (allocated(elev)) deallocate(elev)
  if (allocated(temp)) deallocate(temp)
  if (allocated(salt)) deallocate(salt)
  if (allocated(uu)) deallocate(uu)
  if (allocated(vv)) deallocate(vv)
  if (allocated(ww)) deallocate(ww)
  allocate(elev(npa),temp(nvrt,npa),salt(nvrt,npa),uu(nvrt,npa),vv(nvrt,npa),ww(nvrt,npa))
!  Assign state_p to vars 
  itot=0
  do i=1,npa
     itot=itot+1
     elev(i)=state_p(itot)
  enddo !i
  do i=1,npa
     do k=1,nvrt
        itot=itot+1
        temp(k,i)=state_p(itot)
     end do
     do k=1,kbp(i)-1 ! extrapolation
        temp(k,i)=temp(kbp(i),i)
     end do
  end do
  do i=1,npa
     do k=1,nvrt
        itot=itot+1
        salt(k,i)=state_p(itot)
     end do
     do k=1,kbp(i)-1 ! extrapolation
        salt(k,i)=salt(kbp(i),i)
     end do
  end do
  do i=1,npa
     do k=1,nvrt
        itot=itot+1
        uu(k,i)=state_p(itot)
     end do
     do k=1,kbp(i)-1 ! extrapolation
        uu(k,i)=0.d0
     end do
  end do
  do i=1,npa
     do k=1,nvrt
        itot=itot+1
        vv(k,i)=state_p(itot)
     end do
     do k=1,kbp(i)-1 ! extrapolation
        vv(k,i)=0.d0
     end do
  end do
  do i=1,npa
     do k=1,nvrt
        itot=itot+1
        ww(k,i)=state_p(itot)
     end do
     do k=1,kbp(i)-1 ! extrapolation
        ww(k,i)=0.d0
     end do
  end do




! new28: can we collect m_state_p base on model values,  not using values from
!        state_p by obs index, we need to do interpolation, so index method may
!        not good enough!
!        If m_state_p=-999, We can remove obs_p record here?!

! Get process-local observed state
  ALLOCATE(m_state_p(dim_obs_p))
! write(*,*) 'in obs_op_f_pdaf, dim_obs_p', dim_obs_p,mype_world,task_id,filterpe 
! write(*,*) 'obs_coords_p',maxval(obs_coords_p)

  do i=1,dim_obs_p
     ie=iep_obs_mod(i)
     if ((obstype_mod(i).eq.'z').or.(obstype_mod(i).eq.'Z')) then
!       swild2(1,1:i34(ie))=eta2(elnode(1:i34(ie),ie))
        swild2(1,1:i34(ie))=elev(elnode(1:i34(ie),ie))
     elseif ((obstype_mod(i).eq.'t').or.(obstype_mod(i).eq.'T')) then
!       swild2(1:nvrt,1:i34(ie))=tr_nd(1,1:nvrt,elnode(1:i34(ie),ie))
        swild2(1:nvrt,1:i34(ie))=temp(1:nvrt,elnode(1:i34(ie),ie))
     elseif ((obstype_mod(i).eq.'s').or.(obstype_mod(i).eq.'S')) then
!       swild2(1:nvrt,1:i34(ie))=tr_nd(2,1:nvrt,elnode(1:i34(ie),ie))
        swild2(1:nvrt,1:i34(ie))=salt(1:nvrt,elnode(1:i34(ie),ie))
     elseif ((obstype_mod(i).eq.'u').or.(obstype_mod(i).eq.'U')) then
!       swild2(1:nvrt,1:i34(ie))=uu2(1:nvrt,elnode(1:i34(ie),ie))
        swild2(1:nvrt,1:i34(ie))=uu(1:nvrt,elnode(1:i34(ie),ie))
     elseif ((obstype_mod(i).eq.'v').or.(obstype_mod(i).eq.'V')) then
!       swild2(1:nvrt,1:i34(ie))=vv2(1:nvrt,elnode(1:i34(ie),ie))
        swild2(1:nvrt,1:i34(ie))=vv(1:nvrt,elnode(1:i34(ie),ie))
     else
        call parallel_abort('PDAF: unknown DA data input')
     end if
!    Start to interpolation
     if ((obstype_mod(i).eq.'z').or.(obstype_mod(i).eq.'Z')) then
        m_state_p(i)=sum(arco_obs_mod(i,1:i34(ie))*swild2(1,1:i34(ie)))
     else !3D vars
        if(idry_e(ie)==1) then !dry
            m_state_p(i)=-999.d0
        else !wet
            do m=1,i34(ie) !wet nodes
               nd=elnode(m,ie)
               !Vertical interplation
               if (ics==2) iz=5
               if(obs_coords_p(iz,i)<=znl(kbp(nd),nd)) then
                  k0=kbp(nd); zrat=0.d0
               else if(obs_coords_p(iz,i)>=znl(nvrt,nd)) then
                  k0=nvrt-1; zrat=1.d0
               else
                  k0=0
                  do k=kbp(nd),nvrt-1
                     if(obs_coords_p(iz,i)>=znl(k,nd).and.obs_coords_p(iz,i)<=znl(k+1,nd)) then
                       k0=k
                       zrat=(obs_coords_p(iz,i)-znl(k,nd))/(znl(k+1,nd)-znl(k,nd))
                       exit
                     endif
                  enddo !k
                  if(k0==0) then
                    write(errmsg,*)'PDAF: DA data depth error',i,obs_coords_p(iz,i)
                    call parallel_abort(errmsg)
                  endif
               endif !obs_coords_p(iz,i)
               swild(m)=swild2(k0,m)*(1.d0-zrat)+swild2(k0+1,m)*zrat
            enddo !m

            !Horizonal interplation
            m_state_p(i)=sum(arco_obs_mod(i,1:i34(ie))*swild(1:i34(ie)))
!           write(*,*) 'Vert itp check',i,m_state_p(i),arco_obs_mod(i,1:i34(ie)),swild(1:i34(ie)),task_id,filterpe
!           write(*,*) 'Vert itp check',i,m_state_p(i),task_id,filterpe,mype_filter
        end if
     end if !itp if
  end do   

! Check -999 for m_state_p & delete obs_p record if neceressary
! Do this later
  ibad=0
  do i=1,dim_obs_p
     if (m_state_p(i)<-900.) ibad=ibad+1
  end do
  if (ibad.ne.0) call parallel_abort('PDAF: Some model values at obs pts are at dry region!')

! Get full observed state vector
  CALL PDAF_gather_obs_f(m_state_p, m_state_f, istat)

! Clean up
  DEALLOCATE(m_state_p)
! Deallocate var
  deallocate(elev)
  deallocate(temp)
  deallocate(salt)
  deallocate(uu)
  deallocate(vv)
  deallocate(ww)



END SUBROUTINE obs_op_f_pdaf
