!$Id: init_dim_obs_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_all(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! SCHISM module
  use schism_glbl, only: nea,ne,ics,dp,elnode,rearth_eq,rearth_pole,xctr,yctr,zctr,eframe,i34,small2,pi,idry_e,rkind,nsteps_from_cold,dt,xnd,ynd,xlon,ylat
  use schism_msgp, only: parallel_abort
! PDAF user define
! new28 add in mod_assimilation, add in some schism_interpolation required here
  use mod_assimilation, only: obs_p,iep_obs_mod,obstype_mod,arco_obs_mod,obs_coords_p,rms_type,rms_obs,rms_obs2,Zdepth_limit,dim_obstype_mod,obstype_index_mod,min_MSL_acDay
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! Local vars
  character(len=17) fnDA
  character(len=1), allocatable :: obstype(:) ! z/s/t/u/v
  character(len=1), allocatable :: obstype2(:) ! z/s/t/u/v
  real(rkind), allocatable :: xobs(:),yobs(:),zobs(:),zzobs(:),obsval(:),iep_obs(:),arco_obs(:,:)
  real(rkind), allocatable :: xobs2(:),yobs2(:),zobs2(:),zzobs2(:),obsval2(:),rmsval2(:)
  integer, allocatable :: idx_select(:)
  real(rkind), allocatable :: rmsval(:)
  integer nobs,i,l,itmp,ifl,iobs,istat,j,nd,nobs2,ic
  real(rkind) tmp,xtmp,ytmp,xobsl,yobsl,zcomp
  real(rkind) xndmax,xndmin,yndmax,yndmin
  logical fexist
  integer max_index_obstype,id_type(6) !counter for each type obs

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
! Called by: PDAF_netf_analysis
!EOP

! write(*,*) 'In init_dim_obs_pdaf, check!',mype_world,task_id,filterpe

! ******************************************************************
! *** Initialize observation dimension for PE-local model domain ***
! ******************************************************************


! Data inputs
! Specify i8 as steps
  write(fnDA,'(a5,i8.8,a4)') 'data_',step,'.dat'

  inquire(file='./DA_data/'//fnDA,exist=fexist)
  if (fexist) then ! file exist
     open(31,file='./DA_data/'//fnDA,status='old')
     read(31,*) nobs
  else
     nobs=0
  end if
  allocate(xobs(nobs),yobs(nobs),zobs(nobs),obsval(nobs),obstype(nobs),&
          &iep_obs(nobs),arco_obs(nobs,4),rmsval(nobs),idx_select(nobs),stat=istat)
  if (ics==2) allocate(zzobs(nobs))
  if(istat/=0) call parallel_abort('PDAF: observation allocation failure')

  do i=1,nobs
     if ((rms_type==1).or.(rms_type==2)) then
         read(31,*) obstype(i),xobs(i),yobs(i),zobs(i),obsval(i) 
     else if (rms_type==3) then
         read(31,*) obstype(i),xobs(i),yobs(i),zobs(i),obsval(i),rmsval(i)
     end if
     if (rms_type==2) then
        if ((obstype(i).eq.'a').or.(obstype(i).eq.'A')) rmsval(i)=rms_obs2(1) !SSH-A
        if ((obstype(i).eq.'z').or.(obstype(i).eq.'Z')) rmsval(i)=rms_obs2(1) !Z
        if ((obstype(i).eq.'t').or.(obstype(i).eq.'T')) rmsval(i)=rms_obs2(2) !T
        if ((obstype(i).eq.'s').or.(obstype(i).eq.'S')) rmsval(i)=rms_obs2(3) !S
        if ((obstype(i).eq.'u').or.(obstype(i).eq.'U')) rmsval(i)=rms_obs2(4) !U
        if ((obstype(i).eq.'v').or.(obstype(i).eq.'V')) rmsval(i)=rms_obs2(5) !V
     end if
     if (rms_type==1) rmsval(i)=rms_obs
     zobs(i)=0.-zobs(i) !negtive, Input zobs is Positive, this is for znl itp
     if(ics==2) then
        zzobs(i)=zobs(i) ! save for skip code
        xtmp=xobs(i)/180.d0*pi
        ytmp=yobs(i)/180.d0*pi
        xobs(i)=rearth_eq*cos(ytmp)*cos(xtmp)
        yobs(i)=rearth_eq*cos(ytmp)*sin(xtmp)
        zobs(i)=rearth_pole*sin(ytmp)
     endif !ics
  end do
  if (fexist) close(31) ! file exist

  !Limit obs in subdomain to speed-up searching
! if (ics.eq.2) then !xnd/ynd are already converted
!    xndmin=minval(xlon(:))
!    xndmax=maxval(xlon(:))
!    yndmin=minval(ylat(:))
!    yndmax=maxval(ylat(:))
! else
     xndmin=minval(xnd(:))
     xndmax=maxval(xnd(:))
     yndmin=minval(ynd(:))
     yndmax=maxval(ynd(:))
! end if
  nobs2=0
  idx_select=0
  do i=1,nobs
     if ((xobs(i).gt.xndmin).and.(xobs(i).lt.xndmax).and.(yobs(i).gt.yndmin).and.(yobs(i).lt.yndmax)) then
        nobs2=nobs2+1
        idx_select(i)=1
     end if
  end do
  allocate(xobs2(nobs2),yobs2(nobs2),zobs2(nobs2),obsval2(nobs2),obstype2(nobs2),rmsval2(nobs2))
  if (ics==2) allocate(zzobs2(nobs2))
  ic=0
  do i=1,nobs
     if (idx_select(i).eq.1) then
         ic=ic+1
         xobs2(ic)=xobs(i)
         yobs2(ic)=yobs(i)
         zobs2(ic)=zobs(i)
         obsval2(ic)=obsval(i)
         obstype2(ic)=obstype(i)
         rmsval2(ic)=rmsval(i)
         if (ics==2) zzobs2(ic)=zzobs(i)
      end if
  end do
  !Swap to original vars
!  write(*,*) 'Check obs dim,',mype_world,nobs,nobs2
   nobs=nobs2
   xobs(1:nobs)=xobs2(:)
   yobs(1:nobs)=yobs2(:)
   zobs(1:nobs)=zobs2(:)
   obsval(1:nobs)=obsval2(:)
   rmsval(1:nobs)=rmsval2(:)
   obstype(1:nobs)=obstype2(:)
   if (ics==2) zzobs(1:nobs)=zzobs2(:)
   !Deallocate
   deallocate(xobs2,yobs2,zobs2,obsval2,obstype2,rmsval2)
   if (ics==2) deallocate(zzobs2)
   !End of limit obs
   
! do i=1,ne
!    do j=1,i34(i)
!       if ((i.eq.1).and.(j.eq.1)) then ! initialize
!          xndmax=xnd(elnode(j,i))
!          xndmin=xnd(elnode(j,i))
!          yndmax=ynd(elnode(j,i))
!          yndmin=ynd(elnode(j,i))
!       end if
!       if (xnd(elnode(j,i)).gt.xndmax) xndmax=xnd(elnode(j,i))
!       if (xnd(elnode(j,i)).lt.xndmin) xndmin=xnd(elnode(j,i))
!       if (ynd(elnode(j,i)).gt.yndmax) yndmax=ynd(elnode(j,i))
!       if (ynd(elnode(j,i)).lt.yndmin) yndmin=ynd(elnode(j,i))
!    end do
! end do
  !Pre-select to speed up
! do l=1,nobs
!     xobsl=xobs(l)
!     yobsl=yobs(l)
!     if ((xobsl.gt.xndmin).and.(xobsl.lt.xndmax).and.&
!        &(yobsl.gt.yndmin).and.(yobsl.lt.yndmax)) iep_obs(l)=1
! end do
     
! Find parent elements in argumented
  iep_obs=0 !flag for no-parent
  do i=1,ne ! search in resident domain to avoid overlap use of observations
     if(idry_e(i)==1) cycle ! skip dry points 
     do l=1,nobs
          if(iep_obs(l)/=0) cycle
!         if(iep_obs(l)==0) cycle ! skip to speedup searching

          if(ics==1) then
             xobsl=xobs(l)
             yobsl=yobs(l)
          else !to eframe
             call project_pt('g2l',xobs(l),yobs(l),zobs(l), &
             &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xobsl,yobsl,tmp)
          endif

          if(i34(i)==3) then
             call area_coord(0,i,(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xobsl,yobsl,arco_obs(l,1:3))
             tmp=minval(arco_obs(l,1:3))
             if(tmp>-small2) then
                iep_obs(l)=i
                if(tmp<0) call area_coord(1,i,(/xctr(i),yctr(i),zctr(i)/), &
                  &eframe(:,:,i),xobsl,yobsl,arco_obs(l,1:3)) !fix A.C.
                !skip data > dp, any node depth in this element > dp, then skip it
                 if (ics==1) then
                     zcomp=abs(zobs(l))
                 else
                     zcomp=abs(zzobs(l))
                 end if
                 do j=1,i34(i)
                    nd=elnode(j,i)
                    if (zcomp.gt.dp(nd)) iep_obs(l)=0
                 end do
                !skip SSH & SSH-A data if it locates in depth < Zdepth_limit (on shelf)
                 if ((obstype(l).eq.'a').or.(obstype(l).eq.'A').or. &
                    &(obstype(l).eq.'z').or.(obstype(l).eq.'Z')) then
                    do j=1,i34(i)
                       nd=elnode(j,i)
                       if (dp(nd).lt.Zdepth_limit) iep_obs(l)=0
                    end do
                 end if         
             endif
          else !quad
             call quad_shape(0,0,i,xobsl,yobsl,itmp,arco_obs(l,1:4)) !arco_sta are 4 shape functions
             if(itmp/=0) iep_obs(l)=i
               !skip data > dp, any node depth in this element > dp, then skip it
               if (ics==1) then
                   zcomp=abs(zobs(l))
               else
                   zcomp=abs(zzobs(l))
               end if
               do j=1,i34(i)
                  nd=elnode(j,i)
                  if (zcomp.gt.dp(nd)) iep_obs(l)=0
               end do
               !skip SSH & SSH-A data if it locates in depth < Zdepth_limit (on shelf)
               if ((obstype(l).eq.'a').or.(obstype(l).eq.'A').or. &
                  &(obstype(l).eq.'z').or.(obstype(l).eq.'Z')) then
                  do j=1,i34(i)
                     nd=elnode(j,i)
                     if (dp(nd).lt.Zdepth_limit) iep_obs(l)=0
                  end do
               end if         
          endif !i34
          !Check nsteps_from_cold for SSH-A, skip if nsteps_from_cold too small
          if ((obstype(l).eq.'a').or.(obstype(l).eq.'A')) then
             if ((nsteps_from_cold*dt/86400.d0).lt.min_MSL_acDay) iep_obs(l)=0 ! set 10 days
          end if
          
     enddo !l; build pts

     ifl=0 !flag
     do l=1,nobs
        if(iep_obs(l)==0) then
           ifl=1
           exit
        endif
     end do !l
     if(ifl==0) exit
   enddo !i=1,ne

!  Allocate dim_obstype_mod, currently sepcify with 6 (a/z/s/t/u/v)
   if (allocated(dim_obstype_mod)) deallocate(dim_obstype_mod)
   allocate(dim_obstype_mod(6))
   dim_obstype_mod=0
!  Count obs points in each sub domain
   iobs=0
   do l=1,nobs
      if (iep_obs(l).ne.0) then
         iobs=iobs+1
         if ((obstype(l).eq.'a').or.(obstype(l).eq.'A')) dim_obstype_mod(1)=dim_obstype_mod(1)+1
         if ((obstype(l).eq.'z').or.(obstype(l).eq.'Z')) dim_obstype_mod(2)=dim_obstype_mod(2)+1
         if ((obstype(l).eq.'s').or.(obstype(l).eq.'S')) dim_obstype_mod(3)=dim_obstype_mod(3)+1
         if ((obstype(l).eq.'t').or.(obstype(l).eq.'T')) dim_obstype_mod(4)=dim_obstype_mod(4)+1
         if ((obstype(l).eq.'u').or.(obstype(l).eq.'U')) dim_obstype_mod(5)=dim_obstype_mod(5)+1
         if ((obstype(l).eq.'v').or.(obstype(l).eq.'V')) dim_obstype_mod(6)=dim_obstype_mod(6)+1
      end if
   end do
   dim_obs_p=iobs
!  allocate obstype_index_mod
   if (allocated(obstype_index_mod)) deallocate(obstype_index_mod)
   max_index_obstype=maxval(dim_obstype_mod)
   allocate(obstype_index_mod(6,max_index_obstype))

!  Allocate mod assimilation vars (obs_p,iep_obs_mod,obstype_mod,arco_obs_mod) 
   if (allocated(obs_p)) deallocate(obs_p)
   if (allocated(iep_obs_mod)) deallocate(iep_obs_mod)
   if (allocated(obstype_mod)) deallocate(obstype_mod)
   if (allocated(arco_obs_mod)) deallocate(arco_obs_mod)
   if (allocated(obs_coords_p)) deallocate(obs_coords_p)
   allocate(obs_p(iobs))
   allocate(iep_obs_mod(iobs))
   allocate(obstype_mod(iobs))
   allocate(arco_obs_mod(iobs,4))
   allocate(obs_coords_p(5,iobs)) !5 dim, store x,y,z,rms_obs_vec,zzobs(ics=2)
   
!  We can Put simple QA/QC here, do this later
!   call qaqc(obsval,obstype) 

!  Fill assimilation vars with local, and pass them to obs_op to do vertical interpolation
   iobs=0
   id_type=0
   do l=1,nobs
      if (iep_obs(l).ne.0) then
         iobs=iobs+1
         obs_p(iobs)=obsval(l)
         iep_obs_mod(iobs)=iep_obs(l)
         obstype_mod(iobs)=obstype(l)
         arco_obs_mod(iobs,:)=arco_obs(l,:)
         obs_coords_p(1,iobs)=xobs(l)
         obs_coords_p(2,iobs)=yobs(l)
         obs_coords_p(3,iobs)=zobs(l)
         obs_coords_p(4,iobs)=rmsval(l)
         if (ics==2) obs_coords_p(5,iobs)=zzobs(l)
!        write(*,*) 'check arco_obs_mod',arco_obs_mod(iobs,:),arco_obs(l,:),iep_obs(l)
         ! Store corresponding index for each type obs
         if ((obstype(l).eq.'a').or.(obstype(l).eq.'A')) then
             id_type(1)=id_type(1)+1
             obstype_index_mod(1,id_type(1))=iobs
         end if
         if ((obstype(l).eq.'z').or.(obstype(l).eq.'Z')) then
             id_type(2)=id_type(2)+1
             obstype_index_mod(2,id_type(2))=iobs
         end if
         if ((obstype(l).eq.'s').or.(obstype(l).eq.'S')) then
             id_type(3)=id_type(3)+1
             obstype_index_mod(3,id_type(3))=iobs
         end if
         if ((obstype(l).eq.'t').or.(obstype(l).eq.'T')) then
             id_type(4)=id_type(4)+1
             obstype_index_mod(4,id_type(4))=iobs
         end if
         if ((obstype(l).eq.'u').or.(obstype(l).eq.'U')) then
             id_type(5)=id_type(5)+1
             obstype_index_mod(5,id_type(5))=iobs
         end if
         if ((obstype(l).eq.'v').or.(obstype(l).eq.'V')) then
             id_type(6)=id_type(6)+1
             obstype_index_mod(6,id_type(6))=iobs
         end if
      end if
   end do

!  Clean-up
   deallocate(xobs,yobs,zobs,obsval,obstype,iep_obs,arco_obs,rmsval)
   if (ics==2) deallocate(zzobs)

END SUBROUTINE init_dim_obs_all

