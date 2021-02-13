!$Id: init_dim_obs_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

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
  use schism_glbl, only: nea,ics,dp,elnode,rearth_eq,rearth_pole,xctr,yctr,zctr,eframe,i34,small2,pi,idry_e,rkind
  use schism_msgp, only: parallel_abort
! PDAF user define
! new28 add in mod_assimilation, add in some schism_interpolation required here
  use mod_assimilation, only: obs_p,iep_obs_mod,obstype_mod,arco_obs_mod,obs_coords_p
! Check only
  use mod_parallel_pdaf, only: mype_world,task_id,filterpe


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! Local vars
  character(len=17) fnDA
  character(len=1), allocatable :: obstype(:) ! z/s/t/u/v
  real(rkind), allocatable :: xobs(:),yobs(:),zobs(:),zzobs(:),obsval(:),iep_obs(:),arco_obs(:,:)
  integer nobs,i,l,itmp,ifl,iobs,istat,j,nd
  real(rkind) tmp,xtmp,ytmp,xobsl,yobsl,zcomp
  logical fexist

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
          &iep_obs(nobs),arco_obs(nobs,4),stat=istat)
  if (ics==2) allocate(zzobs(nobs))
  if(istat/=0) call parallel_abort('PDAF: observation allocation failure')

  do i=1,nobs
     read(31,*) obstype(i),xobs(i),yobs(i),zobs(i),obsval(i) 
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

! Find parent elements in argumented
  iep_obs=0 !flag for no-parent
  do i=1,nea ! search in argumented
     if(idry_e(i)==1) cycle ! skip dry points 
     do l=1,nobs
          if(iep_obs(l)/=0) cycle

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
          endif !i34
      enddo !l; build pts

      ifl=0 !flag
      do l=1,nobs
         if(iep_obs(l)==0) then
            ifl=1
            exit
         endif
      end do !l
      if(ifl==0) exit
   enddo !i=1,nea

!  Count obs points in each sub domain
   iobs=0
   do l=1,nobs
      if (iep_obs(l).ne.0) iobs=iobs+1
   end do
   dim_obs_p=iobs

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
   allocate(obs_coords_p(3,iobs)) !3 dim, store x,y,z
   
!  We can Put simple QA/QC here, do this later
!   call qaqc(obsval,obstype) 

!  Fill assimilation vars with local, and pass them to obs_op to do vertical interpolation
!  new28:  if model shows -9999, can we omit this values in obs_op?
   iobs=0
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
!        write(*,*) 'check arco_obs_mod',arco_obs_mod(iobs,:),arco_obs(l,:),iep_obs(l)
      end if
   end do

!  Clean-up
   deallocate(xobs,yobs,zobs,zzobs,obsval,obstype,iep_obs,arco_obs)

END SUBROUTINE init_dim_obs_pdaf

