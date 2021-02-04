!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!===============================================================================
!===============================================================================
! SCHISM FILE I/O SUBROUTINES
!
! subroutine write_netcdf_pdaf
! subroutine writeout_nc_PDAF
! subroutine fill_nc_header_PDAF

!===============================================================================
!===============================================================================

    module output_schism_pdaf
    use schism_glbl, only: nea,nsa,npa,nvrt,idry,idry_e,idry_s,znl,id_out_var,kbp,rkind,&
                   & np,ne,ns,time_stamp
    use schism_msgp, only: myrank,parallel_abort
    use netcdf
    use mod_assimilation, only: ihfskip_PDAF,nspool_PDAF
    implicit none
!    include 'netcdf.inc'
    private

    integer,save :: node_dim,nele_dim,nedge_dim,four_dim,nv_dim, &
    &one_dim,two_dim,time_dim,time_dims(1),itime_id,ele_dims(2),x_dims(1), &
    &y_dims(1),z_dims(1),var2d_dims(2),var3d_dims(3),var4d_dims(4),dummy_dim(1), &
    &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
    &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4)

    character(len=12),save :: ifile_char_PDAF
    
    integer,save :: ncid_PDAF_io,ifile_PDAF,ifile_len_PDAF,ncid_a,ncid_f,ifile_PDAF_a,ifile_PDAF_f
    character(len=1000),save :: out_dir 
    integer,save :: len_out_dir,it_main_PDAF !specify it_main from step
    logical,save :: ifirst

!   public :: writeout_nc_PDAF
!   public :: fill_nc_header_PDAF
    public :: write_netcdf_pdaf

    contains

!===============================================================================
!===============================================================================
!
      subroutine write_netcdf_pdaf( step, dim_p, state_p)
      implicit none
!     character(len=1),intent(in) :: typestr
      integer,intent(in) :: step,dim_p
      real(rkind),intent(in) :: state_p(dim_p)
!     local var
      integer :: itot,i,k
      real,allocatable :: elev(:),salt(:,:),temp(:,:),uu(:,:),vv(:,:),ww(:,:)
      character(len=1) :: typestr

!     Allocate var
      if (allocated(elev)) deallocate(elev)
      if (allocated(temp)) deallocate(temp)
      if (allocated(salt)) deallocate(salt)
      if (allocated(uu)) deallocate(uu)
      if (allocated(vv)) deallocate(vv)
      if (allocated(ww)) deallocate(ww)
      allocate(elev(npa),temp(nvrt,npa),salt(nvrt,npa),uu(nvrt,npa),vv(nvrt,npa),ww(nvrt,npa))

!     Assign state_p to vars for outputs
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

!     Define typestr to avoid other rank has weird value
      if (step==0) typestr='i'
      if (step>0) typestr='a'
      if (step<0) typestr='f'

!     Specify it_main from step to avoid weird error
      it_main_PDAF=abs(step)
      
!     Define out_dir,len_out_dir,ifile_PDAF
      out_dir=adjustl('./DA_output/'//typestr//'/')
      if (myrank==0) write(*,*) 'typestr=',typestr,myrank,step !check
      len_out_dir=len_trim(out_dir)
      if (step.eq.0) then
         if (myrank==0) write(*,*) 'open DA_output file'
         ifile_PDAF=1
         ifile_PDAF_a=1
         ifile_PDAF_f=1
         ! open f & a at time=0
         ifirst=.True. !For 1st time
         out_dir=adjustl('./DA_output/'//'f/')
         len_out_dir=len_trim(out_dir)
         call fill_nc_header_PDAF(0,'f')
         out_dir=adjustl('./DA_output/'//'a/')
         len_out_dir=len_trim(out_dir)
         call fill_nc_header_PDAF(0,'a')
         ifirst=.False. ! after setting ncid_a/f
      end if

!     Output local_global map
!     Copy from outputs?

!     output Hydro only
      if (step.ne.0) then ! only output f/a file when step nq 0
         if(mod(abs(step),nspool_PDAF)==0) then
            if (myrank==0) write(*,*) 'writenc in PDAF with',step, typestr !, out_dir
!            if (typestr=='f') ncid_PDAF_io=ncid_f
!            if (typestr=='a') ncid_PDAF_io=ncid_a
!                  
!   Noted:   idry,idry_e,idry_s,& znl are followed with filter member, usually the 1st
!            we didn't use ens-avg
             call writeout_nc_PDAF(typestr,id_out_var(1),'wetdry_node',1,1,npa,dble(idry))
             call writeout_nc_PDAF(typestr,id_out_var(2),'wetdry_elem',4,1,nea,dble(idry_e))
             call writeout_nc_PDAF(typestr,id_out_var(3),'wetdry_side',7,1,nsa,dble(idry_s))
             !zcor MUST be 1st 3D var output for combine scripts to work!
             call writeout_nc_PDAF(typestr,id_out_var(4),'zcor',2,nvrt,npa,znl(:,:))

!            output ssh,T,S,u,v,w
             call writeout_nc_PDAF(typestr,id_out_var(5),'elev',1,1,npa,elev)
             call writeout_nc_PDAF(typestr,id_out_var(22),'temp',2,nvrt,npa,temp)
             call writeout_nc_PDAF(typestr,id_out_var(23),'salt',2,nvrt,npa,salt)
             call writeout_nc_PDAF(typestr,id_out_var(29),'hvel',2,nvrt,npa,uu,vv)
             call writeout_nc_PDAF(typestr,id_out_var(21),'vertical_velocity',2,nvrt,npa,ww)
         end if


!     Close init file
!     if (step==0) call fill_nc_header_PDAF(1)
!      if (step==0) 
!     Close f & a file at specify step
!     if (step.ne.0) then ! only close after init
         if (mod(abs(step),ihfskip_PDAF)==0) then
            if (typestr=='a') ifile_PDAF_a=ifile_PDAF_a+1 
            if (typestr=='f') ifile_PDAF_f=ifile_PDAF_f+1 
            call fill_nc_header_PDAF(1,typestr)
         end if
      end if

!     Deallocate var
      deallocate(elev)
      deallocate(temp)
      deallocate(salt)
      deallocate(uu)
      deallocate(vv)
      deallocate(ww)

      end subroutine write_netcdf_pdaf

      subroutine writeout_nc_PDAF(typestr,varid,var_nm,i23d,idim1,idim2,outvar1,outvar2)
!-------------------------------------------------------------------------------
!     Netcdf outputs for global arrays. Can be called from any routine, but make sure that
!     the calling routine is called inside the main time loop 
!     exactly ONCE per step! 
!
!     Inputs:
!            var_nm: name of the output variable (to appear in nc file). 
!            i23d: indicates location where outputs are defined. 1:3 - node 2D/3D whole/3D half level
!                  4:6 - elem 2D/3D whole/half levels; 7:9 - side 2D/3D whole/half levels
!            idim1,idim2: dimensions of output array(s) in the driving routine. 
!                         For 2D variables (e.g., bottom
!                         stress), idim1 must be 1; for 3D variables, idim1 must be nvrt.
!                         idim2 must be consistent with the type of output as given by
!                         i23d (e.g., idim2=nea or ne for i23d=4);
!            outvar[1,2](idim1,idim2): output array. outvar2 is optional [for vectors]
!     In/out: varid: 1st call will generate variable ID, which is used later
!-------------------------------------------------------------------------------
      implicit none
   
      character(len=*),intent(in) :: var_nm
      character(len=1), intent(in) :: typestr
      integer,intent(in) :: i23d,idim1,idim2
      real(rkind),intent(in) :: outvar1(idim1,idim2)
      real(rkind),optional,intent(in) :: outvar2(idim1,idim2)
      integer,intent(inout) :: varid
 
      !character(len=3) :: sfix
      character(len=1000) :: var_nm2
      logical :: lex1,lex2
      integer :: i,k,iret,irec,len_var,idim2p,iret2,ivs
      real*4 :: a1d(1) 
      
!     Return if not output step
      if(mod(it_main_PDAF,nspool_PDAF)/=0) return

!     Replace ncid_PDAF_io    
      if (typestr=='f') then
          ncid_PDAF_io=ncid_f
          ifile_PDAF=ifile_PDAF_f
      end if
      if (typestr=='a') then
          ncid_PDAF_io=ncid_a
          ifile_PDAF=ifile_PDAF_a
      end if

      ivs=1
      if(present(outvar2)) ivs=2

      irec=(it_main_PDAF-(ifile_PDAF-1)*ihfskip_PDAF)/nspool_PDAF !time recod #
      if(irec<=0) call parallel_abort('writeout_nc_PDAF: irec<=0')
      var_nm2=var_nm
      var_nm2=adjustl(var_nm2); len_var=len_trim(var_nm2)

      !Define dim/vars
      !nf90_put_var(ncid_PDAF_io,varid,values,start,count,stride)
      !values can be of any type, (optional) start, count, stride are of same dim
      !as values. e.g., to write to 1st entry in an array, start=count=1

      !Dump time
      !Note: using scalar directly won't work; must use array
      a1d(1)=real(time_stamp)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(ncid_PDAF_io,itime_id,a1d,data_start_1d,data_count_1d)

      !Use original dim order in nc
      if(i23d<=3) then !node
        var2d_dims(1)=node_dim
        var3d_dims(2)=node_dim
        var4d_dims(3)=node_dim
        idim2p=np !final output dim
      else if(i23d<=6) then !elem
        var2d_dims(1)=nele_dim
        var3d_dims(2)=nele_dim
        var4d_dims(3)=nele_dim
        idim2p=ne
      else if(i23d<=9) then !side
        var2d_dims(1)=nedge_dim
        var3d_dims(2)=nedge_dim
        var4d_dims(3)=nedge_dim
        idim2p=ns
      else
        call parallel_abort('writeout_nc: unknown i23d')       
      endif

      iret2=nf90_inq_varid(ncid_PDAF_io,var_nm2(1:len_var),i)

      if(mod(i23d-1,3)==0) then !2D var (2D array in nc that has time dim)
        if(iret2/=NF90_NOERR) then !not defined yet
          iret=nf90_redef(ncid_PDAF_io)
          if(ivs==1) then
            var2d_dims(2)=time_dim
            iret=nf90_def_var(ncid_PDAF_io,var_nm2(1:len_var),NF90_FLOAT,var2d_dims,varid)
          else
            var3d_dims(1)=two_dim; var3d_dims(3)=time_dim
            iret=nf90_def_var(ncid_PDAF_io,var_nm2(1:len_var),NF90_FLOAT,var3d_dims,varid)
          endif !ivs
          iret=nf90_put_att(ncid_PDAF_io,varid,'i23d',i23d)
          iret=nf90_put_att(ncid_PDAF_io,varid,'ivs',ivs)
          iret=nf90_enddef(ncid_PDAF_io)
        endif !iret

        if(ivs==1) then
          data_start_2d(1)=1; data_start_2d(2)=irec
          data_count_2d(1)=idim2p; data_count_2d(2)=1
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar1(1,1:idim2p)),data_start_2d,data_count_2d)
        else !vector
          data_start_3d(1)=1; data_start_3d(2)=1; data_start_3d(3)=irec
          data_count_3d(1)=1; data_count_3d(2)=idim2p; data_count_3d(3)=1
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar1(1,1:idim2p)),data_start_3d,data_count_3d)
          data_start_3d(1)=2
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar2(1,1:idim2p)),data_start_3d,data_count_3d)
        endif !ivs
        !write(12,*)'2D:',it_main,varid,var_nm2(1:len_var),iret2
      else !3D
        if(iret2/=NF90_NOERR) then !not defined yet
          iret=nf90_redef(ncid_PDAF_io)
          if(ivs==1) then
            var3d_dims(1)=nv_dim; var3d_dims(3)=time_dim
            iret=nf90_def_var(ncid_PDAF_io,var_nm2(1:len_var),NF90_FLOAT,var3d_dims,varid)
          else
            var4d_dims(1)=two_dim; var4d_dims(2)=nv_dim; var4d_dims(4)=time_dim
            iret=nf90_def_var(ncid_PDAF_io,var_nm2(1:len_var),NF90_FLOAT,var4d_dims,varid)
          endif !ivs
          !write(12,*)'3D def:',var3d_dims,varid,var_nm2(1:len_var),iret
!Add chunking option as well?
          iret=nf90_def_var_deflate(ncid_PDAF_io,varid,0,1,4)
          iret=nf90_put_att(ncid_PDAF_io,varid,'i23d',i23d)
          iret=nf90_put_att(ncid_PDAF_io,varid,'ivs',ivs)
          iret=nf90_enddef(ncid_PDAF_io)
        endif !iret

        if(ivs==1) then
          data_start_3d(1)=1; data_start_3d(2)=1; data_start_3d(3)=irec
          data_count_3d(1)=nvrt; data_count_3d(2)=idim2p; data_count_3d(3)=1
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar1(:,1:idim2p)),data_start_3d,data_count_3d)
        else !vector
          data_start_4d(1:3)=1; data_start_4d(4)=irec
          data_count_4d(1)=1; data_count_4d(2)=nvrt; data_count_4d(3)=idim2p; data_count_4d(4)=1
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar1(:,1:idim2p)),data_start_4d,data_count_4d)
          data_start_4d(1)=2
          iret=nf90_put_var(ncid_PDAF_io,varid,real(outvar2(:,1:idim2p)),data_start_4d,data_count_4d)
        endif !ivs
        !write(12,*)'3D:',it_main,varid,var_nm2(1:len_var),iret2,NF90_NOERR
      endif !2/3D
 
      end subroutine writeout_nc_PDAF

!===============================================================================
      subroutine fill_nc_header_PDAF(iopen,typestr)
!-------------------------------------------------------------------------------
!     Create nc file and define dimension and static info for netcdf output
!     Input: iopen=1: close nc file handle. 0: do not close
!-------------------------------------------------------------------------------
      implicit none

      integer, intent(in) :: iopen
      character(len=1), intent(in) :: typestr
      character(len=140) :: fname
      character(len=4) :: fgb

      integer :: iret

!     !replace ncid_PDAF_io after 2nd time
      if (.not.ifirst) then
          if (typestr=='f') then
              ncid_PDAF_io=ncid_f
              ifile_PDAF=ifile_PDAF_f
          end if
          if (typestr=='a') then
              ncid_PDAF_io=ncid_a
              ifile_PDAF=ifile_PDAF_a
          end if
      end if

      write(ifile_char_PDAF,'(i12)') ifile_PDAF !convert ifile to a string
      ifile_char_PDAF=adjustl(ifile_char_PDAF)  !place blanks at end
      ifile_len_PDAF=len_trim(ifile_char_PDAF)  !length without trailing blanks
      fgb='0000' 
      write(fgb,'(i4.4)') myrank
      fname=out_dir(1:len_out_dir)//('schout_'//fgb//'_'//ifile_char_PDAF(1:ifile_len_PDAF)//'.nc')
!'
      write(*,*) trim(adjustl(fname)),ifile_PDAF !check

      if(iopen==1) iret=nf90_close(ncid_PDAF_io)
      iret=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_PDAF_io)
!     Specify ncid_a/f 
      if (ifirst) then !1st time
          if (typestr=='f') ncid_f=ncid_PDAF_io
          if (typestr=='a') ncid_a=ncid_PDAF_io
      end if
      iret=nf90_def_dim(ncid_PDAF_io,'nSCHISM_hgrid_node',np,node_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'nSCHISM_hgrid_face',ne,nele_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'nSCHISM_hgrid_edge',ns,nedge_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'nMaxSCHISM_hgrid_face_nodes',4, four_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'nSCHISM_vgrid_layers',nvrt,nv_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'one',1,one_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'two',2,two_dim)
      iret=nf90_def_dim(ncid_PDAF_io,'time', NF90_UNLIMITED,time_dim)

      time_dims(1)=time_dim
      iret=nf90_def_var(ncid_PDAF_io,'time',NF90_FLOAT,time_dims,itime_id)
      if(iret.ne.NF90_NOERR) call parallel_abort('fill_nc_header_PDAF: time dim')
!'
      iret=nf90_put_att(ncid_PDAF_io,itime_id,'i23d',0) !set i23d flag
      iret=nf90_enddef(ncid_PDAF_io)

      end subroutine fill_nc_header_PDAF
!===============================================================================
!===============================================================================
! END FILE I/O module
!===============================================================================
!===============================================================================
    end module output_schism_pdaf
