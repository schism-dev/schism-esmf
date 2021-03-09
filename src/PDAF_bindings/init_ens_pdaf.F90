!$Id: init_ens_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens_pdaf --- Initialize ensemble for filter
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! If only a single filter algorithm is used, the 
! ensemble initialization can be performed directly
! in this routine. If a single filter is implemented,
! one can perform the initialization directly here.
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the particular
! ensemble initialization routine for the selected filter.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Check only
  use mod_parallel_pdaf, only: mype_filter,COMM_filter,mype_model!,task_id,filterpe
! new28
! The following is only use for lock-exchange experiment, delete them after done
  use mod_assimilation, only: offset_field_p,varscale,ens_init
  use schism_glbl, only: npa,nvrt,np_global,ipgl,errmsg
  use schism_msgp, only: myrank,rtype,ierr,parallel_abort
  IMPLICIT NONE

! include 'mpif.h' ! for MPI

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! Local variables
  INTEGER ::  member,i,j,k,kk,col,row,is,ie  ! Counters
  real :: increment_init,ic
! Read eof vars
  integer :: rank,dim_state,rcheck,offset(6),offset_p(6),ir,idp,idg
  real,allocatable :: svals(:),meanstate(:)
  real,allocatable :: eof(:),eof_p(:,:),omega(:,:),state_p2(:),ens(:)
  real :: fac

! !CALLING SEQUENCE:
! Called by: PDAF_init       (as U_init_ens)
! Calls: init_seik
! Calls: init_seek
! Calls: init_enkf
!EOP

! write(*,*) 'In init_ens_pdaf, check!',mype_world,task_id,filterpe
! varscale=0.2

! Global offset
  dim_state=np_global*(1+5*nvrt)
  offset(1) = 0  ! elev
  offset(2) = np_global ! temp
  offset(3) = np_global*(1+nvrt) ! salt
  offset(4) = np_global*(1+2*nvrt) ! uu2
  offset(5) = np_global*(1+3*nvrt) ! vv2
  offset(6) = np_global*(1+4*nvrt) ! ww2
! local offset, because we only have 5 fields in offset_field_p
  offset_p(1) = offset_field_p(1) !elev
  offset_p(2) = offset_field_p(2) !tr_nd(1),temp
  offset_p(3) = offset_field_p(2)+npa*nvrt !tr_nd(2), salt
  offset_p(4) = offset_field_p(3) !uu
  offset_p(5) = offset_field_p(4) !vv
  offset_p(6) = offset_field_p(5) !ww

! Do following with ens_init = 1 or 2
  if (ens_init<3) then

  open(10,file='eofs.dat')
  read(10,*) rank
  close(10)
  allocate(svals(rank),eof(dim_state)) !,rank))
  allocate(eof_p(dim_p,rank),omega(dim_ens, dim_ens-1))
  if (ens_init==2) allocate(state_p2(dim_p),meanstate(dim_state))
! open analysis binary file
! do this read file in all filter rank
! IF (mype_filter==0) THEN
     open(12,file='eofs.bin',form='unformatted')
     open(13,file='svd.bin',form='unformatted')
     read(13) rcheck
     if (rcheck.ne.rank) call parallel_abort('Please check svd & eofs binary inputs! ') ! check rank
     read(13) svals
     close(13)
     if (ens_init==2) then
        open(11,file='meanstate.bin',form='unformatted')
        read(11) meanstate
        close(11)
     end if
!    Check
!    IF (mype_filter==0) THEN
!        write(*,*) 'meanstate T666=',meanstate(np_global+665*nvrt+1:np_global+666*nvrt)
!        write(*,*) 'meanstate T2=',meanstate(np_global+nvrt+1:np_global+2*nvrt)
!    end if
!    read(12) eof
!    close(12)
! IF (mype_filter==0) THEN
!    WRITE (*,'(a,8x,a)') 'SCHISM-PDAF','--- Done reading eofs and svd'
! end if
! Bcast to other filter member
! CALL MPI_Bcast(svals, rank, rtype, 0, COMM_filter, ierr)
! CALL MPI_Bcast(meanstate, dim_state, rtype, 0, COMM_filter, ierr)
! CALL MPI_Bcast(eof, dim_state*rank, rtype, 0, COMM_filter, ierr)

! Start to assign local from global values for each filter rank     
! ic=0 !check ipgl
  do ir=1,rank ! read eof for each rank to reduce memory usage
     read(12) eof
     do j=1,6 !nfields
        do i=1,np_global
           if(ipgl(i)%rank==myrank) then
!          if(ipgl(i)%rank==mype_filter) then
!          if(ipgl(i)%rank==mype_model) then
!             adding offset
              if (j.eq.1) then
                 kk=1 !elev
              else
                 kk=nvrt
              end if
              do k=1,kk
                 if (j.ge.2) then
!                   idp=ipgl(i)%id+offset_p(j)+(ipgl(i)%id-1)*nvrt+k
                    idp=offset_p(j)+(ipgl(i)%id-1)*nvrt+k
!                   idg=i+offset(j)+(i-1)*nvrt+k
                    idg=offset(j)+(i-1)*nvrt+k
                 else
                    idp=ipgl(i)%id+offset_p(j)
                    idg=i+offset(j)
                 end if
                ! check idx
                 if (myrank.eq.0) then
                    if (idp.gt.dim_p) write(*,*) 'idp',ir,i,j,k,idp
                    if (idg.gt.dim_state) write(*,*) 'idg',ir,i,j,k,idg
                 end if
!                if ((myrank.eq.0).and.(j.eq.2).and.(ir.eq.1)) then
!                    write(*,*) 'rank',ipgl(i)%rank
!                    write(*,*) i,k,idp,idg,npa,np_global,dim_p,dim_state
!                end if
                 if (ens_init==2) state_p2(idp)=meanstate(idg)
                 eof_p(idp,ir)=eof(idg) !,:) !reduce 1 dimension for eof
!                ic=ic+1
              end do
           end if
        end do!i
     end do !j
  end do !ir
  close(12)
  IF (mype_filter==0) THEN
     WRITE (*,'(a,8x,a)') 'SCHISM-PDAF','--- Done reading eofs and svd'
!        write(*,*) 'state T1410=',state_p2(npa+1409*nvrt+1:npa+1410*nvrt)
!        write(*,*) 'state T2=',state_p2(npa+nvrt+1:npa+2*nvrt)
!        write(*,*) 'state elev=',state_p2(npa-10:npa)
  end if
! check ic identical to dim_p
!  write(*,*) 'check count:',myrank,ic,dim_p
! check eof_p
! if (myrank==0) then
!    write(*,*) 'eof:',np_global,eof(:,1)
!    write(*,*) 'eof_p:',dim_p,eof_p(:,1)
! end if

! *** Read ensemble 
  CALL collect_state_pdaf(dim_p, state_p)


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

  IF (dim_ens>1) THEN
     ! Only initialize Omega if ensemble size > 0

     IF (mype_filter==0) THEN

        WRITE (*,'(a,8x,a)') 'SCHISM-PDAF','--- generate state ensemble'

        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_omega(dim_ens-1, Omega, 1, 1)

        ! ***      Generate ensemble of states         ***
        ! *** x_i = x + sqrt(FAC) eofV (Omega C^(-1))t ***

        ! A = Omega C^(-1)
        DO col = 1, dim_ens-1
           DO row = 1, dim_ens
              Omega(row, col) = Omega(row,col) * svals(col)
           END DO
        END DO
     END IF
     CALL MPI_Bcast(Omega, dim_ens*(dim_ens-1), rtype, 0, COMM_filter, ierr)
  END IF

  ! state_ens = state+ sqrt(dim_ens-1) eofV A^T

  DO col = 1,dim_ens
!    new28
     if (ens_init==1) then
        ens_p(1:dim_p,col) = state_p(1:dim_p) ! from hotstart
     else ! ens_init=2
        ens_p(1:dim_p,col) = state_p2(1:dim_p) ! from meanstate.bin
     end if
  END DO

  IF (dim_ens>1) THEN
     ! Only add perturbations if ensemble size > 0

     fac = varscale * SQRT(REAL(dim_ens-1))

     CALL DGEMM('n', 't', dim_p, dim_ens, dim_ens-1, &
          fac, eof_p, dim_p, Omega, dim_ens, 1.0, ens_p, dim_p)
  END IF

! Adding meanstate as initial error
! DO col = 1,dim_ens
!    ens_p(1:dim_p,col) = ens_p(1:dim_p,col) + state_p2(1:dim_p)
! END DO
! is=offset_field_p(2)+1
! ie=offset_field_p(2)+nvrt*npa
! ens_p(is:ie,:)=ens_p(is:ie,:)-1.0d0 ! shift 1 deg as initial error
! is=offset_field_p(2)+nvrt*npa+1
! ie=offset_field_p(2)+nvrt*npa*2
! ens_p(is:ie,:)=ens_p(is:ie,:)-0.5d0 ! shift 0.5 psu as initial error

! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eof_p, omega,eof)
  if (ens_init==2) DEALLOCATE(state_p2,meanstate)

  else ! ens_init=3
      open(10,file='ens.dat')
      read(10,*) rank
      close(10)
      if (rank.ne.dim_ens) then
          write(errmsg,*) 'Please specify right ens.dat & ens.bin, now is ', rank
          call parallel_abort(errmsg)
      end if
      open(12,file='ens.bin',form='unformatted')
      allocate(ens(dim_state))
      do ir=1,rank ! read eof for each rank to reduce memory usage
         read(12) ens
         do j=1,6 !nfields
            do i=1,np_global
               if(ipgl(i)%rank==myrank) then
                  if (j.eq.1) then
                     kk=1 !elev
                  else
                     kk=nvrt
                  end if
                  do k=1,kk
                     if (j.ge.2) then
                        idp=offset_p(j)+(ipgl(i)%id-1)*nvrt+k
                        idg=offset(j)+(i-1)*nvrt+k
                     else
                        idp=ipgl(i)%id+offset_p(j)
                        idg=i+offset(j)
                     end if
                     ! check idx
                     if (myrank.eq.0) then
                        if (idp.gt.dim_p) write(*,*) 'idp',ir,i,j,k,idp
                        if (idg.gt.dim_state) write(*,*) 'idg',ir,i,j,k,idg
                     end if
                     ens_p(idp,ir)=ens(idg) !assign ens_p 
                  end do !k
               end if !ipgl
            end do!i
         end do !j
      end do !ir
      close(12)
      IF (mype_filter==0) THEN
         WRITE (*,'(a,8x,a)') 'SCHISM-PDAF','--- Done reading restart ens '
      end if
      DEALLOCATE(ens)

  end if !ens_init


! is=offset_field_p(2)+1
! ie=offset_field_p(2)+nvrt*npa
! Debug
! write(*,'(a,6f8.2,i4)') 'state_p in ens',state_p(is-1:is+1),state_p(ie-1:ie+1),mype_model
! write(*,'(a,f6.2,2i4,l2)') 'In init_ens_pdaf, state_p(max)',maxval(state_p),kind(state_p),mype_world,task_id,filterpe
! *** Initialize ens_p
! state_p(is:ie)=state_p(is:ie)-0.25d0 !this is just for lock-exchange test, will remove after test done!
! increment_init= 0.5d0/real(dim_ens-1,8)
! DO member = 1,dim_ens
!       ens_p(:,member) = state_p(:) !+ float(member-1)*1.d-2 ! temporary, add fesom example later
!       state_p(is:ie)=state_p(is:ie)+increment_init !this is just for lock-exchange test, will remove after test done!
!       ens_p(:,member) = state_p(:)  
! END DO
! write(*,*) 'in init_ens_pdaf',state_p(1:3)
! FESOM example is followed by init_seik, need to gen_cov first, then use its
! output to add pertubation into ensembles



! Disable, just collect
! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

! IF (filtertype == 0) THEN
!    ! EOF initialization for SEEK
!    CALL init_seek(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! ELSE IF (filtertype == 2) THEN
!    ! Use random sampling initialization
!    CALL init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! ELSE
!    ! Use 2nd-order exact sampling
!    CALL init_seik(filtertype, dim_p, dim_ens, state_p, Uinv, &
!         ens_p, flag)
! END IF

END SUBROUTINE init_ens_pdaf
