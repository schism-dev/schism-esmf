!     This file is used to generate PDAF cov for ensemble generation
!
!     Inputs: dim.dat,z/t/s/u/v/w.bin
!     Outputs: eofs.bin, eofs.dat, svd.bin, meanstate.bin
!
!
!     real,allocatable :: elev(:,:),t(:,:,:),s(:,:,:),u(:,:,:),v(:,:,:),w(:,:,:)
      use PDAF 
      real,allocatable :: temp3d(:,:),temp2d(:)
      real,allocatable :: states(:,:),svals(:),svdU(:,:),meanstate(:),run_meanstate(:,:)
      integer np,nvrt,ntimes,offset(6),dim_fields(6),nfields,status
      real stddev(6)
      integer dim_state
! Controls for PDAF_eofcovar
      INTEGER :: remove_mean  ! (1) Let PDAF_eofcovar compute and subtract the mean state;
                              ! (0) mean already removed from trajectory 
      INTEGER :: do_mv        ! (1) to perform multivariate normalization; (0) no normalization


      open(10,file='dim.dat')
      read(10,*) np
      read(10,*) nvrt
      read(10,*) ntimes
      close(10)
! set output rank to reduce output
      ioutrank=ntimes-1

! set some parameters
      nfields=6
      remove_mean=1
      do_mv=1
! Lower limit for eigenvalue
      limit = 1.0e-12
! Half time window for running mean (2*irange+1)
      hwindow = 24



!     ncl will dump reverse dimension, just follow state_p in PDAF
!     allocate(elev(np,ntimes),t(nvrt,np,ntimes),s(nvrt,np,ntimes),u(nvrt,np,ntimes),v(nvrt,np,ntimes),w(nvrt,np,ntimes))
      allocate(temp2d(np),temp3d(nvrt,np))
!     allocate(dump(np,nvrt))

!     allocate state
      dim_state=np*(1+5*nvrt)
      allocate(states(dim_state,ntimes))

!     set offset
      offset(1) = 0   ! elev
      offset(2) = np ! temp
      offset(3) = np*(1+nvrt) ! salt
      offset(4) = np*(1+2*nvrt) ! uu2
      offset(5) = np*(1+3*nvrt) ! vv2
      offset(6) = np*(1+4*nvrt) ! ww2

!     set dim_fields
      dim_fields(1)=np
      dim_fields(2)=np*nvrt
      dim_fields(3)=np*nvrt
      dim_fields(4)=np*nvrt
      dim_fields(5)=np*nvrt
      dim_fields(6)=np*nvrt
 

!     Start to read files
      open(11,file='z.bin',form="unformatted")
      open(12,file='t.bin',form="unformatted")
      open(13,file='s.bin',form="unformatted")
      open(14,file='u.bin',form="unformatted")
      open(15,file='v.bin',form="unformatted")
      open(16,file='w.bin',form="unformatted")


!     do nt=1,ntimes
!        read(11) elev(:,nt)
!        read(12) t(:,:,nt)
!        read(13) s(:,:,nt)
!        read(14) u(:,:,nt)
!        read(15) v(:,:,nt)
!        read(16) w(:,:,nt)
!        Check
!        if (nt==10) then
!           write(*,*) t(:,1,nt),s(:,10,nt)
!        end if
!     end do

!     Feed states
      do nt=1,ntimes
         read(11) temp2d(:)
         states(1:np,nt)= temp2d(:) !elev(:,nt)
         ic=np ! start from np
         !temp
         read(12) temp3d(:,:)
         !just check  order
         !if (nt.lt.4) then
         !   write(*,*) nt,1,temp3d(:,100)
         !   write(*,*) nt,2,temp3d(:,200)
         !end if
         do i=1,np
            ! Fill NaN to bottom value as model
            ikk=-1
            do k=1,nvrt
               if (abs(temp3d(k,i)).gt.9999) ikk=k
            end do
           !write(*,*) nt,'T',i,ikk,temp3d(ikk,i)
            if (ikk.eq.nvrt) then
                temp3d(1:ikk,i)=-9999.    
            else
                if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            end if
            do k=1,nvrt
               ic=ic+1
               states(ic,nt)= temp3d(k,i) !t(k,i,nt)
            end do
         end do
!        if (nt.lt.4) then
!           write(*,*) nt,1,temp3d(:,1)
!           write(*,*) nt,2,temp3d(:,2)
!        end if
         !salt
         read(13) temp3d(:,:)
         do i=1,np
            ! Fill NaN to bottom value as model
            ikk=-1
            do k=1,nvrt
               if (abs(temp3d(k,i)).gt.9999) ikk=k
            end do
           !write(*,*) nt,'S',i,ikk,temp3d(ikk,i)
            if (ikk.eq.nvrt) then
                temp3d(1:ikk,i)=-9999.    
            else
                if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            end if
           !if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            do k=1,nvrt
               ic=ic+1
               states(ic,nt)= temp3d(k,i) !s(k,i,nt)
            end do
         end do
!        if (nt.lt.4) then
!           write(*,*) nt,1,temp3d(:,1)
!           write(*,*) nt,2,temp3d(:,2)
!        end if
         !u
         read(14) temp3d(:,:)
         do i=1,np
            ! Fill NaN to bottom value as model
            ikk=-1
            do k=1,nvrt
               if (abs(temp3d(k,i)).gt.9999) ikk=k
            end do
           !write(*,*) nt,'U',i,ikk,temp3d(ikk,i)
            if (ikk.eq.nvrt) then
                temp3d(1:ikk,i)=-9999.    
            else
                if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            end if
           !if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            do k=1,nvrt
               ic=ic+1
               states(ic,nt)= temp3d(k,i) !u(k,i,nt)
            end do
         end do
         !v
         read(15) temp3d(:,:)
         do i=1,np
            ! Fill NaN to bottom value as model
            ikk=-1
            do k=1,nvrt
               if (abs(temp3d(k,i)).gt.9999) ikk=k
            end do
           !write(*,*) nt,'V',i,ikk,temp3d(ikk,i)
            if (ikk.eq.nvrt) then
                temp3d(1:ikk,i)=-9999.    
            else
                if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            end if
           !if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            do k=1,nvrt
               ic=ic+1
               states(ic,nt)= temp3d(k,i) !v(k,i,nt)
            end do
         end do
         !w
         read(16) temp3d(:,:)
         temp3d=0.001 !manualy specify due to lack of data
         do i=1,np
            ! Fill NaN to bottom value as model
            ikk=-1
            do k=1,nvrt
               if (abs(temp3d(k,i)).gt.9999) ikk=k
            end do
           !write(*,*) nt,'W',i,ikk,temp3d(ikk,i)
            if (ikk.eq.nvrt) then
                temp3d(1:ikk,i)=-9999.    
            else
                if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            end if
           !if (ikk.gt.0) temp3d(1:ikk,i)=temp3d(ikk+1,i)
            do k=1,nvrt
               ic=ic+1
               states(ic,nt)= temp3d(k,i) !w(k,i,nt)
            end do
         end do
         !write(*,*) nt,ic
         if (ic.ne.dim_state) then
            write(*,*) 'states array counter is wrong!'
            stop
         end if
      end do !nt

!     Allocate arrays for singular values and vectors
      ALLOCATE(svals(ntimes))
      ALLOCATE(svdU(dim_state, ntimes))

!     Allocate mean_state
      ALLOCATE(meanstate(dim_state))
      ALLOCATE(run_meanstate(dim_state, ntimes))

!     Compute running mean
      if (remove_mean.eq.0) then
      DO i = 1, ntimes
         meanstate = 0.0

         IF (i > hwindow .AND. i <= (ntimes - hwindow)) THEN
            Do j = i - hwindow, i + hwindow
               meanstate = meanstate + states (:, j)
            END DO

         ELSE IF (i <= hwindow) THEN
            DO j = 1, i + hwindow
               meanstate = meanstate + states (:, j)
            END DO
            DO j = ntimes - hwindow +i, ntimes
               meanstate = meanstate + states (:, j)
            END DO

         ELSE IF (i > ntimes - hwindow) THEN
            DO j = i - hwindow, ntimes
               meanstate = meanstate + states (:, j)
            END DO
            DO j = 1, i + hwindow - ntimes
               meanstate = meanstate + states (:, j)
            END DO
          END IF

          meanstate = meanstate / (2 * hwindow + 1)
          run_meanstate (:, i) = meanstate

      END DO

!     get residual
          states = states - run_meanstate
      end if


!     Call PDAF eofcov
      CALL PDAF_eofcovar(dim_state, ntimes, nfields, dim_fields, offset, &
            remove_mean, do_mv, states, stddev, svals, svdU, meanstate, 1, status)

! *** Determine rank to write ***

      getlimit: DO i = 1, ntimes
      IF (svals(i) >= limit) THEN
          rank = i
       ELSE
          EXIT getlimit
      END IF
      END DO getlimit
      IF (rank < ntimes) THEN
         WRITE (*,'(1x,a,i6,a,es10.2)') &
          'Use maximum of ', rank, ' eigenvectors due to eigenvalue-limit of ',limit
      END IF
      IF (rank == ntimes) THEN
         rank = ntimes - 1
         WRITE (*,'(5x,a,i4)') '++ reset rank to ',rank
      END IF

      WRITE (*,'(5x,a)') 'Singular values: '
      DO i = 1, rank
         WRITE (*, '(10x, i4, es12.3)') i, svals(i)
      END DO


      WRITE (*,'(/1x,a)') '------- Write decomposed covariance matrix -------------'


      open(21,file='meanstate.bin',form="unformatted",status="replace")
      open(22,file='eofs.bin',form="unformatted",status="replace")
      open(23,file='svd.bin',form="unformatted",status="replace")
      open(24,file='eofs.dat',status="replace")

!     write(21) run_meanstate(:,ioutrank) ! here only write last rank, try meanstate
      write(21) run_meanstate(:,1) ! here only 1st
!     write(*,*) 'ck',run_meanstate(np+1:np+nvrt,1)
!     write(*,*) 'ck',run_meanstate(np+nvrt+1:np+2*nvrt,1)
!     write(21) meanstate
!     write(22) svdU(:,1:ioutrank) !change format to reduce memory usage  in init_ens_pdaf
      do i=1,ioutrank
         write(22) svdU(:,i)
      end do
      do i=1,ioutrank
         write(*,*) i,svals(i)
      end do

      write(23) ioutrank
      write(23) svals(1:ioutrank)
      write(24,*) int(ioutrank)


! ********************
! *** Finishing up ***
! ********************

      DEALLOCATE(states, meanstate)
      DEALLOCATE(svals, svdU)

      WRITE (*,'(/1x,a/)') '------- END -------------'


      end
