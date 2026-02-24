!$Id: init_n_domains_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
!
! !INTERFACE:
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to set the number of local analysis 
! domains for the PE-local domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use schism_glbl, only: npa,nvrt,ntracers,xnd,ynd,xlon,ylat,znl,ics,pi,kbp,idry
  use mod_assimilation, only: mdl_coords_p,idx_domain_p,domain_limit_depth,Vlocal_opt,vidx_domain_p,local_range
  use mod_parallel_pdaf, only: mype_model
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains

  !Local var
  integer :: ic,i,j,k,ic2
! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_n_domains)
! Called by: PDAF_lestkf_update  (as U_init_n_domains)
! Called by: PDAF_letkf_update   (as U_init_n_domains)
! Called by: PDAF_lnetf_update   (as U_init_n_domains)
!EOP


! ************************************
! *** Initialize number of domains ***
! ************************************
 
  ! Template reminder - delete when implementing functionality
! WRITE (*,*) 'TEMPLATE init_n_domains_pdaf.F90: Set number of local analysis domains here!'

  !n_domains_p = npa ! here we choose npa
  !n_domains_p = npa*(1+nvrt*ntracers+3*nvrt) ! here we choose all for 2D+1D local

  !if (Vlocal_opt.eq.0) then !2D (horizontal only)
     n_domains_p = npa ! here we choose npa
     !Move below to init_pdaf since these are rarely change in the deeper area
     !
     !Derive vidx_domain_p
     !vidx_domain_p=1 !Reset every DA cycle
     !do i=1,npa
     !   do k=nvrt,kbp(i),-1 !Reverse search vert column from surface, since obs are located at surface mostly
     !      if (znl(k,i).gt.(0.d0-local_range(3))) vidx_domain_p(i)=k
     !   end do
     !end do
     !Search with obs_p !Not implement yet 

  if (Vlocal_opt.eq.1) then !2D +1D
  !else !2D+1D (expensive)
    idx_domain_p=0 !reset every DA cycle
    !Set mdl_coords_p for 2D+1D local
    ic=0
    ic2=0
    do i=1,npa !elev
       ic=ic+1
       mdl_coords_p(3,ic)=0.d0
       if (idry(i).eq.0) then 
           idx_domain_p(1,ic)=1  !on/off
           ic2=ic2+1
           idx_domain_p(2,ic2)=ic !index
       end if
    end do

    !Only update znl
    do j=1,ntracers+3 ! + 3 = U,V,W
      do i=1,npa
        do k=1,nvrt
           ic=ic+1
           if (k.le.(kbp(i)-1)) then
              mdl_coords_p(3,ic)=-99999.d0 ! set as large value to avoid localization
              !if ((mype_model.eq.0).and.(i.eq.10).and.(j.eq.1)) write(*,'(a,i,f12.2)') 'In init_n_domains_pdaf, mdl_coords_p: ',k,mdl_coords_p(3,ic)
           else
              mdl_coords_p(3,ic)=znl(k,i)
              !if ((mype_model.eq.0).and.(i.eq.10).and.(j.eq.1)) write(*,'(a,i,f12.2)') 'In init_n_domains_pdaf, mdl_coords_p: ',k,mdl_coords_p(3,ic)
              if ((idry(i).eq.0).and.(mdl_coords_p(3,ic).gt.(0.d0-domain_limit_depth))) then
                      idx_domain_p(1,ic)=1 !Skip bottom & domain_limit
                 ic2=ic2+1
                 idx_domain_p(2,ic2)=ic !index
              end if
           end if
        end do !k
      end do !i 
    end do !j

     n_domains_p=sum(idx_domain_p(1,:)) !overwrite with wet points only
    !if (mype_model.eq.0) write(*,'(a,i)') 'In init_n_domains_pdaf,n_domains_p=',n_domains_p

  end if
END SUBROUTINE init_n_domains_pdaf
