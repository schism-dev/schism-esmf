!$Id: init_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf(schismCount,ierr)

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the online mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_model, &        ! Model variables
!        ONLY: nx, ny
  use schism_glbl, only: errmsg,nea,nsa,npa,nvrt,ntracers
  use schism_msgp, only: parallel_abort
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe !, abort_parallel
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       rms_obs, incremental, covartype, type_forget, forget, &
       rank_analysis_enkf, locweight, local_range, srange, &
       filename, type_trans, type_sqrt, delt_obs,offset_field_p,varscale, &
       ihfskip_PDAF,nspool_PDAF,outf, nhot_PDAF, nhot_write_PDAF, &
       rms_type,rms_obs2,ens_init
! use PDAF_mod_filter, only: dim_p,state !just for check

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
! Calls: PDAF_get_state
!EOP

  integer, intent(in) :: schismCount
  integer, intent(out) :: ierr !return error code

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps,i,lfdb     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation
  character(len=72) :: indir
  integer :: j

! temp-add
! real :: state_p(dim_p)

  ! External subroutines
  EXTERNAL :: init_ens_pdaf         ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens            ! User supplied pre/poststep routine
  
  NAMELIST /pdaf_nml/ screen, filtertype, subtype, &
           delt_obs, rms_obs, &
           type_forget, forget, type_trans, type_sqrt, &
           locweight, local_range, srange,varscale, &
           ihfskip_PDAF,nspool_PDAF,outf,dim_ens, &
           nhot_PDAF, nhot_write_PDAF, &
           rms_type,rms_obs2,ens_init


! ***************************
! ***   Initialize PDAF   ***
! ***************************

  ierr=0

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

! WRITE (*,*) 'TEMPLATE init_pdaf.F90: Initialize state dimension here!'

  ! *** Define state dimension (state var is a long 1D array)
!  dim_state = ?

  !Order of arrays:
  !old setting
  !idry_e,we,tr_el, idry_s,su2,sv2, idry,eta2,tr_nd,tr_nd0,q2,xl,dfv,dfh,dfq1,dfq2
  !dim_state_p=nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+6*nvrt)

  !New setting
  !eta2,tr_nd,uu2,vv2,ww2 --> choose npa vars, and convert them back in distribute routine
  dim_state_p=npa*(1+nvrt*ntracers+3*nvrt)

! Define process-local offsets in local state vector
  allocate(offset_field_p(5)) ! Set 5
  offset_field_p(1) = 0   ! eta2
  offset_field_p(2) = npa ! tr_nd
  offset_field_p(3) = npa*(1+nvrt*ntracers) ! uu2
  offset_field_p(4) = npa*(1+nvrt*ntracers+nvrt) ! vv2
  offset_field_p(5) = npa*(1+nvrt*ntracers+2*nvrt) ! ww2

! The followings are old setting for all hotstart vars
! offset_field_p(1) = 0                                ! idry_e
! offset_field_p(2) = nea                              ! we
! offset_field_p(3) = nea*(1+nvrt)                     ! tr_el
! offset_field_p(4) = nea*(1+nvrt+nvrt*ntracers)              ! idry_s
! offset_field_p(5) = nea*(1+nvrt+nvrt*ntracers)+nsa          ! su2
! offset_field_p(6) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+nvrt) ! sv2
! offset_field_p(7) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)         ! idry
! offset_field_p(8) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa     ! eta2
! offset_field_p(9) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*2   ! tr_nd
! offset_field_p(10) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+nvrt*ntracers)   ! tr_nd0
! offset_field_p(11) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers)        ! q2
! offset_field_p(12) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+nvrt)   ! xl
! offset_field_p(13) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+2*nvrt) ! dfv
! offset_field_p(14) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+3*nvrt) ! dfh
! offset_field_p(15) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+4*nvrt) ! dfq1
! offset_field_p(16) = nea*(1+nvrt+nvrt*ntracers)+nsa*(1+2*nvrt)+npa*(2+2*nvrt*ntracers+5*nvrt) ! dfq2


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 3  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 6    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (8) localized EnKF
                    !   (9) NETF
                    !  (10) LNETF
  !dim_ens = 8       ! Size of ensemble for all ensemble filters
  dim_ens =schismCount ! Size of ensemble for all ensemble filters (all tasks)
                    ! Number of EOFs to be used for SEEK
  subtype = 0       ! subtype of filter: 
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    !   (0) for dim_ens^-1 (old SEIK)
                    !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
                    !   This parameter has also to be set internally in PDAF_init.
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 36     ! Number of time steps between analysis/assimilation steps

! *** specifications for observations ***
  rms_obs = 0.5    ! Observation error standard deviation
                   ! for the Gaussian distribution 
! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  local_range = 1  ! Range in grid points for observation domain in local filters
  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting
  varscale = 1.0    ! Init ensemble variance
  outf = 1          ! output file switch
  nhot_PDAF = 0     ! hotstart rank output switch

! *** File names
  filename = 'output.dat'

! Read pdaf.nml to control filtertype & other parameters
  OPEN (500,file='pdaf.nml')
  READ (500,NML=pdaf_nml)
  CLOSE (500)
  srange = local_range  ! Support range for 5th-order polynomial

! Check pdaf.nml, stop if parameters are wrongly set. 
  if (mod(ihfskip_PDAF,delt_obs).ne.0) then
     write(errmsg,*) 'ihfskip_PDAF has to be multiple of delt_obs, please change accordingly!'
     CALL parallel_abort(errmsg)
  end if
  if (delt_obs.ge.nspool_PDAF) then
     if (mod(delt_obs,nspool_PDAF).ne.0) then
        write(errmsg,*) 'nspool_PDAF and delt_obs must be multiple each other, please change accordingly!'
        CALL parallel_abort(errmsg)
     end if
  else
     if (mod(nspool_PDAF,delt_obs).ne.0) then
        write(errmsg,*) 'nspool_PDAF and delt_obs must be multiple each other, please change accordingly!'
        CALL parallel_abort(errmsg)
     end if
  end if
  if(nhot_PDAF/=0.and.nhot_PDAF/=1.or.nhot_PDAF*mod(nhot_write_PDAF,ihfskip_PDAF)/=0) then
     write(errmsg,*)'Unknown hotout or hotout_write is not multiple of ihfskip in PDAF',nhot_PDAF,ihfskip_PDAF
     call parallel_abort(errmsg)
  endif
  if ((rms_type==2).and.(sum(rms_obs2(:)).eq.0.d0)) then
     write(errmsg,*) 'Please specify rms_obs2, ', (rms_obs2(j),j=1,5)
     call parallel_abort(errmsg)
  end if
  if ((rms_type<1).or.(rms_type>3)) then
     write(errmsg,*) 'Please specify right rms_type (1~3), now is ', rms_type
     call parallel_abort(errmsg)
  end if
  if ((ens_init<1).or.(ens_init>3)) then
     write(errmsg,*) 'Please specify right ens_init (1~3), now is ', ens_init
     call parallel_abort(errmsg)
  end if


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

!  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

! init_ens_pdaf() inside will init state_p(:) (different filters have different
! init routines)
 whichinit: IF (filtertype == 2) THEN
    ! *** EnKF with Monte Carlo init ***
    filter_param_i(1) = dim_state_p ! State dimension
    filter_param_i(2) = dim_ens     ! Size of ensemble
    filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
    filter_param_i(4) = incremental ! Whether to perform incremental analysis
    filter_param_i(5) = 0           ! Smoother lag (not implemented here)
    filter_param_r(1) = forget      ! Forgetting factor
    
    CALL PDAF_init(filtertype, subtype, 0, &
         filter_param_i, 6,&
         filter_param_r, 2, &
         COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, filterpe, init_ens_pdaf, &
         screen, status_pdaf)
 ELSE
    ! *** All other filters                       ***
    ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
    filter_param_i(1) = dim_state_p ! State dimension
    filter_param_i(2) = dim_ens     ! Size of ensemble
    filter_param_i(3) = 0           ! Smoother lag (not implemented here)
    filter_param_i(4) = incremental ! Whether to perform incremental analysis
    filter_param_i(5) = type_forget ! Type of forgetting factor
    filter_param_i(6) = type_trans  ! Type of ensemble transformation
    filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
    filter_param_r(1) = forget      ! Forgetting factor
    
    CALL PDAF_init(filtertype, subtype, 0, &
         filter_param_i, 7,&
         filter_param_r, 2, &
         COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, filterpe, init_ens_pdaf, &
         screen, status_pdaf)
 END IF whichinit

! Call collect_state_pdaf to initialize state_p for each member
! call collect_state_pdaf(dim_p,state)
! state_p=state
! write(*,*) 'in init_pdaf, state_p:',maxval(state_p),mype_world,task_id,filterpe

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     ierr=1
     WRITE (errmsg,*)'init_pdaf error ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL parallel_abort(errmsg)
  END IF


! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***
!new28: This is mainly to init obs
! This must to open for flexible mode initialization

!  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
!       distribute_state_pdaf, prepoststep_ens, status_pdaf)
 
!  IF (status_pdaf /= 0) THEN
!     ierr=1
!     WRITE (errmsg,*)'init_pdaf error ', status_pdaf, &
!          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
!     CALL parallel_abort(errmsg)
!  END IF


END SUBROUTINE init_pdaf
