! This code is a main driver for a coupled SCHISM and PDAF
! program for running multiple schism components concurrently in flexible mode
! (i.e. multiple tasks can share same PET list). The coupling model interface is
! schism_cmi_esmf.F90 (and interfaces are defined in schism_bmi.F90)
!
! @copyright (C) 2018, 2019, 2020-2021 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
! @author Y Joseph Zhang <yjzhang@vims.edu>
!
! @license Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! 		http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! This version accounts for halo(ghost) zone, because ESMF by default
! partitions among nodes instead of elements

!TODO: search for new28

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_pdaf.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
! use schism_msgp, only: parallel_abort,myrank
! use schism_glbl, only: errmsg,tr_el
! USE mod_assimilation, &      ! Variables for assimilation
!      ONLY: filtertype
! use PDAF_mod_filter, only: dim_p,state

  implicit none

  include 'mpif.h'

  type(ESMF_GridComp), allocatable :: schism_components(:)
  type(ESMF_State), allocatable    ::  importStateList(:), exportStateList(:)

  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time
  type(ESMF_Clock)        :: clock, childClock

  type(ESMF_Vm)           :: vm
  type(ESMF_Log)          :: log

  integer(ESMF_KIND_I4)       :: petCountLocal, schismCount=8
  integer(ESMF_KIND_I4)       :: rc, petCount, i, j, inum,localrc, ii
  integer(ESMF_KIND_I4), allocatable    :: petlist(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename='schism_pdaf.cfg', message, message2
  type(ESMF_Config)           :: config
  type(ESMF_Config), allocatable :: configList(:)

  integer(ESMF_KIND_I4)       :: concurrentCount, ncohort,cohortIndex
  integer(ESMF_KIND_I4)       :: sequenceIndex

  integer(ESMF_KIND_I4)         :: alarmCount=0, ringingAlarmCount=0
  type(ESMF_Alarm), allocatable :: alarmList(:), ringingAlarmList(:)
  logical                       :: hasAlarmRung = .false.
  character(len=ESMF_MAXSTR)    :: alarmName

  integer(ESMF_KIND_I4) :: schism_dt,num_schism_dt_in_couple,runhours,num_obs_steps,it,unit,next_obs_step
  integer(ESMF_KIND_I4), allocatable :: list_obs_steps(:)
  real(ESMF_KIND_R8), allocatable :: list_obs_times(:)

! PDAF vars & external routines
! integer doexit,status_pdaf ! PDAF vars
! real time_PDAF
! EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time,
!                                      ! and dimension of next observation
!      distribute_state_pdaf, &        ! Routine to distribute a state vector to  model fields
!      prepoststep_ens                 ! User supplied pre/poststep routine

  namelist /obs_info/list_obs_times

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo for performance reasons, this should be disabled,
  !> preferentially with an #ifdef DEBUG or cmake variable
  call ESMF_LogSet(flush=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! PETs for SCHISM

  call ESMF_VMGetGlobal(vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! PET is persistent globally
  call ESMF_VMGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Get config for number of schism instances and size of cohorts
  ! A cohort is defined as the set of instances that run
  ! concurrently, i.e. on different petLists.
  inquire(file=filename, exist=isPresent)
  if (isPresent) then
    config = ESMF_ConfigCreate(rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ConfigLoadFile(config, filename=filename, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ConfigGetAttribute(config, value=schismCount, label='schism_count:', &
      default=min(petCount,8), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !Components in a cohort are run concurrently on differen PETs. Some components in different
    !cohorts share same PETs
    call ESMF_ConfigGetAttribute(config, value=concurrentCount, label='concurrent_count:', &
      default=max(petCount/schismCount,1), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
    write(message,*)'Please supply .cfg'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if (schismCount>999 .or. schismCount<1) then
    write(message, '(A,I3,A)') 'Number of instances ', &
      schismCount, 'must be in the range [1,999]'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  elseif (concurrentCount > schismCount .or. concurrentCount < 1) then
    write(message, '(A,I3,A,I3,A)') 'Number of cohorts ', &
      concurrentCount, 'must be in the range [1,', schismCount, ']'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !elseif (petCount<concurrentCount) then
  elseif (mod(petCount,concurrentCount)/=0) then !PDAF requires this
    write(message, '(A,I3,A,I3,A,I6,A)') 'Cannot run ', &
      concurrentCount, ' instances on ', petCount, ' PET'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
    write(message, '(A,I3,A,I3,A)') 'Running ', schismCount, &
    ' schism instances in ', concurrentCount, ' concurrent tasks '
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Create all components on their respective parallel
  ! environment provided by each petList
  allocate(schism_components(schismCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(configList(schismCount), stat=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Here, we partition the schismCount on the concurrentCount
  !> concurrent runs. As an example, we assume
  !> petCount=48, schismCount=14, concurrentCount=5
  !> We integer divide petCount by concurrentCount
  !> and obtain the petLists {0..8} {9..17} {18..26}
  !> {27..35} {36..44}, i.e. 9 cores for concurrnt tasks (1,4,7,10,13 etc).
  !> NOTE that # of cores must be
  !> equal as we do not want to repartition the grid etc.

  !> We distribute the schismCount instances on the concurrentCount
  !> to obtain the number of cohorts that run sequentially,
  !> and obtain ncohort=3 in our example
  ncohort = ceiling(real(schismCount)/concurrentCount)
  petCountLocal = petCount/concurrentCount

  if(ncohort<1) then
    write(message,*) 'ncohort<1:',ncohort
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  !> Thus, we obtain the sets of instances that run concurrently: {1 4 7 10 13}
  !etc folowing column major

  !write(0,*) 'schism count, cohort count, maxpercohort ', schismCount, concurrentCount, ncohort
  do i = 1, schismCount

    ! Determine the sequence  and concurrent index of each
    ! instance
!    sequenceIndex = mod(i-1, concurrentCount) !local index in a cohort (0- based); task ID-1 in PDAF
    sequenceIndex=(i-1)/ncohort

    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! Instances within the same cohort use different PET
    do j=1, petCountLocal
      petList(j)=sequenceIndex*petCountLocal+ j-1 !points to global PET #
    end do !j

!   write(0,*) 'Instance ', i, ' uses ', petCountLocal, ' PET in sequence  ', sequenceIndex, 'list is ',  petlist

    write(message2, '(A,I3.3)') 'schism_', i
!   write(0, *) trim(message2), 'list=', petList, 'petCount=', petCount, petCountLocal

    schism_components(i) = ESMF_GridCompCreate(name=trim(adjustl(message2)), &
      petList=petlist, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A,I3.3,A,I3.3,A,I6.6,A,I6.6,A,I6.6)') 'Instance schism_', i, &
      ' with sequence number ', sequenceIndex, ' created on ', petCountLocal, &
      ' PET  ', petList(1), '..', petList(petCountLocal) !,message2
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

    !configList(i) = ESMF_ConfigCreate(rc=localrc)

    !> @todo  the ESMF_ConfigSetAttribute implementation is buggy and thus not enabled, we
    !> use for now the attribute of the component.
    !>
    !call ESMF_ConfigSetAttribute(configList(i), value=i, &
    !  label='schismInstance:', rc=localrc)

    !Put input dir name into attribute to pass onto init P1 etc
    call ESMF_AttributeSet(schism_components(i), name='input_directory', value=trim(message2), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    deallocate(petList)

    !call ESMF_GridCompSet(schism_components(i), config=configList(i), rc=localrc)
    !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do ! loop over schismCount

  write(message, '(A,I3,A)') 'Created ',schismCount,' instances of SCHISM'

  allocate(exportStateList(schismCount))
  allocate(importStateList(schismCount))

  ! Create export states for all instances and register them
  register_loop:  do i = 1, schismCount

    call ESMF_GridCompSetServices(schism_components(i), &
      schismSetServices, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A,I3.3)') 'schism_', i

    exportStateList(i) = ESMF_StateCreate( &
      name='schismExport_'//trim(adjustl(message)), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    importStateList(i) = ESMF_StateCreate( &
      name='schismImport_'//trim(adjustl(message)), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo register_loop

  !Get info on simulation period
!  filename = './global.nml'
!  clock = clockCreateFrmParam(filename, localrc)
  call clockCreateFrmParam(clock,schism_dt,num_schism_dt_in_couple,runhours,num_obs_steps)

  !Read in obs times
  allocate(list_obs_steps(num_obs_steps),list_obs_times(num_obs_steps))
  call ESMF_UtilIOUnitGet(unit, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  open(unit, file='global.nml', iostat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  read(unit, nml=obs_info, iostat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  close(unit)

  !Make sure the list is ascending
  do i=1,num_obs_steps-1
    if(list_obs_times(i+1)<=list_obs_times(i).or.list_obs_times(i)<=0) then
      write(message,*) 'Check obs times:',i,list_obs_times
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
      localrc = ESMF_RC_VAL_OUTOFRANGE
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    endif
  enddo !i
  !Calculate SCHISM time steps for DA
  list_obs_steps=nint(list_obs_times/schism_dt)
  where(list_obs_steps<1) list_obs_steps=1

  !Make sure obs steps are multiples of num_schism_dt_in_couple
  do i=1,num_obs_steps
    if(mod(list_obs_steps(i),num_schism_dt_in_couple)/=0) then
      write(message,*) 'Obs steps must be divisible by num_schism_dt_in_couple:',i,list_obs_steps(i)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
      localrc = ESMF_RC_VAL_OUTOFRANGE
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    endif
  enddo !i

  write(message,*)'List of obs steps:',list_obs_steps
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! Initialize phase 0 and set attribute for calling init
! write(0, *) 'Before init0_loop, ESMF_GridCompInitialize'
  init0_loop: do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=0, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !Put cohort index (0-based) into attribute to pass onto init P1 etc
!   cohortIndex=(i-1)/concurrentCount
    cohortIndex=mod(i-1,ncohort) ! correct this to match PDAF task-ID
    call ESMF_AttributeSet(schism_components(i), 'cohort_index', cohortIndex)

    if(cohortIndex>ncohort-1) then
      write(message,*) 'cohortIndex>ncohort-1:',cohortIndex,ncohort
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
      localrc = ESMF_RC_VAL_OUTOFRANGE
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    endif

    call ESMF_AttributeSet(schism_components(i), name='ncohort', value=ncohort, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_AttributeSet(schism_components(i), name='runhours',value=runhours, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_AttributeSet(schism_components(i), name='schism_dt2',value=schism_dt, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo init0_loop

  ! Read some information on observation data availability and create
  ! a list of alarms for the times that new data is available and should be
  ! assimilated into the model.
  ! filename = './alarmlist.csv'
  ! call alarmListCreateFrmFile(filename, clock)
  ! this function internally calls ESMF_AlarmCreate(clock,ringTime, &
  !  ringInterval, name, rc=localrc)
  ! for each alarm time present in the file
!  if (allocated(alarmList)) then
!    alarmCount=ubound(alarmList,1)
!  endif

  ! Init phase 1: assuming all ensemble members share same parameters and i.c.
  ! and use same # of PETs. The PETs seem to 'block' when doing a task until
  ! it's done so the latter tasks that use same PETs will wait
! write(0, *) 'Before schism_init, ESMF_GridCompInitialize'
  init1_loop: do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=1, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo init1_loop

! write(0, *) 'Before ESMF_StateReconcile'
  reconcile_loop: do i = 1, schismCount

    call ESMF_StateReconcile( importStateList(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateReconcile(exportStateList(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo reconcile_loop

! Init PDAF env
#ifdef USE_PDAF
! write(0, *) 'Before init_parelle_pdaf'
  call init_parallel_pdaf(0,1,schismCount,petCountLocal,concurrentCount)
! write(0, *) 'Before init_pdaf'
  call init_pdaf(schismCount,j)
#endif
  if(j/=0) then
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! Loop over coupling timesteps until stopTime
  it=0
  next_obs_step=1 !init next obs step
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !new28: this is mostly to init analysis routines
    !call PDAF_get_state
    !call PDAF_get_state(it,time_PDAF, doexit, next_observation_pdaf, distribute_state_pdaf, prepoststep_ens, status_pdaf)

    it=it+num_schism_dt_in_couple !SCHISM step #

    do i = 1, schismCount

      !Check if it's analysis step in PDAF
      if (it==list_obs_steps(next_obs_step)) then !DA step
        !Let Run know it's analysis step
        call ESMF_AttributeSet(schism_components(i), name='analysis_step', &
        value=1, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      else
        call ESMF_AttributeSet(schism_components(i), name='analysis_step', &
        value=0, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      endif !DA step

      call ESMF_GridCompRun(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    enddo !i

    if (it==list_obs_steps(next_obs_step)) next_obs_step=min(num_obs_steps,next_obs_step+1)

    call MPI_barrier(MPI_COMM_WORLD,ii)
    if(ii/=MPI_SUCCESS) call MPI_abort(MPI_COMM_WORLD,0,j)

!   Move PDAF_put_state here to match flexible mode
!   Disable local filter for dev
!   IF (filtertype == 4) THEN
!      CALL PDAF_put_state_etkf_si(status_pdaf)
!   ELSEIF (filtertype == 5) THEN
!      CALL PDAF_put_state_letkf_si(status_pdaf)
!   ELSEIF (filtertype == 6) THEN
!      CALL PDAF_put_state_estkf_si(status_pdaf)
!   ELSEIF (filtertype == 7) THEN
!      CALL PDAF_put_state_lestkf_si(status_pdaf)
!   ELSE
!      WRITE (errmsg,*) 'PDAF Filtertype only accept 4,6, please specify right one!'
!      CALL parallel_abort(errmsg)
!   END IF

!   IF (status_pdaf /= 0) THEN
!      WRITE (errmsg,*) &
!           'ERROR ', status_pdaf, &
!           ' in PDAF_put_state(assimilate_pdaf) - stopping! '
!      CALL parallel_abort(errmsg)
!   END IF
!   new28!
!   Can we add "call schism_save_state(cohortIndex)" here to update state_p?
!   write(*,*) 'in schism_pdaf:',tr_el(1,1,1),myrank


!   DA step
!   if (it==list_obs_steps(next_obs_step)) then ! DA step (update ens field)  
!       time_PDAF=list_obs_steps(next_obs_step)
!       write(*,*) 'Before PDAF_get_state in schism_pdaf!', it,next_obs_step,num_obs_steps,time_PDAF,doexit
!       call PDAF_get_state(it,time_PDAF, doexit, next_observation_pdaf, distribute_state_pdaf, prepoststep_ens, status_pdaf)
!       write(*,*) 'After PDAF_get_state in schism_pdaf!', it,next_obs_step,num_obs_steps,time_PDAF,doexit
!       next_obs_step=min(num_obs_steps,next_obs_step+1)
!   endif ! DA step


    ! Advance coupling clock to next timestep
    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!    ! Advance coupling clock to next timestep, any alarms that are going off
!    ! during this timestep will be set to the ringing state
!    call ESMF_ClockAdvance(clock, ringingAlarmCount=ringingAlarmCount, rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!    if (ringingAlarmCount > 0) then
!
!      if (allocated(ringingAlarmList)) deallocate(ringingAlarmList)
!      allocate(ringingAlarmList(ringingAlarmCount))
!
!      !call ESMF_ClockPrint(clock, options="currTime string", message, rc=localrc)
!!      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!      write(message,'(A,I4)') 'Number of ringing alarms = ', ringingAlarmCount
!!      print *, trim(message)
!      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
!
!      call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_RINGING, &
!        alarmList=ringingAlarmList, rc=localrc)
!      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!      do i = 1, ringingAlarmCount
!
!        call ESMF_AlarmGet(ringingAlarmList(i), name=alarmName, rc=localrc)
!        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!        print *, trim(alarmName), " is ringing!"
!        hasAlarmRung=.true.
!
!        call ESMF_AlarmRingerOff(ringingAlarmList(i), rc=localrc)
!        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!      enddo ! i = 1, ringinAlarmCount
!    endif !ringingAlarmCount > 0
!
!    if (allocated(ringingAlarmList)) deallocate(ringingAlarmList)
!
!    if (hasAlarmRung) then
!
!        ! Do something with PDAF, be careful that each instances' clock
!        ! have already advanced
!        ! @todo this is where we have to consider how we want to implement this,
!        ! pdaf usually does this in the timestep routine of the instance, but this
!        ! is more complex with sequential runs (which we have here).
!
!        hasAlarmRung = .false.
!    endif

  end do !do while

  ! @todo also destroy the deep alarm objects
!  if (allocated(alarmList)) deallocate(alarmList)

#ifdef USE_PDAF
! PDAF finalize
    call finalize_pdaf()
#endif

    do i = 1,schismCount
      call ESMF_GridCompFinalize(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), clock=clock, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_StateDestroy( importStateList(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_StateDestroy(exportStateList(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_GridCompDestroy(schism_components(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo !i

  deallocate(exportStateList)
  deallocate( importStateList)
  deallocate(schism_components)
!  deallocate(configList)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize(rc=localrc)

end program main

subroutine clockCreateFrmParam(clock,schism_dt,num_schism_dt_in_couple,runhours,num_obs_steps)
  use esmf
  implicit none

!  character(len=ESMF_MAXSTR), intent(in) :: filename
!  integer(ESMF_KIND_I4), intent(out)     :: rc
  type(ESMF_Clock), intent(out)          :: clock
  !SCHISM dt (sec) must be int
  integer(ESMF_KIND_I4), intent(out) :: schism_dt,num_schism_dt_in_couple,runhours,num_obs_steps

  logical               :: isPresent
  integer(ESMF_KIND_I4) :: unit, localrc,rc
  type(ESMF_Time)       :: stopTime, startTime
  type(ESMF_TimeInterval) :: timeStep

  integer(ESMF_KIND_I4) :: start_year=2000, start_month=1, start_day=1
  integer(ESMF_KIND_I4) :: start_hour=0, itmp
  namelist /sim_time/ start_year,start_month,start_day,start_hour,runhours, &
 &schism_dt,num_schism_dt_in_couple,num_obs_steps

!  inquire(file='global.nml', exist=isPresent)
!  if (isPresent) then
  call ESMF_UtilIOUnitGet(unit, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  open(unit, file='global.nml', iostat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  read(unit, nml=sim_time, iostat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  close(unit)
!  endif

  ! Set day as timestep temporarily to count later to stop time
  call ESMF_TimeSet(startTime, yy=start_year, mm=start_month, dd=start_day, &
    h=start_hour, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !'h': hour
  call ESMF_TimeIntervalSet(timeStep, h=runhours, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  stopTime = startTime + timeStep

  ! Only now define the coupling timestep as fraction of full timeStep (24
  ! coupling steps in total here)
!  timeStep = timeStep / 24

  !Set time interval in sec, used in ESMF main stepping
  itmp=schism_dt*num_schism_dt_in_couple !has to be int
  call ESMF_TimeIntervalSet(timeStep, s=itmp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  clock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, &
    name='main clock', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine clockCreateFrmParam
