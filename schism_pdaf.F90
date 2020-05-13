! This code is a main driver for a coupled SCHISM and PDAF
! program for running multiple schism components concurrently in flexible mode
! (i.e. multiple tasks can share same PET list). The coupling model interface is
! schism_cmi_esmf.F90 (and interfaces are defined in schism_bmi.F90)
!
! @copyright (C) 2018, 2019, 2020 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister
! @author Carsten Lemmen <carsten.lemmen@hzg.de>
! @author Y Joseph Zhang <yjzhang@vims.edu>
!
! @license under the Apache License, Version 2.0 (the "License");
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

  implicit none

  include 'mpif.h'

  !> @todo use this routine from schism_esmf_util, delete local one
  interface
    function clockCreateFrmParam(filename, rc)
      use esmf
        character(len=ESMF_MAXSTR), intent(in) :: filename
        integer(ESMF_KIND_I4), intent(out)     :: rc
        type(ESMF_Clock)                       :: clockCreateFrmParam
    end function clockCreateFrmParam
  end interface

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

  integer(ESMF_KIND_I4)       :: concurrentCount, ncohort
  integer(ESMF_KIND_I4)       :: sequenceIndex

  integer(ESMF_KIND_I4)         :: alarmCount=0, ringingAlarmCount=0
  type(ESMF_Alarm), allocatable :: alarmList(:), ringingAlarmList(:)
  logical                       :: hasAlarmRung = .false.
  character(len=ESMF_MAXSTR)    :: alarmName

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
  elseif (petCount<concurrentCount) then
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
  !> concurrent runs.  As an example, we assume
  !> petCount=48, schismCount=14, concurrentCount=5
  !> We integer divide petCount by concurrentCount
  !> and obtain the petLists {0..8} {9..17} {18..26}
  !> {27..35} {36..44}, i.e. 9 cores for each task.
  !> NOTE that # of cores must be
  !> equal as we do not want to repartition the grid etc.

  !> We distribute the schismCount instances on the concurrentCount
  !> to obtain the number of cohorts that run sequentially,
  !> and obtain ncohort=3 in our example
  ncohort = ceiling(real(schismCount)/concurrentCount)
  petCountLocal = petCount/concurrentCount

  !> Thus, we obtain the sets of instances that run concurrently: {1..5} {6..10}
  !> {11..14}

  !write(0,*) 'schism count, cohort count, maxpercohort ', schismCount, concurrentCount, ncohort
  do i = 1, schismCount

    ! Determine the sequence  and concurrent index of each
    ! instance
    sequenceIndex = mod(i-1, concurrentCount) !local index in a cohort (0- based)

    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! Instances within the same cohort use different PET
    do j=1, petCountLocal
      petList(j)=sequenceIndex*petCountLocal+ j-1 !points to global PET #
    end do !j

    !write(0,*) 'Instance ', i, ' uses ', petCountLocal, ' PET in sequence  ', sequenceIndex, 'list is ',  petlist

    write(message2, '(A,I3.3)') 'schism_', i
    !write(0, *) trim(message2), 'list=', petList, 'petCount=', petCount, petCountLocal

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

    !Put input dir name into attribute to pass onto init P1
    call ESMF_AttributeSet(schism_components(i), name='input_directory', &
      value=trim(message2), rc=localrc)
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

  ! Initialize phase 0 and set attribute for calling init
  init0_loop: do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=0, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !Put sequence_index into attribute to pass onto init P1
    call ESMF_AttributeSet(schism_components(i), 'cohort_index', &
      (i-1)/concurrentCount)
  enddo init0_loop

  !Get info on simulation period
  filename = './global.nml'
  clock = clockCreateFrmParam(filename, localrc)

  ! Read some information on observation data availability and create
  ! a list of alarms for the times that new data is available and should be
  ! assimilated into the model.
  ! filename = './alarmlist.csv'
  ! call alarmListCreateFrmFile(filename, clock)
  ! this function internally calls ESMF_AlarmCreate(clock,ringTime, &
  !  ringInterval, name, rc=localrc)
  ! for each alarm time present in the file
  if (allocated(alarmList)) then
    alarmCount=ubound(alarmList,1)
  endif

  ! Init phase 1: assuming all ensemble members share same parameters and i.c.
  ! and use same # of PETs. The PETs seem to 'block' when doing a task until
  ! it's done so the latter tasks that use same PETs will wait
  init1_loop: do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=1, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo init1_loop

  reconcile_loop: do i = 1, schismCount

    call ESMF_StateReconcile( importStateList(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateReconcile(exportStateList(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo reconcile_loop

! Init PDAF env
  call init_parallel_pdaf(0,1,schismCount,petCountLocal,concurrentCount)
  call init_pdaf(schismCount,j)
  if(j/=0) then
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    do i = 1, schismCount
      call ESMF_GridCompRun(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo

    call MPI_barrier(MPI_COMM_WORLD,ii)
    if(ii/=MPI_SUCCESS) call MPI_abort(MPI_COMM_WORLD,0,j)

    ! Advance coupling clock to next timestep, any alarms that are going off
    ! during this timestep will be set to the ringing state
    call ESMF_ClockAdvance(clock, ringingAlarmCount=ringingAlarmCount, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    if (ringingAlarmCount > 0) then

      if (allocated(ringingAlarmList)) deallocate(ringingAlarmList)
      allocate(ringingAlarmList(ringingAlarmCount))

      !call ESMF_ClockPrint(clock, options="currTime string", message, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      write(message,'(A,I4)') 'Number of ringing alarms = ', ringingAlarmCount
      print *, trim(message)

      call ESMF_ClockGetAlarmList(clock, ESMF_ALARMLIST_RINGING, &
        alarmList=ringingAlarmList, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      do i = 1, ringingAlarmCount

        call ESMF_AlarmGet(ringingAlarmList(i), name=alarmName, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

        print *, trim(alarmName), " is ringing!"
        hasAlarmRung=.true.

        call ESMF_AlarmRingerOff(ringingAlarmList(i), rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      enddo ! i = 1, ringinAlarmCount
    endif
    if (allocated(ringingAlarmList)) deallocate(ringingAlarmList)

    if (hasAlarmRung) then

        ! Do something with PDAF, be careful that each instances' clock
        ! have already advanced
        ! @todo this is where we have to consider how we want to implement this,
        ! pdaf usually does this in the timestep routine of the instance, but this
        ! is more complex with sequential runs (which we have here).

        hasAlarmRung = .false.
    endif

  end do !do while

  ! @todo also destroy the deep alarm objects
  if (allocated(alarmList)) deallocate(alarmList)

!  esmf_loop5: do j=1,ncohort
!    do ii = 1,concurrentCount
    do i = 1,schismCount
!      i=(j-1)*ncohort+ii !component #
!      if(i>schismCount) exit esmf_loop5

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
!  end do esmf_loop5 !j

  deallocate(exportStateList)
  deallocate( importStateList)
  deallocate(schism_components)
!  deallocate(configList)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize(rc=localrc)

end program main

function clockCreateFrmParam(filename, rc) result(clock)

  use esmf
  implicit none

  character(len=ESMF_MAXSTR), intent(in) :: filename
  integer(ESMF_KIND_I4), intent(out)     :: rc
  type(ESMF_Clock)                       :: clock

  logical               :: isPresent
  integer(ESMF_KIND_I4) :: unit, localrc
  type(ESMF_Time)       :: stopTime, startTime
  type(ESMF_TimeInterval) :: timeStep

  integer(ESMF_KIND_I4) :: start_year=2000, start_month=1, start_day=1
  integer(ESMF_KIND_I4) :: start_hour=0, runhours=2, rnday=2 ! rnday not used
  namelist /global/ start_year, start_month, start_day, start_hour, runhours, rnday

  inquire(file=filename, exist=isPresent)
  if (isPresent) then

    call ESMF_UtilIOUnitGet(unit, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    open(unit, file=filename, iostat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    read(unit, nml=global, iostat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    close(unit)
  endif

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
  timeStep = timeStep / 24

  clock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, &
    name='main clock', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end function clockCreateFrmParam
