! This code is a main driver for coupled SCHISM and PDAF
! program for running multiple schism components concurrently
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

  integer(ESMF_KIND_I4)       :: concurrentCount, maxLocalSchismCount
  integer(ESMF_KIND_I4)       :: sequenceIndex, concurrentIndex, maxLocalPetCount

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
  elseif (schismCount>petCount*concurrentCount) then
    write(message, '(A,I3,A,I3,A,I6,A)') 'Cannot run ', &
      schismCount, ' instances as ', concurrentCount, &
      ' on ', petCount, ' PET'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
    write(message, '(A,I3,A,I3,A)') 'Running ', schismCount, &
    ' schism instances in ', concurrentCount, ' concurrent cohorts'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Create all components on their respective parallel
  ! environment provided by each petList
  allocate(schism_components(schismCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(configList(schismCount), stat=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Here, we partition the schismCount on the concurrentCount
  !> concurrent environments.  As an example, we assume
  !> petCount=48, schismCount=14, concurrentCount=5
  !> We integer divide petCount by concurrentCount with remainder
  !> and obtain the petLists {0..9} {10..19} {20..29}
  !> {30..39} {40..47}, with respective petCountLocal [10,10,10,10,8]

  !> We distribute the schismCount instances on the concurrentCount
  !> to obtain the maximum number of cohorts that run sequentially,
  !> and obtain maxLocalSchismCount=3 in our example, and
  !> maxLocalPetCount=10
  maxLocalSchismCount = ceiling(schismCount*1.0/concurrentCount)
  maxLocalPetCount = ceiling(petCount*1.0/concurrentCount)

  !> Thus, we obtain the sets of instances that run sequentially
  !> within their set and concurrently between sets {0..2} {3..5}
  !> {6..8} {9..11} {12..13} with localSchismCounts [3,3,3,3,2]

  !write(0,*) 'schism count, cohort count, maxpercohort ', schismCount, concurrentCount, maxLocalSchismCount
  do i = 1, schismCount

    ! Determine the sequence  and concurrent index of each
    ! instance
    sequenceIndex = mod(i-1, concurrentCount)
    concurrentIndex = (i-1) / maxLocalSchismCount
    petCountLocal = min(maxlocalPetCount, petCount -concurrentIndex*maxLocalPetCount)

    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! Instances within the same cohort use the same PET
    do j=1, petCountLocal
      petList(j)=concurrentIndex*maxLocalPetCount + j-1
    end do !j

    !write(0,*) 'Instance ', i, ' uses ', petCountLocal, ' PET in sequence  ', sequenceIndex, 'list is ',  petlist

    write(message, '(A,I3.3)') 'schism_', i
    !write(0, *) trim(message), 'list=', petList, 'petCount=', petCount, petCountLocal

    schism_components(i) = ESMF_GridCompCreate(name=trim(adjustl(message)), &
      petList=petlist, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A,I3.3,A,I3.3,A,I6.6,A,I6.6,A,I6.6)') 'Instance schism_', i, &
      ' with sequence number ', sequenceIndex, ' created on ', petCountLocal, &
      ' PET  ', petList(1), '..', petList(petCountLocal)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

    !configList(i) = ESMF_ConfigCreate(rc=localrc)

    !> @todo  the ESMF_ConfigSetAttribute implementation is buggy and thus not enabled, we
    !> use for now the attribute of the component.
    !>
    !call ESMF_ConfigSetAttribute(configList(i), value=i, &
    !  label='schismInstance:', rc=localrc)

    call ESMF_AttributeSet(schism_components(i), name='input_directory', &
      value=trim(message), rc=localrc)
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

  ! Initialize phase 1 and set attribute for calling init
  init0_loop: do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=0, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_AttributeSet(schism_components(i), 'sequence_index', &
      mod(i-1, maxLocalSchismCount))
  enddo init0_loop

  !Get info on simulation period
  filename = './global.nml'
  clock = clockCreateFrmParam(filename, localrc)

  ! Init phase 1: assuming all ensemble members share same parameters and i.c.
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

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Save the time step of the current component
    call ESMF_ClockGet(clock, timeStep=timeStep, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do i = 1, schismCount

      call ESMF_GridCompRun(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    enddo

    ! Do something with PDAF, be careful that each instances' clock
    ! have already advanced

    ! Advance coupling clock to next timestep
    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do !do while

  esmf_loop5: do ii = 1,concurrentCount
    do j=1,maxLocalSchismCount
      i=(ii-1)*maxLocalSchismCount+j !component #
      if(i>schismCount) exit esmf_loop5

      call ESMF_GridCompFinalize(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), clock=clock, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_StateDestroy( importStateList(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_StateDestroy(exportStateList(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_GridCompDestroy(schism_components(i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo !j
  end do esmf_loop5 !i

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
  integer(ESMF_KIND_I4) :: start_hour=0, runhours=2, rnday=2 ! not used
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
