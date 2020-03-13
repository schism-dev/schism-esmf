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

  integer(ESMF_KIND_I4)       :: petCountLocal, schismCount
  integer(ESMF_KIND_I4)       :: rc, petCount, i, j, inum,localrc,ii
  integer(ESMF_KIND_I4)       :: ncohort, maxCountperCohort
  integer(ESMF_KIND_I4), allocatable    :: petlist(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename='schism_pdaf.cfg', message, message2
  type(ESMF_Config)           :: config
  type(ESMF_Config), allocatable :: configList(:)

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogSet(flush=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! PETs for SCHISM

  call ESMF_VMGetGlobal(vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !PET is persistent globally
  call ESMF_VMGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Get config info (# of schism instances). Always specify
  inquire(file=filename, exist=isPresent)
  if (isPresent) then
    config = ESMF_ConfigCreate(rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ConfigLoadFile(config, filename=filename, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ConfigGetAttribute(config, value=schismCount, label='schism_instances:', &
      !default: value used if label is not found
      default=min(petCount,8), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ConfigGetAttribute(config, value=ncohort, label='ncohort:', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if(schismCount>999 .or. schismCount<1 .or. schismCount>petCount*ncohort &
    .or. schismCount<ncohort) then
    write(message, '(A,I3,I6,A)') 'Number of instances, cohorts ', &
      schismCount, ncohort, &
      'must be in the range [1,999]'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  elseif (schismCount > petCount) then
!    write(message, '(A)') 'Requested number of instances exceeds available PETs'
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
!    write(message, '(A,I3,A,I3)') 'Reduced from ',schismCount, ' to ', petCount
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
!    schismCount = petCount
  else
    write(message, '(A,I3,I6)') 'Number of SCHISM instances & cohort= ',schismCount,ncohort
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Create all components on their respective parallel
  ! environment provided by each petList
  allocate(schism_components(schismCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(configList(schismCount), stat=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !max # of instances per corhort (may be smaller)
  maxCountperCohort=schismCount/ncohort
  do i = 1, schismCount

    petCountLocal = petCount*ncohort/schismCount
    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !Multiple cohort may use same PET
    do j=1, petCountLocal
      petList(j)=mod(i-1,maxCountperCohort)*petCountLocal + j - 1 ! PET #; 0 based
      if(petList(j)>petCount-1) then
        write(message2,*) 'Exceeded max PET:',petList(j),petCount-1
        call ESMF_LogWrite(trim(message2), ESMF_LOGMSG_ERROR)
        localrc = ESMF_RC_VAL_OUTOFRANGE
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      endif
    end do !j

    write(message, '(A,I3.3)') 'schism_', i
    !write(0, *) trim(message), 'list=', petList, 'petCount=', petCount, petCountLocal

    schism_components(i) = ESMF_GridCompCreate(name=trim(adjustl(message)), &
      petList=petlist, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !configList(i) = ESMF_ConfigCreate(rc=localrc)

    !> @todo  the SetAttribute implementation is buggy and thus not enabled, we
    !> use for now the attribute of the component.
    !>
    !call ESMF_ConfigSetAttribute(configList(i), value=i, &
    !  label='schismInstance:', rc=localrc)

    call ESMF_AttributeSet(schism_components(i), name='input_directory', &
      value=trim(message), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !call ESMF_GridCompSet(schism_components(i), config=configList(i), rc=localrc)
    !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    deallocate(petList)
  end do ! loop over schismCount

  write(message, '(A,I3,A)') 'Created ',schismCount,' instances of SCHISM'

  allocate(exportStateList(schismCount))
  allocate(importStateList(schismCount))

  esmf_loop1: do ii = 1,ncohort !schismCount
    do j=1,maxCountperCohort
      i=(ii-1)*maxCountperCohort+j !component #
      if(i>schismCount) exit esmf_loop1

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
    enddo !j
  end do esmf_loop1 !ii

  !Init phase 0
  esmf_loop2: do ii = 1,ncohort
    do j=1,maxCountperCohort
      i=(ii-1)*maxCountperCohort+j !component #
      if(i>schismCount) exit esmf_loop2

      call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
        exportState=exportStateList(i), phase=0, clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo !j
  end do esmf_loop2 !ii

  !Get info on simulation period
  filename = './global.nml'
  clock = clockCreateFrmParam(filename, localrc)

  !Init phase 1: assuming all ensemble members share same parameters and i.c.
  !Only do 1st cohort to avoid reallocate array error
  esmf_loop3: do ii = 1,1 !ncohort
    do j=1,maxCountperCohort
      i=(ii-1)*maxCountperCohort+j !component #
      if(i>schismCount) exit esmf_loop3

      call ESMF_GridCompInitialize(schism_components(i), importState= importStateList(i), &
      exportState=exportStateList(i), phase=1, clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo !j
  end do esmf_loop3 !ii

  esmf_loop4: do ii = 1,ncohort
    do j=1,maxCountperCohort
      i=(ii-1)*maxCountperCohort+j !component #
      if(i>schismCount) exit esmf_loop4

      call ESMF_StateReconcile( importStateList(i), vm=vm, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      call ESMF_StateReconcile(exportStateList(i), vm=vm, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo !j
  end do esmf_loop4 !ii

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Save the time step of the current component
    call ESMF_ClockGet(clock, timeStep=timeStep, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !Advance model to next coupling step
    esmf_loop6: do ii = 1,ncohort
      do j=1,maxCountperCohort
        i=(ii-1)*maxCountperCohort+j !component #
        if(i>schismCount) exit esmf_loop6

        call ESMF_GridCompRun(schism_components(i), importState= importStateList(i), &
         exportState=exportStateList(i), clock=clock, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

        !> I am not so sure what you want to do here, but let me just play
        !> around:
        !> Reset each of these clocks
        call ESMF_GridCompGet(schism_components(i), clock=childClock, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

        call ESMF_ClockAdvance(childClock, timeStep=-timeStep, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      enddo !j

      !> call something with pdaf on this cohort ???, if not on this cohort
      !> then the whole clock resetting business should be moved to outer loop.

      do j=1,maxCountperCohort
        i=(ii-1)*maxCountperCohort+j

        !> Reset each of these clocks
        call ESMF_GridCompGet(schism_components(i), clock=childClock, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

        call ESMF_ClockAdvance(childClock, timeStep=timeStep, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      enddo !j

    end do esmf_loop6 !ii

    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do !do while

  esmf_loop5: do ii = 1,ncohort
    do j=1,maxCountperCohort
      i=(ii-1)*maxCountperCohort+j !component #
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
  integer(ESMF_KIND_I4) :: start_hour=0, runhours=2
  namelist /global/ start_year, start_month, start_day, start_hour, runhours

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
