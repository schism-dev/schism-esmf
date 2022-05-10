! This code is a main driver
! program for running multiple schism components concurrently
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

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "multi_schism.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
  use schism_esmf_util, only: clockCreateFrmParam

  implicit none

  type(ESMF_GridComp), allocatable :: schism_components(:)
  type(ESMF_State), allocatable    :: schism_imports(:), schism_exports(:)

  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time
  type(ESMF_Clock)        :: clock

  type(ESMF_Vm)           :: vm
  type(ESMF_Log)          :: log

  integer(ESMF_KIND_I4)       :: petCountLocal, schismCount=8
  integer(ESMF_KIND_I4)       :: rc, petCount, i, j, inum, localrc
  integer(ESMF_KIND_I4), allocatable    :: petlist(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename='multi_schism.cfg', message
  type(ESMF_Config)           :: config
  type(ESMF_Config), allocatable :: configList(:)

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogSet(flush=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! PETs for SCHISM.

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

    call ESMF_ConfigGetAttribute(config, value=schismCount, label='count:', &
      !default: value used if label is not found
      default=min(petCount,8), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if(schismCount>999 .or. schismCount<1) then
    write(message, '(A,I3,A)') 'Number of instances ',schismCount, &
      'must be in the range [1,999]'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  elseif (schismCount > petCount) then
    write(message, '(A)') 'Requested number of instances exceeds available PETs'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    write(message, '(A,I3,A,I3)') 'Reduced from ',schismCount, ' to ', petCount
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    schismCount = petCount
  else
    write(message, '(A,I3)') 'Number of SCHISM instances ',schismCount
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Create all components on their respective parallel
  ! environment provided by each petList
  allocate(schism_components(schismCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(configList(schismCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i = 1, schismCount

    petCountLocal = petCount/schismCount
    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do j=1, petCountLocal
      petList(j)=(i-1)*petCountLocal + j - 1 ! PET #; 0 based
    end do

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

  allocate(schism_exports(schismCount))
  allocate(schism_imports(schismCount))

  do i = 1, schismCount

    call ESMF_GridCompSetServices(schism_components(i), &
      schismSetServices, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A,I3.3)') 'schism_', i

    schism_exports(i) = ESMF_StateCreate( &
      name='schismExport_'//trim(adjustl(message)), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    schism_imports(i) = ESMF_StateCreate( &
      name='schismImport_'//trim(adjustl(message)), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  do i = 1, schismCount
    !Init phase 0
    call ESMF_GridCompInitialize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), phase=0, clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  !Get info on simulation period
  filename = './global.nml'
  clock = clockCreateFrmParam(filename, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), phase=1, clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  do i = 1, schismCount

    call ESMF_StateReconcile(schism_imports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateReconcile(schism_exports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    do i=1, schismCount
      call ESMF_GridCompRun(schism_components(i), importState=schism_imports(i), &
        exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo

    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do

  do i = 1, schismCount

    call ESMF_GridCompFinalize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateDestroy(schism_imports(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateDestroy(schism_exports(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompDestroy(schism_components(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do

  deallocate(schism_exports)
  deallocate(schism_imports)
  deallocate(schism_components)
  deallocate(configList)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize(rc=localrc)

end program main
