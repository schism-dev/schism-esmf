! This code is a main driver
! program for running multiple schism components concurrently 
!
! @copyright (C) 2018, 2019, 2020 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister richard.hofmeister@hzg.de
! @author Carsten Lemmen carsten.lemmen@hzg.de
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
#define ESMF_FILENAME "multi_schism.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
!  use atmosphere_cmi_esmf,  only: atmosSetServices => SetServices

  implicit none

  interface
    function clockCreateFrmParam(filename, rc)
      use esmf
        character(len=ESMF_MAXSTR), intent(in) :: filename
        integer(ESMF_KIND_I4), intent(out)     :: rc
        type(ESMF_Clock)                       :: clockCreateFrmParam
    end function clockCreateFrmParam
  end interface

  type(ESMF_GridComp), allocatable :: schism_components(:)
!  type(ESMF_GridComp)     :: atmos_component

  type(ESMF_State), allocatable   :: schism_imports(:), schism_exports(:)
!  type(ESMF_State)        :: atmos_import, atmos_export

  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time
  type(ESMF_Clock)        :: clock

  type(ESMF_Field)               :: field, field_in
  type(ESMF_Field), allocatable  :: fields_out(:)

  type(ESMF_RouteHandle)  :: routehandle_sea2air
  type(ESMF_RouteHandle), allocatable  :: routehandles_air2sea(:)
  type(ESMF_Vm)           :: vm

  integer(ESMF_KIND_I4)       :: petCountLocal, schismCount=8
  integer(ESMF_KIND_I4)       :: rc, petCount, i, j, inum, localrc
  integer(ESMF_KIND_I4), allocatable    :: petlist(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename, message

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! PETs for SCHISM

  call ESMF_VMGetGlobal(vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create all components on their respective parallel
  ! environment provided by each petList
  if(schismCount>999.or.schismCount<1) then
!Error: better way to crash out?
    write(*,*)'Check # of components:',schismCount 
    call ESMF_Finalize(rc=localrc)
    !call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO,rc=localrc)
    !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  allocate(schism_components(schismCount))
  do i = 1, schismCount

    petCountLocal = petCount/schismCount
    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do j=1, petCountLocal
      petList(j)=(i-1)*petCountLocal + j - 1 ! 0 based
    end do


    write(message, '(A,I3.3)') 'schism_', i
    !write(0, *) trim(message), 'list=', petList, 'petCount=', petCount, petCountLocal

    schism_components(i) = ESMF_GridCompCreate(name=trim(adjustl(message)), &
      petList=petlist, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    deallocate(petList)
  end do !i

  allocate(schism_exports(schismCount))
  allocate(schism_imports(schismCount))

  do i = 1, schismCount

    call ESMF_GridCompSetServices(schism_components(i), &
      schismSetServices, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    schism_exports(i) = ESMF_StateCreate(name='schism export state', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    schism_imports(i) = ESMF_StateCreate(name='schism import state', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

!  atmos_component = ESMF_GridCompCreate(name='atmosphere_1', &
!    petList=(/petCount - 1/), rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  atmos_export = ESMF_StateCreate(name='atmos export state', rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  atmos_import = ESMF_StateCreate(name='atmos import state', rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  call ESMF_GridCompSetServices(atmos_component, &
!    atmosSetServices, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  filename = './global.nml'
  clock = clockCreateFrmParam(filename, localrc)

  do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

!  call ESMF_GridCompInitialize(atmos_component, importState=atmos_import, &
!    exportState=atmos_export, clock=clock, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Make sure that all states are reconciled across the
  ! entire VM
  do i = 1, schismCount

    call ESMF_StateReconcile(schism_imports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateReconcile(schism_exports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

!  call ESMF_StateReconcile(atmos_export, vm=vm, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Within the schism component, the following fields are
  ! defined for import and export (y-component not implemented yet)
!  call ESMF_StateGet(atmos_export,'wind_x-velocity', &
!    field=field_in, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  allocate(fields_out(schismCount))
!  allocate(routehandles_air2sea(schismCount))
!
!  do i = 1, schismCount
!    call ESMF_StateGet(schism_imports(i), 'wind_x-velocity_in_10m_height', &
!      field=fields_out(i), rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!    ! Precompute the weights
!    call ESMF_FieldRegridStore(field_in, fields_out(i), &
!      routehandle=routehandles_air2sea(i), rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  enddo

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Directly manipulate the fields from import and export
    !> states.  In less basic applications, this should be
    !> handled by a mediator component.
!    do i=1, schismCount
!      call ESMF_FieldRegrid(field_in, fields_out(i), &
!          routeHandle=routehandles_air2sea(i), rc=localrc)
!        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!    enddo

!    call ESMF_GridCompRun(atmos_component, importState=atmos_import, &
!      exportState=atmos_export, clock=clock, rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do i=1, schismCount
      call ESMF_GridCompRun(schism_components(i), importState=schism_imports(i), &
        exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo

    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do

  ! Clean up

!  deallocate(fields_out)
!  deallocate(routehandles_air2sea)

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

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  call ESMF_StateDestroy(atmos_import, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  call ESMF_StateDestroy(atmos_export, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  call ESMF_GridCompDestroy(atmos_component, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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
  integer(ESMF_KIND_I4) :: start_hour=0, rnday=2 !rnday in hours
  namelist /global/ start_year, start_month, start_day, start_hour, rnday

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

  call ESMF_TimeIntervalSet(timeStep, h=rnday, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  stopTime = startTime + timeStep

  ! Only now define the coupling timestep as fraction of full timeStep (24
  ! coupling steps in total here)
  timeStep = timeStep / 24

  clock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, &
    name='main clock', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end function clockCreateFrmParam
