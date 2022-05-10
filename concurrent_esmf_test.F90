! This code is part of the SCHISM-ESMF interface. It is a main
! program for running the SCHISM component concurrent to a
! dummy component (in src/model/atmosphere_cmi_esmf.F90).
!
! @copyright (C) 2018, 2019, 2020-2021 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister richard.hofmeister@hereon.de
! @author Carsten Lemmen carsten.lemmen@hereon.de
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
#define ESMF_FILENAME "concurrent_esmf_test.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
  use atmosphere_cmi_esmf,  only: atmosSetServices => SetServices

  implicit none

  interface
    function clockCreateFrmParam(filename, rc)
      use esmf
        character(len=ESMF_MAXSTR), intent(in) :: filename
        integer(ESMF_KIND_I4), intent(out)     :: rc
        type(ESMF_Clock)                       :: clockCreateFrmParam
    end function clockCreateFrmParam
  end interface

  type(ESMF_GridComp)     :: schism_component
  type(ESMF_GridComp)     :: atmos_component

  type(ESMF_State)        :: schism_import, schism_export
  type(ESMF_State)        :: atmos_import, atmos_export

  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time
  type(ESMF_Clock)        :: clock

  type(ESMF_Field)        :: field,field_in,field_out
  type(ESMF_RouteHandle)  :: routehandle_air2sea, routehandle_sea2air
  type(ESMF_Vm)           :: vm

  integer(ESMF_KIND_I4)       :: rc, petCount, i, inum, localrc
  integer, allocatable        :: petlist_schism(:),petlist_atmos(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! all but one PET for SCHISM, and last core for the
  ! dummy atmosphere component

  call ESMF_VMGetGlobal(vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(petlist_schism(max(1,petcount-1)), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1, max(1, petcount-1)
    petlist_schism(i)=i-1 !0 based
  end do
  allocate(petlist_atmos(1))
  petlist_atmos(1) = petcount-1

  ! Create both components on their respective parallel
  ! environment provided by each petList, then register
  ! the components entry points.

  schism_component = ESMF_GridCompCreate(name='schism_component', &
    petList=petlist_schism, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  atmos_component = ESMF_GridCompCreate(name='atmosphere_component', &
    petList=petlist_atmos, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(schism_component, &
    schismSetServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(atmos_component, &
    atmosSetServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create states for exchange of information between
  ! components.
  schism_export = ESMF_StateCreate(name='schism export state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  schism_import = ESMF_StateCreate(name='schism import state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  atmos_export = ESMF_StateCreate(name='atmos export state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  atmos_import = ESMF_StateCreate(name='atmos import state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create a clock, this should ideally be done with a configuration
  ! file but is hardcoded here for now. The timestep defined here
  ! is the coupling timestep.

  ! call ESMF_TimeIntervalSet(timestep, s=3600, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! call ESMF_TimeSet(start_time, yy=2018, mm=4, dd=10, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! stop_time = start_time + 2*timestep
  !
  ! clock = ESMF_ClockCreate(timestep, start_time, stopTime=stop_time, &
  !   name='main clock',rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  filename = './global.nml'
  clock = clockCreateFrmParam(filename, localrc)

  ! Initialize SCHISM
  call ESMF_GridCompInitialize(schism_component, &
    importState=schism_import, &
    exportState=schism_export, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompInitialize(atmos_component, importState=atmos_import, &
    exportState=atmos_export, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Make sure that all states are reconciled across the
  ! entire VM
  call ESMF_StateReconcile(schism_import, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateReconcile(atmos_export, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Within the schism component, the following fields are
  ! defined for import and export (y-component not implemented yet)
  call ESMF_StateGet(schism_import, 'wind_x-velocity_in_10m_height', &
    field=field_out, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(atmos_export,'wind_x-velocity', &
    field=field_in, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Precompute the weights
  call ESMF_FieldRegridStore(field_in, field_out, &
    regridMethod =  ESMF_REGRIDMETHOD_BILINEAR, &
    routehandle=routehandle_air2sea, unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Directly manipulate the fields from import and export
    !> states.  In less basic applications, this should be
    !> handled by a mediator component.
    call ESMF_FieldRegrid(field_in, field_out, &
      routeHandle=routehandle_air2sea, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompRun(atmos_component, importState=atmos_import, &
      exportState=atmos_export, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompRun(schism_component, importState=schism_import, &
      exportState=schism_export, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do

  !> Clean up
  call ESMF_GridCompFinalize(schism_component, importState=schism_import, &
    exportState=schism_export, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(schism_import, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(schism_export, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(atmos_import, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(atmos_export, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompDestroy(schism_component, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompDestroy(atmos_component, rc=localrc)
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
