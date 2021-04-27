! This code is part of the SCHISM-ESMF interface. It is a main
! program for running the SCHISM component concurrent to a
! dummy component.
!
! @copyright (C) 2018-2020-2021 Helmholtz-Zentrum Geesthacht
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
#define ESMF_FILENAME "toplevel_schism_atm.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module toplevel_schism_atm

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
  use atmosphere_cmi_esmf,  only: atmosSetServices => SetServices

  implicit none

  public SetServices

  !> @todo this needs to be moved to internal state or dedicated coupler
  type(ESMF_RouteHandle) :: routehandle_air2sea, routehandle_sea2air
  type(ESMF_Field)       :: field_in, field_out
  type(ESMF_GridComp)    :: atmos_component, schism_component
  type(ESMF_State)       :: schism_import, schism_export
  type(ESMF_State)       :: atmos_import, atmos_export

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)                :: comp
  integer(ESMF_KIND_I4), intent(out) :: rc

  integer(ESMF_KIND_I4)              :: localrc

  rc = ESMF_SUCCESS

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phase=0, userRoutine=InitializeP0, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phase=1, userRoutine=InitializeP1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
    userRoutine=Run, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
    userRoutine=Finalize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef  ESMF_METHOD
#define ESMF_METHOD "InitializeP0"
  subroutine InitializeP0(gridComp, importState, exportState, &
    parentClock, rc)

    implicit none

    type(ESMF_GridComp)         :: gridComp
    type(ESMF_State)            :: importState
    type(ESMF_State)            :: exportState
    type(ESMF_Clock)            :: parentClock
    integer, intent(out)        :: rc

    character(len=10)           :: InitializePhaseMap(1)
    character(len=ESMF_MAXSTR)  :: myName
    type(ESMF_Time)             :: currTime
    integer                     :: localrc

    rc=ESMF_SUCCESS

    InitializePhaseMap(1) = "IPDv00p1=1"

    call ESMF_AttributeAdd(gridComp, convention="NUOPC", &
      purpose="General", &
      attrList=(/"InitializePhaseMap"/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_AttributeSet(gridComp, name="InitializePhaseMap", valueList=InitializePhaseMap, &
      convention="NUOPC", purpose="General", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateValidate(importState, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateValidate(exportState, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end subroutine InitializeP0

#undef  ESMF_METHOD
#define ESMF_METHOD "InitializeP1"
subroutine InitializeP1(comp, importState, exportState, clock, rc)


  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: clock
  integer, intent(out)    :: rc

  integer(ESMF_KIND_I4)   :: localrc
  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time

  type(ESMF_Field)        :: field
  type(ESMF_Vm)           :: vm

  integer(ESMF_KIND_I4)       :: petCount, i, inum
  integer, allocatable        :: petList(:)


  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VmGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Use all but one pets in this petList for schism and the last
  ! PET for the dummy atmosphere.  This runs the models concurrently.
  ! For petCount=1, both run sequentially on the same PET  type(ESMF_Field)        :: field, field_in, field_out

  allocate(petList(max(1, petCount-1)))
  do i=1, max(1, petCount-1)
    petList(i) = i-1
  enddo

  ! Create both components on their respective parallel
  ! environment provided by each petList, then register
  ! the components entry points.

  schism_component = ESMF_GridCompCreate(name='schismComponent', &
    petList=petList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  deallocate(petList)

  atmos_component = ESMF_GridCompCreate(name='atmosComponent', &
    petList=(/petCount-1/), rc=localrc)
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
  ! defined for import and export
  call ESMF_StateGet(schism_import, 'wind_x-velocity_in_10m_height', &
    field=field_out, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(atmos_export,'wind_x-velocity', &
    field=field_in, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Precompute the weights
  call ESMF_FieldRegridStore(field_in, field_out, &
    routehandle=routehandle_air2sea, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine InitializeP1

#undef ESMF_METHOD
#define ESMF_METHOD "Run"
subroutine Run(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: clock
  integer, intent(out)    :: rc

  integer(ESMF_KIND_I4)   :: localrc


  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Directly manipulate the fields from import and export
    !> states.  In less basic applications, this should be
    !> handled by a mediator component.
    call ESMF_FieldRegrid(field_in, field_out, &
      routeHandle=routehandle_air2sea,rc=localrc)
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
end subroutine Run

#undef ESMF_METHOD
#define ESMF_METHOD "Finalize"
subroutine Finalize(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: clock
  integer, intent(out)    :: rc

  integer(ESMF_KIND_I4)   :: localrc

  !> Clean up
  call ESMF_GridCompFinalize(schism_component, importState=schism_import, &
    exportState=schism_export, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompFinalize(atmos_component, importState=atmos_import, &
    exportState=atmos_export, clock=clock, rc=localrc)
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

end subroutine Finalize
end module toplevel_schism_atm
