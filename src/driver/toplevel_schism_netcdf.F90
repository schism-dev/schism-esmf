! This code is part of the SCHISM-ESMF interface. It is a main
! program for running the SCHISM component sequential with netcdf output
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
! Unless required by applicable law or agreed to in writESMF_Componenting, software
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
#define ESMF_FILENAME "toplevel_schism_netcdf.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module toplevel_schism_netcdf

  use esmf
  use schism_cmi_esmf, only: schismSetServices => SetServices
  use netcdf_component,  only: netcdfSetServices => SetServices

  implicit none

  public SetServices

  !> @todo this needs to be moved to internal state
  type(ESMF_GridComp)    :: netcdfComponent, schismComponent

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

  integer(ESMF_KIND_I4)   :: petCount, i, inum
  integer, allocatable    :: petList(:)
  type(ESMF_State)        :: schismExport, schismImport, netcdfExport, netcdfImport

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VmGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(petList(petCount))
  do i=1, petCount
    petList(i) = i-1
  enddo

  schismComponent = ESMF_GridCompCreate(name='schism', &
    petList=petList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  netcdfComponent = ESMF_GridCompCreate(name='netcdf', &
    petList=petList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(schismComponent, &
    schismSetServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(netcdfComponent, &
    netcdfSetServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create states for exchange of information between
  ! components.
  schismExport = ESMF_StateCreate(name='schismExport', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  schismImport = ESMF_StateCreate(name='schismImport', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  netcdfExport = ESMF_StateCreate(name='netcdfExport', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  netcdfImport = ESMF_StateCreate(name='netcdfImport', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompInitialize(schismComponent, &
    importState=schismImport, &
    exportState=schismExport, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompInitialize(netcdfComponent, importState=netcdfImport, &
    exportState=netcdfExport, clock=clock, rc=localrc)
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
  type(ESMF_State)        :: schismExport, schismImport, netcdfExport, netcdfImport

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(schismComponent, exportState=schismExport, &
    importState=schismImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(netcdfComponent, exportState=netcdfExport, &
    importState=netcdfImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Directly manipulate the fields from import and exportschism_component
    !> states.  In less basic applications, this should be
    !> handled by a mediator component.
    call ESMF_GridCompRun(schismComponent, importState=schismImport, &
      exportState=schismExport, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompRun(netcdfComponent, importState=schismExport, &
      exportState=netcdfExport, clock=clock, rc=localrc)
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
  type(ESMF_State)        :: schismExport, schismImport, netcdfExport, netcdfImport

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(schismComponent, exportState=schismExport, &
    importState=schismImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(netcdfComponent, exportState=netcdfExport, &
    importState=netcdfImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Clean up
  call ESMF_GridCompFinalize(schismComponent, importState=schismImport, &
    exportState=schismExport, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompFinalize(netcdfComponent, importState=netcdfImport, &
    exportState=netcdfExport, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(schismImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(schismExport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(netcdfImport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(netcdfExport, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompDestroy(schismComponent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompDestroy(netcdfComponent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine Finalize
end module toplevel_schism_netcdf
