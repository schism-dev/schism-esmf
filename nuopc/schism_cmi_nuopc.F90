! This code is part of the SCHISM-ESMF interface.  It defines
! the schism component for a NUOPC coupled system
!
! @copyright (C) 2020 Helmholtz-Zentrum Geesthacht
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

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_cmi_nuopc.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism

  use esmf
  use nuopc
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance

  use schism_bmi

  implicit none

  private
  public SetServices

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4) :: localrc

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_SetClock, &
    specRoutine=SetClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeP1"
subroutine InitializeP1(comp, importState, exportState, clock, rc)

  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)       :: localrc, mpiCommunicator, mpiCommDuplicate
  integer(ESMF_KIND_I4)       :: ntime=0, iths=0
  character(len=ESMF_MAXSTR)  :: message, compName

  type(ESMF_Vm)               :: vm

  rc = ESMF_SUCCESS

  !> @todo replace by NUOPC_CompGet()
  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Get VM for this component
  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, mpiCommunicator=mpiCommunicator, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifndef ESMF_MPIUNI
  call MPI_Comm_dup(mpiCommunicator, mpiCommDuplicate, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing parallel environment ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  call schism_parallel_init(communicator=mpiCommDuplicate)

  write(message, '(A)') trim(compName)//' initialized parallel environment.'
#endif

  write(message, '(A)') trim(compName)//' initializing science model ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_UtilIOMkDir ('./outputs',  relaxedFlag=.true., rc=localrc)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call schism_init('./', iths, ntime)
  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> @todo define more import fields, i.e wind x,y fields
  call NUOPC_Advertise(importState, &
    StandardName="air_pressure_at_sea_level", name="pmsl", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(importState, &
    StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo define more export fields, i.e wind x,y fields
  call NUOPC_Advertise(exportState, &
    StandardName="sea_surface_temperature", name="temperature_at_water_surface", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionaryAddEntry('mesh_topology','1', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionaryAddEntry('mesh_global_node_id','1', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionaryAddEntry('mesh_global_element_id','1', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionaryAddEntry('mesh_element_node_connectivity','1', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="mesh_topology", name="mesh_topology", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="mesh_global_node_id", name="mesh_global_node_id", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="mesh_global_element_id", name="mesh_global_element_id", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="mesh_element_node_connectivity", name="mesh_element_node_connectivity", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InititalizeP2"
subroutine InitializeP2(comp, importState, exportState, clock, rc)

  use schism_esmf_util, only : addSchismMesh
  !> @todo move all use statements of schism into schism_bmi
  use schism_glbl, only: np, pr2, airt2
  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_TimeInterval) :: stabilityTimeStep
  type(ESMF_Field)        :: field
  integer(ESMF_KIND_I4)   :: localrc, i, itemCount

  type(ESMF_CoordSys_Flag) :: coordsys
  type(ESMF_Mesh)          :: mesh2d

  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()

  rc = ESMF_SUCCESS

  call addSchismMesh(comp, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! call ESMF_GridCompGet(comp, exportState=exportState, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! call ESMF_StateGet(exportState,  itemCount=itemCount, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! allocate(itemTypeList(itemCount))
  ! allocate(itemNameList(itemCount))
  !
  ! call ESMF_StateGet(exportState,  itemTypeList=itemTypeList, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! do i=1, itemCount
  !   if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle
  !
  !   call ESMF_StateGet(exportState, trim(itemNameList(i)), field=field, rc=localrc )
  !   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  !   call NUOPC_Realize(exportState, field=field, rc=localrc)
  !   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ! enddo
  ! deallocate(itemTypeList)
  ! deallocate(itemNameList)

  call ESMF_GridCompGet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtr1 => pr2(1:np)
  field = ESMF_FieldCreate(name="pmsl", mesh=mesh2d, &
    farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  field = ESMF_FieldCreate(name="rsns", mesh=mesh2d, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! exportable field: sea_surface_temperature
  field = ESMF_FieldCreate(name="temperature_at_water_surface", mesh=mesh2d, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

subroutine SetClock(comp, rc)

  use schism_bmi, only : schismTimeStep

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)              :: clock
  type(ESMF_TimeInterval)       :: timeStep
  integer(ESMF_KIND_I4)         :: localrc
  real(ESMF_KIND_R8)            :: seconds

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, modelClock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call schismTimeStep(seconds)
  call ESMF_TimeIntervalSet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetClock(comp, clock, timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

subroutine ModelAdvance(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)            :: clock
  type(ESMF_State)            :: importState, exportState
  type(ESMF_Time)             :: currTime
  type(ESMF_TimeInterval)     :: timeStep
  character(len=160)          :: message, compName
  integer(ESMF_KIND_I4)       :: localrc
  integer(ESMF_KIND_I8)       :: advanceCount
  integer, save               :: it=1

  rc = ESMF_SUCCESS

  call ESMF_TraceRegionEnter("schism:ModelAdvance", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, &
    exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(clock, advanceCount=advanceCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call schism_step(it)
  it = it + 1

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

  call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing schism from: ", unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", unit=message, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TraceRegionExit("schism:ModelAdvance")
end subroutine

end module schism
