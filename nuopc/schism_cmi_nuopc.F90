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

  call NUOPC_Advertise(importState, &
    StandardName="air_pressure_at_sea_level", name="pmsl", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(importState, &
    StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="sea_surface_temperature", name="sst", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InititalizeP2"
subroutine InitializeP2(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_TimeInterval) :: stabilityTimeStep
  type(ESMF_Field)        :: field
  type(ESMF_Grid)         :: gridIn
  type(ESMF_Grid)         :: gridOut
  integer(ESMF_KIND_I4) :: localrc

  rc = ESMF_SUCCESS

  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/100, 10/), &
      minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
      maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
      coordSys=ESMF_COORDSYS_CART, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
      rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  gridOut = gridIn ! for now out same as in

  field = ESMF_FieldCreate(name="pmsl", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  field = ESMF_FieldCreate(name="rsns", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Realize(importState, field=field, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! exportable field: sea_surface_temperature
    field = ESMF_FieldCreate(name="sst", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    call NUOPC_Realize(exportState, field=field, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

subroutine SetClock(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)              :: clock
  type(ESMF_TimeInterval)       :: stabilityTimeStep
  integer(ESMF_KIND_I4) :: localrc

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, modelClock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_TimeIntervalSet(stabilityTimeStep, m=5, rc=localrc) ! 5 minute steps
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetClock(comp, clock, stabilityTimeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

subroutine ModelAdvance(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)            :: clock
  type(ESMF_State)            :: importState, exportState
  type(ESMF_Time)             :: currTime
  type(ESMF_TimeInterval)     :: timeStep
  character(len=160)          :: msgString
  integer(ESMF_KIND_I4) :: localrc

  call ESMF_TraceRegionEnter("schism:ModelAdvance")

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, &
    exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

  call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing schism from: ", unit=msgString, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", unit=msgString, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TraceRegionExit("schism:ModelAdvance")
end subroutine

end module schism
