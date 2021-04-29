! This code is part of the SCHISM-ESMF interface.  It defines
! a driver for a NUOPC coupled system of SCHISM with a
! dummy atmosphere without a mediator.
!
! @copyright (C) 2021 Helmholtz-Zentrum Hereon
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
!
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

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "driver_schism_atm.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "driver_schism_atm"
module driver_schism_atm

  use esmf
  use nuopc

  ! We specialize this toplevel as a NUOPC Driver instance
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices, &
    driver_label_SetModelServices => label_SetModelServices

  use atmosphere_cmi_nuopc, only: atmosphereSS => SetServices
  use schism_cmi_nuopc, only: schismSS => SetServices

  use NUOPC_Connector, only: cplSS => SetServices

  implicit none

  private

  public SetServices

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(driver, rc)

  type(ESMF_GridComp)  :: driver
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4) :: localrc

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(driver, driver_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
    specRoutine=SetModelServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(driver, name="Verbosity", value="high", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetModelServices"
subroutine SetModelServices(driver, rc)

  type(ESMF_GridComp)  :: driver
  integer, intent(out) :: rc

  type(ESMF_Grid)               :: grid
  type(ESMF_Field)              :: field
  type(ESMF_Time)               :: startTime, stopTime
  type(ESMF_TimeInterval)       :: timeStep
  type(ESMF_Clock)              :: internalClock
  type(ESMF_GridComp)           :: child
  type(ESMF_CplComp)            :: connector
  integer(ESMF_KIND_I4)         :: localrc, petCount, slotCount, i
  type(ESMF_Vm)                 :: vm
  integer(ESMF_KIND_I4), allocatable :: petList(:)
  type(NUOPC_FreeFormat)        :: freeFormat

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(driver, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VmGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Use all but one pets in this petList for schism and the last
  ! PET for the dummy atmosphere.  This runs the models concurrently.
  ! For petCount=1, both run sequentially on the same PET
  allocate(petList(max(1, petCount-1)))
  do i=1, max(1, petCount-1)
    petList(i) = i-1
  enddo

  ! SetServices for dummy atmosphereosphere and for schism ocean
  ! NUOPC_DriverAddGridComp(driver, compLabel, &
  !   compSetServicesRoutine, compSetVMRoutine, petList, info, driver, rc)
  !> @todo this next statement throws an error without showing this in the log,
  !> but simply exiting to error at main:97
  call NUOPC_DriverAddComp(driver=driver, compLabel="atmosphere", &
    compSetServicesroutine=atmosphereSS, petList=(/petCount-1/), comp=child, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(child, name="Verbosity", value="low", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_DriverAddComp(driver, "schism", schismSS, petList=petList, &
    comp=child, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Doesn't seem to work for now
  !call NUOPC_GridCompAttributeEgest(child, freeFormat, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call NUOPC_FreeFormatLog(freeFormat, iofmt=ESMF_IOFMT_YAML, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! SetServices for connectors atmosphere--ocean and vice versa
  call NUOPC_DriverAddComp(driver, srcCompLabel="atmosphere",   dstCompLabel="schism", &
    compSetServicesRoutine=cplSS, comp=connector, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(connector, name="Verbosity", value="high", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Let's skip the feedback for now
  ! call NUOPC_DriverAddComp(driver, srcCompLabel="schism", dstCompLabel="atmosphere", &
  !   compSetServicesRoutine=cplSS, comp=connector, rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! call NUOPC_CompAttributeSet(connector, name="Verbosity", value="high", rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  deallocate(petList)

  ! call NUOPC_DriverGet(driver, slotCount=slotCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_DriverEgestRunSequence(driver, freeFormat=freeFormat, &
     rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FreeFormatLog(freeFormat, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionaryEgest(freeFormat, iofmt=ESMF_IOFMT_YAML, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FreeFormatLog(freeFormat, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FreeFormatDestroy(freeFormat, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

!NUOPC_StateWrite(state, fieldNameList, fileNamePrefix, overwrite, &
!     status, timeslice, relaxedflag, rc)


end module driver_schism_atm
