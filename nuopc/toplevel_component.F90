! This code is part of the SCHISM-ESMF interface.  It defines
! a toplevel component for a NUOPC coupled system of SCHISM with a
! dummy atmosphere without a mediator.
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
#define ESMF_FILENAME "toplevel_component.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "toplevel"
module toplevel

  use esmf
  use nuopc

  ! We specialize this toplevel as a NUOPC Driver instance
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices, &
    driver_label_SetModelServices => label_SetModelServices

  use atmosphere, only: atmosphereSS => SetServices
  use schism, only: schismSS => SetServices

  use NUOPC_Connector, only: cplSS => SetServices

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

  call NUOPC_CompDerive(comp, driver_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=driver_label_SetModelServices, &
    specRoutine=SetModelServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(comp, name="Verbosity", value="high", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetModelServices"
subroutine SetModelServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Grid)               :: grid
  type(ESMF_Field)              :: field
  type(ESMF_Time)               :: startTime, stopTime
  type(ESMF_TimeInterval)       :: timeStep
  type(ESMF_Clock)              :: internalClock
  type(ESMF_GridComp)           :: child
  type(ESMF_CplComp)            :: connector
  integer(ESMF_KIND_I4)         :: localrc

  rc = ESMF_SUCCESS

  ! SetServices for dummy atmosphereosphere and for schism ocean
  call NUOPC_DriverAddComp(comp, "atmosphere", atmosphereSS, comp=child, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_DriverAddComp(comp, "schism", schismSS, comp=child, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_CompAttributeSet(child, name="Verbosity", value="high", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! SetServices for connectores atmosphere--ocean and vice versa
    call NUOPC_DriverAddComp(comp, srcCompLabel="atmosphere", dstCompLabel="schism", &
      compSetServicesRoutine=cplSS, comp=connector, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_CompAttributeSet(connector, name="Verbosity", value="high", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_DriverAddComp(comp, srcCompLabel="schism", dstCompLabel="atmosphere", &
      compSetServicesRoutine=cplSS, comp=connector, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_CompAttributeSet(connector, name="Verbosity", value="high", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TimeIntervalSet(timeStep, m=15, rc=localrc) ! 15 minute steps
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    !> @todo inherit from main()
    call ESMF_TimeSet(startTime, yy=2010, mm=6, dd=1, h=0, m=0, &
      calkindflag=ESMF_CALKIND_GREGORIAN, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TimeSet(stopTime, yy=2010, mm=6, dd=1, h=1, m=0, &
      calkindflag=ESMF_CALKIND_GREGORIAN, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    internalClock = ESMF_ClockCreate(name="Application Clock", &
      timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompSet(comp, clock=internalClock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end subroutine

end module toplevel
