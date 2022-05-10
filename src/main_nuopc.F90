! This code is part of the SCHISM-ESMF interface.  It defines a main() program
! for a NUOPC coupled system.
!
! @copyright (C) 2021 Helmholtz-Zentrum Hereon
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
!
! @author Carsten Lemmen <carsten.lemmen@hzg.de>
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
#define ESMF_FILENAME "main_nuopc.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use nuopc
  use driver_schism_atm, only: driverSetServices => SetServices
  use schism_esmf_util, only : clockCreateFrmParam

  implicit none

  integer                     :: localrc, userRc, rc
  type(ESMF_GridComp)         :: topComp
  character(len=ESMF_MAXSTR)  :: filename
  type(ESMF_Clock)            :: clock
  character(len=ESMF_MAXSTR)  :: message, string
  type(NUOPC_FreeFormat)      :: freeFormat

  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
    defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo sugarglider crashes on the call with filename argument
  !call NUOPC_FieldDictionarySetup("field_dictionary.yaml", rc=localrc)
  call NUOPC_FieldDictionarySetup(rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !
  ! @todo find out why this does not work
  call NUOPC_FieldDictionaryEgest(freeFormat, iofmt=ESMF_IOFMT_YAML, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FreeFormatLog(freeFormat, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Read clock parameters from global.nml, be relaxed about this
  filename = './global.nml'
  clock = clockCreateFrmParam(filename, relaxedFlag=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create the top level component and register its services along with
  ! profiling attributes
  topComp = ESMF_GridCompCreate(name="toplevel", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(topComp, driverSetServices, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)


  call NUOPC_CompAttributeSet(topComp, name="Profiling", value="0", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockPrint(clock, options="startTime", &
      preString='main starts at ', unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockPrint(clock, options="stopTime", &
      preString='main will stop at ', unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Initialize, run, and finalize top level
  call ESMF_GridCompInitialize(topComp, clock=clock, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  call ESMF_GridCompRun(topComp, clock=clock, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  call ESMF_GridCompFinalize(topComp, clock=clock, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  ! Clean up and finish
  call ESMF_GridCompDestroy(topComp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize()

end program main
