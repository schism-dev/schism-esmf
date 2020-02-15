! This code is part of the SCHISM-ESMF interface.  It defines
! a main() program for a NUOPC coupled system.
!
! @copyright (C) 2020 Helmholtz-Zentrum Geesthacht
! @author Carsten Lemmen <carsten.lemmen@hzg.de>
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
#define ESMF_FILENAME "main_nuopc.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use nuopc

  use toplevel, only: toplevelSetServices => SetServices

  implicit none

  integer                 :: localrc, userRc, rc
  type(ESMF_GridComp)     :: topComp

  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldDictionarySetup("field_dictionary.yaml", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! @todo find out why this does not work
  !call FieldDictionaryLog("field_dictionary.yaml.out", iofmt=ESMF_IOFMT_YAML, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create the top level component and register its services along with
  ! profiling attributes
  topComp = ESMF_GridCompCreate(name="toplevel", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(topComp, toplevelSetServices, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  call NUOPC_CompAttributeSet(topComp, name="Profiling", value="0", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Initialize, run, and finalize top level
  call ESMF_GridCompInitialize(topComp, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  call ESMF_GridCompRun(topComp, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  call ESMF_GridCompFinalize(topComp, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  ! Clean up and finish
  call ESMF_GridCompDestroy(topComp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize()

end program main
