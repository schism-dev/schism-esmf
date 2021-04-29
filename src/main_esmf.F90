! This code is part of the SCHISM-ESMF interface.  It defines
! a main() program for an ESMF coupled system.
!
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
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
#define ESMF_FILENAME "main_esmf.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#undef ESMF_METHOD
#define ESMF_METHOD "main"
program main

  use esmf
  !use toplevel_schism_atm, only: toplevelSetServices => SetServices
  use toplevel_schism_netcdf, only: toplevelSetServices => SetServices
  use schism_esmf_util, only: clockCreateFrmParam

  implicit none

  integer                 :: localrc, userRc, rc
  type(ESMF_GridComp)     :: topComp
  character(len=ESMF_MAXSTR)  :: filename
  type(ESMF_Clock)        :: clock

  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
    defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogSet(flush=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create the top level component and register its services
  topComp = ESMF_GridCompCreate(name="toplevel", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(topComp, toplevelSetServices, userRc=userRc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(rc)

  !> Read clock parameters from global.nml

  filename = './global.nml'
  clock = clockCreateFrmParam(filename, relaxedFlag=.true., rc=localrc)

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
