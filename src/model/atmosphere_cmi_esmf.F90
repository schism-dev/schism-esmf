! This code is part of the SCHISM-ESMF interface.  It defines
! a dummy atmosphere component for an ESMF coupled system.
!
! @copyright (C) 2017--2020-2021 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister
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
#define ESMF_FILENAME "atmosphere_cmi_esmf.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module atmosphere_cmi_esmf

  use esmf

  implicit none
  integer               :: iths=0,ntime=0

  private
  public SetServices

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

end subroutine InitializeP0

#undef  ESMF_METHOD
#define ESMF_METHOD "InitializeP1"
subroutine InitializeP1(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: clock
  integer, intent(out)  :: rc

  type(ESMF_Field)      :: field
  type(ESMF_VM)         :: vm
  type(ESMF_DistGrid)   :: distgrid
  type(ESMF_Grid)       :: grid
  type(ESMF_Array)      :: xarray,yarray,dataarray
  real(ESMF_KIND_R8), dimension(:,:),pointer :: ptr2d
  integer               :: localrc
  integer               :: i, n, petcount

  ! get communicator
  call ESMF_GridCompGet(comp,vm=vm,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm,petCount=petcount,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! create a northwestern hemisphere grid
  grid = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/100, 10/), &
    minCornerCoord=(/-120._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
    maxCornerCoord=(/60._ESMF_KIND_R8, 60._ESMF_KIND_R8/), &
    coordSys=ESMF_COORDSYS_SPH_DEG, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! define fields for export
  field = ESMF_FieldCreate(grid, name='wind_x-velocity', &
                           typekind=ESMF_TYPEKIND_R8, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldGet(field,farrayPtr=ptr2d,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ptr2d(1,1)=0.0
  ptr2d(1,2)=3.0
  ptr2d(2,1)=-3.0
  ptr2d(2,2)=1.0

  call ESMF_StateAddReplace(exportstate, (/field/), rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  rc = ESMF_SUCCESS
end subroutine InitializeP1



  subroutine Run(comp, importState, exportState, clock, rc)
  implicit none
  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: clock
  integer, intent(out)    :: rc


  rc = ESMF_SUCCESS
  end subroutine Run



  subroutine Finalize(comp, importState, exportState, clock, rc)
  implicit none
  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: clock
  type(ESMF_Grid)       :: grid
  type(ESMF_Field)      :: field
  integer, intent(out)  :: rc

  call ESMF_StateGet(importState,'wind_x-velocity',field=field,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldGet(field,grid=grid,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldDestroy(field,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_GridDestroy(grid,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  rc = ESMF_SUCCESS
  end subroutine Finalize


end module atmosphere_cmi_esmf
