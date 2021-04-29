! This code is part of the SCHISM-ESMF interface.  It defines
! a dummy atmosphere component for a NUOPC coupled system.
!
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
! @author Richard Hofmeister
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

module atmosphere

  use esmf
  use nuopc

  ! This is a dummy science model, so derive it from NUOPC model
  use NUOPC_Model, &
    model_routine_SS    => SetServices, &
    model_label_Advance => label_Advance

  implicit none

  private
  public SetServices

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4) :: localRc

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeP1"
subroutine InitializeP1(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4) :: localRc

  rc = ESMF_SUCCESS

  ! This component imports SST and exports SLP and SWFLUX

  !call NUOPC_Advertise(importState, &
  !  StandardName="surface_temperature", name="sst", rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="surface_air_pressure", name="pmsl", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
      StandardName="surface_downwelling_photosynthetic_radiative_flux", name="rsns", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo expand fieldDict to correctly label this
  call NUOPC_Advertise(exportState, &
      StandardName="x_velocity_at_10m_above_sea_surface", name="x_velocity_at_10m_above_sea_surface", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
      StandardName="y_velocity_at_10m_above_sea_surface", name="y_velocity_at_10m_above_sea_surface", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine InitializeP1

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeP2"
subroutine InitializeP2(comp, importState, exportState, clock, rc)

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_Field)        :: field
  type(ESMF_Grid)         :: gridIn, gridOut
  integer(ESMF_KIND_I4)   :: localRc, ubnd(2), lbnd(2), i, j

  real(ESMF_KIND_R8), pointer  :: farrayPtr2(:,:) => null()

  rc = ESMF_SUCCESS

  ! create a Grid object for Fields
  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 100/), &
    minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
    maxCornerCoord=(/2000._ESMF_KIND_R8, 20000._ESMF_KIND_R8/), &
    coordSys=ESMF_COORDSYS_CART, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  gridOut = gridIn ! for now out same as in

  field = ESMF_FieldCreate(name="sst", grid=gridIn, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Disabled for now as we disabled coupling OCN->ATM
  !call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifdef CREATE_AND_REALIZE
  ! This branch shows the standard procedure of creating a complete field
  ! with Grid and memory allocation, and then calling Realize() for it.

  field = ESMF_FieldCreate(name="pmsl", grid=gridOut, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#else
  ! This branch shows the alternative way of "realizing" an advertised
  ! field, it accesses the empty field that was created during advertise,
  ! and finishes it, setting a grid, and then calling FieldEmptyComplete().

  call ESMF_StateGet(exportState, field=field, itemName="pmsl", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldEmptySet(field, grid=gridOut, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! There is not need to formally call Realize() when completing the
  ! adverised field directly. However, calling Realize() also works.
  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
#endif

  ! Other export fields
  field = ESMF_FieldCreate(name="rsns", grid=gridOut, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  field = ESMF_FieldCreate(name="x_velocity_at_10m_above_sea_surface", grid=gridOut, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr2, &
    exclusiveUBound=ubnd, exclusiveLBound=lbnd, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do j=lbnd(2),ubnd(2)
    do i=lbnd(1), ubnd(1)
      farrayPtr2(i,j) = 3*sin(i/ubnd(1)*3.14/180.) &
        + 3*cos(j/ubnd(2)*3.14/180.)
    end do
  enddo

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  field = ESMF_FieldCreate(name="y_velocity_at_10m_above_sea_surface", grid=gridOut, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "ModelAdvance"
subroutine ModelAdvance(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    character(len=ESMF_MAXSTR)  :: message, name
    integer(ESMF_KIND_I4) :: localRc

    rc = ESMF_SUCCESS

    call ESMF_TraceRegionEnter("atmosphere:ModelAdvance", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, &
      exportState=exportState, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompGet(comp, name=name, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ClockPrint(clock, options="currTime", &
      preString=trim(name)//" advances from ", unit=message, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_ClockPrint(clock, options="stopTime", &
      preString=trim(name)//" advances to ", unit=message, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_TraceRegionExit("atmosphere:ModelAdvance", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end subroutine

end module atmosphere
