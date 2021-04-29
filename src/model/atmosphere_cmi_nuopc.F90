! This code is part of the SCHISM-ESMF interface.  It defines
! a dummy atmosphere component for a NUOPC coupled system.
!
! @copyright (C) 2021 Helmholtz-Zentrum Hereon
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
!
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
#define ESMF_FILENAME "atmosphere_cmi_nuopc.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module atmosphere_cmi_nuopc

  use esmf
  use nuopc
  use schism_esmf_util, only : SCHISM_FieldRealize

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

  integer(ESMF_KIND_I4)      :: localRc
  character(len=ESMF_MAXSTR) :: message, compName

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (.not.ESMF_StateIsCreated(importState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Import', stateintent= &
      ESMF_STATEINTENT_IMPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if (.not.ESMF_StateIsCreated(exportState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Export', stateintent= &
      ESMF_STATEINTENT_EXPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! This component imports SST and exports SLP and SWFLUX

  if (.not.NUOPC_FieldDictionaryHasEntry("surface_air_pressure", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("surface_air_pressure", "N m-2", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="surface_air_pressure", name="air_pressure_at_water_surface", &
    SharePolicyField='share', SharePolicyGeomObject='share', &
    TransferOfferGeomObject='will provide',  rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.NUOPC_FieldDictionaryHasEntry("surface_downwelling_photosynthetic_radiative_flux", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("surface_downwelling_photosynthetic_radiative_flux", "W m-2 s-1", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="surface_downwelling_photosynthetic_radiative_flux", &
    name="downwelling_short_photosynthetic_radiation_at_water_surface", &
    SharePolicyField='share', SharePolicyGeomObject='share', &
    TransferOfferGeomObject='will provide',  rc=localrc)

  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.NUOPC_FieldDictionaryHasEntry("x_velocity_at_10m_above_sea_surface", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("x_velocity_at_10m_above_sea_surface", "m s-1", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
      StandardName="x_velocity_at_10m_above_sea_surface", &
      name="x_velocity_at_10m_above_sea_surface", &
      SharePolicyField='share', SharePolicyGeomObject='share', &
      TransferOfferGeomObject='will provide',  rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.NUOPC_FieldDictionaryHasEntry("y_velocity_at_10m_above_sea_surface", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("y_velocity_at_10m_above_sea_surface", "m s-1", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="y_velocity_at_10m_above_sea_surface", &
    name="y_velocity_at_10m_above_sea_surface", &
    SharePolicyField='share', SharePolicyGeomObject='share', &
    TransferOfferGeomObject='will provide',  rc=localrc)

  if (.not.NUOPC_FieldDictionaryHasEntry("air_temperature_at_water_surface", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("air_temperature_at_water_surface", "K", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="air_temperature_at_water_surface", &
    name="air_temperature_at_water_surface", &
    SharePolicyField='share', SharePolicyGeomObject='share', &
    TransferOfferGeomObject='will provide',  rc=localrc)

  !> Create fields that feed back from the ocean to the atmosphere.
  !> @todo Make sure that these remain optional in case there is no coupling
  !> @todo check unit for temperature
  if (.not.NUOPC_FieldDictionaryHasEntry("temperature_at_water_surface", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("temperature_at_water_surface", "deg C", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! call NUOPC_Advertise(importState, &
  !   StandardName="sea_surface_temperature", &
  !   name="temperature_at_water_surface", &
  !   SharePolicyField='share', SharePolicyGeomObject='share', &
  !   TransferOfferGeomObject='will provide',  rc=localrc)

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
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  integer(ESMF_KIND_I4)                   :: itemCount

  rc = ESMF_SUCCESS

  ! create a northwestern hemisphere grid
  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/100, 10/), &
    minCornerCoord=(/-150._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
    maxCornerCoord=(/30._ESMF_KIND_R8, 80._ESMF_KIND_R8/), &
    coordSys=ESMF_COORDSYS_SPH_DEG, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  gridOut = gridIn ! for now out same as in

  field = ESMF_FieldCreate(name="temperature_at_water_surface", grid=gridIn, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Disabled for now as we disabled coupling OCN->ATM
  !call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  !> Realize all export fields using the utility function from schism_esmf_util
  call ESMF_StateGet(exportState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(exportState, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    call SCHISM_FieldRealize(exportState, itemNameList(i), &
      grid=gridOut, typeKind=ESMF_TYPEKIND_R8, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo

  !> @todo Provide values for the fields

  call ESMF_StateGet(exportState, field=field, &
    itemName="x_velocity_at_10m_above_sea_surface", rc=localrc)
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

  call ESMF_StateGet(exportState, field=field, &
    itemName="y_velocity_at_10m_above_sea_surface", rc=localrc)
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

  call ESMF_StateGet(exportState, field=field, &
    itemName="air_pressure_at_water_surface", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr2, &
    exclusiveUBound=ubnd, exclusiveLBound=lbnd, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do j=lbnd(2),ubnd(2)
    do i=lbnd(1), ubnd(1)
      farrayPtr2(i,j) = 100000 + 2000*sin(i/ubnd(1)*3.14/180.) &
        + 2000*cos(j/ubnd(2)*3.14/180.)
    end do
  enddo

  call ESMF_StateGet(exportState, field=field, &
    itemName="air_temperature_at_water_surface", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr2, &
    exclusiveUBound=ubnd, exclusiveLBound=lbnd, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do j=lbnd(2),ubnd(2)
    do i=lbnd(1), ubnd(1)
      farrayPtr2(i,j) = 298 + 10*sin(i/ubnd(1)*3.14/180.) &
        - 10*cos(j/ubnd(2)*3.14/180.)
    end do
  enddo

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

end module atmosphere_cmi_nuopc
