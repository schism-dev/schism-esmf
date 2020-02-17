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
    phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeAdvertise, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeRealize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_SetClock, &
    specRoutine=SetClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Do we need a specialization of Finalize, by adding a label?
  !call NUOPC_CompSpecialize(comp, specLabel=model_label_Finalize, &
  !  specRoutine=Finalize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeAdvertise"
!> @description The purpose of this phase is for your model to advertise its import and export fields. This means that your model announces which model variables it is capable of exporting (e.g., an atmosphere might export air pressure at sea level) and which model variables it requires (e.g., an atmosphere might require sea surface temperature as a boundary condition). The reason there is an explicit advertise phase is because NUOPC dynamically matches fields among all the models participating in a coupled simulation during runtime. So, we need to collect the list of possible input and output fields from all the models during their initialization.
! Note that NUOPC does not allocate memory for fields during the advertise phase or when NUOPC_Advertise is called. Instead, this is simply a way for models to communicate the standard names of fields. During a later phase, only those fields that are connected (e.g., a field exported from one model that is imported by another) need to have memory allocated. Also, since ESMF will accept pointers to pre-allocated memory, it is usually not necessary to change how memory is allocated for your model’s variables.
!> Should be mapped to IPDv00p1, IPDv01p1, IPDv02p1, IPDv03p1, IPDv04p1, IPDv05p1;
!> its implementation is required
subroutine InitializeAdvertise(comp, importState, exportState, clock, rc)

  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)       :: localrc, mpiCommunicator, mpiCommDuplicate
  integer(ESMF_KIND_I4)       :: ntime=0, iths=0, i
  character(len=ESMF_MAXSTR)  :: message, compName
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)

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
    StandardName="surface_air_pressure", name="pmsl", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(importState, &
    StandardName="surface_downwelling_photosynthetic_radiative_flux", name="rsns", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo define more export fields, i.e wind x,y fields
  ! NUOPC_AdvertiseField(state, StandardName, Units, &
  !   LongName, ShortName, name, TransferOfferGeomObject, SharePolicyField, &
  !   SharePolicyGeomObject, vm, field, rc)
  call NUOPC_Advertise(exportState, &
    StandardName="surface_temperature", name="temperature_at_water_surface", &
    TransferOfferGeomObject="will provide", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> The mesh information is usually not in CF standard and therefore needs
  !> to be added to the FieldDictionary before advertising
  allocate(itemNameList(4))
  itemNameList=(/ 'mesh_topology                 ', &
                  'mesh_global_node_id           ', &
                  'mesh_global_element_id        ', &
                  'mesh_element_node_connectivity'/)

  do i=1,4
    call NUOPC_FieldDictionaryAddEntry(trim(itemNameList(i)),'1', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Advertise(exportState, StandardName=trim(itemNameList(i)), &
      SharePolicyField='share', SharePolicyGeomObject='not share', &
      TransferOfferGeomObject='will provide', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeRealize"
!> @description During this phase, fields that were previously advertised should now be realized. Realizing a field means that an ESMF_Field object is created and it is added to the appropriate ESMF_State, either import or export. In order to create an ESMF_Field, you’ll first need to create one of the ESMF geometric types, ESMF_Grid, ESMF_Mesh, or ESMF_LocStream. For 2D and 3D logically rectangular grids (such as a lat-lon grid), the typical choice is ESMF_Grid. For unstructured grids, use an ESMF_Mesh. Fields are put into import or export States by calling NUOPC_Realize.
!> Should be mapped to IPDv00p2, IPDv01p3, IPDv02p3, IPDv03p3, IPDv04p3, IPDv05p4;
!> Required if providing any geometry.  For higher phases IPDv03p5, IPDv04p5, IPDv05p6
!> a separate routine should be provided for accepting fields.
subroutine InitializeRealize(comp, importState, exportState, clock, rc)

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

  ! exportable field: surface_temperature
  field = ESMF_FieldCreate(name="temperature_at_water_surface", mesh=mesh2d, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(exportState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetClock"
!> @description A model's clock copies start/stop time and timestep from its
!> parent's clock.  If a model has a timestep that is different (smaller) than
!> the parent's it needs to be set here.
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

  call NUOPC_CompSetClock(comp, externalClock=clock, &
    stabilityTimeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetRunClock"
!> @description If the timestep of the parent is dynamic, then there might
!> be mismatches between intergral timesteps of a model and the parent
!> timestep.  This is reported as an error unless dealt with in this label
!> @todo (not registered and fully implemented yet, so not used)
subroutine SetRunClock(comp, rc)

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

  !> @todo find integral partitioning of timesteps and
  !> set schism's timestep temporarily
  !> NUOPC_AdjustClock(clock, maxTimestep, rc)

  call NUOPC_CompSetClock(comp, clock, timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine



#undef ESMF_METHOD
#define ESMF_METHOD "ModelAdvance"
!> @description As described in the section 4.2, the subroutine ModelAdvance (shown below) has been registered to the special- ization point with the label model_label_Advance in the SetServices subroutine. This specialization point subroutine is called within the generic NUOPC_Model run phase in order to request that your model take a timestep forward. The code to do this is model dependent, so it does not appear in the subroutine below. Each NUOPC component maintains its own clock (an ESMF_Clock object). The clock is used to indicate the current model time and the timestep size. When the subroutine finishes, your model should be moved ahead in time from the current time by one timestep. NUOPC will automatically advance the clock for you, so there is no explicit call to do that here.
!> Because the import/export states and the clock do not come in through the parameter list, they must be accessed via a call to NUOPC_ModelGet
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