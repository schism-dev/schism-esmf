! This code is part of the SCHISM-ESMF interface.  It defines
! the schism component for a NUOPC coupled system
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
#define ESMF_FILENAME "schism_cmi_nuopc.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_cmi_nuopc

  use esmf
  use nuopc
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance

  use schism_bmi
  use schism_esmf_util

  implicit none

  private
  public SetServices

! The internal state saves data across ESMF phases and is
! persistent throught the lifetime of an instance.  Here, we
! only provide a boilerplate implementation of an empty internal state

  type type_InternalStateStruct
  end type

  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type

  character(len=ESMF_MAXSTR), parameter :: &
    label_InternalState = 'InternalState'

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4) :: localrc
  type(type_InternalState)   :: internalState

  rc = ESMF_SUCCESS

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_UserCompSetInternalState(comp, label_InternalState, &
    internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
  !  phaseLabelList=(/"IPDv00p0"/), userRoutine=InitializeP0, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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
#define ESMF_METHOD "InitializeP0"
subroutine InitializeP0(comp, importState, exportState, parentClock, rc)

!  phase 0: Set Initialize Phase Definition Version (REQUIRED, NUOPC PROVIDED)
!  – Initialize the InitializePhaseMap Attribute according to the NUOPC Initialize Phase Definition (IPD) version 00
! ∗ IPDv00p1 = 1: (REQUIRED, IMPLEMENTOR PROVIDED) · Advertise Fields in import and export States.
! ∗ IPDv00p2 = 2: (REQUIRED, IMPLEMENTOR PROVIDED)
! · Realize the advertised Fields in import and export States.
! ∗ IPDv00p3 = 3: (REQUIRED, NUOPC PROVIDED)
! · Check compatibility of the Fields’ Connected status.
! ∗ IPDv00p4 = 4: (REQUIRED, NUOPC PROVIDED)
! · Handle Field data initialization. Time stamp the export Fields.

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: parentClock
  integer, intent(out)  :: rc

  character(len=10)           :: InitializePhaseMap(4)
  integer                     :: localrc
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: configFileName, compName
  character(len=ESMF_MAXSTR)  :: message
  type(ESMF_Config)           :: config

  rc=ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing (p=0) component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  InitializePhaseMap = (/"IPDv00p1=1","IPDv00p2=2", &
    "IPDv00p3=3","IPDv00p4=4"/)

  call ESMF_AttributeAdd(comp, convention="NUOPC", &
    purpose="General", &
    attrList=(/"InitializePhaseMap"/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, name="InitializePhaseMap", valueList=InitializePhaseMap, &
    convention="NUOPC", purpose="General", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Read the configuration for this component from file if not
  ! already present in the component
  call ESMF_GridCompGet(comp, configIsPresent=isPresent, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (isPresent) then
    call ESMF_GridCompGet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A)') trim(compName)//' uses internal configuration'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  else
    configfilename=trim(compName)//'.cfg'
    inquire(file=trim(configfilename), exist=isPresent)

    config = ESMF_ConfigCreate(rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    if (isPresent) then
      call ESMF_ConfigLoadFile(config, trim(configfilename), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      write(message,'(A)')  trim(compName)//' read configuration from '// trim(configFileName)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    else
      write(message,'(A)')  trim(compName)//' has no configuration; use global config'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompSet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

end subroutine InitializeP0

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
  logical                     :: isPresent

  type(ESMF_Vm)               :: vm

  rc = ESMF_SUCCESS

  !> @todo replace by NUOPC_CompGet()
!  NUOPC_GridCompGet(comp, name, verbosity, profiling, diagnostic, rc)
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

  call ESMF_UtilIOMkDir ('./outputs',  relaxedFlag=.true., rc=localrc)
  write(message, '(A)') trim(compName)//' writes results to directory "./outputs".'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(message, '(A)') trim(compName)//' initializing science model ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  inquire(file='param.nml', exist=isPresent)
  if (.not.isPresent) then
    write(message, '(A)') trim(compName)//' cannot find requred file "param.nml"'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_FILE_OPEN
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  call schism_init(0, './', iths, ntime)
  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (.not.NUOPC_FieldDictionaryHasEntry("surface_air_pressure", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("surface_air_pressure", "N m-2", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! call NUOPC_Advertise(importState, &
  !   StandardName="surface_air_pressure", name="air_pressure_at_water_surface", &
  !   SharePolicyField='share', SharePolicyGeomObject='not share', &
  !   TransferOfferGeomObject='will provide',  rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.NUOPC_FieldDictionaryHasEntry("surface_downwelling_photosynthetic_radiative_flux", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("surface_downwelling_photosynthetic_radiative_flux", "W m-2 s-1", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! call NUOPC_Advertise(importState, &
  !   StandardName="surface_downwelling_photosynthetic_radiative_flux", &
  !   name="downwelling_short_photosynthetic_radiation_at_water_surface", &
  !   SharePolicyField='share', SharePolicyGeomObject='not share', &
  !   TransferOfferGeomObject='will provide',  rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.NUOPC_FieldDictionaryHasEntry("surface_temperature", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("surface_temperature", "degree C", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! call NUOPC_Advertise(importState, &
  !   StandardName="surface_temperature", name="air_temperature_at_water_surface", &
  !   SharePolicyField='share', SharePolicyGeomObject='not share', &
  !   TransferOfferGeomObject='will provide',  rc=localrc)
  ! _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(importState, &
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

  !> The mesh information is usually not in CF standard and therefore needs
  !> to be added to the FieldDictionary before advertising
  allocate(itemNameList(4))
  itemNameList=(/ 'mesh_topology                 ', &
                  'mesh_global_node_id           ', &
                  'mesh_global_element_id        ', &
                  'mesh_element_node_connectivity'/)

  do i=1,4
    if (.not.NUOPC_FieldDictionaryHasEntry(trim(itemNameList(i)), rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry(trim(itemNameList(i)),'1', rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    endif
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Advertise(exportState, StandardName=trim(itemNameList(i)), &
      SharePolicyField='share', SharePolicyGeomObject='share', &
      TransferOfferGeomObject='will provide', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

  if (.not.NUOPC_FieldDictionaryHasEntry("sea_surface_temperature", rc=localrc)) then
      call NUOPC_FieldDictionaryAddEntry("sea_surface_temperature", "degree C", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Advertise(exportState, &
    StandardName="sea_surface_temperature", name="temperature_at_water_surface", &
    SharePolicyField='share', SharePolicyGeomObject='not share', &
    TransferOfferGeomObject='will provide',  rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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
  integer(ESMF_KIND_I4)   :: localrc, i

  type(ESMF_CoordSys_Flag) :: coordsys
  type(ESMF_Mesh)          :: mesh2d

  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  integer(ESMF_KIND_I4)                   :: itemCount

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

  !> @todo change variable here
  farrayPtr1 => pr2(1:np)
  field = ESMF_FieldCreate(name="x_velocity_at_10m_above_sea_surface", mesh=mesh2d, &
    farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo Disabled until we fix the coupling
  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtr1 => pr2(1:np)
  field = ESMF_FieldCreate(name="air_pressure_at_water_surface", mesh=mesh2d, &
    farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo Disabled until we fix the coupling
  !call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  field = ESMF_FieldCreate(name="downwelling_short_photosynthetic_radiation_at_water_surface", mesh=mesh2d, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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
      mesh=mesh2d, typeKind=ESMF_TYPEKIND_R8, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  enddo

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

end module schism_cmi_nuopc
