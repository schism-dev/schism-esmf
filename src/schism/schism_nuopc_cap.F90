! This code is part of the SCHISM-ESMF interface.  It defines
! the schism component for a NUOPC coupled system
!
! @copyright (C) 2021-2022 Helmholtz-Zentrum Hereon
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
#define ESMF_FILENAME "schism_nuopc_cap.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_nuopc_cap

  use esmf
  use nuopc
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance

  use schism_bmi
  use schism_esmf_util
  use schism_nuopc_util

  implicit none

  private
  public SetServices

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)             :: localrc

  character(len=ESMF_MAXSTR), parameter :: label_InternalState = 'InternalState'
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_UserCompSetInternalState(comp, label_InternalState, &
    internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! NUOPC automatically has an entry point for InitializeP0, so do not
  ! call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
  !  phaseLabelList=(/"IPDv00p0"/), userRoutine=InitializeP0, rc=localrc)

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

  !call ESMF_AttributeAdd(comp, convention="NUOPC", &
  !  purpose="General", &
  !  attrList=(/"InitializePhaseMap"/), rc=localrc)
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

  type(ESMF_VM)               :: vm
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  !> @todo replace by NUOPC_CompGet()
!  NUOPC_GridCompGet(comp, name, verbosity, profiling, diagnostic, rc)
  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap
  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0
  
  if (.not.ESMF_StateIsCreated(importState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Import', stateintent= &
      ESMF_STATEINTENT_IMPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Import" for import'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  endif

  if (.not.ESMF_StateIsCreated(exportState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Export', stateintent= &
      ESMF_STATEINTENT_EXPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Import" for export'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
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
    write(message, '(A)') trim(compName)//' cannot start without required file "param.nml".'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_FILE_OPEN
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  call schism_init(0, './', iths, ntime)
  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! Deal with setting up the Field dictionary here and advertising the
  ! variables.
  ! @todo this should be customized by a field dictionary-like configuration
  ! file,  for uncopuled applications, we cannot advertise as we get a NUOPC not
  ! connected error message

  !call NUOPC_FieldDictionaryAddIfNeeded("surface_air_pressure", "N m-2", localrc)
  call NUOPC_FieldDictionaryAddIfNeeded("air_pressure_at_sea_level", "N m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("surface_downwelling_photosynthetic_radiative_flux", "W m-2 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("surface_temperature", "degree C", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_merid_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_zonal_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! for coupling to ATMESH
  call NUOPC_Advertise(importState, "air_pressure_at_sea_level", rc=localrc)
  call NUOPC_Advertise(importState, "inst_zonal_wind_height10m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_merid_wind_height10m", rc=localrc)

  ! for coupling to WW3DATA
  call NUOPC_Advertise(importState, "eastward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_northward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "northward_wave_radiation_stress", rc=localrc)

  ! call NUOPC_Advertise(importState, &
  !   StandardName="surface_temperature", name="air_temperature_at_water_surface", &


  !> The mesh information is usually not in CF standard and therefore needs
  !> to be added to the FieldDictionary before advertising
  allocate(itemNameList(4))
  itemNameList=(/ 'mesh_topology                 ', &
                  'mesh_global_node_id           ', &
                  'mesh_global_element_id        ', &
                  'mesh_element_node_connectivity'/)

  do i=1,4
    call NUOPC_FieldDictionaryAddIfNeeded(trim(itemNameList(i)), "1", localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Advertise(exportState, StandardName=trim(itemNameList(i)), &
      SharePolicyField='share', SharePolicyGeomObject='share', &
      TransferOfferGeomObject='will provide', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

  ! @todo the following fails since the canonical unit for SST is K
  !call NUOPC_FieldAdvertise(exportState, "sea_surface_temperature", "degree C", localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeRealize"
!> @description During this phase, fields that were previously advertised should now be realized. Realizing a field means that an ESMF_Field object is created and it is added to the appropriate ESMF_State, either import or export. In order to create an ESMF_Field, you’ll first need to create one of the ESMF geometric types, ESMF_Grid, ESMF_Mesh, or ESMF_LocStream. For 2D and 3D logically rectangular grids (such as a lat-lon grid), the typical choice is ESMF_Grid. For unstructured grids, use an ESMF_Mesh. Fields are put into import or export States by calling NUOPC_Realize.
!> Should be mapped to IPDv00p2, IPDv01p3, IPDv02p3, IPDv03p3, IPDv04p3, IPDv05p4;
!> Required if providing any geometry.  For higher phases IPDv03p5, IPDv04p5, IPDv05p6
!> a separate routine should be provided for accepting fields.
subroutine InitializeRealize(comp, importState, exportState, clock, rc)

!  use schism_esmf_util, only : addSchismMesh
  !> @todo move all use statements of schism into schism_bmi
  use schism_glbl, only: np, pr2, windx2, windy2, srad, nws, rkind
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

  type(ESMF_DistGrid)                :: nodalDistgrid
  type(ESMF_Array)                   :: array

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()
!  real(ESMF_KIND_R8), target, allocatable :: farray1(:)

  !> @todo move to internal state
  real(ESMF_KIND_R8), target, allocatable :: eastward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: eastward_northward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: northward_wave_radiation_stress(:)

  rc = ESMF_SUCCESS
  localrc= ESMF_SUCCESS

  call ESMF_LogWrite("before addSchismMesh",ESMF_LOGMSG_WARNING)
  !call addSchismMesh(comp, ownedNodes=ownedNodes, foreignNodes=foreignNodes, rc=localrc)
  call addSchismMesh(comp, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call ESMF_LogWrite("after addSchismMesh",ESMF_LOGMSG_WARNING)

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

  call ESMF_MeshGet(mesh2d, nodalDistgrid=nodalDistgrid, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="inst_zonal_wind_height10m", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="inst_zonal_wind_height10m", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !farrayPtr1(1:npe) = windx2(1:npe)

  ! The postfix 2 on windx, windy, pr denotes the information at the next
  ! timestep 
!  if(.not.allocated(farray1)) allocate(farray1(np))
!  farray1=windx2(1:np)

!  farrayPtr1 =>windx2(1:np)
!  field = ESMF_FieldCreate(name="inst_zonal_wind_height10m", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, meshloc=ESMF_MESHLOC_NODE, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (NUOPC_IsConnected(importState, "inst_zonal_wind_height10m", rc=localrc) & 
    .and. nws /= 3) then
    call ESMF_LogWrite("Connected zonal wind will not be used if nws /=3", ESMF_LOGMSG_WARNING)
  endif

  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="inst_merid_wind_height10m", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="inst_merid_wind_height10m", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !farrayPtr1(1:np) = windy2(1:np)

!  farrayPtr1 => windy2(1:np)
!  field = ESMF_FieldCreate(name="inst_merid_wind_height10m", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (NUOPC_IsConnected(importState, "inst_merid_wind_height10m", rc=localrc) & 
    .and. nws /= 3) then
    call ESMF_LogWrite("Connected meridional wind will not be used if nws /=3", ESMF_LOGMSG_WARNING)
  endif

  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="air_pressure_at_sea_level", rc=localrc)
  !farrayPtr1 => pr2(1:np)

  !write(0,*) 'farrayptr1 ',ubound(farrayPtr1,1), ubound(pr2,1)
  !array = ESMF_ArrayCreate(nodalDistgrid, farrayPtr=farrayPtr1, &
  !  name="air_pressure_at_sea_level", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="air_pressure_at_sea_level", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !farrayPtr1(1:np) = pr2(1:np)

!  farrayPtr1 => pr2(1:np)
!  field = ESMF_FieldCreate(name="air_pressure_at_sea_level", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (NUOPC_IsConnected(importState, "air_pressure_at_sea_level", rc=localrc) & 
    .and. nws /= 3) then
    call ESMF_LogWrite("Connected air presure will not be used if nws /=3", ESMF_LOGMSG_WARNING)
  endif

#if 0
!@TODO: add other air vars: airt2, shum2, hradd, fluxprc
  farrayPtr1 => srad(1:np)
  field = ESMF_FieldCreate(name="downwelling_short_photosynthetic_radiation_at_water_surface", mesh=mesh2d, &
    typekind=ESMF_TYPEKIND_R8, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (NUOPC_IsConnected(importState, "downwelling_short_photosynthetic_radiation_at_water_surface", rc=localrc) & 
    .and. nws /= 3) then
    call ESMF_LogWrite("Connected downwelling par will not be used if nws /=3", ESMF_LOGMSG_WARNING)
  endif
#endif

  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="eastward_wave_radiation_stress", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="eastward_wave_radiation_stress", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(eastward_wave_radiation_stress(np))
!  farrayPtr1 => eastward_wave_radiation_stress(1:np)
!  field = ESMF_FieldCreate(name="eastward_wave_radiation_stress", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="eastward_northward_wave_radiation_stress", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="eastward_northward_wave_radiation_stress", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(eastward_northward_wave_radiation_stress(np))
!  farrayPtr1 => eastward_northward_wave_radiation_stress(1:np)
!  field = ESMF_FieldCreate(name="eastward_northward_wave_radiation_stress", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
!  allocate(northward_wave_radiation_stress(np))

  array = ESMF_ArrayCreate(nodalDistgrid, typekind=ESMF_TYPEKIND_R8, &
    name="northward_wave_radiation_stress", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  field = ESMF_FieldCreate(name="northward_wave_radiation_stress", mesh=mesh2d, array=array, &
     meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  farrayPtr1 => northward_wave_radiation_stress(1:np)
!  field = ESMF_FieldCreate(name="northward_wave_radiation_stress", mesh=mesh2d, &
!    farrayPtr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_Realize(importState, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


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

  !> As IPDv01p3 (NUOPC Provided) fails when fields are not connected, we here
  !> remove all unconnected fields from the import and export stabilityTimeStep
  call SCHISM_RemoveUnconnectedFields(importState, rc=localrc)
  call SCHISM_RemoveUnconnectedFields(exportState, rc=localrc)

  !@TODO: should we destroy field & array vars?
end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetClock"
!> @description A model's clock copies start/stop time and timestep from its
!> parent's clock.  If a model has a timestep that is different (smaller) than
!> the parent's it needs to be set here.
subroutine SetClock(comp, rc)

  use schism_bmi, only : schismTimeStep
  use schism_glbl, only : wtiminc

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)              :: clock
  type(ESMF_TimeInterval)       :: timeStep
  integer(ESMF_KIND_I4)         :: localrc
  real(ESMF_KIND_R8)            :: seconds
  character(len=ESMF_MAXSTR)    :: message

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, modelClock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Check that wtiminc, i.e. the time between two new atmospheric inputs
  !>  corresponds to the parent (coupling) time step
  call ESMF_ClockGet(clock, timeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_TimeIntervalGet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (abs(wtiminc - seconds) > 1e-6) then 
    write(message, '(A,I7,A,I7)') 'Check setting of wtiminc = ', int(wtiminc), &
      ' in param.nml! Auto-resetting to wtiminc = ', int(seconds)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  wtiminc = seconds

  !> Now get schism's internal timestep and synchronize that with the 
  !> component's time step
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

  use schism_glbl, only: wtiminc,windx2,windy2,pr2,np

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
  integer                     :: npe ! number of exclusive nodes <= np

  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: farrayPtr1(:)
  type(ESMF_StateItem_Flag) :: itemType

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

  !> Obtain radiation tensor from wave component and calculate the wave stress
  call SCHISM_StateImportWaveTensor(importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(importState, itemname='inst_zonal_wind_height10m', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(importState, itemname='inst_zonal_wind_height10m', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_FieldGet(field, farrayptr=farrayPtr1, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! The following assumes that halo exchange on nodes is done inside SCHISM
    ! and we do not have to worry about this here
    npe = ubound(farrayPtr1,1)
    windx2(1:npe)=farrayPtr1(1:npe)
  endif 

  call ESMF_StateGet(importState, itemname='inst_merid_wind_height10m', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(importState, itemname='inst_merid_wind_height10m', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_FieldGet(field, farrayptr=farrayPtr1, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
 
    npe = ubound(farrayPtr1,1)
    windy2(1:npe)=farrayPtr1(1:npe)
  endif 

  call ESMF_StateGet(importState, itemname='air_pressure_at_sea_level', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(importState, itemname='air_pressure_at_sea_level', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_FieldGet(field, farrayptr=farrayPtr1, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    
    npe = ubound(farrayPtr1,1)
    pr2(1:npe)=farrayPtr1(1:npe)
  endif 

!  call ESMF_StateGet(importState, itemname='inst_merid_wind_height10m', field=field, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  call ESMF_FieldGet(field, farrayptr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  windy2(1:np)=farrayPtr1(1:np)

!  call ESMF_StateGet(importState, itemname='air_pressure_at_sea_level', field=field, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  call ESMF_FieldGet(field, farrayptr=farrayPtr1, rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  pr2(1:np)=farrayPtr1(1:np)

  call schism_step(it)
  it = it + 1

    call ESMF_TraceRegionExit("schism:ModelAdvance")
end subroutine

subroutine SCHISM_RemoveUnconnectedFields(state, rc)

  implicit none

  type(ESMF_State), intent(inout) :: state
  integer(kind=ESMF_KIND_I4), optional :: rc

  integer(kind=ESMF_KIND_I4)              :: rc_, localrc, itemCount, i
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  type(ESMF_Field)                        :: field

  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  allocate(itemTypeList(itemCount))
  allocate(itemNameList(itemCount))

  call ESMF_StateGet(state, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    if (.not.NUOPC_IsConnected(state, trim(itemNameList(i)), rc=localrc)) then
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

      call ESMF_StateRemove(state, itemNameList(i:i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    endif

  enddo
end subroutine SCHISM_RemoveUnconnectedFields

#undef  ESMF_METHOD
#define ESMF_METHOD "addSchismMesh"
subroutine addSchismMesh(comp, rc)
! Define ESMF domain partition
  !> @todo apply only filter to 'use schism_glbl'
  use schism_esmf_util, only: type_InternalStateWrapper, type_InternalState
  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  use schism_glbl, only: ylat, xlon, npa, np, nea, ne, ics
  use schism_glbl, only:  nvrt

  implicit none

  type(ESMF_GridComp)  :: comp
  integer(ESMF_KIND_I4), intent(out) :: rc

  type(ESMF_Mesh)          :: mesh2d
  type(ESMF_DistGrid)      :: elementDistgrid, distgrid
  type(ESMF_CoordSys_Flag) :: coordsys
  integer, dimension(:), allocatable            :: nodeids, elementids, nv
  real(ESMF_KIND_R8), dimension(:), allocatable :: nodecoords2d
  real(ESMF_KIND_R8), dimension(:), allocatable :: elementcoords2d
  integer, dimension(:), allocatable            :: nodeowners, elementtypes
  integer, dimension(:), allocatable            :: nodemask, elementmask
  integer, dimension(:), allocatable            :: tmpIdx, tmpIdx2, localNodes, nodeHaloIdx
  integer, dimension(:), allocatable            :: schismTolocalNodes,testids
  integer, dimension(1:4)                       :: elLocalNode
  integer               :: numNodeHaloIdx
  integer               :: i,n,nvcount
  integer               :: ii,ip,ie, localrc
  integer               :: mynp,myne,rank2
  type(llist_type),pointer :: nextp=>null()

  integer(ESMF_KIND_I4), pointer, dimension(:)  :: farrayPtrI41 => null()
  integer(ESMF_KIND_I4), pointer, dimension(:,:):: farrayPtrI42 => null()

  real(ESMF_KIND_R8), parameter :: rad2deg=180.0d0/pi

  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  character(len=ESMF_MAXSTR)           :: compName, message, fieldName
  type(ESMF_State)                     :: exportState
  type(ESMF_Field)                     :: field

  integer(ESMF_KIND_I4) :: numOwnedNodes, numForeignNodes, localPet
  integer(ESMF_KIND_I4) :: ownedCount, foreignCount

  integer(ESMF_KIND_I4), allocatable, target :: ownedNodeIdx(:)
  integer(ESMF_KIND_I4), allocatable, target :: foreignNodeIdx(:)

  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

!  write(0,*)'__LINE__ inside addSchismMesh'
!  call ESMF_Finalize() 

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo all non-ESMF stuff should be outsourced to schism_bmi.F90
  ! prepare mesh
  ! a) take local elements and get number of nodes for definition
  ! b) get number of augmented nodes, which are not connected to local elements
  ! c) allocate node arrays such that nodeids gives all nods belonging to local
  !    elements, all other nodes as part of the augmented domain are outside the
  !    exclusive region
  ! d) allocate element arrays such that local elements (ne) are in array and
  !    augmented elements are defined in the computational domain outside the
  !    exclusive domain
!  numNodeHaloIdx=0
!  npa=npa

!  allocate(localNodes(npa), stat=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  do i=1,npa
!    localNodes(i)=i
!  end do

  ! get inverse matching
!  allocate(schismToLocalNodes(npa))
!  schismToLocalNodes(:) = -1 ! initialize to -1
!  do i=1,npa
!    schismToLocalNodes(localNodes(i))=i
!  end do

  ! define mesh
  !Points to global node #
  allocate(nodeids(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Global coordinates
  allocate(nodecoords2d(2*np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !A node is owned by same rank across PETs; interface nodes are owned by min rank
  allocate(nodeowners(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  nodeowners(:) = -1

  allocate(nodemask(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Points to global elem #
  allocate(elementids(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementtypes(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementmask(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementcoords2d(2*ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! nv (elemConn): 1D array for connectivity (packed from 2D array elnode).
  ! Outputs local node # 
  allocate(nv(sum(i34(1:ne))), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! set ESMF coordSys type
  if (ics==2) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    write(message, '(A)') trim(compName)//' uses a cartesian coordinate system'// &
       ' which may not be suitable for coupling'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    coordsys=ESMF_COORDSYS_CART
  endif

  do ip=1, np
    nodeids(ip)=iplg(ip) !global node #
    if (ics==2) then
      ! if geographical coordinates present
      nodecoords2d(2*ip-1) = rad2deg*xlon(ip)
      nodecoords2d(2*ip)   = rad2deg*ylat(ip)
    else
      ! use cartesian coordinates
      nodecoords2d(2*ip-1) = xnd(ip)
      nodecoords2d(2*ip)   = ynd(ip)
    end if

    !nodeowners must be unique cross all PETs
    rank2=ipgl(iplg(ip))%rank
    nodeowners(ip) =rank2 !init 
    if(associated(ipgl(iplg(ip))%next)) then !interface or ghost node
      if(ipgl(iplg(ip))%next%rank<rank2) then
        nodeowners(ip) =ipgl(iplg(ip))%next%rank
      endif
    endif

!    if (ip<=np) then
!      ! if iplg() is resident, ipgl(iplg(i))%rank=myrank
!      nodeowners(ip) = ipgl(iplg(ip))%rank
!    else !not executed at the moment
!      ! get owner of foreign node, the SCHISM manual tells us to not use
!      ! ipgl for foreign nodes ???
!      ! ipgl%next%next%next.... is the linked list, with ranks in ascending order.  Unless ipgb is an interface node (i.e., resident in more than 1 process), the list has only 1 entry
!      nextp => ipgl(iplg(ip))
!      ! We here advance to the end of the list, Joseph suggested to just take the next element
!      do while (associated(nextp%next))
!        nextp => nextp%next
!      end do
!      nodeowners(ip) = nextp%rank
!    end if
    nodemask(ip)         = idry(ip)
  end do !ip

  ! As the list of owned and non-owned nodes is not preserved in the ESMF_Mesh
  ! structure, we need to save this information to an internal state, for later 
  ! use in Array/Field creation.
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0
  do ip=1,np
    if (nodeowners(ip) == localPet) then 
      isDataPtr%numOwnedNodes = isDataPtr%numOwnedNodes + 1
    else 
      isDataPtr%numForeignNodes = isDataPtr%numForeignNodes + 1
    endif
  enddo

  if (isDataPtr%numForeignNodes + isDataPtr%numOwnedNodes /= np) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
      'mismatching number of resident np=',np,', owned=',isDataPtr%numOwnedNodes, &
      ' and foreign=', isDataPtr%numForeignNodes,' nodes'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  write(message,'(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
    'matching number of resident np=',np,', owned=',isDataPtr%numOwnedNodes, &
    ' and foreign=', isDataPtr%numForeignNodes,' nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  
  allocate(isDataPtr%ownedNodeIds(isDataPtr%numOwnedNodes), stat=localrc)
  allocate(isDataPtr%foreignNodeIds(isDataPtr%numForeignNodes), stat=localrc)

  ownedCount = 0
  foreignCount = 0
  do ip=1,np
    write(0,*) np, ip, ownedCount, foreignCount!, nodeowners(ip)
    if (nodeowners(ip) == localPet) then
      ownedCount=ownedCount + 1
      isDataPtr%ownedNodeIds(ownedCount) = ip
    else
      foreignCount=foreignCount + 1
      isDataPtr%foreignNodeIds(foreignCount) = ip
    endif
  enddo

  write(message, '(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
    'number of resident np=',np,', owned=',isDataPtr%numOwnedNodes, &
    ' and foreign=', isDataPtr%numForeignNodes,' nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(*,*) 'Owned nodes on PET ',localPet,isDataPtr%ownedNodeIds
  write(*,*) 'Foreign nodes on PET ',localPet,isDataPtr%foreignNodeIds

  nvcount=0
  do i=1,ne
    elementids(i)=ielg(i)
    if(i34(i)==3) then
      elementtypes(i)=ESMF_MESHELEMTYPE_TRI
    else
      elementtypes(i)=ESMF_MESHELEMTYPE_QUAD
    endif
    elementmask(i)=idry_e(i)
    do ii=1,i34(i)
      elLocalNode(ii)=elnode(ii,i)
      nvcount = nvcount+1
      nv(nvcount) =elnode(ii,i) 
    end do

    ! Element coord is the sum of nodes divided by number of nodes
    elementcoords2d(2*i-1)=sum(nodecoords2d(2*elLocalNode(1:i34(i))-1))/i34(i)
    elementcoords2d(2*i)=sum(nodecoords2d(2*elLocalNode(1:i34(i))))/i34(i)
  end do !i

  if(ubound(nv,1)/=nvcount) then
    localrc=ESMF_RC_ARG_SIZE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc) 
  endif

#if 0
  write(0,*) 'Local Nodes, np=',np
  do i=1,npa
   write(0,*) 'localid, globalid, nodeowner',i,nodeids(i),nodeowners(i)
  end do
  nvcount=1
  do i=1,ne
    write(0,*) 'ne, conn',i,nv(nvcount:nvcount+i34(i)-1)
    nvcount = nvcount+i34(i)
  end do
#endif

  ! create element distgrid (distribute)
  elementDistgrid = ESMF_DistgridCreate(elementids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  mesh2d = ESMF_MeshCreate(parametricDim=2,spatialdim=2,nodeIds=nodeids, &
             nodeCoords=nodecoords2d,nodeOwners=nodeowners, &
             coordSys=coordsys, &
             nodeMask=nodemask,elementMask=elementmask, &
             elementIds=elementids, elementTypes=elementtypes, &
             elementCoords=elementcoords2d, &
             elementDistgrid=elementDistgrid, &
             elementConn=nv, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshGet(mesh2d, numOwnedNodes=mynp, numOwnedElements=myne, elementDistgrid=distgrid, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(testids(myne))
  call ESMF_DistGridGet(distgrid,localDE=0,seqIndexList=testids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,I3.3,A,I3.3,A)') trim(compName)//' created mesh from "', np, &
    'resident nodes and ', myne, ' resident elements in SCHISM'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  write(message, '(A,I3.3,A,I3.3,A)') trim(compName)//' created mesh with "', mynp, &
    'owned nodes and ', myne, ' owned elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  !> @todo the following might overflow the message buffer easily ...
  !write(message,*) 'elementIds:',elementIds
  !call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  !write(message,*) 'distgridElementIds:',testids
  !call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  deallocate(testids)

  call ESMF_GridCompSet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo the following steps don't work in the NUOPC cap yet

  !> Create fields for export to describe mesh (this information is not yet
  !> accessible with ESMF_MeshGet calls)
  !> @todo remove this part of the code once there is a suitable ESMF implementation

  !> Create a dummy field to satisfy ugrid conventions
  field = ESMF_FieldEmptyCreate(name='mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'cf_role', 'mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'topology_dimension', 2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'node_coordinates', &
   'mesh_node_lon mesh_node_lat', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'face_node_connectivity', 'mesh_element_node_connectivity', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  fieldName = 'mesh_global_node_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_NODE, typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtrI41 = nodeIds(1:np)

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  fieldName = 'mesh_global_element_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_ELEMENT, typeKind=ESMF_TYPEKIND_I4, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtrI41 = elementIds(1:ne)

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  nullify(farrayPtrI41)

  fieldName = 'mesh_element_node_connectivity'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
    meshloc=ESMF_MESHLOC_ELEMENT, ungriddedLBound=(/1/), ungriddedUBound=(/4/), &
    typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI42, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1,ne
    do n=1,i34(i)
      farrayPtrI42(i,n) = iplg(elnode(n,i))
    end do
  end do

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  nullify(farrayPtrI42)

  ! clean up
  deallocate(nodeids, stat=localrc)
  deallocate(nodecoords2d, stat=localrc)
  deallocate(nodeowners, stat=localrc)
  deallocate(nodemask, stat=localrc)
  deallocate(elementids, stat=localrc)
  deallocate(elementmask, stat=localrc)
  deallocate(elementtypes, stat=localrc)
  if (allocated(elementCoords2d)) deallocate(elementCoords2d, stat=localrc)
  deallocate(nv, stat=localrc)

  write(message, '(A)') trim(compName)//' created 2D mesh"'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

end subroutine addSchismMesh

end module schism_nuopc_cap
